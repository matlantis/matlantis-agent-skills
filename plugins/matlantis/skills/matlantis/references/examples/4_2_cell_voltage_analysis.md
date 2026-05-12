# 置換構造探索を用いた平衡電位計算の解析

置換構造探索の結果を利用してLiCoO2の平衡電位を計算します。


```python
%pip install -U mtcsp[matlantis]
%pip install pandas
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain

from pathlib import Path

import ase.io
import ase.visualize
import pandas as pd
import plotly.express as px
from pymatgen.analysis.phase_diagram import CompoundPhaseDiagram
from pymatgen.core import Composition

from mtcsp.analysis import MTCSPEntry
from mtcsp.analysis import select_min_energy_crystal_structures_by_composition
from mtcsp.atoms import FileSystemAtomsStore
from mtcsp.crystal_structure import CrystalStructureID
from mtcsp.experiment import FrozenExperiment
```


## 探索結果の読み込み

解析用の Experiment である FrozenExperiment を読み込みます。
Experiment にはメタデータのみが入っているため結晶構造の実体を取り出すための AtomsStore も作成します。
`db_file` と `atoms_store_dir` は探索時の設定と同じものを設定してください。
ここでは、ディレクトリ直下の最新のdb_fileに対応するexperimentを読み込みます。


```python
db_file = sorted([path.name for path in Path(".").glob("cell_voltage_*.journal")])[
    -1
]  # Choose the latest journal file
experiment_name = Path(db_file).stem  # Get the experiment name from the journal file
print(f"{experiment_name=}")

atoms_store_dir = "atoms_store"
```


```python
atoms_store = FileSystemAtomsStore(atoms_store_dir)
experiment = FrozenExperiment(
    experiment_name=experiment_name,
    storage=db_file,
    atoms_store=atoms_store,
)
```


Matlantis CSP では一つの構造のメタデータを CrystalStructure というクラスにまとめて保持しています。
全メタデータは `Experiment.completed_crystal_structure` で取得できます。
探索中の全構造を対象に解析するのは数が多すぎるので、ここでは各組成で最小のエネルギーを持つ CrystalStructure を
スクリーニングする関数 `select_min_energy_crystal_structures_by_composition` を利用します。
また、 `Experiment.get_crystal_structure_info_by_id_list` を利用すると結晶構造の安定性、対称性、Cell の情報などが表として取得できます。


```python
cs_list = [
    cs for cs in experiment.completed_crystal_structures if len(set(cs.chemical_symbols)) &gt; 1
]  # remove simple structures
min_energy_cs_list = select_min_energy_crystal_structures_by_composition(cs_list)
```


```python
data = experiment.get_crystal_structure_info_by_id_list(
    id_list=[cs.id for cs in min_energy_cs_list],
)
df = pd.DataFrame(data).fillna(0).set_index("id")
df.head(4)
```


## 平衡電位計算

置換探索の結果をもとに以下の反応式の平衡電位を計算します。

$
\mathrm{Li}_{x_1}\mathrm{Co}\mathrm{O}_{2} \rightarrow \mathrm{Li}_{x_2}\mathrm{Co}\mathrm{O}_{2} + (x_1 - x_2) \mathrm{Li}
\quad (\mathrm{with}\, (x_1 - x_2) e^-)
$

まず、各構造に対してLi サイトにおける空孔の比率 (`x_vacancy`) を計算します。
この計算方法は母構造や置換方法によって異なるため以下を参考にケースに合わせて計算してください。
なお、CrystalStructure が unit cell を何倍して作られているかは `CrystalStructure.cell_size` で取得できます。


```python
base_composition = Composition("LiCoO2")

num_Li_in_unit_cell = base_composition.get_el_amt_dict()["Li"]

data = []
for cs in min_energy_cs_list:
    assert cs.chemical_symbols is not None
    assert cs.cell_size is not None
    n_Li = cs.chemical_symbols.count("Li") or 0
    x_Li = n_Li / cs.cell_size
    x_vacancy = num_Li_in_unit_cell - x_Li
    entry = MTCSPEntry.from_crystal_structure(cs)

    data.append(
        {
            "id": cs.id,
            "x_vacancy": x_vacancy,
            "cell_size": cs.cell_size,
            "entry": entry,
            # Number of atoms in the formula, that is, 4-x for Li_{1-x}CoO2
            "num_atoms_in_formula": entry.composition.num_atoms / cs.cell_size,
        }
    )

df = df.merge(pd.DataFrame(data), on="id").set_index("id")
```


次に、置換構造だけのエネルギー凸包を考えたときに、各構造が凸包上の安定な構造であるかどうかを判定しておきます。


```python
# Endpoints for derivative-structure-only convex hull
endpoint_x_vacancy_0: MTCSPEntry = df[df["x_vacancy"] == 0].iloc[0]["entry"]
endpoint_x_vacancy_1: MTCSPEntry = df[df["x_vacancy"] == 1].iloc[0]["entry"]

derivative_structure_phase_diagram = CompoundPhaseDiagram(
    entries=[entry for entry in df["entry"]],
    terminal_compositions=[endpoint_x_vacancy_0.composition, endpoint_x_vacancy_1.composition],
)

pseudo_stable = {}
for i, (index, row) in enumerate(df.iterrows()):
    e_above_hull = derivative_structure_phase_diagram.get_e_above_hull(
        derivative_structure_phase_diagram.entries[i]
    )
    assert e_above_hull is not None
    pseudo_stable[index] = e_above_hull &lt; 1e-8
df["pseudo_stable"] = pseudo_stable
```


```python
fig = derivative_structure_phase_diagram.get_plot()
fig.update_layout(
    title="Compound Phase Diagram only for derivative structures",
    paper_bgcolor="white",
)
fig
```


置換構造間の形成エネルギー差から平衡電位を計算します。
原子あたりの形成エネルギーを組成あたりの形成エネルギーに変換する必要がある点に注意ください。


```python
df_pseudo_stable = df[df["pseudo_stable"]].sort_values("x_vacancy").reset_index()

data = []
for i, (_, row) in enumerate(df_pseudo_stable.iterrows()):
    if i == len(df_pseudo_stable) - 1:
        continue
    next_row = df_pseudo_stable.iloc[i + 1]

    voltage = (
        next_row["formation_energy_eV_per_atom"] * next_row["num_atoms_in_formula"]
        - row["formation_energy_eV_per_atom"] * row["num_atoms_in_formula"]
    ) / (next_row["x_vacancy"] - row["x_vacancy"])
    for r in [row, next_row]:
        data.append(
            {
                "x_vacancy": r["x_vacancy"],
                "voltage_eV": voltage,
                "formula": r["formula"],
                "id": r["id"],
            }
        )

df_voltage = pd.DataFrame(data)
```


空孔の比率と平衡電位をプロットします。
図の点にマウスオーバーすると CrystalStructure の ID が表示されます。
このID を指定して結晶構造を取得できます。


```python
fig = px.line(
    df_voltage,
    x="x_vacancy",
    y="voltage_eV",
    hover_name="id",
    hover_data=["formula"],
    labels={
        "voltage_eV": "Voltage (V)",
        "x_vacancy": "x in Li_{1-x}CoO2",
    },
    markers=True,
)
fig.update_traces(
    marker=dict(size=12),
)
fig.update_layout(
    title="Equilibrium voltage of LiCoO2",
    yaxis_range=(0, 5),
)
fig
```


```python
cs_id = CrystalStructureID(83)
atoms = experiment.get_atoms_by_id_list([cs_id])[0]
ase.visualize.view(atoms.repeat(2), viewer="x3d")
```
