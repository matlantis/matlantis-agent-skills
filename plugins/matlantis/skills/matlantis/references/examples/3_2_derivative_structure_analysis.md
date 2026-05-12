# 指定した母構造からの置換構造探索の解析

Matlantis CSP で行った置換構造探索の解析を行います。


```python
%pip install -U mtcsp[matlantis]
%pip install pandas plotly
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain

from pathlib import Path

from ase.visualize import view
from pandas import DataFrame
from plotly import express as px
from pymatgen.analysis.phase_diagram import PDPlotter

from mtcsp.analysis import MTCSPEntry
from mtcsp.analysis import MTCSPPhaseDiagram
from mtcsp.analysis import select_min_energy_crystal_structures_by_composition
from mtcsp.atoms import FileSystemAtomsStore
from mtcsp.experiment import FrozenExperiment
```


## Experimentの読み込み

解析用の Experiment である FrozenExperiment を読み込みます。
Experiment にはメタデータのみが入っているため結晶構造の実体を取り出すための AtomsStore も作成します。
`db_file` と `atoms_store_dir` は探索時の設定と同じものを設定してください。
ここでは、ディレクトリ直下の最新のdb_fileに対応するexperimentを読み込みます。


```python
db_file = sorted([path.name for path in Path(".").glob("derivative_structure_search_*.journal")])[
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


## Energy above Hull 計算

安定性を評価するためには convex hull からのエネルギー差 energy above hull を確認することが重要となり、そのためにはまず hull を定義する必要があります。hull を定義するには以下の複数の方法があります。 
1. Materials Project のようなデータベースを利用する
2. 自身で所有しているデータを入力する。
3. Matlantis CSP の Convex hull search を利用して探索する。

1 や 2 の方法はよく知られている系や、知見の蓄積がある系で利用できます。
あまり性質が知られていない系では 3 を利用するとより良い解析が可能です。

Materials Project を利用するには API Key を取得し、設定する必要があります。
取得した API Key は、pmg コマンドで登録してください。(cf. [pymatgen](https://pymatgen.org/usage.html#setting-the-pmg_mapi_key-in-the-config-file )) 
設定後に notebook の kernel の restart が必要な場合があります。 

以下のセルをアンコメントすると探索した系で Materials Project の hull 上の構造をダウンロードして置換探索と同条件で構造最適化の後エネルギーを評価することが出来ます。


```python
"""
from tqdm import tqdm

from mtcsp.localopt import Relaxer
from mtcsp.materials_project import download_mp_atoms


elements = experiment.elements
atoms_list, _, id_list = download_mp_atoms(elements=elements, max_atoms=100, max_e_above_hull=1e-8)

mp_entries = []
relaxer = experiment.get_relaxer()
for atoms, id in zip(tqdm(atoms_list, desc="Optimizing MP atoms"), id_list):
    atoms_opt, e_pot = relaxer(atoms)
    mp_entries.append(MTCSPEntry(atoms.get_chemical_formula(), e_pot, id))
"""
```


ここでは API Key をまだ設定していない場合のためにこの系に関してあらかじめ取得した結果を使います。


```python
mp_hull_data = [
    {"composition": "Cu3Au1", "energy": -14.602007030317928, "id": "mp-2258"},
    {"composition": "Cu2Au2", "energy": -14.153859548053688, "id": "mp-522"},
]
mp_entries = [
    MTCSPEntry(**data)  # ty: ignore[invalid-argument-type]
    for data in mp_hull_data
]
```


比較対象となる convex hull を描画します。


```python
phase_diagram = MTCSPPhaseDiagram(
    mp_entries, elements=experiment.elements, reference_entries=experiment.reference_entries
)
# e_above_hullが`show_unstable`(eV/atom 単位)以下の構造を表示
plotter = PDPlotter(phase_diagram, show_unstable=0.01)
plotter.show()
```


比較対象となる convex hull が用意できたので、探索結果の情報を取得します。
Matlantis CSP では一つの構造のメタデータを CrystalStructure という class にまとめて保持しています。
全メタデータの list は `Experiment.completed_crystal_structure` で取得できます。
探索中の全構造を対象に解析するのは数が多すぎるのでここでは各組成で最小のエネルギーを持つ CrystalStructure を
スクリーニングする関数 `select_min_energy_crystal_structures_by_composition` を利用します。
また、 `Experiment.get_crystal_structure_info_by_id_list` を利用すると結晶構造の安定性、対称性、Cell の情報などが表として取得できます。
オプション変数である `phase_diagram` を渡すことで Energy above hull を含めた表にすることが出来ます。


```python
cs_list = experiment.completed_crystal_structures
cs_list = [cs for cs in cs_list if len(set(cs.chemical_symbols)) &gt; 1]  # remove simple structures
min_energy_cs_list = select_min_energy_crystal_structures_by_composition(cs_list)
min_energy_cs_id_list = [cs.id for cs in min_energy_cs_list]
```


```python
data = experiment.get_crystal_structure_info_by_id_list(
    min_energy_cs_id_list, phase_diagram=phase_diagram
)

df = DataFrame(data).fillna(0).sort_values("energy_above_hull_eV_per_atom", ignore_index=True)
df.head(5)
# df.to_csv("hull_cs_list.csv", index=False)
```


## 結晶構造の可視化

表や相図から可視化したい CrystalStructure の ID を取得してその構造を可視化します。
cif や xyz 形式に出力して他の解析に利用することも出来ます。


```python
cs_id = df.id[0]
atoms = experiment.get_atoms_by_id_list([cs_id])[0]
atoms.wrap()
view(atoms.repeat(2), viewer="x3d")

# import ase.io
# ase.io.write(f"{atoms.get_chemical_formula()}.cif", atoms)
```


## 置換した組成空間での相安定性の可視化

この例ではCu サイトを CuもしくはAuで置換しているので、これらの比率を軸として可視化し、有望な構造を探すことが出来ます。
まずはそのために各構造に対して、Cu, Au (`x_Cu`, `x_Au`) を計算します。
この計算方法は母構造や置換方法によって異なるため以下を参考にケースに合わせて計算してください。

以下の例では元のCuサイトは unit cell あたり1箇所であることからそれぞれの比率を計算しています。
また、各 CrystalStructure が unit cell を何倍して作られているかは `CrystalStructure.cell_size` で取得できます。


```python
data = []
num_Cu_in_unit_cell = 1
for cs in min_energy_cs_list:
    assert cs.cell_size is not None
    assert cs.chemical_symbols is not None
    n_Cu = cs.chemical_symbols.count("Cu") or 0
    n_Au = cs.chemical_symbols.count("Au") or 0
    x_Cu = n_Cu / (cs.cell_size * num_Cu_in_unit_cell)
    x_Au = n_Au / (cs.cell_size * num_Cu_in_unit_cell)
    data.append(
        {
            "id": cs.id,
            "x_Cu": x_Cu,
            "x_Au": x_Au,
        }
    )

df_substitution = DataFrame(data).fillna(0)

df = df.merge(df_substitution, on="id")
```


`x_Cu`, `x_Au`は足して1になる変数であることから、以下のように可視化が出来ます。
図の点にマウスオーバーすると CrystalStructure の ID が表示されます。
有望そうな構造があったら ID を指定して結晶構造を取得できます。


```python
fig = px.scatter(
    df,
    x="x_Cu",
    y="energy_above_hull_eV_per_atom",
    hover_name="id",
    hover_data=["formula", "formation_energy_eV_per_atom", "spacegroup_symbol"],
    labels={
        "x_Cu": "Ratio of Cu",
        "energy_above_hull_eV_per_atom": "Energy above hull&lt;br&gt;[eV/atom]",
    },
)
fig.update_traces(marker=dict(size=12))
fig.update_layout(
    title="Energy above hull of the min-energy crystal structures",
    paper_bgcolor="white",
    width=768,
    height=512,
)
fig.show()
```


```python
from mtcsp.crystal_structure import CrystalStructureID


cs_id = CrystalStructureID(129)
atoms = experiment.get_atoms_by_id_list([cs_id])[0]
atoms.wrap()
view(atoms.repeat(2), viewer="x3d")

# import ase.io
# ase.io.write(f"{atoms.get_chemical_formula()}.cif", atoms)
```


## Experimentのエクスポート

探索結果を保存用に出力することが出来ます。

以下の2つのファイルが出力されます。
- &lt;experiment_name&gt;.sqlite: Experiment に含まれるのメタデータ情報が保存された sqlite ファイル
- &lt;experiment_name&gt;.zip: 結晶構造の情報が含まれるファイル


```python
from mtcsp.export import export_experiment


export_experiment(experiment_name, db_file, atoms_store)
```


以下のように export したファイルから FrozenExperiment を読み出すことが出来ます。


```python
from mtcsp.atoms import ZipFileAtomsStore


atoms_store_zipfile = experiment_name + ".zip"
experiment_db_file = experiment_name + ".sqlite"
atoms_store = ZipFileAtomsStore(atoms_store_zipfile)
experiment = FrozenExperiment(
    experiment_name=experiment_name,
    storage=experiment_db_file,
    atoms_store=atoms_store,
)
```
