# 指定した元素系全体の結晶構造探索の解析

Matlantis CSP で行った元素系全体の結晶構造探索の解析を行います。


```python
%pip install -U mtcsp[matlantis]
%pip install pandas
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain

from pathlib import Path

from ase.visualize import view
from pandas import DataFrame
from pymatgen.analysis.phase_diagram import PDPlotter

from mtcsp.analysis import MTCSPEntry
from mtcsp.analysis import MTCSPPhaseDiagram
from mtcsp.atoms import FileSystemAtomsStore
from mtcsp.experiment import FrozenExperiment
```


## Experimentの読み込み

解析用の Experiment である FrozenExperiment を読み込みます。
Experiment にはメタデータのみが入っているため結晶構造の実体を取り出すための AtomsStore も作成します。
`db_file` と `atoms_store_dir` は探索時の設定と同じものを設定してください。
ここでは、ディレクトリ直下の最新のdb_fileに対応するexperimentを読み込みます。


```python
db_file = sorted([path.name for path in Path(".").glob("convex_hull_search_*.journal")])[
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
crystal_structures = experiment.completed_crystal_structures
print(f"Number of completed crystal structures: {len(crystal_structures)}")
```


## 状態図に基づく解析

探索結果の可視化、解析には pymatgen を利用します。
`MTCSPEntry.from_crystal_structure` を利用すると pymatgen の [PDEntry](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.phase_diagram.PDEntry) 互換のオブジェクトとして探索結果が取得できます。
描画機能の詳細は [pymatgen](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.phase_diagram.PhaseDiagram) を参照してください。
相図上の構造をマウスオーバーすると今回の探索における構造の ID が表示されます。


```python
csp_entries = [MTCSPEntry.from_crystal_structure(cs) for cs in crystal_structures]
phase_diagram = MTCSPPhaseDiagram(
    csp_entries,
    elements=experiment.elements,
    reference_entries=experiment.reference_entries,
)
# Display structures with e_above_hull less than `show_unstable` (in eV/atom)
plotter = PDPlotter(phase_diagram, show_unstable=0.01)
plotter.show()
```


### エネルギー凸包上の結晶構造

Matlantis CSP では一つの構造のメタデータを CrystalStructure という class にまとめて保持しています。
Convex hull 上に位置する CrystalStructure を取り出し、それらの結晶構造情報を取り出すには以下のようにします。
また、探索中のすべての CrystalStructure は `experiment.completed_crystal_structures` で取得出来ます。


```python
hull_cs_list = experiment.get_final_stable_crystal_structures()
hull_cs_id_list = [cs.id for cs in hull_cs_list]
```


Convex hull 上の構造の情報を表にして表示します。
`FrozenExperiment.get_crystal_structure_info_by_id_list` に CrystalStructure の ID の list を与えると Cell や対称性の情報が集約され、pandas の DataFrame に変換して表示・解析することが出来ます。
ここでは生成エネルギーでソートしています。


```python
data = experiment.get_crystal_structure_info_by_id_list(hull_cs_id_list)

df = DataFrame(data).sort_values("formation_energy_eV_per_atom", ignore_index=True)
df
# df.to_csv("hull_cs_list.csv", index=False)
```


### 結晶構造の可視化

表や相図から可視化したい CrystalStructure の ID を取得してその構造を可視化します。
実際の結晶構造はメタデータである CrystalStructure class とは別に保存されています。
`FrozenExperiment.get_atoms_by_id_list` に ID の list を与えて ase.Atoms の list を取得します。
cif などの形式で保存すれば別のソフトウェアでの解析に用いることも出来ます。


```python
cs_id = df.id[0]
atoms = experiment.get_atoms_by_id_list([cs_id])[0]
view(atoms.repeat(2), viewer="x3d")


# import ase.io
# ase.io.write("atoms.cif", atoms)
```


### Materials Projectの凸包との比較

CSP 探索結果を Materials Project と比較します。
Materials Project の API Key を取得し、pmg コマンドで登録してください（cf. [pymatgen](https://pymatgen.org/usage.html#setting-the-pmg_mapi_key-in-the-config-file ))。 
設定後に notebook の kernel の restart が必要な場合があります。 


以下のコメントアウトされているセルでは Materials Project の hull 上の構造をダウンロードして置換探索と同条件で構造最適化の後エネルギーを評価します。
ここでは API Key をまだ設定していない場合のために Sr-Ti-O 系に関してあらかじめ取得した結果を使います。


```python
# from tqdm import tqdm
#
# from mtcsp.localopt import Relaxer
# from mtcsp.materials_project import download_mp_atoms
#
#
# elements = experiment.elements
# atoms_list, _, id_list = download_mp_atoms(elements=elements, max_atoms=100, max_e_above_hull=1e-8)
#
# mp_entries = []
# relaxer = Relaxer()
# for atoms, id in zip(tqdm(atoms_list, desc="Optimizing MP atoms"), id_list):
#     atoms_opt, e_pot = relaxer(atoms)
#     mp_entries.append(MTCSPEntry(atoms.get_chemical_formula(), e_pot, id))
```


```python
mp_hull_data = [
    {"composition": "OSr", "energy": -12.737429328359836, "id": "mp-2472"},
    {"composition": "O2Sr", "energy": -18.12319359096691, "id": "mp-2697"},
    {"composition": "O4Sr2Ti", "energy": -54.98943117557418, "id": "mp-5532"},
    {"composition": "O7Sr3Ti2", "energy": -97.08787997210938, "id": "mp-3349"},
    {"composition": "O40Sr2Ti22", "energy": -591.9580541421976, "id": "mp-28740"},
    {"composition": "O6Sr2Ti2", "energy": -84.11562420837957, "id": "mp-4651"},
    {"composition": "OTi2", "energy": -26.659255059617834, "id": "mp-1215"},
    {"composition": "O6Ti4", "energy": -93.85556433286406, "id": "mp-458"},
    {"composition": "O4Ti12", "energy": -138.6395509621912, "id": "mp-2591"},
    {"composition": "O10Ti6", "energy": -150.42471796982232, "id": "mp-8057"},
    {"composition": "O2Ti12", "energy": -116.16734132231868, "id": "mp-554098"},
    {"composition": "O3Ti3", "energy": -55.04663688850377, "id": "mp-1071163"},
    {"composition": "O4Ti2", "energy": -56.40386631620808, "id": "mp-390"},
]


mp_entries = [
    MTCSPEntry(**data)  # ty: ignore[invalid-argument-type]
    for data in mp_hull_data
]
```


```python
# Add entries for elemental phases.
simple_entries = [
    MTCSPEntry.from_crystal_structure(cs)
    for cs in experiment.completed_crystal_structures[: len(experiment.elements)]
]
mp_entries = mp_entries + simple_entries
```


PFP でエネルギーを評価した場合の Materials Project の相図を描画します。


```python
phase_diagram = MTCSPPhaseDiagram(
    mp_entries, elements=experiment.elements, reference_entries=experiment.reference_entries
)
# e_above_hullが`show_unstable`(eV/atom 単位)以下の構造を表示
plotter = PDPlotter(phase_diagram, show_unstable=0.01)
plotter.show()
```


`FrozenExperiment.get_crystal_structure_info_by_id_list` PhaseDiagram のオブジェクトを渡すと Energy above Hull (eV/atom) も計算し、列を追加できます。
Materials Project の PhaseDiagram オブジェクト渡せば Materials Project の hull と比較した場合のエネルギーを計算することが出来ます。


```python
data = experiment.get_crystal_structure_info_by_id_list(
    hull_cs_id_list, phase_diagram=phase_diagram
)

df = DataFrame(data).sort_values("energy_above_hull_eV_per_atom", ignore_index=True)
df
# df.to_csv("hull_cs_list.csv", index=False)
```


### 探索経過の表示

世代ごとのエネルギー凸包の体積をプロットします。
横軸は世代数(1世代あたりに`population_size`個の構造を含む)を、縦軸はエネルギー凸包の体積(eV/atom単位)を表しています。
これによって、探索数(`n_crystal_structures`)に対して系の探索がどの程度進んでいるかを確認できます。
エネルギー凸包の体積がほぼ一定になって変化しなくなれば、探索数が十分であると判断できます。


```python
from mtcsp.visualization import plot_hull_hypervolume_history_by_generation
```


```python
fig = plot_hull_hypervolume_history_by_generation(experiments=[experiment])
fig
```


## 凸包付近の構造の取得と結晶構造の解析


```python
from mtcsp.analysis import get_structures_and_e_above_hull
from mtcsp.analysis import unique_computed_structure_entries
from mtcsp.export import export_as_computed_structure_entry
```


得られた凸包から測ってenergy above hull が`max_e_above_hull`以下の構造を取得し、結晶構造の解析を行います。


```python
near_hull_data = get_structures_and_e_above_hull(experiment, max_e_above_hull=0.03)
for entry, e_above_hull in near_hull_data:
    print(
        f"ID {entry.crystal_structure_id}: {entry.composition} - E above hull: {e_above_hull:.4f} eV/atom"
    )
```


PFPで評価したエネルギーと結晶構造をpymatgenの`ComputedStructureEntry`オブジェクトに変換します。


```python
near_hull_crystal_structures = experiment.get_crystal_structure_by_id_list(
    [entry.crystal_structure_id for entry, _ in near_hull_data]
)

near_hull_computed_structure_entries = []
for cs in near_hull_crystal_structures:
    cse = export_as_computed_structure_entry(cs, experiment.atoms_store)
    if cse is None:
        continue
    near_hull_computed_structure_entries.append(cse)
```


`unique_computed_structure_entries`で重複する結晶構造を除去します。


```python
uniqued_near_hull_computed_structure_entries = unique_computed_structure_entries(
    near_hull_computed_structure_entries
)
```


Materials Project API を設定している場合は、以下をアンコメントすることで得られた構造をMaterials Project掲載構造と照合することができます。


```python
# from mtcsp.analysis import MaterialsProjectMatcher
# mpm = MaterialsProjectMatcher(elements=experiment.elements)
#
# for cse in uniqued_near_hull_computed_structure_entries:
#     result = mpm.fit_structure(cse.structure)
#     if result is not None:
#         print(f"{cse.composition} matched with {result}")
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
