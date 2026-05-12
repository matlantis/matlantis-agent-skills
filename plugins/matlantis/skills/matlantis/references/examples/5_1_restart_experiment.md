# 指定した元素系全体の結晶構造探索(hull search)の再開

このノートブックでは、一度実行していた探索を再開する方法を紹介します。


```python
%pip install -U mtcsp[matlantis]
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain


from mtcsp.atoms import FileSystemAtomsStore
from mtcsp.convex_hull_search import initialize_hull_search
from mtcsp.convex_hull_search import perform_hull_search
from mtcsp.convex_hull_search import SearchConfig
from mtcsp.convex_hull_search import SystemConfig
from mtcsp.experiment import FrozenExperiment
```


## Experimentの作成

探索の再開の例のために、まずは探索を実行します。


```python
experiment_name = "csp-Ti-O-examle"

search_config = SearchConfig(
    elements=("Ti", "O"),
    experiment_name=experiment_name,
    # CSP search,
    population_size=4,
    max_atoms=8,
    add_mp_crystals_to_initial_population=False,
    # PFP
    pfp_version=None,  # Use default PFP version
)
system_config = SystemConfig.from_experiment_name(experiment_name)
```


```python
initialize_hull_search(search_config=search_config, system_config=system_config)
```


```python
first_n_crystal_structures = 8

perform_hull_search(
    search_config=search_config,
    system_config=system_config,
    n_crystal_structures=first_n_crystal_structures,
    parallelism=6,
)
```


## 探索の再実行

探索の再開方法を説明するために、まずmtcspが探索結果を内部的に保存している方式を説明します。
下に概要図を示します。
mtcspでは[Optuna](https://optuna.org/)を用いて、探索結果をデータベース上に保存しています。
`experiment_name="csp-Ti-O-example"`のように探索名を指定すると、対応するテーブルデータが`db_file`で指定したデータベース上に保存されます。
データベース上のエントリは`mtcsp.crystal_structure.CrystalStructure`に対応しており、組成式や形成エネルギー、結晶構造を保存しているファイルへのパスなどが保存されています。
具体的な結晶構造はJSON形式で`atoms_store_dir`で指定したディレクトリ以下に保存されています。

![mtcsp_data_layout](./assets/5_1_restart_experiment/mtcsp_data_layout.png)


`SystemConfig` (`db_file`と`atoms_store_dir`)と`experiment_name`を指定することで、探索結果を復元することができます。
探索を再開するためには、加えて探索アルゴリズムのパラメータ(`SearchConfig`)を指定する必要があります。
以上から、`SystemConfig`と`SearchConfig`を指定すれば、Experimentを再読み込みして探索を再開することができます。

ここでは、探索する結晶構造の累計を増やして追加で探索を行う例を示します。
既に`first_n_crystal_structure = 8`個の結晶構造を探索しているので、追加で`second_n_crystal_structures - first_n_crystal_structures = 12`個の結晶構造を探索、同一のデータベース(`db_file`)に保存します。


```python
second_n_crystal_structures = 20

perform_hull_search(
    search_config=search_config,
    system_config=system_config,
    n_crystal_structures=second_n_crystal_structures,
    parallelism=6,
)
```


探索結果を読み込むと、探索完了した構造`FrozenExperiment.completed_crystal_structures`の数が、後から指定した探索数分だけ増えていることが確認できます。


```python
experiment = FrozenExperiment(
    experiment_name=search_config.experiment_name,
    storage=str(system_config.db_file),
    atoms_store=FileSystemAtomsStore(base_path=system_config.atoms_store_dir),
)
```


```python
assert len(experiment.completed_crystal_structures) &gt;= second_n_crystal_structures + len(
    search_config.elements
)
```
