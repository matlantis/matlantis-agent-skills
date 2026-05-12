# 既存の探索結果を初期構造に追加する

このノートブックでは、既存の探索結果を初期構造に追加して探索を新たに行う方法を紹介します。


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

例示のために、まずは探索を実行します。


```python
presearch_experiment_name = "csp-Ti-O-presearch"

presearch_config = SearchConfig(
    elements=("Ti", "O"),
    experiment_name=presearch_experiment_name,
    # CSP search,
    population_size=2,
    max_atoms=8,
    add_mp_crystals_to_initial_population=False,
    # PFP
    pfp_version=None,  # Use default PFP version
)
presearch_system_config = SystemConfig.from_experiment_name(presearch_experiment_name)
```


```python
initialize_hull_search(search_config=presearch_config, system_config=presearch_system_config)
```


```python
perform_hull_search(
    search_config=presearch_config,
    system_config=presearch_system_config,
    n_crystal_structures=4,
)
```


## 初期構造に探索結果を追加したExperimentを作成する

まず、上で実行した探索を`FrozenExperiment`として読み込みます。


```python
presearch_experiment = FrozenExperiment(
    experiment_name=presearch_experiment_name,
    storage=str(presearch_system_config.db_file),
    atoms_store=FileSystemAtomsStore(presearch_system_config.atoms_store_dir),
)
```


新たに実行する`Experiment`を作成します。


```python
experiment_name = "csp-Sr-Ti-O-examle"

search_config = SearchConfig(
    elements=("Sr", "Ti", "O"),
    experiment_name=experiment_name,
    # CSP search,
    population_size=2,
    max_atoms=8,
    add_mp_crystals_to_initial_population=False,
    # PFP
    pfp_version=None,
)
system_config = SystemConfig.from_experiment_name(experiment_name)
```


```python
experiment = initialize_hull_search(search_config=search_config, system_config=system_config)
```


初期化後の`Experiment`の`add_crystal_structures_from_experiment`メソッドを用いて指定したenergy above hull 以下の構造を追加します。　


```python
experiment.add_crystal_structures_from_experiment(
    other_experiment=presearch_experiment,
    e_above_hull=0.1,  # in eV/atom
)
```


初期構造に追加した後は、通常通り探索を行うことができます。


```python
perform_hull_search(
    search_config=search_config,
    system_config=system_config,
    n_crystal_structures=4,
)
```


# 初期構造を直接Experimentに追加する

`Experiment.add_unoptimized_atoms`を用いて`ase.Atoms`のリストを直接`Experiment`に追加することもできます。

```python
list_atoms: list[ase.Atoms]
experiment.add_unoptimized_atoms(list_atoms)
```
