# クイックスタート

このクイックスタートチュートリアルでは、MTCSPの基本的な使用方法を説明し、凸包探索の実行と結果の分析を行います。
このチュートリアルでは、利用可能なすべてのオプションと機能について説明していません。詳細については、後続のチュートリアルとドキュメントを参照してください。

## インストール

Matlantis CSPで探索を行うには、まず`mtcsp`パッケージをインストールします。


```python
%pip install -U mtcsp[matlantis]
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain

from datetime import datetime
from pathlib import Path

from pymatgen.analysis.phase_diagram import PDPlotter

from mtcsp.analysis import MTCSPEntry
from mtcsp.analysis import MTCSPPhaseDiagram
from mtcsp.atoms import FileSystemAtomsStore
from mtcsp.convex_hull_search import initialize_hull_search
from mtcsp.convex_hull_search import perform_hull_search
from mtcsp.convex_hull_search import SearchConfig
from mtcsp.convex_hull_search import SystemConfig
from mtcsp.experiment import FrozenExperiment
```


## 探索パラメータの設定

まず、探索の設定(元素系の指定を含む)を行います。


```python
elements = ["Cu", "Au"]

experiment_name = f"quickstart_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
search_config = SearchConfig(
    elements=elements,
    experiment_name=experiment_name,
    population_size=8,  # Small value for quickstart
)
system_config = SystemConfig(
    db_file=Path(f"{experiment_name}.journal"),
    atoms_store_dir=Path("atoms_store"),
    log_file=None,
)
```


## 探索の初期化と実行

定義した設定を使用して探索の初期化を行います。


```python
initialize_hull_search(search_config=search_config, system_config=system_config)
```


初期化した探索を実行して、安定な結晶構造を探索します。
このセルの実行には数分かかります。


```python
perform_hull_search(
    search_config=search_config,
    system_config=system_config,
    n_crystal_structures=16,  # Small value for quickstart
    parallelism=6,
)
```


## 探索結果の解析

探索が完了したら、まず先程設定した `db_file` と `atoms_store_dir` から結果(`mtcsp`では`FrozenExperiment`と呼ばれます)を読み込みます。


```python
experiment = FrozenExperiment(
    experiment_name=search_config.experiment_name,
    storage=str(system_config.db_file),
    atoms_store=FileSystemAtomsStore(system_config.atoms_store_dir),
)
entries = [MTCSPEntry.from_crystal_structure(cs) for cs in experiment.completed_crystal_structures]
print(f"Number of completed crystal structures: {len(entries)}")
```


その後、探索結果から相図を構築し、可視化することができます。


```python
phase_diagram = MTCSPPhaseDiagram(
    entries,
    elements=experiment.elements,
    reference_entries=experiment.reference_entries,
)
plotter = PDPlotter(phase_diagram)
plotter.show()
```
