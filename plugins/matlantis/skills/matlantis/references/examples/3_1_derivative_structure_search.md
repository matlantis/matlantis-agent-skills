# 指定した母構造からの置換構造探索

この notebook では指定した母構造から置換構造を生成して探索を行う例を示します。
長時間の探索を行う場合には backgroud notebook として実行してください。

## インストール

Matlantis CSP では Matlantis 上で提供されている mtcsp という python library を用いて探索を行うため下記コマンドでインストールを行います。


```python
%pip install -U mtcsp[matlantis]
%pip install pandas
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain

from datetime import datetime
from pathlib import Path

import ase
import ase.build
import ase.data
from ase.data import atomic_numbers
from ase.visualize import view
from pandas import DataFrame

from mtcsp.analysis import standardize_atoms
from mtcsp.convex_hull_search import SystemConfig
from mtcsp.derivative_structure import append_derivative_structures
from mtcsp.derivative_structure import DerivativeStructureSearchConfig
from mtcsp.derivative_structure import initialize_derivative_structure_search
from mtcsp.derivative_structure import list_num_derivative_structures
from mtcsp.derivative_structure import perform_derivative_structure_search
from mtcsp.derivative_structure import prepare_base_structure_from_atoms
from mtcsp.localopt import Compatibility
```


## 探索

### 探索のシステム設定

探索名(`experiment_name`)を決め、以下を設定した`SystemConfig`クラスを作成します。
- `SystemConfig.db_file`: メタデータを保存するデータベースファイル、拡張子は`.journal`としてください。
- `SystemConfig.atoms_store_dir`: 原子構造の情報を保存するディレクトリ
- `SystemConfig.log_file`: ログファイル、`None`を指定すると標準エラー出力にログが出力されます。

探索名が同じ場合、同じファイルを使用しようとするため、異なる実験を行いたい場合には探索名を使い回さないよう注意してください。

以下ではユーティリティ関数`start_or_continue`で指定した`elements`とタイムスタンプから探索名を生成しています。
この関数は既存の探索があれば自動でそれを選択し、なければ新規に作成します。


```python
symbols = ["Cu"]
```


```python
def start_or_continue(elements: tuple[str, ...], prefix: str):
    experiment_base_name = prefix + "_" + "-".join(sorted(elements))
    files = list(p for p in Path.cwd().glob("*.journal"))
    files_filtered = list(filter(lambda p: p.is_file() and experiment_base_name in p.name, files))
    files_sorted = sorted(files_filtered, key=lambda p: p.stat().st_mtime)

    # Create new experiment name if no previous experiments were found
    if len(files_sorted) &lt; 1:
        experiment_name = experiment_base_name + "_" + datetime.now().strftime("%Y%m%d%H%M%S")
        print(f"Starting new experiment '{experiment_name}'")
    else:
        # Take most recently modified experiment.
        experiment_name = files_sorted[0].stem
        print(f"Continuing from '{experiment_name}'")
    return experiment_name


experiment_name = start_or_continue(
    elements=tuple(sorted(set(symbols))), prefix="derivative_structure_search"
)

# Alternatively, enter custom_experiment_name to override
# experiment_name = "my_experiment_name"
```


`SystemConfig.from_experiment_name`クラスメソッドを用いると、探索名から上記の設定を自動で行うことができます。


```python
system_config = SystemConfig.from_experiment_name(experiment_name)
```


### 母構造(base structure)の設定

母構造の設定を行います。
この例では fcc Cu を母構造として設定します。
cif や xyz などのファイルから読み込んで設定することも可能です。
置換探索では母構造を primitive cell で取る必要があります。
この例ではあえて repeat してAtomsを作成し、primitive cell が取り直されている様子を可視化しています。


```python
parent_atoms = ase.build.bulk("".join(symbols), crystalstructure="fcc", cubic=True)
view(parent_atoms, viewer="x3d")
```


```python
primitive_atoms = standardize_atoms(parent_atoms, primitive=True)
view(primitive_atoms, viewer="x3d")
```


### 置換元素の設定

元素の置換方法を設定します。
置換を行わない元素も含めて母構造に存在するすべての元素について記述してください。
空孔を設定したい場合は `VACANCY_NUMBER` を指定してください。


```python
feasible_number_mapping = {
    atomic_numbers["Cu"]: [
        atomic_numbers["Cu"],
        atomic_numbers["Au"],
    ],
}

search_config = DerivativeStructureSearchConfig.from_feasible_number_mapping(
    feasible_number_mapping=feasible_number_mapping,
    experiment_name=experiment_name,
    pfp_version="v8.0.0",
    pfp_compatibility=Compatibility.PFP2023COMPATIBILITY,
)
```


Primitive cell の何倍までの supercell を取るかを示す `max_index` 設定します。
同じN倍 cell でも repeat させ方は様々な方法がありますが、Matlantis CSP ではこれらで表現される結晶構造を重複無く列挙することが出来ます。
(cf. Fig.2 of [Hart and Forcade 2008, Phys. Rev. B 77, 224115](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.224115))
`max_index` を大きくすると置換構造数が指数関数的に大きくなってしまうので設定には注意が必要です。
以下では index をいくつにすると合計の構造数 `total_num_structures` がいくつになるかを示しています。
目安としては10万構造以下になるように設定すると現実的な時間(1x notebook で１週間程度、4x notebook で数日程度)で探索が終了します。

置換構造数が膨大となってしまう場合には、[組成とスーパーセルを限定した置換構造探索](3_3_derivative_structure_search_with_composition_and_supercell)に示す、組成とスーパーセルを限定して置換構造をランダムサンプリングする方法を検討してください。


```python
# Please provide the composition as an integer ratio for each element as shown below.
# composition = {
#     atomic_numbers['Cu']: 1,
#     atomic_numbers['Au']: 1,
# }

base_structure = prepare_base_structure_from_atoms(
    atoms=primitive_atoms,
    feasible_number_mapping=feasible_number_mapping,
    # composition=composition,
)
data = list_num_derivative_structures(base_structure=base_structure)
df = DataFrame(data).set_index("index")
df
```


`total_num_structures` をみて max_index を設定します。
この例では数分で探索が終了するように200構造程度を評価するように設定します。


```python
max_index = 6
```


### 探索の初期化

Matlantis CSP では一つの探索は Experiment と呼びます。
Experiment にはシステム、アルゴリズム、探索設定などの情報が集約されています。

Experiment の初期化を実行します。
ここで構造を列挙し追加していくので構造数が多い場合には少し時間がかかります。


```python
experiment = initialize_derivative_structure_search(
    search_config=search_config,
    system_config=system_config,
)
append_derivative_structures(
    experiment=experiment,
    base_structure=base_structure,
    max_index=max_index,
)
```


### 探索の実行

探索を行います。
並列で探索するプロセス数を`parallelism`で指定します。
通常notebookでは 6、2x notebook では 12、4x notebookでは24程度に設定してください。

`SystemConfig.log_file`を指定している場合、探索の進捗ログは指定されたファイルに追記されます。


```python
perform_derivative_structure_search(
    search_config=search_config,
    system_config=system_config,
    parallelism=6,
)
```


以上で探索は終了です。
探索結果の結果の解析は`3_2_derivative_structure_analysis.ipynb`で行います。
