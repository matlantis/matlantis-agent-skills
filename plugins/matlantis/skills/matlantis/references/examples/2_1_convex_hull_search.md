# 指定した元素系全体の結晶構造探索

この notebook では元素系を指定して相図全体を探索し、データベースにはない未知の安定結晶候補を探す例を示します。
長時間の探索を行う場合には backgroud notebook として実行してください。

## インストール
Matlantis CSP では Matlantis 上で提供されている mtcsp という python library を用いて探索を行うため下記コマンドでインストールを行います。
Python 環境のバージョンは python 3.13 以上を用いてください。


```python
%pip install -U mtcsp[matlantis]
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain

from datetime import datetime
from pathlib import Path

from mtcsp.convex_hull_search import initialize_hull_search
from mtcsp.convex_hull_search import perform_hull_search
from mtcsp.convex_hull_search import SearchConfig
from mtcsp.convex_hull_search import SystemConfig
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
elements = ("Sr", "Ti", "O")
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


experiment_name = start_or_continue(elements, prefix="convex_hull_search")

# Alternatively, enter custom_experiment_name to override
# experiment_name = "my_experiment_name"
```


`SystemConfig.from_experiment_name`クラスメソッドを用いると、探索名から上記の設定を自動で行うことができます。


```python
system_config = SystemConfig.from_experiment_name(experiment_name)
```


### 探索のパラメータ設定

探索のパラメータを設定します。
設定すべきパラメータは元素系に関するパラメータ、探索アルゴリズムに関するパラメータなどです。
特に探索に合わせてユーザーが変更する部分に関して以下で説明します。

### `SearchConfig`の代表的なパラメータとその説明

- `elements`
    - 探索の対象となる元素系。無償版でサポートされているのは3元素系までです。
- `population_size`
    - 遺伝的アルゴリズムの1世代で生成する結晶構造の数。多様性を確保するために64-128程度に設定するのが望ましいです。
- `max_atoms`
    - 結晶構造をランダム生成する際の最大原子数。32-64 程度に設定するのが標準的です。
- `add_mp_crystals_to_initial_population`
    - 初期 population に Materials Project の構造を追加し、効率的に探索することが出来ます。
    - Materials Project の API Key を[こちら](https://next-gen.materialsproject.org/api)で取得し、pmg コマンドで登録する必要があります（cf. [pymatgen](https://pymatgen.org/usage.html#setting-the-pmg_mapi_key-in-the-config-file ))。 設定後に notebook の kernel の restart が必要な場合があります。
- `pfp_version`
    - PFPのバージョン。指定しない場合、pfp-api-clientのデフォルトが使われます。
- `pfp_compatibility`
    - 組成空間全体で統一的に`EstimatorCalcMode`を選ぶスキーム。
      各オプションでの振る舞いはドキュメントのFundamental-&gt;About MTCSPを参照してください。
- `seed`
    - 乱数の seed 値。指定した場合は同じseed値での探索を再現することが出来ます。

**以下の例ではテスト用に各パラメータを小さく設定しています。プロダクトランでは表にしたがってパラメータを設定してください。**


```python
search_config = SearchConfig(
    elements=elements,
    experiment_name=experiment_name,
    # CSP search,
    population_size=16,
    max_atoms=32,
    add_mp_crystals_to_initial_population=False,
    # PFP,
    pfp_version="v8.0.0",
)
```


### 探索の初期化

Matlantis CSP では一つの探索は Experiment と呼びます。
Experiment にはシステム、アルゴリズム、探索設定などの情報が集約されています。

Experiment の初期化を実行します。


```python
initialize_hull_search(search_config=search_config, system_config=system_config)
```


### 探索の実行

探索する構造数の総数を`n_crystal_structures`で指定します。
以下ではテスト用に小さい値を指定していますが、プロダクトランでは10000 以上に設定するのが望ましく、3元素系では40000程度が推奨されます。

探索の並列数を`parallelism`で指定します。
通常notebookでは 6、2x notebook では 12、4x notebookでは24程度に設定してください。
通常notebookの場合、40000 構造をサンプルする探索は通常2-3日程度かかります。
notebookスペックに対して並列数を増やしすぎるとCPUが常に最大まで使われて処理が滞り、かえって探索が遅くなる場合があることにご注意ください。

`SystemConfig.log_file`を指定している場合、探索の進捗ログは指定されたファイルに追記されます。


```python
perform_hull_search(
    search_config=search_config,
    system_config=system_config,
    n_crystal_structures=50,
    parallelism=6,
)
```


以上で探索は終了です。
探索結果の結果の解析は`2_2_convex_hull_analysis.ipynb`で行います。
