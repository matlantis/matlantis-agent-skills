---
name: mt-setup
description: >
  Calculator初期化、PFP Estimator設定、計算モード(calc_mode)選択、構造ファイルI/O、
  バッチ実行設定を扱うセットアップスキルです。
  pfp_api_client, ASECalculator, Estimator, EstimatorCalcMode, estimator_fn, pfp_estimator_fn,
  max_retries, calc_mode選択, PBE, PBE_U, R2SCAN, WB97XD, D3,
  ase.io.read, ase.io.write, .cif, .xyz, POSCAR, .traj,
  run_jobs, ResourceAwareJobScheduler, calculator.reset,
  RetriesExceeded, MultiCalculatorUseDetected, ConcurrentUseDetected,
  PFPAPIError, Too many atoms, Too many neighbors
  に関するコード生成時に使用してください。
---
# Calculator 初期化・構造 I/O

## 概要

Matlantis でシミュレーションを行うための基盤となるガイドです。PFP (PreFerred Potential) の Calculator 初期化、計算モード (calc_mode) の選択、原子構造ファイルの入出力、バッチ実行の設定、および共通エラーへの対処を包括的に扱います。

PFP は PFN が開発したグラフニューラルネットワークベースの汎用原子間ポテンシャルで、DFT に近い精度を桁違いに短い時間で実現します。Matlantis 上のすべてのシミュレーションは、Estimator と ASECalculator の初期化から始まります。

## Environment Preflight

Python を使う前に、必ず Python が使えるか確認してください。

```bash
python --version > /dev/null 2>&1
```

- `python --version > /dev/null 2>&1` が成功する場合: 通常どおりセットアップを進めてください。
- `python --version > /dev/null 2>&1` の終了コードが 1 の場合: その場で停止し、`Please exit Claude, run use_venv, and relaunch.` と表示して `use_venv` 実行を推奨してください。

上記の条件を満たさない状態では、Python 実行や Python パッケージ前提の手順に進まないでください。

## Notebook Kernel Policy

Matlantis 向けに `.ipynb` を新規作成・更新する場合は、Notebook metadata の `kernelspec` に **既存の kernel 名** を設定してください。`python3` のような汎用名を新規に書き込んではいけません。background job で `kernel 'python3' is not available` になりやすいためです。

- 既知の既存 kernel がある場合: その名前を `metadata.kernelspec.name` に設定してください（例: `python313`, `python311`）
- どの kernel を使うか明確な場合: `display_name` も対応する既存 kernel に揃えてください
- 既存 Notebook を編集する場合: 元の `metadata.kernelspec` を保持してください
- 新規 Notebook を作る場合: 利用可能な既存 kernel を選んで metadata に明示してください

## ワークフロー

```
1. Import     - ase, pfp_api_client 等の必須ライブラリをインポート
2. Select     - 計算対象に合わせて calc_mode を選択
3. Initialize - Estimator / ASECalculator を生成
4. Load       - 構造ファイルを読み込み ase.Atoms オブジェクトを生成
5. Attach     - atoms.calc = calculator で計算準備を完了
```

### 計算モード選択フローチャート

```
計算対象は？
├── 有機分子・単分子系
│   ├── 9 元素のみ（H C N O F P S Cl Br） → WB97XD
│   └── それ以外の元素を含む → R2SCAN_PLUS_D3
├── 遷移金属酸化物・強相関系 → PBE_U（分散力が重要なら PBE_U_PLUS_D3）
├── 分子結晶・液体・吸着系
│   ├── 実験値との定量比較が重要 → R2SCAN_PLUS_D3
│   └── 定性的・スクリーニング用途 → PBE_PLUS_D3
├── バルク結晶の形成エネルギー
│   ├── 実験値と比較 → R2SCAN
│   └── Materials Project と比較 → PBE（または PBE_U）
└── 迷ったとき → PBE
```

## 実装パターン

### パターン A: Estimator / ASECalculator の基本初期化

最も基本的な初期化パターンです。`model_version` を常に明示して再現性を確保します。

```python
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

# Estimator を作成（モデルバージョン・計算モードを指定）
estimator = Estimator(
    model_version="v8.0.0",
    calc_mode=EstimatorCalcMode.PBE,
    max_retries=15  # サーバー混雑時のリトライ（10-15 推奨）
)

# ASE Calculator にラップ
calculator = ASECalculator(estimator)

# Atoms オブジェクトに設定
atoms.calc = calculator
```

リトライ付きファクトリ関数として定義する場合:

```python
def get_calculator(calc_mode: str = "PBE", max_retries: int = 15) -> ASECalculator:
    """PFP Calculator を初期化する。"""
    estimator = Estimator(
        calc_mode=calc_mode,
        model_version="v8.0.0",
        max_retries=max_retries
    )
    calculator = ASECalculator(estimator)
    return calculator
```

### パターン B: matlantis-features 用の estimator_fn ファクトリパターン

matlantis-features の各 Feature には `estimator_fn`（Estimator を毎回新規生成するファクトリ関数）を渡します。Estimator インスタンスを直接渡してはいけません。

```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode

# 組み込みファクトリ
estimator_fn = pfp_estimator_fn(
    model_version="v8.0.0",
    calc_mode=EstimatorCalcMode.PBE,
)

# またはカスタムファクトリ
def my_estimator_fn():
    return Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.R2SCAN)
```

### パターン C: エネルギー・力・応力の取得

```python
from ase.build import bulk

atoms = bulk("Si")
atoms.calc = calculator

energy = atoms.get_potential_energy()  # eV
forces = atoms.get_forces()            # eV/A
stress = atoms.get_stress()            # eV/A^3（周期系のみ）
```

stress は周期境界条件 (`pbc=True`) が設定された系でのみ意味を持ちます。

### パターン D: 構造ファイルの読み書き

ASE の `read` / `write` は拡張子からフォーマットを自動判別します。

```python
from ase.io import read, write

# 読み込み（最後のフレーム）
atoms = read("structure.cif")

# 読み込み（最初のフレーム）
atoms = read("structure.cif", index=0)

# 全フレーム読み込み（トラジェクトリ等）
frames = read("opt.traj", index=":")

# 書き出し（拡張子で自動判定）
write("output.cif", atoms)
write("output.xyz", atoms)
write("POSCAR", atoms, format="vasp")
```

対応フォーマット:

| 拡張子 | 形式 | 備考 |
|--------|------|------|
| `.cif` | Crystallographic Information File | 結晶構造の標準形式 |
| `.xyz` | XYZ | シンプルな座標形式 |
| `.vasp` / `POSCAR` | VASP | VASP 入出力形式 |
| `.traj` | ASE Trajectory | エネルギー・力・応力も保持するバイナリ形式 |

### パターン E: バッチ実行と資源管理

`pfcc_extras.job_scheduler` を使って Notebook のバッチ実行を管理します。

```python
from pfcc_extras.job_scheduler.runner import run_jobs
```

実行形態の使い分け:

| スケジューラ | 用途 |
|-------------|------|
| `ResourceAwareJobScheduler` | 並列投入（資源制約を監視しながら） |
| `QueueScheduler` | 逐次完了待ち（フラット配列は 1 件ずつ、二重配列はグループ単位） |
| `PapermillScheduler` | パラメータ掃引（Papermill 連携） |

資源制約パラメータ:

| パラメータ | 説明 |
|-----------|------|
| `my_limit` | 個人トークン率の上限 |
| `tenant_limit` | テナントトークン率の上限 |
| `mem_limit_mb` | メモリ使用量の上限 (MB) |
| `cpu_limit_percent` | CPU 使用率の上限 (%) |
| `post_submit_wait_sec` | 投入直後の待機時間（過投入防止） |

### パターン F: .traj 変換と外部解析ツール連携

ASE の `.traj` ファイルを MDAnalysis / MDTraj で解析するための変換:

```python
from pfcc_extras.structure.ase_traj_converter import convert_traj
```

### パターン G: 入力サイズ制限の確認

計算後に `calc_stats` を確認することで、入力構造がサーバー制限にどの程度近いかを把握できます。

```python
atoms.calc = calculator
energy = atoms.get_potential_energy()
stats = atoms.calc.results['calc_stats']
```

### パターン H: in-place 修正後の Calculator リセット

Atoms オブジェクトを in-place で変更した場合、Calculator のキャッシュをリセットする必要があります。

```python
atoms.positions[0] += [0.1, 0.0, 0.0]  # in-place 修正
calculator.reset()  # キャッシュをクリア
energy = atoms.get_potential_energy()  # 再計算
```

### パターン I: サポート元素の確認

モデルバージョン・計算モードごとに計算可能な元素を確認します。

```python
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator(model_version="v8.0.0", calc_mode="PBE")
elements = estimator.supported_elements(model_version="v8.0.0", calc_mode="PBE")
print(f"Supported elements: {elements}")
print(f"Number of supported elements: {len(elements)}")
```

サポート元素はモードによって異なります。特に WB97XD は H, C, N, O, F, P, S, Cl, Br の 9 元素のみです。

### パターン J: 可視化による構造確認

計算投入前に `pfcc_extras.show_gui` で構造を確認します。

```python
from pfcc_extras import show_gui

# 単一構造の可視化
show_gui(atoms)

# トラジェクトリ（構造のリスト）の可視化
frames = read("opt.traj", index=":")
show_gui(frames)
```

周期境界がおかしい、原子が重なっている、意図しない構造になっているといった問題を可視化で即座に発見できます。

### パターン K: 出力ディレクトリの退避

既存成果物の上書き事故を防ぐため、出力ディレクトリが既に存在する場合はタイムスタンプを付けて退避します。

```python
import os
import shutil
from datetime import datetime

def ensure_output_dir(path: str) -> str:
    """出力ディレクトリを安全に作成する。既存の場合はタイムスタンプ付きで退避。"""
    if os.path.exists(path):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = f"{path}_backup_{timestamp}"
        shutil.move(path, backup)
        print(f"Existing directory moved to: {backup}")
    os.makedirs(path, exist_ok=True)
    return path
```

## ベストプラクティス

### 計算モード一覧

| モード | 汎関数 | 対応元素数 | 主な用途 |
|--------|--------|-----------|---------|
| `PBE` | GGA-PBE | 72 | 汎用（デフォルト） |
| `PBE_U` | GGA-PBE + Hubbard U | 72 | 遷移金属酸化物・強相関系 |
| `PBE_PLUS_D3` | GGA-PBE + D3 分散補正 | 72 | 分子結晶・液体（定性） |
| `PBE_U_PLUS_D3` | GGA-PBE + U + D3 | 72 | 酸化物 + 有機配位子 |
| `R2SCAN` | meta-GGA R2SCAN | 72 | バルク形成エネルギー |
| `R2SCAN_PLUS_D3` | meta-GGA R2SCAN + D3 | 72 | 分子・吸着（定量） |
| `WB97XD` | hybrid wB97X-D | 9 | 有機分子（単分子系専用） |

### 重要な注意事項

1. **model_version は常に明示する**: 指定しないと最新版が使われ、将来的に計算結果が変動するリスクがあります。`model_version="v8.0.0"` のようにバージョンを固定してください。

2. **Estimator は共有しない**: 1 つの Estimator を複数の ASECalculator で共有してはいけません。Calculator ごとに個別の Estimator を生成してください。

3. **異なるモード間のエネルギー比較は不可**: 各モードは汎関数が異なるためエネルギー基準点の絶対値が一致しません。エネルギー差を議論する際は同一モードで揃えてください。

4. **WB97XD は単分子系専用**: 液体・固体など分子間相互作用が重要な系では精度が大幅低下します。

5. **MD の各ステップでトークンを消費**: 長時間 MD の前に短いテストランでトークン消費量を見積もってください。

6. **PBC の設定**: 結晶・スラブ系では `atoms.pbc = True` を設定してください。設定を忘れると物理的に無意味な結果になります。

7. **可視化による確認**: 計算投入前に構造を確認してください。周期境界の誤りや原子の重なりは可視化で即座に発見できます。

8. **Platform Recovery**: 環境更新後に不安定な場合は、Dashboard の `Restart Application` または `Force Stop Application` で復旧してください。

## よくあるエラーと対処

### pfp-api-client エラー

| エラー | 原因 | 対処 |
|--------|------|------|
| `RuntimeError: Atoms object has no calculator` | Calculator が未設定のまま `get_potential_energy` 等を呼び出した | `atoms.calc = calculator` を実行してください |
| `RetriesExceeded` | サーバーアクセス集中でタイムアウト | `max_retries=15` 程度に増やしてください |
| `ConcurrentUseDetected` | 1 つの ASECalculator を複数スレッドで使い回した | 各スレッドで個別に ASECalculator を生成してください |
| `MultiCalculatorUseDetected` | 1 つの Estimator を複数の ASECalculator で共有した | Calculator ごとに Estimator を新規作成してください |

### PFPAPIError 詳細

| メッセージ | 原因 | 対処 |
|-----------|------|------|
| `Atoms are far away from the primitive cell` | 原子とセルの距離が遠すぎる | セルと原子座標を確認してください |
| `Cell is too small` | セルが PBC に対して小さすぎる | スーパーセル化してセルサイズを大きくしてください |
| `Illegal atomic number was detected` | サポートされていない元素 | モード別サポート元素を確認してください |
| `Internal server error` | サーバー側の予期せぬエラー | 時間をおいて再試行し、継続する場合はサポートに連絡してください |
| `Operation incomplete due to timeout` | タイムアウト | しばらく待って再実行してください |
| `Too big for GPU memory` | GPU メモリ不足 | 入力構造を削減するか、LightPFP の利用を検討してください |
| `Too many atoms` | 原子数が上限超過 | 原子数を減らしてください |
| `Too many atoms considering periodic boundary conditions` | PBC でのゴースト原子含む原子数が上限超過 | 原子数を減らしてください |
| `Too many neighbors` | 近傍数が上限超過（密度が高すぎる系） | 原子間距離を確認し、重なりを解消してください |

### matlantis-features エラー

| エラー | 原因 | 対処 |
|--------|------|------|
| `MatlantisError` | サーバー側でエラー発生 | エラーメッセージ詳細を確認してください |
| `RetriesExceeded` | 再実行数が `max_retries` を超えた | サーバー利用状況を確認するか `max_retries` を増やしてください |
| `PermutationIndexError` | フォノン symmetry での原子位置順列検出失敗 | 入力構造の対称性を確認してください |
| `PFPAPIError` | 汎用サーバーエラー | 一時的な接続不良の可能性があるため再試行してください |

## 関連ガイド

- **ASE 基礎操作** (ase-basics/SKILL.md): Atoms 生成、構造編集
- **構造モデリング** (modeling/SKILL.md): 表面、分子、置換
- **構造最適化** (optimization/SKILL.md): 構造最適化
- **分子動力学** (dynamics/SKILL.md): 分子動力学
