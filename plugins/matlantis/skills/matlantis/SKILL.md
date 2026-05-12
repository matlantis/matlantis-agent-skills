---
name: matlantis
description: >
  Guide for writing Python code using Matlantis, PFN's cloud-based atomistic simulation platform
  powered by neural network potentials (PFP). Use this skill whenever the user asks about
  Matlantis, PFP (Preferred Potential), matlantis-features, pfp-api-client, MTCSP,
  or any atomistic/molecular simulation on the Matlantis platform. Also trigger when code involves
  imports like pfp_api_client, matlantis_features, or mtcsp_api_client.
---

## What is Matlantis

Matlantis は Preferred Networks (PFN) が提供するクラウドベースの原子シミュレーションプラットフォームです。
コアとなるのは PFP (Preferred Potential) で、グラフニューラルネットワークベースの汎用原子間ポテンシャルです。
DFT に近い精度を低コストで実現し、大規模（10,000原子以上）・長時間（数百万ステップ）の分子動力学が可能です。

Matlantis はクラウド上の JupyterLab で動作し、全ての計算は PFP の API サーバーを経由します。

## Environment Preflight

Python を使うコード生成・実行・デバッグに入る前に、必ず Python が使えるか確認してください。

```bash
python --version > /dev/null 2>&1
```

- `python --version > /dev/null 2>&1` が成功する場合: そのまま Python 関連の作業を進めてください。
- `python --version > /dev/null 2>&1` の終了コードが 1 の場合: Python はまだ使えない前提で処理を止め、`Please exit Claude, run use_venv, and relaunch.` と案内して `use_venv` 実行を推奨してください。

上記の条件を満たさない状態では、Python コードの実行や Python 前提の追加調査を続けてはいけません。


## Library Architecture

```
pfp-api-client          <- Core: PFP server API (energy, forces, stress)
  +-- ASECalculator     <- ASE integration (the standard interface)

matlantis-features      <- High-level: property computation (optimization, MD, elasticity, phonon, etc.)

mtcsp-api-client        <- Crystal structure prediction (convex hull, derivative structures)

matlantis-group-drive-client <- Shared file storage
```

全てが ASE (Atomic Simulation Environment) の上に構築されています。ASE を知っていれば Matlantis を使えます。


## PFP: 計算モードの選択

計算モードは最も重要な設定です。以下のフローチャートに従って選択してください。

```
計算対象は？
+-- 有機分子・単分子系
|   +-- 9元素のみ (H,C,N,O,F,P,S,Cl,Br) -> WB97XD
|   +-- それ以外 -> R2SCAN_PLUS_D3
+-- 遷移金属酸化物・強相関系 -> PBE_U (分散力が重要なら PBE_U_PLUS_D3)
+-- 分子結晶・液体・吸着系
|   +-- 実験値との定量比較が重要 -> R2SCAN_PLUS_D3
|   +-- 定性的・スクリーニング用途 -> PBE_PLUS_D3
+-- バルク結晶の形成エネルギー
|   +-- 実験値と比較 -> R2SCAN
|   +-- Materials Project と比較 -> PBE (または PBE_U)
+-- 迷ったとき -> PBE
```

### モード一覧

| モード | 汎関数 | 対応元素 | 主な用途 |
|--------|--------|----------|----------|
| `PBE` | GGA-PBE | 72元素 | 汎用（デフォルト） |
| `PBE_U` | GGA-PBE + Hubbard U | 72元素 | 遷移金属酸化物・強相関系 |
| `PBE_PLUS_D3` | GGA-PBE + D3分散補正 | 72元素 | 分子結晶・液体（定性） |
| `PBE_U_PLUS_D3` | GGA-PBE + U + D3 | 72元素 | 酸化物+有機配位子 |
| `R2SCAN` | meta-GGA R2SCAN | 72元素 | バルク形成エネルギー |
| `R2SCAN_PLUS_D3` | meta-GGA R2SCAN + D3 | 72元素 | 分子・吸着（定量） |
| `WB97XD` | hybrid omegaB97X-D | 9元素のみ | 有機分子 |

### 対応元素

- PBE/PBE_U/R2SCAN: 最大96元素 (H〜Pu)、バージョン依存
- WB97XD: 9元素のみ (H, C, N, O, F, P, S, Cl, Br)
- 実行時確認: `estimator.supported_elements()`

### PFP の適用範囲外

PFP は以下を扱えません:
- 励起状態
- 帯電系
- 電子物性 (HOMO/LUMO, バンドギャップ)
- 粗視化モデル
- 磁気秩序 (Hubbard U 以上の精度)


## Calculator の基本セットアップ

```python
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

# 基本セットアップ
estimator = Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.PBE)
calculator = ASECalculator(estimator)

# Atoms にアタッチして計算
from ase.build import bulk
atoms = bulk("Si")
atoms.calc = calculator
energy = atoms.get_potential_energy()  # eV
forces = atoms.get_forces()            # eV/A
stress = atoms.get_stress()            # eV/A^3
```

### matlantis-features 使用時の estimator_fn パターン

matlantis-features の Feature は `estimator_fn`（ファクトリ関数）を受け取ります。
Estimator インスタンスを直接渡してはいけません。

```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode

estimator_fn = pfp_estimator_fn(
    model_version="v8.0.0",
    calc_mode=EstimatorCalcMode.PBE,
)
```


## Critical Rules & Common Pitfalls

### 必守ルール

1. **Estimator の共有禁止**: Estimator インスタンスを複数の Calculator や Feature で共有しないでください。必ず個別に生成します。`MultiCalculatorUseDetected` エラーの原因になります。
2. **estimator_fn の使用**: matlantis-features 使用時は `estimator_fn`（ファクトリ関数）を渡してください。Estimator インスタンスを直接渡すと予期しない動作になります。
3. **calculator.reset()**: `Atoms` オブジェクトをインプレースで変更した後は `calculator.reset()` を呼んでください。Calculator は結果をキャッシュしています。
4. **PBC の設定**: 結晶・スラブ系では `atoms.pbc = True` を設定してください。PBC なしでは周期系の計算結果が不正になります。

### Notebook 実行ポリシー

Notebook 実行を伴う依頼では、**foreground 実行**と**background 実行**の両方が選べることを明示してください。

| 実行方式 | 使いどころ | 代表手段 |
|----------|------------|----------|
| foreground 実行 | 短時間の確認、対話的に様子を見ながら進めたい場合 | `jupyter nbconvert --to notebook --execute ...` |
| background 実行 | 長時間計算、離席、ブラウザを閉じても継続したい場合 | `mtl-bg-job run input.ipynb output.ipynb` |

以下を必ず守ってください。

1. **選択肢の明示**: Notebook を実行する提案では、foreground と background の両方の選択肢があることを明示してください。
2. **background 指定時は `mtl-bg-job` を優先**: ユーザーが「バックグラウンドジョブで」「background で」「Run in Background」などと指定した場合、`mtl-bg-job run` を使ってください。
3. **`nbconvert` は foreground 専用**: `jupyter nbconvert --to notebook --execute --inplace ...` は foreground 実行の選択肢としてのみ扱ってください。background 指定時には使わないでください。
4. **長時間計算では background を推奨**: MD、phonon、QHA、大量スクリーニング、長い post-processing では background 実行を推奨してください。
5. **方式未指定でも background を隠さない**: ユーザーが単に「実行してください」と言った場合も、foreground だけに決め打ちせず、background 実行の選択肢があることを伝えてください。
6. **backgroundで実行する場合にipynb形式にわざわざ変換しない**: background実行する場合、既に `.py` 形式でコードが生成されている場合は、わざわざ `.ipynb` に変換しないでください。
7. **Notebook の kernel は既存名を使う**: `.ipynb` を新規作成・更新する場合、`metadata.kernelspec.name` には `python313` や `python311` のような既存 kernel 名を設定し、`python3` のような汎用名を新規に書き込まないでください。
8. **`mtl-bg-job` の kernel エラーは即フォールバック**: `mtl-bg-job run` が `kernel '...' is not available` で失敗した場合は、利用可能な既存 kernel を選び、`mtl-bg-job run --kernel <existing-kernel> ...` で再実行してください。

### よくある落とし穴

| 落とし穴 | 症状 | 対処 |
|----------|------|------|
| 計算モードの不一致 | 結果が不正（エラーにならない） | フローチャートに従い系に合ったモードを選択 |
| 異なるモード間のエネルギー比較 | エネルギー差が無意味 | 同一モードで統一 |
| WB97XD を凝縮相に使用 | 精度が大幅低下 | 液体・固体には R2SCAN_PLUS_D3 等を使用 |
| D3 補正の未適用 | 分子間相互作用の過小評価 | 分子結晶・液体・吸着・層状物質には `*_PLUS_D3` を使用 |
| セル最適化の欠落 | 格子定数が未最適化 | 結晶には `filter=True` または `UnitCellASEFilter` を使用 |
| fmax が不適切 | 収束不足または計算過多 | 本番: 0.01 eV/A、簡易: 0.05 eV/A |
| フォノンのスーパーセル不足 | 結果が不正確 | `supercell=4-6` から開始、必要に応じ増加 |
| 入力サイズ超過 | `Too many atoms` エラー | `atoms.calc.results['calc_stats']` で確認 |
| PBE_U を非酸化物に適用 | 不適切な補正 | 酸化物以外の遷移金属は PBE を使用 |

### PFP エラーリファレンス

| 例外名 | 原因 | 対処 |
|--------|------|------|
| `RetriesExceeded` | サーバー再試行回数超過 | `max_retries` 増加またはサーバー状況確認 |
| `ConcurrentUseDetected` | ASECalculator への並列アクセス | スレッドごとに個別インスタンス生成 |
| `MultiCalculatorUseDetected` | Estimator の複数 Calculator 共有 | Calculator ごとに個別 Estimator 生成 |
| `Too many atoms` | 原子数上限超過 | 原子数を削減 |
| `Too many neighbors` | 近傍数上限超過（高密度系） | 原子数削減または密度低下 |
| `Cell is too small` | PBC に対しセルが小さすぎる | セルサイズ確認 |
| `Too big for GPU memory` | GPU メモリ不足 | 待機後再実行、入力構造削減 |


## 関連スキル

以下の個別スキルがタスクに応じて動的に読み込まれます。
各スキルは独立した SKILL.md を持ち、トリガー条件に基づいて自動的にコンテキストに追加されます。

| スキル | 内容 | ユースケース |
|--------|------|-------------|
| `setup/` | Calculator初期化、構造I/O | 全ての計算の出発点 |
| `ase-basics/` | ASE Atoms基本操作 | 構造の作成・編集の基礎 |
| `modeling/` | 構造モデリング | 表面/分子/液体/HEA等の構造構築 |
| `optimization/` | 構造最適化 | 安定構造の探索 |
| `dynamics/` | 分子動力学 | 有限温度での時間発展 |
| `reaction/` | 反応経路探索 | NEB/ReactionString/遷移状態 |
| `property-analysis/` | 物性計算 | 弾性/フォノン/拡散/粘度/熱伝導 |
| `visualization/` | 可視化 | 構造・軌跡の表示 |
| `external-data/` | 外部データ連携 | PubChem/Materials Project/VASP/Gaussian |
| `storage/` | Group Driveストレージ | データの永続保存・チーム共有 |
| `background-job/` | バックグラウンドジョブ | 長時間計算の実行・管理 |
| `pfcc-extras/` | pfcc_extrasユーティリティ | ジョブスケジューリング・堆積・軌跡変換 |
| `integrated-workflows/` | 統合ワークフロー | スクリーニング/吸着/欠陥パイプライン |
| `ssh/` | SSH接続 | ローカルPCからMatlantisへのSSH接続 |


## Reference Navigation Guide

全ての詳細ドキュメントは `references/` ディレクトリにあります。

### Getting Started
- `references/guidebook/` - プラットフォーム利用ガイド
- `references/tutorials/` - ASE基礎からMDまでの7章チュートリアル

### Model Theory
- `references/theory/pfp-description/` - PFPアーキテクチャ、訓練データ、対応元素
- `references/theory/d3-description/` - D3分散補正の詳細
- `references/theory/pfp-validation/` - PFP精度ベンチマーク
- `references/theory/d3-validation/` - D3パラメータ最適化結果
- `references/theory/mtcsp/` - MTCSP結晶構造予測の概要

### API Reference
- `references/api/pfp-api-client/` - pfp_api_client モジュールドキュメント
- `references/api/matlantis-features/` - matlantis_features モジュールドキュメント（最大セクション）
- `references/api/mtcsp-api-reference/` - MTCSP API ドキュメント
- `references/api/matlantis-group-drive-client/` - Group Drive クライアントドキュメント

### Examples
- `references/examples/` - 公式例（全主要機能のランナブルノートブック）
- `references/community-examples/` - コミュニティ貢献例（触媒/電池/ポリマー/潤滑剤等）

### API ドキュメントの検索方法

`references/api/matlantis-features/` 内の `generated/` ディレクトリにクラス別APIドキュメントがあります。
クラス名でファイルを検索してください（例: `matlantis_features.features.common.opt.BFGSASEOptFeature.md`）。
`tips/` ディレクトリには高度な使用パターンがあります。


## サポート・お問い合わせ

スキルファイルに関する、ご要望・ご質問がある場合は、Matlantis サポート窓口へお問い合わせください: https://matlantis.zendesk.com/hc/ja
