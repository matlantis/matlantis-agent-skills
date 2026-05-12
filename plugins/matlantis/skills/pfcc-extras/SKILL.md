---
name: mt-pfcc-extras
description: >
  pfcc_extras ユーティリティツール群を扱うスキルです。
  pfcc_extras, show_gui, view_ngl, run_jobs,
  ResourceAwareJobScheduler, QueueScheduler, PapermillScheduler,
  DepositionScheduler, DeleteMoleculeScheduler, TemperatureScaleScheduler,
  ApplyUniformEfield, ase_traj_converter,
  PartialOccupancy, wrap_molecule, get_mol_list,
  ジョブスケジューリング, 堆積シミュレーション, 軌跡変換
  に関するコード生成時に使用してください。
---
# pfcc_extras ユーティリティ

## 概要

`pfcc_extras` は、Matlantis / ASE ワークフローで不足しやすい「可視化」「ジョブ実行制御」「MD 補助」「構造加工」を補完する実務ユーティリティ群です。個別の ASE 操作を組み合わせる際に生じるボイラープレートコードを削減し、再現可能な計算ワークフローの構築を支援します。

### 機能レイヤー

| レイヤー | 役割 | 主要 API |
|---|---|---|
| Viewer Layer | 構造・軌跡の即時確認 | `show_gui`, `view_ngl`, `AddEditor` |
| Structure Layer | 入力構造の整備・加工 | `PartialOccupancy`, `wrap_molecule`, `get_mol_list` |
| Dynamics Layer | MD 実験条件の制御 | `DepositionScheduler`, `ApplyUniformEfield` |
| Execution Layer | Notebook ジョブの投入・監視 | `run_jobs`, `ResourceAwareJobScheduler`, `QueueScheduler` |

## ワークフロー

```
1. ユーティリティカテゴリを選択する（可視化 / 構造 / MD / 実行制御）
   ↓
2. 対象 API を設定・パラメータを調整する
   ↓
3. 処理を実行する
   ↓
4. 結果を検証する（構造の整合性、ジョブステータスなど）
```

## 実装パターン

### パターン A: 可視化（Viewer Layer）

解析前に構造の破綻・衝突・ラップ不整合を短時間で検出するためのエントリーポイントです。

#### show_gui による可視化

```python
from pfcc_extras import show_gui
from ase.io import read

atoms = read("input.traj")
show_gui(atoms, panel=True, axes=True, atom_index=True)
```

**主要パラメータ**:
- `panel=True`: 操作パネルを表示
- `ball_size=0.3`: 原子の表示サイズ
- `axes=True`: 座標軸を表示
- `atom_index=True`: 原子インデックスを表示

#### view_ngl による NGL ビューア表示

```python
from pfcc_extras import view_ngl

# 単一構造の表示
view_ngl(atoms)

# 軌跡の表示（view_ngl_traj）
from pfcc_extras import view_ngl_traj
view_ngl_traj(trajectory)
```

**バリエーション**:
- `view_ngl(atoms_or_traj)`: 汎用表示
- `view_ngl_atoms(atoms)`: 単一構造の表示
- `view_ngl_traj(traj)`: 軌跡のアニメーション表示

#### AddEditor による構造編集

```python
from pfcc_extras import AddEditor

editor = AddEditor(atoms, fallback_calc_mode="pbe", fallback_model_version="v8.0.0")
```

GUI 上で原子の追加・削除・移動などをインタラクティブに行えます。

### パターン B: 構造加工（Structure Layer）

#### 部分占有構造の離散化（PartialOccupancy）

部分占有構造（結晶学データで各サイトに複数元素が確率的に配置される構造）を、計算可能な離散配置へ変換します。

```python
from pfcc_extras.structure import PartialOccupancy

po = PartialOccupancy(
    input_structure=atoms,
    csv_path="occupancy.csv",
    occupancy_header="occupancy",
    random_seed=42,
)
po.set_parameters(cell_repeat=(2, 2, 2))
realized_atoms = po.assign_atom_positions()
```

**実務要点**:
- `cell_repeat` を大きくすると占有率の整数比近似が改善されます
- 占有率が小数で割り切れない場合は `cell_repeat` を増やして離散化誤差を低減してください
- `charge_values` を与えると、丸め後の電荷整合性を検証しやすくなります

#### 分子ラップと分子リスト

PBC（周期境界条件）処理を分子単位で適用し、分子が断裂するのを防ぎます。

```python
from pfcc_extras.structure import wrap_molecule, get_mol_list

# 分子単位のリストを取得
mol_list = get_mol_list(atoms)

# 分子単位で PBC ラップ
wrapped_atoms = wrap_molecule(atoms)
```

**ポイント**:
- PBC 処理は原子単位ではなく、`get_mol_list` / `wrap_molecule` で分子単位に適用してください
- 原子単位で wrap すると分子が分断される場合があります

#### 近傍・結合解析

```python
from pfcc_extras.structure import get_neighbors, get_connectivity_matrix

# 近傍原子の取得
neighbors = get_neighbors(atoms, cutoff=3.0)

# 結合行列の取得
connectivity = get_connectivity_matrix(atoms, cutoff=1.8)
```

**ポイント**:
- `cutoff` の閾値はワークフロー全体で統一し、結合判定の揺らぎを防いでください
- 系に応じて適切な `cutoff` 値を設定してください（過結合・未結合の両方を防ぐ）

#### 表面・欠陥処理

表面スラブの構築や欠陥構造の生成に関するユーティリティも提供されています。吸着探索では初期力閾値（`TH_max_f`, `TH_min_f`）と結合長変化閾値（`tol`）を使って不適切な初期配置を速やかに除外できます。

### パターン C: MD イベントスケジューリング（Dynamics Layer）

#### DepositionScheduler（分子堆積）

時間軸に沿って分子を堆積させる非平衡 MD を構築します。

```python
from pfcc_extras.molecular_dynamics.scheduler import (
    DepositionScheduler,
    convert_kinetic_energy_to_velocity,
)

# エネルギーから速度への変換（比較可能な条件設定のため）
v = convert_kinetic_energy_to_velocity(1.0, 18.0)  # 1 eV, 18 amu

scheduler = DepositionScheduler(
    molecules=[molecule_atoms],
    fractions=[1.0],
    incident_velocities=[[0.0, 0.0, -v]],
    initial_height=5.0,
    axis=2,
    seed=42,
)
```

**主要パラメータ**:
- `molecules`: 堆積する分子の ASE Atoms リスト
- `fractions`: 各分子の堆積比率（正規化される）
- `incident_velocities` または `incident_energy`: 入射速度またはエネルギー
- `initial_height`: 初期高さ
- `axis`: 堆積方向の軸（0=x, 1=y, 2=z）
- `seed`: 再現性のための乱数シード
- `temperature_K`: 初期温度

**ガードレール**:
- `fractions` の総和は自動的に正規化されます
- 分子サイズに応じて `distance_tol` を設定してください
- `calculate_default_distance_tol(atoms, ratio=...)` で適切な距離閾値を算出できます

**ベストプラクティス**:
- 速度の直接指定ではなく「入射エネルギー → 速度変換」を共通化すると、ケース間の比較が容易になります

#### ApplyUniformEfield（外部電場）

MD シミュレーション中に一様電場を印加します。

```python
from pfcc_extras.molecular_dynamics.constraints import ApplyUniformEfield

efield_constraint = ApplyUniformEfield(
    efield=(0.0, 0.0, 0.01),  # (Ex, Ey, Ez) in V/A
    charges=[...]              # 各原子の電荷リスト
)
atoms.set_constraint(efield_constraint)
```

**ガードレール**:
- `charges` の長さは原子数と一致させてください
- 時間依存電場を適用する場合は、MD ステップごとに `efield` を更新してください
- `ApplyUniformEfield` は `FixConstraint` ベースの実装です

### パターン D: ジョブ実行制御（Execution Layer）

#### run_jobs による並列/逐次実行

複数の Notebook を自動実行します。

```python
from pfcc_extras.job_scheduler.runner import run_jobs
from pfcc_extras.job_scheduler.core import Job

jobs = [Job("a.ipynb"), Job("b.ipynb")]
results = run_jobs(
    jobs,
    scheduler="parallel",
    max_workers=2,
    monitor_interval=10,
    wait_on_finish=True,
)
```

**`scheduler` パラメータの選択肢**:
- `"parallel"`: 並列実行（`max_workers` で同時実行数を制御）
- `"sequential"`: 逐次実行
- `"papermill"`: Papermill ベースの実行

**ジョブステータス**:

| ステータス | 意味 |
|---|---|
| `READY` | 実行準備完了 |
| `WAITING` | 実行待ち |
| `FAILED` | 実行失敗 |
| `SKIPPED` | スキップ（出力 Notebook が既に存在する場合など） |

**ポイント**:
- `JobResult.status` を終端状態まで確認し、`FAILED` と `SKIPPED` を分離処理してください
- 出力 Notebook が既に存在する場合に `SKIPPED` 扱いになる挙動を前提に、再実行可能なバッチを設計してください

#### ParameterizedJob

パラメータを変えて同じ Notebook を複数回実行する場合に使用します。

```python
from pfcc_extras.job_scheduler.core import ParameterizedJob

jobs = [
    ParameterizedJob("template.ipynb", parameters={"temperature": 300}),
    ParameterizedJob("template.ipynb", parameters={"temperature": 500}),
]
results = run_jobs(jobs, scheduler="parallel", max_workers=2)
```

#### ResourceAwareJobScheduler

CPU / メモリの制約を見ながら自動的にジョブ投入を制御します。

```python
from pfcc_extras.job_scheduler.scheduler import ResourceAwareJobScheduler

scheduler = ResourceAwareJobScheduler(
    resource_limits={"cpu": 80, "memory": 70},  # 使用率の上限 (%)
    monitor_interval=30,
)
```

**ポイント**:
- 並列数を固定するのではなく、`ResourceAwareJobScheduler` を使って CPU / メモリ閾値に基づく投入制御を行ってください
- リソースのボトルネックを検出し、過負荷を防止します

#### QueueScheduler

同時実行数を制限したキューベースのスケジューラです。

```python
from pfcc_extras.job_scheduler.scheduler import QueueScheduler

scheduler = QueueScheduler(
    max_concurrent=3,
    monitor_interval=10,
)
```

#### PapermillScheduler

Papermill をバックエンドとして Notebook を実行します。

```python
from pfcc_extras.job_scheduler.scheduler import PapermillScheduler
```

### パターン E: 軌跡の相互変換（Trajectory Interoperability）

ASE の軌跡を MD 解析ライブラリ（MDTraj / MDAnalysis）に変換します。

```python
from pfcc_extras.trajectory import asetraj_to_mdtraj, asetraj_to_mdanalysis

# MDTraj 形式に変換
md_traj = asetraj_to_mdtraj(traj, topology=None, set_PBC=True, bond_cutoff=1.2)

# MDAnalysis 形式に変換
universe = asetraj_to_mdanalysis(traj, topology=None, set_PBC=True, bond_cutoff=1.2)
```

**ガードレール**:
- 全フレームで原子数・元素順を一致させてください
- `bond_cutoff` は過結合 / 未結合の両方を防ぐよう、系ごとに調整してください
- 変換前に PBC・セル・原子順を固定し、フォーマット変換で意味が変わらないことを確認してください

## ベストプラクティス

### ジョブ実行の冪等化

出力 Notebook が既に存在する場合に `SKIPPED` 扱いになる挙動を前提に、再実行可能なバッチを設計してください。冪等なジョブ設計により、途中で失敗しても安全に再実行できます。

### 資源待ちの明示化

並列数を固定するのではなく、`ResourceAwareJobScheduler` を使い、CPU / メモリ閾値を超える投入を抑制してください。これにより、計算ノードの過負荷を防止しつつ、利用可能なリソースを最大限活用できます。

### 堆積条件の比較可能化

速度を直接指定するのではなく、`convert_kinetic_energy_to_velocity` による「入射エネルギー → 速度変換」を共通化してください。これにより、異なるケース間での条件比較が容易になります。

### 分子断裂の予防

PBC 処理は原子単位ではなく、`get_mol_list` / `wrap_molecule` で分子単位に適用してください。原子単位のラップでは、分子がセル境界で分断される場合があります。

### 近傍探索の一貫性

`get_neighbors` / `get_connectivity_matrix` の閾値をワークフロー全体で統一し、結合判定の揺らぎを防いでください。異なる処理段階で異なる `cutoff` を使うと、結合の有無が矛盾する結果になります。

## よくあるエラーと対処

| エラー / 問題 | 原因 | 対処 |
|---|---|---|
| `run_jobs` でジョブが `FAILED` になる | Notebook 内のエラー | `JobResult` の詳細を確認し、Notebook 単体で実行してデバッグしてください |
| `SKIPPED` と `FAILED` の混同 | 出力ファイルが既に存在 | `SKIPPED` は正常動作（冪等性）です。`FAILED` とは区別して処理してください |
| `DepositionScheduler` で分子が重なる | `distance_tol` が小さすぎる | `calculate_default_distance_tol(atoms, ratio=...)` で適切な距離閾値を算出してください |
| `ApplyUniformEfield` でエラー | `charges` の長さが原子数と不一致 | `len(charges) == len(atoms)` を確認してください |
| 部分占有の離散化誤差が大きい | `cell_repeat` が小さい | `cell_repeat` を増やして離散化誤差を低減してください |
| 軌跡変換でフレーム数が合わない | 原子数・元素順の不一致 | 全フレームで原子数と元素順が一致しているか確認してください |
| `bond_cutoff` で過結合が発生 | 閾値が大きすぎる | 系に応じて `bond_cutoff` を調整してください |
| 分子が PBC 境界で分断される | 原子単位の wrap | `wrap_molecule` で分子単位に wrap してください |
| `NotImplementedError` が発生 | 未実装モード | 一部の欠陥探索モード等は未実装です。運用分岐を設けてください |

## 関連ガイド

- **統合ワークフロー** (`integrated-workflows/SKILL.md`): pfcc_extras を使ったエンドツーエンドのパイプライン構築
- **バックグラウンドジョブ** (`background-job/SKILL.md`): 長時間ジョブの実行管理
