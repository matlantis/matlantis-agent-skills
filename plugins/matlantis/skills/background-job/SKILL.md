---
name: mt-background-job
description: >
  バックグラウンドジョブの実行・管理を扱うスキルです。
  mtl-bg-job, background job, バックグラウンドジョブ, 長時間計算,
  チェックポイント, checkpoint, キュー管理, queue,
  Notebook分離実行, ブラウザ閉じても計算継続
  に関するコード生成時に使用してください。
---
# バックグラウンドジョブ

## 概要

Matlantis のバックグラウンドジョブ機能は、長時間計算を Notebook セッションから切り離して実行するための仕組みです。ブラウザを閉じても計算が継続されるため、数時間から数日に及ぶ MD シミュレーション、大量の構造最適化、phonon / QHA 計算などに適しています。

バックグラウンドジョブでは、投入時に Notebook のコピーが作成され、コピー側が実行されます。元ファイルはジョブ実行中でも自由に編集できますが、Python コードで読み込む外部ファイル（データファイルなど）はコピーされない点に注意してください。

このスキルは、ユーザーが background 実行を明示した場合、または Notebook 実行時に foreground / background の選択肢を提示する必要がある場合に優先して参照してください。CLI から background 実行する標準手段は `mtl-bg-job` です。

### 主な制約

- デフォルト最大並列実行数: **3 ジョブ**
- ディスク上限: **100 GiB**（トラジェクトリファイルは肥大化しやすい）
- Notebook を停止すると実行中のジョブも**全て停止**される
- バックグラウンドジョブで作成した出力ファイルでは 3D NGL Viewer が表示されない

## ワークフロー

バックグラウンドジョブの利用は、以下のフローで進めます。

```
1. バックグラウンドジョブが必要か判断する
   ↓
2. Notebook にチェックポイント処理を追加する
   ↓
3. ジョブを投入する（GUI またはコマンドライン）
   ↓
4. ジョブを監視する（ジョブ一覧画面 / mtl-bg-job list）
   ↓
5. 完了後に結果を取得・確認する
```

### いつバックグラウンドジョブを使うか

以下の判断表を参考にしてください。

| バックグラウンドジョブを使う | 通常実行で良い |
|---|---|
| 数時間〜数日の MD シミュレーション | 数分で終わる計算 |
| 大量の構造最適化（ハイスループットスクリーニング） | インタラクティブに確認しながら進めたい計算 |
| phonon / QHA / 長い post-processing | 小さな分子の単発最適化 |
| ブラウザを閉じても継続させたい計算 | 短時間ジョブを大量投入したい場合（非効率） |

短時間計算を大量に投入するケースでは、バックグラウンドジョブではなく通常実行でループ処理する方が効率的です。キュー待ちが発生し、並列実行数の上限（3）がボトルネックになります。

### 実行方式の選び方

Notebook 実行を提案する際は、以下のルールで方式を選んでください。

1. ユーザーが方式を未指定なら、foreground と background の両方を明示する
2. ユーザーが background を指定したら、`mtl-bg-job run` を使う
3. `jupyter nbconvert --to notebook --execute ...` は foreground 実行としてのみ扱う
4. 長時間計算では background を推奨する

出力 Notebook 名が未指定なら、`input.ipynb` に対して `input_results.ipynb` や `input_bg.ipynb` のような派生名を提案して構いません。ただし、実行方式は `mtl-bg-job` から変更しないでください。

Notebook を自動生成・編集してから background 実行する場合は、Notebook metadata の `metadata.kernelspec.name` が **既存 kernel** を指していることを前提にしてください。`python3` のような汎用名を新規に入れてはいけません。

## 実装パターン

### パターン A: GUI からのジョブ投入

Notebook ツールバーの**再生ボタン**をクリックして投入します。

1. Notebook を開く
2. ツールバーの再生ボタンをクリック
3. 出力先ファイル名を設定する
4. 優先度を設定する（高優先度 / 低優先度）
5. 投入する

ツールバーの再生ボタン右隣のボタンから**ジョブ一覧画面**にアクセスできます。Active タブでは実行中・待機中のジョブを確認でき、ドラッグ&ドロップでキューの並び替えが可能です。

### パターン B: コマンドラインからのジョブ投入

`mtl-bg-job` CLI を使用します。ユーザーが「バックグラウンドジョブで実行して」と言った場合の第一候補はこの方法です。

```shell
# ジョブ投入
mtl-bg-job run phonon.ipynb phonon_results.ipynb

# 実行中のジョブ一覧
mtl-bg-job list --phase active

# JSON 形式で取得（スクリプト連携用）
mtl-bg-job list --phase active --format json | jq '.jobs | map(.status) | unique'

# 詳細ヘルプ
mtl-bg-job --help
```

kernel 周りでは次のルールを守ってください。

1. Notebook 側に既存 kernel (`python313`, `python311` など) が設定されている場合は、その kernel を前提に `mtl-bg-job run` してください。
2. `mtl-bg-job run` が `kernel '...' is not available` で失敗した場合は、エラーメッセージに出た利用可能 kernel のいずれか既存のものを使って即座に再実行してください。
3. 再実行時は `mtl-bg-job run --kernel <existing-kernel> input.ipynb output.ipynb` を使ってください。
4. Notebook をこれから書き直す場合は、再発防止のため Notebook metadata の `metadata.kernelspec.name` も同じ既存 kernel に合わせてください。

```shell
# Notebook metadata の kernel が使えない場合は既存 kernel を明示する
mtl-bg-job run --kernel python313 phonon.ipynb phonon_results.ipynb
```

foreground 実行が必要な場合のみ、別手段として `jupyter nbconvert --to notebook --execute ...` を使ってください。background 実行の代替として `nbconvert` を選んではいけません。

### パターン C: MD のチェックポイント設計

長時間 MD では、再起動時に最初からやり直しにならないよう `Trajectory` の append モード (`"a"`) を使います。

```python
from ase.io import read, Trajectory
from ase.md.langevin import Langevin
from ase import units
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

atoms = read('initial.cif')
atoms.calc = ASECalculator(
    Estimator(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
)

# "a" = append モード: 再起動後も続きから書き込まれる
traj = Trajectory('md.traj', 'a', atoms)
dyn = Langevin(atoms, timestep=1.0*units.fs, temperature_K=300, friction=0.01)
dyn.attach(traj.write, interval=100)  # 100 ステップごとに保存

dyn.run(100000)
```

**ポイント**:
- `Trajectory` の第 2 引数に `"a"` を指定すると、ファイルが既に存在する場合は末尾に追記されます
- `interval=100` のように保存間隔を調整し、ディスク消費を抑えてください
- `.traj` ファイルは大きくなりやすいため、ディスク上限 100 GiB に注意が必要です

### パターン D: ループ計算のチェックポイント設計

大量の構造最適化など、ループで回す計算では JSON ファイルに進捗を逐次保存します。

```python
import json, os
from ase.optimize import BFGS
from ase.io import read
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

structures = read('structures.xyz', ':')
results_file = 'optimization_results.json'

# チェックポイントからの復帰
if os.path.exists(results_file):
    with open(results_file) as f:
        results = json.load(f)
    start_idx = len(results)
else:
    results = []
    start_idx = 0

for i, atoms in enumerate(structures[start_idx:], start=start_idx):
    atoms.calc = ASECalculator(
        Estimator(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
    )
    opt = BFGS(atoms, trajectory=f'opt_{i}.traj')
    opt.run(fmax=0.05)

    results.append({
        'index': i,
        'energy': atoms.get_potential_energy(),
        'converged': opt.converged()
    })

    # 毎回保存することで、中断時も途中まで の結果が残る
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Completed {i+1}/{len(structures)}")
```

**ポイント**:
- `start_idx = len(results)` で、既に完了した計算をスキップします
- `json.dump` を毎イテレーション実行することで、中断しても途中結果が失われません
- リストに結果を溜め続けるとメモリが増大するため、大規模な場合はファイルへの逐次書き出しを検討してください

### パターン E: ログ管理

バックグラウンドジョブでは、進捗をログとして出力することが重要です。

```python
from datetime import datetime

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

log("Starting calculation...")
for i in range(total_steps):
    if i % 1000 == 0:
        log(f"Step {i}/{total_steps} completed")
log("Calculation finished")
```

**ポイント**:
- `flush=True` で出力を即座に書き込みます（バッファリングによる遅延を防止）
- 全ステップでログを出力するとファイルが肥大化するため、一定間隔で出力してください

### パターン F: メモリ最適化

長時間計算ではメモリ消費に注意が必要です。

```python
# NG: リストに結果を溜め続ける（メモリが増大する）
energies = []
for i in range(100000):
    energies.append(atoms.get_potential_energy())

# OK: ファイルに逐次書き出す
with open('energies.txt', 'w') as f:
    for i in range(100000):
        energy = atoms.get_potential_energy()
        f.write(f"{i} {energy}\n")
        f.flush()
```

## ベストプラクティス

### pfcc-extras を優先する

バックグラウンドジョブのコード実装においても、ASE と pfcc-extras の両方で実現できる機能は pfcc-extras を優先して使用してください。特に以下の点に注意してください。
- **可視化**: `ase.visualize.view` ではなく `pfcc_extras.show_gui` / `view_ngl` を使用
- **PBC ラッピング**: `atoms.wrap()` ではなく `pfcc_extras.structure.wrap_molecule` を使用（分子断裂を防止）
- **衝突検出**: `pfcc_extras.structure.collision.CollisionDetector` を使用
- **ジョブ実行制御**: 複数 Notebook の並列・逐次実行には `pfcc_extras.job_scheduler.runner.run_jobs` を使用
- **軌跡変換**: MDTraj / MDAnalysis への変換には `pfcc_extras.trajectory` のユーティリティを使用

### 優先度の設定

| 優先度 | 用途 |
|---|---|
| 高優先度（デフォルト） | 急ぎの計算 |
| 低優先度 | 時間的余裕のある計算。他ユーザーへの影響を抑え、チーム全体のトークン消費を平準化 |

優先度は投入ダイアログ、またはジョブ一覧画面で待機中ジョブに対して変更可能です。

### 並列数の管理

- デフォルトの最大並列実行数は **3** です
- 長時間計算を 3 つまで同時に実行できます（数時間以上のジョブ向き）
- 短時間計算を大量に投入するのはキュー待ちが発生するため非効率です

### チーム運用のルール

1. **大規模計算は事前共有**: 長時間かかる計算はチーム内で事前に共有してください
2. **低優先度の活用**: 締切に余裕がある計算は低優先度で投入し、リソースを譲り合ってください
3. **完了後は Group Drive に保存**: 計算結果を Group Drive に置いてチームで再利用可能にしてください

### メール通知

Matlantis 設定ダイアログで通知先メールアドレスを設定すると、ジョブ完了時にメール通知を受け取れます。長時間ジョブの完了を見逃さないために設定を推奨します。

### チェックポイント設計の原則

| 計算タイプ | チェックポイント方式 | キーポイント |
|---|---|---|
| MD シミュレーション | `Trajectory("file.traj", "a", atoms)` | append モードで再起動後も続きから |
| ループ計算（最適化バッチ等） | JSON ファイルに毎回保存 | `start_idx = len(results)` で復帰位置を計算 |
| 大規模データ出力 | ファイルへの逐次書き出し | メモリに溜めず `f.write()` + `f.flush()` |

## よくあるエラーと対処

| エラー / 問題 | 原因 | 対処 |
|---|---|---|
| ジョブが突然停止した | Notebook が停止された | Notebook を停止するとジョブも全て停止します。Notebook を再起動してジョブを再投入してください |
| ディスク容量不足 | トラジェクトリファイルの肥大化 | `interval` を大きくして保存頻度を下げる。不要な `.traj` ファイルを削除する。上限は 100 GiB |
| 再起動後に最初からやり直しになる | チェックポイントが未設定 | MD は append モード (`"a"`)、ループ計算は JSON チェックポイントを実装してください |
| NGL Viewer が表示されない | バックグラウンドジョブの制約 | バックグラウンドジョブの出力ファイルでは 3D Viewer が利用できません。結果確認は通常の Notebook で行ってください |
| ジョブがキュー待ちのまま進まない | 並列数上限（3）に達している | ジョブ一覧画面で優先度を調整するか、不要なジョブをキャンセルしてください |
| メモリ不足（OOM） | リストへの結果蓄積 | ファイルへの逐次書き出しに切り替えてください |
| ログが出力されない | バッファリング | `print(..., flush=True)` を使用してください |
| 外部ファイルが見つからない | ジョブ投入時にコピーされない | Notebook がコピーされますが、参照する外部ファイルはコピーされません。絶対パスまたは共通の場所にファイルを配置してください |

## 関連ガイド

- **統合ワークフロー** (`integrated-workflows/SKILL.md`): 大量の構造最適化をバックグラウンドジョブで実行するパターン
- **pfcc_extras ユーティリティ** (`pfcc-extras/SKILL.md`): `run_jobs` によるNotebookジョブの自動実行
- **SSH 接続** (`ssh/SKILL.md`): SSH 経由でのジョブ監視・ファイル転送
