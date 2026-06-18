---
name: mt-integrated-workflows
description: >
  統合ワークフロー（パイプライン）テンプレートを扱うスキルです。
  構造スクリーニング, structure screening, 吸着構造探索, adsorption structure search,
  欠陥解析, defect analysis pipeline, バッチ処理, batch processing,
  パイプライン, pipeline, ワークフロー, workflow,
  エネルギーランキング, energy ranking, 早期棄却フィルタ, early rejection filter,
  2段階緩和, two-stage relaxation
  に関するコード生成時に使用してください。
---
# 統合ワークフロー

## 概要

統合ワークフローは、構造セットアップ、最適化、MD、反応解析、結果の保存などの個別処理を**実務でそのまま使える一連のパイプライン**として接続するためのプレイブックです。本ガイドでは、構造スクリーニング、吸着構造探索、欠陥解析の 3 つの代表的なワークフローを提供します。

### 設計原則

すべてのワークフローに共通する設計原則は以下のとおりです。

| 原則 | 説明 |
|---|---|
| 入力検証 | `Path.exists()` とファイル拡張子の事前検証を必須化する |
| 早期棄却 | 強すぎる力・原子の重なり構造は最適化前に除外し、計算資源の無駄を防ぐ |
| 例外処理 | 例外で全体停止せず、失敗した計算をサマリに残して次の候補へ継続する |
| 再実行性 | 出力は `runs/<workflow_name>/` 配下に分け、タイムスタンプ付きで保存する |
| サマリ再生成 | サマリ CSV は毎回再生成する（中間ファイル再利用時でも最後に再ソート） |
| 比較の一貫性 | 同一の `calc_mode` / `model_version` / 収束条件で横比較する |

## ワークフロー

統合ワークフローの全体フローは以下のとおりです。

```
1. ワークフローを定義する（対象構造、パラメータ、出力先）
   ↓
2. 入力を検証する（ファイル存在確認、拡張子チェック）
   ↓
3. パイプラインを実行する（生成 → 棄却 → 緩和 → 保存）
   ↓
4. エラーを処理する（失敗をログに記録、次候補へ継続）
   ↓
5. サマリを生成する（ランキング CSV、エラーログ）
```

## 実装パターン

### パターン A: 構造スクリーニングパイプライン

**目的**: 複数の構造ファイル（`.cif` / `.traj`）を一括で最適化し、エネルギーでランキングします。

**入力**: 構造ファイル一覧（`*.cif` / `*.traj`）
**出力**: 最適化済み構造群 + エネルギーランキング CSV

```python
from pathlib import Path
import csv
from typing import List, Tuple

from ase.io import read, write
from ase.optimize import FIRE
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator


def build_calculator(calc_mode: str = "pbe", model_version: str = "latest") -> ASECalculator:
    """統一された calculator を生成する。全構造で同一条件を保証。"""
    estimator = Estimator(calc_mode=calc_mode, model_version=model_version)
    return ASECalculator(estimator)


def optimize_structure(input_path: Path, output_dir: Path, calculator: ASECalculator) -> Tuple[str, float, str]:
    """単一構造の最適化。例外は呼び出し側でハンドリング。"""
    atoms = read(str(input_path))
    atoms.calc = calculator

    opt = FIRE(atoms, logfile=None)
    converged = opt.run(fmax=0.05, steps=500)

    energy = atoms.get_potential_energy()
    status = "converged" if converged else "max_steps"

    optimized_path = output_dir / f"{input_path.stem}_opt.traj"
    write(str(optimized_path), atoms)
    return input_path.name, energy, status


def run_screening(input_dir: str = "inputs", run_root: str = "runs/screening"):
    """スクリーニングパイプラインのメイン関数。"""
    input_dir_path = Path(input_dir)
    run_root_path = Path(run_root)
    run_root_path.mkdir(parents=True, exist_ok=True)

    result_dir = run_root_path / "optimized_structures"
    result_dir.mkdir(exist_ok=True)

    csv_path = run_root_path / "energy_ranking.csv"
    log_path = run_root_path / "errors.log"

    calculator = build_calculator()
    records: List[Tuple[str, float, str]] = []

    input_files = sorted(list(input_dir_path.glob("*.cif")) + list(input_dir_path.glob("*.traj")))
    for input_file in input_files:
        try:
            record = optimize_structure(input_file, result_dir, calculator)
            records.append(record)
        except Exception as exc:
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(f"{input_file.name}: {exc}\n")

    # エネルギーでソートしてランキング CSV を出力
    records.sort(key=lambda x: x[1])
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "energy_eV", "status"])
        writer.writerows(records)

    print(f"Done. ranked={len(records)} errors_log={log_path}")


if __name__ == "__main__":
    run_screening()
```

**ポイント**:
- `build_calculator` で統一された calculator を生成し、全構造で同一条件を保証します
- 例外は `try/except` で捕捉し、`errors.log` に記録して次の構造に進みます
- ランキング CSV はエネルギー順にソートして毎回再生成します

### パターン B: 吸着構造探索パイプライン

**目的**: 表面構造に吸着種をランダムに配置し、初期配置の力に基づく早期棄却 → 最適化 → ランキングの流れで安定吸着構造を探索します。

**入力**: 表面構造 + 吸着種 + 配置候補数
**出力**: フィルタリング後の最適化構造 + ワークフロー結果 CSV

```python
from pathlib import Path
import csv
from copy import deepcopy
from typing import List, Dict

import numpy as np
from ase.io import read, write
from ase.build import add_adsorbate
from ase.optimize import FIRE
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator


def setup_calc() -> ASECalculator:
    return ASECalculator(Estimator(calc_mode="pbe_plus_d3", model_version="latest"))


def generate_candidates(slab, adsorbate, n_candidates: int = 20) -> List:
    """ランダムに吸着構造候補を生成する。"""
    candidates = []
    rng = np.random.default_rng(42)
    for _ in range(n_candidates):
        trial = slab.copy()
        x_frac, y_frac = rng.random(), rng.random()
        xy = np.matmul([x_frac, y_frac, 0.0], trial.cell)[:2]
        height = 1.8 + 1.2 * rng.random()
        add_adsorbate(trial, deepcopy(adsorbate), height=height, position=xy)
        candidates.append(trial)
    return candidates


def filter_by_initial_force(candidates: List, calculator: ASECalculator, max_force_threshold: float = 8.0) -> List:
    """初期の力が閾値以下の候補のみを残す（早期棄却フィルタ）。"""
    survived = []
    for atoms in candidates:
        atoms.calc = calculator
        max_force = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
        if max_force <= max_force_threshold:
            survived.append(atoms)
    return survived


def run_workflow(
    slab_path: str = "inputs/slab.traj",
    adsorbate_path: str = "inputs/adsorbate.traj",
    output_root: str = "runs/adsorption_workflow",
):
    out_root = Path(output_root)
    out_root.mkdir(parents=True, exist_ok=True)
    traj_dir = out_root / "survivor_structures"
    traj_dir.mkdir(exist_ok=True)

    slab = read(slab_path)
    ads = read(adsorbate_path)
    calc = setup_calc()

    # 候補生成 → 早期棄却 → 最適化
    candidates = generate_candidates(slab, ads, n_candidates=30)
    survivors = filter_by_initial_force(candidates, calc, max_force_threshold=8.0)

    rows: List[Dict] = []
    for idx, atoms in enumerate(survivors):
        atoms.calc = calc
        opt = FIRE(atoms, logfile=None)
        status = "converged" if opt.run(fmax=0.08, steps=300) else "max_steps"
        energy = atoms.get_potential_energy()

        file_path = traj_dir / f"candidate_{idx:03d}.traj"
        write(str(file_path), atoms)
        rows.append({"candidate": idx, "energy_eV": energy, "status": status, "path": str(file_path)})

    # エネルギーでソートしてサマリ CSV を出力
    rows = sorted(rows, key=lambda x: x["energy_eV"])
    with open(out_root / "adsorption_summary.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["candidate", "energy_eV", "status", "path"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"candidates={len(candidates)} survivors={len(survivors)}")


if __name__ == "__main__":
    run_workflow()
```

**ポイント**:
- `filter_by_initial_force` で初期の力が大きい候補を早期に除外します。これにより、不適切な初期配置での無駄な最適化を回避し、全体処理時間を大幅に短縮できます
- 吸着計算では `calc_mode="pbe_plus_d3"` のように D3 分散補正を含むモードを使用してください（分子間相互作用が重要なため）
- 乱数シード (`seed=42`) を固定することで再現性を確保しています
- 吸着構造探索では `TH_max_f`, `TH_min_f` による初期の力の閾値と結合長変化の閾値 `tol` も活用できます
- 単一系での吸着エネルギー計算やサイト探索の基礎は `adsorption/SKILL.md` を参照してください

### パターン C: 欠陥パイプライン

**目的**: 母相構造から空孔欠陥を生成し、2 段階緩和を経て欠陥エネルギーを比較します。

**入力**: 母相構造 + 欠陥候補インデックス
**出力**: 各欠陥の 2 段階緩和結果 + 比較 CSV

```python
from pathlib import Path
import csv
from typing import List, Dict

from ase.io import read, write
from ase.optimize import FIRE
from ase.constraints import FixAtoms
from ase.geometry import get_distances
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator


def create_calc() -> ASECalculator:
    return ASECalculator(Estimator(calc_mode="pbe", model_version="latest"))


def make_vacancy_candidate(base_atoms, remove_index: int):
    """空孔候補を生成する。除去した原子の座標（欠陥中心）も返す。"""
    atoms = base_atoms.copy()
    defect_center = atoms.positions[remove_index].copy()
    del atoms[remove_index]
    return atoms, defect_center


def relax_two_stage(atoms, calculator: ASECalculator, defect_center, r_local: float = 4.0):
    """2段階緩和: まず欠陥近傍のみ局所緩和、次にセル固定で内部座標を緩和。

    欠陥中心からの距離で固定領域を決める。セルは固定し内部座標のみ最適化する。
    """
    atoms.calc = calculator

    # Stage 1: 欠陥中心から r_local 以内のみ動かす（遠方を固定して局所緩和）
    _, dlen = get_distances(atoms.positions, defect_center.reshape(1, 3),
                            cell=atoms.cell, pbc=atoms.pbc)
    dists = dlen.reshape(-1)
    far_indices = [i for i, d in enumerate(dists) if d > r_local]
    atoms.set_constraint(FixAtoms(indices=far_indices))
    stage1 = FIRE(atoms, logfile=None)
    stage1.run(fmax=0.10, steps=200)

    # Stage 2: 拘束を解除し、セル固定のまま内部座標を全体緩和
    atoms.set_constraint()
    stage2 = FIRE(atoms, logfile=None)
    converged = stage2.run(fmax=0.05, steps=400)
    return converged, atoms.get_potential_energy()


def run_defect_workflow(
    bulk_path: str = "inputs/bulk.traj",
    defect_indices: List[int] = None,
    output_root: str = "runs/defect_workflow",
):
    if defect_indices is None:
        defect_indices = [0, 1, 2]

    out_root = Path(output_root)
    out_root.mkdir(parents=True, exist_ok=True)
    structure_dir = out_root / "defect_structures"
    structure_dir.mkdir(exist_ok=True)

    base = read(bulk_path)
    calc = create_calc()

    summary: List[Dict] = []
    for idx in defect_indices:
        try:
            defect_atoms, defect_center = make_vacancy_candidate(base, idx)
            converged, energy = relax_two_stage(defect_atoms, calc, defect_center)

            out_file = structure_dir / f"vacancy_{idx:03d}.traj"
            write(str(out_file), defect_atoms)

            summary.append({
                "vacancy_index": idx,
                "energy_eV": energy,
                "status": "converged" if converged else "max_steps",
                "path": str(out_file),
            })
        except Exception as exc:
            summary.append({
                "vacancy_index": idx,
                "energy_eV": None,
                "status": f"failed: {exc}",
                "path": "",
            })

    # 成功行と失敗行を分離してソート
    valid_rows = [r for r in summary if isinstance(r["energy_eV"], (int, float))]
    valid_rows.sort(key=lambda x: x["energy_eV"])
    invalid_rows = [r for r in summary if r not in valid_rows]

    with open(out_root / "defect_summary.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["vacancy_index", "energy_eV", "status", "path"])
        writer.writeheader()
        writer.writerows(valid_rows + invalid_rows)

    print(f"done valid={len(valid_rows)} failed={len(invalid_rows)}")


if __name__ == "__main__":
    run_defect_workflow()
```

**ポイント**:
- 2 段階緩和により、欠陥周辺のみを先に緩和してから全体を緩和します。これにより収束が安定しやすくなります
- 対象はバルクの結晶です。欠陥中心からの距離 `r_local` で固定領域を決めます
- Stage 1 では欠陥中心から `r_local`（既定 4.0 Å）より遠い原子を `FixAtoms` で固定し、`fmax=0.10` で近傍のみ粗く緩和します
- Stage 2 では拘束を解除し、セルは固定したまま内部座標のみを `fmax=0.05` で精密に全体緩和します
- 失敗した候補は `energy_eV=None`, `status=f"failed: {exc}"` としてサマリに残し、全体を停止しません

## ベストプラクティス

### 実行順序を固定する

「生成 → 棄却 → 緩和 → 保存 → 集計」の順序をテンプレート化すると、失敗が減ります。各パイプラインで共通のパターンを維持してください。

### ログと成果物を分離する

`errors.log` と `*_summary.csv` を別管理し、再実行時の差分確認を容易にしてください。エラーログには失敗した構造名と例外メッセージを記録し、サマリ CSV には成功した計算結果のランキングを出力します。

### ランキングは常に再計算する

中間ファイルを再利用する場合でも、最後にサマリ CSV を再ソートして出力してください。古いランキングが残っていると混乱の原因になります。

### 段階緩和を使い分ける

欠陥・吸着では局所緩和 → 全体緩和の 2 段階が安定しやすいです。単一ステージの緩和では、遠方の原子が不必要に動いて収束が遅くなる場合があります。

### 資源制約を前提に設計する

不適切な候補を最初に落とすだけで、全体処理時間を大きく短縮できます。吸着パイプラインの `filter_by_initial_force` のように、安価な事前チェックで明らかに不適切な候補を除外してください。

### 統一された計算条件

同一ワークフロー内では `calc_mode`, `model_version`, 収束条件（`fmax`, `steps`）を統一してください。条件が異なると横比較ができなくなります。

## よくあるエラーと対処

| エラー / 問題 | 原因 | 対処 |
|---|---|---|
| 入力ファイルが見つからない | パスの誤り / ファイル未配置 | `Path.exists()` で事前検証してください。ファイル拡張子も確認してください |
| 最適化が `max_steps` で終了する | 収束しない構造 | `steps` を増やすか、初期構造を見直してください。不安定な構造は早期棄却の対象にしてください |
| メモリ不足 | 大量の構造を同時にメモリに保持 | 構造を逐次処理し、結果はファイルに逐次保存してください |
| CSV に `None` が含まれる | 計算失敗 | 正常動作です。`valid_rows` と `invalid_rows` を分離して処理してください |
| 異なる条件で比較してしまう | calculator の設定不統一 | `build_calculator` / `setup_calc` / `create_calc` を一箇所で定義し、全体で使い回してください |
| 再実行時に結果が上書きされる | 出力先が同じ | `runs/<workflow_name>_<timestamp>/` のようにタイムスタンプ付きで出力先を分けてください |

## 関連ガイド

- **吸着** (`adsorption/SKILL.md`): 吸着構造の構築・吸着エネルギー計算・サイト探索の基礎
- **バックグラウンドジョブ** (`background-job/SKILL.md`): 大量のスクリーニングをバックグラウンドで実行する方法
- **pfcc_extras** (`pfcc-extras/SKILL.md`): `run_jobs` によるパイプラインの並列実行、構造の可視化・加工
- **SSH 接続** (`ssh/SKILL.md`): リモートからのジョブ監視・ファイル転送
