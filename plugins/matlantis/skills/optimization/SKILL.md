---
name: mt-optimization
description: >
  構造最適化（構造緩和）を扱うスキルです。
  LBFGS, BFGS, FIRE, fmax, 構造最適化, 構造緩和, geometry optimization,
  FrechetCellFilter, ExpCellFilter, UnitCellASEFilter, セル最適化, cell optimization,
  FixSymmetry, 対称性保持, maxstep,
  BFGSASEOptFeature, FireASEOptFeature, LBFGSASEOptFeature, filter=True,
  多段階最適化, 収束基準, trajectory
  に関するコード生成時に使用してください。
---
# 構造最適化

## 概要

原子にかかる力が閾値 (`fmax`) 以下になるまで原子位置およびセル形状を動かし、ポテンシャルエネルギー曲面上の極小点（安定構造）を探索するプロセスです。このガイドでは、固定セル最適化、可変セル最適化、対称性保持最適化、局所緩和、matlantis-features の高レベル API、およびオプティマイザの選択基準を統合的に解説します。

構造最適化はほぼすべてのシミュレーションワークフローの前提工程です。モデリングで生成した構造は平衡状態にないため、物性値の計算前に必ず最適化を行ってください。

## ワークフロー

```
1. Prepare    - 初期構造 (atoms) と Calculator (calc) を用意
2. Choose     - 目的に応じてオプティマイザとフィルタを選択
3. Constrain  - 必要に応じて拘束条件を設定（固定原子、対称性保持など）
4. Optimize   - 緩和計算を実行
5. Validate   - 収束確認、構造の可視化チェック
6. Save       - 最適化済み構造とトラジェクトリを保存
```

## 実装パターン

### パターン A: 固定セル最適化 (原子位置のみ)

格子定数を固定し、原子位置のみを最適化します。分子（孤立系）、表面モデル（スラブ）、欠陥を含むスーパーセル、MD 前の初期緩和に適しています。

```python
from ase.optimize import LBFGS
from ase.io import write

def run_fixed_cell_optimization(
    atoms,
    fmax: float = 0.05,
    steps: int = 1000,
    trajectory_file: str = "opt_fixed.traj"
):
    """格子定数を固定して原子位置のみを最適化する。"""
    if atoms.calc is None:
        raise RuntimeError("Calculator is not attached to atoms.")

    opt = LBFGS(atoms, trajectory=trajectory_file, logfile="opt_fixed.log")
    opt.run(fmax=fmax, steps=steps)

    epot = atoms.get_potential_energy()
    print(f"Optimization finished. Energy: {epot:.4f} eV")

    write("optimized_structure.cif", atoms)
    return atoms
```

分子の最小構成:

```python
from ase.optimize import LBFGS

atoms.calc = calculator
opt = LBFGS(atoms, trajectory="opt.traj")
opt.run(fmax=0.01)
```

### パターン B: 可変セル最適化 (原子位置 + セル)

原子位置と格子ベクトル（セル形状・体積）を同時に最適化します。バルク結晶の格子定数決定、相転移の探索に適しています。

ASE の Filter を使って Atoms をラップすることで、Optimizer がセル変形を扱えるようになります。

#### FrechetCellFilter を使用する場合

```python
from ase.filters import FrechetCellFilter
from ase.optimize import LBFGS
from ase.io import write

def run_variable_cell_optimization(
    atoms,
    fmax: float = 0.05,
    steps: int = 1000,
    trajectory_file: str = "opt_cell.traj"
):
    """原子位置と格子定数 (セル) を同時に最適化する。"""
    if atoms.calc is None:
        raise RuntimeError("Calculator is not attached to atoms.")

    atoms_filter = FrechetCellFilter(atoms)
    opt = LBFGS(atoms_filter, trajectory=trajectory_file, logfile="opt_cell.log")
    opt.run(fmax=fmax, steps=steps)

    epot = atoms.get_potential_energy()
    cell = atoms.get_cell_lengths_and_angles()
    print(f"Energy: {epot:.4f} eV")
    print(f"Cell: {cell}")

    write("optimized_cell_structure.cif", atoms)
    return atoms
```

#### ExpCellFilter を使用する場合（簡潔パターン）

```python
from ase.filters import ExpCellFilter
from ase.optimize import FIRE

atoms.calc = calculator
filtered = ExpCellFilter(atoms)
opt = FIRE(filtered, trajectory="opt.traj")
opt.run(fmax=0.01)
```

### パターン C: 対称性保持最適化

結晶の空間群対称性を維持したまま最適化を行います。

```python
from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter
from ase.optimize import FIRE

atoms.calc = calculator
atoms.set_constraint(FixSymmetry(atoms))
filtered = ExpCellFilter(atoms)
opt = FIRE(filtered, trajectory="opt.traj")
opt.run(fmax=0.01)
```

`FixSymmetry` は対称性操作を拘束として適用し、最適化中に対称性が破れることを防ぎます。

### パターン D: 多段階最適化

収束が困難な系では、粗い閾値で大まかに最適化してから精密な最適化を行います。

```python
from ase.optimize import LBFGS

# Stage 1: 粗い最適化
opt1 = LBFGS(atoms, trajectory="opt_stage1.traj", logfile="stage1.log")
opt1.run(fmax=0.1, steps=500)
print("Stage 1 complete")

# Stage 2: 精密な最適化
opt2 = LBFGS(atoms, trajectory="opt_stage2.traj", logfile="stage2.log")
opt2.run(fmax=0.05, steps=1000)
print("Stage 2 complete")
```

これは初期構造が平衡から遠い場合や、第一ステップで大きな力が働く場合に特に有効です。

### パターン E: maxstep による爆発防止

初期構造が不安定な場合、1 ステップの最大移動量を制限して構造の崩壊を防ぎます。

```python
from ase.optimize import LBFGS

# maxstep で 1 ステップあたりの最大変位を制限
opt = LBFGS(atoms, trajectory="opt.traj", maxstep=0.2)
opt.run(fmax=0.05)
```

`maxstep` のデフォルトは通常 0.04A 程度ですが、不安定な系では 0.1-0.2 に明示的に設定すると安全です。

### パターン F: 局所クラスタ最適化

大規模系で欠陥周辺のみを緩和する際に、計算コストを大幅に削減できます。

```python
from ase.constraints import FixAtoms

# 中心 (center) から r_fix 以内の原子のみ可動、外側は固定
center_pos = atoms.positions[defect_index]
distances = atoms.get_distances(defect_index, range(len(atoms)), mic=True)
indices_fix = [i for i, d in enumerate(distances) if d > r_fix]

atoms.set_constraint(FixAtoms(indices=indices_fix))

opt = LBFGS(atoms, trajectory="opt_local.traj")
opt.run(fmax=0.05)
```

運用指針:
- 最終報告構造: `full` 最適化（全原子緩和）
- 探索中の多数試行: `cluster` / `in_sphere` 最適化（局所緩和）

### パターン G: matlantis-features による高レベル最適化

matlantis-features は構造最適化のための高レベル API を提供します。

```python
from matlantis_features.features.common.opt import BFGSASEOptFeature
from matlantis_features.utils.calculators import pfp_estimator_fn

estimator_fn = pfp_estimator_fn(model_version="v8.0.0")

# 固定セル最適化
opt = BFGSASEOptFeature(
    fmax=0.01,
    n_run=1000,
    estimator_fn=estimator_fn,
)
result = opt(atoms)
optimized_atoms = result.atoms.ase_atoms
print(f"Converged: {result.converged}")
```

利用可能なオプティマイザ:

| Feature クラス | アルゴリズム |
|---------------|------------|
| `BFGSASEOptFeature` | BFGS |
| `LBFGSASEOptFeature` | L-BFGS |
| `FireASEOptFeature` | FIRE |
| `FireLBFGSASEOptFeature` | FIRE + L-BFGS 切り替え |

#### セル最適化 (matlantis-features)

結晶のセル最適化には `filter` パラメータを使用します。

```python
from matlantis_features.features.common.opt import BFGSASEOptFeature
from matlantis_features.filters import UnitCellASEFilter

opt = BFGSASEOptFeature(
    fmax=0.01,
    filter=UnitCellASEFilter(),  # または filter=True
    estimator_fn=estimator_fn,
)
result = opt(atoms)
```

`filter=True` でもセル最適化が有効になります。より詳細な制御には `UnitCellASEFilter()` を直接渡します。

### パターン H: サーバーエラーからの再開

`PFPAPIError` や `RetriesExceeded` が発生しても、Atoms オブジェクトは最後の座標を保持しています。

```python
from ase.optimize import LBFGS

opt = LBFGS(atoms, trajectory="opt.traj")

try:
    opt.run(fmax=0.05, steps=1000)
except Exception as e:
    print(f"Error: {e}")
    print("Retrying from last position...")
    # Optimizer を再作成せず、そのまま再実行（前回のステップから再開）
    opt.run(fmax=0.05, steps=1000)
```

## ベストプラクティス

### オプティマイザの選択基準

| オプティマイザ | 推奨用途 | 備考 |
|--------------|---------|------|
| `LBFGS` | 分子の最適化、初心者向け | メモリ効率が良い。デフォルト選択 |
| `BFGS` | 一般的な最適化 | 安定性が高い |
| `FIRE` | 結晶のセル最適化、NEB | 大きな力に強い |

1 ページで複数のオプティマイザを並べすぎないでください。用途に応じて 1 つを選択するのが基本です。

### 収束判定基準 (fmax)

| 用途 | 推奨 fmax (eV/A) |
|------|-----------------|
| 粗い最適化（探索用） | 0.1 |
| 通常用途 | 0.05 |
| 高精度（弾性定数、格子定数） | 0.01 |
| フォノン計算の前処理 | 0.001 - 0.05 |
| matlantis-features のデフォルト | 0.01 |

注意: fmax を厳しくしすぎると数値誤差により収束しない場合があります。`0.001` 以下は通常不要です。

### その他の推奨事項

1. **多段階最適化**: 収束が困難な系では、`fmax=0.1` で粗く最適化してから `fmax=0.05` で精密化してください。

2. **対称性の保持**: PFP 計算自体は対称性を厳密に保持しません。対称性維持が重要な場合は `FixSymmetry` 拘束を使用してください。

3. **可視化による確認**: 最適化前後は必ず構造を可視化してください。意図しない構造変化（結合の開裂、アモルファス化）がないか確認が必要です。

4. **冪等な再実行**: `opt.run()` は途中で中断しても、同じ Atoms オブジェクトで再実行すれば前回の位置から継続します。

5. **フォノン計算の前提**: `ForceConstantFeature` を使用する前に、`fmax < 0.05` まで最適化してください。

## よくあるエラーと対処

| エラー・問題 | 原因 | 対処 |
|-------------|------|------|
| `Too many neighbors` | 初期構造の原子間距離が近すぎる（重なり） | モデリングを見直すか、微小ランダム変位 (Rattle) を加えて重なりを解消してください |
| `Cell is too small for book keeping` | セルサイズが PFP のカットオフ半径に対して不足 | スーパーセルを作成してください (`atoms *= (2,2,2)`) |
| `PFPAPIError` / `RetriesExceeded` | サーバー通信エラー | `opt.run()` をそのまま再実行してください（前回の位置から再開されます） |
| 構造が崩れる (Explosion) | 初期構造が物理的にありえない配置 | `maxstep=0.2` で 1 ステップの最大移動量を制限してください |
| 収束しない | fmax が厳しすぎる、または局所最小に陥っている | 多段階最適化を試すか、fmax を緩めてください |
| `RuntimeError: Calculator is not attached` | Calculator が未設定 | `atoms.calc = calculator` を実行してください |
| 対称性が破れる | PFP は対称性を厳密に保持しない | `FixSymmetry` 拘束を追加してください |

## 関連ガイド

- **Calculator 初期化** (setup/SKILL.md): 最適化に必須の前工程
- **構造モデリング** (modeling/SKILL.md): 最適化対象の構造作成
- **分子動力学** (dynamics/SKILL.md): 最適化後の動的シミュレーション
