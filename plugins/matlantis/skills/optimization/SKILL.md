---
name: mt-optimization
description: >
  構造最適化（構造緩和）を扱うスキルです。
  LBFGS, LBFGSLineSearch, BFGS, BFGSLineSearch, QuasiNewton, FIRE, FIRE2, ABC-FIRE,
  GoodOldQuasiNewton, CellAwareBFGS, RFO, GPMin, MDMin,
  PreconLBFGS, PreconFIRE, SciPyFminBFGS, SciPyFminCG,
  BasinHopping, MinimaHopping,
  fmax, 構造最適化, 構造緩和, geometry optimization,
  FrechetCellFilter, ExpCellFilter, UnitCellASEFilter, セル最適化, cell optimization,
  FixSymmetry, 対称性保持, maxstep,
  ASEOptFeature, BFGSASEOptFeature, BFGSLineSearchASEOptFeature,
  LBFGSASEOptFeature, LBFGSLineSearchASEOptFeature, FireASEOptFeature,
  FireLBFGSASEOptFeature, MDMinASEOptFeature, filter=True,
  多段階最適化, 収束基準, trajectory, preconditioned optimizer,
  グローバル最適化, global optimization
  に関するコード生成時に使用してください。
---
# 構造最適化

## 概要

原子にかかる力が閾値 (`fmax`) 以下になるまで原子位置およびセル形状を動かし、ポテンシャルエネルギー曲面上の極小点（安定構造）を探索するプロセスです。構造最適化はほぼすべてのシミュレーションワークフローの前提工程です。モデリングで生成した構造は平衡状態にないため、物性値の計算前に必ず最適化を行ってください。以下のサンプルコード中で事前に生成していないオブジェクトに関しては事前に生成しておくようにすること。

## ワークフロー

```
1. Prepare    - 初期構造 (atoms) と Calculator (calc) を用意
2. Choose     - 目的に応じてオプティマイザとフィルタを選択
3. Constrain  - 必要に応じて拘束条件を設定（固定原子、対称性保持など）
4. Optimize   - 緩和計算を実行
5. Validate   - 収束確認、構造の可視化チェック
6. Save       - 最適化済み構造とトラジェクトリを保存
```

## 使用例

### 固定セル最適化（原子位置のみ）

格子定数を固定し、原子位置のみを最適化します。分子、スラブ、欠陥スーパーセル、MD 前の初期緩和に適しています。

```python
from ase.optimize import LBFGS
from ase.io import write

atoms.calc = calculator
opt = LBFGS(atoms, trajectory="opt.traj", logfile="opt.log")
opt.run(fmax=0.05, steps=1000)

print(f"Energy: {atoms.get_potential_energy():.4f} eV")
write("optimized.cif", atoms)
```

### 可変セル最適化（原子位置 + セル）

原子位置と格子ベクトルを同時に最適化します。バルク結晶の格子定数決定に適しています。Filter で Atoms をラップすることでオプティマイザがセル変形を扱えるようになります。

```python
from ase.filters import FrechetCellFilter
from ase.optimize import LBFGS
from ase.io import write

atoms.calc = calculator
filtered = FrechetCellFilter(atoms)
opt = LBFGS(filtered, trajectory="opt.traj", logfile="opt.log")
opt.run(fmax=0.01, steps=1000)

print(f"Energy: {atoms.get_potential_energy():.4f} eV")
print(f"Cell: {atoms.cell.cellpar()}")   # [a, b, c, alpha, beta, gamma]
write("optimized.cif", atoms)
```

`ExpCellFilter` を使う場合（※新しい ASE では `ExpCellFilter` は deprecated で、セル変数に対する収束性の良い `FrechetCellFilter` の使用が推奨されます）:

```python
from ase.filters import ExpCellFilter
from ase.optimize import FIRE

atoms.calc = calculator
filtered = ExpCellFilter(atoms)
opt = FIRE(filtered, trajectory="opt.traj")
opt.run(fmax=0.01)
```

### 対称性保持最適化

`FixSymmetry` 拘束を加えることで、空間群対称性を維持したまま最適化できます。PFP 自体は対称性を厳密に保持しないため、必要な場合は明示的に設定してください。

```python
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import LBFGS

atoms.calc = calculator
atoms.set_constraint(FixSymmetry(atoms))
filtered = FrechetCellFilter(atoms)
opt = LBFGS(filtered, trajectory="opt.traj")
opt.run(fmax=0.01)
```

### 多段階最適化

初期構造が平衡から遠い場合や大きな力が働く場合は、粗い閾値で大まかに最適化してから精密化します。

```python
from ase.optimize import FIRE, LBFGS

opt = FIRE(atoms, trajectory="opt_stage1.traj", logfile="stage1.log")
opt.run(fmax=0.1, steps=500)

opt = LBFGS(atoms, trajectory="opt_stage2.traj", logfile="stage2.log")
opt.run(fmax=0.05, steps=1000)
```

### maxstep による爆発防止

初期構造が不安定な場合、`maxstep` で 1 ステップの最大変位（Å）を制限して構造の崩壊を防ぎます。デフォルトは 0.2 Åなので、崩壊を防ぐにはそれより小さい値（例: 0.05 Å 程度）を指定します。

```python
from ase.optimize import LBFGS

opt = LBFGS(atoms, trajectory="opt.traj", maxstep=0.05)  # デフォルト0.2より小さく制限
opt.run(fmax=0.05)
```

### 局所クラスタ最適化（欠陥周辺のみ緩和）

大規模系で欠陥周辺のみを緩和し、計算コストを削減します。

```python
from ase.constraints import FixAtoms
from ase.optimize import LBFGS

distances = atoms.get_distances(defect_index, range(len(atoms)), mic=True)
indices_fix = [i for i, d in enumerate(distances) if d > r_fix]
atoms.set_constraint(FixAtoms(indices=indices_fix))

opt = LBFGS(atoms, trajectory="opt_local.traj")
opt.run(fmax=0.05)
```

最終構造は全原子緩和、探索中の多数試行は局所緩和が推奨です。

### matlantis-features による高レベル最適化

matlantis-features は構造最適化のための高レベル API を提供します。

```python
from matlantis_features.features.common.opt import BFGSASEOptFeature
from matlantis_features.utils.calculators import pfp_estimator_fn

estimator_fn = pfp_estimator_fn(model_version="v8.0.0")

opt = BFGSASEOptFeature(fmax=0.01, n_run=1000, estimator_fn=estimator_fn)
result = opt(atoms)
optimized_atoms = result.atoms.ase_atoms
print(f"Converged: {result.converged}")
```

セル最適化には `filter` パラメータを使用します。

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

利用可能なオプティマイザ（すべて `matlantis_features.features.common.opt` にあり、共通引数 `n_run`（最大ステップ数、デフォルト 200）、`fmax`、`filter`、`estimator_fn`、`attach_methods`、`show_progress_bar`、`show_logger` を取ります）:

| Feature クラス | アルゴリズム |
|---------------|------------|
| `BFGSASEOptFeature` | BFGS |
| `BFGSLineSearchASEOptFeature` | BFGS + ラインサーチ |
| `LBFGSASEOptFeature` | L-BFGS |
| `LBFGSLineSearchASEOptFeature` | L-BFGS + ラインサーチ |
| `FireASEOptFeature` | FIRE |
| `FireLBFGSASEOptFeature` | FIRE + L-BFGS 切り替え |
| `MDMinASEOptFeature` | MDMin |
| `ASEOptFeature` | 任意の ASE Optimizer をラップ（`get_optimizer` を指定） |

`ASEOptFeature` は上記以外の ASE オプティマイザ（`GPMin`, `CellAwareBFGS`, `PreconLBFGS` など）を matlantis-features として使うための汎用ラッパーです。`get_optimizer`（`ASEAtoms` を受け取り `Optimizer` を返す関数）を渡します。

```python
from matlantis_features.features.common.opt import ASEOptFeature
from matlantis_features.utils.calculators import pfp_estimator_fn
from ase.optimize import GPMin

estimator_fn = pfp_estimator_fn(model_version="v8.0.0")

def get_optimizer(atoms):
    return GPMin(atoms, logfile=None)

opt = ASEOptFeature(get_optimizer=get_optimizer, fmax=0.01, n_run=1000, estimator_fn=estimator_fn)
result = opt(atoms)
```

> 注: 各サブクラス（`BFGSASEOptFeature` など）の `fmax` デフォルトは 0.01、基底の `ASEOptFeature` のみ `fmax` デフォルトが 0.001 です。

### サーバーエラーからの再開

`PFPAPIError` や `RetriesExceeded` が発生しても Atoms オブジェクトは最後の座標を保持しています。Optimizer を再作成せずそのまま再実行すれば前回の位置から継続します。

```python
from ase.optimize import LBFGS

opt = LBFGS(atoms, trajectory="opt.traj")
try:
    opt.run(fmax=0.05, steps=1000)
except Exception as e:
    print(f"Error: {e}. Retrying from last position...")
    opt.run(fmax=0.05, steps=1000)
```

## オプティマイザ別リファレンス

### LBFGSLineSearch / BFGSLineSearch（QuasiNewton）

Wolfe 条件（エネルギー単調減少 + 十分な勾配減少）を満たすようにステップ長を決定するラインサーチ付きオプティマイザです。標準版より収束安定性が向上しますが、1 ステップあたりの Calculator 呼び出しが増えることがあります。

```python
from ase.optimize import LBFGSLineSearch

# 方向のみ Hessian で決め、ステップ長は力から決定
opt = LBFGSLineSearch(atoms, trajectory="opt.traj", logfile="opt.log")
opt.run(fmax=0.05)
```

```python
from ase.optimize import BFGSLineSearch  # QuasiNewton という別名でもインポート可

# Wolfe 条件ラインサーチ付き BFGS。Hessian は pickle で再起動可能
opt = BFGSLineSearch(atoms, trajectory="opt.traj", restart="opt.pckl")
opt.run(fmax=0.05)
```

⚠ `BFGSLineSearch` は NEB 計算には使用不可。NEB には `LBFGS` または `FIRE` を使うこと。

### FIRE2 / ABC-FIRE

FIRE の改良版。同じ `FIRE2` クラスで両方が利用でき、`use_abc=True` で ABC-FIRE（Accelerated Bias-Corrected FIRE）が有効になります。

```python
from ase.optimize import FIRE2

opt = FIRE2(atoms, trajectory="opt.traj")           # FIRE2.0
opt.run(fmax=0.05)

opt = FIRE2(atoms, trajectory="opt.traj", use_abc=True)  # ABC-FIRE
opt.run(fmax=0.05)
```

### GPMin（ガウス過程最小化）

ガウス過程回帰でポテンシャルエネルギー曲面を構築し、局所最小化を高速化します。

```python
from ase.optimize import GPMin

opt = GPMin(atoms, trajectory="opt.traj", logfile="opt.log")
opt.run(fmax=0.05)
```

⚠ メモリ使用量が O(n²N²)（n=ステップ数、N=原子数）でスケールするため、**100 原子以下の系**に限定してください。超過時は警告が表示されます。

### MDMin

速度 Verlet MD を改変した最適化アルゴリズム。運動量ベクトルと力ベクトルの内積がゼロ（局所最小通過）になった時点で全運動量をリセットします。物理的に自然な緩和が得られますが、二次収束領域では `LBFGS` への切り替えが効率的です。

```python
from ase.optimize import MDMin

opt = MDMin(atoms, trajectory="opt.traj", logfile="opt.log")
opt.run(fmax=0.05)
```

### GoodOldQuasiNewton

ASE の旧来の準ニュートン法。ステップ数が少なく収束することが多く、現代的なオプティマイザと遜色ない場面があります。

```python
from ase.optimize import GoodOldQuasiNewton

opt = GoodOldQuasiNewton(atoms, trajectory="opt.traj")
opt.run(fmax=0.05)
```

### CellAwareBFGS

バルク弾性率やポアソン比の情報をセル最適化に組み込める特殊 BFGS。`FrechetCellFilter` と組み合わせて使用します。

```python
from ase.optimize import CellAwareBFGS
from ase.filters import FrechetCellFilter

filtered = FrechetCellFilter(atoms)
opt = CellAwareBFGS(filtered, trajectory="opt.traj")
opt.run(fmax=0.05)
```

### RFO（有理関数オプティマイザ）

BFGS 型 Hessian 更新を使いながら、二次領域では準ニュートンステップ、それ以外では減衰した有理関数ステップを取ります。`damping`（座標スケールの逆数、単位 Å⁻¹）で領域切り替えの閾値を制御します（デフォルトは `None` で内部既定値が使われます）。

```python
from ase.optimize import RFO

opt = RFO(atoms, trajectory="opt.traj", damping=1.0)
opt.run(fmax=0.05)
```

⚠ `RFO` は比較的新しい ASE で追加されたクラスです。`ImportError` が出る場合は ASE のバージョンを確認してください。

### PreconLBFGS / PreconFIRE（Preconditioned オプティマイザ）

局所結合トポロジーの情報を座標変換に組み込み、収束を高速化します。大規模固体系で標準 LBFGS の最大 **10 倍** の高速化が報告されています。100 原子未満の系では効果が小さいため、標準 LBFGS/BFGS を推奨します。

```python
from ase.optimize.precon import PreconLBFGS, Exp

# 固体・金属系：Exp プリコンディショナー + Armijo ラインサーチ
opt = PreconLBFGS(atoms, precon=Exp(A=3), use_armijo=True, trajectory="opt.traj")
opt.run(fmax=1e-3)
```

```python
from ase.optimize.precon import PreconFIRE, Exp

opt = PreconFIRE(atoms, precon=Exp(A=3), trajectory="opt.traj")
opt.run(fmax=0.05)
```

利用可能なプリコンディショナー:

| クラス | 用途 |
|--------|------|
| `Exp` | 固体・金属系（汎用） |
| `FF` | 気相分子系 |
| `Exp_FF` | 分子結晶（`Exp` + `FF` の合成） |

### SciPy オプティマイザ

SciPy のオプティマイザを ASE インターフェース経由で利用します。`alpha` は Hessian の初期推定値（デフォルト 70.0）。

```python
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG

opt = SciPyFminBFGS(atoms, trajectory="opt.traj", alpha=70.0)  # BFGS
opt.run(fmax=0.05)

opt = SciPyFminCG(atoms, trajectory="opt.traj", alpha=70.0)    # 非線形共役勾配法
opt.run(fmax=0.05)
```

## グローバル最適化

### Basin Hopping

局所最適化と温度制御したモンテカルロサンプリングを組み合わせてグローバル最小を探索します。

```python
from ase.optimize.basin import BasinHopping
from ase.optimize import LBFGS
from ase.units import kB

bh = BasinHopping(
    atoms=atoms,
    temperature=100 * kB,  # 障壁を越えるための「温度」
    dr=0.5,                # 最大ステップ幅 (Å)
    optimizer=LBFGS,
    fmax=0.1,
)
bh.run(steps=50)
```

### Minima Hopping

NVE MD と局所最適化を交互に繰り返します。温度と受容閾値 (`Ediff`) を動的に調整します。

```python
from ase.optimize.minimahopping import MinimaHopping
from ase.optimize import BFGSLineSearch

opt = MinimaHopping(
    atoms=atoms,
    T0=1000.,       # 初期 MD 温度 (K)
    Ediff0=0.5,     # 初期エネルギー受容閾値 (eV)
    fmax=0.05,
    optimizer=BFGSLineSearch,
)
opt(totalsteps=10)
```

## ベストプラクティス

### オプティマイザの選択基準

| オプティマイザ | 推奨用途 | 備考 |
|--------------|---------|------|
| `LBFGS` | 分子・スラブ、初心者向け | メモリ効率が良い。デフォルト選択 |
| `LBFGSLineSearch` | LBFGS より収束安定性が欲しい場合 | Wolfe 条件。方向は Hessian、ステップ長は力 |
| `BFGS` | 一般的な最適化 | 安定性が高い |
| `BFGSLineSearch` (`QuasiNewton`) | 準ニュートン法の最良選択 | Wolfe 条件。NEB には不可 |
| `FIRE` | セル最適化、NEB、大変形系 | 大きな力に強い |
| `FIRE2` | FIRE より改善された収束性が必要な場合 | `use_abc=True` で ABC-FIRE |
| `MDMin` | シンプルな MD ベース緩和 | 物理的に自然。二次収束域では LBFGS に切替推奨 |
| `GoodOldQuasiNewton` | ステップ数を少なくしたい場合 | 旧来の準ニュートン法 |
| `CellAwareBFGS` | バルク弾性率・ポアソン比が既知のセル最適化 | `FrechetCellFilter` と組み合わせ必須 |
| `RFO` | 安定性を保ちつつ大変位を扱いたい場合 | `damping` で領域切り替え閾値を制御 |
| `GPMin` | 小系（≤100原子）の高効率最適化 | ガウス過程。メモリが O(n²N²) でスケール |
| `PreconLBFGS` | 大規模固体（100原子超）の高速化 | 最大 10x 高速化。100原子未満は標準 LBFGS |
| `SciPyFminBFGS` / `SciPyFminCG` | SciPy との連携が必要な場合 | ASE ネイティブより遅いことが多い |

### 収束判定基準（fmax）

| 用途 | 推奨 fmax (eV/Å) |
|------|-----------------|
| 粗い最適化（探索用） | 0.1 |
| 通常用途 | 0.05 |
| 高精度（弾性定数、格子定数） | 0.01 |
| フォノン計算の前処理 | 0.001 - 0.05 |
| matlantis-features のデフォルト | 0.01 |

fmax を厳しくしすぎると数値誤差により収束しない場合があります。`0.001` 以下は通常不要です。

### 推奨事項

- **多段階最適化**: 収束が困難な系では `fmax=0.1` で粗く最適化してから `fmax=0.05` で精密化してください。
- **対称性の保持**: PFP 計算自体は対称性を厳密に保持しません。対称性維持が重要な場合は `FixSymmetry` 拘束を使用してください。
- **可視化による確認**: 最適化前後は必ず構造を可視化してください。結合の開裂やアモルファス化がないか確認が必要です。
- **冪等な再実行**: `opt.run()` は途中で中断しても、同じ Atoms オブジェクトで再実行すれば前回の位置から継続します。
- **フォノン計算の前提**: `ForceConstantFeature` を使用する前に `fmax < 0.05` まで最適化してください。

## よくあるエラーと対処

| エラー・問題 | 原因 | 対処 |
|-------------|------|------|
| `Too many neighbors` | 原子間距離が近すぎる（重なり） | モデリングを見直すか、微小ランダム変位 (Rattle) を加えてください |
| `Cell is too small for book keeping` | セルサイズが PFP のカットオフ半径に対して不足 | スーパーセルを作成してください (`atoms *= (2,2,2)`) |
| `PFPAPIError` / `RetriesExceeded` | サーバー通信エラー | `opt.run()` をそのまま再実行してください（前回の位置から再開されます） |
| 構造が崩れる (Explosion) | 初期構造が物理的にありえない配置 | `maxstep=0.2` で 1 ステップの最大移動量を制限してください |
| 収束しない | fmax が厳しすぎる、または局所最小に陥っている | 多段階最適化を試すか、fmax を緩めてください |
| `RuntimeError: Calculator is not attached` | Calculator が未設定 | `atoms.calc = calculator` を実行してください |
| 対称性が破れる | PFP は対称性を厳密に保持しない | `FixSymmetry` 拘束を追加してください |
| GPMin でメモリエラー | 原子数が多い（>100） | `LBFGS` または `LBFGSLineSearch` に切り替えてください |
| `BFGSLineSearch` で NEB が動かない | アルゴリズムの互換性問題 | NEB には `LBFGS` または `FIRE` を使用してください |

## 関連ガイド

- **Calculator 初期化** (setup/SKILL.md): 最適化に必須の前工程
- **構造モデリング** (modeling/SKILL.md): 最適化対象の構造作成
- **分子動力学** (dynamics/SKILL.md): 最適化後の動的シミュレーション
