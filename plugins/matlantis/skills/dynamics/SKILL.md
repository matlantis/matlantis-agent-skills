---
name: mt-dynamics
description: >
  分子動力学(MD)シミュレーションを扱うスキルです。
  Langevin, VelocityVerlet, NVTBerendsen, NPTBerendsen, NVE, NVT, NPT,
  MaxwellBoltzmannDistribution, MDLogger, timestep, friction, temperature_K,
  MDFeature, MDFeatureResult, ASEMDSystem, LangevinIntegrator, NPTIntegrator,
  チェックポイント, checkpoint,
  PostMDDiffusionFeature, PostMDSpecificHeatFeature,
  PostEMDViscosityFeature, PostNEMDViscosityFeature, PostNEMDThermalConductivityFeature,
  PLUMED, メタダイナミクス, metadynamics, FES, collective variable,
  熱分解, pyrolysis, get_mol_list, フラグメント解析, fragment analysis,
  DepositionScheduler, ApplyUniformEfield, ElasticVirtualWall,
  Matlantis-LAMMPS, LAMMPS, OpenMM, pair_style pfp_api,
  分子動力学, MD, 拡散係数, 粘度, 熱伝導率
  に関するコード生成時に使用してください。
---
# 分子動力学

## 概要

有限温度における原子・分子の時間発展をシミュレーションするガイドです。

MD 計算は **「準備 → 実行 → 解析」の 3 段**で構成されます。実行には 2 つのルート（A: ASE 標準 MD、B: matlantis-features の `MDFeature`）があり、**解析コードは使用したルートによって異なります**（A は `.traj` ファイルの直接解析、B は `MDFeatureResult` と Post-MD Feature）。

分子動力学 (MD) は拡散係数、粘性、熱伝導率の計算、アモルファス構造の作成（Melt-Quench 法）、相転移のシミュレーションなどに使用されます。MD の各ステップでトークンを消費するため、本番計算の前に短いテストランで消費量を見積もることを推奨します。

長時間の MD Notebook 実行では、foreground 実行だけに決め打ちせず background 実行も常に選択肢として提示してください。数時間以上の計算や離席を伴う場合は `background-job/SKILL.md` に従い `mtl-bg-job` を優先します。

本ガイドに記載のない熱浴・バロスタットや細かな MD 設定が必要な場合は、Matlantis-LAMMPS（`pair_style pfp_api` で PFP を LAMMPS から利用）や OpenMM との連携を検討してください（例えば SLLOD 法によるせん断 NEMD、rRESPA 等の多重時間刻み、グランドカノニカル MC、SHAKE/SETTLE による結合・剛体拘束、Monte Carlo バロスタットなどは ASE 標準 MD にはなく、LAMMPS / OpenMM 側の機能です）。LAMMPS 連携の実例は `matlantis/references/community-examples/uniaxial_tensile_lammps/` および `plumed_well_tempered_metadynamics/` を参照してください。

## ワークフロー

```
1. 準備段階
   構造最適化済みの atoms を用意 → Calculator 設定 → 初期速度の付与 → ロガー/トラジェクトリ設定
   ↓
2. 平衡化（equilibration）
   NVT / NPT で目標温度（NPT では密度も）に到達させ、温度・エネルギーのプラトーを確認
   ↓
3. 本計算（production; ルート A / B のいずれかを選択）
   ルート A: ASE 標準 MD（NVE / NVT / NPT を直接構築）
   ルート B: matlantis-features の MDFeature（Integrator + ASEMDSystem）
   ↓
4. 解析（使用したルートに対応; 本計算のデータのみを解析する）
   A → .traj ファイルの読み込みと解析
   B → MDFeatureResult + Post-MD Feature（輸送特性）
```

### アンサンブル選択ガイド

| 目的 | 推奨アンサンブル |
| :--- | :--- |
| 密度の緩和、有限温度での格子定数の決定 | NPT |
| 平衡状態のサンプリング、アニーリング | NVT（Langevin） |
| 拡散係数など動的性質の評価 | 平衡化（NVT）後に NVE、または弱い Langevin（friction~0.002） |
| 急速な平衡化 | NVTBerendsen / 強めの Langevin（friction~0.02） |
| 真空を含む系（slab・分子・クラスター） | NVE / NVT（**NPT は使用不可**。後述の警告参照） |

### ルート選択ガイド

| ルート | 特徴 | 向くケース |
| :--- | :--- | :--- |
| **A: ASE 標準 MD** | ASE の MD クラスを直接構築。attach による柔軟なカスタマイズ | 単発・短時間の MD、独自のロギングや制御を組み込みたい場合 |
| **B: MDFeature** | チェックポイント/リスタート対応、Post-MD Feature（輸送特性）と直結 | 長時間 MD、拡散係数・粘度・熱伝導率の算出が目的の場合 |

## 準備段階

MD を開始する前に、構造最適化済みの atoms（optimization/SKILL.md 参照）に Calculator を設定し、初期温度に応じた速度分布を付与します。

```python
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.io.trajectory import Trajectory
from ase import units

# Calculator を設定（構造最適化済みの atoms を使用する）
atoms.calc = calculator

# 初期速度の付与
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
Stationary(atoms)       # 重心運動（全体の並進）を除去

# 孤立系（分子・クラスター）の場合のみ、全体の回転も除去する
# 周期系（バルク・スラブ）では全体回転は定義されないため適用しない
# ZeroRotation(atoms)

# トラジェクトリ保存
traj = Trajectory("md.traj", "w", atoms)
```

**初期温度の低下に注意**: 構造最適化済み（ポテンシャル極小）の構造に `MaxwellBoltzmannDistribution(temperature_K=T)` で速度を与えると、エネルギー等分配により運動エネルギーの約半分がポテンシャルエネルギーに流れ、温度は T/2 程度まで低下します。NVT / NPT では熱浴が回復させますが、NVE では低温のまま進行します。NVE で目標温度にしたい場合は 2T で初期化するか、先に NVT で平衡化してから NVE に切り替えてください。

**再現性が必要な場合**は乱数生成器を明示的に渡します（速度の初期化と Langevin の乱雑力の両方が乱数を使います）。

```python
import numpy as np

rng = np.random.default_rng(42)
MaxwellBoltzmannDistribution(atoms, temperature_K=300, rng=rng)
# Langevin(atoms, ..., rng=rng) も同様に指定可能
```

ログ出力には `MDLogger` が利用できます（公式チュートリアル準拠）。温度・エネルギーの時系列をファイルに記録します。

```python
from ase.md import MDLogger

# dyn は後続の MD オブジェクト（VelocityVerlet, Langevin 等）
dyn.attach(
    MDLogger(dyn, atoms, "md.log", header=True, stress=False, peratom=False, mode="w"),
    interval=100
)
```

## MD 実行

### 平衡化 → 本計算の 2 段階プロトコル

実務の MD は 1 本の計算ではなく、**平衡化（equilibration）と本計算（production）の 2 段階**で実行します。平衡化が完了していないデータを解析に含めると、物性値に系統誤差が入ります。

1. **平衡化**: 強めの熱浴で目標温度（NPT では密度も）に到達させる。温度・ポテンシャルエネルギー（NPT では体積）の時系列がプラトーに達し、ドリフトがないことを確認する（解析セクションの時系列 plot を使用）
2. **本計算**: 平衡化済みの構造から、解析目的に合ったアンサンブル・弱い熱浴で実行する。**解析には本計算のデータのみを使用する**

ここで「強めの熱浴」とは、系と熱浴のカップリングを密にする（Langevin の `friction` を大きく、Berendsen の `taut` を小さくする）ことを指します。カップリングを強めると熱浴とのエネルギー交換が速くなり、運動温度が目標値へ収束する時定数（Langevin で ~1/friction、Berendsen で taut）が短くなるため、初期速度の付与や最適化構造から生じた余剰エネルギー（「初期温度の低下に注意」参照）が速やかに散逸し、温度のドリフトやオーバーシュートが抑えられます。ただし速くなるのは運動温度（熱）の緩和であって、拡散・密度・構造といった遅い自由度の緩和は加速されない点に注意してください（詳細は「サーモスタットの強弱と物性への影響」）。

```python
from ase import units
from ase.io import Trajectory
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

# --- 1. 平衡化: 強めの熱浴で目標温度に到達させる ---
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
Stationary(atoms)

eq_dyn = Langevin(
    atoms, timestep=1.0 * units.fs, temperature_K=300,
    friction=0.02,                       # 平衡化では強めの熱浴
    trajectory="equilibration.traj", loginterval=100,
)
eq_dyn.run(10000)  # 10 ps 程度から。プラトー未到達なら延長する

# --- 2. 本計算: 平衡化済み構造から開始し、このデータのみを解析する ---
prod_dyn = Langevin(
    atoms, timestep=1.0 * units.fs, temperature_K=300,
    friction=0.002,                      # 本計算では弱い熱浴（動的性質を乱さない）
    trajectory="production.traj", loginterval=10,
)
prod_dyn.run(100000)  # 必要な長さはベストプラクティスの目安表を参照
```

拡散係数など動的性質を厳密に評価する場合は、本計算を NVE（`VelocityVerlet`）に切り替えるとサーモスタットの影響を完全に排除できます。

### ルート A: ASE 標準 MD

#### NVE アンサンブル（VelocityVerlet）

体積 (V) とエネルギー (E) を一定に保つ、最も基本的なアンサンブルです。サーモスタットのバイアスがないため、動的性質（拡散係数など）の正確な評価に適しています。

```python
from ase.md.verlet import VelocityVerlet

dyn = VelocityVerlet(atoms, timestep=1.0 * units.fs)
dyn.attach(traj.write, interval=10)  # 10 ステップごとに保存
dyn.run(1000)
traj.close()
```

NVE ではエネルギー保存が正しいかを確認してください。保存が悪い場合はタイムステップを小さくします。

#### NVT アンサンブル（Langevin）

体積 (V) と温度 (T) を一定に保つシミュレーションです。最も安定しており、平衡状態のサンプリングや高温での構造探索（アニール）に適しています。

```python
from ase.md.langevin import Langevin
from ase import units
from ase.io import Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

def run_nvt_md(
    atoms,
    temp_k: float = 300.0,
    steps: int = 5000,
    dt_fs: float = 1.0,
    friction: float = 0.002,
    traj_file: str = "md_nvt.traj"
):
    """
    Langevin サーモスタットを用いた NVT 計算。

    Args:
        temp_k: 設定温度 (K)
        dt_fs: タイムステップ (fs). 水素含む系は 0.5-1.0, それ以外は 1.0-2.0 推奨
        friction: 摩擦係数 (ASE 時間単位の逆数). 0.002-0.02 が一般的（公式チュートリアルは 0.002）
                  大きいほど温度制御が強いが、ダイナミクスに粘性が入る
    """
    MaxwellBoltzmannDistribution(atoms, temperature_K=temp_k)
    Stationary(atoms)

    traj = Trajectory(traj_file, 'w', atoms)

    dyn = Langevin(
        atoms,
        timestep=dt_fs * units.fs,
        temperature_K=temp_k,
        friction=friction,
        trajectory=traj,
        loginterval=10
    )

    print(f"Starting NVT-MD: {steps} steps @ {temp_k}K")
    dyn.run(steps)
    print("MD Finished.")
    return atoms
```

#### NVT アンサンブル（NVTBerendsen）

Berendsen サーモスタットによる NVT です。温度の収束が速く平衡化に向きますが、正準分布を厳密には再現しないため、本番サンプリングには Langevin や NVE 切り替えを検討してください。

```python
from ase.md.nvtberendsen import NVTBerendsen

dyn = NVTBerendsen(
    atoms,
    timestep=1.0 * units.fs,
    temperature_K=300,
    # taut: 温度カップリングの時定数。小さいほど強い緩和（公式チュートリアルは 1 fs）。
    # 急速な平衡化には小さく、平衡後のサンプリングでは大きめ（数十-数百 fs）にする
    taut=1.0 * units.fs,
)
dyn.attach(traj.write, interval=10)
dyn.run(5000)
traj.close()
```

#### NPT アンサンブル（NPT / Nose-Hoover）

粒子数 (N)、圧力 (P)、温度 (T) を一定に保ちます。結晶の格子定数を有限温度で決定する場合や、密度を緩和させる場合に使用します。

**警告: 真空を含む系（slab・孤立分子・クラスター）に NPT を使用しないでください。** 圧力制御はセル全体を変形させるため、真空領域ごとセルが潰れます。表面系は NVT を使用してください。表面系で面内格子のみ緩和したい場合は、`mask` でセルの自由度を制限します。

```python
import numpy as np

# z 方向（真空方向）のセルを固定し、面内 (x, y) のみ可変にする例
dyn = NPT(atoms, ..., mask=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]]))
```

**注意**: `ase.md.npt.NPT` はセルが**上三角行列**であることを要求します。`bulk("Cu")` のようなプリミティブセルをそのまま渡すと実行時に `NotImplementedError` になります（よくあるエラー参照）。また ASE 3.28 以降では `ase.md.melchionna.MelchionnaNPT` への移行が告知されています（`NPT` も FutureWarning 付きで動作します）。

```python
from ase.io import Trajectory
from ase.md.npt import NPT
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

def run_npt_md(
    atoms,
    temp_k: float = 300.0,
    pressure_bar: float = 1.01325,
    steps: int = 5000,
    dt_fs: float = 1.0,
    ttime: float = 25.0,
    ptime: float = 160.0,
    pfactor: float = None,
    traj_file: str = "md_npt.traj"
):
    """
    Nose-Hoover thermostat + Parrinello-Rahman barostat を用いた NPT 計算。

    Args:
        ttime: 温度制御の特性時間 (fs). 通常 25*dt 程度（公式チュートリアルは 20 fs）
        ptime: 圧力制御の特性時間 (fs). 温度制御より十分長くとる（100-200 fs 程度）
        pfactor: 圧力制御の係数 (bulk modulus に関連).
                 None の場合 pfactor = B * ptime^2 で自動設定（B = 75 GPa, 金属想定）。
                 公式チュートリアルの典型値は 2e6 GPa*fs^2
    """
    external_stress = pressure_bar * units.bar

    if pfactor is None:
        bulk_modulus_guess_gpa = 75.0
        pfactor = (bulk_modulus_guess_gpa * units.GPa) * (ptime * units.fs)**2

    MaxwellBoltzmannDistribution(atoms, temperature_K=temp_k)
    Stationary(atoms)

    traj = Trajectory(traj_file, 'w', atoms)

    dyn = NPT(
        atoms,
        timestep=dt_fs * units.fs,
        temperature_K=temp_k,
        externalstress=external_stress,
        ttime=ttime * units.fs,
        pfactor=pfactor,
        trajectory=traj,
        loginterval=10
    )

    print(f"Starting NPT-MD: {steps} steps @ {temp_k}K, {pressure_bar}bar")
    dyn.run(steps)
    return atoms
```

#### NPT アンサンブル（NPTBerendsen）

Berendsen バロスタットによる簡潔な NPT です。`compressibility_au`（等温圧縮率）は材料依存のパラメータで、公式チュートリアルでも明示指定されています。

```python
from ase.md.nptberendsen import NPTBerendsen

dyn = NPTBerendsen(
    atoms,
    timestep=1.0 * units.fs,
    temperature_K=300,
    pressure_au=1.01325 * units.bar,          # 1 atm
    taut=100 * units.fs,                       # 温度カップリング時定数
    taup=1000 * units.fs,                      # 圧力カップリング時定数
    compressibility_au=4.57e-5 / units.bar,    # 等温圧縮率（水の値。材料に応じて設定）
)
dyn.attach(traj.write, interval=10)
dyn.run(10000)
traj.close()
```

### ルート B: matlantis-features MDFeature

matlantis-features の `MDFeature` は、Integrator（アンサンブル指定）と `ASEMDSystem`（MD 系の状態管理）を組み合わせた高レベル MD API です。チェックポイント/リスタートに対応し、Post-MD Feature（輸送特性計算）と直結します。

```python
import numpy as np
from ase import units
from matlantis_features.features.md import (
    MDFeature,
    ASEMDSystem,
    NPTIntegrator,
)
from matlantis_features.utils.calculators import pfp_estimator_fn

estimator_fn = pfp_estimator_fn(model_version="v8.0.0", calc_mode="PBE")

# MD 系の作成と初期温度の付与
mdsys = ASEMDSystem(atoms.copy())
mdsys.init_temperature(300)

# Integrator の選択（ここでは NPT: Nose-Hoover + Parrinello-Rahman）
# timestep / ttime は fs 単位の数値で指定する（公式例準拠）
integrator = NPTIntegrator(
    timestep=1.0,
    temperature=300,
    pressure=101325 * units.Pascal,
    ttime=100,
    pfactor=100 * units.fs,
    mask=np.eye(3),
)

md = MDFeature(
    integrator=integrator,
    n_run=20000,                  # MD ステップ数
    traj_file_name="md.traj",
    traj_freq=100,                # トラジェクトリ保存間隔
    show_progress_bar=True,
    show_logger=True,
    logger_interval=1000,
    estimator_fn=estimator_fn,
)

md_result = md(mdsys)
# md_result.traj_path: トラジェクトリの保存先
# md_result.checkpoint_path: チェックポイントの保存先（checkpoint_freq で間隔指定）
```

#### Integrator の一覧

| Integrator | アンサンブル | 主な引数 |
| :--- | :--- | :--- |
| `VelocityVerletIntegrator` | NVE | `timestep` |
| `LangevinIntegrator` | NVT | `timestep`, `temperature`, `friction` |
| `AndersenIntegrator` | NVT | `timestep`, `temperature`, `collision_frequency` |
| `NVTBerendsenIntegrator` | NVT | `timestep`, `temperature`, `taut` |
| `NPTBerendsenIntegrator` | NPT | `timestep`, `temperature`, `pressure`, `taut`, `taup`, `compressibility` |
| `NPTIntegrator` | NPT | `timestep`, `temperature`, `pressure`, `ttime`, `pfactor`, `mask` |

チェックポイントからのリスタートや拡張機能（`TemperatureScheduler` による温度スケジュール等）の詳細は、公式例 `matlantis/references/examples/MD_adv_usage.md` を参照してください。

## 解析（使用したルートに対応）

### ルート A の解析: .traj ファイルの読み込み

ASE 標準 MD のトラジェクトリは `.traj` ファイルとして保存されます。温度・エネルギーの時系列を確認し、平衡化の完了と異常の有無を判断します。

```python
import matplotlib.pyplot as plt
from ase.io import read

images = read("md_nvt.traj", index=":")

temps = [atoms.get_temperature() for atoms in images]
energies = [atoms.get_potential_energy() for atoms in images]

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(temps)
axes[0].set_xlabel("frame")
axes[0].set_ylabel("temperature (K)")
axes[0].grid(True)
axes[1].plot(energies)
axes[1].set_xlabel("frame")
axes[1].set_ylabel("potential energy (eV)")
axes[1].grid(True)
plt.show()
```

MSD・RDF・輸送特性（拡散係数、粘度、熱伝導率）の算出方法は **物性計算（property-analysis/SKILL.md）** を参照してください。

### ルート B の解析: MDFeatureResult と Post-MD Feature

`MDFeature` の結果は `MDFeatureResult` で返されます。保存済みトラジェクトリから復元することもできます。輸送特性は Post-MD Feature に `MDFeatureResult` を渡して計算します。

```python
import numpy as np
from matlantis_features.features.md import MDFeatureResult, PostMDDiffusionFeature

# 実行直後の md_result をそのまま使うか、保存済み traj から復元する
md_result = MDFeatureResult.from_traj_obj("md.traj")

# 拡散係数の計算（MSD から）
diffusion = PostMDDiffusionFeature()
diffusion_result = diffusion(
    md_result, init_time=1000.0, stride=10,
    number_of_segments=1, direction=np.array([[1, 1, 1]]),
    method="segment", effective_msd_range=(0.1, 0.9),
)
```

| Feature | 計算する物性 |
|---------|------------|
| `PostMDDiffusionFeature` | 拡散係数 (MSD から) |
| `PostMDSpecificHeatFeature` | 比熱 |
| `PostMDThermalExpansionFeature` | 熱膨張係数 |
| `PostEMDViscosityFeature` | 粘度 (平衡 MD, Green-Kubo 法) |
| `PostNEMDViscosityFeature` | 粘度 (非平衡 MD, NEMD 法) |
| `PostNEMDThermalConductivityFeature` | 熱伝導率 (非平衡 MD) |

EMD (Equilibrium MD) と NEMD (Non-Equilibrium MD) の使い分け:
- `PostEMDViscosityFeature`: 平衡 MD から Green-Kubo 関係で計算。小さい系でも適用可能
- `PostNEMDViscosityFeature`: 非平衡 MD（せん断流）で直接計算。大きな系で高精度

各 Feature の詳細な使い方・パラメータは **物性計算（property-analysis/SKILL.md）** を参照してください。

## 発展的な MD 制御

### 温度スケジュール（アニーリング・Melt-Quench）

昇温 → 保持 → 降温のように温度を段階的に変える場合は、温度を変えながら MD を順に実行します。アモルファス構造の作成（Melt-Quench 法）や、準安定構造からの脱出（アニーリング）に使用します。

```python
from ase import units
from ase.md.langevin import Langevin

# 昇温 → 高温保持 → 降温（Melt-Quench の基本形）
temperatures = [300, 600, 1200, 1200, 600, 300]
for i, temp in enumerate(temperatures):
    dyn = Langevin(
        atoms, timestep=1.0 * units.fs,
        temperature_K=temp, friction=0.02,
        trajectory=f"anneal_{i:02d}_{temp}K.traj", loginterval=100,
    )
    dyn.run(5000)  # 各温度 5 ps（急冷速度は段数とステップ数で調整）
```

急冷速度（降温の刻みと各温度の保持時間）はアモルファス構造の品質に直結します。Melt-Quench によるアモルファス作成の指針は modeling/SKILL.md を、ルート B での温度スケジュール（`TemperatureScheduler`）は公式例 `matlantis/references/examples/MD_adv_usage.md` を参照してください。

### 熱分解（パイロリシス）シミュレーション

熱分解のような多数の逐次反応は、NEB で 1 つずつ経路を追うのではなく、**高温 MD（2000〜4000 K 程度に昇温して加速）で反応を直接観測する**のが定石です。結果は昇温速度に依存するため複数条件での比較が必要で、加速用の高温では実験条件と反応機構が異なる可能性にも注意してください。

凝集構造の作成 → NPT/NVT 平衡化 → 昇温条件の比較 → `get_mol_list` によるフラグメント解析（分解生成物のカウント・プロット）までの完全な実装が、公式コミュニティ例 `matlantis/references/community-examples/epoxy_pyrolysis/`（エポキシ樹脂の熱分解）にあります。熱分解シミュレーションを実装する際はこの例を参照してください。

### メタダイナミクス（PLUMED）

PLUMED を使った Well-Tempered Metadynamics で自由エネルギー面 (FES) を探索します。

**CV (Collective Variable, 反応座標)** は系の状態変化を代表する低次元の指標（原子間距離、結合角、配位数など）です。メタダイナミクスではこの CV 空間にガウス型バイアスを蓄積して遷移を促進します。

#### PLUMED 入力ファイルの例

```text
UNITS LENGTH=A ENERGY=eV
dist: DISTANCE ATOMS=1,2
uwall: UPPER_WALLS ARG=dist AT=3.5 KAPPA=15.0 EXP=2 EPS=1 OFFSET=0
METAD ARG=dist SIGMA=0.2 HEIGHT=0.05 PACE=200 BIASFACTOR=10.0 TEMP=300.0 LABEL=metad FILE=output/HILLS
```

#### パラメータ設定ガイドライン

| パラメータ | 初期値の目安 | 説明 |
|-----------|------------|------|
| `HEIGHT` | kT 程度 (300K で約 0.025 eV) | ガウス高さ。高すぎると不安定化 |
| `SIGMA` | CV 可動範囲の 5-10% | ガウス幅。距離 CV が 1-5A なら約 0.2A |
| `BIASFACTOR` | 10-15 から開始 | Well-Tempered 係数。大きいほど探索が進むが収束判定が困難 |
| `UPPER_WALLS` / `LOWER_WALLS` | 物理的に意味のある範囲 | CV の探索範囲制限。結合解離系で特に有効 |

#### 実行環境のセットアップ

Matlantis 上での PLUMED の実行（PFP との接続）は LAMMPS 経由で行います（`fix plumed`）。PLUMED のインストール、環境変数の設定、LAMMPS との連携、`sum_hills` による FES の後処理まで含めた完全な手順は、公式コミュニティ例 `matlantis/references/community-examples/plumed_well_tempered_metadynamics/` を参照してください。

### スケジューラ・外場・仮想壁（pfcc-extras）

堆積（`DepositionScheduler`）、分子削除（`DeleteMoleculeScheduler`）、温度スケジュール（`TemperatureScaleScheduler`）、外部電場（`ApplyUniformEfield`）、境界制御（`ElasticVirtualWall`）などの高度な MD 制御は、**pfcc_extras ユーティリティ（pfcc-extras/SKILL.md）** が実装パターンをカバーしています。これらが必要な場合は ASE の低レベル API を直接操作せず、pfcc-extras のスケジューラ・制約 API を優先してください。

## ベストプラクティス

### パラメータ設定の指針

| パラメータ | 推奨範囲 | 備考 |
|-----------|---------|------|
| タイムステップ (dt) | 0.5-2.0 fs | 水素含む系: 0.5-1.0 fs, 重原子のみ: 1.0-2.0 fs |
| Friction (NVT) | 0.002-0.02 | 弱い: 拡散係数の正確な評価向き。強い: 温度安定性重視 |
| pfactor (NPT) | bulk modulus 依存 | 柔らかい材料: 小さめ、硬い材料: 大きめ |
| 保存間隔 | 10-100 ステップ | 解析精度と保存容量のバランス |

**水素を含む系でタイムステップを大きくするのは非推奨です。** 水素 (H) は最も軽い元素で振動周期が短い（X-H 伸縮振動は周期 ~10 fs）ため、タイムステップが大きいと 1 ステップあたりの変位が振動を追跡できないほど大きくなり、水素が飛んで行ったり他の原子にめり込んだりする物理的に不適当な挙動を示すことがあります。水素を含む系では 1.0 fs 以下（安全には 0.5 fs）を使用してください。

### シミュレーション長と系サイズの目安

必要なシミュレーション長は解析対象の物性で決まります。例示コードのステップ数（数千〜数万）は動作確認用であり、本番計算では以下を目安に設計してください。

| 目的 | シミュレーション長の目安 | 備考 |
| :--- | :--- | :--- |
| 平衡化 | 10-100 ps | 温度・エネルギー（NPT では体積）のプラトー到達で判断 |
| RDF・構造解析 | 数十 ps | 平衡化後のデータで十分収束しやすい |
| 拡散係数（液体） | 数百 ps〜数 ns | MSD の線形領域が十分に取れる長さが必要 |
| 粘度（Green-Kubo） | 数 ns | 相関関数の収束が遅い。複数の独立ランの平均も有効 |

系サイズは数百〜数千原子を目安にしてください。セルが小さすぎると周期イメージ間の人工的な相互作用（有限サイズ効果）が物性値に影響します。拡散係数や粘度は系サイズ依存性が知られているため、可能ならサイズを変えて収束を確認してください。

### 温度のゆらぎは異常ではない

カノニカルアンサンブルでは温度の統計ゆらぎの標準偏差は σ_T ≈ T·√(2/(3N)) で、系が小さいほど大きくなります（例: N=100 の系では 300 K で ±25 K 程度が正常）。瞬時温度の振れを「温度制御の失敗」と誤認して摩擦係数を強めたりせず、**平均温度**が設定値に一致しているかで判断してください。

### 平衡化が困難な系への対処

ガラス・アモルファス、長鎖ポリマー、低温の固体内拡散・合金の規則化、核生成などは、構造緩和時間が MD の到達可能時間（~ns）を超えるため、計算時間を延ばすだけでは平衡化できません（化学反応と同じレアイベント問題が構造緩和に現れたものです）。以下の対処を優先順に検討してください。

1. **平衡に近い初期構造を作ってから MD を始める**（最重要・最安価）: MD に平衡化を任せず、液体は `LiquidGenerator`、長鎖ポリマーは EMC でほぼ平衡な凝集構造を生成する（modeling/SKILL.md）
2. **高温平衡化 → 段階降温**: 緩和時間は温度に対して指数的に短くなる。高温で平衡化してから温度スケジュール（発展的な MD 制御参照）で目標温度へ段階的に冷却する。降温が速すぎると高温構造が凍結される
3. **独立レプリカによる検証**: 初期構造・乱数シードを変えた複数の独立ランで物性値が同じ値に収束するか確認する。収束しなければ未平衡であり、結果に初期条件依存のバイアスがある
4. **拡張サンプリング法**: メタダイナミクス（本ガイドの PLUMED 参照）で遅い自由度にバイアスをかける。レプリカ交換 MD（REMD）は ASE 標準にはなく、Matlantis-LAMMPS（`temper` コマンド）の利用が現実的
5. **モンテカルロの併用**: 合金の規則化・偏析など原子の入れ替わりが遅い系では、拡散を MD で待たず原子スワップ MC（LAMMPS の `fix atom/swap` 等）で配置の平衡化を先に済ませる
6. **Light-PFP で到達時間を延ばす**: 特定材料向けに学習した Light-PFP カスタムモデル（PFP の 15-40 倍高速）で 1-2 桁長いシミュレーションが可能。緩和時間が「あと 1 桁」の系に有効
7. **高温で測って外挿する**: 拡散係数などは複数の高温点で計算し Arrhenius 則で目標温度へ外挿する。機構が温度で変わらないことの確認が前提

### サーモスタットの強弱と物性への影響

- **弱いサーモスタット (friction: 0.002)**: 拡散係数、粘度など動的性質の評価に適する。サーモスタットのバイアスが小さい
- **強いサーモスタット (friction: 0.02)**: 運動温度の平衡化を速め、過渡の温度安定性が高い（平衡化向き）。ただし速まるのは熱的緩和のみで、拡散・構造などの遅い緩和は加速せず、強すぎると過減衰でかえって遅くなる
- **NVE に切り替え**: サーモスタットの影響を完全に排除したい場合。平衡化後に NVT -> NVE へ切り替える手法もある

### ZeroRotation は孤立系のみ

`Stationary`（全体並進の除去）はすべての系で適用しますが、`ZeroRotation`（全体回転の除去）は孤立系（分子・クラスター）のみに適用してください。周期系では全体回転が定義されないため不適切です。

### ルート選択の指針

チェックポイント/リスタートが必要な長時間 MD や、輸送特性の算出が目的の場合はルート B（MDFeature）を選択してください。Post-MD Feature は `MDFeatureResult` を入力に取るため、ルート B との組み合わせが最も自然です。

### その他の推奨事項

1. **平衡化 (Equilibration)**: 必ず「平衡化 → 本計算」の 2 段階プロトコル（MD 実行セクション参照）に従い、解析には本計算のデータのみを使用してください。

2. **トラジェクトリの可視化**: 計算終了後、必ずトラジェクトリを可視化して、原子が不自然に飛び出したりクラスタ化したりしていないか確認してください（visualization/SKILL.md 参照）。

3. **チェックポイント/リスタート**: 長時間計算では途中経過を保存し再開できるようにしてください。ルート A では `.traj` から `read('last.traj')` で再開、ルート B では checkpoint 機能を使用します。

4. **長時間 MD**: 長時間 MD はバックグラウンドジョブとして実行することを推奨します（background-job/SKILL.md 参照）。

5. **トークン消費の見積もり**: 本番計算の前に短いテストラン (100-1000 ステップ) でトークン消費量を確認してください。

## よくあるエラーと対処

| エラー・問題 | 原因 | 対処 |
|-------------|------|------|
| Simulation Explosion (原子が飛び散る) | 初期構造の原子間距離が近すぎる、温度上昇が急激、タイムステップが大きすぎる（特に水素を含む系） | 事前に構造最適化を行い、タイムステップを小さくし (0.1 fs)、低温 (10K) から徐々に昇温してください |
| `Too many neighbors` | NPT でセルが潰れて原子密度が極端に高くなった | `pfactor` を調整するか、初期密度が妥当か確認してください |
| `NotImplementedError: Can (so far) only operate on lists of atoms where the computational box is a triangular matrix` | `ase.md.npt.NPT` は上三角セル必須。プリミティブセル等の非上三角セルを渡した | セルを上三角形式に変換する（例: `ase.build.bulk(..., cubic=True)` を使う、または cell を標準化）。`matlantis-contrib` の NPT example も参照 |
| 温度が設定値より高い | `FixAtoms` 拘束使用時の自由度計算ずれ | 自由度補正を行うか、拘束なしで計算してください |
| 温度が大きく揺れる | 小さい系の正常な統計ゆらぎ（σ_T ≈ T·√(2/(3N))） | 異常ではない。平均温度で判断する。ゆらぎ自体を小さくしたい場合は系を大きくする |
| NVE の温度が設定値の半分程度になる | 最適化済み構造への速度付与後、運動エネルギーの約半分がポテンシャルへ流れる（エネルギー等分配） | 2 倍の温度で初期化するか、NVT で平衡化してから NVE に切り替える |
| NPT でセルが真空ごと潰れる | slab・分子など真空を含む系に NPT を適用した | NVT を使用する。面内のみ緩和したい場合は `mask` でセル自由度を制限する |
| エネルギー保存が悪い (NVE) | タイムステップが大きすぎる | タイムステップを小さくしてください (0.5 fs 以下) |
| `FutureWarning: NPT thermostat has been moved/renamed` | ASE 3.28 以降で `NPT` が `MelchionnaNPT` に改名された | 警告のみで動作は継続。ASE 3.28 以降では `ase.md.melchionna.MelchionnaNPT` への移行を検討 |
| 化学反応が起こらない | 活性化エネルギーの大きな化学反応はレアイベントであるため、MD のタイムスケール（~ns）では観測困難 | NEB などの反応経路探索手法（reaction/SKILL.md）、加速 MD、メタダイナミクス（本ガイドの PLUMED 参照）など、化学反応の計算に適した手法を使用する。熱分解のような逐次反応は高温 MD で加速して直接観測する（本ガイドの熱分解参照） |
| 計算が途中で止まる | サーバータイムアウトやリソース不足 | `max_retries` を増やし、チェックポイントからリスタートしてください |

## 関連ガイド

- **Calculator 初期化** (setup/SKILL.md): Calculator 初期化
- **構造モデリング** (modeling/SKILL.md): MD 用初期構造の作成
- **構造最適化** (optimization/SKILL.md): MD 前の必須工程
- **物性計算** (property-analysis/SKILL.md): MSD/RDF・拡散係数・粘度・熱伝導率の算出（MD 後処理）
- **反応経路探索** (reaction/SKILL.md): MD で観測困難な化学反応（レアイベント）の NEB / 遷移状態探索
- **バックグラウンドジョブ** (background-job/SKILL.md): 長時間 MD のバックグラウンド実行
- **pfcc_extras ユーティリティ** (pfcc-extras/SKILL.md): スケジューラ・外場・仮想壁による高度な MD 制御
- **可視化** (visualization/SKILL.md): トラジェクトリの確認
