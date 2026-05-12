---
name: mt-dynamics
description: >
  分子動力学(MD)シミュレーションを扱うスキルです。
  Langevin, VelocityVerlet, NPTBerendsen, NVE, NVT, NPT,
  MaxwellBoltzmannDistribution, timestep, friction, temperature_K,
  MDFeature, PostMDDiffusionFeature, PostMDSpecificHeatFeature,
  PostEMDViscosityFeature, PostNEMDViscosityFeature, PostNEMDThermalConductivityFeature,
  PLUMED, メタダイナミクス, metadynamics, FES, collective variable,
  DepositionScheduler, ApplyUniformEfield, ElasticVirtualWall,
  分子動力学, MD, 拡散係数, 粘度, 熱伝導率
  に関するコード生成時に使用してください。
---
# 分子動力学

## 概要

有限温度における原子・分子の時間発展をシミュレーションするガイドです。NVE / NVT / NPT アンサンブルの基本実装、PLUMED メタダイナミクス、matlantis-features の高レベル MD API、トラジェクトリの後処理による輸送特性計算、およびスケジューラを使った高度な MD 制御を包括的に扱います。

分子動力学 (MD) は拡散係数、粘性、熱伝導率の計算、アモルファス構造の作成（Melt-Quench 法）、相転移のシミュレーションなどに使用されます。MD の各ステップでトークンを消費するため、本番計算の前に短いテストランで消費量を見積もることを推奨します。

長時間の MD Notebook 実行では、foreground 実行だけに決め打ちせず background 実行も常に選択肢として提示してください。数時間以上の計算や離席を伴う場合は `background-job/SKILL.md` に従い `mtl-bg-job` を優先します。

## ワークフロー

```
1. Prepare     - 構造最適化済みの atoms を用意し、Calculator をセット
2. Initialize  - 初期温度に応じた速度分布 (MaxwellBoltzmann) を付与
3. Choose      - 目的に応じたアンサンブル (NVE / NVT / NPT) を選択
4. Configure   - タイムステップ、摩擦係数、圧力パラメータ等を設定
5. Run         - 時間発展計算を実行、トラジェクトリを保存
6. Analyze     - 保存された .traj ファイルの解析（MSD, RDF, 輸送特性等）
```

## 実装パターン

### パターン A: 共通セットアップ

すべての MD 計算で共通する初期速度の付与とトラジェクトリ保存の設定です。

```python
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.io.trajectory import Trajectory
from ase import units

# Calculator を設定
atoms.calc = calculator

# 初期速度の付与
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
Stationary(atoms)       # 重心運動（全体の並進）を除去
ZeroRotation(atoms)     # 全体の回転を除去

# トラジェクトリ保存
traj = Trajectory("md.traj", "w", atoms)
```

`Stationary` と `ZeroRotation` をセットで適用することで、全体の並進・回転の寄生成分を確実に除去します。

### パターン B: NVE アンサンブル (VelocityVerlet)

体積 (V) とエネルギー (E) を一定に保つ、最も基本的なアンサンブルです。サーモスタットのバイアスがないため、動的性質（拡散係数など）の正確な評価に適しています。

```python
from ase.md.verlet import VelocityVerlet

dyn = VelocityVerlet(atoms, timestep=1.0 * units.fs)
dyn.attach(traj.write, interval=10)  # 10 ステップごとに保存
dyn.run(1000)
traj.close()
```

NVE ではエネルギー保存が正しいかを確認してください。保存が悪い場合はタイムステップを小さくします。

### パターン C: NVT アンサンブル (Langevin)

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
    friction: float = 0.01,
    traj_file: str = "md_nvt.traj"
):
    """
    Langevin サーモスタットを用いた NVT 計算。

    Args:
        temp_k: 設定温度 (K)
        dt_fs: タイムステップ (fs). 水素含む系は 0.5-1.0, それ以外は 1.0-2.0 推奨
        friction: 摩擦係数 (atomic units). 0.002-0.02 が一般的
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

簡潔な最小構成:

```python
from ase.md.langevin import Langevin

dyn = Langevin(atoms, timestep=1.0 * units.fs, temperature_K=300, friction=0.01)
dyn.attach(traj.write, interval=10)
dyn.run(5000)
traj.close()
```

### パターン D: NPT アンサンブル (Nose-Hoover)

粒子数 (N)、圧力 (P)、温度 (T) を一定に保ちます。結晶の格子定数を有限温度で決定する場合や、密度を緩和させる場合に使用します。

```python
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
    pfactor: float = None,
    traj_file: str = "md_npt.traj"
):
    """
    Nose-Hoover barostat を用いた NPT 計算。

    Args:
        ttime: 温度制御の特性時間 (fs). 通常 25*dt 程度
        pfactor: 圧力制御の係数 (bulk modulus に関連).
                 None の場合、75 GPa (金属想定) で自動設定
    """
    external_stress = pressure_bar * units.bar

    if pfactor is None:
        bulk_modulus_guess_gpa = 75.0
        pfactor = (bulk_modulus_guess_gpa * units.GPa) * (ttime * units.fs)**2

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

#### NPTBerendsen（簡潔パターン）

```python
from ase.md.nptberendsen import NPTBerendsen

dyn = NPTBerendsen(
    atoms,
    timestep=1.0 * units.fs,
    temperature_K=300,
    pressure_au=0.0,
    taut=100 * units.fs,
    taup=1000 * units.fs,
)
dyn.attach(traj.write, interval=10)
dyn.run(10000)
traj.close()
```

### パターン E: PLUMED メタダイナミクス

PLUMED を使った Well-Tempered Metadynamics で自由エネルギー面 (FES) を探索します。

**CV (Collective Variable, 反応座標)** は系の状態変化を代表する低次元の指標（原子間距離、結合角、配位数など）です。メタダイナミクスではこの CV 空間にガウス型バイアスを蓄積して遷移を促進します。

#### PLUMED 入力ファイルの例

```text
UNITS LENGTH=A ENERGY=eV
dist: DISTANCE ATOMS=1,2
uwall: UPPER_WALLS ARG=dist AT=3.5 KAPPA=15.0 EXP=2 EPS=1 OFFSET=0
METAD ARG=dist SIGMA=0.2 HEIGHT=0.05 PACE=200 BIASFACTOR=80.0 TEMP=300.0 LABEL=metad FILE=output/HILLS
```

#### パラメータ設定ガイドライン

| パラメータ | 初期値の目安 | 説明 |
|-----------|------------|------|
| `HEIGHT` | kT 程度 (300K で約 0.025 eV) | ガウス高さ。高すぎると不安定化 |
| `SIGMA` | CV 可動範囲の 5-10% | ガウス幅。距離 CV が 1-5A なら約 0.2A |
| `BIASFACTOR` | 10-15 から開始 | Well-Tempered 係数。大きいほど探索が進むが収束判定が困難 |
| `UPPER_WALLS` / `LOWER_WALLS` | 物理的に意味のある範囲 | CV の探索範囲制限。結合解離系で特に有効 |

### パターン F: matlantis-features による MD

matlantis-features の `MDFeature` は高レベルな MD API を提供します。

```python
from matlantis_features.features.md import MDFeature

md = MDFeature(
    estimator_fn=estimator_fn,
    # ... MD パラメータ
)
result = md(atoms)
```

#### MD 後処理 Feature

MD トラジェクトリから輸送特性を計算する後処理 Feature が用意されています。

| Feature | 計算する物性 |
|---------|------------|
| `PostMDDiffusionFeature` | 拡散係数 (MSD から) |
| `PostMDSpecificHeatFeature` | 比熱 |
| `PostMDThermalExpansionFeature` | 熱膨張係数 |
| `PostEMDViscosityFeature` | 粘度 (平衡 MD, Green-Kubo 法) |
| `PostNEMDViscosityFeature` | 粘度 (非平衡 MD, NEMD 法) |
| `PostNEMDThermalConductivityFeature` | 熱伝導率 (非平衡 MD) |

```python
from matlantis_features.features.md import (
    PostMDDiffusionFeature,
    PostMDSpecificHeatFeature,
    PostMDThermalExpansionFeature,
    PostEMDViscosityFeature,
    PostNEMDViscosityFeature,
    PostNEMDThermalConductivityFeature,
)
```

EMD (Equilibrium MD) と NEMD (Non-Equilibrium MD) の使い分け:
- `PostEMDViscosityFeature`: 平衡 MD から Green-Kubo 関係で計算。小さい系でも適用可能
- `PostNEMDViscosityFeature`: 非平衡 MD（せん断流）で直接計算。大きな系で高精度

### パターン G: スケジューラによる高度な MD 制御

`pfcc_extras` のスケジューラを MD に組み込むことで、堆積・削除・温調を自動制御します。

```python
# 1 つの dyn に複数スケジューラを attach
# DepositionScheduler: 分子の逐次堆積
# DeleteMoleculeScheduler: 特定条件で分子を削除
# TemperatureScaleScheduler: 温度の段階的変更
```

#### 入射エネルギーから速度への変換

堆積計算では入射エネルギー (eV) から速度ベクトルへの変換が条件比較に便利です。

#### ElasticVirtualWall による境界制御

蒸発や飛散が支配的な系で計算領域を制御します。

```python
# ElasticVirtualWall: axis, wall_direction (top/bottom), wall_position で指定
```

削除境界は固定値より `cell_z - margin` 形式で定義し、セルサイズ変更時の設定破綻を防ぎます。

### パターン H: 時間依存外場 (AC 電場)

一様電場をステップごとに更新して AC 場応答を追跡します。

```python
# E(t) = E0 * cos(2 * pi * f * t)
# 各ステップで電場の値を再評価して拘束を更新
```

静的電場ではなく時間依存の AC 場応答を得る場合に使用します。

## ベストプラクティス

### pfcc-extras のスケジューラ・制約を優先する

MD シミュレーションにおいて、堆積（`DepositionScheduler`）、外部電場（`ApplyUniformEfield`）、境界制御（`ElasticVirtualWall`）、温度スケジュール（`TemperatureScaleScheduler`）、分子削除（`DeleteMoleculeScheduler`）などの高度な制御には、ASE の低レベル API を直接操作するのではなく pfcc-extras のスケジューラ・制約 API を優先して使用してください。入射エネルギーから速度への変換には `convert_kinetic_energy_to_velocity` を使用し、条件比較を容易にしてください。

### パラメータ設定の指針

| パラメータ | 推奨範囲 | 備考 |
|-----------|---------|------|
| タイムステップ (dt) | 0.5-2.0 fs | 水素含む系: 0.5-1.0 fs, 重原子のみ: 1.0-2.0 fs |
| Friction (NVT) | 0.002-0.02 | 弱い: 拡散係数の正確な評価向き。強い: 温度安定性重視 |
| pfactor (NPT) | bulk modulus 依存 | 柔らかい材料: 小さめ、硬い材料: 大きめ |
| 保存間隔 | 10-100 ステップ | 解析精度と保存容量のバランス |

### サーモスタットの強弱と物性への影響

- **弱いサーモスタット (friction: 0.002)**: 拡散係数、粘度など動的性質の評価に適する。サーモスタットのバイアスが小さい
- **強いサーモスタット (friction: 0.02)**: 温度安定性が高い。構造探索や平衡化に適する
- **NVE に切り替え**: サーモスタットの影響を完全に排除したい場合。平衡化後に NVT -> NVE へ切り替える手法もある

### その他の推奨事項

1. **平衡化 (Equilibration)**: 計算開始直後は系が安定していないため、最初の数千ステップはデータ解析から除外してください。

2. **トラジェクトリの可視化**: 計算終了後、必ずトラジェクトリを可視化して、原子が不自然に飛び出したりクラスタ化したりしていないか確認してください。

3. **チェックポイント/リスタート**: 長時間計算では途中経過を `.traj` に保存し、`read('last.traj')` から再開できるようにしてください。

4. **長時間 MD**: 長時間 MD はバックグラウンドジョブとして実行することを推奨します。

5. **トークン消費の見積もり**: 本番計算の前に短いテストラン (100-1000 ステップ) でトークン消費量を確認してください。

## よくあるエラーと対処

| エラー・問題 | 原因 | 対処 |
|-------------|------|------|
| Simulation Explosion (原子が飛び散る) | 初期構造の原子間距離が近すぎる、温度上昇が急激 | 事前に構造最適化を行い、タイムステップを小さくし (0.1 fs)、低温 (10K) から徐々に昇温してください |
| `Too many neighbors` | NPT でセルが潰れて原子密度が極端に高くなった | `pfactor` を調整するか、初期密度が妥当か確認してください |
| `NotImplementedError` (convert_atoms_to_upper) | ASE の NPT モジュールのバージョン依存バグ (triclinic cell) | セルを直交化に近い形に変形するか、`matlantis-contrib` の NPT example を参照してください |
| 温度が設定値より高い | `FixAtoms` 拘束使用時の自由度計算ずれ | 自由度補正を行うか、拘束なしで計算してください |
| エネルギー保存が悪い (NVE) | タイムステップが大きすぎる | タイムステップを小さくしてください (0.5 fs 以下) |
| 計算が途中で止まる | サーバータイムアウトやリソース不足 | `max_retries` を増やし、チェックポイントからリスタートしてください |

## 関連ガイド

- **Calculator 初期化** (setup/SKILL.md): Calculator 初期化
- **構造モデリング** (modeling/SKILL.md): MD 用初期構造の作成
- **構造最適化** (optimization/SKILL.md): MD 前の必須工程
