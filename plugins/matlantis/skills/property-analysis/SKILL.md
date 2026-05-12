---
name: mt-property-analysis
description: >
  物性計算とMD後処理を扱うスキルです。
  ElasticTensorFeature, PostElasticPropertiesFeature, 弾性定数, elastic, Young率, バルク弾性率,
  ForceConstantFeature, PostPhononBandFeature, PostPhononDOSFeature, PostPhononThermochemistryFeature,
  フォノン, phonon, 力定数, QHA, 準調和近似,
  VibrationFeature, 振動解析,
  PostMDDiffusionFeature, 拡散係数, MSD, 平均二乗変位,
  PostEMDViscosityFeature, PostNEMDViscosityFeature, 粘度, Green-Kubo,
  PostNEMDThermalConductivityFeature, 熱伝導率,
  PostMDSpecificHeatFeature, 比熱, PostMDThermalExpansionFeature, 熱膨張,
  GasStandardFormationEnthalpyFeature, PostVibrationGasThermoFeature, ガス熱化学,
  RDF, 動径分布関数, 振動状態密度
  に関するコード生成時に使用してください。
---

# 物性計算

## 概要

Matlantis を用いた物性計算と MD 後処理のガイドです。PFP（Preferred Potential）の高速な原子間ポテンシャル計算を活用して、弾性定数、フォノン、準調和近似（QHA）、拡散係数、粘度、熱伝導率、動径分布関数（RDF）、振動状態密度など、幅広い物性を評価できます。

matlantis-features ライブラリは、これらの物性計算を高レベル API として提供しています。各 Feature は callable なオブジェクトであり、`result = feature(atoms)` のパターンで統一的に利用できます。多くの物性計算は 2 段階パイプライン（計算 Feature + 後処理 PostFeature）で構成されています。

本ガイドでは、理論的な背景と実装を分離し、MD trajectory の生成コードと後処理コードを明確に区別して説明します。

## ワークフロー

物性計算の標準的なワークフローは以下の通りです。

1. **構造の準備**: 対象系の構造を作成し、構造最適化を完了する
2. **計算手法の選択**: 目的の物性に応じた Feature クラスを選択する
3. **一次計算の実行**: ForceConstantFeature、ElasticTensorFeature、MDFeature 等を実行する
4. **後処理の実行**: Post-Feature で結果を解析・変換する
5. **結果の解析**: 数値の妥当性を検証し、可視化する

### Feature 選択の判断基準

| 計算対象 | 一次 Feature | 後処理 Feature |
| :--- | :--- | :--- |
| 弾性定数 | `ElasticTensorFeature` | `PostElasticPropertiesFeature` |
| フォノンバンド | `ForceConstantFeature` | `PostPhononBandFeature` |
| フォノン DOS | `ForceConstantFeature` | `PostPhononDOSFeature` |
| 熱力学量 | `ForceConstantFeature` | `PostPhononThermochemistryFeature` |
| QHA | `ForceConstantFeature` | `PostPhononQHAGibbsFeature` |
| 拡散係数 | `MDFeature` | `PostMDDiffusionFeature` |
| 比熱 | `MDFeature` | `PostMDSpecificHeatFeature` |
| 熱膨張係数 | `MDFeature` | `PostMDThermalExpansionFeature` |
| 粘度（EMD） | `MDFeature` | `PostEMDViscosityFeature` |
| 粘度（NEMD） | `MDFeature` | `PostNEMDViscosityFeature` |
| 熱伝導率 | `MDFeature` | `PostNEMDThermalConductivityFeature` |
| 振動解析 | `VibrationFeature` | --- |
| 気体熱力学 | `VibrationFeature` | `PostVibrationGasThermoFeature` |
| 標準生成エンタルピー | --- | `GasStandardFormationEnthalpyFeature` |

## 実装パターン

### パターン A: 弾性定数

最適化済み結晶に小ひずみを付加し、応力応答から弾性テンソルを算出します。Young 率、剪断弾性率、体積弾性率、Poisson 比を得ることができます。

#### 理論

弾性定数は、結晶に微小ひずみを付加した際の応力応答から求めます。

1. 最適化済みの結晶構造を用意する
2. 6 つの独立な小ひずみ（Voigt 表記）を付加する
3. 各ひずみ下で原子位置のみ再緩和する
4. `get_stress()` から弾性定数テンソルを計算する

#### 実装

```python
from matlantis_features.features.elasticity import (
    ElasticTensorFeature,
    PostElasticPropertiesFeature,
)
from matlantis_features.utils.calculators import pfp_estimator_fn

estimator_fn = pfp_estimator_fn(model_version="v8.0.0", calc_mode="PBE")

# Step 1: 弾性テンソルの計算
elastic = ElasticTensorFeature(estimator_fn=estimator_fn)
elastic_result = elastic(optimized_atoms)

# Step 2: 物性値の算出
props = PostElasticPropertiesFeature()
props_result = props(elastic_result)

# 結果の取得
# Young's modulus, shear modulus, bulk modulus, Poisson's ratio
print(f"Young's modulus: {props_result.youngs_modulus} GPa")
print(f"Shear modulus: {props_result.shear_modulus} GPa")
print(f"Bulk modulus: {props_result.bulk_modulus} GPa")
print(f"Poisson's ratio: {props_result.poissons_ratio}")
```

### パターン B: フォノン計算

フォノンバンド構造、状態密度（DOS）、および熱力学量を計算します。力定数の計算を起点として、段階的に各種物性を算出します。

#### 理論

フォノン計算は有限変位法（Finite Displacement Method）に基づきます。

1. primitive cell を使用する（計算効率とバンド構造の正確性のため）
2. 十分なサイズのスーパーセルを構築する
3. 各原子を微小変位させてフォルダー定数を計算する
4. フォノンバンド、DOS、熱力学量を段階的に算出する

#### 実装

```python
from matlantis_features.features.phonon import (
    ForceConstantFeature,
    PostPhononBandFeature,
    PostPhononDOSFeature,
    PostPhononThermochemistryFeature,
)

estimator_fn = pfp_estimator_fn(model_version="v8.0.0", calc_mode="PBE")

# Step 1: 力定数の計算
# supercell: スーパーセルのサイズ（大きいほど正確だが遅い。4-6で開始推奨）
# delta: 変位量 (Angstrom)
fc = ForceConstantFeature(
    supercell=6,
    delta=0.01,
    estimator_fn=estimator_fn,
)
fc_result = fc(primitive_atoms)

# Step 2a: バンド構造
band = PostPhononBandFeature()
band_result = band(fc_result)

# Step 2b: 状態密度
dos = PostPhononDOSFeature(mesh=[20, 20, 20])
dos_result = dos(fc_result)

# Step 2c: 熱力学量（自由エネルギー、エントロピー、比熱など）
thermo = PostPhononThermochemistryFeature()
thermo_result = thermo(fc_result)
```

### パターン C: 準調和近似（QHA）

QHA は、異なる体積でのフォノン計算結果を統合し、有限温度での熱力学量を精密に評価します。

```python
from matlantis_features.features.phonon import (
    ForceConstantFeature,
    PostPhononQHAGibbsFeature,
)

# 複数体積でのForceConstant計算結果を用いてQHA Gibbsエネルギーを算出
fc = ForceConstantFeature(supercell=6, delta=0.01, estimator_fn=estimator_fn)

# 体積を変えた複数構造でfc_resultを取得
# ...

qha = PostPhononQHAGibbsFeature()
qha_result = qha(fc_results)
```

### パターン D: 拡散係数

MD trajectory から平均二乗変位（MSD）を計算し、拡散係数を算出します。

#### 理論

1. 十分な長さの MD trajectory を生成する
2. MSD を計算する（周期境界を考慮した unwrapping が必要）
3. MSD の線形領域を特定してフィッティングする
4. Einstein の関係式から拡散係数を求める

#### 実装

```python
from matlantis_features.features.md import (
    MDFeature,
    PostMDDiffusionFeature,
)

estimator_fn = pfp_estimator_fn(model_version="v8.0.0", calc_mode="PBE")

# Step 1: MD trajectoryの生成
md = MDFeature(
    estimator_fn=estimator_fn,
    # MD パラメータ（温度、ステップ数、アンサンブル等）
)
md_result = md(atoms)

# Step 2: 拡散係数の計算
diffusion = PostMDDiffusionFeature()
diffusion_result = diffusion(md_result)
```

#### MSD 解析の注意点

- **Unwrapping**: 周期境界条件下では原子がセル境界を越えると座標が折り返されるため、MSD 計算前に座標を unwrap する必要があります
- **線形領域の特定**: MSD が時間に対して線形になる領域（バリスティック領域の後、統計誤差が大きくなる前）を選択する必要があります
- **統計精度**: 複数の時間起点から MSD を平均化し、統計誤差を低減させます

### パターン E: 粘度

平衡 MD（EMD）の Green-Kubo 法、または非平衡 MD（NEMD）による粘度計算です。

#### 理論（Green-Kubo 法）

- 平衡 MD で応力テンソルの自己相関関数を計算する
- Green-Kubo 積分により粘度を求める
- 統計誤差の評価と平均化が重要

#### 実装

```python
from matlantis_features.features.md import (
    MDFeature,
    PostEMDViscosityFeature,
    PostNEMDViscosityFeature,
)

# EMD（Green-Kubo法）
md_result = md(atoms)
emd_viscosity = PostEMDViscosityFeature()
emd_result = emd_viscosity(md_result)

# NEMD
nemd_viscosity = PostNEMDViscosityFeature()
nemd_result = nemd_viscosity(md_result)
```

### パターン F: 熱伝導率

```python
from matlantis_features.features.md import (
    MDFeature,
    PostNEMDThermalConductivityFeature,
)

md_result = md(atoms)
thermal = PostNEMDThermalConductivityFeature()
thermal_result = thermal(md_result)
```

統計誤差や平均化の必要性に注意してください。単一の MD run からの結果は統計的に不十分な場合があります。

### パターン G: 比熱・熱膨張係数

```python
from matlantis_features.features.md import (
    PostMDSpecificHeatFeature,
    PostMDThermalExpansionFeature,
)

# 比熱
specific_heat = PostMDSpecificHeatFeature()
sh_result = specific_heat(md_result)

# 熱膨張係数
thermal_expansion = PostMDThermalExpansionFeature()
te_result = thermal_expansion(md_result)
```

### パターン H: 動径分布関数（RDF）

MD trajectory から RDF を計算し、構造の特徴を解析します。

```python
# MD trajectoryからRDFを計算
# ase.geometry.analysis または matlantis-features の後処理機能を使用
from ase.geometry.analysis import Analysis

analysis = Analysis(trajectory)
rdf = analysis.get_rdf(rmax=10.0, nbins=200)
```

### パターン I: 振動解析と気体熱力学

分子の振動モードと、気体状態の熱力学量を計算します。

```python
from matlantis_features.features.common.vibration import VibrationFeature
from matlantis_features.features.common.gas_thermo import PostVibrationGasThermoFeature
from matlantis_features.features.common.gas_formation_enthalpy import (
    GasStandardFormationEnthalpyFeature,
)

estimator_fn = pfp_estimator_fn(model_version="v8.0.0", calc_mode="PBE")

# 振動解析
vib = VibrationFeature(delta=0.01, estimator_fn=estimator_fn)
vib_result = vib(molecule_atoms)

# 気体熱力学量（理想気体近似）
gas_thermo = PostVibrationGasThermoFeature()
gas_result = gas_thermo(vib_result)

# 標準生成エンタルピー
formation = GasStandardFormationEnthalpyFeature(estimator_fn=estimator_fn)
formation_result = formation(molecule_atoms)
```

### パターン J: 振動状態密度（VDOS from MD）

MD trajectory から振動状態密度を計算します。速度自己相関関数のフーリエ変換に基づきます。

```python
# MD trajectoryから速度自己相関関数を計算し、
# フーリエ変換で振動状態密度を得る
# 実装は解析ライブラリに依存
```

### パターン K: Dipole 計算（gRPC）

双極子モーメントの計算には gRPC 経由のインターフェースを使用します。

```python
# Dipole feature は gRPC 経由で利用可能
# 詳細は matlantis-features の API リファレンスを参照
```

## ベストプラクティス

### 構造最適化の徹底

物性計算の入力構造は必ず事前に最適化してください。結晶系では `filter=True` または `UnitCellASEFilter` を使ってセル形状も最適化する必要があります。

### Primitive Cell の使用（フォノン計算）

フォノン計算では primitive cell を使用してください。Conventional cell を使用すると、不要なバンドフォールディングが発生し、結果の解釈が困難になります。

### スーパーセルサイズの選択（フォノン計算）

`ForceConstantFeature(supercell=N)` の N は大きいほど正確ですが計算コストが増大します。N=4-6 から開始し、結果が収束しない場合に増やしてください。

### 理論と実装の分離

物性計算のコードを書く際は、理論的な説明（なぜこの手法を使うか）と実装（具体的なコード）を明確に分離してください。

### MD 生成と後処理の分離

MD trajectory の生成コードと後処理（MSD 計算、Green-Kubo 積分等）のコードは分離してください。これにより、同じ trajectory に対して異なる解析を適用したり、trajectory の生成条件を変更したりすることが容易になります。

### 統計誤差への配慮

粘度、熱伝導率などの輸送特性の計算では、統計誤差が大きくなりやすいです。以下の点に注意してください。

- 十分な長さの MD trajectory を生成する
- 複数の独立した MD run から結果を平均化する
- 相関時間を考慮したブロック平均を使用する

### estimator_fn の使用

matlantis-features の Feature クラスには、Estimator インスタンスではなく `estimator_fn`（ファクトリ関数）を渡してください。Feature 内部で必要に応じて新しい Estimator が生成され、インスタンス共有による問題を防ぎます。

## よくあるエラーと対処

| エラー・症状 | 原因 | 対処法 |
| :--- | :--- | :--- |
| フォノンバンドに虚の周波数（負の値）が出る | 構造が十分に最適化されていない | より厳密な fmax（例: 0.001 eV/A）で再最適化する |
| フォノン計算が非常に遅い | supercell が大きすぎる | supercell を小さくして結果の収束を確認する |
| MSD が線形にならない | MD が平衡化していない、または trajectory が短すぎる | 平衡化フェーズを長くする、trajectory を延長する |
| 粘度の値が不安定 | 統計サンプル不足 | MD run を延長する、複数 run を平均化する |
| `MultiCalculatorUseDetected` | Estimator を共有した | `estimator_fn` パターンを使用する |
| 弾性定数が対称でない | 構造が十分に最適化されていない | セル最適化（filter=True）を含む最適化を実行する |
| Conventional cell で不要なバンドが出る | primitive cell を使っていない | `spglib` 等で primitive cell に変換する |

## 関連ガイド

- 構造最適化: 物性計算の前段階として必須
- **反応経路探索** (reaction/SKILL.md): 振動解析を遷移状態検証に使用
- **可視化** (visualization/SKILL.md): 計算結果の構造確認
- **Group Drive ストレージ** (storage/SKILL.md): 大規模計算結果の保存
