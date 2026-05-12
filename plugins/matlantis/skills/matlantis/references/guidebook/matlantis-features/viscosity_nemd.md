# PostNEMDViscosityFeature ​

## はじめに ​

「粘度」は液体の変形に対する抵抗の強さを表す特徴量です。 粘度は複数の方法で求めることができますが、 [PostNEMDViscosityFeature](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeature.html#matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeature) では、 reverse non-equilibrium molecular dynamics (rNEMD) 法を用いて計算を行います。

## 計算手法 ​

rNEMD法では、 系に人工的な速度勾配を発生させたときの速度プロファイルを測定します。

具体的には、以下の手順を繰り返すことで、系全体での温度を変えることなくシミュレーション領域の中央部と下部の間で速度勾配を発生させます。

  1. シミュレーション対象の領域ををある方向（通常はz方向）に沿ってN個のスラブに分割します。
  2. 一番下のスラブから+x方向に最も大きい運動量をもつ原子を見つけ出します。
  3. その原子と同じ元素で、なおかつ-x方向に最も大きい運動量をもつ原子を中央のスラブから見つけ出します。
  4. この2つの原子の速度を交換します。



この方法により、底のスラブが-x方向、中央のスラブが+x方向の速度をもつような速度勾配を発生させることができます。 交換された運動量は次のようになります。

pexchanged=px,bottom−px,middle

ニュートンの粘性法則によると、運動量束は粘性と速度勾配の積となります。

∑pexchanged2tSxy=η∂vx∂z

ここで、 η は粘性、 t は時間、 Sxy はxy面の面積、vx は速度のx方向成分を表します。

このため、系が定常状態に達した後の速度勾配を測定します。 各スラブにおけるx方向の平均速度が計算されます。 以下の図に示すように、 vx と z の間に明確な線形関係が見いだせます。

∂vx∂z はフィッティングした際の傾きとなります。

この計算を実行するには、MD計算そのものに対する改変が必要です。 これは、 [MDFeature](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.md.md.MDFeature.html#matlantis_features.features.md.md.MDFeature) に [RNEMDExtension](/api/docs/ja/matlantis-features/latest/generated/matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension.html#matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension) を渡すことで実行可能です。 `PostNEMDViscosityFeature` によりrNEMDによるシミュレーション結果を加工し、粘性を取得することができます。

以下は `PostNEMDViscosityFeature` を使用した例です。 `assets/viscosity/n-decane.xyz` を入手するには、 ランチャーのExamplesにある Viscosity rNEMD の “Copy files to the current directory” を選択してください。 これによりファイルが生成されます。

default
    
    
    from matlantis_features.atoms import MatlantisAtoms
    from matlantis_features.features.md import MDFeature, ASEMDSystem
    from matlantis_features.features.md.post_features.rnemd_extension import RNEMDExtension
    from matlantis_features.features.md import NVTBerendsenIntegrator
    from matlantis_features.features.md import PostNEMDViscosityFeature
    
    from ase import units
    
    atoms = MatlantisAtoms.from_file("assets/viscosity/n-decane.xyz")
    system = ASEMDSystem(atoms)
    system.init_temperature(293.15)
    integ = NVTBerendsenIntegrator(timestep=1.0, temperature=293.15)
    
    n_slab = 20
    rnemd_interval = 10
    n_run=1000
    md = MDFeature(integrator=integ, n_run=n_run)
    rnemd_extension = RNEMDExtension(rnemd_type="viscosity", n_slab=n_slab)
    nemd_result = md(system, extensions=[(rnemd_extension, rnemd_interval)])
    
    post_feature = PostNEMDViscosityFeature()
    post_feature(nemd_result)

## Tips ​

  * スラブ数 `n_slab` は偶数である必要があります。
  * `PostNEMDThermalConductivityFeature` での分析を実行するには運動量の交換を記録した情報が必要なため、 `RNEMDExtension` を適用した `MDFeature` の出力トラジェクトリを使用する必要があります。
  * 速度勾配の線形性が見られるかどうか確認してください。 もし線形性が見られない場合、MDシミュレーション時間を長くするか、速度交換の頻度をより高く設定してみてください。



## Reference ​

  * Bordat, P. & Müller-Plathe, F. The shear viscosity of molecular fluids: A calculation by reverse nonequilibrium molecular dynamics. J. Chem. Phys. 116, 3362–3369 (2002).


