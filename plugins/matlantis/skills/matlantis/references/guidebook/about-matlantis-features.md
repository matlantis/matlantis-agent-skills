# matlantis-featuresとは ​

## はじめに ​

`matlantis-features` では、さまざまな物性値を計算するための高度な関数、 `Features` を提供しています。

`matlantis-features` を用いることで、短くシンプルなコードで物性値を計算することができます。 皆さんが計算アルゴリズムの詳細に精通していなくても使用でき、複雑な実装も必要ありません。

Launcherの `Examples` 内で `matlantis-features` の実用例を紹介しているのでご参照ください。

## Features ​

`Features` は、その処理内容からcalculation featuresとpost-calculation featuresの2種類に分類できます。 calculation featuresは内部的にPFP API serverを呼び出し、PFPから出力されたenergy、force、stressなどの値を用いて原子シミュレーションを行います。 `MDFeatures` などがこれに該当します。

一方、post-calculation featuresは計算後の結果を分析するための `Features` であり、内部でのPFP API serverの呼び出しはありません。

名前が `Complex` から始まる `Features` は、いくつかの `Features` を複合させたものです。 たとえば、 `ComplexElasticFeature` は `FireASEOptFeature`、 `ElasticFeature` および `PostElasticPropertiesFeature` を複合させた `Feature` です。

`matlantis-features` はさまざまな物性値やシミュレーション手法をカバーしています。 `matlantis-features` で実装されている `Features` および取得できる物性値は以下の表の通りです。

#### 表: matlantis-featuresに実装されているすべてのfeature ​

Categories| Features| Physical Properties  
---|---|---  
Calculation Features| Post Calculation Features|   
structure optimization| BFGSASEOptFeature|   
FireASEOptFeature| |   
FireLBFGSASEOptFeature| |   
MD related features| MDFeatures| [PostMDDiffusionFeature](./matlantis-features/diffusion.html)  
[PostEMDViscosityFeature](./matlantis-features/viscosity_emd.html)| Viscosity|   
[PostNEMDViscosityFeature](./matlantis-features/viscosity_nemd.html)| Viscosity|   
[PostNEMDThermalConductivityFeature](./matlantis-features/thermal_conduct.html)| Thermal conductivity|   
PostMDSpecificHeatFeature| Specific heat|   
PostMDThermalExpansionFeature| Thermal expansion|   
Elastic related features| [ElasticTensorFeature](./matlantis-features/elastic_tensor.html)| PostElasticPropertiesFeature  
ComplexElasticFeature| | Same as above  
Vibration| [VibrationFeature](./matlantis-features/vibration.html)| [PostVibrationGasThermoFeature](./matlantis-features/vibration_thermo.html)  
ComplexGasFormationEnthalpy| | Gas formation enthalpy  
GasStandardFormationEnthalpy| | Gas standard formation enthalpy  
Phonon related features| [ForceConstantFeature](./matlantis-features/force_constant.html)| [PostPhononBandFeature](./matlantis-features/phonon_band.html)  
[PostPhononDOSFeature](./matlantis-features/phonon_dos.html)| Phonon density of states|   
PostPhononThermochemistryFeature| \- specific heat  
\- entropy  
\- internal energy  
\- helmholtz free energy|   
PostPhononModeFeature| Atomic vibration in a phonon mode|   
Reaction pathway| NEBFeature|   
ReactionStringFeature| [ReactionStringFeature](./matlantis-features/reaction_string.html)|   
Scan| RestScanFeature| [RestScanFeature](./matlantis-features/rest_scan.html)  
  
## ASEとの関係 ​

Atomic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/)) は原子シミュレーションにおいて広く使われているパッケージです。 `matlantis-features` は内部的にASEに依存している箇所はあるものの、多くの相違点があります。

  * 設計思想: 
    * ASEは様々なシミュレーション手法の呼び出しを重視して設計されていますが、 `matlantis-features` は物性値を求めることを重視して設計されています。
  * Calculator: 
    * ASEは任意のcalculatorをセットできますが、 `matlantis-features` はPFP calculator専用に設計されています。
  * コーディングスタイル: 
    * `matlantis-features` は深層学習ライブラリ風のコーディングスタイルを採用しています。 ASEでは物性値を計算するたびに `Atoms` インスタンスを作成する必要がありますが、 `matlantis-features` では1つの `Atoms` インスタンスに対して複数の物理量を計算することができ、 ハイスループット計算に適しています。 
      * ASEの場合:

python
            
            dynamics = Dynamics(atoms)
            result = dynamics()

      * `matlantis-features` の場合:

python
            
            feature = Feature()
            result = feature(atoms)

    * `matlantis-features` では、複雑なタスクをこなしたいときなどに独自のfeatureを定義することができます。 以下に、大量の候補物質に対しスクリーニングすることを想定した独自featureを定義する例を示します。

python
          
          # featureの定義: 構造最適化を行った後、弾性特性および熱化学特性を計算します。
          class MyTaskFeature(FeatureBase):
              def __init__(self):
                  super().__init__()
                  with self.init_scope():
                      self.optimization = FireASEOptFeature()
                      self.elastic_tensor = ElasticTensorFeature()
                      self.elasticity = PostElasticPropertiesFeature()
                      self.force_constant = ForceConstantFeature()
                      self.thermal = PostPhononThermochemistryFeature()
              def __call__(self, input_atoms):
                  opt_result = self.optimization(input_atoms)
                  elastic_tensor_result = self.elastic_tensor(opt_result.atoms)
                  elasticity_result = self.elasticity(elastic_tensor_result)
                  force_constant_result = self.force_constant(opt_result.atoms)
                  thermal_result = self.thermal(force_constant_result)
                  return opt_result, elastic_tensor_result, \
                      elasticity_result, force_constant_result, thermal_result
          
          # 独自定義したfeatureの初期化
          my_task_feature = MyTaskFeature()
          
          # 独自定義したfeatureを用いて、候補物質に対してハイスループット計算を実行
          results = [my_task_feature(atoms) for atoms in material_list]




## matlantis-featuresで使われるestimatorの設定方法 ​

matlantis-featuresで使われるestimatorの設定(PFPの計算モードやモデルバージョンなど)は変更可能です。 環境変数を用いる方法と `estimator_fn` オプションを指定する方法があります。

### 1\. 環境変数を用いてPFPの計算モードやモデルバージョンを変更する方法 ​

#### 注意 ​

matlantis-features v0.7.0 および v0.8.0 には `MATLANTIS_PFP_CALC_MODE` の読み込みにバグがあります。 環境変数を用いて計算モードの変更をしている場合 v0.8.1 かそれより新しい matlantis-features へのアップグレードをお願いします。 [2022/08/15 のリリースノート](/api/docs/ja/release-notes/#_2022-08-15) も併せてご参照ください。

環境変数 `MATLANTIS_PFP_CALC_MODE` および `MATLANTIS_PFP_MODEL_VERSION` を指定することで、 それぞれ計算モードとモデルバージョンを指定することができます。

python
    
    
    from ase.build import molecule
    import os
    from pfp_api_client.pfp.estimator import Estimator
    from pfp_api_client.utils.messages import MessageEnum
    
    from matlantis_features.features.common.opt import FireASEOptFeature
    
    import os
    os.environ["MATLANTIS_PFP_MODEL_VERSION"] = "v2.0.0"
    os.environ["MATLANTIS_PFP_CALC_MODE"] = "pbe_u_plus_d3"
    
    atoms = molecule("CH3OH")
    
    opt = FireASEOptFeature()
    result = opt(atoms)

### 2\. estimator_fnを指定することでPFPの計算モードやモデルバージョンを変更する方法 ​

環境変数を用いる代わりに `estimator_fn` を指定することでも設定を変更できます。 `estimator_fn` はfeaturesの内部で使われるestimatorを定義するfactory methodです。

PFPの計算モードとバージョンを指定するだけであれば、 `pfp_estimator_fn` を用いるのが便利です。 `pfp_estimator_fn` によって作られた `estimator_fn` を各featureのイニシャライザに渡すことで使用します。

python
    
    
    from ase.build import molecule
    
    from pfp_api_client.pfp.estimator import EstimatorCalcMode
    
    from matlantis_features.features.common.opt import FireASEOptFeature
    from matlantis_features.utils.calculators import pfp_estimator_fn
    
    
    atoms = molecule("CH3OH")
    
    opt = FireASEOptFeature(
        estimator_fn=pfp_estimator_fn(
            model_version="v2.0.0",
            calc_mode=EstimatorCalcMode.PBE_U_PLUS_D3,
        ),
    )
    result = opt(atoms)

さらに、`estimator_fn` の定義を独自に行うことでestimatorを詳しく設定することも可能です。 以下の例では、estimatorによるwarning messageを非表示にしています。

python
    
    
    from ase.build import molecule
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    from pfp_api_client.utils.messages import MessageEnum
    
    from matlantis_features.features.common.opt import FireASEOptFeature
    
    
    def estimator_fn() -> Estimator:
        # change the model version.
        estimator = Estimator(
            model_version="v2.0.0",
            calc_mode=EstimatorCalcMode.PBE_U_PLUS_D3,
        )
    
        # disable the specific type of warning.
        estimator.set_message_status(
            message=MessageEnum.ExperimentalElementWarning,
            message_enable=False,
        )
    
        return estimator
    
    atoms = molecule("CH3OH")
    opt = FireASEOptFeature(estimator_fn=estimator_fn)
    result = opt(atoms)

### 注意点 ​

#### estimatorの共有 ​

複数のfeaturesにまたがってestimatorを共有することはできません。 `estimator_fn` を定義する場合は、その内側でestimatorのインスタンスを作るようにしてください。

##### GOOD ​

python
    
    
    def estimator_fn() -> Estimator:
        return Estimator()

##### BAD ​

python
    
    
    estimator = Estimator()
    
    def estimator_fn() -> Estimator:
        return estimator

#### estimator_fnと環境変数の優先順位 ​

環境変数と `estimator_fn` が同時に定義されている場合、 `estimator_fn` の設定が優先して使用されます。

## Index ​

  * [ElasticTensorFeature](./matlantis-features/elastic_tensor.html)
  * [VibrationFeature](./matlantis-features/vibration.html)
  * [PostVibrationGasThermoFeature](./matlantis-features/vibration_thermo.html)
  * [PostMDDiffusionFeature](./matlantis-features/diffusion.html)
  * [PostEMDViscosityFeature](./matlantis-features/viscosity_emd.html)
  * [PostNEMDViscosityFeature](./matlantis-features/viscosity_nemd.html)
  * [ForceConstantFeature](./matlantis-features/force_constant.html)
  * [PostPhononBandFeature](./matlantis-features/phonon_band.html)
  * [PostPhononDOSFeature](./matlantis-features/phonon_dos.html)
  * [PostNEMDThermalConductivityFeature](./matlantis-features/thermal_conduct.html)
  * [ReactionStringFeature](./matlantis-features/reaction_string.html)
  * [RestScanFeature](./matlantis-features/rest_scan.html)

