# Matlantis-features: 熱伝導率 (rNEMD)


## 概要
物質の熱の伝わりやすさを表す基本的な物理量、それが熱伝導率です。

熱伝導率を原子シミュレーションで計算する方法として、実時間でのフォノン計算、グリーン-久保公式、非平衡MDなどの方法があります。
``matlantis-features`` では、 reverse non-equilibrium molecular dynamics (rNEMD) 法を提供しています。
rNEMD 法の詳細については、&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/thermal_conduct.html" target="_blank"&gt;Matlantis-Guidebook&lt;/a&gt; を参照してください。

この例では：
* reverse non-equilibrium MD (rNEMD)という手法により熱伝導率を求めるフィーチャー`ComplexNEMDThermalConductivityFeature`の紹介です。


## 初期設定

* 必要なライブラリ等の準備を行います。


```python
# If you have already installed matlantis-features, you can skip this.
#!pip install 'matlantis-features&gt;=0.8.1'
```


```python
import matlantis_features
import pfp_api_client
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator()
print(f"matlantis_features: {matlantis_features.__version__}")
print(f"pfp_api_client: {pfp_api_client.__version__}")
print(f"current pfp model version: {estimator.model_version}")
```

    matlantis_features: 1.0.1
    pfp_api_client: 2.0.1
    current pfp model version: v8.0.0


### `estimator_fn` による計算モードとモデルバージョンの指定

* Feature に `estimator_fn` 引数を使うことで、matlantis-features で使われる PFP の計算モードとモデルバージョンを指定することができます。

* `estimator_fn` は、Estimator オブジェクトを作る factory method です。計算モードとモデルバージョンのみを指定する場合は、以下のように `pfp_estimator_fn` を利用できます。より細かい指定が必要な場合は、自身で factory method を定義できます。詳細に関しては [Matlantis Guidebook](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html#estimator-fnpfp) をご参照ください。

* `estimator_fn` が指定されない場合は、環境変数で指定された値が使用されます。`estimator_fn` も環境変数も指定されない場合は、デフォルトのモデルバージョンと計算モード（`PBE`）が使用されます。

* 環境変数と `estimator_fn` が同時に定義されている場合、 `estimator_fn` の設定が優先して使用されます。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode

estimator_fn = pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
```


## 材料の用意

`Thermal conductivity (water)` の例では、液体の熱伝導率の計算に rNEMD 法を使用できることを示します。 ここでは、固体材料アモルファス SiO2 に使用します。

アモルファス SiO2 の構造は、melt quenching法、つまり結晶構造を非常に高温に加熱し、その後システムを急冷して無秩序な構造を得る方法から得ることができます。
この例では、事前に準備された構造体をロードするだけです。

**Note**:シミュレーション領域はフォノンの平均自由行程よりも十分広くなければいけません。そうでないと、熱伝導率が大幅に過小評価されてしまいます。
結晶中のフォノンは非常に大きい平均自由行程をもつので、大きなシミュレーション領域が必要になります。
しかし、PFPで計算できる系のサイズの限界から、完全な結晶やアモルファスのような無秩序結晶を除く結晶での熱伝導率計算は困難です。


```python
import logging
import pathlib
import numpy as np
from plotly.offline import iplot
from IPython.display import clear_output

from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.md import ASEMDSystem, NVTBerendsenIntegrator
from matlantis_features.features.md import (
    ComplexNEMDThermalConductivityFeature
)

try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()
    
logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

```


```python
atoms_SiO2 = MatlantisAtoms.from_file(str(dir_path/"assets/thermal_conductivity/a_sio2.xyz"))
```


## rNEMD

rNEMD法では、原子の速度を交換することで人為的な温度勾配を発生させます。 ``PostNEMDViscosityFeature`` と類似の手法を用いています。

rNEMD法は以下の手順を繰り返すことで、系全体での温度を変えることなく、シミュレーション領域の中央部と下部の間で温度勾配を発生させます。

1. シミュレーション対象の領域ををある方向（通常はz方向）に沿ってN個のスラブに分割します。
2. 一番下のスラブから最も温度が高い原子（運動エネルギーが最も大きい原子）を見つけ出します。
3. その原子と同じ元素で、なおかつ最小の運動エネルギーをもつ原子を中央のスラブから見つけ出します。
4. この2つの原子を交換することで、高温原子を下部から中央のスラブに移動させます。

rNEMD法を実行するには、MDの手順を変更する必要があります。
これは、 ``RNEMDExtension`` を用いることで可能です。


```python
integ = NVTBerendsenIntegrator(timestep=1.0, temperature=293.15)
mdsys = ASEMDSystem(atoms_SiO2)
mdsys.init_temperature(293.15)

feature = ComplexNEMDThermalConductivityFeature(
    integrator=integ,
    n_run=20000,
    rnemd_n_slab=20,
    rnemd_interval=50,
    init_time=2000, 
    md_feature_args={
        "traj_file_name": "md_thermal_conductivity2.traj",
        "traj_freq": 100,
        "show_progress_bar": True,
        "show_logger": True,
        "logger_interval": 100,
        "estimator_fn": estimator_fn,
    }
)
result = feature(mdsys)

clear_output(wait=True)
```


      0%|          | 0/20000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/md_thermal_conductivity2.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-6.23 eV/atom  total energy: -6.19 eV/atom  temperature: 292.24 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   100  energy：-6.23 eV/atom  total energy: -6.19 eV/atom  temperature: 351.23 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   200  energy：-6.23 eV/atom  total energy: -6.19 eV/atom  temperature: 338.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   300  energy：-6.24 eV/atom  total energy: -6.19 eV/atom  temperature: 347.32 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   400  energy：-6.24 eV/atom  total energy: -6.19 eV/atom  temperature: 340.38 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   500  energy：-6.24 eV/atom  total energy: -6.19 eV/atom  temperature: 327.06 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   600  energy：-6.24 eV/atom  total energy: -6.19 eV/atom  temperature: 322.08 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 324.10 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 319.23 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:   900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 317.16 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 322.80 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 306.07 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 306.60 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 313.45 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 312.53 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 303.83 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.49 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 302.45 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 302.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  1900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.29 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.53 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.23 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.10 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.73 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 303.29 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.81 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.08 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.71 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  2900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.37 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 302.21 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.51 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.51 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.18 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.30 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.45 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.16 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  3900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.34 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.04 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.05 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.80 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 285.03 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.08 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.61 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.71 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.42 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  4900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.54 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.13 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.34 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.31 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.72 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.46 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 285.17 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.65 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.04 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  5900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.21 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.21 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.48 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.46 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.80 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.89 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 286.09 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.42 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.41 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  6900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.09 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.68 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.92 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.03 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.59 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.88 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.74 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.92 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.62 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.43 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  7900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.30 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.59 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.14 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.19 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.76 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.19 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.81 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.99 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.82 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  8900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.29 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.66 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.83 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.38 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.23 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.67 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.66 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.95 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.63 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.80 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps:  9900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.54 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.30 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.10 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.19 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.80 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.68 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.72 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.83 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.55 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.68 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 10900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.52 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.05 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.66 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.99 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.87 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.39 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.47 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 305.20 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.33 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.52 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 11900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 282.05 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.78 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.09 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.53 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.09 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.44 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.08 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.53 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.76 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 12900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.41 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.06 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.38 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.17 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.00 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 284.80 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.99 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.24 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.15 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 13900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.88 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.59 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.27 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.96 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.62 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.85 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.97 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.00 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.87 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.35 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 14900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.48 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.01 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.90 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.36 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.48 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.57 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.11 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.35 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.55 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.42 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 15900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.87 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.10 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.18 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.27 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.61 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 286.82 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.63 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.36 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 16900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.81 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.86 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.82 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.67 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.65 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.32 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.25 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.82 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.97 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.93 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 17900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.79 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.82 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 285.56 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.86 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 304.18 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.12 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 304.67 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.59 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 285.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.41 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 18900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.74 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 282.10 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.85 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.89 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.43 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.85 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.64 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.68 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.16 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 286.70 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 19900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.99 K  volume: 26428 Ang^3  density: 2.174 g/cm^3
    steps: 20000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 286.63 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-6.23 eV/atom  total energy: -6.19 eV/atom  temperature: 296.70 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   100  energy：-6.24 eV/atom  total energy: -6.19 eV/atom  temperature: 358.13 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   200  energy：-6.23 eV/atom  total energy: -6.19 eV/atom  temperature: 332.89 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   300  energy：-6.24 eV/atom  total energy: -6.19 eV/atom  temperature: 341.83 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   400  energy：-6.23 eV/atom  total energy: -6.19 eV/atom  temperature: 324.00 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   500  energy：-6.24 eV/atom  total energy: -6.19 eV/atom  temperature: 332.56 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   600  energy：-6.24 eV/atom  total energy: -6.19 eV/atom  temperature: 334.42 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 321.50 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 320.09 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:   900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 313.57 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 318.40 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 321.23 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 304.38 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 306.00 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 315.04 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.07 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 305.47 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.61 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 305.73 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  1900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 302.23 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.00 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 305.07 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.47 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 303.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.64 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.69 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.87 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.54 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.70 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  2900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.59 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.05 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.27 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.00 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 302.95 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.29 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.16 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.53 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.41 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  3900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.51 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.47 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.13 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.85 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.46 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.38 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.34 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.38 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.00 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.23 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  4900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.36 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.87 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 302.26 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.06 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.77 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.53 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.07 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.55 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.62 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.03 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  5900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.26 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.30 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.02 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.92 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.42 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.21 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.19 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 302.90 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.83 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.61 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  6900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.52 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.72 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.46 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.25 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.55 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.76 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.69 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.43 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.78 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.89 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  7900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.95 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.60 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 306.05 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.15 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.10 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 284.16 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.56 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.38 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.57 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.41 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  8900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 284.19 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.70 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.33 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.06 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.02 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 286.78 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.39 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.44 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.85 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.39 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps:  9900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.86 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.61 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.54 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.33 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.57 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.04 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.05 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.83 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.03 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.11 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 10900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.59 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.76 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.07 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.89 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.82 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.44 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.92 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.81 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.86 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.35 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 11900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.33 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.57 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.50 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.14 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.25 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 303.11 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.86 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.93 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.63 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.06 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 12900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.61 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.90 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 286.94 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.91 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.46 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.49 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.38 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.05 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.70 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.33 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 13900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.50 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.09 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.28 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.19 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.13 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.64 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.84 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.44 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.23 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.10 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 14900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 282.93 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.91 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.08 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.62 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.59 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.88 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.30 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.33 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.20 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 281.82 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 15900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.16 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.25 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.27 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 286.29 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.00 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.56 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 286.59 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.89 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.10 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.56 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 16900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.62 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 300.91 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.07 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.85 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 299.11 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.80 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.81 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 302.41 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 288.44 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.42 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 17900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 289.82 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 297.25 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 296.39 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.72 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 298.39 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 301.42 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.54 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.51 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.71 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 291.87 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 18900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 304.49 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.89 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19100  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 290.25 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19200  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.28 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19300  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 287.39 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19400  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 292.30 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19500  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.97 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19600  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.26 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19700  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 294.86 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19800  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.86 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 19900  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 295.17 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


    steps: 20000  energy：-6.24 eV/atom  total energy: -6.20 eV/atom  temperature: 293.88 K  volume: 26428 Ang^3  density: 2.174 g/cm^3


## 熱伝導率の計算

rNEMD計算を行った後、 ``PostNEMDThermalConductivityFeature`` を使うことでrNEMDのトラジェクトリを分析できます。
熱伝導率は以下のように計算されます:

$$\lambda = - \frac{\sum E_{kinetic, hot} - E_{kinetic, cold}}{2 t S_{xy} &lt;\frac{\partial T}{\partial z}&gt;} $$

$\lambda$ は熱伝導率、$t$ は時間、 $S_{xy}$ はxy面の面積、 $\frac{\partial T}{\partial z}$ は $z$ 方向に沿った温度勾配を表します。


```python
result.post_nemd_thermal_conductivity_result.thermal_conductivity
```


    {'eV/fs/Angstrom/K': [5.505405321251307e-07],
     'J/s/m/K': [0.8820631693736758],
     'W/m/K': [0.8820631693736758]}


## 温度プロファイルのプロット

熱伝導率は、z 方向に沿った温度勾配に直接関係しています。 ここでは、温度プロファイルをプロットします。


```python
fig = result.post_nemd_thermal_conductivity_result.plot()
fig.update_layout(title="temperature profile in rNEMD")
fig.update_layout(width=800, height=500)
iplot(fig)
```


## 実験値との比較

* The thermal conductivity of amorphous SiO2 at 293.15K is about 1.2~1.3 W/(m K). (1. Callard, S., Tallarida, G., Borghesi, A. &amp; Zanotti, L. Journal of non-crystalline solids 245, 203–209 (1999). 2. Cao, S., He, H. &amp; Zhu, W. Aip Advances 7, 015038 (2017).)
