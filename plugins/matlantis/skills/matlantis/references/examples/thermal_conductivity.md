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
!pip install 'matlantis-features&gt;=0.8.1'
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

rNEMD 法は、液体の熱伝導率の計算に適しています。 ここでは例として水を使用します。

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
atoms_water = MatlantisAtoms.from_file(str(dir_path/"assets/thermal_conductivity/water.xyz"))
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
mdsys = ASEMDSystem(atoms_water)
mdsys.init_temperature(293.15)

feature = ComplexNEMDThermalConductivityFeature(
    integrator=integ,
    n_run=10000,
    rnemd_n_slab=20,
    rnemd_interval=50,
    init_time=2000,
    md_feature_args={
        "traj_file_name": "md_thermal_conductivity.traj",
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


      0%|          | 0/10000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/md_thermal_conductivity.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-3.48 eV/atom  total energy: -3.44 eV/atom  temperature: 304.33 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   100  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 372.25 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   200  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 346.81 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   300  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 348.29 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   400  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 373.63 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   500  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 366.32 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   600  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 350.63 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   700  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 377.13 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   800  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 375.17 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:   900  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 376.00 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1000  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 368.35 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1100  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 360.32 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1200  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 342.37 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1300  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 372.38 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1400  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 361.41 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1500  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 394.00 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1600  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 398.81 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1700  energy：-3.48 eV/atom  total energy: -3.44 eV/atom  temperature: 339.92 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1800  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 401.13 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  1900  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 380.12 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2000  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 364.83 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2100  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 348.01 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2200  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 364.62 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2300  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 374.32 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2400  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 365.28 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2500  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 313.99 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2600  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 326.70 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2700  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 341.48 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2800  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 332.40 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  2900  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 338.27 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3000  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 317.76 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3100  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 316.81 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3200  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 331.40 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3300  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 331.36 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3400  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 359.07 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3500  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 320.29 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3600  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 326.15 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3700  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 349.49 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3800  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 380.87 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  3900  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 385.74 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4000  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 352.17 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4100  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 386.72 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4200  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 377.37 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4300  energy：-3.50 eV/atom  total energy: -3.44 eV/atom  temperature: 396.03 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4400  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 362.23 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4500  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 372.32 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4600  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 310.25 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4700  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 372.58 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4800  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 345.88 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  4900  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 326.63 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5000  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 322.80 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5100  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 338.32 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5200  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 347.36 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5300  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 311.64 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5400  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 319.20 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5500  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 344.13 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5600  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 362.46 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5700  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 339.73 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5800  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 351.44 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  5900  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 342.61 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6000  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 365.92 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6100  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 355.85 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6200  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 345.86 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6300  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 346.31 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6400  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 331.34 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6500  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 368.10 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6600  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 330.88 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6700  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 325.90 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6800  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 327.32 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  6900  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 326.54 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7000  energy：-3.50 eV/atom  total energy: -3.44 eV/atom  temperature: 400.49 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7100  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 372.14 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7200  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 364.05 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7300  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 331.32 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7400  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 337.43 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7500  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 342.67 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7600  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 343.55 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7700  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 335.56 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7800  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 340.48 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  7900  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 334.78 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8000  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 346.13 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8100  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 352.18 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8200  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 344.00 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8300  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 362.80 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8400  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 349.02 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8500  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 339.84 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8600  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 325.40 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8700  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 325.16 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8800  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 337.78 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  8900  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 363.01 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9000  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 332.53 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9100  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 323.33 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9200  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 328.67 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9300  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 334.08 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9400  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 392.23 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9500  energy：-3.48 eV/atom  total energy: -3.44 eV/atom  temperature: 372.39 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9600  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 400.36 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9700  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 404.54 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9800  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 393.33 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps:  9900  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 410.35 K  volume:  2429 Ang^3  density: 0.998 g/cm^3
    steps: 10000  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 398.49 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-3.48 eV/atom  total energy: -3.44 eV/atom  temperature: 320.69 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   100  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 393.55 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   200  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 403.98 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   300  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 381.61 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   400  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 364.23 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   500  energy：-3.49 eV/atom  total energy: -3.43 eV/atom  temperature: 397.71 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   600  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 393.46 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   700  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 396.74 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   800  energy：-3.49 eV/atom  total energy: -3.43 eV/atom  temperature: 431.97 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:   900  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 387.89 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1000  energy：-3.48 eV/atom  total energy: -3.42 eV/atom  temperature: 398.66 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1100  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 414.13 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1200  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 420.35 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1300  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 405.11 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1400  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 358.27 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1500  energy：-3.48 eV/atom  total energy: -3.44 eV/atom  temperature: 367.79 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1600  energy：-3.48 eV/atom  total energy: -3.42 eV/atom  temperature: 399.89 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1700  energy：-3.48 eV/atom  total energy: -3.42 eV/atom  temperature: 432.97 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1800  energy：-3.48 eV/atom  total energy: -3.42 eV/atom  temperature: 423.37 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  1900  energy：-3.48 eV/atom  total energy: -3.42 eV/atom  temperature: 401.71 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2000  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 399.49 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2100  energy：-3.48 eV/atom  total energy: -3.42 eV/atom  temperature: 428.07 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2200  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 384.05 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2300  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 383.79 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2400  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 378.35 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2500  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 376.80 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2600  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 362.91 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2700  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 377.51 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2800  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 394.11 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  2900  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 393.06 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3000  energy：-3.48 eV/atom  total energy: -3.43 eV/atom  temperature: 378.51 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3100  energy：-3.48 eV/atom  total energy: -3.44 eV/atom  temperature: 377.78 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3200  energy：-3.49 eV/atom  total energy: -3.43 eV/atom  temperature: 410.50 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3300  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 408.89 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3400  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 364.32 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3500  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 375.34 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3600  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 387.40 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3700  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 359.12 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3800  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 356.72 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  3900  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 381.34 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4000  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 357.50 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4100  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 312.62 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4200  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 333.84 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4300  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 318.62 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4400  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 338.79 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4500  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 330.46 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4600  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 323.93 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4700  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 339.72 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4800  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 353.90 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  4900  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 351.71 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5000  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 387.73 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5100  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 363.30 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5200  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 333.94 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5300  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 336.17 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5400  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 356.15 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5500  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 310.01 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5600  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 343.46 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5700  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 330.40 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5800  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 339.57 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  5900  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 359.45 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6000  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 357.54 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6100  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 404.93 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6200  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 355.96 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6300  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 378.12 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6400  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 348.56 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6500  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 351.25 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6600  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 341.06 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6700  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 345.44 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6800  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 361.19 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  6900  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 339.08 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7000  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 359.43 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7100  energy：-3.49 eV/atom  total energy: -3.44 eV/atom  temperature: 375.97 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7200  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 359.45 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7300  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 333.09 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7400  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 354.96 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7500  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 365.00 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7600  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 357.58 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7700  energy：-3.49 eV/atom  total energy: -3.45 eV/atom  temperature: 326.37 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7800  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 368.43 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  7900  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 367.61 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8000  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 364.97 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8100  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 334.06 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8200  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 372.21 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8300  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 338.13 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8400  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 348.43 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8500  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 344.33 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8600  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 351.16 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8700  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 370.93 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8800  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 333.77 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  8900  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 338.34 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9000  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 351.51 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9100  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 334.25 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9200  energy：-3.50 eV/atom  total energy: -3.46 eV/atom  temperature: 367.39 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9300  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 370.96 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9400  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 337.06 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9500  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 389.80 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9600  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 348.91 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9700  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 354.24 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9800  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 376.86 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps:  9900  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 352.07 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


    steps: 10000  energy：-3.50 eV/atom  total energy: -3.45 eV/atom  temperature: 329.53 K  volume:  2429 Ang^3  density: 0.998 g/cm^3


## 熱伝導率の計算

rNEMD計算を行った後、 ``PostNEMDThermalConductivityFeature`` を使うことでrNEMDのトラジェクトリを分析できます。
熱伝導率は以下のように計算されます:

$$\lambda = - \frac{\sum E_{kinetic, hot} - E_{kinetic, cold}}{2 t S_{xy} &lt;\frac{\partial T}{\partial z}&gt;} $$

$\lambda$ は熱伝導率、$t$ は時間、 $S_{xy}$ はxy面の面積、 $\frac{\partial T}{\partial z}$ は $z$ 方向に沿った温度勾配を表します。


```python
result.post_nemd_thermal_conductivity_result.thermal_conductivity
```


    {'eV/fs/Angstrom/K': [5.420904143957603e-07],
     'J/s/m/K': [0.8685245883046709],
     'W/m/K': [0.8685245883046709]}


## 温度プロファイルのプロット

熱伝導率は、z 方向に沿った温度勾配に直接関係しています。 ここでは、温度プロファイルをプロットします。


```python
fig = result.post_nemd_thermal_conductivity_result.plot()
fig.update_layout(title="temperature profile in rNEMD")
fig.update_layout(width=800, height=500)
iplot(fig)
```


## 実験値との比較

* The water thermal conductivity at 293.15K is about 0.6 W/(m K). (Maria L. V. Ramires et al. Journal of Physical and Chemical Reference Data 24, 1377 (1995))
