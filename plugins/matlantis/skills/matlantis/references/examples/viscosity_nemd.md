# Matlantis-features: 粘性係数 (rNEMD)


## 概要

* reverse non-equilibrium MD (rNEMD)という手法により粘性係数を求めるフィーチャー`ComplexNEMDViscosityFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/viscosity_nemd.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;) の紹介です。標準的な分子動力学計算を利用して粘性係数を求める手法 (`viscosity.ipynb`) を参考にしてください。

* rNEMDでは、系に人工的な速度勾配をつけて、そこからの変化の度合いを測定することで粘性を計算する手法です。速度勾配の強さという計算に必要なパラメータが増えますが、一般により少ない計算時間で計算結果の統計的なばらつきを抑えることが可能となります。

* 計算手法の特徴として、`viscosity.ipynb`とは異なり、MD計算自体に操作を加える必要があります。そのためrNEMDではMD計算とポスト処理をまとめた`ComplexNEMDViscosityFeature`を提供しています。


## 初期設定

* 必要なライブラリ等の準備を行います。


```python
# If you have already installed matlantis-features, you can skip this.
# Version 0.8.1 or later is required to run this notebook.
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


```python
import logging
import pathlib
import numpy as np
from plotly.offline import iplot
from IPython.display import clear_output

from ase.io import read
from matlantis_features.atoms import MatlantisAtoms

from matlantis_features.features.md import ASEMDSystem, NVTBerendsenIntegrator
from matlantis_features.features.md import MDFeature
from matlantis_features.features.md import (
    ComplexNEMDViscosityFeature,
)

try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

```


## 材料の用意

* これ以降は、`viscosity.ipynb`と同様の系に対して`ComplexNEMDViscosityFeature`を使用した例を紹介します。


```python
atoms = read(str(dir_path/"assets/viscosity/n-decane.xyz")) * (1,1,2)
```


## 分子動力学計算

* ここでは`MDFeature`の代わりに`ComplexNEMDViscosityFeature`を実行することで動力学計算の結果を得ています。


```python
integ = NVTBerendsenIntegrator(timestep=1.0, temperature=480.0)
mdsys = ASEMDSystem(atoms)
mdsys.init_temperature(480.0)

feature = ComplexNEMDViscosityFeature(
    integrator=integ, n_run=10000, rnemd_n_slab=20, rnemd_interval=50, init_time=2000,
    md_feature_args={
        "traj_file_name": "md_viscosity_nemd.traj",
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


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/md_viscosity_nemd.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 463.23 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   100  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 487.80 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   200  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 479.44 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   300  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 482.20 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   400  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 471.71 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   500  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 471.30 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   600  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 478.56 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   700  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 476.92 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 471.70 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:   900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 487.96 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 480.32 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 481.93 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 473.18 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 490.28 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1400  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 475.09 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 500.37 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 494.95 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 486.72 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 482.82 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  1900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 488.34 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 492.30 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2100  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 479.31 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2200  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 481.79 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2300  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 474.35 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 510.39 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 510.55 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2600  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 478.91 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 498.12 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 516.06 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  2900  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 496.61 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 491.00 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3100  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 487.52 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 498.28 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3300  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 484.13 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 502.29 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 521.05 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 502.83 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3700  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 496.89 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 514.50 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  3900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 502.97 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 493.53 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4100  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 494.37 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4200  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 493.31 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 509.34 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 503.12 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 490.80 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4600  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 491.53 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 500.18 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4800  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 480.09 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  4900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 497.64 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 498.68 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5100  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 484.98 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 501.07 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 499.89 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 505.01 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 500.08 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 497.30 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 517.34 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 482.69 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  5900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 484.39 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 496.27 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 484.24 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 504.92 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 489.95 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 480.87 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6500  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 468.08 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 478.08 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.43 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 511.42 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  6900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 496.60 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 503.91 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 487.19 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 505.01 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 497.51 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 510.36 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 491.71 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 499.56 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 496.44 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 532.73 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  7900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.62 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 496.71 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 471.61 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 504.66 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 497.40 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8400  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 455.26 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.59 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 507.85 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8700  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 492.48 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 513.42 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  8900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.57 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 503.77 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 477.25 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 514.69 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 481.83 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 507.98 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 489.25 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 485.92 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.86 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 501.28 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps:  9900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 471.14 K  volume:  7701 Ang^3  density: 0.614 g/cm^3
    steps: 10000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 490.13 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:     0  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 483.84 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 484.51 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 507.87 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   300  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 485.31 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 491.41 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   500  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 478.02 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 495.54 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 504.27 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   800  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 470.94 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:   900  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 480.55 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 482.89 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1100  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 476.23 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 500.98 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.89 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 497.24 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 489.05 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 490.82 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 483.38 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 505.04 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  1900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 505.49 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 508.62 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 502.29 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 489.95 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 478.16 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 507.40 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.55 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 503.86 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.92 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 492.91 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  2900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 497.96 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 490.24 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 492.28 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 491.65 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.96 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 517.33 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 490.55 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 479.72 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 488.38 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 506.85 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  3900  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 475.53 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 491.99 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 488.84 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4200  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 476.56 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 492.20 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4400  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 484.79 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 497.35 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4600  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 480.62 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4700  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 483.01 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 512.05 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  4900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 492.76 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 477.55 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5100  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 467.51 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 489.92 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 491.63 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 501.92 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 500.59 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 475.74 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 500.88 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 477.77 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  5900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 475.42 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 488.54 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 479.09 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 491.93 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 505.04 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 481.54 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6500  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 474.62 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6600  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 478.82 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6700  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 472.65 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 507.66 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  6900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 485.31 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 483.48 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 500.06 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 485.31 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 486.19 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7400  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 494.30 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7500  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 492.43 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 508.75 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 508.45 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7800  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 480.65 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  7900  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 482.22 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 502.04 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 510.83 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8200  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 486.51 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 491.91 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 484.72 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8500  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 498.80 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 500.71 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 489.93 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 498.09 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  8900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 488.19 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 490.93 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9100  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 482.30 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9200  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 485.44 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9300  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 489.75 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9400  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 483.23 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9500  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 466.13 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9600  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 488.16 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9700  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 499.90 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9800  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.21 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps:  9900  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 488.36 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


    steps: 10000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 508.46 K  volume:  7701 Ang^3  density: 0.614 g/cm^3


## 粘性係数の計算

* 必要に応じて、通常の分子動力学計算を利用して粘性係数を求める`viscosity.ipynb`も参照してください。


```python
result.post_nemd_viscosity_result.viscosity
```


    {'amu/A/fs': [0.011794349326067996],
     'Pa s': [0.000195849775073336],
     'mPa s': [0.19584977507333595]}


## 速度プロファイルのプロット

* rNEMDによって空間的に速度勾配がついています。そのことを確認してみましょう。

* 初回起動時などで図が表示されない場合、ブラウザでページをリロードすると表示されることがあります。


```python
fig = result.post_nemd_viscosity_result.plot()
fig.update_layout(title="velocity profile in rNEMD")
fig.update_layout(width=800, height=500)
iplot(fig)
```


