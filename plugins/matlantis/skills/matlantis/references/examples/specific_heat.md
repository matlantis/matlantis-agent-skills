# Matlantis-features: 比熱


## 概要

* この Notebook では、 Si 結晶を対象として比熱の計算方法を紹介します。
* 複数の温度で生成された軌跡を入力として各温度での比熱を計算するフィーチャーである PostMDSpecificHeatFeature を紹介します。この Example では、MDFeature を用いて軌跡を生成します。

## 初期設定

* はじめに、必要なライブラリ等の準備を行います。


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


```python
import logging
import pathlib
import numpy as np
import plotly.graph_objects as go
from plotly.offline import iplot
from IPython.display import clear_output

from ase import units
from ase.build import bulk

from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.md import ASEMDSystem, MDFeature, NPTBerendsenIntegrator
from matlantis_features.features.md.post_features import PostMDSpecificHeatFeature

try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

```


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


## 関数の用意
* MDFeature を利用して、分子動力学計算の実行と対応する軌跡ファイルの生成を行う関数を定義しておきます。


```python
def run_md(temperature):
    atoms = MatlantisAtoms(bulk("Si", cubic=True) * (2, 2, 2))

    npt = NPTBerendsenIntegrator(
        timestep=1.0, temperature=temperature, pressure=101315*units.Pascal, taut=50, taup=50
    )
    mdsys = ASEMDSystem(atoms)
    mdsys.init_temperature(temperature)

    md = MDFeature(
        integrator=npt,
        traj_file_name=f"output/md_specific_heat_{str(temperature).zfill(4)}K.traj",
        n_run=2000, show_progress_bar=True,
        show_logger=True,
        logger_interval=100,
        estimator_fn=estimator_fn,
    )
    md_res = md(mdsys)
    return md, md_res

```


```python
md_features = []
res = []

for i, t in enumerate([200, 300, 400, 500, 600, 700]):
    logger.info(f"Running MD T={t}")
    md, md_res = run_md(t)
    md_features.append(md)
    res.append(md_res)

```


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0200K.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.52 eV/atom  temperature: 230.14 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:   100  energy：-4.53 eV/atom  total energy: -4.51 eV/atom  temperature: 193.35 K  volume:  1308 Ang^3  density: 2.281 g/cm^3
    steps:   200  energy：-4.52 eV/atom  total energy: -4.50 eV/atom  temperature: 133.98 K  volume:  1306 Ang^3  density: 2.285 g/cm^3
    steps:   300  energy：-4.52 eV/atom  total energy: -4.50 eV/atom  temperature: 160.09 K  volume:  1306 Ang^3  density: 2.285 g/cm^3
    steps:   400  energy：-4.52 eV/atom  total energy: -4.50 eV/atom  temperature: 182.11 K  volume:  1308 Ang^3  density: 2.282 g/cm^3
    steps:   500  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 197.45 K  volume:  1309 Ang^3  density: 2.281 g/cm^3
    steps:   600  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 239.77 K  volume:  1309 Ang^3  density: 2.281 g/cm^3
    steps:   700  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 191.80 K  volume:  1308 Ang^3  density: 2.282 g/cm^3
    steps:   800  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 197.45 K  volume:  1308 Ang^3  density: 2.283 g/cm^3
    steps:   900  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 222.76 K  volume:  1309 Ang^3  density: 2.280 g/cm^3
    steps:  1000  energy：-4.52 eV/atom  total energy: -4.50 eV/atom  temperature: 173.50 K  volume:  1309 Ang^3  density: 2.281 g/cm^3
    steps:  1100  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 189.54 K  volume:  1306 Ang^3  density: 2.286 g/cm^3
    steps:  1200  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 213.03 K  volume:  1307 Ang^3  density: 2.283 g/cm^3
    steps:  1300  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 247.48 K  volume:  1309 Ang^3  density: 2.280 g/cm^3
    steps:  1400  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 208.10 K  volume:  1308 Ang^3  density: 2.281 g/cm^3
    steps:  1500  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 219.85 K  volume:  1308 Ang^3  density: 2.282 g/cm^3
    steps:  1600  energy：-4.52 eV/atom  total energy: -4.50 eV/atom  temperature: 145.52 K  volume:  1308 Ang^3  density: 2.281 g/cm^3
    steps:  1700  energy：-4.52 eV/atom  total energy: -4.50 eV/atom  temperature: 161.55 K  volume:  1307 Ang^3  density: 2.283 g/cm^3
    steps:  1800  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 201.76 K  volume:  1308 Ang^3  density: 2.282 g/cm^3
    steps:  1900  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 208.61 K  volume:  1310 Ang^3  density: 2.279 g/cm^3
    steps:  2000  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 204.20 K  volume:  1308 Ang^3  density: 2.282 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0300K.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.51 eV/atom  temperature: 316.26 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:   100  energy：-4.52 eV/atom  total energy: -4.49 eV/atom  temperature: 264.82 K  volume:  1310 Ang^3  density: 2.278 g/cm^3
    steps:   200  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 199.36 K  volume:  1308 Ang^3  density: 2.282 g/cm^3
    steps:   300  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 257.27 K  volume:  1307 Ang^3  density: 2.283 g/cm^3
    steps:   400  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 251.67 K  volume:  1310 Ang^3  density: 2.279 g/cm^3
    steps:   500  energy：-4.52 eV/atom  total energy: -4.48 eV/atom  temperature: 320.59 K  volume:  1310 Ang^3  density: 2.278 g/cm^3
    steps:   600  energy：-4.52 eV/atom  total energy: -4.48 eV/atom  temperature: 345.37 K  volume:  1312 Ang^3  density: 2.274 g/cm^3
    steps:   700  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 318.03 K  volume:  1311 Ang^3  density: 2.276 g/cm^3
    steps:   800  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 270.97 K  volume:  1310 Ang^3  density: 2.279 g/cm^3
    steps:   900  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 269.70 K  volume:  1310 Ang^3  density: 2.278 g/cm^3
    steps:  1000  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 263.67 K  volume:  1310 Ang^3  density: 2.278 g/cm^3
    steps:  1100  energy：-4.52 eV/atom  total energy: -4.48 eV/atom  temperature: 335.77 K  volume:  1307 Ang^3  density: 2.283 g/cm^3
    steps:  1200  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 294.76 K  volume:  1311 Ang^3  density: 2.277 g/cm^3
    steps:  1300  energy：-4.52 eV/atom  total energy: -4.47 eV/atom  temperature: 330.57 K  volume:  1312 Ang^3  density: 2.274 g/cm^3
    steps:  1400  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 282.29 K  volume:  1310 Ang^3  density: 2.278 g/cm^3
    steps:  1500  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 281.21 K  volume:  1310 Ang^3  density: 2.278 g/cm^3
    steps:  1600  energy：-4.52 eV/atom  total energy: -4.48 eV/atom  temperature: 348.55 K  volume:  1310 Ang^3  density: 2.278 g/cm^3
    steps:  1700  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 287.57 K  volume:  1309 Ang^3  density: 2.280 g/cm^3
    steps:  1800  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 273.16 K  volume:  1310 Ang^3  density: 2.279 g/cm^3
    steps:  1900  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 289.73 K  volume:  1311 Ang^3  density: 2.277 g/cm^3
    steps:  2000  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 318.69 K  volume:  1309 Ang^3  density: 2.280 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0400K.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.49 eV/atom  temperature: 444.13 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:   100  energy：-4.51 eV/atom  total energy: -4.46 eV/atom  temperature: 349.19 K  volume:  1313 Ang^3  density: 2.273 g/cm^3
    steps:   200  energy：-4.49 eV/atom  total energy: -4.46 eV/atom  temperature: 271.05 K  volume:  1309 Ang^3  density: 2.280 g/cm^3
    steps:   300  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 333.98 K  volume:  1308 Ang^3  density: 2.283 g/cm^3
    steps:   400  energy：-4.49 eV/atom  total energy: -4.45 eV/atom  temperature: 326.55 K  volume:  1312 Ang^3  density: 2.276 g/cm^3
    steps:   500  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 405.08 K  volume:  1314 Ang^3  density: 2.271 g/cm^3
    steps:   600  energy：-4.51 eV/atom  total energy: -4.45 eV/atom  temperature: 452.14 K  volume:  1315 Ang^3  density: 2.270 g/cm^3
    steps:   700  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 383.16 K  volume:  1314 Ang^3  density: 2.272 g/cm^3
    steps:   800  energy：-4.49 eV/atom  total energy: -4.45 eV/atom  temperature: 353.14 K  volume:  1312 Ang^3  density: 2.276 g/cm^3
    steps:   900  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 390.13 K  volume:  1313 Ang^3  density: 2.274 g/cm^3
    steps:  1000  energy：-4.49 eV/atom  total energy: -4.45 eV/atom  temperature: 359.37 K  volume:  1312 Ang^3  density: 2.274 g/cm^3
    steps:  1100  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 426.61 K  volume:  1309 Ang^3  density: 2.280 g/cm^3
    steps:  1200  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 384.09 K  volume:  1312 Ang^3  density: 2.276 g/cm^3
    steps:  1300  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 440.04 K  volume:  1314 Ang^3  density: 2.272 g/cm^3
    steps:  1400  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 449.55 K  volume:  1313 Ang^3  density: 2.274 g/cm^3
    steps:  1500  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 386.37 K  volume:  1313 Ang^3  density: 2.273 g/cm^3
    steps:  1600  energy：-4.51 eV/atom  total energy: -4.45 eV/atom  temperature: 432.45 K  volume:  1313 Ang^3  density: 2.274 g/cm^3
    steps:  1700  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 366.23 K  volume:  1312 Ang^3  density: 2.276 g/cm^3
    steps:  1800  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 397.64 K  volume:  1312 Ang^3  density: 2.276 g/cm^3
    steps:  1900  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 408.58 K  volume:  1312 Ang^3  density: 2.274 g/cm^3
    steps:  2000  energy：-4.49 eV/atom  total energy: -4.45 eV/atom  temperature: 371.88 K  volume:  1312 Ang^3  density: 2.276 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0500K.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.47 eV/atom  temperature: 598.43 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:   100  energy：-4.49 eV/atom  total energy: -4.44 eV/atom  temperature: 418.49 K  volume:  1315 Ang^3  density: 2.270 g/cm^3
    steps:   200  energy：-4.47 eV/atom  total energy: -4.43 eV/atom  temperature: 312.68 K  volume:  1309 Ang^3  density: 2.280 g/cm^3
    steps:   300  energy：-4.49 eV/atom  total energy: -4.43 eV/atom  temperature: 502.01 K  volume:  1307 Ang^3  density: 2.283 g/cm^3
    steps:   400  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 435.63 K  volume:  1313 Ang^3  density: 2.273 g/cm^3
    steps:   500  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 545.67 K  volume:  1314 Ang^3  density: 2.272 g/cm^3
    steps:   600  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 538.09 K  volume:  1315 Ang^3  density: 2.270 g/cm^3
    steps:   700  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 472.07 K  volume:  1311 Ang^3  density: 2.276 g/cm^3
    steps:   800  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 449.84 K  volume:  1313 Ang^3  density: 2.273 g/cm^3
    steps:   900  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 547.05 K  volume:  1314 Ang^3  density: 2.272 g/cm^3
    steps:  1000  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 461.41 K  volume:  1310 Ang^3  density: 2.279 g/cm^3
    steps:  1100  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 419.31 K  volume:  1309 Ang^3  density: 2.281 g/cm^3
    steps:  1200  energy：-4.49 eV/atom  total energy: -4.43 eV/atom  temperature: 527.43 K  volume:  1313 Ang^3  density: 2.274 g/cm^3
    steps:  1300  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 541.05 K  volume:  1315 Ang^3  density: 2.270 g/cm^3
    steps:  1400  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 504.85 K  volume:  1311 Ang^3  density: 2.276 g/cm^3
    steps:  1500  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 494.88 K  volume:  1311 Ang^3  density: 2.277 g/cm^3
    steps:  1600  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 533.33 K  volume:  1314 Ang^3  density: 2.271 g/cm^3
    steps:  1700  energy：-4.50 eV/atom  total energy: -4.42 eV/atom  temperature: 554.94 K  volume:  1313 Ang^3  density: 2.273 g/cm^3
    steps:  1800  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 487.15 K  volume:  1312 Ang^3  density: 2.276 g/cm^3
    steps:  1900  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 448.19 K  volume:  1312 Ang^3  density: 2.274 g/cm^3
    steps:  2000  energy：-4.47 eV/atom  total energy: -4.42 eV/atom  temperature: 413.24 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0600K.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.48 eV/atom  temperature: 560.97 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:   100  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 505.95 K  volume:  1313 Ang^3  density: 2.273 g/cm^3
    steps:   200  energy：-4.47 eV/atom  total energy: -4.41 eV/atom  temperature: 434.24 K  volume:  1309 Ang^3  density: 2.281 g/cm^3
    steps:   300  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 565.37 K  volume:  1311 Ang^3  density: 2.277 g/cm^3
    steps:   400  energy：-4.46 eV/atom  total energy: -4.40 eV/atom  temperature: 464.04 K  volume:  1312 Ang^3  density: 2.275 g/cm^3
    steps:   500  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 659.64 K  volume:  1315 Ang^3  density: 2.269 g/cm^3
    steps:   600  energy：-4.48 eV/atom  total energy: -4.39 eV/atom  temperature: 652.72 K  volume:  1319 Ang^3  density: 2.262 g/cm^3
    steps:   700  energy：-4.47 eV/atom  total energy: -4.39 eV/atom  temperature: 553.70 K  volume:  1315 Ang^3  density: 2.269 g/cm^3
    steps:   800  energy：-4.47 eV/atom  total energy: -4.40 eV/atom  temperature: 573.43 K  volume:  1314 Ang^3  density: 2.271 g/cm^3
    steps:   900  energy：-4.47 eV/atom  total energy: -4.40 eV/atom  temperature: 585.41 K  volume:  1316 Ang^3  density: 2.269 g/cm^3
    steps:  1000  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 620.03 K  volume:  1314 Ang^3  density: 2.272 g/cm^3
    steps:  1100  energy：-4.47 eV/atom  total energy: -4.40 eV/atom  temperature: 585.53 K  volume:  1318 Ang^3  density: 2.264 g/cm^3
    steps:  1200  energy：-4.49 eV/atom  total energy: -4.40 eV/atom  temperature: 716.76 K  volume:  1317 Ang^3  density: 2.266 g/cm^3
    steps:  1300  energy：-4.46 eV/atom  total energy: -4.39 eV/atom  temperature: 540.29 K  volume:  1315 Ang^3  density: 2.269 g/cm^3
    steps:  1400  energy：-4.46 eV/atom  total energy: -4.39 eV/atom  temperature: 517.13 K  volume:  1315 Ang^3  density: 2.270 g/cm^3
    steps:  1500  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 651.68 K  volume:  1315 Ang^3  density: 2.270 g/cm^3
    steps:  1600  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 619.10 K  volume:  1320 Ang^3  density: 2.261 g/cm^3
    steps:  1700  energy：-4.49 eV/atom  total energy: -4.40 eV/atom  temperature: 696.78 K  volume:  1317 Ang^3  density: 2.267 g/cm^3
    steps:  1800  energy：-4.47 eV/atom  total energy: -4.40 eV/atom  temperature: 554.73 K  volume:  1316 Ang^3  density: 2.269 g/cm^3
    steps:  1900  energy：-4.45 eV/atom  total energy: -4.40 eV/atom  temperature: 451.48 K  volume:  1315 Ang^3  density: 2.269 g/cm^3
    steps:  2000  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 609.24 K  volume:  1312 Ang^3  density: 2.274 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0700K.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.46 eV/atom  temperature: 668.43 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:   100  energy：-4.47 eV/atom  total energy: -4.40 eV/atom  temperature: 594.47 K  volume:  1315 Ang^3  density: 2.270 g/cm^3
    steps:   200  energy：-4.45 eV/atom  total energy: -4.38 eV/atom  temperature: 495.82 K  volume:  1309 Ang^3  density: 2.281 g/cm^3
    steps:   300  energy：-4.47 eV/atom  total energy: -4.38 eV/atom  temperature: 717.39 K  volume:  1305 Ang^3  density: 2.286 g/cm^3
    steps:   400  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 581.79 K  volume:  1311 Ang^3  density: 2.277 g/cm^3
    steps:   500  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 714.31 K  volume:  1314 Ang^3  density: 2.271 g/cm^3
    steps:   600  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 751.43 K  volume:  1318 Ang^3  density: 2.264 g/cm^3
    steps:   700  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 744.04 K  volume:  1316 Ang^3  density: 2.268 g/cm^3
    steps:   800  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 775.69 K  volume:  1317 Ang^3  density: 2.266 g/cm^3
    steps:   900  energy：-4.46 eV/atom  total energy: -4.36 eV/atom  temperature: 700.28 K  volume:  1317 Ang^3  density: 2.267 g/cm^3
    steps:  1000  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 635.78 K  volume:  1308 Ang^3  density: 2.281 g/cm^3
    steps:  1100  energy：-4.44 eV/atom  total energy: -4.37 eV/atom  temperature: 526.22 K  volume:  1311 Ang^3  density: 2.277 g/cm^3
    steps:  1200  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 669.63 K  volume:  1315 Ang^3  density: 2.269 g/cm^3
    steps:  1300  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 679.05 K  volume:  1315 Ang^3  density: 2.270 g/cm^3
    steps:  1400  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 632.59 K  volume:  1313 Ang^3  density: 2.273 g/cm^3
    steps:  1500  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 678.20 K  volume:  1319 Ang^3  density: 2.263 g/cm^3
    steps:  1600  energy：-4.48 eV/atom  total energy: -4.37 eV/atom  temperature: 852.16 K  volume:  1318 Ang^3  density: 2.265 g/cm^3
    steps:  1700  energy：-4.46 eV/atom  total energy: -4.36 eV/atom  temperature: 767.38 K  volume:  1313 Ang^3  density: 2.273 g/cm^3
    steps:  1800  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 682.98 K  volume:  1311 Ang^3  density: 2.277 g/cm^3
    steps:  1900  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 595.88 K  volume:  1309 Ang^3  density: 2.280 g/cm^3
    steps:  2000  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 631.18 K  volume:  1310 Ang^3  density: 2.279 g/cm^3


    steps:   500  energy：-4.52 eV/atom  total energy: -4.48 eV/atom  temperature: 334.19 K  volume:  1310 Ang^3  density: 2.279 g/cm^3


    steps:   600  energy：-4.52 eV/atom  total energy: -4.48 eV/atom  temperature: 372.68 K  volume:  1312 Ang^3  density: 2.275 g/cm^3


    steps:   700  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 280.35 K  volume:  1309 Ang^3  density: 2.279 g/cm^3


    steps:   800  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 252.24 K  volume:  1309 Ang^3  density: 2.281 g/cm^3


    steps:   900  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 284.32 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:  1000  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 284.27 K  volume:  1310 Ang^3  density: 2.279 g/cm^3


    steps:  1100  energy：-4.52 eV/atom  total energy: -4.48 eV/atom  temperature: 346.20 K  volume:  1307 Ang^3  density: 2.283 g/cm^3


    steps:  1200  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 284.54 K  volume:  1309 Ang^3  density: 2.279 g/cm^3


    steps:  1300  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 313.12 K  volume:  1311 Ang^3  density: 2.277 g/cm^3


    steps:  1400  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 302.57 K  volume:  1309 Ang^3  density: 2.280 g/cm^3


    steps:  1500  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 296.08 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:  1600  energy：-4.52 eV/atom  total energy: -4.48 eV/atom  temperature: 334.50 K  volume:  1310 Ang^3  density: 2.279 g/cm^3


    steps:  1700  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 281.91 K  volume:  1309 Ang^3  density: 2.281 g/cm^3


    steps:  1800  energy：-4.51 eV/atom  total energy: -4.48 eV/atom  temperature: 273.04 K  volume:  1309 Ang^3  density: 2.280 g/cm^3


    steps:  1900  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 283.92 K  volume:  1310 Ang^3  density: 2.279 g/cm^3


    steps:  2000  energy：-4.51 eV/atom  total energy: -4.47 eV/atom  temperature: 323.08 K  volume:  1308 Ang^3  density: 2.281 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0400K.traj.


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-4.55 eV/atom  total energy: -4.50 eV/atom  temperature: 422.40 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:   100  energy：-4.51 eV/atom  total energy: -4.46 eV/atom  temperature: 354.15 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:   200  energy：-4.49 eV/atom  total energy: -4.45 eV/atom  temperature: 254.77 K  volume:  1307 Ang^3  density: 2.284 g/cm^3


    steps:   300  energy：-4.49 eV/atom  total energy: -4.45 eV/atom  temperature: 310.38 K  volume:  1306 Ang^3  density: 2.285 g/cm^3


    steps:   400  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 364.65 K  volume:  1308 Ang^3  density: 2.282 g/cm^3


    steps:   500  energy：-4.51 eV/atom  total energy: -4.45 eV/atom  temperature: 468.88 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:   600  energy：-4.51 eV/atom  total energy: -4.45 eV/atom  temperature: 459.47 K  volume:  1313 Ang^3  density: 2.274 g/cm^3


    steps:   700  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 376.19 K  volume:  1311 Ang^3  density: 2.276 g/cm^3


    steps:   800  energy：-4.49 eV/atom  total energy: -4.45 eV/atom  temperature: 341.04 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:   900  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 410.58 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:  1000  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 376.01 K  volume:  1310 Ang^3  density: 2.279 g/cm^3


    steps:  1100  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 411.15 K  volume:  1309 Ang^3  density: 2.280 g/cm^3


    steps:  1200  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 372.38 K  volume:  1309 Ang^3  density: 2.280 g/cm^3


    steps:  1300  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 400.88 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:  1400  energy：-4.51 eV/atom  total energy: -4.45 eV/atom  temperature: 454.37 K  volume:  1309 Ang^3  density: 2.280 g/cm^3


    steps:  1500  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 380.64 K  volume:  1311 Ang^3  density: 2.277 g/cm^3


    steps:  1600  energy：-4.51 eV/atom  total energy: -4.45 eV/atom  temperature: 469.78 K  volume:  1312 Ang^3  density: 2.275 g/cm^3


    steps:  1700  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 358.30 K  volume:  1311 Ang^3  density: 2.276 g/cm^3


    steps:  1800  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 387.66 K  volume:  1309 Ang^3  density: 2.281 g/cm^3


    steps:  1900  energy：-4.50 eV/atom  total energy: -4.45 eV/atom  temperature: 410.54 K  volume:  1308 Ang^3  density: 2.282 g/cm^3


    steps:  2000  energy：-4.49 eV/atom  total energy: -4.45 eV/atom  temperature: 331.71 K  volume:  1309 Ang^3  density: 2.279 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0500K.traj.


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-4.55 eV/atom  total energy: -4.49 eV/atom  temperature: 492.33 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:   100  energy：-4.50 eV/atom  total energy: -4.44 eV/atom  temperature: 454.70 K  volume:  1312 Ang^3  density: 2.275 g/cm^3


    steps:   200  energy：-4.48 eV/atom  total energy: -4.43 eV/atom  temperature: 344.98 K  volume:  1307 Ang^3  density: 2.284 g/cm^3


    steps:   300  energy：-4.48 eV/atom  total energy: -4.43 eV/atom  temperature: 419.99 K  volume:  1305 Ang^3  density: 2.287 g/cm^3


    steps:   400  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 419.10 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:   500  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 551.87 K  volume:  1311 Ang^3  density: 2.277 g/cm^3


    steps:   600  energy：-4.50 eV/atom  total energy: -4.42 eV/atom  temperature: 583.28 K  volume:  1313 Ang^3  density: 2.272 g/cm^3


    steps:   700  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 526.09 K  volume:  1313 Ang^3  density: 2.273 g/cm^3


    steps:   800  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 454.31 K  volume:  1313 Ang^3  density: 2.274 g/cm^3


    steps:   900  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 510.94 K  volume:  1313 Ang^3  density: 2.274 g/cm^3


    steps:  1000  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 482.36 K  volume:  1309 Ang^3  density: 2.280 g/cm^3


    steps:  1100  energy：-4.47 eV/atom  total energy: -4.42 eV/atom  temperature: 422.04 K  volume:  1308 Ang^3  density: 2.281 g/cm^3


    steps:  1200  energy：-4.49 eV/atom  total energy: -4.43 eV/atom  temperature: 459.62 K  volume:  1311 Ang^3  density: 2.276 g/cm^3


    steps:  1300  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 545.17 K  volume:  1315 Ang^3  density: 2.271 g/cm^3


    steps:  1400  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 501.13 K  volume:  1313 Ang^3  density: 2.273 g/cm^3


    steps:  1500  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 486.21 K  volume:  1313 Ang^3  density: 2.274 g/cm^3


    steps:  1600  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 537.26 K  volume:  1313 Ang^3  density: 2.273 g/cm^3


    steps:  1700  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 509.51 K  volume:  1314 Ang^3  density: 2.272 g/cm^3


    steps:  1800  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 485.59 K  volume:  1312 Ang^3  density: 2.276 g/cm^3


    steps:  1900  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 513.69 K  volume:  1309 Ang^3  density: 2.280 g/cm^3


    steps:  2000  energy：-4.48 eV/atom  total energy: -4.42 eV/atom  temperature: 455.48 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0600K.traj.


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-4.55 eV/atom  total energy: -4.48 eV/atom  temperature: 559.45 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:   100  energy：-4.49 eV/atom  total energy: -4.42 eV/atom  temperature: 498.64 K  volume:  1316 Ang^3  density: 2.269 g/cm^3


    steps:   200  energy：-4.46 eV/atom  total energy: -4.40 eV/atom  temperature: 466.24 K  volume:  1313 Ang^3  density: 2.272 g/cm^3


    steps:   300  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 612.91 K  volume:  1311 Ang^3  density: 2.277 g/cm^3


    steps:   400  energy：-4.46 eV/atom  total energy: -4.40 eV/atom  temperature: 497.54 K  volume:  1315 Ang^3  density: 2.270 g/cm^3


    steps:   500  energy：-4.49 eV/atom  total energy: -4.40 eV/atom  temperature: 719.54 K  volume:  1317 Ang^3  density: 2.267 g/cm^3


    steps:   600  energy：-4.49 eV/atom  total energy: -4.40 eV/atom  temperature: 692.29 K  volume:  1321 Ang^3  density: 2.260 g/cm^3


    steps:   700  energy：-4.47 eV/atom  total energy: -4.39 eV/atom  temperature: 604.49 K  volume:  1317 Ang^3  density: 2.266 g/cm^3


    steps:   800  energy：-4.47 eV/atom  total energy: -4.40 eV/atom  temperature: 553.93 K  volume:  1317 Ang^3  density: 2.266 g/cm^3


    steps:   900  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 615.88 K  volume:  1316 Ang^3  density: 2.268 g/cm^3


    steps:  1000  energy：-4.47 eV/atom  total energy: -4.39 eV/atom  temperature: 565.43 K  volume:  1315 Ang^3  density: 2.269 g/cm^3


    steps:  1100  energy：-4.46 eV/atom  total energy: -4.40 eV/atom  temperature: 464.16 K  volume:  1319 Ang^3  density: 2.262 g/cm^3


    steps:  1200  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 660.29 K  volume:  1318 Ang^3  density: 2.264 g/cm^3


    steps:  1300  energy：-4.47 eV/atom  total energy: -4.39 eV/atom  temperature: 562.11 K  volume:  1318 Ang^3  density: 2.264 g/cm^3


    steps:  1400  energy：-4.47 eV/atom  total energy: -4.39 eV/atom  temperature: 581.03 K  volume:  1316 Ang^3  density: 2.268 g/cm^3


    steps:  1500  energy：-4.48 eV/atom  total energy: -4.40 eV/atom  temperature: 651.96 K  volume:  1317 Ang^3  density: 2.266 g/cm^3


    steps:  1600  energy：-4.47 eV/atom  total energy: -4.39 eV/atom  temperature: 583.42 K  volume:  1319 Ang^3  density: 2.262 g/cm^3


    steps:  1700  energy：-4.48 eV/atom  total energy: -4.39 eV/atom  temperature: 694.34 K  volume:  1315 Ang^3  density: 2.270 g/cm^3


    steps:  1800  energy：-4.47 eV/atom  total energy: -4.40 eV/atom  temperature: 567.03 K  volume:  1316 Ang^3  density: 2.267 g/cm^3


    steps:  1900  energy：-4.47 eV/atom  total energy: -4.39 eV/atom  temperature: 570.64 K  volume:  1315 Ang^3  density: 2.270 g/cm^3


    steps:  2000  energy：-4.47 eV/atom  total energy: -4.40 eV/atom  temperature: 572.76 K  volume:  1313 Ang^3  density: 2.274 g/cm^3


      0%|          | 0/2000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/output/md_specific_heat_0700K.traj.


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-4.55 eV/atom  total energy: -4.46 eV/atom  temperature: 706.41 K  volume:  1281 Ang^3  density: 2.330 g/cm^3




    steps:   100  energy：-4.47 eV/atom  total energy: -4.39 eV/atom  temperature: 588.76 K  volume:  1314 Ang^3  density: 2.271 g/cm^3


    steps:   200  energy：-4.43 eV/atom  total energy: -4.38 eV/atom  temperature: 417.19 K  volume:  1308 Ang^3  density: 2.282 g/cm^3


    steps:   300  energy：-4.46 eV/atom  total energy: -4.38 eV/atom  temperature: 666.89 K  volume:  1305 Ang^3  density: 2.288 g/cm^3


    steps:   400  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 567.00 K  volume:  1311 Ang^3  density: 2.277 g/cm^3


    steps:   500  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 782.66 K  volume:  1315 Ang^3  density: 2.270 g/cm^3


    steps:   600  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 772.01 K  volume:  1318 Ang^3  density: 2.264 g/cm^3


    steps:   700  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 683.28 K  volume:  1315 Ang^3  density: 2.269 g/cm^3


    steps:   800  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 746.52 K  volume:  1314 Ang^3  density: 2.272 g/cm^3


    steps:   900  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 769.91 K  volume:  1314 Ang^3  density: 2.271 g/cm^3


    steps:  1000  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 628.23 K  volume:  1314 Ang^3  density: 2.272 g/cm^3


    steps:  1100  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 653.46 K  volume:  1312 Ang^3  density: 2.275 g/cm^3


    steps:  1200  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 746.56 K  volume:  1315 Ang^3  density: 2.270 g/cm^3


    steps:  1300  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 694.03 K  volume:  1315 Ang^3  density: 2.269 g/cm^3


    steps:  1400  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 657.43 K  volume:  1310 Ang^3  density: 2.278 g/cm^3


    steps:  1500  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 663.93 K  volume:  1317 Ang^3  density: 2.267 g/cm^3


    steps:  1600  energy：-4.47 eV/atom  total energy: -4.37 eV/atom  temperature: 827.25 K  volume:  1317 Ang^3  density: 2.267 g/cm^3


    steps:  1700  energy：-4.47 eV/atom  total energy: -4.36 eV/atom  temperature: 780.73 K  volume:  1312 Ang^3  density: 2.275 g/cm^3


    steps:  1800  energy：-4.45 eV/atom  total energy: -4.37 eV/atom  temperature: 647.19 K  volume:  1313 Ang^3  density: 2.273 g/cm^3


    steps:  1900  energy：-4.44 eV/atom  total energy: -4.37 eV/atom  temperature: 565.12 K  volume:  1311 Ang^3  density: 2.277 g/cm^3


    steps:  2000  energy：-4.46 eV/atom  total energy: -4.37 eV/atom  temperature: 705.61 K  volume:  1312 Ang^3  density: 2.275 g/cm^3


## 比熱の計算
* PostMDSpecificHeatFeature を利用して、以下のように比熱の計算とプロットの生成を行います。
* 初回起動時などで図が表示されない場合、ブラウザでページをリロードすると表示されることがあります。


```python
specific_heat_feature = PostMDSpecificHeatFeature()
specific_heat = specific_heat_feature(res, [300, 350, 400, 450, 500, 550, 600], init_time=400)
specific_heat_fig = specific_heat.plot()
iplot(specific_heat_fig)
```


## 実験値との比較
NIST Chemistry WebBook (https://webbook.nist.gov/chemistry) による Si 結晶の比熱の実験値は以下となります。

| 温度 [K] | C&lt;sub&gt;p&lt;/sub&gt; [J/(mol K)] |
|-----------------|----------------|
|300              |20.05           |
|400              |22.15           |
|500              |23.34           |
|600              |24.15           |

以下のコードを実行することで Matlantis による計算値と実験値の比較を可視化するプロットを生成できます。
初回起動時などで図が表示されない場合、ブラウザでページをリロードすると表示されることがあります。


```python
T = [300., 400., 500., 600.]
Cp = [20.05, 22.15, 23.34, 24.15]

specific_heat_fig.add_trace(
    go.Scatter(
        x=T, y=[i/28 for i in Cp],
        name="expt"
    ),
    row=1,
    col=2,
)
iplot(specific_heat_fig)
```


## 結論
1. 分子動力学計算に基づいている Matlantis-features の PostMDSpecificHeatFeature を利用して、 Si 結晶の比熱を計算しました。
2. Matlantis-features によって得られた比熱の計算値（0.9 J/K/g）が、Dulong–Petit の法則による推定値（0.89 J/K/g）と近い値をとることがわかりました。
3. 低温での比熱の計算値は実験値との誤差が比較的大きいことがわかりました。原因としては、分子動力学計算には量子効果が考慮されていないことが考えられます。
