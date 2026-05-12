# Matlantis-features: 粘性係数 (EMD)


## 概要

* 粘性係数は、主に流体におけるねばりの強さを表す量です。

* Matlantis-featuresでは、分子動力学計算の軌跡に対して粘性係数を計算するフィーチャーとして `PostEMDViscosityFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/viscosity_emd.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;) を提供しています。このフィーチャーでは、拡散係数と同様に、分子動力学計算の軌跡を処理することにより求めます。

* 別途、reverse non-equilibrium MD (rNEMD)という手法により粘性係数を求める手法があります。そちらは粘性係数(rNEMD)のExample (`viscosity_rnemd.ipynb`) を参照してください。


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
from matlantis_features.features.md import PostEMDViscosityFeature

try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

```


## 材料の用意

* ここでは、有機分子のデカンに対する粘性係数を計算してみます。


```python
atoms_n_decane = MatlantisAtoms.from_file(str(dir_path/"assets/viscosity/n-decane.xyz"))
atoms = atoms_n_decane.ase_atoms
```


## 分子動力学計算

* 特定の温度での分子動力学シミュレーションを行い、後の粘性係数の計算に使う軌跡を計算しています。シミュレーションには一定の時間がかかります。

* 以下は、NVT(Berendsen)アンサンブルにおける温度一定のMD計算を行う例です。


```python
integ = NVTBerendsenIntegrator(timestep=1.0, temperature=480.0)
mdsys = ASEMDSystem(atoms)
mdsys.init_temperature(480.0)
md = MDFeature(
    integrator=integ,
    n_run=20000,
    traj_file_name="md_viscosity_emd.traj",
    traj_freq=100,
    show_progress_bar=True,
    show_logger=True,
    logger_interval=1000,
    estimator_fn=estimator_fn,
)
md_res = md(mdsys)

clear_output(wait=True)

```


      0%|          | 0/20000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/md_viscosity_emd.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 490.72 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  1000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 494.19 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  2000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 475.15 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  3000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 482.53 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  4000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 488.64 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  5000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 507.27 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  6000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 479.75 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  7000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 470.62 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  8000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 496.35 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps:  9000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 486.41 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 10000  energy：-4.13 eV/atom  total energy: -4.08 eV/atom  temperature: 451.06 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 11000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 496.97 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 12000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 467.51 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 13000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 478.58 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 14000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 514.92 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 15000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 507.79 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 16000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 491.48 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 17000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 480.52 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 18000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 509.46 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 19000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 468.49 K  volume:  3851 Ang^3  density: 0.614 g/cm^3
    steps: 20000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 476.74 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 498.78 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  1000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 459.12 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  2000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 480.36 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  3000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 482.61 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  4000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 437.86 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  5000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 454.85 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  6000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 472.98 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  7000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 476.44 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  8000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 474.96 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps:  9000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 470.53 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 10000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 486.59 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 11000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 476.11 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 12000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 481.05 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 13000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 504.55 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 14000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 466.72 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 15000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 473.95 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 16000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 490.12 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 17000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 469.64 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 18000  energy：-4.13 eV/atom  total energy: -4.07 eV/atom  temperature: 448.90 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 19000  energy：-4.14 eV/atom  total energy: -4.08 eV/atom  temperature: 475.39 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


    steps: 20000  energy：-4.14 eV/atom  total energy: -4.07 eV/atom  temperature: 493.19 K  volume:  3851 Ang^3  density: 0.614 g/cm^3


## 粘性係数の計算

* MD計算の軌跡から、粘性係数を計算します。


```python
vis = PostEMDViscosityFeature()
vis_res = vis(md_res, init_time=1000, number_of_segments=4)
print(vis_res.viscosity)
```

    {'amu/A/fs': [0.04364935211339314], 'Pa s': [0.0007248145325499582], 'mPa s': [0.7248145325499581]}


## 結果のプロット

* n-decaneの粘性係数としては 0.178 mPa s という値が推定されています。 (S. T. Cui et al., Molecular Physics, 93 (1998) 117.)

* 初回起動時などで図が表示されない場合、ブラウザでページをリロードすると表示されることがあります。


```python
fig = vis_res.plot()
fig.update_layout(width=1000, height=400)
iplot(fig)
```


