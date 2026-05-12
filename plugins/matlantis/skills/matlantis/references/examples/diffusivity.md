# Matlantis-features: 拡散係数


## 概要

* 拡散係数は、主に流体などにおける原子の拡散の速さに相当する量です。分子動力学計算では通常、原子の平均二乗変位から計算されます。

* Matlantis-featuresでは、分子動力学計算の軌跡に対して拡散係数を計算するフィーチャーとして `PostMDDiffusionFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/diffusion.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;) を提供しています。

* このExampleは、分子動力学(MD)計算の例にもなっています。以下ではMatlantis-featuresのMD計算ツールを利用した例をあわせて紹介します。


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
import pandas as pd
import plotly.graph_objects as go
from ase.io import read
from plotly.offline import iplot
from IPython.display import clear_output

from matlantis_features.atoms import MatlantisAtoms
from ase import units
from matlantis_features.features.md import NPTIntegrator
from matlantis_features.features.md import ASEMDSystem
from matlantis_features.features.md import MDFeature
from matlantis_features.features.md import PostMDDiffusionFeature, OldPostMDDiffusionFeature

try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

```


## 材料の用意

* ここでは、液体のアルミニウムを例として拡散係数を計算してみます。


```python
ase_atoms_Al = read(str(dir_path/"assets/diffusivity/liquid_Al.xyz"))
atoms = MatlantisAtoms.from_ase_atoms(ase_atoms_Al)
```


## 分子動力学計算

* 特定の温度での分子動力学シミュレーションを行い、後の拡散係数の計算に使う軌跡を計算しています。分子動力学計算には少なくとも数十分程度の時間がかかります。

* 以下は、NPTアンサンブルにおける温度一定のMD計算を行う例です。Nose-Hoover thermostat を使用して温度を制御し、Parrinello-Rahman barostat を使用して圧力を制御します。


```python
mdsys = ASEMDSystem(atoms.ase_atoms.copy())
mdsys.init_temperature(1020)

integrator = NPTIntegrator(
    timestep=1.0,
    temperature=1020,
    pressure=101325 * units.Pascal,
    ttime=100,
    pfactor=100 * units.fs,
    mask=np.eye(3)
)

md = MDFeature(
    integrator=integrator,
    n_run=20000,
    traj_file_name="md_diffusion.traj",
    traj_freq=100,
    show_progress_bar=True,
    show_logger=True,
    logger_interval=1000,
    estimator_fn=estimator_fn,
)


md_results = md(mdsys)

clear_output(wait=True)
```


      0%|          | 0/20000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/md_diffusion.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    WARNING: NPT: Setting the center-of-mass momentum to zero (was -1.2828 10.6755 4.39525)
    steps:     0  energy：-3.17 eV/atom  total energy: -3.05 eV/atom  temperature: 914.08 K  volume:   611 Ang^3  density: 2.345 g/cm^3
    steps:     0  energy：-3.17 eV/atom  total energy: -3.05 eV/atom  temperature: 914.08 K  volume:   611 Ang^3  density: 2.345 g/cm^3
    steps:  1000  energy：-3.17 eV/atom  total energy: -3.07 eV/atom  temperature: 777.05 K  volume:   617 Ang^3  density: 2.325 g/cm^3
    steps:  2000  energy：-3.21 eV/atom  total energy: -3.07 eV/atom  temperature: 1135.39 K  volume:   602 Ang^3  density: 2.380 g/cm^3
    steps:  3000  energy：-3.16 eV/atom  total energy: -3.01 eV/atom  temperature: 1159.69 K  volume:   612 Ang^3  density: 2.342 g/cm^3
    steps:  4000  energy：-3.20 eV/atom  total energy: -3.07 eV/atom  temperature: 971.40 K  volume:   594 Ang^3  density: 2.413 g/cm^3
    steps:  5000  energy：-3.18 eV/atom  total energy: -3.04 eV/atom  temperature: 1046.83 K  volume:   597 Ang^3  density: 2.402 g/cm^3
    steps:  6000  energy：-3.20 eV/atom  total energy: -3.08 eV/atom  temperature: 973.97 K  volume:   600 Ang^3  density: 2.391 g/cm^3
    steps:  7000  energy：-3.22 eV/atom  total energy: -3.10 eV/atom  temperature: 914.02 K  volume:   598 Ang^3  density: 2.398 g/cm^3
    steps:  8000  energy：-3.17 eV/atom  total energy: -3.05 eV/atom  temperature: 930.18 K  volume:   621 Ang^3  density: 2.309 g/cm^3
    steps:  9000  energy：-3.15 eV/atom  total energy: -3.04 eV/atom  temperature: 856.97 K  volume:   601 Ang^3  density: 2.385 g/cm^3
    steps: 10000  energy：-3.22 eV/atom  total energy: -3.08 eV/atom  temperature: 1039.03 K  volume:   609 Ang^3  density: 2.352 g/cm^3
    steps: 11000  energy：-3.17 eV/atom  total energy: -3.05 eV/atom  temperature: 943.10 K  volume:   636 Ang^3  density: 2.256 g/cm^3
    steps: 12000  energy：-3.15 eV/atom  total energy: -3.04 eV/atom  temperature: 849.87 K  volume:   625 Ang^3  density: 2.293 g/cm^3
    steps: 13000  energy：-3.20 eV/atom  total energy: -3.07 eV/atom  temperature: 980.83 K  volume:   606 Ang^3  density: 2.365 g/cm^3
    steps: 14000  energy：-3.21 eV/atom  total energy: -3.06 eV/atom  temperature: 1176.16 K  volume:   600 Ang^3  density: 2.392 g/cm^3
    steps: 15000  energy：-3.18 eV/atom  total energy: -3.05 eV/atom  temperature: 1016.14 K  volume:   614 Ang^3  density: 2.335 g/cm^3
    steps: 16000  energy：-3.19 eV/atom  total energy: -3.09 eV/atom  temperature: 837.43 K  volume:   609 Ang^3  density: 2.353 g/cm^3
    steps: 17000  energy：-3.20 eV/atom  total energy: -3.07 eV/atom  temperature: 981.60 K  volume:   619 Ang^3  density: 2.316 g/cm^3
    steps: 18000  energy：-3.16 eV/atom  total energy: -3.04 eV/atom  temperature: 917.83 K  volume:   626 Ang^3  density: 2.292 g/cm^3
    steps: 19000  energy：-3.15 eV/atom  total energy: -3.03 eV/atom  temperature: 940.76 K  volume:   630 Ang^3  density: 2.275 g/cm^3
    steps: 20000  energy：-3.19 eV/atom  total energy: -3.08 eV/atom  temperature: 832.68 K  volume:   594 Ang^3  density: 2.415 g/cm^3


    Note: The max disk size of /home/jovyan is about 99G.


    WARNING: NPT: Setting the center-of-mass momentum to zero (was 3.32533 -10.7961 -7.1486)


    steps:     0  energy：-3.17 eV/atom  total energy: -3.03 eV/atom  temperature: 1018.75 K  volume:   611 Ang^3  density: 2.345 g/cm^3


    steps:     0  energy：-3.17 eV/atom  total energy: -3.03 eV/atom  temperature: 1018.75 K  volume:   611 Ang^3  density: 2.345 g/cm^3


    steps:  1000  energy：-3.14 eV/atom  total energy: -3.02 eV/atom  temperature: 972.19 K  volume:   646 Ang^3  density: 2.221 g/cm^3


    steps:  2000  energy：-3.20 eV/atom  total energy: -3.08 eV/atom  temperature: 949.68 K  volume:   631 Ang^3  density: 2.273 g/cm^3


    steps:  3000  energy：-3.15 eV/atom  total energy: -2.99 eV/atom  temperature: 1219.74 K  volume:   611 Ang^3  density: 2.348 g/cm^3


    steps:  4000  energy：-3.18 eV/atom  total energy: -3.05 eV/atom  temperature: 1032.28 K  volume:   589 Ang^3  density: 2.435 g/cm^3


    steps:  5000  energy：-3.14 eV/atom  total energy: -3.03 eV/atom  temperature: 845.19 K  volume:   623 Ang^3  density: 2.302 g/cm^3


    steps:  6000  energy：-3.17 eV/atom  total energy: -3.06 eV/atom  temperature: 908.51 K  volume:   606 Ang^3  density: 2.367 g/cm^3


    steps:  7000  energy：-3.17 eV/atom  total energy: -3.04 eV/atom  temperature: 1030.17 K  volume:   634 Ang^3  density: 2.262 g/cm^3


    steps:  8000  energy：-3.15 eV/atom  total energy: -3.03 eV/atom  temperature: 865.76 K  volume:   643 Ang^3  density: 2.228 g/cm^3


    steps:  9000  energy：-3.14 eV/atom  total energy: -3.00 eV/atom  temperature: 1044.19 K  volume:   656 Ang^3  density: 2.184 g/cm^3


    steps: 10000  energy：-3.20 eV/atom  total energy: -3.09 eV/atom  temperature: 830.77 K  volume:   614 Ang^3  density: 2.336 g/cm^3


    steps: 11000  energy：-3.21 eV/atom  total energy: -3.07 eV/atom  temperature: 1046.73 K  volume:   596 Ang^3  density: 2.405 g/cm^3


    steps: 12000  energy：-3.19 eV/atom  total energy: -3.08 eV/atom  temperature: 882.14 K  volume:   609 Ang^3  density: 2.353 g/cm^3


    steps: 13000  energy：-3.19 eV/atom  total energy: -3.06 eV/atom  temperature: 989.15 K  volume:   599 Ang^3  density: 2.393 g/cm^3


    steps: 14000  energy：-3.18 eV/atom  total energy: -3.05 eV/atom  temperature: 1076.19 K  volume:   620 Ang^3  density: 2.311 g/cm^3


    steps: 15000  energy：-3.22 eV/atom  total energy: -3.10 eV/atom  temperature: 885.19 K  volume:   608 Ang^3  density: 2.357 g/cm^3


    steps: 16000  energy：-3.18 eV/atom  total energy: -3.05 eV/atom  temperature: 989.96 K  volume:   597 Ang^3  density: 2.401 g/cm^3


    steps: 17000  energy：-3.22 eV/atom  total energy: -3.09 eV/atom  temperature: 994.41 K  volume:   604 Ang^3  density: 2.376 g/cm^3


    steps: 18000  energy：-3.18 eV/atom  total energy: -3.05 eV/atom  temperature: 1018.40 K  volume:   611 Ang^3  density: 2.348 g/cm^3


    steps: 19000  energy：-3.14 eV/atom  total energy: -3.04 eV/atom  temperature: 777.14 K  volume:   626 Ang^3  density: 2.292 g/cm^3


    steps: 20000  energy：-3.23 eV/atom  total energy: -3.08 eV/atom  temperature: 1138.76 K  volume:   597 Ang^3  density: 2.403 g/cm^3


## 拡散係数の計算

* MD計算の軌跡から、拡散係数を計算します。
* MDの開始時に系が平衡状態にないと想定される場合、`init_time` パラメータを指定することでMDの初期のトラジェクトリを破棄します。
* `direction` パラメータは、原子の拡散係数を x, y, z 以外の特定の方向に沿って計算する際に指定します。
* 拡散係数の計算手法には `normal` と `segment` の2つの方法があります。一般的には `segment` 法がより正確ですが、計算時間が長くなります。
* `effective_msd_range` パラメータにより、平均二乗変位のどの部分(MDステップの範囲)を用いて拡散係数を計算するかを指定できます。


```python
from matlantis_features.features.md import MDFeatureResult

md_results = MDFeatureResult.from_traj_obj("md_diffusion.traj")

diffusion = PostMDDiffusionFeature()
diffusion_results = diffusion(
    md_results, init_time=1000.0, stride=10,
    number_of_segments=1, direction=np.array([[1,1,1]]),
    method="segment", effective_msd_range=(0.1, 0.9),
)
```


* Al 原子の拡散係数は `PostMDDiffusionFeatureResult` に保存されます。加えて、x, y, z 軸方向およびユーザが指定した方向 (i.e. [1,1,1]) に沿った拡散係数も提供されます。


```python
diffusion_results.diffusion_coefficient["cm^2/s"]
```


    {'Al': np.float64(5.859912127029276e-05),
     'Al_x': np.float64(8.197017043291903e-05),
     'Al_y': np.float64(6.0334088517563805e-05),
     'Al_z': np.float64(3.3493104860395424e-05),
     'Al_direc_0': np.float64(6.391028492696725e-05)}


## 平均二乗変位のプロット

* 初回起動時などで図が表示されない場合、ブラウザでページをリロードすると表示されることがあります。
* 図の上で平均二乗変位の変化が直線的になっているかどうかにご注意ください。基本的に、平均二乗変位の変化は直線的になるべきですが、そうでない場合は拡散係数の算出に問題が生じることがあります。原因として、(1)MDシミュレーションが短すぎる、(2)拡散係数が低すぎて原子の長距離変位が観測されない、などの理由が考えられます。
* 平均二乗変位は `PostMDDiffusionFeatureResult` の `plot` 関数でプロットすることができます。図では、フィッティングされた直線が破線、拡散係数のフィッティングに用いた有効範囲が灰色の長方形で示されます。


```python
fig = diffusion_results.plot(show_effective_range=True, show_fit_line=True)
fig.update_layout(title="MSD of liquid Al @ 1020 K")
fig.update_layout(width=800, height=500)
iplot(fig)
```


## 実験値との比較
* 液体アルミニウムの自己拡散係数の実験値は、incoherent quasielastic neutron scatteringによる測定により1020 Kで 7.9x10-5 cm^2/s という値が得られています。 (1. Kargl, F., Weis, H., Unruh, T. &amp; Meyer, A. Self diffusion in liquid aluminium. J. Phys. Conf. Ser. 340, 6–11 (2012).)


```python
print("Matlantis: {:.1e} cm^2/s".format(diffusion_results.dict["diffusion_coefficient"]["cm^2/s"]["Al"]))
print("Experiment: {:.1e} cm^2/s".format(7.9e-5))
```

    Matlantis: 5.9e-05 cm^2/s
    Experiment: 7.9e-05 cm^2/s


## 補足: 以前のバージョンの PostMDDiffusionFeature ("matlantis_features&lt;0.8.0") の利用

`PostMDDiffusionFeature` にはバージョン 0.8.0 にて破壊的変更が加えられ、以前のバージョンの `PostMDDiffusionFeature` は `OldPostMDDiffusionFeature` にリネームされています。`OldPostMDDiffusionFeature` の利用方法を以下に示します。


```python
diffusion_ase = OldPostMDDiffusionFeature()
diffusion_ase_results = diffusion(
    md_results, init_time=0.0, stride=10,
    number_of_segments=1
)
```

    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/md/post_features/diffusion.py:548: FutureWarning:

    The OldPostMDDiffusionFeature will be deprecated in the future. Please use PostMDDiffusionFeature instead.



```python
fig = diffusion_ase_results.plot()
fig.update_layout(title="MSD of liquid Al @ 1020 K")
fig.update_layout(width=800, height=500)
iplot(fig)
```


```python
print("Matlantis: {:.1e}".format(diffusion_ase_results.dict["diffusion_coefficient"]["cm^2/s"]["Al"]))
print("Experiment: {:.1e}".format(7.9e-5))
```

    Matlantis: 5.1e-05
    Experiment: 7.9e-05

