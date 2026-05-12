# Simulation Visualization

このExampleでは、物質シミュレーションでよく用いられる計算結果の可視化の例を紹介します。こちらのコードをスニペットとして用いることで図の作成を簡単に行うことができます。


## 構造最適化

* [Structural Optimization](#opt_1_jp)

&lt;img src="assets/visualization/opt_1.png" width=25%&gt;

&lt;img src="assets/visualization/opt_2.png" width=25%&gt;


## 分子動力学法(MD)

* [MD 基本情報](#md_2_jp)

&lt;img src="assets/visualization/md_plt.png" width=25%&gt;

&lt;img src="assets/visualization/md_plotly.png" width=25%&gt;


* [MD 平均二乗変位](#md_4_jp)

&lt;img src="assets/visualization/msd.png" width=25%&gt;

* [MD 平均二乗変位 (matlantis_features)](#md_5_jp)

&lt;img src="assets/visualization/diffusion_mf.png" width=25%&gt;

* [MD 応力 ACF と粘性係数 (matlantis_features)](#md_6_jp)

&lt;img src="assets/visualization/viscosity_mf.png" width=25%&gt;


## フォノン

* [フォノン 分散と状態密度](#phonon_1_jp)

&lt;img src="assets/visualization/phonon.png" width=25%&gt;


* [フォノン 熱力学特性](#phonon_2_jp)

&lt;img src="assets/visualization/phonon_thermo.png" width=25%&gt;


## Nudged elastic band

* [反応経路 エネルギープロファイル](#neb_1_jp)

&lt;img src="assets/visualization/neb.png" width=25% &gt;


# ソースコード


&lt;a id='opt_1_jp'&gt;&lt;/a&gt;
## 構造最適化
このexampleでは、構造最適化中のエネルギーと力の変化をプロットします。まず、構造最適化を実行し、エネルギーと力の変化をトラジェクトリから抽出します。


```python
%pip install "matlantis_features&gt;=0.14.1" "pfp_api_client&gt;=1.7.0"
```



```python
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

import numpy as np
from ase.io import Trajectory
from ase.build import molecule
from ase.optimize import LBFGS

estimator = Estimator(calc_mode=EstimatorCalcMode.PBE_U)
calc = ASECalculator(estimator)

atoms = molecule("CH3CH2OH")
atoms.calc = calc

opt = LBFGS(atoms, logfile=None, trajectory="opt.traj")
opt.run(fmax=0.001)

energy_profile = [a.get_potential_energy()/len(a) for a in Trajectory("opt.traj")]
forces_profile = [np.linalg.norm(a.get_forces(), axis=1).max() for a in Trajectory("opt.traj")]
```


`matplotlib` を使ってエネルギーと力をプロットします。


```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,2, figsize=(10,4))

ax[0].plot(energy_profile)
ax[0].set_xlabel("steps")
ax[0].set_ylabel("energy (eV/atom)")

ax[1].plot(forces_profile)
ax[1].set_xlabel("steps")
ax[1].set_ylabel("max force (eV/A)")
```


    Text(0, 0.5, 'max force (eV/A)')



![png](output_10_1.png)



`plotly` を使ってエネルギーと力をプロットします。


```python
import plotly.subplots
import plotly.graph_objects as go


steps = np.arange(len(energy_profile))

fig = plotly.subplots.make_subplots(1, 2)
fig.add_trace(go.Scatter(x=steps, y=energy_profile, mode="lines", showlegend=False), row=1, col=1)
fig.add_trace(go.Scatter(x=steps, y=forces_profile, mode="lines", showlegend=False), row=1, col=2)

fig.update_xaxes(title="steps", row=1, col=1)
fig.update_yaxes(title="energy (eV/atom)", row=1, col=1)
fig.update_xaxes(title="steps", row=1, col=2)
fig.update_yaxes(title="max force (eV/atom)", row=1, col=2)

fig.update_layout(height=500, width=1000)
```


&lt;a id='md_2_jp'&gt;&lt;/a&gt;
### MD 基本情報
このexampleでは、MDシミュレーションにおける系のポテンシャルエネルギー、温度、体積、エンタルピーなどの基本的な情報をプロットします。

MDのトラジェクトリから各種情報を取得します。ここでは、事前に実行済みのMDのトラジェクトリを使用します。実際のMDの実行方法やトラジェクトリの保存方法については `Welcome` tutorialをご覧ください。


```python
from ase import units
from ase.io import Trajectory

traj = Trajectory("assets/visualization/md_example.traj")

step = list(range(len(traj)))
E = [a.get_potential_energy() / len(a) for a in traj]
T = [a.get_temperature() for a in traj]
V = [a.get_volume() / len(a) for a in traj]
H = [
    (a.get_potential_energy() + a.get_stress(voigt=False, include_ideal_gas=True).trace()/3.0 * a.get_volume())/len(a)
    for a in traj
]
```


まずは `matplotlib` を用いてプロットしてみましょう。


```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 2, figsize=(10, 6))

ax[0][0].plot(step, E)
ax[0][1].plot(step, T)
ax[1][0].plot(step, V)
ax[1][1].plot(step, H)
ax[0][0].set_xlabel("MD steps")
ax[0][1].set_xlabel("MD steps")
ax[1][0].set_xlabel("MD steps")
ax[1][1].set_xlabel("MD steps")
ax[0][0].set_ylabel("potential energy [eV/atom]")
ax[0][1].set_ylabel("temperature [K]")
ax[1][0].set_ylabel("volume [A^3/atom]")
ax[1][1].set_ylabel("enthalpy [eV/atom]")
```


    Text(0, 0.5, 'enthalpy [eV/atom]')



![png](output_16_1.png)



次に、 `plotly` を用いてプロットしてみましょう。


```python
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots

fig = make_subplots(rows=2, cols=2)
fig.add_trace(go.Scatter(x=step, y=E, mode="lines", showlegend=False), row=1, col=1)
fig.add_trace(go.Scatter(x=step, y=T, mode="lines", showlegend=False), row=1, col=2)
fig.add_trace(go.Scatter(x=step, y=V, mode="lines", showlegend=False), row=2, col=1)
fig.add_trace(go.Scatter(x=step, y=H, mode="lines", showlegend=False), row=2, col=2)
fig.update_xaxes(title="MD steps", row=1, col=1)
fig.update_xaxes(title="MD steps", row=1, col=2)
fig.update_xaxes(title="MD steps", row=2, col=1)
fig.update_xaxes(title="MD steps", row=2, col=2)
fig.update_yaxes(title="potential energy [eV/atom]", row=1, col=1)
fig.update_yaxes(title="temperature [K]", row=1, col=2)
fig.update_yaxes(title="volume [A^3/atom]", row=2, col=1)
fig.update_yaxes(title="enthalpy [eV/atom]", row=2, col=2)
fig.update_layout(width=1000, height=600)
```


&lt;a id='md_4_jp'&gt;&lt;/a&gt;
### MD 平均二乗変位

平均二乗変位 (Mean squared displacement; MSD) により系の中の粒子の移動しやすさ表しています。 拡散係数はMSDから計算することができます。この例ではMSDのプロット方法を紹介します。

ここでは、事前に実行済みのMDのトラジェクトリを使用します。実際のMDの実行方法やトラジェクトリの保存方法については `Welcome` tutorialをご覧ください。

`ase.md.analysis.DiffusionCoefficient` はMSDを計算するためのクラスです。
- `timestep`: トラジェクトリの2つの隣接イメージ間の時間幅。
- `ignore_n_images`: トラジェクトリの始点から無視したいイメージの数。系が平衡化するまでのイメージを取り除くために指定します。


```python
from ase import units
from ase.io import Trajectory
from ase.md.analysis import DiffusionCoefficient

traj = Trajectory("assets/visualization/md_example.traj")[::4]
diffusion = DiffusionCoefficient(traj, timestep=400*units.fs)
diffusion.calculate(ignore_n_images=0)
```


MSDを線形フィッティングした際の傾きから、$MSD(t) = 2nDt$の関係を用いることで拡散係数$D$を求められます。ここで、$n$は系の次元です。


```python
%pip install joblib
```





```python
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)
diffusion.plot(ax=ax)
```



![png](output_23_0.png)



&lt;a id='md_5_jp'&gt;&lt;/a&gt;
### MD 平均二乗変位 (matlantis_features)
`matlantis-features` では、ASEよりも高速で動作する mean square displacement (MSD) 計算機能を提供しています。

ここでは、予め保存しておいた `PostMDDiffusionFeature` 形式の計算結果を用います。この結果の取得方法については、`Diffusivity` example および guidebookの`PostMDDiffusionFeature` をご参照ください。

`PostMDDiffusionFeature` がプロット機能を持っており(`PostMDDiffusionFeature.plot`)、`plotly.Figure`形式のオブジェクトを生成します。以下の例のようにプロットされる図をカスタマイズすることができます。


```python
%pip install joblib
```





```python
from joblib import load

diffusion_results = load("assets/visualization/diffusion_results.pk")

fig = diffusion_results.plot(show_effective_range=True, show_fit_line=True)
fig.update_layout(title="MSE of ethanol at 300 K")
fig.update_layout(width=800, height=500)
fig.add_annotation(
    text=(
        "Diffusion coefficient &lt;br&gt;"
        f"H: {diffusion_results.diffusion_coefficient['cm^2/s']['H']:.3e} cm^2/s &lt;br&gt;"
        f"C: {diffusion_results.diffusion_coefficient['cm^2/s']['C']:.3e} cm^2/s &lt;br&gt;"
        f"O: {diffusion_results.diffusion_coefficient['cm^2/s']['O']:.3e} cm^2/s &lt;br&gt;"
    ),
    xref="paper", yref="paper",
    x=0.1, y=0.9, showarrow=False
)
fig.show()
```


&lt;a id='md_6_jp'&gt;&lt;/a&gt;
### MD 応力 ACF と粘性係数 (matlantis_features)
この例では、応力のautocorrelation function (ACF)をプロットし、粘性係数を計算する方法を示します。 `matlantis-features` の `PostMDViscosityFeature` を用いることで応力ACFおよび粘性係数を計算することができます。

ここでは、予め保存しておいた計算結果を読み込みます。 この結果の取得方法については、`Viscosity` example および guidebookの`PostEMDViscosityFeature` をご参照ください。
`PostMDViscosityFeature` がプロット機能を持っており(`PostMDViscosityFeature.plot`)、`plotly.Figure`形式のオブジェクトを生成します。


```python
%pip install joblib
```





```python
from joblib import load

viscosity_results = load("assets/visualization/viscosity_results.pk")

fig = viscosity_results.plot()
fig.update_layout(title="Viscosity of ethanol at 300 K")
fig.update_layout(width=800, height=500)

fig.add_annotation(
    text=(
        f"Viscosity: {viscosity_results.viscosity['mPa s'][0]:.3f} mPa s &lt;br&gt;"
    ),
    xref="paper", yref="paper",
    x=0.85, y=0.95, showarrow=False,
    font={"size":16}
)

fig.show()
```


&lt;a id='phonon_1_jp'&gt;&lt;/a&gt;
### フォノン 分散と状態密度 (matlantis_features)
このexampleでは、フォノン分散とDOSを `matlantis-features` の機能を用いてプロットします。

ここでは、予め保存しておいた力定数をファイルから読み込みます。 力定数の計算方法については、`Phonon` example およびguidebookの `ForceConstantFeature` をご参照ください。

フォノン分散およびDOSはそれぞれ`PostPhononBandFeature`, `PostPhononDOSFeature`によってプロットできます。


```python
%pip install joblib
```





```python
from joblib import load
from matlantis_features.features.phonon import (
    PostPhononBandFeature,
    PostPhononDOSFeature,
    PostPhononThermochemistryFeature,
    plot_band_dos,
)

force_constant = load("assets/visualization/Si_fc.pk")

band = PostPhononBandFeature()
band_results = band(force_constant)

dos = PostPhononDOSFeature()
dos_results = dos(force_constant, kpts = [20,20,20])
```


```python
fig = plot_band_dos(band_results, dos_results)
fig.update_layout(width=800, height=500)
fig.show()
```


&lt;a id='phonon_2_jp'&gt;&lt;/a&gt;
### フォノン 熱力学特性 (matlantis_features)
このexampleでは、フォノンから計算される熱力学特性のプロット方法を示します。

ここでは、予め保存しておいた力定数をファイルから読み込みます。 力定数の計算方法については、`Phonon` example およびguidebookの `ForceConstantFeature` をご参照ください。

プロットには`PostPhononThermochemistryFeature` の `plot` 関数を用います。


```python
%pip install joblib
```





```python
from joblib import load
from matlantis_features.features.phonon import PostPhononThermochemistryFeature

force_constant = load("assets/visualization/Si_fc.pk")

thermo = PostPhononThermochemistryFeature()
thermo_results = thermo(
    force_constant,
    kpts = [20, 20, 20],
)
```


```python
fig = thermo_results.plot()
fig.update_layout(width=800, height=500)
fig.show()
```


&lt;a id='neb_1_jp'&gt;&lt;/a&gt;
### 反応経路 エネルギープロファイル

このexampleでは、反応経路のエネルギープロファイルをプロットする方法を示します。

反応経路はnudged elastic band (NEB) 法によって計算することができます。ここでは、予め計算しておいた反応経路をファイルから読み込みます。 NEB計算の詳細については `NEB` tutorial をご参照ください。。

エネルギープロファイルのプロットに関しては、ASEが強力なプロット機能`ase.mep.NEBTools`を提供しています。


```python
from ase.io import Trajectory

neb_traj = Trajectory("assets/visualization/neb.traj")
```


```python
from ase.mep import NEBTools

fig = NEBTools(neb_traj).plot_band()
```



![png](output_40_0.png)


