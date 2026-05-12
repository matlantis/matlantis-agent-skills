# Matlantis-features: Phonon QHA, Gibbs Free Energy &amp; T-P Phase diagram


## 概要

このexampleの前半では、シリコン結晶を例に`matlantis-features`を用いたGibbs自由エネルギーの計算方法を紹介します。後半ではGibbs自由エネルギーの計算結果を用いてシリコン結晶の温度(T)-圧力(P)相図を図示します。

Gibbs自由エネルギーの計算は **準調和近似(quasi-harmonic approximation; QHA)** に基づきます。一般に、フォノンの計算は調和振動子模型に基づいており、原子間の相互作用は調和振動のみで記述されます。この仮定の下では熱膨張率やGibbs自由エネルギーなどの体積が関係する物理量を表現することができません。そこで、準調和近似では特定の温度、圧力下でのGibbs自由エネルギーを以下の手続きで計算します。

* 体積(V)を初期結晶構造から $V_0$, $V_1$, ... $V_n$ へ変化させる。
* それぞれの体積でのフォノン分散関係を計算し、それぞれの温度 ($T_1$, $T_2$, ... $T_m$ K)におけるHelmholtz自由エネルギー (F)を得る。
* 温度 $T_m$ Kにおける $F$-$V$ の関係式を得る。一般にこの関係式はBirch-Murnaghan状態方程式またはその他の状態方程式によって近似できる。
* 次のステップで必要となる、温度 $T_i$, 圧力 $P_j$における結晶の体積を得たい。
  熱力学によると、圧力はHelmholtz自由エネルギー(F)の体積(V)による偏微分で定義される。
$$
P=- \frac {\partial F} {\partial V}
$$
  そこで、温度 $T_i$における$F$-$V$ の関係式（Birch-Murnaghan状態方程式）の微分と負の圧力　$-P_j$ が等しくなる体積を探すことで、求めたい体積 $V_{P=P_j}$ を得ることができる。
* 以上により、次の公式からGibbs自由エネルギーを計算することができる。
$$
G(T_i, P_j) = F(T_i, V_{P=P_j}) + P_j V_{P=P_j}
$$
なお、この例で使用する構造ファイルは [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) の下で Materials Project から入手したものです。
* Dimond Si (mp-149, Primitive Cell)
* β-Sn Si (mp-92, Primitive Cell)
* Clathrate Si (mp-16220, Primitive Cell)
* Simple hexagonal Si (mp-34, Primitive Cell)


## 初期設定

* 必要なライブラリ等の準備を行います。


```python
# If you have already installed matlantis-features, you can skip this.
!pip install 'matlantis-features&gt;=0.16.1'
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
import pathlib
try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()
```


```python
import logging

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)
```


## シリコン結晶のGibbs自由エネルギー


### 初期設定

* 必要なライブラリ等の準備を行います。


```python
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.offline import iplot
from tqdm.auto import tqdm

from ase import units
from ase.io import read

from matlantis_features.features.common.opt import BFGSASEOptFeature
from matlantis_features.features.phonon import ForceConstantFeature
from matlantis_features.features.phonon import PostPhononBandFeature
from matlantis_features.features.phonon import PostPhononQHAGibbsFeature
from matlantis_features.functions.eos import BirchMurnaghanEOS
from matlantis_features.utils.calculators import get_calculator
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


### 状態方程式の計算
状態方程式(equation of states; EOS)は材料の基本的な特性の一つです。EOSから体積とポテンシャルエネルギーの関係を知ることができます。
自由エネルギー計算に直接に使用する訳ではありませんが、最初にシリコン結晶のEOSを見ることにします。

ダイヤモンド構造のシリコン結晶構造を読み込んで体積(V)を80%~140%の範囲で変化させた後、各体積でのエネルギー(E)を計算してE-Vの関係を描画するとBirch-Murnaghan状態方程式に良く従っていることが分かります。


```python
atoms = read(str(dir_path/"assets/phonon_qha/Si_diamond.xyz"))

bfgs = BFGSASEOptFeature(fmax=0.001, filter=True, show_progress_bar=True, estimator_fn=estimator_fn)
opt_res = bfgs(atoms)

atoms_opt = opt_res.atoms.ase_atoms
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


```python
calculator = get_calculator()

deforms = np.arange(0.8, 1.4, 0.04)
V_eos = []
E_eos = []
for i, deform in enumerate(deforms):
    atoms_deformed = atoms_opt.copy()
    atoms_deformed.set_cell(atoms_opt.get_cell()[:] * deform ** (1 / 3.0), scale_atoms=True)
    atoms_deformed.calc = calculator
    V_eos.append(atoms_deformed.get_volume()/len(atoms_opt))
    E_eos.append(atoms_deformed.get_potential_energy()/len(atoms_opt))

eos = BirchMurnaghanEOS.from_fit(V_eos, E_eos)
V = np.linspace(V_eos[0], V_eos[-1], 101, endpoint=True)
E_pred = eos(V)
plt.figure()
plt.scatter(V_eos, E_eos)
plt.plot(V, E_pred, c="k", label="Birch-Murnaghan EOS")
plt.xlabel("Volume (Å)")
plt.ylabel("Potential energy (eV/atom)")
plt.legend()

print(f"Bulk modulus: {eos.bulk_modulus/units.GPa:3.3f} GPa")

```

    Bulk modulus: 89.137 GPa



![png](output_14_1.png)



### Gibbs自由エネルギーの計算
`ForceConstantFeature`を用いて各体積におけるシリコン結晶の力定数を計算します。以下の例の中で、力定数計算前に構造最適化をしていることに注意してください。


```python
force_constant_list = []
opt = BFGSASEOptFeature(n_run=1000, fmax=0.001, estimator_fn=estimator_fn)
force_constant_feature = ForceConstantFeature(supercell=(10, 10, 10), delta=0.2, estimator_fn=estimator_fn)
for i, deform in tqdm(enumerate(deforms), total=len(deforms)):
    atoms_deformed = atoms_opt.copy()
    atoms_deformed.set_cell(atoms_opt.get_cell()[:] * deform ** (1 / 3.0))
    opt_result = opt(atoms_deformed)
    force_constant_list.append(force_constant_feature(opt_result.atoms))

```


      0%|          | 0/15 [00:00&lt;?, ?it/s]


各構造に対するフォノン分散関係は次図のようになります。


```python
band = PostPhononBandFeature()
plt.figure()
cmap = plt.get_cmap("tab20")
for i, (fc, deform) in enumerate(zip(force_constant_list, deforms)):
    band_res = band(fc)
    x = band_res.coords
    for j, freq in enumerate(band_res.frequency.T):
        if (i == 0 or i == (len(force_constant_list)-1)) and j == 0:
            plt.plot(x, freq, c=cmap(i), label=f"deform: {deform:.2f}")
        else:
            plt.plot(x, freq, c=cmap(i))

plt.legend()
```


    &lt;matplotlib.legend.Legend at 0x7f03fbec7ed0&gt;



![png](output_18_1.png)



各体積におけるシリコン結晶の力定数を`PostPhononQHAGibbsFeature`に入力し、温度と圧力を指定することでGibbs自由エネルギー(G)を計算することができます。


```python
temperature = list(range(0, 1500, 100))
pressure = [0.0, 10*units.GPa]

gibbs = PostPhononQHAGibbsFeature()
gibbs_results = gibbs(
    force_constant_list,
    kpts=[10, 10, 10],
    temperatures=temperature,
    pressures=pressure,
)
```


温度0, 100, 200 ... 1400K、圧力0GPa, 10GPaにおけるダイヤモンド結晶シリコンのGibbs自由エネルギーは次のようになります。単位は`eV/atom`です。


```python
gibbs_results.gibbs_free_energy
```


    array([[-4.49031244, -3.27308245],
           [-4.49159294, -3.27479806],
           [-4.50009673, -3.28386445],
           [-4.51677878, -3.30062214],
           [-4.54042589, -3.32403876],
           [-4.56981894, -3.35302774],
           [-4.60402184, -3.38672035],
           [-4.64233399, -3.42445205],
           [-4.68422202, -3.4657105 ],
           [-4.72926996, -3.51009306],
           [-4.77714591, -3.55727675],
           [-4.82757957, -3.60699763],
           [-4.88034692, -3.65903637],
           [-4.93525947, -3.71320811],
           [-4.99215655, -3.76935502]])


QHAの計算結果を可視化します。左図の水平方向の曲線は各温度におけるF-Vの関係を、垂直方向の曲線は圧力0GPa, 10GPaの下での各温度の平衡体積を表します。右図の2曲線はそれぞれ圧力0GPa, 10GPaの下でのGibbs自由エネルギーと温度の関係を表します。


```python
fig = gibbs_results.plot()
fig.update_layout(width=1200, height=500)
iplot(fig)
```


```python
fig2 = gibbs_results.plot_3d()
fig2.update_layout(width=800, height=700, scene_camera_eye={"x": 1.25, "y": 1.25, "z": 0.9})
```


## Si T-P 相図


### 初期設定

* 必要なライブラリ等の準備を行います。


```python
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.offline import iplot
from tqdm.auto import tqdm

from ase import units
from ase.io import read

from matlantis_features.features.base import FeatureBase
from matlantis_features.features.common.opt import BFGSASEOptFeature
from matlantis_features.features.phonon import ForceConstantFeature
from matlantis_features.features.phonon import PostPhononQHAGibbsFeature
```


後半ではGibbs自由エネルギーを簡単に計算できるFeatureを実装していきます。計算方法は前半と同様です。ここでは`BFGSASEOptFeature`, `ForceConstantFeature`, `PostPhononQHAGibbsFeature`を組み合わせることで実装します。


```python
class ComplexGibbsFreeEnergyFeature(FeatureBase):
    def __init__(self, supercell, delta):
        super(ComplexGibbsFreeEnergyFeature, self).__init__()
        self.opt_pre = BFGSASEOptFeature(n_run=1000, fmax=0.001, filter=True, estimator_fn=estimator_fn)
        self.opt = BFGSASEOptFeature(n_run=1000, fmax=0.001, estimator_fn=estimator_fn)
        self.fc = ForceConstantFeature(supercell=supercell, delta=delta, estimator_fn=estimator_fn)
        self.gibbs = PostPhononQHAGibbsFeature()
    def __call__(self, atoms, deforms, temperature, pressure):
        opt_res = self.opt_pre(atoms)
        atoms_opt = opt_res.atoms.ase_atoms
        force_constant_list = []
        for i, deform in tqdm(enumerate(deforms), total=len(deforms)):
            atoms_deformed = atoms_opt.copy()
            atoms_deformed.set_cell(atoms_opt.get_cell()[:] * deform ** (1 / 3.0), scale_atoms=True)
            opt_result = self.opt(atoms_deformed)
            force_constant_list.append(self.fc(opt_result.atoms))
        gibbs_results = self.gibbs(
            force_constant_list,
            kpts=[10, 10, 10],
            temperatures=temperature,
            pressures=pressure,
        )
        return gibbs_results
```


上で定義したFeatureを使っていくつかの有名なシリコン結晶構造についてGibbs自由エネルギーを計算します。


```python
label = [
    "diamond",
    "beta_Sn",
    "clathrate",
    "sh",
]
supercell = [
    (10, 10, 10),
    (10, 10, 10),
    (4, 4, 4),
    (12, 12, 12),
]
deforms = [
    list(np.arange(0.8, 1.4, 0.04)),
    list(np.arange(0.8, 1.4, 0.04)),
    list(np.arange(0.84, 1.5, 0.04)),
    list(np.arange(0.8, 1.3, 0.04)),
]

gibbs_results = []
for l, s, d in zip(label, supercell, deforms):
    atoms = read(str(dir_path/f"assets/phonon_qha/Si_{l}.xyz"))
    gibbs = ComplexGibbsFreeEnergyFeature(supercell=s, delta=0.2)
    gibbs_results.append(
        gibbs(
            atoms,
            deforms=d,
            temperature = list(range(0, 1500, 25)),
            pressure = list(np.arange(-15, 20, 0.5)*units.GPa),
        )
    )
```


      0%|          | 0/15 [00:00&lt;?, ?it/s]


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 0.00 K -15.00 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 0.00 K -14.50 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 0.00 K -14.00 GPa
    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 0.00 K -13.50 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 25.00 K -15.00 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 25.00 K -14.50 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 25.00 K -14.00 GPa
    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 25.00 K -13.50 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 50.00 K -15.00 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 50.00 K -14.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 50.00 K -14.00 GPa
    The equlibrium volume 28.3 exceed valid range (16.3 27.8) at 50.00 K -13.50 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 75.00 K -15.00 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 75.00 K -14.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 75.00 K -14.00 GPa
    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 75.00 K -13.50 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 100.00 K -15.00 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 100.00 K -14.50 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 100.00 K -14.00 GPa
    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 100.00 K -13.50 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 125.00 K -15.00 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 125.00 K -14.50 GPa
    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 125.00 K -14.00 GPa
    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 125.00 K -13.50 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 150.00 K -15.00 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 150.00 K -14.50 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 150.00 K -14.00 GPa
    The equlibrium volume 28.5 exceed valid range (16.3 27.8) at 150.00 K -13.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 175.00 K -15.00 GPa
    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 175.00 K -14.50 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 175.00 K -14.00 GPa
    The equlibrium volume 28.6 exceed valid range (16.3 27.8) at 175.00 K -13.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 200.00 K -15.00 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 200.00 K -14.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 200.00 K -14.00 GPa
    The equlibrium volume 28.7 exceed valid range (16.3 27.8) at 200.00 K -13.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 225.00 K -15.00 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 225.00 K -14.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 225.00 K -14.00 GPa
    The equlibrium volume 28.9 exceed valid range (16.3 27.8) at 225.00 K -13.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 250.00 K -15.00 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 250.00 K -14.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 250.00 K -14.00 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 250.00 K -13.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 275.00 K -15.00 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 275.00 K -14.50 GPa
    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 275.00 K -14.00 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 275.00 K -13.50 GPa
    The equlibrium volume 27.8 exceed valid range (16.3 27.8) at 275.00 K -13.00 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 300.00 K -15.00 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 300.00 K -14.50 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 300.00 K -14.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 300.00 K -13.50 GPa
    The equlibrium volume 27.9 exceed valid range (16.3 27.8) at 300.00 K -13.00 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 325.00 K -15.00 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 325.00 K -14.50 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 325.00 K -14.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 325.00 K -13.50 GPa
    The equlibrium volume 28.0 exceed valid range (16.3 27.8) at 325.00 K -13.00 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 350.00 K -15.00 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 350.00 K -14.50 GPa
    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 350.00 K -14.00 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 350.00 K -13.50 GPa
    The equlibrium volume 28.2 exceed valid range (16.3 27.8) at 350.00 K -13.00 GPa
    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 375.00 K -15.00 GPa
    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 375.00 K -14.50 GPa
    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 375.00 K -14.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 375.00 K -13.50 GPa
    The equlibrium volume 28.3 exceed valid range (16.3 27.8) at 375.00 K -13.00 GPa
    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 400.00 K -15.00 GPa
    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 400.00 K -14.50 GPa
    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 400.00 K -14.00 GPa
    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 400.00 K -13.50 GPa
    The equlibrium volume 28.5 exceed valid range (16.3 27.8) at 400.00 K -13.00 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 425.00 K -15.00 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 425.00 K -14.50 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 425.00 K -14.00 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 425.00 K -13.50 GPa
    The equlibrium volume 28.6 exceed valid range (16.3 27.8) at 425.00 K -13.00 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 450.00 K -15.00 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 450.00 K -14.50 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 450.00 K -14.00 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 450.00 K -13.50 GPa
    The equlibrium volume 28.8 exceed valid range (16.3 27.8) at 450.00 K -13.00 GPa
    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 475.00 K -15.00 GPa
    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 475.00 K -14.50 GPa
    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 475.00 K -14.00 GPa
    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 475.00 K -13.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 475.00 K -13.00 GPa
    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 500.00 K -15.00 GPa
    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 500.00 K -14.50 GPa
    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 500.00 K -14.00 GPa
    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 500.00 K -13.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 500.00 K -13.00 GPa
    The equlibrium volume 27.8 exceed valid range (16.3 27.8) at 500.00 K -12.50 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 525.00 K -15.00 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 525.00 K -14.50 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 525.00 K -14.00 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 525.00 K -13.50 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 525.00 K -13.00 GPa
    The equlibrium volume 27.9 exceed valid range (16.3 27.8) at 525.00 K -12.50 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 550.00 K -15.00 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 550.00 K -14.50 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 550.00 K -14.00 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 550.00 K -13.50 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 550.00 K -13.00 GPa
    The equlibrium volume 28.0 exceed valid range (16.3 27.8) at 550.00 K -12.50 GPa
    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 575.00 K -15.00 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 575.00 K -14.50 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 575.00 K -14.00 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 575.00 K -13.50 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 575.00 K -13.00 GPa
    The equlibrium volume 28.2 exceed valid range (16.3 27.8) at 575.00 K -12.50 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -15.00 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -14.50 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -14.00 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -13.50 GPa
    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -13.00 GPa
    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 600.00 K -12.50 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -15.00 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -14.50 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -14.00 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -13.50 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -13.00 GPa
    The equlibrium volume 28.6 exceed valid range (16.3 27.8) at 625.00 K -12.50 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -15.00 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -14.50 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -14.00 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -13.50 GPa
    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -13.00 GPa
    The equlibrium volume 28.9 exceed valid range (16.3 27.8) at 650.00 K -12.50 GPa
    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -15.00 GPa
    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -14.50 GPa
    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -14.00 GPa
    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -13.50 GPa
    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -13.00 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 675.00 K -12.50 GPa
    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 700.00 K -15.00 GPa
    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 700.00 K -14.50 GPa
    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 700.00 K -14.00 GPa
    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 700.00 K -13.50 GPa
    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 700.00 K -13.00 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 700.00 K -12.50 GPa
    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -15.00 GPa
    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -14.50 GPa
    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -14.00 GPa
    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -13.50 GPa
    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -13.00 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 725.00 K -12.50 GPa
    The equlibrium volume 27.8 exceed valid range (16.3 27.8) at 725.00 K -12.00 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -15.00 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -14.50 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -14.00 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -13.50 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -13.00 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -12.50 GPa
    The equlibrium volume 28.0 exceed valid range (16.3 27.8) at 750.00 K -12.00 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -15.00 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -14.50 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -14.00 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -13.50 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -13.00 GPa
    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -12.50 GPa
    The equlibrium volume 28.2 exceed valid range (16.3 27.8) at 775.00 K -12.00 GPa
    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -15.00 GPa
    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -14.50 GPa
    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -14.00 GPa
    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -13.50 GPa
    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -13.00 GPa
    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -12.50 GPa
    The equlibrium volume 28.5 exceed valid range (16.3 27.8) at 800.00 K -12.00 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -15.00 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -14.50 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -14.00 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -13.50 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -13.00 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -12.50 GPa
    The equlibrium volume 28.8 exceed valid range (16.3 27.8) at 825.00 K -12.00 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -15.00 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -14.50 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -14.00 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -13.50 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -13.00 GPa
    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -12.50 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 850.00 K -12.00 GPa
    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -15.00 GPa
    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -14.50 GPa
    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -14.00 GPa
    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -13.50 GPa
    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -13.00 GPa
    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -12.50 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 875.00 K -12.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -15.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -14.50 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -14.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -13.50 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -13.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -12.50 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -12.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -15.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -14.50 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -14.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -13.50 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -13.00 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -12.50 GPa
    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -12.00 GPa
    The equlibrium volume 27.9 exceed valid range (16.3 27.8) at 925.00 K -11.50 GPa
    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -15.00 GPa
    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -14.50 GPa
    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -14.00 GPa
    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -13.50 GPa
    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -13.00 GPa
    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -12.50 GPa
    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -12.00 GPa
    The equlibrium volume 28.1 exceed valid range (16.3 27.8) at 950.00 K -11.50 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -15.00 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -14.50 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -14.00 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -13.50 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -13.00 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -12.50 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -12.00 GPa
    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 975.00 K -11.50 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -15.00 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -14.50 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -14.00 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -13.50 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -13.00 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -12.50 GPa
    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -12.00 GPa
    The equlibrium volume 28.7 exceed valid range (16.3 27.8) at 1000.00 K -11.50 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -15.00 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -14.50 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -14.00 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -13.50 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -13.00 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -12.50 GPa
    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -12.00 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1025.00 K -11.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -15.00 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -14.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -14.00 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -13.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -13.00 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -12.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -12.00 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -11.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -15.00 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -14.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -14.00 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1075.00 K -13.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -13.00 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -12.50 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -12.00 GPa
    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -11.50 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -15.00 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -14.50 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -14.00 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -13.50 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -13.00 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -12.50 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -12.00 GPa
    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -11.50 GPa
    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -15.00 GPa
    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -14.50 GPa
    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -14.00 GPa
    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -13.50 GPa
    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -13.00 GPa
    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -12.50 GPa
    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -12.00 GPa
    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -11.50 GPa
    The equlibrium volume 28.0 exceed valid range (16.3 27.8) at 1125.00 K -11.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -15.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -14.50 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -14.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -13.50 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -13.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -12.50 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -12.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -11.50 GPa
    The equlibrium volume 28.3 exceed valid range (16.3 27.8) at 1150.00 K -11.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -15.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -14.50 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -14.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -13.50 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -13.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -12.50 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -12.00 GPa
    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -11.50 GPa
    The equlibrium volume 28.8 exceed valid range (16.3 27.8) at 1175.00 K -11.00 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -15.00 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -14.50 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -14.00 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -13.50 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -13.00 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -12.50 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -12.00 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -11.50 GPa
    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -11.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -15.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -14.50 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -14.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -13.50 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -13.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -12.50 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -12.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -11.50 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -11.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -15.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -14.50 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1250.00 K -14.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -13.50 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -13.00 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -12.50 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -12.00 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1250.00 K -11.50 GPa
    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -11.00 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -15.00 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -14.50 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -14.00 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -13.50 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -13.00 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -12.50 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -12.00 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -11.50 GPa
    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -11.00 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -15.00 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -14.50 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -14.00 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -13.50 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -13.00 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -12.50 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -12.00 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -11.50 GPa
    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -11.00 GPa
    The equlibrium volume 27.9 exceed valid range (16.3 27.8) at 1300.00 K -10.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -15.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -14.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -14.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -13.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -13.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -12.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -12.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -11.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -11.00 GPa
    The equlibrium volume 28.3 exceed valid range (16.3 27.8) at 1325.00 K -10.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -15.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -14.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -14.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -13.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -13.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -12.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -12.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -11.50 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -11.00 GPa
    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -10.50 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -15.00 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -14.50 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -14.00 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -13.50 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -13.00 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -12.50 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -12.00 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -11.50 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -11.00 GPa
    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -10.50 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -15.00 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -14.50 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -14.00 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -13.50 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -13.00 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -12.50 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -12.00 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -11.50 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -11.00 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -10.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -15.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -14.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -14.00 GPa
    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1425.00 K -13.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -13.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -12.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -12.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -11.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -11.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -10.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -15.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -14.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -14.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -13.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -13.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -12.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -12.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -11.50 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -11.00 GPa
    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -10.50 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -15.00 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -14.50 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -14.00 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -13.50 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -13.00 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -12.50 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -12.00 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -11.50 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -11.00 GPa
    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -10.50 GPa


      0%|          | 0/15 [00:00&lt;?, ?it/s]


    The equlibrium volume 21.1 exceed valid range (12.4 21.0) at 0.00 K -15.00 GPa
    The equlibrium volume 21.1 exceed valid range (12.4 21.0) at 25.00 K -15.00 GPa
    The equlibrium volume 21.2 exceed valid range (12.4 21.0) at 50.00 K -15.00 GPa
    The equlibrium volume 21.3 exceed valid range (12.4 21.0) at 75.00 K -15.00 GPa
    The equlibrium volume 21.4 exceed valid range (12.4 21.0) at 100.00 K -15.00 GPa
    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 125.00 K -15.00 GPa
    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 150.00 K -15.00 GPa
    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 175.00 K -15.00 GPa
    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 200.00 K -15.00 GPa
    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 225.00 K -15.00 GPa
    The equlibrium volume 21.4 exceed valid range (12.4 21.0) at 250.00 K -15.00 GPa
    The equlibrium volume 21.3 exceed valid range (12.4 21.0) at 275.00 K -15.00 GPa
    The equlibrium volume 21.2 exceed valid range (12.4 21.0) at 300.00 K -15.00 GPa
    The equlibrium volume 21.1 exceed valid range (12.4 21.0) at 325.00 K -15.00 GPa
    The equlibrium volume 21.0 exceed valid range (12.4 21.0) at 350.00 K -15.00 GPa


      0%|          | 0/17 [00:00&lt;?, ?it/s]


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 0.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 25.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 50.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 75.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 100.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 125.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 150.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 175.00 K 19.00 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 175.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 200.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 200.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 225.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 225.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 250.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 250.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 275.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 275.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 300.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 300.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 325.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 325.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 350.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 350.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 375.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 375.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 400.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 400.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 425.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 425.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 450.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 450.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 475.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 475.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 500.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 500.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 525.00 K 18.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 525.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 525.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 550.00 K 18.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 550.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 550.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 575.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 575.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 575.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 600.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 600.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 600.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 625.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 625.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 625.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 650.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 650.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 650.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 675.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 675.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 675.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 700.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 700.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 700.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 725.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 725.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 725.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 750.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 750.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 750.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 775.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 775.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 775.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 800.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 800.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 800.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 825.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 825.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 825.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 850.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 850.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 850.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 875.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 875.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 875.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 900.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 900.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 900.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 925.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 925.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 925.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 950.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 950.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 950.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 975.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 975.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 975.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1000.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1000.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1000.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1025.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1025.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1025.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1050.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1050.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1050.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1075.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1075.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1075.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1100.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1100.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1100.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1125.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1125.00 K 19.00 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1125.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1150.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1150.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1150.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1175.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1175.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1175.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1200.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1200.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1200.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1225.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1225.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1225.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1250.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1250.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1250.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1275.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1275.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1275.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1300.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1300.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1300.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1325.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1325.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1325.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1350.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1350.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1350.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1375.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1375.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1375.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1400.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1400.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1400.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1425.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1425.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1425.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1450.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1450.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1450.00 K 19.50 GPa
    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1475.00 K 18.50 GPa
    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1475.00 K 19.00 GPa
    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1475.00 K 19.50 GPa


      0%|          | 0/13 [00:00&lt;?, ?it/s]


    The equlibrium volume 22.0 exceed valid range (12.1 19.4) at 0.00 K -15.00 GPa
    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 0.00 K -14.50 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 0.00 K -14.00 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 0.00 K -13.50 GPa
    The equlibrium volume 22.1 exceed valid range (12.1 19.4) at 25.00 K -15.00 GPa
    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 25.00 K -14.50 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 25.00 K -14.00 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 25.00 K -13.50 GPa
    The equlibrium volume 22.6 exceed valid range (12.1 19.4) at 50.00 K -15.00 GPa
    The equlibrium volume 20.9 exceed valid range (12.1 19.4) at 50.00 K -14.50 GPa
    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 50.00 K -14.00 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 50.00 K -13.50 GPa
    The equlibrium volume 23.0 exceed valid range (12.1 19.4) at 75.00 K -15.00 GPa
    The equlibrium volume 21.0 exceed valid range (12.1 19.4) at 75.00 K -14.50 GPa
    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 75.00 K -14.00 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 75.00 K -13.50 GPa
    The equlibrium volume 23.1 exceed valid range (12.1 19.4) at 100.00 K -15.00 GPa
    The equlibrium volume 21.2 exceed valid range (12.1 19.4) at 100.00 K -14.50 GPa
    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 100.00 K -14.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 100.00 K -13.50 GPa
    The equlibrium volume 23.1 exceed valid range (12.1 19.4) at 125.00 K -15.00 GPa
    The equlibrium volume 21.4 exceed valid range (12.1 19.4) at 125.00 K -14.50 GPa
    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 125.00 K -14.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 125.00 K -13.50 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 125.00 K -13.00 GPa
    The equlibrium volume 23.2 exceed valid range (12.1 19.4) at 150.00 K -15.00 GPa
    The equlibrium volume 21.5 exceed valid range (12.1 19.4) at 150.00 K -14.50 GPa
    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 150.00 K -14.00 GPa
    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 150.00 K -13.50 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 150.00 K -13.00 GPa
    The equlibrium volume 23.3 exceed valid range (12.1 19.4) at 175.00 K -15.00 GPa
    The equlibrium volume 21.6 exceed valid range (12.1 19.4) at 175.00 K -14.50 GPa
    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 175.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 175.00 K -13.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 175.00 K -13.00 GPa
    The equlibrium volume 23.4 exceed valid range (12.1 19.4) at 200.00 K -15.00 GPa
    The equlibrium volume 21.6 exceed valid range (12.1 19.4) at 200.00 K -14.50 GPa
    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 200.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 200.00 K -13.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 200.00 K -13.00 GPa
    The equlibrium volume 23.5 exceed valid range (12.1 19.4) at 225.00 K -15.00 GPa
    The equlibrium volume 21.6 exceed valid range (12.1 19.4) at 225.00 K -14.50 GPa
    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 225.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 225.00 K -13.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 225.00 K -13.00 GPa
    The equlibrium volume 23.7 exceed valid range (12.1 19.4) at 250.00 K -15.00 GPa
    The equlibrium volume 21.5 exceed valid range (12.1 19.4) at 250.00 K -14.50 GPa
    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 250.00 K -14.00 GPa
    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 250.00 K -13.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 250.00 K -13.00 GPa
    The equlibrium volume 23.3 exceed valid range (12.1 19.4) at 275.00 K -15.00 GPa
    The equlibrium volume 21.4 exceed valid range (12.1 19.4) at 275.00 K -14.50 GPa
    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 275.00 K -14.00 GPa
    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 275.00 K -13.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 275.00 K -13.00 GPa
    The equlibrium volume 22.6 exceed valid range (12.1 19.4) at 300.00 K -15.00 GPa
    The equlibrium volume 21.3 exceed valid range (12.1 19.4) at 300.00 K -14.50 GPa
    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 300.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 300.00 K -13.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 300.00 K -13.00 GPa
    The equlibrium volume 22.3 exceed valid range (12.1 19.4) at 325.00 K -15.00 GPa
    The equlibrium volume 21.2 exceed valid range (12.1 19.4) at 325.00 K -14.50 GPa
    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 325.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 325.00 K -13.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 325.00 K -13.00 GPa
    The equlibrium volume 22.0 exceed valid range (12.1 19.4) at 350.00 K -15.00 GPa
    The equlibrium volume 21.1 exceed valid range (12.1 19.4) at 350.00 K -14.50 GPa
    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 350.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 350.00 K -13.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 350.00 K -13.00 GPa
    The equlibrium volume 21.8 exceed valid range (12.1 19.4) at 375.00 K -15.00 GPa
    The equlibrium volume 21.0 exceed valid range (12.1 19.4) at 375.00 K -14.50 GPa
    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 375.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 375.00 K -13.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 375.00 K -13.00 GPa
    The equlibrium volume 21.6 exceed valid range (12.1 19.4) at 400.00 K -15.00 GPa
    The equlibrium volume 20.9 exceed valid range (12.1 19.4) at 400.00 K -14.50 GPa
    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 400.00 K -14.00 GPa
    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 400.00 K -13.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 400.00 K -13.00 GPa
    The equlibrium volume 21.4 exceed valid range (12.1 19.4) at 425.00 K -15.00 GPa
    The equlibrium volume 20.8 exceed valid range (12.1 19.4) at 425.00 K -14.50 GPa
    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 425.00 K -14.00 GPa
    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 425.00 K -13.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 425.00 K -13.00 GPa
    The equlibrium volume 21.2 exceed valid range (12.1 19.4) at 450.00 K -15.00 GPa
    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 450.00 K -14.50 GPa
    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 450.00 K -14.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 450.00 K -13.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 450.00 K -13.00 GPa
    The equlibrium volume 21.1 exceed valid range (12.1 19.4) at 475.00 K -15.00 GPa
    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 475.00 K -14.50 GPa
    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 475.00 K -14.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 475.00 K -13.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 475.00 K -13.00 GPa
    The equlibrium volume 21.0 exceed valid range (12.1 19.4) at 500.00 K -15.00 GPa
    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 500.00 K -14.50 GPa
    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 500.00 K -14.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 500.00 K -13.50 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 500.00 K -13.00 GPa
    The equlibrium volume 20.8 exceed valid range (12.1 19.4) at 525.00 K -15.00 GPa
    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 525.00 K -14.50 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 525.00 K -14.00 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 525.00 K -13.50 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 525.00 K -13.00 GPa
    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 550.00 K -15.00 GPa
    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 550.00 K -14.50 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 550.00 K -14.00 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 550.00 K -13.50 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 550.00 K -13.00 GPa
    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 575.00 K -15.00 GPa
    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 575.00 K -14.50 GPa
    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 575.00 K -14.00 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 575.00 K -13.50 GPa
    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 600.00 K -15.00 GPa
    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 600.00 K -14.50 GPa
    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 600.00 K -14.00 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 600.00 K -13.50 GPa
    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 625.00 K -15.00 GPa
    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 625.00 K -14.50 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 625.00 K -14.00 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 625.00 K -13.50 GPa
    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 650.00 K -15.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 650.00 K -14.50 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 650.00 K -14.00 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 650.00 K -13.50 GPa
    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 675.00 K -15.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 675.00 K -14.50 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 675.00 K -14.00 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 675.00 K -13.50 GPa
    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 700.00 K -15.00 GPa
    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 700.00 K -14.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 700.00 K -14.00 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 700.00 K -13.50 GPa
    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 725.00 K -15.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 725.00 K -14.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 725.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 750.00 K -15.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 750.00 K -14.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 750.00 K -14.00 GPa
    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 775.00 K -15.00 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 775.00 K -14.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 775.00 K -14.00 GPa
    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 800.00 K -15.00 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 800.00 K -14.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 800.00 K -14.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 825.00 K -15.00 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 825.00 K -14.50 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 825.00 K -14.00 GPa
    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 850.00 K -15.00 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 850.00 K -14.50 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 875.00 K -15.00 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 875.00 K -14.50 GPa
    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 900.00 K -15.00 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 900.00 K -14.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 925.00 K -15.00 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 925.00 K -14.50 GPa
    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 950.00 K -15.00 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 950.00 K -14.50 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 975.00 K -15.00 GPa
    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 1000.00 K -15.00 GPa
    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 1025.00 K -15.00 GPa


    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 125.00 K -14.50 GPa


    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 125.00 K -14.00 GPa


    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 125.00 K -13.50 GPa


    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 150.00 K -15.00 GPa


    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 150.00 K -14.50 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 150.00 K -14.00 GPa


    The equlibrium volume 28.5 exceed valid range (16.3 27.8) at 150.00 K -13.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 175.00 K -15.00 GPa


    The equlibrium volume 32.1 exceed valid range (16.3 27.8) at 175.00 K -14.50 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 175.00 K -14.00 GPa


    The equlibrium volume 28.6 exceed valid range (16.3 27.8) at 175.00 K -13.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 200.00 K -15.00 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 200.00 K -14.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 200.00 K -14.00 GPa


    The equlibrium volume 28.7 exceed valid range (16.3 27.8) at 200.00 K -13.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 225.00 K -15.00 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 225.00 K -14.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 225.00 K -14.00 GPa


    The equlibrium volume 28.9 exceed valid range (16.3 27.8) at 225.00 K -13.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 250.00 K -15.00 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 250.00 K -14.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 250.00 K -14.00 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 250.00 K -13.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 275.00 K -15.00 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 275.00 K -14.50 GPa


    The equlibrium volume 32.0 exceed valid range (16.3 27.8) at 275.00 K -14.00 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 275.00 K -13.50 GPa


    The equlibrium volume 27.8 exceed valid range (16.3 27.8) at 275.00 K -13.00 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 300.00 K -15.00 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 300.00 K -14.50 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 300.00 K -14.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 300.00 K -13.50 GPa


    The equlibrium volume 27.9 exceed valid range (16.3 27.8) at 300.00 K -13.00 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 325.00 K -15.00 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 325.00 K -14.50 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 325.00 K -14.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 325.00 K -13.50 GPa


    The equlibrium volume 28.0 exceed valid range (16.3 27.8) at 325.00 K -13.00 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 350.00 K -15.00 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 350.00 K -14.50 GPa


    The equlibrium volume 31.9 exceed valid range (16.3 27.8) at 350.00 K -14.00 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 350.00 K -13.50 GPa


    The equlibrium volume 28.2 exceed valid range (16.3 27.8) at 350.00 K -13.00 GPa


    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 375.00 K -15.00 GPa


    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 375.00 K -14.50 GPa


    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 375.00 K -14.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 375.00 K -13.50 GPa


    The equlibrium volume 28.3 exceed valid range (16.3 27.8) at 375.00 K -13.00 GPa


    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 400.00 K -15.00 GPa


    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 400.00 K -14.50 GPa


    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 400.00 K -14.00 GPa


    The equlibrium volume 31.8 exceed valid range (16.3 27.8) at 400.00 K -13.50 GPa


    The equlibrium volume 28.5 exceed valid range (16.3 27.8) at 400.00 K -13.00 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 425.00 K -15.00 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 425.00 K -14.50 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 425.00 K -14.00 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 425.00 K -13.50 GPa


    The equlibrium volume 28.6 exceed valid range (16.3 27.8) at 425.00 K -13.00 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 450.00 K -15.00 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 450.00 K -14.50 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 450.00 K -14.00 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 450.00 K -13.50 GPa


    The equlibrium volume 28.8 exceed valid range (16.3 27.8) at 450.00 K -13.00 GPa


    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 475.00 K -15.00 GPa


    The equlibrium volume 31.7 exceed valid range (16.3 27.8) at 475.00 K -14.50 GPa


    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 475.00 K -14.00 GPa


    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 475.00 K -13.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 475.00 K -13.00 GPa


    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 500.00 K -15.00 GPa


    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 500.00 K -14.50 GPa


    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 500.00 K -14.00 GPa


    The equlibrium volume 31.6 exceed valid range (16.3 27.8) at 500.00 K -13.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 500.00 K -13.00 GPa


    The equlibrium volume 27.8 exceed valid range (16.3 27.8) at 500.00 K -12.50 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 525.00 K -15.00 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 525.00 K -14.50 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 525.00 K -14.00 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 525.00 K -13.50 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 525.00 K -13.00 GPa


    The equlibrium volume 27.9 exceed valid range (16.3 27.8) at 525.00 K -12.50 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 550.00 K -15.00 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 550.00 K -14.50 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 550.00 K -14.00 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 550.00 K -13.50 GPa


    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 550.00 K -13.00 GPa


    The equlibrium volume 28.0 exceed valid range (16.3 27.8) at 550.00 K -12.50 GPa


    The equlibrium volume 31.5 exceed valid range (16.3 27.8) at 575.00 K -15.00 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 575.00 K -14.50 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 575.00 K -14.00 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 575.00 K -13.50 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 575.00 K -13.00 GPa


    The equlibrium volume 28.2 exceed valid range (16.3 27.8) at 575.00 K -12.50 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -15.00 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -14.50 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -14.00 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -13.50 GPa


    The equlibrium volume 31.4 exceed valid range (16.3 27.8) at 600.00 K -13.00 GPa


    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 600.00 K -12.50 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -15.00 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -14.50 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -14.00 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -13.50 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 625.00 K -13.00 GPa


    The equlibrium volume 28.6 exceed valid range (16.3 27.8) at 625.00 K -12.50 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -15.00 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -14.50 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -14.00 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -13.50 GPa


    The equlibrium volume 31.3 exceed valid range (16.3 27.8) at 650.00 K -13.00 GPa


    The equlibrium volume 28.9 exceed valid range (16.3 27.8) at 650.00 K -12.50 GPa


    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -15.00 GPa


    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -14.50 GPa


    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -14.00 GPa


    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -13.50 GPa


    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 675.00 K -13.00 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 675.00 K -12.50 GPa


    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 700.00 K -15.00 GPa


    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 700.00 K -14.50 GPa


    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 700.00 K -14.00 GPa


    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 700.00 K -13.50 GPa


    The equlibrium volume 31.2 exceed valid range (16.3 27.8) at 700.00 K -13.00 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 700.00 K -12.50 GPa


    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -15.00 GPa


    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -14.50 GPa


    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -14.00 GPa


    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -13.50 GPa


    The equlibrium volume 31.1 exceed valid range (16.3 27.8) at 725.00 K -13.00 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 725.00 K -12.50 GPa


    The equlibrium volume 27.8 exceed valid range (16.3 27.8) at 725.00 K -12.00 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -15.00 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -14.50 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -14.00 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -13.50 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -13.00 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 750.00 K -12.50 GPa


    The equlibrium volume 28.0 exceed valid range (16.3 27.8) at 750.00 K -12.00 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -15.00 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -14.50 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -14.00 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -13.50 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -13.00 GPa


    The equlibrium volume 31.0 exceed valid range (16.3 27.8) at 775.00 K -12.50 GPa


    The equlibrium volume 28.2 exceed valid range (16.3 27.8) at 775.00 K -12.00 GPa


    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -15.00 GPa


    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -14.50 GPa


    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -14.00 GPa


    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -13.50 GPa


    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -13.00 GPa


    The equlibrium volume 30.9 exceed valid range (16.3 27.8) at 800.00 K -12.50 GPa


    The equlibrium volume 28.5 exceed valid range (16.3 27.8) at 800.00 K -12.00 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -15.00 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -14.50 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -14.00 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -13.50 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -13.00 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 825.00 K -12.50 GPa


    The equlibrium volume 28.8 exceed valid range (16.3 27.8) at 825.00 K -12.00 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -15.00 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -14.50 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -14.00 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -13.50 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -13.00 GPa


    The equlibrium volume 30.8 exceed valid range (16.3 27.8) at 850.00 K -12.50 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 850.00 K -12.00 GPa


    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -15.00 GPa


    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -14.50 GPa


    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -14.00 GPa


    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -13.50 GPa


    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -13.00 GPa


    The equlibrium volume 30.7 exceed valid range (16.3 27.8) at 875.00 K -12.50 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 875.00 K -12.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -15.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -14.50 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -14.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -13.50 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -13.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -12.50 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 900.00 K -12.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -15.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -14.50 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -14.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -13.50 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -13.00 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -12.50 GPa


    The equlibrium volume 30.6 exceed valid range (16.3 27.8) at 925.00 K -12.00 GPa


    The equlibrium volume 27.9 exceed valid range (16.3 27.8) at 925.00 K -11.50 GPa


    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -15.00 GPa


    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -14.50 GPa


    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -14.00 GPa


    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -13.50 GPa


    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -13.00 GPa


    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -12.50 GPa


    The equlibrium volume 30.5 exceed valid range (16.3 27.8) at 950.00 K -12.00 GPa


    The equlibrium volume 28.1 exceed valid range (16.3 27.8) at 950.00 K -11.50 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -15.00 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -14.50 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -14.00 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -13.50 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -13.00 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -12.50 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 975.00 K -12.00 GPa


    The equlibrium volume 28.4 exceed valid range (16.3 27.8) at 975.00 K -11.50 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -15.00 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -14.50 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -14.00 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -13.50 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -13.00 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -12.50 GPa


    The equlibrium volume 30.4 exceed valid range (16.3 27.8) at 1000.00 K -12.00 GPa


    The equlibrium volume 28.7 exceed valid range (16.3 27.8) at 1000.00 K -11.50 GPa


    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -15.00 GPa


    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -14.50 GPa


    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -14.00 GPa


    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -13.50 GPa


    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -13.00 GPa


    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -12.50 GPa


    The equlibrium volume 30.3 exceed valid range (16.3 27.8) at 1025.00 K -12.00 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1025.00 K -11.50 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -15.00 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -14.50 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -14.00 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -13.50 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -13.00 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -12.50 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -12.00 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1050.00 K -11.50 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -15.00 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -14.50 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -14.00 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1075.00 K -13.50 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -13.00 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -12.50 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -12.00 GPa


    The equlibrium volume 30.2 exceed valid range (16.3 27.8) at 1075.00 K -11.50 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -15.00 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -14.50 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -14.00 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -13.50 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -13.00 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -12.50 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -12.00 GPa


    The equlibrium volume 30.1 exceed valid range (16.3 27.8) at 1100.00 K -11.50 GPa


    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -15.00 GPa


    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -14.50 GPa


    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -14.00 GPa


    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -13.50 GPa


    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -13.00 GPa


    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -12.50 GPa


    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -12.00 GPa


    The equlibrium volume 30.0 exceed valid range (16.3 27.8) at 1125.00 K -11.50 GPa


    The equlibrium volume 28.0 exceed valid range (16.3 27.8) at 1125.00 K -11.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -15.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -14.50 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -14.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -13.50 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -13.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -12.50 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -12.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1150.00 K -11.50 GPa


    The equlibrium volume 28.3 exceed valid range (16.3 27.8) at 1150.00 K -11.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -15.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -14.50 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -14.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -13.50 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -13.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -12.50 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -12.00 GPa


    The equlibrium volume 29.9 exceed valid range (16.3 27.8) at 1175.00 K -11.50 GPa


    The equlibrium volume 28.8 exceed valid range (16.3 27.8) at 1175.00 K -11.00 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -15.00 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -14.50 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -14.00 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -13.50 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -13.00 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -12.50 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -12.00 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -11.50 GPa


    The equlibrium volume 29.8 exceed valid range (16.3 27.8) at 1200.00 K -11.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -15.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -14.50 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -14.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -13.50 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -13.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -12.50 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -12.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -11.50 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1225.00 K -11.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -15.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -14.50 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1250.00 K -14.00 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -13.50 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -13.00 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1250.00 K -12.50 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -12.00 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1250.00 K -11.50 GPa


    The equlibrium volume 29.7 exceed valid range (16.3 27.8) at 1250.00 K -11.00 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -15.00 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -14.50 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -14.00 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -13.50 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -13.00 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -12.50 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -12.00 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -11.50 GPa


    The equlibrium volume 29.6 exceed valid range (16.3 27.8) at 1275.00 K -11.00 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -15.00 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -14.50 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -14.00 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -13.50 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -13.00 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -12.50 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -12.00 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -11.50 GPa


    The equlibrium volume 29.5 exceed valid range (16.3 27.8) at 1300.00 K -11.00 GPa


    The equlibrium volume 27.9 exceed valid range (16.3 27.8) at 1300.00 K -10.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -15.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -14.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -14.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -13.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -13.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -12.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -12.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -11.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1325.00 K -11.00 GPa


    The equlibrium volume 28.3 exceed valid range (16.3 27.8) at 1325.00 K -10.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -15.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -14.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -14.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -13.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -13.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -12.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -12.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -11.50 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -11.00 GPa


    The equlibrium volume 29.4 exceed valid range (16.3 27.8) at 1350.00 K -10.50 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -15.00 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -14.50 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -14.00 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -13.50 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -13.00 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -12.50 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -12.00 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -11.50 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -11.00 GPa


    The equlibrium volume 29.3 exceed valid range (16.3 27.8) at 1375.00 K -10.50 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -15.00 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -14.50 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -14.00 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -13.50 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -13.00 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -12.50 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -12.00 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -11.50 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -11.00 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1400.00 K -10.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -15.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -14.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -14.00 GPa


    The equlibrium volume 29.2 exceed valid range (16.3 27.8) at 1425.00 K -13.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -13.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -12.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -12.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -11.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -11.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1425.00 K -10.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -15.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -14.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -14.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -13.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -13.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -12.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -12.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -11.50 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -11.00 GPa


    The equlibrium volume 29.1 exceed valid range (16.3 27.8) at 1450.00 K -10.50 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -15.00 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -14.50 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -14.00 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -13.50 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -13.00 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -12.50 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -12.00 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -11.50 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -11.00 GPa


    The equlibrium volume 29.0 exceed valid range (16.3 27.8) at 1475.00 K -10.50 GPa


      0%|          | 0/15 [00:00&lt;?, ?it/s]


    The equlibrium volume 21.1 exceed valid range (12.4 21.0) at 0.00 K -15.00 GPa


    The equlibrium volume 21.1 exceed valid range (12.4 21.0) at 25.00 K -15.00 GPa


    The equlibrium volume 21.2 exceed valid range (12.4 21.0) at 50.00 K -15.00 GPa


    The equlibrium volume 21.3 exceed valid range (12.4 21.0) at 75.00 K -15.00 GPa


    The equlibrium volume 21.4 exceed valid range (12.4 21.0) at 100.00 K -15.00 GPa


    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 125.00 K -15.00 GPa


    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 150.00 K -15.00 GPa


    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 175.00 K -15.00 GPa


    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 200.00 K -15.00 GPa


    The equlibrium volume 21.5 exceed valid range (12.4 21.0) at 225.00 K -15.00 GPa


    The equlibrium volume 21.4 exceed valid range (12.4 21.0) at 250.00 K -15.00 GPa


    The equlibrium volume 21.3 exceed valid range (12.4 21.0) at 275.00 K -15.00 GPa


    The equlibrium volume 21.2 exceed valid range (12.4 21.0) at 300.00 K -15.00 GPa


    The equlibrium volume 21.1 exceed valid range (12.4 21.0) at 325.00 K -15.00 GPa


      0%|          | 0/17 [00:00&lt;?, ?it/s]


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 0.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 25.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 50.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 75.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 100.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 125.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 150.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 175.00 K 19.00 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 175.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 200.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 200.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 225.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 225.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 250.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 250.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 275.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 275.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 300.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 300.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 325.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 325.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 350.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 350.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 375.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 375.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 400.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 400.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 425.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 425.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 450.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 450.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 475.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 475.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 500.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 500.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 525.00 K 18.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 525.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 525.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 550.00 K 18.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 550.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 550.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 575.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 575.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 575.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 600.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 600.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 600.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 625.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 625.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 625.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 650.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 650.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 650.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 675.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 675.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 675.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 700.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 700.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 700.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 725.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 725.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 725.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 750.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 750.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 750.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 775.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 775.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 775.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 800.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 800.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 800.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 825.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 825.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 825.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 850.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 850.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 850.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 875.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 875.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 875.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 900.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 900.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 900.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 925.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 925.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 925.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 950.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 950.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 950.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 975.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 975.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 975.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1000.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1000.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1000.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1025.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1025.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1025.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1050.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1050.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1050.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1075.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1075.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1075.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1100.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1100.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1100.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1125.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1125.00 K 19.00 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1125.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1150.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1150.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1150.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1175.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1175.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1175.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1200.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1200.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1200.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1225.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1225.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1225.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1250.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1250.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1250.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1275.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1275.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1275.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1300.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1300.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1300.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1325.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1325.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1325.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1350.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1350.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1350.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1375.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1375.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1375.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1400.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1400.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1400.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1425.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1425.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1425.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1450.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1450.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1450.00 K 19.50 GPa


    The equlibrium volume 20.0 exceed valid range (20.0 35.3) at 1475.00 K 18.50 GPa


    The equlibrium volume 19.9 exceed valid range (20.0 35.3) at 1475.00 K 19.00 GPa


    The equlibrium volume 19.8 exceed valid range (20.0 35.3) at 1475.00 K 19.50 GPa


      0%|          | 0/13 [00:00&lt;?, ?it/s]


    The equlibrium volume 22.0 exceed valid range (12.1 19.4) at 0.00 K -15.00 GPa


    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 0.00 K -14.50 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 0.00 K -14.00 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 0.00 K -13.50 GPa


    The equlibrium volume 22.1 exceed valid range (12.1 19.4) at 25.00 K -15.00 GPa


    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 25.00 K -14.50 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 25.00 K -14.00 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 25.00 K -13.50 GPa


    The equlibrium volume 22.7 exceed valid range (12.1 19.4) at 50.00 K -15.00 GPa


    The equlibrium volume 20.9 exceed valid range (12.1 19.4) at 50.00 K -14.50 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 50.00 K -14.00 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 50.00 K -13.50 GPa


    The equlibrium volume 23.0 exceed valid range (12.1 19.4) at 75.00 K -15.00 GPa


    The equlibrium volume 21.0 exceed valid range (12.1 19.4) at 75.00 K -14.50 GPa


    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 75.00 K -14.00 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 75.00 K -13.50 GPa


    The equlibrium volume 23.1 exceed valid range (12.1 19.4) at 100.00 K -15.00 GPa


    The equlibrium volume 21.2 exceed valid range (12.1 19.4) at 100.00 K -14.50 GPa


    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 100.00 K -14.00 GPa


    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 100.00 K -13.50 GPa


    The equlibrium volume 23.1 exceed valid range (12.1 19.4) at 125.00 K -15.00 GPa


    The equlibrium volume 21.4 exceed valid range (12.1 19.4) at 125.00 K -14.50 GPa


    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 125.00 K -14.00 GPa


    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 125.00 K -13.50 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 125.00 K -13.00 GPa


    The equlibrium volume 23.2 exceed valid range (12.1 19.4) at 150.00 K -15.00 GPa


    The equlibrium volume 21.6 exceed valid range (12.1 19.4) at 150.00 K -14.50 GPa


    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 150.00 K -14.00 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 150.00 K -13.50 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 150.00 K -13.00 GPa


    The equlibrium volume 23.3 exceed valid range (12.1 19.4) at 175.00 K -15.00 GPa


    The equlibrium volume 21.7 exceed valid range (12.1 19.4) at 175.00 K -14.50 GPa


    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 175.00 K -14.00 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 175.00 K -13.50 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 175.00 K -13.00 GPa


    The equlibrium volume 23.4 exceed valid range (12.1 19.4) at 200.00 K -15.00 GPa


    The equlibrium volume 21.7 exceed valid range (12.1 19.4) at 200.00 K -14.50 GPa


    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 200.00 K -14.00 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 200.00 K -13.50 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 200.00 K -13.00 GPa


    The equlibrium volume 23.5 exceed valid range (12.1 19.4) at 225.00 K -15.00 GPa


    The equlibrium volume 21.6 exceed valid range (12.1 19.4) at 225.00 K -14.50 GPa


    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 225.00 K -14.00 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 225.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 225.00 K -13.00 GPa


    The equlibrium volume 23.6 exceed valid range (12.1 19.4) at 250.00 K -15.00 GPa


    The equlibrium volume 21.6 exceed valid range (12.1 19.4) at 250.00 K -14.50 GPa


    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 250.00 K -14.00 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 250.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 250.00 K -13.00 GPa


    The equlibrium volume 23.8 exceed valid range (12.1 19.4) at 275.00 K -15.00 GPa


    The equlibrium volume 21.5 exceed valid range (12.1 19.4) at 275.00 K -14.50 GPa


    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 275.00 K -14.00 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 275.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 275.00 K -13.00 GPa


    The equlibrium volume 22.8 exceed valid range (12.1 19.4) at 300.00 K -15.00 GPa


    The equlibrium volume 21.4 exceed valid range (12.1 19.4) at 300.00 K -14.50 GPa


    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 300.00 K -14.00 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 300.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 300.00 K -13.00 GPa


    The equlibrium volume 22.4 exceed valid range (12.1 19.4) at 325.00 K -15.00 GPa


    The equlibrium volume 21.3 exceed valid range (12.1 19.4) at 325.00 K -14.50 GPa


    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 325.00 K -14.00 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 325.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 325.00 K -13.00 GPa


    The equlibrium volume 22.1 exceed valid range (12.1 19.4) at 350.00 K -15.00 GPa


    The equlibrium volume 21.1 exceed valid range (12.1 19.4) at 350.00 K -14.50 GPa


    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 350.00 K -14.00 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 350.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 350.00 K -13.00 GPa


    The equlibrium volume 21.8 exceed valid range (12.1 19.4) at 375.00 K -15.00 GPa


    The equlibrium volume 21.0 exceed valid range (12.1 19.4) at 375.00 K -14.50 GPa


    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 375.00 K -14.00 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 375.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 375.00 K -13.00 GPa


    The equlibrium volume 21.6 exceed valid range (12.1 19.4) at 400.00 K -15.00 GPa


    The equlibrium volume 20.9 exceed valid range (12.1 19.4) at 400.00 K -14.50 GPa


    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 400.00 K -14.00 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 400.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 400.00 K -13.00 GPa


    The equlibrium volume 21.5 exceed valid range (12.1 19.4) at 425.00 K -15.00 GPa


    The equlibrium volume 20.8 exceed valid range (12.1 19.4) at 425.00 K -14.50 GPa


    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 425.00 K -14.00 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 425.00 K -13.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 425.00 K -13.00 GPa


    The equlibrium volume 21.3 exceed valid range (12.1 19.4) at 450.00 K -15.00 GPa


    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 450.00 K -14.50 GPa


    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 450.00 K -14.00 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 450.00 K -13.50 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 450.00 K -13.00 GPa


    The equlibrium volume 21.1 exceed valid range (12.1 19.4) at 475.00 K -15.00 GPa


    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 475.00 K -14.50 GPa


    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 475.00 K -14.00 GPa


    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 475.00 K -13.50 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 475.00 K -13.00 GPa


    The equlibrium volume 21.0 exceed valid range (12.1 19.4) at 500.00 K -15.00 GPa


    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 500.00 K -14.50 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 500.00 K -14.00 GPa


    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 500.00 K -13.50 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 500.00 K -13.00 GPa


    The equlibrium volume 20.9 exceed valid range (12.1 19.4) at 525.00 K -15.00 GPa


    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 525.00 K -14.50 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 525.00 K -14.00 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 525.00 K -13.50 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 525.00 K -13.00 GPa


    The equlibrium volume 20.8 exceed valid range (12.1 19.4) at 550.00 K -15.00 GPa


    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 550.00 K -14.50 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 550.00 K -14.00 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 550.00 K -13.50 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 550.00 K -13.00 GPa


    The equlibrium volume 20.7 exceed valid range (12.1 19.4) at 575.00 K -15.00 GPa


    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 575.00 K -14.50 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 575.00 K -14.00 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 575.00 K -13.50 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 575.00 K -13.00 GPa


    The equlibrium volume 20.6 exceed valid range (12.1 19.4) at 600.00 K -15.00 GPa


    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 600.00 K -14.50 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 600.00 K -14.00 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 600.00 K -13.50 GPa


    The equlibrium volume 20.5 exceed valid range (12.1 19.4) at 625.00 K -15.00 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 625.00 K -14.50 GPa


    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 625.00 K -14.00 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 625.00 K -13.50 GPa


    The equlibrium volume 20.4 exceed valid range (12.1 19.4) at 650.00 K -15.00 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 650.00 K -14.50 GPa


    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 650.00 K -14.00 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 650.00 K -13.50 GPa


    The equlibrium volume 20.3 exceed valid range (12.1 19.4) at 675.00 K -15.00 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 675.00 K -14.50 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 675.00 K -14.00 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 675.00 K -13.50 GPa


    The equlibrium volume 20.2 exceed valid range (12.1 19.4) at 700.00 K -15.00 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 700.00 K -14.50 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 700.00 K -14.00 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 700.00 K -13.50 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 725.00 K -15.00 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 725.00 K -14.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 725.00 K -14.00 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 725.00 K -13.50 GPa


    The equlibrium volume 20.1 exceed valid range (12.1 19.4) at 750.00 K -15.00 GPa


    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 750.00 K -14.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 750.00 K -14.00 GPa


    The equlibrium volume 20.0 exceed valid range (12.1 19.4) at 775.00 K -15.00 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 775.00 K -14.50 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 775.00 K -14.00 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 800.00 K -15.00 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 800.00 K -14.50 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 800.00 K -14.00 GPa


    The equlibrium volume 19.9 exceed valid range (12.1 19.4) at 825.00 K -15.00 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 825.00 K -14.50 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 825.00 K -14.00 GPa


    The equlibrium volume 19.8 exceed valid range (12.1 19.4) at 850.00 K -15.00 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 850.00 K -14.50 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 850.00 K -14.00 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 875.00 K -15.00 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 875.00 K -14.50 GPa


    The equlibrium volume 19.7 exceed valid range (12.1 19.4) at 900.00 K -15.00 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 900.00 K -14.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 925.00 K -15.00 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 925.00 K -14.50 GPa


    The equlibrium volume 19.6 exceed valid range (12.1 19.4) at 950.00 K -15.00 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 950.00 K -14.50 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 975.00 K -15.00 GPa


    The equlibrium volume 19.5 exceed valid range (12.1 19.4) at 1000.00 K -15.00 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 1025.00 K -15.00 GPa


    The equlibrium volume 19.4 exceed valid range (12.1 19.4) at 1050.00 K -15.00 GPa


Gibbs自由エネルギーの計算結果から、特定の温度・圧力下でGibbs自由エネルギーが最小になる結晶構造を知ることができます。これにより、相図を描くことができます。`matlantis-features`では相図を描くための関数 `draw_phase_diagram` を提供しています。


```python
from matlantis_features.features.phonon.qha_gibbs_free_energy import draw_phase_diagram

fig, ax = draw_phase_diagram(label, gibbs_results)
```



![png](output_34_0.png)



### 実験値と計算結果の比較
文献から相転移点における温度・圧力を参照します。
* $\triangledown$: ダイヤモンド相 -&gt; $\beta$ -Sn相の境界 [1].
* $\vartriangle$: ダイヤモンド相 -&gt; $\beta$ -Sn相の境界 [2].
* $\times$: 第一原理計算によるクラスレート層 -&gt; ダイヤモンド相の境界 (約 -3GPa)とダイヤモンド層 -&gt; $\beta$ -Sn相の境界 [3].
* 水色の破線: 第一原理計算によるダイヤモンド相 -&gt; $\beta$ -Sn相の境界 [4].
* 紫色の破線: 第一原理計算による $\beta$ -Sn相 -&gt; sh相の境界 [4].

[1] G. A. Voronin, C. Pantea, T. W. Zerda, L. Wang and Y. Zhao, Phys. Rev. B. 68, 020102 (2003).

[2] J. Z. Hu, L. D. Merkle, C. S. Menoni and I. L. Spain, Phys. Rev. B 34, 4679 (1986).

[3] M. Kaczmarski, O. N. Bedoya-Martínez, and E. R. Hernández, Phys. Rev. Lett. 94, 095701 (2005).

[4] C. Li, C. P. Wang, J. J. Han, L. H. Yan, B. Deng &amp; X. J. Liu, J Mater. Sci. 53, 7475 (2018)


```python
# G. A. Voronin, C. Pantea, T. W. Zerda, L. Wang and Y. Zhao, Phys. Rev. B. 68, 020102 (2003).
expt = np.array([
    [10.50, 973],
    [10.60, 773],
    [10.40, 573],
    [10.50, 273],
    [10.80, 973],
    [11.70, 773],
    [12.40, 573],
    [13.40, 273],
    [10.50, 1003],
])

#J. Z. Hu, L. D. Merkle, C. S. Menoni and I. L. Spain, Phys. Rev. B 34, 4679 (1986).
expt2 = np.array([
    [11.40, 77],
    [12.36, 77],
])

# M. Kaczmarski, O. N. Bedoya-Martínez, and E. R. Hernández, Phys. Rev. Lett. 94, 095701 (2005).
calc = np.array([
    [-2.45, 500],
    [-2.38, 1000],
    [12.62, 1000],
    [15.48, 500],
    [10.0, 1231],
])

# C. Li, C. P. Wang, J. J. Han, L. H. Yan, B. Deng &amp; X. J. Liu, J Mater. Sci. 53, 7475 (2018)
calc2 = np.array([
    [10.2, 1362],
    [10.4, 1235],
    [10.6, 1097],
    [10.8, 1005],
    [11.0, 630],
    [11.2, 554],
    [11.4, 438],
    [11.6, 243],
    [11.8, 138],
    [12, 59],
])

calc3 = np.array([
    [15.7, 46],
    [15.0, 158],
    [14.0, 312],
    [13.2, 496],
    [12.6, 657],
    [12.3, 817],
    [12.3, 974],
    [12.8, 1121],
    [13.7, 1216],
    [15.0, 1295],
    [16.6, 1339],
])
```


```python
ax.scatter(expt[:,0], expt[:,1], c="k", marker = "v", s = 15)
ax.scatter(expt2[:,0], expt2[:,1], c="k", marker = "^", s=15)
ax.scatter(calc[:,0], calc[:,1], c="k", marker = "x", s=15)
ax.plot(calc2[:,0], calc2[:,1], linestyle="--", linewidth=2, c="c")
ax.plot(calc3[:,0], calc3[:,1], linestyle="--", linewidth=2, c="m")
fig
```



![png](output_37_0.png)



計算結果はクラスレート相は低圧力帯、 sh相は高圧力帯に存在すると示していますが、これらは実験値や第一原理計算といくつかの違いがあります。
* クラスレート相 -&gt; ダイヤモンド相の転移がある圧力について、このnotebookでの計算結果(\~-10GPa)は第一原理計算の結果(\~-3GPa)よりも小さい値になっています。
* ダイヤモンド相 -&gt; sh相の転移がある圧力は実験値より少し小さい値になっています。
* 第一原理計算によると、低温度(0~800K)・高圧力(&gt;10GPa)ではダイヤモンド相は$\beta$-Sn相を経てsh相へ転移していきますが、このnotebookの計算結果では高圧力帯に$\beta$-Sn相が存在せず、直接sh相へ転移します。
