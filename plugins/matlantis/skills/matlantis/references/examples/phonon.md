# Matlantis-features: フォノン


## 概要

* フォノンは、結晶構造における格子振動に対応する量子です。古典的には、結晶中の弾性波に対応する量です。フォノンは、原子の振動が調和振動的であるという近似に基づいて計算されます。この近似は温度が低いときに特に有効で、結晶に対する種々の熱力学特性を効率的に計算することができます。

* フォノンの計算は、原子の位置の変化に対する力の応答を表現する力定数マトリクスという量を求めることによって行われます。これは実質的にエネルギー曲面の2階微分に相当する量を計算していることになります。

## Matlantis-featuresにおけるフォノン計算

* Matlantis-featuresが提供するフィーチャーには以下のものがあります。
    * `ForceConstantFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/force_constant.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;)
    * `PostPhononBandFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/phonon_band.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;)
    * `PostPhononDOSFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/phonon_dos.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;)
    * `PostPhononModeFeature`
    * `PostPhononThermochemistryFeature`


* このうち、フォノンおよびフォノンから計算される種々の量を求めるのに必要なのは力定数マトリクスで、これは`ForceConstantFeature`によって得ることができます。それ以外の計算は力定数マトリクスを使ったポスト処理と考えることができます。

* 以下の例では、シリコンのフォノン計算を行います。はじめに力定数マトリクスを求めた後、バンド構造、状態密度、熱に関する物性値を計算します。最後に実験値およびDFT計算の値との比較を掲載しています。


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


### `estimator_fn` による計算モードとモデルバージョンの指定

* Feature に `estimator_fn` 引数を使うことで、matlantis-features で使われる PFP の計算モードとモデルバージョンを指定することができます。

* `estimator_fn` は、Estimator オブジェクトを作る factory method です。計算モードとモデルバージョンのみを指定する場合は、以下のように `pfp_estimator_fn` を利用できます。より細かい指定が必要な場合は、自身で factory method を定義できます。詳細に関しては [Matlantis Guidebook](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html#estimator-fnpfp) をご参照ください。

* `estimator_fn` が指定されない場合は、環境変数で指定された値が使用されます。`estimator_fn` も環境変数も指定されない場合は、デフォルトのモデルバージョンと計算モード（`PBE`）が使用されます。

* 環境変数と `estimator_fn` が同時に定義されている場合、 `estimator_fn` の設定が優先して使用されます。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode

estimator_fn = pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE_PLUS_D3)
```


```python
import logging
import pathlib
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from ase.io import read
from plotly.offline import iplot

from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.common import FireASEOptFeature
from matlantis_features.features.phonon import (
    ForceConstantFeature,
    PostPhononBandFeature,
    PostPhononDOSFeature,
    PostPhononModeFeature,
    PostPhononThermochemistryFeature,
    plot_band_dos,
)

try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

```


## 0. 材料の準備

* 構造ファイルを読み込み、構造緩和計算によって安定構造を得ます。なお、Si (mp-149)の構造ファイルは [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) の下で Materials Project から入手したものです。


```python
atoms_Si = MatlantisAtoms.from_file(str(dir_path/"assets/phonon/Si-mp-149.xyz"))

opt = FireASEOptFeature(fmax=0.0001, filter=True, show_progress_bar=True, estimator_fn=estimator_fn)

atoms_opt = opt(atoms_Si).atoms
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


## 1. 力定数マトリクスの計算

* 力定数マトリクスの計算には、k点(`supercell`)および原子の微小変形に相当する`delta`の設定が必要です。k点は多いほど正確ですが、計算時間が増大します。`delta`値は調和近似が成立する範囲で計算誤差を減らしたいという視点から、適切な数値に設定する必要があります。現在のPFPでは、出力の細かいノイズの影響を抑制するため、比較的大きいdelta値(0.1)を使うことが望ましいです。


```python
fc = ForceConstantFeature(
    supercell = (10,10,10),
    delta = 0.1,
    show_logger=True,
    show_progress_bar=True,
    estimator_fn=estimator_fn,
)
force_constant = fc(atoms_opt)

print(force_constant.force_constant.shape)
print(np.round(force_constant.force_constant[0], 2))

```


      0%|          | 0/13 [00:00&lt;?, ?it/s]


    Force calculation (no displacement)
    Force calculation (0th atoms displacement by [-0.100  0.000  0.000])
    Force calculation (0th atoms displacement by [ 0.100  0.000  0.000])
    Force calculation (0th atoms displacement by [ 0.000 -0.100  0.000])
    Force calculation (0th atoms displacement by [ 0.000  0.100  0.000])
    Force calculation (0th atoms displacement by [ 0.000  0.000 -0.100])
    Force calculation (0th atoms displacement by [ 0.000  0.000  0.100])
    Force calculation (1th atoms displacement by [-0.100  0.000  0.000])
    Force calculation (1th atoms displacement by [ 0.100  0.000  0.000])
    Force calculation (1th atoms displacement by [ 0.000 -0.100  0.000])
    Force calculation (1th atoms displacement by [ 0.000  0.100  0.000])
    Force calculation (1th atoms displacement by [ 0.000  0.000 -0.100])
    Force calculation (1th atoms displacement by [ 0.000  0.000  0.100])


    (1000, 6, 6)
    [[13.41 -0.    0.   -3.28 -2.33 -2.33]
     [-0.   13.41  0.   -2.33 -3.28 -2.33]
     [ 0.    0.   13.41 -2.33 -2.33 -3.28]
     [-3.28 -2.33 -2.33 13.41 -0.   -0.  ]
     [-2.33 -3.28 -2.33 -0.   13.41  0.  ]
     [-2.33 -2.33 -3.28 -0.    0.   13.41]]


    Force calculation (0th atoms displacement by [ 0.000  0.000  0.100])


    Force calculation (1th atoms displacement by [-0.100  0.000  0.000])


    Force calculation (1th atoms displacement by [ 0.100  0.000  0.000])


    Force calculation (1th atoms displacement by [ 0.000 -0.100  0.000])


    Force calculation (1th atoms displacement by [ 0.000  0.100  0.000])


    Force calculation (1th atoms displacement by [ 0.000  0.000 -0.100])


    Force calculation (1th atoms displacement by [ 0.000  0.000  0.100])


    (1000, 6, 6)
    [[13.41 -0.   -0.   -3.28 -2.33 -2.33]
     [-0.   13.41  0.   -2.33 -3.28 -2.33]
     [-0.    0.   13.41 -2.33 -2.33 -3.28]
     [-3.28 -2.33 -2.33 13.41  0.   -0.  ]
     [-2.33 -3.28 -2.33  0.   13.41 -0.  ]
     [-2.33 -2.33 -3.28 -0.   -0.   13.41]]


## 2. バンド構造の計算
* labels と special_kpts が指定されていない場合、kpoint path は結晶構造に基づいて自動的に生成されます。


```python
band = PostPhononBandFeature()
band_results = band(
    force_constant,
    labels = ['Gamma', 'X', 'W', 'K', 'Gamma', 'L', 'U', 'W', 'L', 'K', '|', 'U', 'X'],
    special_kpts = np.array([
        [0.000, 0.000, 0.000], # Gamma
        [0.500, 0.000, 0.500], # X
        [0.500, 0.250, 0.750], # W
        [0.375, 0.375, 0.750], # K
        [0.000, 0.000, 0.000], # Gamma
        [0.500, 0.500, 0.500], # L
        [0.625, 0.250, 0.625], # U
        [0.500, 0.250, 0.750], # W
        [0.500, 0.500, 0.500], # L
        [0.375, 0.375, 0.75 ], # K
        [0.625, 0.250, 0.625], # U
        [0.500, 0.000, 0.500], # X
    ]),
    total_n_kpts = 100,
    unit = "meV"
)
```


* バンド構造とk点パスを可視化してみます。

* 初回起動時などで図が表示されない場合、ブラウザでページをリロードすると表示されることがあります。


```python
from ipywidgets import HBox
dispersion = band_results.plot()
dispersion.update_layout(width=600, height=500)

brillouin = band_results.visual_k_path()
brillouin.update_layout(width=600, height=500)

HBox([go.FigureWidget(dispersion), go.FigureWidget(brillouin)])

```


    HBox(children=(FigureWidget({
        'data': [{'marker': {'color': '#636EFA'},
                  'mode': 'lines',
       …


* それぞれのk点におけるフォノンを確認してみましょう。


```python
df = band_results.to_dataframe()
print(df.round(2).loc[:20])
```

          kx    ky    kz  coords labels      0      1      2      3      4      5
    0   0.00  0.00  0.00    0.00  Gamma  -0.00  -0.00  -0.00  62.96  62.96  62.96
    1   0.03  0.00  0.03    0.04   None   2.58   2.70   3.92  62.86  62.86  62.92
    2   0.06  0.00  0.06    0.09   None   5.10   5.21   7.79  62.58  62.58  62.80
    3   0.09  0.00  0.09    0.13   None   7.51   7.53  11.59  62.14  62.14  62.62
    4   0.12  0.00  0.12    0.18   None   9.70   9.74  15.32  61.56  61.57  62.36
    5   0.16  0.00  0.16    0.22   None  11.70  11.72  19.01  60.91  60.92  62.02
    6   0.19  0.00  0.19    0.27   None  13.39  13.39  22.64  60.23  60.23  61.60
    7   0.22  0.00  0.22    0.31   None  14.69  14.69  26.15  59.58  59.58  61.08
    8   0.25  0.00  0.25    0.35   None  15.60  15.60  29.49  59.01  59.01  60.44
    9   0.28  0.00  0.28    0.40   None  16.14  16.15  32.64  58.52  58.52  59.67
    10  0.31  0.00  0.31    0.44   None  16.38  16.38  35.60  58.14  58.14  58.77
    11  0.34  0.00  0.34    0.49   None  16.39  16.39  38.38  57.72  57.84  57.84
    12  0.38  0.00  0.38    0.53   None  16.27  16.27  41.00  56.53  57.61  57.61
    13  0.41  0.00  0.41    0.57   None  16.12  16.12  43.47  55.17  57.43  57.43
    14  0.44  0.00  0.44    0.62   None  15.99  15.99  45.81  53.63  57.30  57.30
    15  0.47  0.00  0.47    0.66   None  15.90  15.90  48.01  51.92  57.22  57.22
    16  0.50  0.00  0.50    0.71      X  15.87  15.87  50.05  50.05  57.20  57.20
    17  0.50  0.03  0.53    0.75   None  16.18  16.18  49.82  49.82  57.26  57.26
    18  0.50  0.06  0.56    0.80   None  17.03  17.03  49.19  49.19  57.43  57.43
    19  0.50  0.09  0.59    0.84   None  18.30  18.30  48.27  48.27  57.63  57.63
    20  0.50  0.12  0.62    0.88   None  19.78  19.78  47.18  47.18  57.82  57.82


## 3. 状態密度(Density of States, DoS)の計算


```python
dos = PostPhononDOSFeature()
dos_results = dos(
    force_constant,
    kpts = [20,20,20],
    freq_min = -5.0,
    freq_max = 70.0,
    freq_bin =  1.0,
    scheme = "mp",
    unit = "meV"
)
```


* 同様に、可視化してみます。


```python
fig = dos_results.plot()
fig.update_layout(width=800, height=500)
iplot(fig)
```


* それぞれの周波数におけるDoSを見てみましょう。


```python
df = dos_results.to_dataframe()
print(df.loc[:20])
```

        energy     total
    0     -4.5  0.000000
    1     -3.5  0.000000
    2     -2.5  0.000000
    3     -1.5  0.000000
    4     -0.5  0.000000
    5      0.5  0.000000
    6      1.5  0.000083
    7      2.5  0.000042
    8      3.5  0.000375
    9      4.5  0.000208
    10     5.5  0.001125
    11     6.5  0.000250
    12     7.5  0.001958
    13     8.5  0.002042
    14     9.5  0.002708
    15    10.5  0.003375
    16    11.5  0.004083
    17    12.5  0.007708
    18    13.5  0.020208
    19    14.5  0.021542
    20    15.5  0.026000


* フォノンの band structure と DOS を同じ図に表示することができます。


```python
fig = plot_band_dos(band_results, dos_results)
fig.update_layout(width=1000, height=500)
iplot(fig)
```


## 4. 熱力学特性の計算


```python
thermo = PostPhononThermochemistryFeature()
thermo_results = thermo(
    force_constant,
    kpts = [20, 20, 20],
    tmin = 0.0,
    tmax = 1000.0,
    tstep = 10.0,
    scheme = "mp",
)
```


```python
fig = thermo_results.plot()
fig.update_layout(width=800, height=500)
iplot(fig)
```


```python
df = thermo_results.to_dataframe()
print(df[["temperature (K)", "specific_heat (J/K/g)", "entropy (J/K/g)", "internal_energy (J/g)", "helmholtz_free_energy (J/g)"]].loc[::10])
```

         temperature (K)  specific_heat (J/K/g)  entropy (J/K/g)  \
    0                0.0               0.000033         0.000008
    10             100.0               0.274813         0.153582
    20             200.0               0.562209         0.440448
    30             300.0               0.708109         0.699748
    40             400.0               0.777787         0.914147
    50             500.0               0.814501         1.092048
    60             600.0               0.835798         1.242593
    70             700.0               0.849133         1.372507
    80             800.0               0.857996         1.486509
    90             900.0               0.864168         1.587944
    100           1000.0               0.868633         1.679236

         internal_energy (J/g)  helmholtz_free_energy (J/g)
    0               210.336332                   210.336293
    10              220.852624                   205.494424
    20              263.916174                   175.826657
    30              328.359370                   118.434882
    40              403.051875                    37.393148
    50              482.844660                   -63.179091
    60              565.448135                  -180.107438
    70              649.742809                  -311.012152
    80              735.127411                  -454.079942
    90              821.253106                  -607.896065
    100             907.904583                  -771.331258


## 5. 実験値およびDFT計算値との比較


### 5-1. フォノンDoSのDFT計算値との比較
* DFT計算の値は、Materials Projectを参照しています。(https://materialsproject.org/materials/mp-149/)


```python
dos_expt = pd.read_csv(f"{dir_path}/assets/phonon/dos.csv")
dos_expt["frequencies"] /= 0.24180 # change THz to meV

fig = dos_results.plot()
fig.add_trace(
    go.Scatter(
        x = dos_expt["frequencies"],
        y = [v*0.045 for v in dos_expt["densities"]],
        name = "DFT",
    )
)
fig.update_layout(width=800, height=500)
iplot(fig)
```


### 5-2. 熱力学特性の実験値およびDFT計算値の比較
* DFT計算の値は、同様にMaterials Projectを参照しています。 (https://materialsproject.org/materials/mp-149/)
* 実験値は、NISTデータベースを参照しています。 (https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440213&amp;Units=SI&amp;Mask=2#Thermo-Condensed)


```python
thermo_dft = pd.read_csv(f"{dir_path}/assets/phonon/thermo.csv")
thermo_expt = pd.read_csv(f"{dir_path}/assets/phonon/thermo_expt.csv")
```


```python
fig = thermo_results.plot()
fig.add_trace(
    go.Scatter(
        x = thermo_dft["T"].loc[:100],
        y = thermo_dft["entropy"].loc[:100] / 28.09 / 2.,
        name = "S DFT(J/K/g)"
    )
)
fig.add_trace(
    go.Scatter(
        x = thermo_dft["T"].loc[:100],
        y = thermo_dft["cv"].loc[:100] / 28.09 / 2.,
        name = "Cv DFT(J/K/g)"
    )
)
fig.add_trace(
    go.Scatter(
        x = thermo_dft["T"].loc[:100],
        y = thermo_dft["internal_energy"].loc[:100] / 28.09 / 2./1000.,
        name = "U DFT (kJ/g)"
    )
)
fig.add_trace(
    go.Scatter(
        x = thermo_dft["T"].loc[:100],
        y = thermo_dft["helmholtz_free_energy"].loc[:100] / 28.09 / 2./1000.,
        name = "F DFT(kJ/g)"
    )
)
fig.add_trace(
    go.Scatter(
        x = thermo_expt["T"],
        y = thermo_expt["cp"],
        mode = "markers",
        name = "Cp expt(J/k/g)"
    )
)
fig.add_trace(
    go.Scatter(
        x = thermo_expt["T"],
        y = thermo_expt["entropy"],
        mode = "markers",
        name = "S expt(J/k/g)"
    )
)
fig.update_layout(width=800, height=500)

```


