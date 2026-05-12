# Matlantis-features: 弾性率


## このexampleについて

* この資料は弾性率計算についてのexampleですが、同時に、matlantis-featuresを利用した計算の流れや一般的な使い方について実例を通して紹介しています。

## 概要

### 弾性率について

* 体積弾性率、ヤング率、ポアソン比をはじめとした種々の弾性特性は、外部応力を受けたときの材料のひずみを表す量であり、力学的特性や熱的特性と密接な関係があります。

* これらの弾性特性は、弾性テンソルという物理量から計算することができます。弾性テンソルの計算は、ひずみを静的に与えたときの応力の応答を使って求めます。計算の流れは以下のようになります。
    1) 入力された構造体の原子位置とセル形状の両方を緩和します。
    2) 格子ベクトルにひずみを与えます。
    3) セル形状を固定したまま原子位置を緩和し、変形による応力を計算します。
    4) 歪みと応力の線形関係をフィッティングすることで、6×6の弾性テンソルを計算します。


* 計算した弾性テンソルから、様々な弾性特性に換算することができます。
    * 以下の物性値を計算できます。
        * 体積弾性率 (bulk_modulus)
        * 剛性率 (shear_modulus)
        * ヤング率 (Young’s_modulus)
        * ポアソン比 (poisson_ratio)
        * 異方性 (anisotropy)
    * また、特定の方向に対する弾性特性も計算できます。

### Matlantis-featuresにおける弾性率計算

* matlantis-featuresが提供する、弾性率に関連したfeatureには以下のものがあります。
    * `matlantis_features.features.elasticity.ElasticFeature`
    * `matlantis_features.features.elasticity.ElasticTensorFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/elastic_tensor.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;)
    * `matlantis_features.features.elasticity.PostElasticPropertiesFeature`


* このうち、下2つの`ElasticTensorFeature`および`PostElasticPropertiesFeature`は、弾性率を計算するための基本的なfeatureです。`ElasticTensorFeature`はある結晶構造を入力として、応力とひずみの関係式から弾性テンソルを計算します。`PostElasticPropertiesFeature`は弾性テンソルを入力として種々の弾性特性を計算します。

* 一方、`ElasticFeature`はこれらのfeatureを組み合わせた計算を自動で実行するfeatureです。もし弾性率について標準的な計算のみを行うのであれば、ユーザーは`ElasticFeature`のみを使うことで目的が達成できます。一方で、独自の計算を行うために弾性テンソルの計算そのものや弾性特性を自前の実装で計算したい場合、`ElasticTensorFeature`や`PostElasticPropertiesFeature`のいずれかにあたる処理を自分で書き、それ以外の部分では既存のfeatureを使う、という方法もとれます。

* 以下の計算例ではSiおよびAl2O3の弾性率の計算を通して、それぞれのfeatureの使い方を紹介します。まず`ElasticTensorFeature`および`PostElasticPropertiesFeature`を順に使用した例、次に`ElasticFeature`を使用した例を紹介します。

* 計算の詳細などの詳しい情報については &lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/elastic_tensor.html" target="_blank"&gt;GuidebookのElasticTensorFeatureのページ&lt;/a&gt; もご参照ください。


## 初期設定

* はじめに、必要なライブラリ等の準備を行います。

### matlantis-featuresのインストール

* matlantis-featuresをはじめて使うときは、matlantis-featuresそのものをインストールする必要があります。

* これを行うには、下の`!pip install 'matlantis-features&gt;=0.8.1'`と書いてあるセルを実行してください。初回のインストールには数分程度かかります。


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

* feature に `estimator_fn` 引数を使うことで、matlantis-features で使われる PFP の計算モードとモデルバージョンを指定することができます。

* `estimator_fn` は、Estimator オブジェクトを作る factory method です。計算モードとモデルバージョンのみを指定する場合は、以下のように `pfp_estimator_fn` を利用できます。より細かい指定が必要な場合は、自身で factory method を定義できます。詳細に関しては [Matlantis Guidebook](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html#estimator-fnpfp) をご参照ください。

* `estimator_fn` を指定しない場合は、環境変数で指定された値が使用されます。`estimator_fn` も環境変数も指定されない場合は、デフォルトのモデルバージョンと計算モード（`PBE`）が使用されます。

* 環境変数と `estimator_fn` が同時に定義されている場合、 `estimator_fn` の設定が優先して使用されます。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode

estimator_fn = pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
```


### ライブラリのインポート

* matlantis-featuresからは、上記で紹介したfeatureの他に`FireASEOptFeature`をインポートしています。

* `FireASEOptFeature`は構造緩和計算に使うオプティマイザに相当するfeatureで、弾性テンソルを計算する`ElasticTensorFeature`の引数として使います。


```python
import logging
import pathlib
import numpy as np
from ase.io import read
from matlantis_features.atoms import MatlantisAtoms

from matlantis_features.features.elasticity.complex_elastic_feature import ComplexElasticFeature
from matlantis_features.features.common import FireASEOptFeature
from matlantis_features.features.elasticity.elastic_properties import PostElasticPropertiesFeature
from matlantis_features.features.elasticity.elastic_tensor import (
    ElasticTensorFeature,
)

try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

```


## 原子構造の用意
* [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) の下で Materials Project から入手した Si (mp-149) および Al2O3 (mp-1143) の構造ファイルを読み込みます。


```python
atoms_Si = MatlantisAtoms.from_file(str(dir_path/"assets/elastic/Si-mp-149.xyz"))
atoms_Al2O3 = MatlantisAtoms.from_file(str(dir_path/"assets/elastic/Al2O3-mp-1143.xyz"))
```


## `ElasticTensorFeature` および `PostElasticPropertiesFeature` を利用した計算

* 弾性テンソルを計算する際のセルの変形の大きさを`diagonal_strains`と`off_diagonal_strains`によって指定します。`diagonal_strains`はセルの伸張と圧縮の大きさに相当します。`off_diagonal_strains`はセルの剪断変形に相当します。

* 必要に応じて，セル変形後の原子位置の最適化に使用するオプティマイザーを指定することもできます．デフォルトでは、`FireASEOptFeature`が使用されます。上記の例では明示的に`FireASEOptFeature`を指定しています。収束判定のようなオプティマイザー内のパラメータも必要に応じて変更可能です。

* 弾性テンソルが得られたら、`PostElasticPropertiesFeature`を使って、ヤング率などの種々の弾性特性を計算できます。


```python
opt = FireASEOptFeature(filter=True, show_progress_bar=True, tqdm_options={"desc":"cell optimization"}, estimator_fn=estimator_fn)

diagonal_strains = [-0.015, -0.01, -0.005, 0.005, 0.01, 0.015]
off_diagonal_strains = [-0.03, -0.02, -0.01, 0.01, 0.02, 0.03]
elastic_tensor_feature = ElasticTensorFeature(
    diagonal_strains=diagonal_strains,
    off_diagonal_strains=off_diagonal_strains,
    optimizer=FireASEOptFeature(estimator_fn=estimator_fn),
    show_logger=True,
    show_progress_bar=True,
    tqdm_options={"desc":"elastic tensor"},
    check_symmetry=True,
    estimator_fn=estimator_fn,
)

post_elastic = PostElasticPropertiesFeature()
```


* 以下はSiの計算例です。このセルの実行には数十秒程度の時間がかかります。


```python
opt_results = opt(atoms_Si)
elastic_tensor = elastic_tensor_feature(opt_results.atoms)
results_Si = post_elastic(elastic_tensor)
```


    cell optimization:   0%|          | 0/201 [00:00&lt;?, ?it/s]


    elastic tensor:   0%|          | 0/36 [00:00&lt;?, ?it/s]


    Calculating stress tensor
    Deform [0.984886 1.000000 1.000000 0.000000 0.000000 0.000000]
    Stress [-1.42e-02 -5.61e-03 -5.61e-03 5.77e-09 -7.20e-10 1.04e-08]
    Calculating stress tensor
    Deform [0.989949 1.000000 1.000000 0.000000 0.000000 0.000000]
    Stress [-9.34e-03 -3.54e-03 -3.54e-03 -9.98e-09 1.07e-09 1.07e-08]
    Calculating stress tensor
    Deform [0.994987 1.000000 1.000000 0.000000 0.000000 0.000000]
    Stress [-4.48e-03 -1.56e-03 -1.56e-03 7.12e-09 -2.43e-09 -2.68e-11]
    Calculating stress tensor
    Deform [1.004988 1.000000 1.000000 0.000000 0.000000 0.000000]
    Stress [5.47e-03 2.41e-03 2.41e-03 -3.30e-09 3.59e-09 4.27e-09]
    Calculating stress tensor
    Deform [1.009950 1.000000 1.000000 0.000000 0.000000 0.000000]
    Stress [1.03e-02 4.35e-03 4.35e-03 -1.62e-09 -7.47e-09 -2.97e-09]
    Calculating stress tensor
    Deform [1.014889 1.000000 1.000000 0.000000 0.000000 0.000000]
    Stress [1.50e-02 6.18e-03 6.18e-03 4.07e-09 -9.79e-09 5.02e-09]
    Calculating stress tensor
    Deform [1.000000 0.984886 1.000000 0.000000 0.000000 0.000000]
    Stress [-5.61e-03 -1.42e-02 -5.61e-03 1.50e-09 -3.21e-09 4.37e-09]
    Calculating stress tensor
    Deform [1.000000 0.989949 1.000000 0.000000 0.000000 0.000000]
    Stress [-3.54e-03 -9.34e-03 -3.54e-03 -6.93e-09 2.31e-09 2.68e-09]
    Calculating stress tensor
    Deform [1.000000 0.994987 1.000000 0.000000 0.000000 0.000000]
    Stress [-1.56e-03 -4.48e-03 -1.56e-03 3.34e-11 -1.73e-09 -6.19e-09]
    Calculating stress tensor
    Deform [1.000000 1.004988 1.000000 0.000000 0.000000 0.000000]
    Stress [2.41e-03 5.47e-03 2.41e-03 -3.29e-10 6.95e-10 -1.47e-09]
    Calculating stress tensor
    Deform [1.000000 1.009950 1.000000 0.000000 0.000000 0.000000]
    Stress [4.35e-03 1.03e-02 4.35e-03 -7.07e-09 -4.64e-09 -5.20e-10]
    Calculating stress tensor
    Deform [1.000000 1.014889 1.000000 0.000000 0.000000 0.000000]
    Stress [6.18e-03 1.50e-02 6.18e-03 -8.19e-09 -5.91e-09 2.39e-09]
    Calculating stress tensor
    Deform [1.000000 1.000000 0.984886 0.000000 0.000000 0.000000]
    Stress [-5.61e-03 -5.61e-03 -1.42e-02 -1.68e-09 -8.79e-09 -1.03e-08]
    Calculating stress tensor
    Deform [1.000000 1.000000 0.989949 0.000000 0.000000 0.000000]
    Stress [-3.54e-03 -3.54e-03 -9.34e-03 -1.09e-08 -1.60e-09 5.88e-09]
    Calculating stress tensor
    Deform [1.000000 1.000000 0.994987 0.000000 0.000000 0.000000]
    Stress [-1.56e-03 -1.56e-03 -4.48e-03 -7.17e-09 -1.05e-08 1.48e-09]
    Calculating stress tensor
    Deform [1.000000 1.000000 1.004988 0.000000 0.000000 0.000000]
    Stress [2.41e-03 2.41e-03 5.47e-03 -2.55e-09 -6.16e-09 8.42e-09]
    Calculating stress tensor
    Deform [1.000000 1.000000 1.009950 0.000000 0.000000 0.000000]
    Stress [4.35e-03 4.35e-03 1.03e-02 -1.30e-08 -1.27e-09 -8.50e-09]
    Calculating stress tensor
    Deform [1.000000 1.000000 1.014889 0.000000 0.000000 0.000000]
    Stress [6.18e-03 6.18e-03 1.50e-02 -4.99e-09 -5.00e-09 -9.50e-10]
    Calculating stress tensor
    Deform [1.000000 0.999549 0.999549 -0.030014 0.000000 0.000000]
    Stress [1.48e-03 -9.33e-04 -9.33e-04 -2.65e-02 4.10e-09 -9.51e-09]
    Calculating stress tensor
    Deform [1.000000 0.999800 0.999800 -0.020004 0.000000 0.000000]
    Stress [1.13e-03 -2.19e-05 -2.19e-05 -1.78e-02 -1.73e-08 -3.53e-09]
    Calculating stress tensor
    Deform [1.000000 0.999950 0.999950 -0.010001 0.000000 0.000000]
    Stress [7.98e-04 5.18e-04 5.18e-04 -8.56e-03 -1.12e-10 -3.94e-10]
    Calculating stress tensor
    Deform [1.000000 0.999950 0.999950 0.010001 0.000000 0.000000]
    Stress [7.98e-04 5.18e-04 5.18e-04 8.56e-03 -2.63e-09 -4.81e-09]
    Calculating stress tensor
    Deform [1.000000 0.999800 0.999800 0.020004 0.000000 0.000000]
    Stress [1.13e-03 -2.20e-05 -2.20e-05 1.78e-02 -2.61e-09 -1.50e-09]
    Calculating stress tensor
    Deform [1.000000 0.999549 0.999549 0.030014 0.000000 0.000000]
    Stress [1.48e-03 -9.33e-04 -9.33e-04 2.65e-02 9.49e-11 -4.00e-09]
    Calculating stress tensor
    Deform [0.999549 1.000000 0.999549 0.000000 -0.030014 0.000000]
    Stress [-9.33e-04 1.48e-03 -9.33e-04 -4.30e-09 -2.65e-02 -1.83e-09]
    Calculating stress tensor
    Deform [0.999800 1.000000 0.999800 0.000000 -0.020004 0.000000]
    Stress [-2.19e-05 1.13e-03 -2.19e-05 6.77e-09 -1.78e-02 1.55e-09]
    Calculating stress tensor
    Deform [0.999950 1.000000 0.999950 0.000000 -0.010001 0.000000]
    Stress [5.18e-04 7.98e-04 5.18e-04 -1.72e-09 -8.56e-03 -7.00e-09]
    Calculating stress tensor
    Deform [0.999950 1.000000 0.999950 0.000000 0.010001 0.000000]
    Stress [5.18e-04 7.98e-04 5.18e-04 -3.09e-09 8.56e-03 7.11e-09]
    Calculating stress tensor
    Deform [0.999800 1.000000 0.999800 0.000000 0.020004 0.000000]
    Stress [-2.20e-05 1.13e-03 -2.20e-05 -7.03e-09 1.78e-02 -1.50e-09]
    Calculating stress tensor
    Deform [0.999549 1.000000 0.999549 0.000000 0.030014 0.000000]
    Stress [-9.33e-04 1.48e-03 -9.33e-04 -8.19e-09 2.65e-02 -8.92e-09]
    Calculating stress tensor
    Deform [0.999549 0.999549 1.000000 0.000000 0.000000 -0.030014]
    Stress [-9.33e-04 -9.33e-04 1.48e-03 -2.35e-09 -7.43e-10 -2.65e-02]
    Calculating stress tensor
    Deform [0.999800 0.999800 1.000000 0.000000 0.000000 -0.020004]
    Stress [-2.19e-05 -2.19e-05 1.13e-03 -1.25e-09 -3.75e-09 -1.78e-02]
    Calculating stress tensor
    Deform [0.999950 0.999950 1.000000 0.000000 0.000000 -0.010001]
    Stress [5.18e-04 5.18e-04 7.98e-04 2.12e-09 -3.80e-10 -8.56e-03]
    Calculating stress tensor
    Deform [0.999950 0.999950 1.000000 0.000000 0.000000 0.010001]
    Stress [5.18e-04 5.18e-04 7.98e-04 1.81e-09 7.42e-09 8.56e-03]
    Calculating stress tensor
    Deform [0.999800 0.999800 1.000000 0.000000 0.000000 0.020004]
    Stress [-2.19e-05 -2.19e-05 1.13e-03 1.17e-08 -1.05e-08 1.78e-02]
    Calculating stress tensor
    Deform [0.999549 0.999549 1.000000 0.000000 0.000000 0.030014]
    Stress [-9.33e-04 -9.33e-04 1.48e-03 2.62e-09 -5.90e-10 2.65e-02]


    Deform [0.984886 1.000000 1.000000 0.000000 0.000000 0.000000]


    Stress [-1.42e-02 -5.61e-03 -5.61e-03 2.08e-09 7.70e-11 2.74e-09]


    Calculating stress tensor


    Deform [0.989949 1.000000 1.000000 0.000000 0.000000 0.000000]


    Stress [-9.34e-03 -3.54e-03 -3.54e-03 -6.85e-09 3.38e-09 1.75e-09]


    Calculating stress tensor


    Deform [0.994987 1.000000 1.000000 0.000000 0.000000 0.000000]


    Stress [-4.48e-03 -1.56e-03 -1.56e-03 -2.53e-10 -7.53e-09 -3.58e-09]


    Calculating stress tensor


    Deform [1.004988 1.000000 1.000000 0.000000 0.000000 0.000000]


    Stress [5.47e-03 2.41e-03 2.41e-03 2.37e-09 -7.25e-09 -1.36e-09]


    Calculating stress tensor


    Deform [1.009950 1.000000 1.000000 0.000000 0.000000 0.000000]


    Stress [1.03e-02 4.35e-03 4.35e-03 1.14e-09 -8.17e-09 1.31e-08]


    Calculating stress tensor


    Deform [1.014889 1.000000 1.000000 0.000000 0.000000 0.000000]


    Stress [1.50e-02 6.18e-03 6.18e-03 5.45e-09 -7.62e-09 3.69e-09]


    Calculating stress tensor


    Deform [1.000000 0.984886 1.000000 0.000000 0.000000 0.000000]


    Stress [-5.61e-03 -1.42e-02 -5.61e-03 1.27e-09 4.26e-09 -4.24e-10]


    Calculating stress tensor


    Deform [1.000000 0.989949 1.000000 0.000000 0.000000 0.000000]


    Stress [-3.54e-03 -9.34e-03 -3.54e-03 -4.09e-09 4.22e-09 4.08e-09]


    Calculating stress tensor


    Deform [1.000000 0.994987 1.000000 0.000000 0.000000 0.000000]


    Stress [-1.56e-03 -4.48e-03 -1.56e-03 -4.71e-09 -7.82e-09 8.24e-10]


    Calculating stress tensor


    Deform [1.000000 1.004988 1.000000 0.000000 0.000000 0.000000]


    Stress [2.41e-03 5.47e-03 2.41e-03 4.52e-10 -1.52e-09 -1.00e-09]


    Calculating stress tensor


    Deform [1.000000 1.009950 1.000000 0.000000 0.000000 0.000000]


    Stress [4.35e-03 1.03e-02 4.35e-03 -1.07e-08 -9.26e-09 -2.00e-09]


    Calculating stress tensor


    Deform [1.000000 1.014889 1.000000 0.000000 0.000000 0.000000]


    Stress [6.18e-03 1.50e-02 6.18e-03 -1.56e-09 -7.80e-09 9.41e-09]


    Calculating stress tensor


    Deform [1.000000 1.000000 0.984886 0.000000 0.000000 0.000000]


    Stress [-5.61e-03 -5.61e-03 -1.42e-02 -8.27e-09 -4.10e-09 5.28e-10]


    Calculating stress tensor


    Deform [1.000000 1.000000 0.989949 0.000000 0.000000 0.000000]


    Stress [-3.54e-03 -3.54e-03 -9.34e-03 -8.28e-09 1.12e-08 -8.94e-09]


    Calculating stress tensor


    Deform [1.000000 1.000000 0.994987 0.000000 0.000000 0.000000]


    Stress [-1.56e-03 -1.56e-03 -4.48e-03 -3.72e-09 -1.85e-08 -1.04e-09]


    Calculating stress tensor


    Deform [1.000000 1.000000 1.004988 0.000000 0.000000 0.000000]


    Stress [2.41e-03 2.41e-03 5.47e-03 -1.70e-09 -9.64e-10 1.26e-08]




    Calculating stress tensor


    Deform [1.000000 1.000000 1.009950 0.000000 0.000000 0.000000]


    Stress [4.35e-03 4.35e-03 1.03e-02 -7.51e-10 -9.59e-09 -2.13e-09]


    Calculating stress tensor


    Deform [1.000000 1.000000 1.014889 0.000000 0.000000 0.000000]


    Stress [6.18e-03 6.18e-03 1.50e-02 -7.90e-10 -1.02e-08 9.91e-09]


    Calculating stress tensor


    Deform [1.000000 0.999549 0.999549 -0.030014 0.000000 0.000000]


    Stress [1.48e-03 -9.32e-04 -9.32e-04 -2.65e-02 -7.49e-09 -6.97e-09]


    Calculating stress tensor


    Deform [1.000000 0.999800 0.999800 -0.020004 0.000000 0.000000]


    Stress [1.13e-03 -2.19e-05 -2.19e-05 -1.78e-02 5.53e-09 4.10e-09]


    Calculating stress tensor


    Deform [1.000000 0.999950 0.999950 -0.010001 0.000000 0.000000]


    Stress [7.98e-04 5.18e-04 5.18e-04 -8.56e-03 4.96e-09 3.44e-09]


    Calculating stress tensor


    Deform [1.000000 0.999950 0.999950 0.010001 0.000000 0.000000]


    Stress [7.98e-04 5.18e-04 5.18e-04 8.56e-03 1.01e-08 -1.39e-08]


    Calculating stress tensor


    Deform [1.000000 0.999800 0.999800 0.020004 0.000000 0.000000]


    Stress [1.13e-03 -2.19e-05 -2.19e-05 1.78e-02 1.32e-08 -1.45e-08]


    Calculating stress tensor


    Deform [1.000000 0.999549 0.999549 0.030014 0.000000 0.000000]


    Stress [1.48e-03 -9.32e-04 -9.32e-04 2.65e-02 2.93e-10 2.91e-09]


    Calculating stress tensor


    Deform [0.999549 1.000000 0.999549 0.000000 -0.030014 0.000000]


    Stress [-9.33e-04 1.48e-03 -9.33e-04 7.12e-09 -2.65e-02 -6.09e-09]


    Calculating stress tensor


    Deform [0.999800 1.000000 0.999800 0.000000 -0.020004 0.000000]


    Stress [-2.19e-05 1.13e-03 -2.19e-05 -2.59e-09 -1.78e-02 -8.01e-09]


    Calculating stress tensor


    Deform [0.999950 1.000000 0.999950 0.000000 -0.010001 0.000000]


    Stress [5.18e-04 7.98e-04 5.18e-04 -1.46e-09 -8.56e-03 -1.20e-09]


    Calculating stress tensor


    Deform [0.999950 1.000000 0.999950 0.000000 0.010001 0.000000]


    Stress [5.18e-04 7.98e-04 5.18e-04 -3.69e-09 8.56e-03 7.91e-09]


    Calculating stress tensor


    Deform [0.999800 1.000000 0.999800 0.000000 0.020004 0.000000]


    Stress [-2.19e-05 1.13e-03 -2.19e-05 4.94e-09 1.78e-02 -4.27e-09]


    Calculating stress tensor


    Deform [0.999549 1.000000 0.999549 0.000000 0.030014 0.000000]


    Stress [-9.33e-04 1.48e-03 -9.33e-04 7.01e-10 2.65e-02 1.32e-08]


    Calculating stress tensor


    Deform [0.999549 0.999549 1.000000 0.000000 0.000000 -0.030014]


    Stress [-9.33e-04 -9.33e-04 1.48e-03 -3.30e-09 -1.41e-09 -2.65e-02]


    Calculating stress tensor


    Deform [0.999800 0.999800 1.000000 0.000000 0.000000 -0.020004]


    Stress [-2.19e-05 -2.19e-05 1.13e-03 -3.57e-10 -9.79e-09 -1.78e-02]


    Calculating stress tensor


    Deform [0.999950 0.999950 1.000000 0.000000 0.000000 -0.010001]


    Stress [5.18e-04 5.18e-04 7.98e-04 4.85e-09 8.36e-10 -8.56e-03]


    Calculating stress tensor


    Deform [0.999950 0.999950 1.000000 0.000000 0.000000 0.010001]


    Stress [5.18e-04 5.18e-04 7.98e-04 6.34e-09 3.24e-09 8.56e-03]


    Calculating stress tensor


    Deform [0.999800 0.999800 1.000000 0.000000 0.000000 0.020004]


    Stress [-2.19e-05 -2.19e-05 1.13e-03 5.54e-09 -7.18e-09 1.78e-02]


    Calculating stress tensor


    Deform [0.999549 0.999549 1.000000 0.000000 0.000000 0.030014]


    Stress [-9.33e-04 -9.33e-04 1.48e-03 5.92e-09 9.19e-10 2.65e-02]


* 計算結果を表示してみます。


```python
print("Elastic tensor")
print(np.round(elastic_tensor.elastic_tensor))
print(f"Bulk Modulus {results_Si.bulk_modulus:3.2f} GPa")
print(f"Young's Modulus {results_Si.young_modulus:3.2f} GPa")
print(f"Shear Modulus {results_Si.shear_modulus:3.2f} GPa")
print(f"Poisson ratio {results_Si.poisson_ratio:3.2f}")

prop = results_Si.get_prop_results(dir1=[1,1,1], dir2=[1,-1,0])
print(f"Young's Modulus along [1 1 1]", prop["young's modulus along [1 1 1]"]["GPa"])
print(f"Shear modulus along [1 1 1] [1 -1 0]", prop["shear modulus along [1 1 1] [1 -1 0]"]["GPa"])
```

    Elastic tensor
    [[157.  63.  63.   0.  -0.  -0.]
     [ 63. 157.  63.  -0.  -0.  -0.]
     [ 63.  63. 157.   0.  -0.   0.]
     [  0.  -0.   0.  71.  -0.   0.]
     [ -0.  -0.  -0.  -0.  71.  -0.]
     [ -0.  -0.   0.   0.  -0.  71.]]
    Bulk Modulus 94.36 GPa
    Young's Modulus 148.49 GPa
    Shear Modulus 59.99 GPa
    Poisson ratio 0.24
    Young's Modulus along [1 1 1] 169.71301717289177
    Shear modulus along [1 1 1] [1 -1 0] 52.84381007603689


## `ElasticFeature`を利用した計算

* `ElasticFeature`は、ここまでに紹介した`ASEOptFeature`、`ElasticTensorFeature`、`PostElasticPropertiesFeature`を組み合わせて一度に計算してくれるfeatureです。今度はこれを使ってAl2O3の計算をしてみましょう。

* なお、計算誤差などの都合で、弾性テンソルの中でしばしば0になるべき値に微少量が入ることがあります。このときにWarningが表示されることがあります。


```python
elastic_feature = ComplexElasticFeature(
    FireASEOptFeature(
        filter=True,
        show_progress_bar=True,
        tqdm_options={"desc": "cell optimization"},
        estimator_fn=estimator_fn,
    ),
    ElasticTensorFeature(
        diagonal_strains,
        off_diagonal_strains,
        FireASEOptFeature(estimator_fn=estimator_fn),
        show_progress_bar=True,
        tqdm_options={"desc": "elastic tensor"},
        check_symmetry=True,
        estimator_fn=estimator_fn,
    ),
    PostElasticPropertiesFeature(),

)
results_Al2O3 = elastic_feature(atoms_Al2O3)
```


    cell optimization:   0%|          | 0/201 [00:00&lt;?, ?it/s]


    elastic tensor:   0%|          | 0/36 [00:00&lt;?, ?it/s]


## 実験結果との比較

* ここまでに計算した弾性率について、実験データと比較してみましょう。


```python
Si_expt = {
    "young's_modulus": 160,
    "bulk_modulus": 100,
    "shear_modulus": 62,
    "poisson_ratio": 0.27,
}
Al2O3_expt = {
    "young's_modulus": 375,
    "bulk_modulus": 228,
    "shear_modulus": 152,
    "poisson_ratio": 0.22,
}

```


```python
from plotly.offline import iplot
import plotly.graph_objs as go
from plotly.subplots import make_subplots

def plot_comparison_fig(matlantis, expt, material=""):
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    keys = ["young's_modulus", "bulk_modulus", "shear_modulus"]
    values = [matlantis[k]["GPa"] for k in keys]
    fig.add_trace(
        go.Bar(
            x=keys, y=values,
            textposition="auto", name="matlantis", marker={"color": "#636EFA"},
        )
    )
    values = [expt[k] for k in keys]
    fig.add_trace(
        go.Bar(
            x=keys, y=values,
            textposition="auto", name="experiment", marker={"color": "#EF553B"},
        )
    )
    fig.add_trace(
        go.Bar(
            x=["poisson_ratio"], y=[matlantis["poisson_ratio"]],
            textposition="auto", showlegend=False, marker={"color": "#636EFA"},
        ),
        secondary_y=True,
    )
    fig.add_trace(
        go.Bar(
            x=["poisson_ratio"], y=[expt["poisson_ratio"]],
            textposition="auto", showlegend=False, marker={"color": "#EF553B"},
        ),
        secondary_y=True,
    )
    fig.update_layout(
        title={"text": material, "x": 0.2},
        barmode="group",
        legend={"xanchor": "right", "yanchor": "top", "y": 1.35},
        width=500,
        height=400,
        margin=go.layout.Margin(l=50, r=25, b=50, t=50, pad=4),
    )
    fig.update_yaxes(title_text="Elasticity [GPa]")
    fig.update_yaxes(title_text="possion ratio", secondary_y=True)
    return fig

```


* Siの結果を以下に示します。

* 初回起動時などで図が表示されない場合、ブラウザでページをリロードすると表示されることがあります。


```python
fig = plot_comparison_fig(results_Si.get_prop_results(), Si_expt, "Si")
iplot(fig)
```


* Al2O3の結果を以下に示します。


```python
fig = plot_comparison_fig(results_Al2O3.elastic_properties.get_prop_results(), Al2O3_expt, "Al2O3")
iplot(fig)
```


