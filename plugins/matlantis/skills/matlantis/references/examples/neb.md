# Matlantis-features: NEB


## 概要
* NEB(Nudged Elastic Band)法は、初期構造と最終構造を離散的に補完した構造を多数用意し、それら複数の構造を位相空間上で互いに結合させつつ構造最適化を行うことで反応経路を探す手法です。

* Matlantis-featuresでは、NEB法を実行するするフィーチャーとして、`NEBFeature`を提供しています。


## 初期設定
* 必要なライブラリ等の準備を行います。


```python
# If you have already installed matlantis-features, you can skip this.
%pip install 'matlantis-features&gt;=0.18.0'
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
from plotly.offline import iplot

import ase.io
from ase.visualize import view

from matlantis_features.utils.calculators import pfp_estimator_fn
from matlantis_features.filters import UnitCellASEFilter
from matlantis_features.features.common import (
    BFGSASEOptFeature, FireASEOptFeature, NEBFeature
)
# For save and reload, please use matlantis_features.symmetry.FixSymmetry,
# rather than ase.spacegroup.symmetrize.FixSymmetry.
from matlantis_features.symmetry import FixSymmetry

try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

```


## 材料の準備
* ここではナトリウム電池材料 NaFePO$_4$中のNa拡散に対する反応経路（拡散パス）を扱います。
* まずは、NaFePO$_4$の構造をcifファイルから読み込み、構造緩和を行います。なお、cifファイルは NaFePO4 (mp-746030) を [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) の下で Materials Project から入手したものです。
* ここでは、元の構造の対称性を保持する`FixSymmetry`を拘束条件として設定し、セルパラメータを含めた最適化を行うため、`UnitCellASEFilter`を用います。


```python
bulk_ = ase.io.read( dir_path / "assets/neb/NaFePO4.cif" )
bulk_.set_constraint([FixSymmetry(bulk_, adjust_cell=True, adjust_positions=True)])
ucf = UnitCellASEFilter()
```


* BFGS法による構造緩和を行う`BFGSASEOptFeature`を用います。`UnitCellASEFilter`の効果を得るために、`ucf`を`filter`パラメータとして渡すことに注意してください。
* 対称性による拘束や、セルパラメータを含めた最適化は以降行わないため、`ase.Atoms`オブジェクトである`bulk`を取り出したのち、拘束条件をリセットします。


## `estimator_fn` による計算モードとモデルバージョンの指定
* 計算対象には遷移金属（Fe）が含まれているので、デフォルトの `PBE` ではなく、Hubbard U 補正が含まれている `PBE_U` モードを指定する必要があります。
* Feature に `estimator_fn` 引数を使うことで、matlantis-features で使われる PFP の計算モードとモデルバージョンを指定することができます。
* `estimator_fn` は、Estimator オブジェクトを作る factory method です。計算モードとモデルバージョンのみを指定する場合は、以下のように `pfp_estimator_fn` を利用できます。より細かい指定が必要な場合は、自身で factory method を定義できます。詳細に関しては [Matlantis Guidebook](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html#estimator-fnpfp) をご参照ください。
* 構造緩和で使う `BFGSASEOptFeature` と、 `NEBFeature` で使う `FireASEOptFeature` に `estimator_fn` を渡します。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode

pfp_crystal_estimator_fn = pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE_U)
```


```python
opt = BFGSASEOptFeature(
    fmax=0.01,
    show_logger=True,
    logger_interval=10,
    show_progress_bar=True,
    filter=ucf,
    estimator_fn=pfp_crystal_estimator_fn,
)
bulk = opt(bulk_).atoms.ase_atoms
bulk.set_constraint([])

display(view(bulk, viewer="ngl"))
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


    steps:     0  energy：-4.47 eV/atom  max_force:  0.17 eV/Ang
    steps:    10  energy：-4.47 eV/atom  max_force:  0.04 eV/Ang
    steps:    20  energy：-4.47 eV/atom  max_force:  0.02 eV/Ang
    steps:    30  energy：-4.47 eV/atom  max_force:  0.01 eV/Ang





    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Na', 'Fe', 'P', 'O'),…


## 初期構造と最終構造の準備
* 初期構造と最終構造を準備します。まずは$2\times2\times2$のスーパーセルを用意し、0番目のNaが3番目のNaの位置に移動する構造を生成します。続いて、`BFGSASEOptFeature`により構造緩和を行います。
* `BFGSASEOptFeature` に `estimator_fn` 引数を渡すことで、計算モードを `PBE_U` モードに設定します。


```python
# initial: remove index 'i' atom
# final: swap positions of index 'i' atom and 'j' atom, then remove index 'i' atoms
def get_initial_and_final(atoms, i, j):
    initial = atoms.copy()
    final = atoms.copy()
    final.positions[[i,j]] = final.positions[[j,i]]
    initial.pop(i)
    final.pop(i)
    return initial, final
```


```python
supercell = bulk * (2,2,2)
initial_, final_ = get_initial_and_final(supercell, 3, 0)

opt_initial = BFGSASEOptFeature(
    fmax=0.01,
    show_logger=True,
    logger_interval=10,
    show_progress_bar=True,
    estimator_fn=pfp_crystal_estimator_fn,
)
initial = opt_initial(initial_).atoms

opt_final = BFGSASEOptFeature(
    fmax=0.01,
    show_logger=True,
    logger_interval=10,
    show_progress_bar=True,
    estimator_fn=pfp_crystal_estimator_fn,
)
final = opt_final(final_).atoms
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


    steps:     0  energy：-4.95 eV/atom  max_force:  1.74 eV/Ang
    steps:    10  energy：-4.95 eV/atom  max_force:  0.22 eV/Ang
    steps:    20  energy：-4.95 eV/atom  max_force:  0.04 eV/Ang


      0%|          | 0/201 [00:00&lt;?, ?it/s]


    steps:     0  energy：-4.95 eV/atom  max_force:  1.74 eV/Ang
    steps:    10  energy：-4.95 eV/atom  max_force:  0.22 eV/Ang
    steps:    20  energy：-4.95 eV/atom  max_force:  0.04 eV/Ang


    steps:    10  energy：-4.95 eV/atom  max_force:  0.22 eV/Ang


    steps:    20  energy：-4.95 eV/atom  max_force:  0.04 eV/Ang


## 初期構造と最終構造の確認
* 準備した、初期構造と最終構造を見てみます。Naの1つが移動していることがわかります。


```python
v = view([initial.ase_atoms, final.ase_atoms], viewer="ngl")
display(v)
```


    HBox(children=(NGLWidget(max_frame=1), VBox(children=(Dropdown(description='Show', options=('All', 'Na', 'Fe',…


## NEBの実行
* `NEBFeature`を用いてNEBを実行します。`FireASEOptFeature`をオプティマイザとして指定します。`FireASEOptFeature` に `estimator_fn` 引数を渡すことで、 計算モードを `PBE_U` モードに設定します。初期構造と最終構造の中間構造の数は`n_images=5`で、5つに設定します。また、遷移状態付近の構造をえるために`climb=True`を指定し、NEB法を改良したCI-NEB法をここでは実行します。
*`idpp=True`とすると、IDPP (Image Dependent Pair Potential)法により中間構造の初期推定を行うことができます。この手法は分子などにおいては有効な一方で、結晶などでは破綻してしまうケースが散見されます。このため、`idpp=False`を指定しています。
* この計算には少なくとも数分以上かかります。


```python
neb = NEBFeature(
    optimizer = FireASEOptFeature(
        fmax=0.05,
        n_run=200,
        show_logger=True,
        logger_interval=10,
        show_progress_bar=True,
        estimator_fn=pfp_crystal_estimator_fn,
    ),
    n_images = 5,
    climb = True, # use CI-NEB to obtain activation energy
    idpp = False, # idpp interpolate may create broken initial NEB image for crystal system
)
neb_res = neb(initial, final)
neb_res.converged
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


    steps:     0  energies: -1104.05 -1103.73 -1102.99 -1102.55 -1102.99 -1103.73 -1104.05 eV  max_force:  3.22 eV/Ang
    steps:    10  energies: -1104.05 -1103.92 -1103.66 -1103.54 -1103.66 -1103.92 -1104.05 eV  max_force:  0.70 eV/Ang
    steps:    20  energies: -1104.05 -1103.94 -1103.75 -1103.68 -1103.75 -1103.94 -1104.05 eV  max_force:  0.42 eV/Ang
    steps:    30  energies: -1104.05 -1103.95 -1103.77 -1103.72 -1103.77 -1103.95 -1104.05 eV  max_force:  0.07 eV/Ang
    steps:    40  energies: -1104.05 -1103.95 -1103.78 -1103.73 -1103.78 -1103.95 -1104.05 eV  max_force:  0.06 eV/Ang


    True


## 結果の確認
* 得られた反応経路のエネルギーは以下のようにプロットすることができます。

* 初回起動時などで図が表示されない場合、ブラウザでページをリロードすると表示されることがあります。


```python
fig = neb_res.plot()
fig.update_layout(width=800, height=500)
iplot(fig)
```


* 活性化エネルギーを確認します。既存のDFT計算での値は 0.35 eV です。(doi:10.1021/acsomega.9b04213)


```python
energies = neb_res.energy_log[-1]
print("delta_Ea = {:.2f} eV".format(np.max(energies) - np.min(energies)))
```

    delta_Ea = 0.32 eV


* 最後に得られた反応経路の構造を確認しましょう。


```python
v = view([atoms.ase_atoms for atoms in neb_res.neb_images], viewer="ngl")
display(v)
```


    HBox(children=(NGLWidget(max_frame=6), VBox(children=(Dropdown(description='Show', options=('All', 'Na', 'Fe',…


#### JSONファイルへの結果とパラメータの保存と読み込み

シミュレーション結果と呼び出し引数はJSONファイルに保存して、後から読み込むことができます。
なお、いくつかの `Feature` の初期化パラメータは読み込みに対応していません。
具体的には、`estimator_fn` はサポートされていません。
これを `Feature` 初期化時に設定したとしても、 `saved_result.feature` に入っている `Feature` オブジェクトのパラメーターはデフォルト値となっています。


```python
# Save the optimization result.
save_success = neb_res.save("neb_result.json")
print(save_success)
```

    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/utils/save_util.py:62: UserWarning:

    The parameter Spacegroup does not support saving and reloading. It will be ignored upon saving and reloading.



    True


```python
# Reload the saved result. Note that UserWarning will be shown because estimator_fn is not supported for savaing/reloading.
from matlantis_features.features.common.neb import NEBFeatureResult

saved_result = NEBFeatureResult.load("neb_result.json")
```

    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/common/opt.py:154: UserWarning:

    Model version is not specified, default 'v8.0.0' will be used.Please specify it by providing an estimator_fn or setting the 'MATLANTIS_PFP_MODEL_VERSION' environment variable.

    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/common/opt.py:154: UserWarning:

    Calculation mode is not specified, default 'EstimatorCalcMode.PBE' will be used.Please specify it by providing an estimator_fn or setting the 'MATLANTIS_PFP_CALC_MODE' environment variable.



```python
# Inspect loaded results

# PFP parameters (model version, calc mode, etc.)
print("Calculator info:", saved_result.calculator_info)

# Package version information
print("Version info:", saved_result.version_info)

# Inspect the feature init parameters (as dict)
print("Feature init params:", saved_result.feature.to_dict())
```

    Calculator info: {'model_version': 'v8.0.0', 'calc_mode': 'pbe_u', 'method_type': 'pfvm_d3_pfvm'}
    Version info: {'matlantis-feature-result-version': '0.0.1', 'ase': '3.26.0', 'pfp-api-client': '2.0.1', 'matlantis-features': '1.0.1', 'python': '3.13'}
    Feature init params: {'cls_path': 'matlantis_features.features.common.neb.NEBFeature', 'optimizer': {'cls_path': 'matlantis_features.features.common.opt.FireASEOptFeature', 'n_run': 200, 'fmax': 0.05, 'show_progress_bar': True, 'tqdm_options': None, 'show_logger': True, 'logger_interval': 10, 'filter': None, 'trajectory': None, 'maxstep': 0.2}, 'n_images': 5, 'k': 0.1, 'climb': True, 'method': 'aseneb', 'idpp': False}


```python
# Rerunning the simulation
# result_rerun = saved_result.feature(**saved_result.call_params)
```
