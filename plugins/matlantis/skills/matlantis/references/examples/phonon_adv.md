# Matlantis-features: フォノン Advanced usage
**Contents**

* [0. 初期設定](#advanced_0)
* [1. 力定数マトリクスの計算で格子対称性を利用する方法は？](#advanced_1)
* [2. 入力構造をprimitive cellに変換する方法は？](#advanced_2)
* [3. バンド構造計算でk点パスを設定する方法は？](#advanced_3)
* [4. Partial DOSを取得する方法は？](#advanced_4)
* [5. フォノンモードを視覚化する方法は？](#advanced_5)


&lt;a id='advanced_0'&gt;&lt;/a&gt;
## 0. 初期設定

* 必要なライブラリ等の準備を行います。


```python
# If you have already installed matlantis-features, you can skip this.
# Version 0.11.0 or later is required to run this notebook.
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


&lt;a id='advanced_1'&gt;&lt;/a&gt;
## 1. 力定数マトリクスの計算で格子対称性を利用する方法は？
* 結晶材料は通常、多くの対称操作を持っています。 格子対称性を利用することで、フォノン計算時の力推定の回数を大幅に減らすことができます。
* Matlantis-featuresでは、ForceConstantFeatureで `use_symmetry = True`を設定するだけで、格子対称性を使用できます。


### 1-1 材料の準備
* 例としてCrystal Siを使用します。


```python
import numpy as np
from ase.build import bulk

from matlantis_features.features.common import FireASEOptFeature
from matlantis_features.features.phonon import ForceConstantFeature

atoms_Si = bulk("Si")
opt = FireASEOptFeature(fmax=0.0001, filter=True, show_progress_bar=True, estimator_fn=estimator_fn)
atoms_opt = opt(atoms_Si).atoms
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


### 1-2対称性が　ある/ない　場合の力定数マトリクス計算


```python
fc_wo_symmetry = ForceConstantFeature(
    supercell = (4,4,4),
    delta = 0.1,
    show_logger=False,
    show_progress_bar=True,
    estimator_fn=estimator_fn,
)
fc_w_symmetry = ForceConstantFeature(
    supercell = (4,4,4),
    delta = 0.1,
    show_logger=False, show_progress_bar=True,
    use_symmetry=True,
    estimator_fn=estimator_fn,
)
```


* 結晶Si（primitive cell）の力定数マトリクスは、格子対称性を考慮しない場合、13回の力推定が必要です。
* これに対して、格子対称性を考慮した場合、力定数の計算に必要な力の推定は1回のみです。 力推定の回数が大幅に削減されます。


```python
force_constant_wo_symmetry = fc_wo_symmetry(atoms_opt)
force_constant_w_symmetry  = fc_w_symmetry(atoms_opt)
```


      0%|          | 0/13 [00:00&lt;?, ?it/s]


    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/phonon/symmetry.py:47: DeprecationWarning: dict interface is deprecated. Use attribute interface instead
      self.num_symmetry_operations = len(self.database["rotations"])


      0%|          | 0/1 [00:00&lt;?, ?it/s]


* 格子対称性を利用しても計算結果は変わりません。


```python
print(np.round(force_constant_wo_symmetry.force_constant[0], 2))
print(np.round(force_constant_w_symmetry.force_constant[0], 2))
```

    [[13.41  0.   -0.   -3.28 -2.33 -2.33]
     [ 0.   13.41  0.   -2.33 -3.28 -2.33]
     [-0.    0.   13.41 -2.33 -2.33 -3.28]
     [-3.28 -2.33 -2.33 13.41  0.   -0.  ]
     [-2.33 -3.28 -2.33  0.   13.41  0.  ]
     [-2.33 -2.33 -3.28 -0.    0.   13.41]]
    [[13.41  0.   -0.   -3.28 -2.33 -2.33]
     [ 0.   13.41 -0.   -2.33 -3.28 -2.33]
     [-0.    0.   13.41 -2.33 -2.33 -3.28]
     [-3.28 -2.33 -2.33 13.41  0.   -0.  ]
     [-2.33 -3.28 -2.33  0.   13.41 -0.  ]
     [-2.33 -2.33 -3.28 -0.   -0.   13.41]]


&lt;a id='advanced_2'&gt;&lt;/a&gt;
## 2. 入力構造をprimitive cellに変換する方法は？

* 場合によっては、入力構造がprimitive cell（原子数が最小の単位）ではなく、conventional cell またはさらに大きな supercell であることがあります。 その場合、以下のような2つの問題が起きます。
    * Force Constant の計算に時間がかかる。
    * フォノンバンド構造により多くのバンドがあり、画像が乱雑になる。

* Matlantis-featuresは、supercellをprimitive cellに変換する方法を提供します。


### 2-1 材料の準備
* 入力構造は、32個の原子を含むCrystal Siの2x2x1 supercellです。
* （注：Crystal Siのprimitive cellには2つの原子しかありません）


```python
import numpy as np

from ase.build import bulk
from ase.visualize import view

from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.common import FireASEOptFeature
from matlantis_features.features.phonon import ForceConstantFeature
from matlantis_features.features.phonon.utils import get_primitive_structure

atoms_Si = bulk("Si", cubic=True) * (2,2,1)
view(atoms_Si, viewer="ngl")
```





    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), value='All'), D…


### 2-2 2x2x1 supercellを使用した力定数マトリクス計算
* 2x2x1 supercellでは、193回の力推定する必要があります。


```python
fc = ForceConstantFeature(
    supercell = (4,4,4),
    delta = 0.1,
    show_logger=False,
    show_progress_bar=True,
    estimator_fn=estimator_fn,
)

force_constant = fc(atoms_Si)
```


      0%|          | 0/193 [00:00&lt;?, ?it/s]


### 2-3 SupercellをPrimitive Cellに変換して力定数マトリクスを計算する
* ForceConstantFeatureの呼び出し時に `primitive_matrix`が指定されている場合、力定数マトリクスを計算する前に、入力構造がPrimitive Cellに変換されます。
* 分率座標 `primitive_matrix`の1行目、2行目、3行目はプリミティブセル** a **、** b **、** c **の基底ベクトル。
* `primitive_matrix`パラメータを使用してprimitive cellを明示的に定義したとき、ForceConstantFeatureの力推定の回数が193から13に減少します。


```python
force_constant = fc(
    atoms_Si,
    primitive_matrix=np.array(
        [
            [0.00, 0.25, 0.25],
            [0.25, 0.00, 0.25],
            [0.50, 0.50, 0.00],
        ]
    ),
)
```


      0%|          | 0/13 [00:00&lt;?, ?it/s]


* phonon計算やその他の目的のために `get_primitive_structure` 関数を直接使用して primitive cell を得ることも可能です。


```python
atoms_Si = bulk("Si", cubic=True) * (2,2,1)
atoms_Si_primitive = get_primitive_structure(
    atoms_Si,
    primitive_matrix=np.array(
        [
            [0.0, 0.25, 0.25],
            [0.25, 0.0, 0.25],
            [0.5, 0.5, 0.0],
        ]
    )
)
view(atoms_Si_primitive, viewer="ngl")
```


    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), value='All'), D…


* `primitive_matrix` が厳密でない場合には ValueErrorが起きることに注意が必要です。


```python
atoms_Si = bulk("Si", cubic=True) * (2,2,1)
try:
    atoms_Si_primitive = get_primitive_structure(
        atoms_Si,
        primitive_matrix=np.array(
            [
                [0.0, 0.1, 0.1],
                [0.1, 0.0, 0.1],
                [0.1, 0.1, 0.0],
            ]
        )
    )
except ValueError as e:
    print("Raise Exception ", e)

```

    Raise Exception


### 2-4 Primitive cell の自動検出
* primitive cell の基底ベクトルが不明な場合、primitive_matrix='auto'を設定することで、primitive cellを自動的に見つけさせることができます。
* この場合、入力構造の対称性を解析することでprimitive cellを生成します。
* 力推定の回数が13に減少したことから、自動検出法を用いて以前と 上記の実験と同じprimitive cellを見つけられたことが示唆されます。


```python
force_constant = fc(atoms_Si, primitive_matrix="auto")
```


      0%|          | 0/13 [00:00&lt;?, ?it/s]


&lt;a id='advanced_3'&gt;&lt;/a&gt;
## 3. バンド構造計算でk点パスを設定する方法は？
* フォノンバンド構造を取得するには、k点パスを定義する必要があります。 matlantis-featuresでは、k点パスを指定する方法がいくつかあります。


```python
import numpy as np
from ase.build import bulk

from matlantis_features.features.common import FireASEOptFeature
from matlantis_features.features.phonon import ForceConstantFeature, PostPhononBandFeature
from plotly.offline import iplot
```


### 3-1 材料と力定数マトリクスの準備
* 例としてCrystalSiを使用しています。


```python
atoms_Si = bulk("Si")
opt = FireASEOptFeature(fmax=0.0001, filter=True, show_progress_bar=True, estimator_fn=estimator_fn)
atoms_opt = opt(atoms_Si).atoms

fc = ForceConstantFeature(
    supercell = (10,10,10),
    delta = 0.1,
    show_logger=True,
    show_progress_bar=True,
    estimator_fn=estimator_fn,
)
force_constant = fc(atoms_opt)
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


      0%|          | 0/13 [00:00&lt;?, ?it/s]


### 3-2 k点パスの自動定義
* 最も簡単なのは、デフォルトのk点パスを使用することです。
* PostPhononBandFeatureを呼び出すときに `labels`および` special_kpts`パラメーターを指定しない場合、k点パスはBravais latticeに従って自動的に生成されます。


```python
band = PostPhononBandFeature()
band_results = band(force_constant)

fig = band_results.plot()
fig.update_layout(width=800, height=500)
iplot(fig)
```


## 3-3 labelsを使用してk点パスを定義
* Matlantis-featuresではいくつかのspecial k-pointsを接続することでkpointパスを自分で定義することもできます。
* special kpoint名のリストを `labels = ['G'、 'X'、 'W']` とした場合、kpointパスは、2つの隣接するspecial kpoint間を線形補間することによって構築されます。
* special kpointの名称と座標を以下に紹介します。
* 符号「|」 は切断のために使用されます。 たとえば、 `labels = ['G'、 'X'、 '|'、 'L'、 'U']`は、GからXへのパス、およびLからUへのパスを意味しますが、XとLは切断されています。


```python
band = PostPhononBandFeature()
band_results = band(force_constant, labels = ['G', 'X', 'W', 'K', 'G', 'L', '|', 'U', 'X'])
fig = band_results.plot()
fig.update_layout(width=650, height=500)
iplot(fig)
```


* special k-pointsの名前と座標はASEで定義され、次の方法で確認することができます。


```python
from ase.lattice import identify_lattice

lattice = identify_lattice(atoms_Si.cell)[0]
print("List of special k-points", lattice.special_point_names)
print("Coordination of X point", lattice.get_special_points()["X"])
```

    List of special k-points ['G', 'K', 'L', 'U', 'W', 'X']
    Coordination of X point [0.5 0.  0.5]


### 3-4 k点パスを自分で定義する
* Matlantis-featuresでは上記の2つの方法で要求を満たすことができない場合に、ユーザー自身でspecial k-pointsを定義することもできます。
* この場合、ユーザーはspecial k-points名のリストを `label`に渡し、対応するspecial k-pointsの座標を `special_kpts`に渡す必要があります。
* 注意：
     * `|` を除いて、各ラベルには対応する座標が必要です。
     * この場合、special k-pointsの名前はASEの定義によって異なる場合があります。


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
)
fig = band_results.plot()
fig.update_layout(width=800, height=500)
iplot(fig)
```


### 3-5 k点の数を設定する
* ユーザーはintのリストを `n_kpts`に渡して、各セグメントの補間数を設定できます。
* この場合、ユーザーは各セグメントにintを指定する必要があります。 たとえば、5つのspecial k-pointsがある場合、ユーザーは `n_kpts`に4つのintを指定する必要があります


```python
band = PostPhononBandFeature()
band_results = band(force_constant, labels = ['G', 'X', 'W', 'K', 'G', 'L', '|', 'U', 'X'], n_kpts=[50, 10, 20, 10, 10, 0, 30])
fig = band_results.plot()
for trace in fig.select_traces():
    trace.mode = "lines+markers"
fig.update_layout(width=650, height=500)
iplot(fig)
```


&lt;a id='advanced_4'&gt;&lt;/a&gt;
## 4. Partial DOSを取得する方法は？
* Partial DOSは、多元素材料の各元素からの寄与として決定されます。 これにより、材料研究者はさまざまな振動モードの性質を理解できます。


```python
import numpy as np
from ase import Atoms
from matlantis_features.features.common import FireASEOptFeature
from matlantis_features.features.phonon import ForceConstantFeature, PostPhononDOSFeature
from plotly.offline import iplot
```


### 4-1 材料と力定数マトリクスの準備
*　Partial DOSを計算するために、例として2成分酸化物（ZnO）を使用します。


```python
atoms_ZnO = Atoms(
        cell=np.array(
            [
                [3.289102, 0.000000, 0.000000],
                [-1.644551, 2.848446, 0.000000],
                [0.000000, 0.000000, 5.306821],
            ]
        ),
        symbols=["Zn", "Zn", "O", "O"],
        scaled_positions=np.array(
            [
                [0.333333, 0.666667, 0.000548],
                [0.666667, 0.333333, 0.500548],
                [0.333333, 0.666667, 0.379762],
                [0.666667, 0.333333, 0.879762],
            ]
        ),
        pbc=True,
    )
opt = FireASEOptFeature(fmax=0.0001, filter=True, show_progress_bar=True, estimator_fn=estimator_fn)
atoms_opt = opt(atoms_ZnO).atoms
fc = ForceConstantFeature(
    supercell = (4,4,4),
    delta = 0.1,
    show_logger=False, show_progress_bar=True,
    use_symmetry=True,
    estimator_fn=estimator_fn,
)
force_constant = fc(atoms_opt)
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/phonon/symmetry.py:47: DeprecationWarning:

    dict interface is deprecated. Use attribute interface instead



      0%|          | 0/4 [00:00&lt;?, ?it/s]


### 4-2 Partial DOS計算
* Partial DOSは、PostPhononDOSFeatureを初期化するときに `partial = True`を設定することで簡単に有効化できます。
* 図ではZnおよびO元素の寄与を示しています。
* Partial DOSは追加の計算を行うことに注意してください。必要に応じて有効にしてください。


```python
dos = PostPhononDOSFeature(partial=True)
dos_results = dos(force_constant, kpts=[5, 5, 5], scheme="gamma", unit="meV")
fig = dos_results.plot()
fig.update_layout(width=800, height=500)
iplot(fig)
```


&lt;a id='advanced_5'&gt;&lt;/a&gt;
## 5. フォノンモードを視覚化する方法は？
* PostPhononModeFeatureを使用して、特定のフォノンモードで原子がどのように移動するかを視覚化できます。


```python
from ase.build import bulk
from ase.io import Trajectory
from ase.visualize import view

from matlantis_features.features.common import FireASEOptFeature
from matlantis_features.features.phonon import ForceConstantFeature, PostPhononModeFeature
```


### 5-1 材料と力定数マトリクスの準備
* 例としてCrystalSiを使用しています。


```python
atoms_Si = bulk("Si")
opt = FireASEOptFeature(fmax=0.0001, filter=True, show_progress_bar=True, estimator_fn=estimator_fn)
atoms_opt = opt(atoms_Si).atoms
fc = ForceConstantFeature(
    supercell = (10,10,10),
    delta = 0.1,
    show_logger=False, show_progress_bar=True,
    estimator_fn=estimator_fn,
)
force_constant = fc(atoms_opt)
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]


      0%|          | 0/13 [00:00&lt;?, ?it/s]


### 5-2 フォノンモードを視覚化する
* PostPhononModeFeatureは、特定のフォノンモードでの実際の原子の動きを取得するために使用されます。 与えられたフォノンモードでの各原子の固有変位を抽出し、それを一連の原子構造に変換します。
* 見たいk-pointを`k_point`に定義します。 次の例では、Xポイント（0.5、0.0、0.5）でフォノンモードが表示されます。
* Crystal SiはX点で6つのフォノン周波数を持っているため、6つの振動軌道が得られ、これは `mode_results.trajectories`に保存されます。
* フォノン周波数の昇順で軌跡を保存します。 例えば mode_results.trajectories[0]は、Xポイントで最も低い周波数を持つフォノンモードです。


```python
mode = PostPhononModeFeature()
mode_results = mode(
    force_constant,
    k_point=(0.5, 0.0, 0.5), # X
    repeat_of_cell=(4, 4, 4),
    n_images=30,
    amplitude=1.0
)
```


```python
v=view(mode_results.trajectories[0], viewer="ngl")
v.view.add_representation("ball+stick")
display(v)
```


    HBox(children=(NGLWidget(max_frame=29), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), val…

