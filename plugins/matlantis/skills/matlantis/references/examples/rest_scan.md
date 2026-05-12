# Restraint Scan

Restscanは、分子などの構造を連続的に変化させ、異なる構造を得るための機能です。Gaussianに実装されているscanに類似しています。この機能を使用することで、単純な構造を逐次的に変化させ、より複雑な構造を得ることができます。また、反応経路のinitial guessを作成するために、initial stateから構造を逐次的に変化させ、反応経路とfinal stateを得ることもできます。

Restscanの特徴は、構造変化の実装に holonomic constraint を用いない点にあります。代わりに、分子の元々の力を一部だけ無効化した上で、条件を満たす方向に力をかけ、力の方向に構造を変化させます。これによって、様々な内部座標方向に構造を変化させることができるだけでなく、複数の構造変化を同時に適用したり、逐次的に適用したりすることができます。

構造を変化させる方向の定義は、ユーザーが定義することもできますが、よく使われる内部座標の定義に基づいて構造を変化させることもできます。事前に実装されている内部座標を使用して、構造を変化させる例を以下に示します。


## 初期設定

* 必要なライブラリ等の準備を行います。


```python
# If you have already installed matlantis-features, you can skip this.
!pip install 'matlantis-features&gt;=0.16.0'
```




```python
import logging
import matlantis_features
import pfp_api_client
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

model_version = "v8.0.0"
calc_mode = EstimatorCalcMode.PBE

estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
calc = ASECalculator(estimator)

print(f"matlantis_features: {matlantis_features.__version__}")
print(f"pfp_api_client: {pfp_api_client.__version__}")
print(f"current pfp model version: {estimator.model_version}")

logger = logging.getLogger("matlantis_features")
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.setLevel(logging.INFO)
logger.addHandler(handler)
```

    matlantis_features: 1.0.1
    pfp_api_client: 2.0.1
    current pfp model version: v8.0.0


```python
import nglview as nv
import numpy as np
from ase import Atoms
from ase.visualize import view
from ase.build import molecule
from ase.io import Trajectory

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

from matlantis_features.features.common.opt import LBFGSASEOptFeature
from matlantis_features.features.reaction.rest_scan import (
    RestScanFeature,
    RestScanFeatureResult,
    ScanDihedral,
    ScanDistance,
    calculate_dihedral,
)
from matlantis_features.features.reaction.reaction_string import ReactionStringFeature
from matlantis_features.utils.calculators import pfp_estimator_fn
```





```python
def view_with_atomindex(atoms: Atoms, label_color: str = "blue", label_scale: float = 1.0):
    v = nv.show_ase(atoms, viewer="ngl")
    v.add_label(
        color=label_color, labelType="text",
        labelText=[atoms[i].symbol + str(i) for i in range(atoms.get_global_number_of_atoms())],
        zOffset=1.0, attachment='middle_center', radius=label_scale
    )
    return v
```


## 二つのscanを直列繋ぎで実行する例

ここでは、エタン分子を例に、Restscanを用いて (i)C-C 結合を切断し、(ii)結合切断後の2つの CH3 グループ間の二面角を回転させる 方法について説明します。
この場合、2つのscanを連続して準備する必要があります。最初のscanは結合を切断し、2つ目のscanは、結合の切断後に二面角を回転させるためのものです。

### 初期構造の作成

まずエタン分子を作り、`LBFGSASEOptFeature` を使って構造最適化を行います。
NGL Viewerを用いて構造と原子の番号を可視化してみます。


```python
atoms = molecule("C2H6")

opt_feature = LBFGSASEOptFeature(
    fmax=0.01, n_run=1000, show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version=model_version, calc_mode=calc_mode),
)
opt_result = opt_feature(atoms)
```


      0%|          | 0/1001 [00:00&lt;?, ?it/s]


```python
view_with_atomindex(atoms)
```


    NGLWidget()


### RestScan 計算

以下の例では原子0と原子1の距離を2Åまで伸ばした後、2-0-1-5で定義される二面角を90度回転させようとします。
各原子の番号は上の図で示されている通りです。

ただし、二面角を90度よりもさらに大きい角度回転させるとより安定になります。この追加の回転も自動で行われます。
幾何的には120度回転させると対称性が高く安定になると予想されますが、この例ではC-C距離を伸ばしてから回転した影響で、120度回転する前にforceが十分小さくなってしまい止まります。

現在、ScanDihedral および ScanDistance の2種類のscanが提供されています。

`ScanDistance` および `ScanDihedral` には、以下の引数が指定できます。

|Parameter|Description|
|---|---|
|indices|スキャンされる原子のインデックス。 `ScanDistance`では2つ、`ScanDihedral`では4つのインデックスを指定。|
|direction|大きな距離or角度に向かってスキャンする場合は1、小さい距離or角度に向かってスキャンする場合は-1に設定。|
|destination|距離(`ScanDistance`)または二面角(`ScanDihedral`)がこの値に達するとスキャンは停止。|
|exceed|このクラスでは通常、スキャン力は目的地に到達するとオフになります。この値によりスキャンが目的地を超えることができる距離を設定できます。オプティマイザが感じる力場を滑らかにするために使われます。|
|fmax|拘束スキャン力。|
|current|(`ScanDihedral`のみ)計算再開時に指定するパラメータ。通常は指定する必要がありません。|


```python
current_dihedral = calculate_dihedral(opt_result.atoms.ase_atoms, [2, 0, 1, 5])
scans = [
    [ScanDistance(indices=[0, 1], direction=1, destination=2.0, exceed=0.1, fmax=0.05)],
    [
        ScanDihedral(
            indices=[2, 0, 1, 5], direction=-1, destination=current_dihedral - 90.0, exceed=1.0, fmax=0.05
        )
    ],
]
```


`RestScanFeature` は `Scan` のリストのリストを入力として受け取ります。 また、Feature に `estimator_fn` 引数を渡すことで、matlantis-features で使用される PFP 計算モードやモデルバージョンを指定できます。初期状態と最終状態の間のイメージを保存する場合は `trajectory` 引数を指定してください。指定しない場合、イメージは保存されません。


```python
scan_feature = RestScanFeature(
    scans=scans,
    trajectory="scan.traj",
    estimator_fn=pfp_estimator_fn(model_version=model_version, calc_mode=calc_mode),
    show_progress_bar=True, show_logger=True, logger_interval=10
)

scan_result = scan_feature(opt_result.atoms.ase_atoms, fmax=0.01, n_run=1000)
```


    Scan steps (ends when converged):   0%|          | 0/1000 [00:00&lt;?, ?it/s]


    Running scan 1/2
    steps:     0  energy：-3.88 eV/atom  max_force:  0.05 eV/Ang
    steps:    10  energy：-3.88 eV/atom  max_force:  0.05 eV/Ang
    steps:    20  energy：-3.83 eV/atom  max_force:  0.10 eV/Ang
    steps:    30  energy：-3.81 eV/atom  max_force:  0.05 eV/Ang
    steps:    40  energy：-3.79 eV/atom  max_force:  0.06 eV/Ang
    steps:    50  energy：-3.75 eV/atom  max_force:  0.05 eV/Ang
    steps:    60  energy：-3.72 eV/atom  max_force:  0.05 eV/Ang
    steps:    70  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    Running scan 2/2
    steps:    80  energy：-3.71 eV/atom  max_force:  0.08 eV/Ang
    steps:    90  energy：-3.71 eV/atom  max_force:  0.09 eV/Ang
    steps:   100  energy：-3.71 eV/atom  max_force:  0.09 eV/Ang
    steps:   110  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang
    steps:   120  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang
    steps:   130  energy：-3.71 eV/atom  max_force:  0.31 eV/Ang
    steps:   140  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang
    steps:   150  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang
    steps:   160  energy：-3.71 eV/atom  max_force:  0.11 eV/Ang
    steps:   170  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang
    steps:   180  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang
    steps:   190  energy：-3.71 eV/atom  max_force:  0.08 eV/Ang
    steps:   200  energy：-3.71 eV/atom  max_force:  0.04 eV/Ang
    steps:   210  energy：-3.71 eV/atom  max_force:  0.02 eV/Ang
    steps:   220  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    steps:   230  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    steps:   240  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    steps:   250  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    steps:   260  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    steps:   270  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    steps:   280  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    steps:   290  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang
    steps:   300  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


    steps:     0  energy：-3.88 eV/atom  max_force:  0.05 eV/Ang


    steps:    10  energy：-3.84 eV/atom  max_force:  0.05 eV/Ang


    steps:    20  energy：-3.82 eV/atom  max_force:  0.05 eV/Ang


    steps:    30  energy：-3.75 eV/atom  max_force:  0.06 eV/Ang


    steps:    40  energy：-3.72 eV/atom  max_force:  0.05 eV/Ang


    Running scan 2/2


    steps:    50  energy：-3.71 eV/atom  max_force:  0.08 eV/Ang


    steps:    60  energy：-3.71 eV/atom  max_force:  0.08 eV/Ang


    steps:    70  energy：-3.71 eV/atom  max_force:  0.09 eV/Ang


    steps:    80  energy：-3.71 eV/atom  max_force:  0.11 eV/Ang


    steps:    90  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang


    steps:   100  energy：-3.71 eV/atom  max_force:  0.11 eV/Ang


    steps:   110  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang


    steps:   120  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang


    steps:   130  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang


    steps:   140  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang


    steps:   150  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang


    steps:   160  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang


    steps:   170  energy：-3.71 eV/atom  max_force:  0.10 eV/Ang


    steps:   180  energy：-3.71 eV/atom  max_force:  0.06 eV/Ang


    steps:   190  energy：-3.71 eV/atom  max_force:  0.02 eV/Ang


    steps:   200  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


    steps:   210  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


    steps:   220  energy：-3.71 eV/atom  max_force:  0.02 eV/Ang


    steps:   230  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


    steps:   240  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


    steps:   250  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


    steps:   260  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


    steps:   270  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


    steps:   280  energy：-3.71 eV/atom  max_force:  0.01 eV/Ang


### 計算結果

`RestScanFeature` の結果は `RestScanFeatureResult` のインスタンスに保存されます。`RestScanFeatureResult` は4つのattributesを持ちます。

|Parameter|Description|
|---|---|
|trajectory_logger|メモリ上に格納されたトラジェクトリロガー。|
|traj_path|トラジェクトリファイルのパス。|
|converged|計算が正常に終わったかどうか。|
|warnflag|計算が収束するまえに終了した場合の理由。 1: 最大ステップ数へ到達した, 2: ExitEnergyに到達した。|


Restscan の計算結果を集約し、初期状態と最終状態の間のイメージを取得するには、 `RestScanFeatureResult.get_images` メソッドを使用できます。このメソッドには、以下の引数が指定できます。

|Parameter|Description|
|---|---|
|kink_dx (float, optional)|最大距離が `dx` 以下のimageは同一とみなし、省略されます。デフォルト値は0.1。|
|local_minima_de (float or None, optional)|`de` が与えられた場合、最終的なlocal minimumを集約して返します。ここで， `de` はlocal minimumとみなす閾値です。`None` を指定すると、生成されたすべてのimageを返します。デフォルト値は `None` です。|


```python
nv.show_asetraj(scan_result.get_images(kink_dx=0.1))
```


    NGLWidget(max_frame=11)


保存されたファイルからトラジェクトリをロードできます。


```python
nv.show_asetraj(Trajectory("scan.traj"))
```


    NGLWidget(max_frame=306)


`RestScanFeatureResult.plot` メソッドを使用して、トラジェクトリに沿ったエネルギープロファイルをプロットすることもできます。


```python
fig = scan_result.plot(plt_name=None, kink_dx=0.1)
fig
```


## scanを最初のlocal minimumまで実行する例

### 初期構造の作成

再び例としてエタンを使用します。


```python
atoms = molecule("C2H6")

opt_feature = LBFGSASEOptFeature(
    fmax=0.01, n_run=1000, show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version=model_version, calc_mode=calc_mode),
)
opt_result = opt_feature(atoms)
```


      0%|          | 0/1001 [00:00&lt;?, ?it/s]


### RestScan計算

よくある事例として、結合の長さや角度をエネルギーの山を越えるまで変化させたいが、何度回転させればいいか or 何Å伸ばせば良いかは不明である、というケースがあります。このような場合、回転させる角度や伸ばす距離を大きめに設定しておき、山を越えた時にscanを抜けるようにすると得たい経路が得られる可能性があります。

以下の例では、とりあえず360度回転するように設定しておき、山を越えたら抜けるように設定しています。


```python
current_dihedral = calculate_dihedral(opt_result.atoms.ase_atoms, [2, 0, 1, 5])
scans = [
    [
        ScanDihedral(
            indices=[2, 0, 1, 5], direction=-1, destination=current_dihedral - 360.0, exceed=1.0, fmax=0.05
        )
    ],
]
```


Restscan を早期終了させるための引数として、以下のものがあります。

|Parameter|Description|
|-|-|
|energy_stop (float or None, optional)|  系のpotential energyがこの値より高い場合に中断します。デフォルトはNoneです。|
|local_minima_de (float, optional)| local maximaに対し、この値以上低いlocal minimumをlocal_minima_nの判定においてカウントする対象とします。|
|local_minima_n (float, optional)| トラジェクトリがn個目のlocal minumumに到達したとき終了します。|

この例では、初めての局所最小値を見つけるために `local_minima_n=1` を設定しました。さらに、`energy_stop` 引数も使用して、現実的でない構造（ポテンシャルエネルギーが非常に高い構造）に変化するのを防ぎました。システムのポテンシャルエネルギーが `energy_stop` を超えた場合は、計算が中止されます。


```python
scan_feature = RestScanFeature(
    scans=scans,
    logfile=None,
    estimator_fn=pfp_estimator_fn(model_version=model_version, calc_mode=calc_mode),
    show_progress_bar=True, show_logger=True, logger_interval=10
)

de = 0.01
energy_stop = opt_result.energy_log[-1] + 2.0

scan_result = scan_feature(
    opt_result.atoms,
    fmax=0.01,
    n_run=1000,
    energy_stop=energy_stop,
    local_minima_de=de,
    local_minima_n=1
)
```


    Scan steps (ends when converged):   0%|          | 0/1000 [00:00&lt;?, ?it/s]


    Running scan 1/1
    steps:     0  energy：-3.88 eV/atom  max_force:  0.05 eV/Ang
    steps:    10  energy：-3.88 eV/atom  max_force:  0.09 eV/Ang
    steps:    20  energy：-3.88 eV/atom  max_force:  0.09 eV/Ang
    steps:    30  energy：-3.88 eV/atom  max_force:  0.13 eV/Ang
    steps:    40  energy：-3.88 eV/atom  max_force:  0.14 eV/Ang
    steps:    50  energy：-3.87 eV/atom  max_force:  0.10 eV/Ang
    steps:    60  energy：-3.87 eV/atom  max_force:  0.11 eV/Ang
    steps:    70  energy：-3.87 eV/atom  max_force:  0.11 eV/Ang
    steps:    80  energy：-3.87 eV/atom  max_force:  0.12 eV/Ang
    steps:    90  energy：-3.87 eV/atom  max_force:  0.12 eV/Ang
    steps:   100  energy：-3.87 eV/atom  max_force:  0.12 eV/Ang
    steps:   110  energy：-3.88 eV/atom  max_force:  0.10 eV/Ang
    steps:   120  energy：-3.88 eV/atom  max_force:  0.10 eV/Ang
    steps:   130  energy：-3.88 eV/atom  max_force:  0.09 eV/Ang
    steps:   140  energy：-3.88 eV/atom  max_force:  0.10 eV/Ang
    1 th local minima reached.


    steps:     0  energy：-3.88 eV/atom  max_force:  0.05 eV/Ang


    steps:    10  energy：-3.88 eV/atom  max_force:  0.09 eV/Ang


    steps:    20  energy：-3.88 eV/atom  max_force:  0.10 eV/Ang


    steps:    30  energy：-3.88 eV/atom  max_force:  0.14 eV/Ang


    steps:    40  energy：-3.88 eV/atom  max_force:  0.10 eV/Ang


    steps:    50  energy：-3.87 eV/atom  max_force:  0.10 eV/Ang


    steps:    60  energy：-3.87 eV/atom  max_force:  0.19 eV/Ang


    steps:    70  energy：-3.87 eV/atom  max_force:  0.11 eV/Ang


    steps:    80  energy：-3.87 eV/atom  max_force:  0.12 eV/Ang


    steps:    90  energy：-3.87 eV/atom  max_force:  0.13 eV/Ang


    steps:   100  energy：-3.87 eV/atom  max_force:  0.12 eV/Ang


    steps:   110  energy：-3.87 eV/atom  max_force:  0.11 eV/Ang


    steps:   120  energy：-3.88 eV/atom  max_force:  0.14 eV/Ang


    steps:   130  energy：-3.88 eV/atom  max_force:  0.11 eV/Ang


    steps:   140  energy：-3.88 eV/atom  max_force:  0.10 eV/Ang


    steps:   150  energy：-3.88 eV/atom  max_force:  0.10 eV/Ang


    steps:   160  energy：-3.88 eV/atom  max_force:  0.10 eV/Ang


    1 th local minima reached.


### 計算結果

初期状態と最終状態の間のトラジェクトリを表示してみましょう。


```python
nv.show_asetraj(scan_result.get_images(kink_dx=0.1))
```


    NGLWidget(max_frame=10)


トラジェクトリに沿ったエネルギープロファイルも見てみましょう。


```python
fig = scan_result.plot(plt_name=None, kink_dx=0.1, local_minima_de=de)
fig.show()
```


## ReactionStringによる反応経路最適化

RestScanで構造を変化させることができますが、反応速度を見積もるために活性化障壁を求めたい場合には、RestScanの結果をそのまま使うことは出来ません。
これはRestScanの結果は多くの場合に最適な反応経路上を通っておらず、ポテンシャルエネルギー面の一次の鞍点となるような遷移状態も通っていないためです。
活性化障壁が必要な場合、ReactionStringのような反応経路と遷移状態を最適化するための手法を使う必要があります。
RestScanを用いて得られた大まかな経路をReactionStringで最適化することで、最小エネルギー経路と遷移状態を得ることができます。


```python
scan_images = [atoms.copy() for atoms in scan_result.get_images()]
```


```python
rp = ReactionStringFeature(
    optimize_is=True,
    optimize_fs=True,
    fmax_rp=0.05,
    fmax_ts=0.05,
    fmax_rd=0.05,
    fmax_eq=0.01,
    estimator_fn=pfp_estimator_fn(model_version=model_version, calc_mode=calc_mode),
)
rp_result = rp([scan_images])
```

    &lt;frozen reactionstring.reparametrize.search&gt;:151: RuntimeWarning:

    Method 'bounded' does not support relative tolerance in x; defaulting to absolute tolerance.



### 計算結果

得られた反応経路のトラジェクトリを表示してみましょう。


```python
nv.show_asetraj(scan_result.get_images(kink_dx=0.1))
```


    NGLWidget(max_frame=10)


反応経路に沿ったエネルギープロファイルも見てみましょう。
ReactionStringで得られる活性化エネルギーは最小エネルギー経路のエネルギーの最大値となっており、計算化学で得られた活性化障壁として用いることが出来ます。


```python
rp_result.plot_energy()
```


