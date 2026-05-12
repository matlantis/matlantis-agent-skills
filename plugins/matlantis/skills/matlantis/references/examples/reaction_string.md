# Reaction string method

`ReactionStringFeature` はCI-NEB等と同分類の反応経路探索アルゴリズムです。String法をベースにTSを高精度かつ自動的に計算するためいくつかの改善を行いました。

この手法の主な利点は以下のようになっています。
・反応経路最適化とTS最適化が自動で進行します。
・反応経路の長さに応じて、計算するイメージの数を自動的に調整します。
・Dimer法やSellaに近い精度でTSの最適化を行います。
・障壁が複数あるような反応経路も適切に取り扱うことが出来ます。
・kink(反応経路のねじれ)の解消機能が実装されており、kinkによる計算失敗を防ぎます。

本exampleではPt表面上でCH3がCH2とHに解離する反応の経路探索を扱います。


```python
# If you have already installed matlantis-features, you can skip this.
# Version 0.18.0 or later is required to run this notebook.
!pip install 'matlantis-features&gt;=0.22.0'
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
import pickle
from datetime import timedelta
from math import inf
from pathlib import Path
from typing import List

import nglview as nv
import numpy as np
from ase import Atoms
from ase.build import add_adsorbate, fcc111, molecule
from ase.constraints import FixAtoms
from ase.io import Trajectory
from matlantis_features.features.reaction import ReactionStringFeature
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator
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


## 初期構造の作成
以下のセルでは、Ptのスラブを作成し、`ini` にCH3、 `fin` にCH2とHを配置し、可視化します。Ptの原子配置は固定します。


```python
ini = fcc111("Pt", (4, 4, 4), vacuum=20.0)
ini.constraints.append(FixAtoms(mask=ini.positions[:, 2] &lt; ini.cell[2, 2] / 2))
position = tuple((ini.cell[0] + ini.cell[1]) / 2)[:2]
add_adsorbate(ini, molecule("CH3"), 2.0, position=position)
```


```python
fin = fcc111("Pt", (4, 4, 4), vacuum=20.0)
fin.constraints.append(FixAtoms(mask=fin.positions[:, 2] &lt; fin.cell[2, 2] / 2))
ads = molecule("CH3")
del ads[3]
position_ch2 = tuple((fin.cell[0] + fin.cell[1]) / 2)[:2]
position_h = tuple((fin.cell[0] + fin.cell[1]) / 4)[:2]
add_adsorbate(fin, ads, 2.0, position=position_ch2)
add_adsorbate(fin, "H", 1.4, position=position_h)
```


```python
def view_ngl(images: List[Atoms]) -&gt; nv.widget.NGLWidget:
    atoms = images[0]
    view = nv.show_asetraj(images)
    view.add_label(
        color="red",
        scale=1.3,
        labelType="text",
        labelText=[str(i) for i in range(len(atoms))],
        zOffset=2.0,
        attachment="middle_center",
    )
    return view
```


```python
view_ngl([ini, fin])
```


    NGLWidget(max_frame=1)


# `ReactionStringFeature` の入力について

`ReactionStringFeature` オブジェクトは `List[List[MatlantisAtoms | AseAtoms]]` 型の入力を受け付けます。これは後述する出力の型と同じであり、複数の障壁からなる反応経路の概形を入力として取ることが出来るようになっています。

入力を真の反応経路に近くすることで計算時間を大きく削減することが出来ますが、ISとFSの線形補完から最適化を開始しても問題がない場合、`[[IS, FS]]` という形でISとFSのみ指定することも出来ます。 本exampleではISとFSの線形補完から開始します。


# 引数について

fmaxから始まる4つの引数で精度と速度のトレードオフを調整します。
- `fmax_rp` はreaction path全体のfmaxを指定します。
- `fmax_ts` はTSのfmaxを指定します。
- `fmax_rd` は最も高いTS (rate determining step）のfmaxを指定します。
- `fmax_eq` は全てのlocal minimaの最適化のfmaxを指定します。

また、それぞれの値は`inf`にすることで収束しなくても良いという条件を指定することも出来ます。


# 出力について
`ReactionStringFeature` の計算結果は `ReactionStringFeatureResult` という型に格納されます。この返り値は4つのattributeを持ちます:

- `num_segments(int)` は障壁の数を示します。
- `reaction_string_images (List[List[MatlantisAtoms]])` 反応経路を示します。反応経路は途中の全てのEQ (stable point)で分けられ、一つの `List[MatlantisAtoms]` は一つ以下のTS (saddle point)を含み、両端がEQであるように分けられます。これを本ノートブックでは sub-RPと呼ぶこととします。ある sub-RP の最後の構造と次の sub-RP の最初の構造は必ず一致します。
- `transition_state_indices (List[Optional[int]])` は各山の頂点、すなわちTSが何番目の `MatlantisAtoms` であるかを示します。経路は一つの山に一つのTSが含まれるように分けられるため、基本的に単一の整数値で山のインデックスを示します。しかし、TSが計算されなかった場合は整数値の代わりに `None` が格納されます。TSが計算されない原因として、fmaxを `inf` にする等で計算をスキップした場合やタイムアウトした場合、経路があまりにも短すぎた場合や障壁が `local_extrema_tol` より小さい場合等がありえます。
- `imaginary_eigenvectors (List[np.ndarray])` は計算されたそれぞれの山の固有ベクトルが格納されます。


## 全てのTSを求める場合
`ReactionStringFeature` オブジェクトを作成してみましょう。以下の例では、全ての経路をfmax=0.05eV/Aの基準で収束させ、全てのTSをfmax=0.01eV/Aの基準で収束させるように設定しています。

作成したインスタンスに `ini`, `fin` を入力することで反応経路探索が開始します。この例では10～20分程度時間がかかります。


```python
%%time
rp = ReactionStringFeature(
    fmax_rp=0.05,
    fmax_ts=0.01,
    fmax_rd=0.01,
    fmax_eq=0.001,
    timeout=None,
    estimator_fn=estimator_fn,
)
rp_result = rp([[ini, fin]])
```

    &lt;frozen reactionstring.reparametrize.search&gt;:151: RuntimeWarning: Method 'bounded' does not support relative tolerance in x; defaulting to absolute tolerance.


    CPU times: user 33.5 s, sys: 1.77 s, total: 35.3 s
    Wall time: 4min 23s


```python
# It visualizes the result. x-axis is an index of a reaction path and y-axis is eV.
# `reaction_string_images` are segmented by stable points (EQ, red points in the figure).
# `plot_energy()` show the entire reaction path connecting all paths by EQ.
rp_result.plot_energy()
```


## JSONファイルへの結果とパラメータの保存と読み込み

シミュレーション結果と呼び出し引数はJSONファイルに保存して、後から読み込むことができます。
なお、いくつかの `Feature` の初期化パラメータは読み込みに対応していません。
具体的には、`ReactionStringFeature` の `interesting` 、`timeout` 、`observers_rp` 、 `observers_eq` 、`estimator_fn` パラメーターはサポートされていません。
これらのパラメーターを `Feature` 初期化時に設定したとしても、 `saved_result.feature` に入っている `Feature` オブジェクトのパラメーターはデフォルト値となっています。


```python
# Save the result to json (some parameters cannot be serialized).
rp_result.save("result.json")
```

    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/utils/save_util.py:62: UserWarning:

    The parameter estimator_fn does not support saving and reloading. It will be ignored upon saving and reloading.



    True


```python
# Reloading the result
# UserWarning will be shown because estimator_fn does not support saving and reloading.
from matlantis_features.features.reaction.reaction_string import ReactionStringFeatureResult
saved_result = ReactionStringFeatureResult.load("result.json")
```

    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/reaction/reaction_string.py:278: UserWarning:

    Model version is not specified, default 'v8.0.0' will be used.Please specify it by providing an estimator_fn or setting the 'MATLANTIS_PFP_MODEL_VERSION' environment variable.

    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/reaction/reaction_string.py:278: UserWarning:

    Calculation mode is not specified, default 'EstimatorCalcMode.PBE' will be used.Please specify it by providing an estimator_fn or setting the 'MATLANTIS_PFP_CALC_MODE' environment variable.



```python
# Inspecting the result
print("Number of segments:", saved_result.num_segments)
print("Is converged:", saved_result.converged)

# PFP parameters (model version, calc mode, etc.)
print("Calculator info: ", saved_result.calculator_info)

# Package version information
print("Version info: ", saved_result.version_info)
```

    Number of segments: 2
    Is converged: True
    Calculator info:  {'model_version': 'v8.0.0', 'calc_mode': 'pbe', 'method_type': 'pfvm_d3_pfvm'}
    Version info:  {'matlantis-feature-result-version': '0.0.1', 'ase': '3.26.0', 'pfp-api-client': '2.0.1', 'matlantis-features': '1.0.1', 'python': '3.13'}


```python
# Rerunning the experiment
# result_rerun = saved_result.feature(**saved_result.call_params)
```


```python
# Saving a result as traj file.
rp_result.to_traj("reacstringCH3onPT.traj")
```


```python
# Visualizing a result as GIF. It takes minutes.
rp_result.to_gif("react.gif")
```



![png](output_23_0.png)



分解された反応ごとのsub-RPのリスト(`list[list[MatlantisAtoms]]`)が `rp_result.reaction_string_images` に格納されます。以下のように`AseAtoms`のリストが取り出せます。

各反応の終状態と次の反応の始状態は同じであることに注意してください。


```python
rs = rp_result.reaction_string_images
atomslist = [rs[0][0].ase_atoms]
atomslist += [matatoms.ase_atoms for reac in rs for matatoms in reac[1:]]
print(len(atomslist))
view_ngl(atomslist)
```

    45


    NGLWidget(max_frame=44)


## 最も高いTSだけ精度を上げたい場合
`fmax_rd` を指定することで、最も高いTS(rate-determining step)だけ精度を上げることができます。以下の例では、最も高いTSのみ収束基準を0.005eV/Aとし、他のすべてのTSの収束基準を0.01eV/Aとしています。


```python
%%time
rp = ReactionStringFeature(
    fmax_rp=0.05,
    fmax_ts=0.01,
    fmax_rd=0.005,
    fmax_eq=0.001,
    timeout=None,
    estimator_fn=estimator_fn,
)
rp_result = rp([[ini, fin]])
```

    &lt;frozen reactionstring.reparametrize.search&gt;:151: RuntimeWarning:

    Method 'bounded' does not support relative tolerance in x; defaulting to absolute tolerance.



    CPU times: user 24.3 s, sys: 2.03 s, total: 26.3 s
    Wall time: 6min 4s


```python
rp_result.plot_energy()
```


## 一番高いTSのみを求める場合
以下の例では、全ての経路をfmax=0.05eV/Aの基準で収束させ、最も高いTSのみfmax=0.01eV/Aの基準で収束させます。最も高いTS以外のTSは一切収束させません。山の数が多い時は少し速くなることもあります。

`ReactionStringFeatureResult.transition_state_indices` には各山のTSに対応するAtomsのindexが格納されますが、一番高い山以外ではTSを計算しないため整数値のindexの代わりにNoneが格納されます。このため、計算がスムーズに進んだ場合は整数が1つだけ格納され、残りはNoneになります。しかし、最適化の途中で山の高さが変化し、最も高い山が切り替わることがあります。このような場合にはTS計算を行なった全ての山の情報が格納されるため、Noneでない数が複数入ることもあります。


```python
%%time
rp = ReactionStringFeature(
    fmax_rp=0.05,
    fmax_ts=inf,
    fmax_rd=0.01,
    fmax_eq=0.001,
    timeout=None,
    estimator_fn=estimator_fn,
)
rp_result = rp([[ini, fin]])
```

    CPU times: user 23 s, sys: 2.05 s, total: 25.1 s
    Wall time: 10min 23s


# TSが必要ない場合
NEBのように大体の経路だけが必要でTSが必要ない場合、`fmax_ts`と`fmax_rd`の両方を`inf`にすることで計算時間が大幅に短縮されることがあります。


```python
%%time
rp = ReactionStringFeature(
    fmax_rp=0.05,
    fmax_ts=inf,
    fmax_rd=inf,
    fmax_eq=0.001,
    timeout=None,
    estimator_fn=estimator_fn,
)
rp_result = rp([[ini, fin]])
```

    CPU times: user 13.2 s, sys: 1.14 s, total: 14.4 s
    Wall time: 6min 3s


## 時間がかかる計算を途中で抜けたい場合
時間がかかる計算を途中で抜けたい場合、timeout引数にtimedeltaで時間を入力します。 ここでは1分で抜ける例を示します。
事前計算などのいくつかの処理は停止する対象に含まれていません。このため、指定した時間に実行が終了しないことがあります。

## 並行数を自分で制御したい場合
`max_workers` 引数で並行実行するスレッドの数を変更できます。デフォルトの並行数の最大値は `max(1, min(50, 10000 // n_atoms))` という計算式によって決定されます。
ここで、 `n_atoms` は始状態の原子数です。

## 途中の探索結果をファイルに出力したい場合
`dump_frequency` および `dump_directory` を設定すると探索中に反応経路を `dump_directory` に出力できます。
`dump_directory` 以下に `2024-10-30-04:41:11` のように dump 処理が実行された時刻を表すフォルダが作成され、i 番目の sub-RP のトラジェクトリが `{i}.traj` という形式で出力されます。

```
dump-reactionstring-results
├── 2024-10-30-04:41:11
│&nbsp;&nbsp; └── 1.traj
├── 2024-10-30-04:42:16
│&nbsp;&nbsp; ├── 1.traj
│&nbsp;&nbsp; └── 2.traj
├── 2024-10-30-04:43:17
│&nbsp;&nbsp; ├── 1.traj
│&nbsp;&nbsp; └── 2.traj
└── 2024-10-30-04:44:17
    ├── 1.traj
    └── 2.traj
```


```python
%%time
rp = ReactionStringFeature(
    fmax_rp=0.05,
    fmax_ts=inf,
    fmax_rd=0.01,
    fmax_eq=0.001,
    timeout=timedelta(minutes=1),
    estimator_fn=estimator_fn,
    max_workers=20,
    dump_directory=Path("dump-reactionstring-results"),
    dump_frequency=timedelta(seconds=20),
)
rp_result = rp([[ini, fin]])
```

    CPU times: user 831 ms, sys: 47.2 ms, total: 879 ms
    Wall time: 1min


    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/reaction/reaction_string.py:466: UserWarning:

    calculator is not attached.



# 欲しい反応の種類を自分で指定したい場合

今回取り扱う反応では、64番目の原子と67番目の原子の解離反応に興味があります。しかし、ISとFSを簡素な本法で作成したために、不要な表面拡散反応を多数含んでしまい計算に時間がかかってしまっています。

注目したい山が経路の中で一番高いことがわかっていれば `fmax_rd` を `math.inf` にするなどで計算時間が長くなることを回避することができます。しかし、エネルギーに関する事前知識がなく、反応の形状だけで判断したい時もあります。このようなケースに対応するため、欲しい反応だけを精密に計算し、残りの反応はほとんど計算しないようにするためのオプションが用意されています。

この例では64番目と67番目の原子の間のCH結合の解離のみを高精度に計算します。注目すべき反応を抽出する部分は、反応の形状に応じて定義する必要があります。したがって、ユーザーによりコードを書いていただく必要がありますが、プログラマブルな形で提供されています。以下の`InterestingBondChange` では、 `i1` 番目の原子と `i2` 番目の原子の距離を監視し、結合距離が `r` をまたいで変化した山に対してのみ `True` を割り当て、それ以外の山に `False` を割り当てるように記述されています。


```python
class InterestingBondChange:
    def __init__(self, i1: int, i2: int, r: float) -&gt; None:
        self.i1 = i1
        self.i2 = i2
        self.r = r

    def __call__(self, images: List[List[Atoms]]) -&gt; List[bool]:
        result = []
        for im in images:
            ini_connected = (
                np.linalg.norm(im[0][self.i1].position - im[0][self.i2].position)
                &lt; self.r
            )
            fin_connected = (
                np.linalg.norm(im[-1][self.i1].position - im[-1][self.i2].position)
                &lt; self.r
            )
            if ini_connected and (not fin_connected):
                result.append(True)
            else:
                result.append(False)
        return result
```


```python
%%time
interesting = InterestingBondChange(i1=64, i2=67, r=1.6)
rp = ReactionStringFeature(
    fmax_rp=0.05,
    fmax_ts=0.01,
    fmax_rd=0.01,
    fmax_eq=0.001,
    timeout=None,
    interesting=interesting,
    estimator_fn=estimator_fn,
)
rp_result = rp([[ini, fin]])
```

    CPU times: user 17.3 s, sys: 1.66 s, total: 18.9 s
    Wall time: 6min 53s

