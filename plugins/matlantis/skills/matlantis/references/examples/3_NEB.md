# 準備

PFPを利用するためには、`pfp_api_client`という専用のパッケージが必要です。

初めてPFPを利用する際には、以下のコマンドを実行してください。

一度実行すれば、再度このノートブックを開いた場合や、他のノートブックを利用する場合に再実行する必要はありません。


```python
!pip install pfp-api-client
```


# NEB法による反応経路解析

次は毛色の違う例として、Nudged Elastic Band (NEB)法による反応経路解析を行ってみましょう。

まずは、今回のチュートリアルで使うライブラリをimportしておきます。


```python
import numpy as np
from IPython.display import Image
import matplotlib.pyplot as plt

import ase
from ase import Atoms
from ase.visualize import view
from ase.optimize import BFGS
from ase.optimize import FIRE
from ase.neb import NEB
from ase.build.rotate import minimize_rotation_and_translation

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator()
calculator = ASECalculator(estimator)
```


今回は、クルチウス転位と呼ばれる有機化学反応の一例を題材に進めていきます。

今回の題材にする構造を設定します。NEB法では、ある注目したい反応についてその前後の構造からスタートします。

今回のチュートリアルでは単一のファイルで完結させるため、直接構造を書き下しておきます。


```python
react = ase.Atoms(
    symbols="C2NON2H3",
    positions = [
        [ 2.12, -0.48,  0.00],
        [ 0.76,  0.16, -0.00],
        [-0.28, -0.79, -0.00],
        [ 0.57,  1.35,  0.00],
        [-1.42, -0.28,  0.00],
        [-2.48,  0.07,  0.00],
        [ 2.67, -0.14, -0.88],
        [ 2.06, -1.57, -0.00],
        [ 2.67, -0.14,  0.88],
    ],
)

prod = ase.Atoms(
    symbols="C2NON2H3",
    positions=[
        [ 2.30, -0.87,  0.00],
        [ 0.68,  0.92, -0.00],
        [ 0.99, -0.25,  0.00],
        [ 0.25,  2.03,  0.00],
        [-2.34,  0.11,  0.00],
        [-3.13, -0.67, -0.00],
        [ 2.87, -0.57, -0.89],
        [ 2.19, -1.96,  0.00],
        [ 2.87, -0.57,  0.89],
    ],
)
```


まずは構造最適化計算を行います。構造最適化計算のやり方は今までと同様です。

また、後の作業のため、反応前後の2つの構造がなるべく近くなるように位置を移動させておきます。


```python
react.calc = calculator
opt = BFGS(react)
opt.run(fmax=0.01)

prod.calc = calculator
opt = BFGS(prod)
opt.run(fmax=0.01)

minimize_rotation_and_translation(react, prod)
```


可視化のスクリプトも下に置いておきます。可視化ウィンドウ内のバーを操作すると反応前後の構造を切り替えできます。


```python
v = view([react, prod], viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


ここからがNEB法の出番です。NEB法は、反応の前後の構造を離散的に補完した構造を多数用意し、それら複数の構造を位相空間上で互いに結合させつつ構造最適化を行うことで反応経路を探す手法です。

中間点の補完やNEBの実行そのもののはASEの組み込みのものを使うことができます。一点注意することとして、NEB 計算を並列に実行するために、各構造ごとに　Calculator のインスタンスを作成します。また、NEBの引数に`allow_shared_calculator=False, parallel=True` を指定します。

NEBは複数の構造を内部で使う手法のため、MD計算と比較してもステップあたりの計算時間は大きくなります。それでもPFPを使えば、この程度の計算であれば数分で終わらせることができます。


```python
images = [react.copy()]
images += [react.copy() for i in range(30)]
images += [prod.copy()]
for image in images:
    estimator = Estimator()
    calculator = ASECalculator(estimator)
    image.calc = calculator
neb = NEB(images, k=0.1, climb=True, allow_shared_calculator=False, parallel=True)
neb.interpolate()
opt = FIRE(neb)
status = opt.run(fmax=0.05, steps=500)
```


もしログを抑制したい場合は、`FIRE(neb, logfile=None)`のような指定をするとデフォルトのoptのログが表示されなくなります。

結果を見てみましょう。Notebook環境は数値データのグラフ可視化などにも向いています。以下はmatplotlibを用いてエネルギーの軌跡を可視化したものです。

エネルギーが一度上がった後で下がる反応が見えていると思います。このエネルギーが最大の点(transition state)と左右の安定点それぞれとのエネルギー差が活性化エネルギーに相当します。


```python
energies = [image.get_total_energy() for image in images]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(energies)
ax.set_xlabel("replica")
ax.set_ylabel("energy")
fig.show()
```


結果を見てみましょう。構造についても同様にNotebookに可視化してもいいですが、画像をファイルとして保存することもできます。今回は試しにそちらを実行してみましょう。

今回はNEBのため、結果は構造1つではなく一組の軌跡として出力されています。ASEからアニメーションGIFとして保存してみましょう。


```python
fig = plt.figure(facecolor="white")
ax = fig.add_subplot()
ase.io.write(
    "curtius_NEB.gif",
    images,
    format="gif",
    ax=ax
)
```


上には最終構造の静止画が表示されていると思いますが、左のファイルビューアから"curtius_NEB.gif"をダブルクリックして開くと直接アニメーションGIFを見ることができます。

正常に計算が実行されていれば、以下のようなアニメーションGIFが見られます。

&lt;img src="assets/tutorial_NEB/curtius_NEB.gif" alt="curtius_NEB" /&gt;

以下のようにしてファイルの中身をここに表示することもできます。


```python
Image(data=open("./curtius_NEB.gif", "rb").read(), format='gif')
```
