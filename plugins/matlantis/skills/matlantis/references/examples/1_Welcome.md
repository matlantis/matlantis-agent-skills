# はじめに

Matlantisへようこそ！ここでは、基礎的な操作を通してMatlantisを利用した原子系シミュレーションの導入のお手伝いをします。

## Jupyter Notebookについて

この操作画面は、Jupyter Notebookを利用して構成されています。既にJupyter Notebookを使用した経験がある方は次の章まで進んでください。

この文章が書かれているエリアがそのまま作業エリアとなります。プログラミング言語はPythonです。

プログラムはセル(四角い枠)に書いていきます。上部の再生ボタン(▶)を押すことで実行されます。
すぐ下のセルにフォーカスを合わせて(クリックして)、上の再生ボタンを押してみましょう。


```python
print("Hello, World!")
```


セルの下に "Hello, World!" と表示されていれば成功です。セルにフォーカスがある状態で`Ctrl+Enter`を入力しても実行されます。`Ctrl+Enter`の代わりに`Shift+Enter`を入力すると、実行した後に次のセルにフォーカスが移動します。(Notebookの末尾の場合、下に新しくセルが作られます。)

セルの追加は上部のセル追加ボタン(+)からも行うことができます。

実行中に宣言された変数や`import`したモジュールなどは、セルの実行が終わった後も保持されます。これによって以前の実行結果を別のセルで使い回すことができます。

次のセルを順番に実行してみましょう。


```python
a = 1
```


```python
a += 1
print(a)
```


実行後に "2" と表示されていれば成功です。(もし上のセルを実行せずに下のセルを実行しようとすると、`a`が未定義である旨のエラーが表示されます。)

また、一度実行したセルを再度実行することもできます。すぐ上の`a += 1`から始まるセルを再度実行すると、aの値が再度加算されて3になっていることが確認できます。

一度まっさらな初期状態に戻したくなったら、上部の再起動ボタン(⟳)を押してください。

このチュートリアルでは、これ以降も解説と並行して順にセルを実行していく形を想定して書いています。

頭にエクスクラメーションマーク`!`をつけるとシェルのコマンドとして実行されます。


```python
!ls -l
```


シェルのコマンドを通して、Python以外のツールの使用やPythonのライブラリの追加などの操作を行うことができます。

なお、時間のかかる処理を実行した場合、処理が終わるまでは上部の"Idle"と表記されている部分が"Busy"という表記に変わります。

Jupyter Notebookの使い方については公式や各所から解説が提供されているため、詳しい使い方についてはそれらを参照されることをおすすめします。

## このファイルについて

このチュートリアルもひとつのJupyter Notebookのファイルとして用意されています。左のファイルビューアに"1_Welcome.ipynb"という名前のファイルができていると思いますが、それがこのチュートリアルの実体です。(説明の部分はPythonではなくMarkdownという文書向けの形式のセルで書かれています。)

通常のNotebookファイルのため好きに改変することができますが、もし後にオリジナルに戻したくなった場合は、左のファイルビューア上部のLauncher(ロケットのアイコン)から再度チュートリアルを選んでください。"Overwrite and Open"を選ぶとチュートリアルの元のファイルで上書きされます。(上書きを避けたい場合、このファイルの名前を変更しておくとよいと思います。)


# シミュレータを使ってみる

## PFPについて

Matlantisで提供する汎用ニューラルネットワークポテンシャルであるPFPを利用するためには、`pfp_api_client`という専用のパッケージが必要です。

初めてPFPを利用する際には、以下のコマンドを実行してください。これによって後述のASEというライブラリも同時にインストールされます。また、PFPをASEから利用する関数もパッケージに同梱されています。


```python
!pip install pfp-api-client
```


一度実行すれば、再度このノートブックを開いた場合や、他のノートブックを利用する場合に再実行する必要はありません。

## ASEについて

PFPを利用できるデフォルトの分子動力学シミュレーションの環境として、[ASE](https://wiki.fysik.dtu.dk/ase/index.html)を提供しています。

ASEに関しても公式にチュートリアル及び各機能のドキュメントが揃っているため、詳しい使用方法についてはそちらを参照することをおすすめします。以下ではASEに沿って説明を進めていきます。

### モジュールを読み込む

上の`!pip install`コマンドを実行していれば、ASEは既に環境にインストール済みです。使うには通常のPythonのモジュールと同じように`import`すれば大丈夫です。

`ModuleNotFoundError: No module named 'ase'`のようなエラーが出た場合は、上部の再起動ボタン(⟳)→`Restart`を押して、この下のセルから再度実行してみてください。


```python
import ase
```


### シミュレーション対象を定義する

ASEでは、シミュレーションを行う対象は`Atoms`クラスのオブジェクトとして管理されます。

ASEでは非常に多くの種類のファイルの入出力に対応しています。詳細は[リファレンス](https://wiki.fysik.dtu.dk/ase/ase/io/io.html)を参照してみてください。

ここでは、簡単のためにビルトインの分子をロードしてみましょう。


```python
from ase.build import molecule
atoms = molecule("CH3CH2OCH3")

print(atoms.get_chemical_symbols())
```


## 可視化してみる

Matlantisではいくつかの可視化ツールをプリインストールしてあります。ここではnglviewerを使ってみましょう。

ASEからnglviewerを起動するのは簡単です。以下のツールを動かすとビューアが起動します。マウスで操作できます。ドラッグで回転、原子をクリックで対象の原子を中央に移動します。


```python
from ase.visualize import view
v = view(atoms, viewer='ngl')
display(v)
```


見た目をカスタマイズすることもできます。以下のコマンドを使うと結合が表示されます。詳細は[リファレンス](https://wiki.fysik.dtu.dk/ase/ase/visualize/visualize.html)を参照してください。


```python
v.view.add_representation("ball+stick")
```


# ポテンシャルを使う

ここまでは標準的なASEの使用方法の解説でした。ここからはいよいよ、PFPを使った分子動力学計算の紹介をしていきます。

## ポテンシャルの読み込みと適用

PFPの本体は、`Estimator`というクラスで定義されています。

また、ASEではポテンシャルを`Calculator`というクラスで定義し、これを`Atoms`に設定することによってポテンシャルの利用が可能になります。`Estimator`をASEで利用するためのクラスが`ASECalculator`です。

スクリプト上でこれらの初期設定を行うのは簡単です。以下のスクリプトを実行するとcalculatorが使える状態として初期化されます。


```python
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator()
calculator = ASECalculator(estimator)
```


あとは、これをAtomsに設定すれば準備完了です。

Atomsにcalculatorを設定すると、エネルギーや力などの各種物性値が計算できるようになります。


```python
atoms.calc = calculator

print("Energy: {}".format(atoms.get_potential_energy()))
print(atoms.get_forces())
```


## 構造最適化計算

calculatorを設定することで、各種MD計算を回すことが可能となります。特殊な計算を回したい場合はユーザ自身で独自の時間発展スクリプトを書きます。

ここではチュートリアルなので、ASEのビルトイン計算を使ってみましょう。基本的な使い方としては、ダイナミクスを定義する一連のクラスがあり、それにAtomsを定義した後に`run()`を呼ぶというスタイルです。

説明よりも先に以下のサンプルを見たほうがわかりやすいかもしれません。ここでは、LBFGS法(準ニュートン法の一種)によって構造最適化を行います。


```python
from ase.optimize import LBFGS

opt = LBFGS(atoms)
opt.run()
```


## 動力学計算

まったく同様にして、分子動力学計算を行うことができます。

すぐ上で構造最適化をしてしまったため、このまま動力学計算をしても原子が動く様子が見えずおもしろくないので、初期速度を加えてみましょう。これらの関数もASEに備わっています。

今回はデモのために100ステップとしたのですぐ終了しますが、一般にMD計算は長時間行うことが多いと思います。また、ASEのMD計算はデフォルトでは計算過程を出力しないため、そのままでは何も出力されない画面のまま待つことになります。ここでは`dyn.attach`を使って計算中の状況を出力する関数を呼ぶようにしています。

さらに、今回はオプションを追加して軌跡も保存してみましょう。Jupyter環境ではファイルの入出力もできます。生成されたファイルは左のファイルビューアに表示されます。


```python
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase.io import Trajectory
from ase import units

# Set the momenta corresponding to T=500K.
MaxwellBoltzmannDistribution(atoms, temperature_K=500.0)
# Sets the center-of-mass momentum to zero.
Stationary(atoms)
# Run MD using the VelocityVerlet algorithm
dyn = VelocityVerlet(atoms, 1.0 * units.fs, trajectory="dyn.traj")

def print_dyn():
    print(f"Dyn  step: {dyn.get_number_of_steps(): &gt;3}, energy: {atoms.get_total_energy():.3f}")

dyn.attach(print_dyn, interval=10)
dyn.run(100)
```


軌跡はAtomsの配列のような情報を持っています。ASEではTrajectoryクラスによって管理されます。

Atomsと同様に可視化することができます。今回は軌跡なので、再生することができます。以下の可視化画面の右下にマウスカーソルを合わせると、再生ボタンが表示されます。クリックするとMD計算の結果を動画で表示します。


```python
traj = Trajectory("dyn.traj")
v = view(traj, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


## 次に進む

ここまで読んでプログラムを問題なく実行できた方は、PFPを使った原子シミュレーションを利用できる環境は整っています。もしJupyter環境およびASEを利用した分子動力学シミュレーションについて既に知識のある方は、さっそく新しいNotebookを開いてシミュレーションを始めることもできます。

もしもっと実例について知りたい方は、Launcher(このファイルを開いたタブ、左上のロケットマークから開けます)から、引き続き"Elasticity"および"NEB"などの具体例を題材としたチュートリアルを触ることができます。

また、Matlantisでは直接MDシミュレーションのプログラムを1行ずつ書く代わりに、特定の物性値を計算するなどの定形計算を短いコマンドでできるよう設計したライブラリ`matlantis-features`を提供しています。`matlantis-features`の使い方や使用例について知りたい人向けに、同じくLauncherからMatlantis-featuresのExampleを提供しています。Exampleの中ではじめて読む例としては"Matlantis-features : Elasticity"がおすすめです。
