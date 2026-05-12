# Materials Projectの利用例

## アクセス方法

Materials Projectのページにアクセスします。

https://materialsproject.org/

Materials Projectの利用にはユーザー登録が必要です。`Sign In or Register`からユーザー登録が行えます。

&lt;img src="assets/tutorial_materials_project/registration.png" alt="registration" /&gt;

登録はメールアドレスやGoogleアカウントなど、いくつかの選択肢が用意されています。例えばメールアドレスによる登録では、記入したアドレスにログイン用のURLが書かれたメールが送信されます。URLへアクセスするとログインの手続きが完了します。

## 検索例

ここでは例としてNa2FePO4F(文献: https://doi.org/10.1038/nmat2007 )の構造を調べたいケースを考えます。

構造の検索の方法はいくつかあって、"Explore Materials `by Elements`"のプルダウンメニューから様々な調べ方が選べます。ここではデフォルトの`by Elements`で進めてみます。

直接キーボードで入力してもよいですが、下の周期表から`Na`, `Fe`, `P`, `O`, `F`をクリックすると検索欄に自動的に入力されます。このうえで`Search`をクリックすると下に該当の構造が出てきます。

https://materialsproject.org/#search/materials/{%22nelements%22%3A5%2C%22elements%22%3A%22Na-Fe-P-O-F%22}

Materials Projectは登録されている構造が多いので、同じ元素の組み合わせでも多数の構造が出てくることがあります。(例えばSi-Oなどで調べると数百件ヒットします。)もし調べたい構造が決まっていて、物性値について情報がある場合には、右の"External Provenance"から検索条件を絞り込むこともできます。

今回は該当のNa2FePO4Fの構造はmp-1194940としてリストに乗っています。今回のように組成までわかっている場合には、検索欄に直接`Na2FePO4F`と書いてしまえば絞り込むことができます。("Explore Materials `by Formula`"に自動的に切り替わります。)

https://materialsproject.org/#search/materials/{%22reduced_cell_formula%22%3A%22Na2FePO4F%22}

### `by Elements` での検索
&lt;img src="assets/tutorial_materials_project/search_by_elements.gif" alt="search_by_elements" /&gt;

### `by Formula` での検索
&lt;img src="assets/tutorial_materials_project/search_by_formula.gif" alt="search_by_formula" /&gt;

## 構造のダウンロード

構造をクリックすると構造に対する詳細ページに飛びます。

https://materialsproject.org/materials/mp-1194940/

可視化された構造および、エネルギーやバンドギャップなどDFT計算から得られる多種の計算結果を見ることができます。(構造はマウスで回転させることができます。)

`CIF`と書かれたダウンロードボタンから構造のCIFファイルをダウンロードすることができます。結晶の対称性やセルの扱いに関してCIFファイルの選択肢が出ますが、ここでは計算に使われた構造に相当する`Computed`を選択してダウンロードしてみましょう。

ダウンロードしたファイル名は手元では`Na2FePO4F_mp-1194940_computed.cif`となっていました。

&lt;img src="assets/tutorial_materials_project/download.gif" alt="download" /&gt;




## Notebookでの利用例

ダウンロードしたファイルをNotebookにアップロードすることができます。Notebookのページ左上の`Upload files`から、あるいは単純にファイルを左のフォルダ欄にドラッグ&amp;ドロップすることでファイルがアップロードできます。

アップロードしたCIFファイルはpythonのASEライブラリで扱うことができます。

&lt;img src="assets/tutorial_materials_project/upload.gif" alt="upload" /&gt;



```python
# If you have installed 'pfp-api-client', you can skip this.
!pip install pfp-api-client
```


```python
import ase.io
from ase.visualize import view

atoms = ase.io.read("Na2FePO4F_mp-1194940_computed.cif")

v = view(atoms, viewer='ngl')
display(v)
```


ここから先は任意のMDシミュレーションに利用することができます。以下は有限温度でのシミュレーションの例です。


```python
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase.io import Trajectory
from ase import units

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator()
calculator = ASECalculator(estimator)
atoms *= (2, 1, 1)
atoms.calc = calculator

MaxwellBoltzmannDistribution(atoms, temperature_K=1000.0)
Stationary(atoms)
dyn = VelocityVerlet(atoms, 1.0 * units.fs, trajectory="dyn.traj")

def print_dyn():
    print(f"Dyn  step: {dyn.get_number_of_steps(): &gt;3}, energy: {atoms.get_total_energy():.3f}")

dyn.attach(print_dyn, interval=10)
dyn.run(100)

traj = Trajectory("dyn.traj")
v = view(traj, viewer='ngl')
display(v)
```
