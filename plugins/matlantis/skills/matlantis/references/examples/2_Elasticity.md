# 準備

PFPを利用するためには、`pfp_api_client`という専用のパッケージが必要です。

初めてPFPを利用する際には、以下のコマンドを実行してください。

一度実行すれば、再度このノートブックを開いた場合や、他のノートブックを利用する場合に再実行する必要はありません。


```python
!pip install pfp-api-client
```


# 結晶構造と分子静力学

ここでは、材料の弾性的性質の解析を題材として、構造最適化計算の実例を紹介します。

この例題では結晶シリコンを対象とします。

まずは関連するモジュールをインポートしてしまいます。


```python
import numpy as np

import ase
from ase import Atoms, units
from ase.constraints import ExpCellFilter
from ase.visualize import view
from ase.optimize import BFGS

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator()
calculator = ASECalculator(estimator)
```


結晶構造の作成方法は複数あります。ASEにある結晶構造作成ツールを使うこともできますし、知られている構造であれば外部のデータベースを参照することもできます。

今回はチュートリアルということなので、直接指定して`Atoms`を作ってみることにします。


```python
lattice_constant = 5.43

atoms = ase.Atoms(
    symbols=["Si"] * 8,
    positions = lattice_constant * np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]),
    cell=[lattice_constant] * 3,
    pbc=True,
)

v = view(atoms, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


格子定数は文献値をあたりました。密度から計算することもできます。positionsではdiamond構造の単位格子を指定しています。

`cell`および`pbc`が指定されていることに注意してください。これにより系に周期境界条件が適用されるため、表面を持たない無限に続く構造のモデルになります。

解析の第一歩は、構造最適化計算によって安定構造を探すことです。ここでは原子の位置に加えてセルサイズも最適化します。これをしないと初期設定の格子定数に引っ張られて不自然な圧力がかかった構造になってしまうためです。

ASEでは、`Atoms`のラッパーとしてこのような用途のために自由度を増やしたクラス(`Filter`)が用意されています。今回はこれを使います。


```python
atoms.calc = calculator
filter = ExpCellFilter(atoms)
opt = BFGS(filter)
opt.run(fmax=0.01)
```


最適化計算によって格子の長さが少し変わっています。また、圧力が緩和されたため系の圧力(正確には応力テンソル)がほぼ0になっているはずです。確認してみましょう。


```python
print(atoms.get_cell().lengths().tolist())
print(atoms.get_stress().tolist())
```


さて、次は結晶の弾性定数を求めてみましょう。弾性定数は、セルのひずみに対する応力の変化の度合いとして定義されています。

x方向に少しひずませてみましょう。+1%ひずませることを考えると、ひずみテンソルと変形勾配テンソルは以下のように定義できます。これを使ってAtomsオブジェクトを変形させることができます。その後、Cell sizeを固定したまま、原子座標のみ最適化を行うことで、セルがひずんだ時の応力を求めます。(値を変更して再計算するときは、上の構造最適化計算から実行し直してください。)


```python
x_strain = 0.01
strain_matrix = np.array([
    [x_strain, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
])
deform_matrix = strain_matrix + np.eye(3)
deformed_cell = atoms.get_cell().dot(deform_matrix)
atoms.set_cell(deformed_cell, scale_atoms=True)

opt = BFGS(atoms)
opt.run(fmax=0.01)
```


このときの応力テンソルをみてみましょう。


```python
print(atoms.get_stress().tolist())
```


ひずみは0.01, ASEの単位系がAngstrom及びeVであることに注意すると、応力をひずみで割った弾性定数C11およびC12は以下のように求まります。


```python
c11 = atoms.get_stress()[0] / x_strain
c12 = atoms.get_stress()[1] / x_strain
print("C11: {:.3f} GPa".format(c11 / units.GPa))
print("C12: {:.3f} GPa".format(c12 / units.GPa))
```


なお実験値はC11が約167 GPa, C12が約65 GPaとされています。
