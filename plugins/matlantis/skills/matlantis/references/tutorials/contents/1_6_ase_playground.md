# Atomsの操作

ここまででASEの基礎を学び Atoms や Calculatorの扱い方を学びました。&lt;br/&gt;
本節では、Atomsを操作していく実例を通してその扱いにより深く慣れていきましょう。


```python
from ase import Atoms
from ase.build import molecule, bulk
```

## コピー

`atoms`を作成したあと、そのコピーを作成するには `copy()` methodを使うことができます。&lt;br/&gt;
原子の元素種・座標値や、セルなどがコピーされます。

以下では`atoms2` が`atoms`と同じ座標値を持ったH2Oのコピーになっています。


```python
atoms = molecule("H2O")
atoms2 = atoms.copy()

print("atoms :", atoms)
print("pos", atoms.positions)
print()
print("atoms2:", atoms2)
print("pos", atoms2.positions)
```

    atoms : Atoms(symbols='OH2', pbc=False)
    pos [[ 0.        0.        0.119262]
     [ 0.        0.763239 -0.477047]
     [ 0.       -0.763239 -0.477047]]
    
    atoms2: Atoms(symbols='OH2', pbc=False)
    pos [[ 0.        0.        0.119262]
     [ 0.        0.763239 -0.477047]
     [ 0.       -0.763239 -0.477047]]


## propertyの書き換え

Atomsは様々なget, set 関数を持っています。これらを通じてAtomsの持つ属性値を変更することが可能です。

 - https://wiki.fysik.dtu.dk/ase/ase/atoms.html#working-with-the-array-methods-of-atoms-objects
 
例えば前節では `set_momenta` 関数を使用して運動量を設定する例がありました。

## 座標値の変更

座標値の変更の場合、`atoms.positions` をそのまま書き換えることができます。


```python
atoms = molecule("CH3CHO")
atoms.positions
```




    array([[ 1.218055,  0.36124 ,  0.      ],
           [ 0.      ,  0.464133,  0.      ],
           [-0.477241,  1.465295,  0.      ],
           [-0.948102, -0.700138,  0.      ],
           [-0.385946, -1.634236,  0.      ],
           [-1.596321, -0.652475,  0.880946],
           [-1.596321, -0.652475, -0.880946]])




```python
from pfcc_extras.visualize.view import view_ngl

view_ngl(atoms, representations=["ball+stick"], w=400, h=300)
```


    





    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'H', 'C', 'O'), value=…



0番目のO原子を `[1.218055, 0.36124, 0.]` から、 `[2.0, 0, 0]`　に動かしてみます。

`atoms.positions` をそのまま書き換えることができます。


```python
# First axis 0 specifies atom index = O
# Second axis with the value 0,1,2 corrensponds to x, y, z axis respectively.
atoms.positions[0] = [2.0, 0, 0]
```


```python
view_ngl(atoms, representations=["ball+stick"], w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'H', 'C', 'O'), value=…



O原子が動いたことが確認できました。

## 平行移動、回転

atoms全体の平行移動や回転は、`translate`, `rotate` 関数を使うことができます。

`translate`関数は (3,) のxyzベクトルを指定するとすべての分子を同じだけ平行移動し、 (n, 3) のベクトルを指定するとn個の原子を別々の量平行移動します。&lt;br/&gt;
以下の例では、すべての原子を `[1.0, 0, 0]` 平行移動しています。 &lt;br/&gt;
実際に座標値を見ると x軸方向だけ平行移動されているのがわかります。


```python
atoms.translate([1.0, 0, 0])
atoms.positions
```




    array([[ 3.      ,  0.      ,  0.      ],
           [ 1.      ,  0.464133,  0.      ],
           [ 0.522759,  1.465295,  0.      ],
           [ 0.051898, -0.700138,  0.      ],
           [ 0.614054, -1.634236,  0.      ],
           [-0.596321, -0.652475,  0.880946],
           [-0.596321, -0.652475, -0.880946]])



`rotate`関数は `atoms`を回転させるための関数です。

以下の例は z軸方向を回転軸として 90°回転させています。&lt;br/&gt;
回転軸は、`v=[0,0,1]`のようにベクトルで指定することもできます。


```python
atoms.rotate(90, v="z")
atoms.positions
```




    array([[ 1.8369702e-16,  3.0000000e+00,  0.0000000e+00],
           [-4.6413300e-01,  1.0000000e+00,  0.0000000e+00],
           [-1.4652950e+00,  5.2275900e-01,  0.0000000e+00],
           [ 7.0013800e-01,  5.1898000e-02,  0.0000000e+00],
           [ 1.6342360e+00,  6.1405400e-01,  0.0000000e+00],
           [ 6.5247500e-01, -5.9632100e-01,  8.8094600e-01],
           [ 6.5247500e-01, -5.9632100e-01, -8.8094600e-01]])




```python
view_ngl(atoms, representations=["ball+stick"], w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'H', 'C', 'O'), value=…



周期構造を持つ系を回転させた場合、以下のように、Cellはそのままで原子の座標のみが回転されるため、原子はCellからはみ出してしまいます。


```python
atoms = bulk("Fe") * (2, 3, 4)

atoms.rotate(90, v=[0, 0, 1])
view_ngl(atoms, representations=["ball+stick"], w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Fe'), value='All'), D…



Cellの座標も同時に回転させるには、 `rotate_cell=True` とします。


```python
atoms = bulk("Fe") * (2, 3, 4)

atoms.rotate(90, v=[0, 0, 1], rotate_cell=True)
view_ngl(atoms, representations=["ball+stick"], w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Fe'), value='All'), D…



## ランダム移動

`rattle` 関数を使うと、それぞれの原子をランダム移動させて少し乱雑な構造を作成することができます。

ここでは、まずSi 結晶を作成します。


```python
atoms = bulk("Si") * (2, 3, 4)
view_ngl(atoms, representations=["ball+stick"], w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), value='All'), D…



この系を `rattle` で少し構造を乱してみます。変位の方向は`seed`の指定を変えることで、変えることができます。


```python
atoms.rattle(stdev=0.2, seed=1)
view_ngl(atoms, w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), value='All'), D…



`stdev`の値を変えることでより大きく乱す事もできます。


```python
atoms.rattle(stdev=0.5)
view_ngl(atoms, w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), value='All'), D…



[Note (Advanced)] ASEの`rattle` 関数の引数である`seed` が指定されていないときは **`seed=42` が使われるという仕様になっており、結果`rattle`の結果はDeterministicになる**ので注意してください。

以下の例では、2回 `rattle`を行った結果が等しくなっていることを確認しています。


```python
import numpy as np


atoms = bulk("Si") * (2, 3, 4)
atoms2 = atoms.copy()

atoms.rattle(stdev=0.2)
atoms2.rattle(stdev=0.2)

print("atoms &amp; atoms2 positions are same? --&gt; ", np.allclose(atoms.positions, atoms2.positions))
```

    atoms &amp; atoms2 positions are same? --&gt;  True


もし、毎回違う結果をランダムに得たい場合は `rng=np.random.RandomState()` を指定してください。


```python
atoms = bulk("Si") * (2, 3, 4)
atoms2 = atoms.copy()

atoms.rattle(stdev=0.2, rng=np.random.RandomState())
atoms2.rattle(stdev=0.2, rng=np.random.RandomState())

print("atoms &amp; atoms2 positions are same? --&gt; ", np.allclose(atoms.positions, atoms2.positions))
```

    atoms &amp; atoms2 positions are same? --&gt;  False


## wrap

上記のSi結晶のように、Cellの外側に原子が飛び出してしまった場合、 `wrap` 関数を使うことで再び対応する周期構造の内側に戻す事ができます。


```python
atoms.wrap()
view_ngl(atoms, w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), value='All'), D…


