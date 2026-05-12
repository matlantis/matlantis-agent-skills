# ASE Basic

Pythonで原子シミュレーションを進めていく際に便利なOSSライブラリとして、Atomic Simulation Environment: ASEがあります。&lt;br/&gt;
本節ではASEの基本的な使い方を紹介していきます。

 - Documentation: https://wiki.fysik.dtu.dk/ase/
 - gitlab: https://gitlab.com/ase/ase

## Atomsについて

ASEでは原子が複数個集まって作られる系を `Atoms` クラスで表現します。

 - https://wiki.fysik.dtu.dk/ase/ase/atoms.html

`Atoms`クラスは前節で説明した、原子シミュレーションに必要な構造を表現するための以下のような`attribute`(変数)を保持しています。

 - 各原子の元素種類
 - 座標値
 - 速度(運動量)
 - セル
 - 周期境界条件など

### Atomsの作成方法： 原子種・座標を直接指定

`Atoms`を作るための一番原始的な方法として、原子種とその座標値を直接指定する方法があります。&lt;br/&gt;
以下は、1つめのHがxyz座標値 `[0, 0, 0]`、2つめのHがxyz座標値 `[1.0, 0, 0]` に存在する水素分子H2を作る例です。


```python
from ase import Atoms

atoms = Atoms("H2", [[0, 0, 0], [1.0, 0, 0]])
```

原子の可視化を行うこともできます。
ここでは、ASEのview 関数を用い、`nglviewer`というライブラリを用いて可視化を行っています。

`nglviewer` を用いると、３次元の原子構造を表示しながら、マウスでインタラクティブに操作することもできます。

可視化についての詳細は[Appendix_1_visualization.ipynb](./Appendix_1_visualization.ipynb)をご参照ください。


```python
from ase.visualize import view
from pfcc_extras.visualize.view import view_ngl

#view(atoms, viewer="ngl")
view_ngl(atoms, representations=["ball+stick"], w=400, h=300)
```


    





    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'H'), value='All'), Dr…



以下のように、ASEを用いてpng静止画を描画することも可能です。


```python
from ase.io import write
from IPython.display import Image

write("output/H2.png", atoms, rotation="0x,0y,0z", scale=100)
Image(url='output/H2.png', width=150)
```




&lt;img src="output/H2.png" width="150"/&gt;



以下は上記のプログラムを１行でできるようにした便利関数です。今後はこれを用います。


```python
from pfcc_extras.visualize.ase import view_ase_atoms

view_ase_atoms(atoms, figsize=(4, 4), title="H2 atoms", scale=100)
```


    
![png](output_9_0.png)
    


元素記号 `symbols`の代わりに、原子番号を`numbers` として指定することも可能です。


```python
co_atoms = Atoms(numbers=[6, 8], positions=[[0, 0, 0], [1.0, 0, 0]])
view_ngl(co_atoms, representations=["ball+stick"], w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'C', 'O'), value='All'…




```python
view_ase_atoms(co_atoms, figsize=(4, 4), title="CO atoms", scale=100)
```


    
![png](output_12_0.png)
    


周期境界条件を含む系を定義したい場合は、`cell` に周期情報を指定し、`pbc` にa軸b軸c軸それぞれに対して周期境界条件を適用するかどうかのON/OFFを指定することができます。


```python
from ase import Atoms


na2_atoms = Atoms(
    symbols="Na2",
    positions=[[0, 0, 0], [2.115, 2.115, 2.115]],
    cell=[4.23, 4.23, 4.23],
    pbc=[True, True, True]
)
```


```python
view_ngl(na2_atoms, w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Na'), value='All'), D…




```python
view_ase_atoms(na2_atoms, rotation="10x,10y,0z", figsize=(3, 3), title="Na2 atoms", scale=30)
```


    
![png](output_16_0.png)
    


### Atomsの持つattribute, method

&lt;figure style="width:200px"&gt;
&lt;img src="../assets/ch1/atoms.png"/&gt;
&lt;/figure&gt;

`atoms` はその系の情報を保持しており、Attribute や `get_XXX` 関数を通して参照することができます。

Attributeとは Atomsクラスの持つ変数のことで、 `atoms.xxx` と直接参照することが出来るものもあります。&lt;br/&gt;
関数を通してアクセスできる情報もあります。

以下は主要なAttributeや関数の返り値(numpy array としてのshapeがどの様になっているか)をまとめたものです。

 - Attribute
   - `symbols`: 元素種とその数をまとめたものを返す
   - `numbers`: 原子番号のarray。shapeは(N,)
   - `positions`: 各原子の座標値を示す。shapeは(N, 3)。 3はxyz座標を表す。
   - `cell`: ASEのCellをあつかうCellクラスで表現される。
     - Cellクラスに関して詳しくは後述する。[三斜晶系](https://ja.wikipedia.org/wiki/%E4%B8%89%E6%96%9C%E6%99%B6%E7%B3%BB)など一般的にはa, b, c軸のベクトルをそれぞれ表す(3, 3) の行列で表されるが、[直方晶](https://ja.wikipedia.org/wiki/%E7%9B%B4%E6%96%B9%E6%99%B6%E7%B3%BB)の場合はa, b, c軸の長さのみが(3,) で表現される。
   - `pbc`: 各方向に周期境界条件(Periodic Boundary Condition) があるかどうかを示す。
 - 関数
   - `get_masses()`: 各原子の質量を返す。
   - `get_momenta()`: 各原子に設定されている運動量を返す。
   - `len(atoms)`: Atoms全体の原子数`N`を返す。


```python
print(f"symbols  : {atoms.symbols}")
print(f"positions: {atoms.positions}")
print(f"cell     : {atoms.cell}")
print(f"pbc      : {atoms.pbc}")

# "Atomic numbers" and "Mass" can be automatically calculated from "Symbol"
print(f"numbers  : {atoms.numbers}")
print(f"massess  : {atoms.get_masses()}")
print(f"momenta  : {atoms.get_momenta()}")
print(f"Numbers of atoms  : {len(atoms)}")
```

    symbols  : H2
    positions: [[0. 0. 0.]
     [1. 0. 0.]]
    cell     : Cell([0.0, 0.0, 0.0])
    pbc      : [False False False]
    numbers  : [1 1]
    massess  : [1.008 1.008]
    momenta  : [[0. 0. 0.]
     [0. 0. 0.]]
    Numbers of atoms  : 2


### 単位系

ASEでは、エネルギーの単位はeV (約 $1.602 \times 10^{-19}$ J)、座標系の単位はÅ ($1 \times 10^{-10}$ m)で扱われる事が多いです。&lt;br/&gt;
力や応力の単位はこれらの複合単位で、力についてはeV/Å、応力についてはeV/Å2です。chargeの単位は電荷素量です。

 - eV: [電子ボルト](https://ja.wikipedia.org/wiki/%E9%9B%BB%E5%AD%90%E3%83%9C%E3%83%AB%E3%83%88)
 - Å: [オングストローム](https://ja.wikipedia.org/wiki/%E3%82%AA%E3%83%B3%E3%82%B0%E3%82%B9%E3%83%88%E3%83%AD%E3%83%BC%E3%83%A0)
 - e: [電気素量](https://ja.wikipedia.org/wiki/%E9%9B%BB%E6%B0%97%E7%B4%A0%E9%87%8F)

[Tips] Jupyter 上で各クラス・関数のより詳細な情報を知りたい場合、　? を後ろにつけることでドキュメントの表示ができます。


```python
Atoms?
```


    [0;31mInit signature:[0m
    [0mAtoms[0m[0;34m([0m[0;34m[0m
    [0;34m[0m    [0msymbols[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mpositions[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mnumbers[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mtags[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmomenta[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmasses[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mmagmoms[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mcharges[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mscaled_positions[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mcell[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mpbc[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mcelldisp[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mconstraint[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mcalculator[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0minfo[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m    [0mvelocities[0m[0;34m=[0m[0;32mNone[0m[0;34m,[0m[0;34m[0m
    [0;34m[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0;31mDocstring:[0m     
    Atoms object.
    
    The Atoms object can represent an isolated molecule, or a
    periodically repeated structure.  It has a unit cell and
    there may be periodic boundary conditions along any of the three
    unit cell axes.
    Information about the atoms (atomic numbers and position) is
    stored in ndarrays.  Optionally, there can be information about
    tags, momenta, masses, magnetic moments and charges.
    
    In order to calculate energies, forces and stresses, a calculator
    object has to attached to the atoms object.
    
    Parameters
    ----------
    symbols : str | list[str] | list[Atom]
        Chemical formula, a list of chemical symbols, or list of
        :class:`~ase.Atom` objects (mutually exclusive with ``numbers``).
    
        - ``'H2O'``
        - ``'COPt12'``
        - ``['H', 'H', 'O']``
        - ``[Atom('Ne', (x, y, z)), ...]``
    
    positions : list[tuple[float, float, float]]
        Atomic positions in Cartesian coordinates
        (mutually exclusive with ``scaled_positions``).
        Anything that can be converted to an ndarray of shape (n, 3) works:
        [(x0, y0, z0), (x1, y1, z1), ...].
    scaled_positions : list[tuple[float, float, float]]
        Atomic positions in units of the unit cell
        (mutually exclusive with ``positions``).
    numbers : list[int]
        Atomic numbers (mutually exclusive with ``symbols``).
    tags : list[int]
        Special purpose tags.
    momenta : list[tuple[float, float, float]]
        Momenta for all atoms in Cartesian coordinates
        (mutually exclusive with ``velocities``).
    velocities : list[tuple[float, float, float]]
        Velocities for all atoms in Cartesian coordinates
        (mutually exclusive with ``momenta``).
    masses : list[float]
        Atomic masses in atomic units.
    magmoms : list[float] | list[tuple[float, float, float]]
        Magnetic moments.  Can be either a single value for each atom
        for collinear calculations or three numbers for each atom for
        non-collinear calculations.
    charges : list[float]
        Initial atomic charges.
    cell : 3x3 matrix or length 3 or 6 vector, default: (0, 0, 0)
        Unit cell vectors.  Can also be given as just three
        numbers for orthorhombic cells, or 6 numbers, where
        first three are lengths of unit cell vectors, and the
        other three are angles between them (in degrees), in following order:
        [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)].
        First vector will lie in x-direction, second in xy-plane,
        and the third one in z-positive subspace.
    celldisp : tuple[float, float, float], default: (0, 0, 0)
        Unit cell displacement vector. To visualize a displaced cell
        around the center of mass of a Systems of atoms.
    pbc : bool | tuple[bool, bool, bool], default: False
        Periodic boundary conditions flags.
    
        - ``True``
        - ``False``
        - ``0``
        - ``1``
        - ``(1, 1, 0)``
        - ``(True, False, False)``
    
    constraint : constraint object(s)
        One or more ASE constraints applied during structure optimization.
    calculator : calculator object
        ASE calculator to obtain energies and atomic forces.
    info : dict | None, default: None
        Dictionary with additional information about the system.
        The following keys may be used by ASE:
    
        - spacegroup: :class:`~ase.spacegroup.Spacegroup` instance
        - unit_cell: 'conventional' | 'primitive' | int | 3 ints
        - adsorbate_info: Information about special adsorption sites
    
        Items in the info attribute survives copy and slicing and can
        be stored in and retrieved from trajectory files given that the
        key is a string, the value is JSON-compatible and, if the value is a
        user-defined object, its base class is importable.  One should
        not make any assumptions about the existence of keys.
    
    Examples
    --------
    &gt;&gt;&gt; from ase import Atom
    
    N2 molecule (These three are equivalent):
    
    &gt;&gt;&gt; d = 1.104  # N2 bondlength
    &gt;&gt;&gt; a = Atoms('N2', [(0, 0, 0), (0, 0, d)])
    &gt;&gt;&gt; a = Atoms(numbers=[7, 7], positions=[(0, 0, 0), (0, 0, d)])
    &gt;&gt;&gt; a = Atoms([Atom('N', (0, 0, 0)), Atom('N', (0, 0, d))])
    
    FCC gold:
    
    &gt;&gt;&gt; a = 4.05  # Gold lattice constant
    &gt;&gt;&gt; b = a / 2
    &gt;&gt;&gt; fcc = Atoms('Au',
    ...             cell=[(0, b, b), (b, 0, b), (b, b, 0)],
    ...             pbc=True)
    
    Hydrogen wire:
    
    &gt;&gt;&gt; d = 0.9  # H-H distance
    &gt;&gt;&gt; h = Atoms('H', positions=[(0, 0, 0)],
    ...           cell=(d, 0, 0),
    ...           pbc=(1, 0, 0))
    [0;31mFile:[0m           ~/.py311/lib/python3.11/site-packages/ase/atoms.py
    [0;31mType:[0m           type
    [0;31mSubclasses:[0m     MSONAtoms, Cluster, Lattice

