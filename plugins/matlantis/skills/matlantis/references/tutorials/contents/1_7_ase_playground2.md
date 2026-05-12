# Atomsの操作2

前節に引き続き、Atomsを操作していく実例を通してその扱いにより深く慣れていきましょう。


```python
from ase import Atoms
from ase.build import molecule, bulk

from pfcc_extras.visualize.view import view_ngl
```


    


## Super cellの作成

前節でもすでに説明したように、周期系の構造は`repeat`または、掛け算 "*" の記法を用いる用いることでCellサイズを大きくすることができます。


```python
atoms = bulk("Si") * (2, 3, 4)
atoms222 = atoms * (2, 2, 2)
view_ngl(atoms222, w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), value='All'), D…



## Cell の取り扱い

周期構造を持つ `atoms`は `cell` propertyを持ちますが、これは `Cell` クラスという特別なクラスで定義されており、いくつか追加の関数を持っています。


```python
# atoms = bulk("Fe", cubic=True)
atoms = bulk("Fe") * (2, 2, 2)
view_ngl(atoms, w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Fe'), value='All'), D…




```python
atoms.cell
```




    Cell([[-2.87, 2.87, 2.87], [2.87, -2.87, 2.87], [2.87, 2.87, -2.87]])



&lt;img src=https://upload.wikimedia.org/wikipedia/commons/5/5e/UnitCell.png width="400"&gt;

&lt;cite&gt;[分率座標](https://ja.wikipedia.org/wiki/%E5%88%86%E7%8E%87%E5%BA%A7%E6%A8%99)より&lt;/cite&gt;

結晶のCellはa軸, b軸, c軸のベクトルをそれぞれ示した3x3 の配列でも表せますが、分率座標を用いると、軸の長さ a, b, cと、軸の角度 α, β, γ で表すこともできます。&lt;br/&gt;
以下の関数はそれぞれ次の値を返します。

 - `cellpar()`: (a, b, c, α, β, γ)
 - `lengths()`: (a, b, c)
 - `angles()`: (α, β, γ)


```python
print("cellpar(): ", atoms.cell.cellpar())
print("lengths(): ", atoms.cell.lengths())
print("angles() : ", atoms.cell.angles())
```

    cellpar():  [  4.97098582   4.97098582   4.97098582 109.47122063 109.47122063
     109.47122063]
    lengths():  [4.97098582 4.97098582 4.97098582]
    angles() :  [109.47122063 109.47122063 109.47122063]


3x3 配列にアクセスするには `array`で numpy array が取得できます。


```python
atoms.cell.array
```




    array([[-2.87,  2.87,  2.87],
           [ 2.87, -2.87,  2.87],
           [ 2.87,  2.87, -2.87]])



`cell` はnumpy の関数を使用することもできます。例えば以下はnumpyの`diagonal`関数で対角成分のみを取得しています。


```python
atoms.cell.diagonal()
```




    array([-2.87, -2.87, -2.87])



`get_bravais_lattice` でブラベ格子や、`bandpath` でバンドパスを得ることもできます。


```python
atoms.cell.get_bravais_lattice()
```




    BCC(a=5.7400000000000002132)




```python
atoms.cell.bandpath()
```




    BandPath(path='GHNGPH,PN', cell=[3x3], special_points={GHNP}, kpts=[25x3])



`standard_form()` methodはCellの形を下三角行列に標準化した場合のCellの形 `rcell` と、そのようにするための回転行列`Q` を計算します。


```python
rcell, Q = atoms.cell.standard_form()
rcell
```




    Cell([[4.970985817722678, 0.0, 0.0], [-1.6569952725742265, 4.686690374525147, 0.0], [-1.6569952725742256, -2.3433451872625737, 4.058792924010784]])



## 元素置換

ある結晶に対して、少量の添加物を加えた場合の構造を作りたい場合、ベースとなる結晶を作成した後に原子を置換する形でモデリングすることができます。

以下の例は、Fe の金属中にVを添加する例です。

※表面構造を生成する関数 `bcc111` では、 `vacuum` を指定することでz軸方向の真空層の長さを指定することができます。&lt;br/&gt;
また、`periodic=True`とすることで、z軸方向に対しても周期境界条件を適用します。&lt;br/&gt;
これを指定しないとCellサイズが定まらず、nglviewerの可視化でエラーが起こるので注意してください。


```python
from ase.build import bcc111

atoms = bcc111("Fe", size=(4, 4, 6), periodic=True, vacuum=10)
view_ngl(atoms, w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Fe'), value='All'), D…



はじめはすべてFe (原子番号26)となっています。


```python
numbers = atoms.get_atomic_numbers()
numbers
```




    array([26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
           26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
           26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
           26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
           26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
           26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26])



今回はランダムに5つの原子をV (原子番号23) に変えてみます。


```python
import numpy as np

replace_index = np.random.choice(np.arange(len(numbers)), size=5, replace=False)
print(f"Replace {replace_index} atom")
numbers[replace_index] = 23
atoms.set_atomic_numbers(numbers)
```

    Replace [83 59 93  7 79] atom



```python
view_ngl(atoms, w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Fe', 'V'), value='All…



ちなみに、それぞれの元素の原子番号を知りたい場合は、 `ase.data.atomic_numbers` を参照することができます。


```python
from ase.data import atomic_numbers
print("atomic_numbers: ", atomic_numbers)

print(f"V's atomic number is", atomic_numbers['V'])
```

    atomic_numbers:  {'X': 0, 'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}
    V's atomic number is 23


逆に、元素番号から元素記号・名前を知りたい場合は、 `chemical_symbols`, `atomic_names` を使用できます。


```python
from ase.data import atomic_names, chemical_symbols
print("chemical_symbols: ", chemical_symbols)
print("atomic_names    : ", atomic_names)

print(f"atomic number 23 is", chemical_symbols[23], atomic_names[23])
```

    chemical_symbols:  ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    atomic_names    :  ['', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron', 'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium', 'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur', 'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium', 'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron', 'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium', 'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium', 'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium', 'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium', 'Promethium', 'Samarium', 'Europium', 'Gadolinium', 'Terbium', 'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium', 'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury', 'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine', 'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium', 'Protactinium', 'Uranium', 'Neptunium', 'Plutonium', 'Americium', 'Curium', 'Berkelium', 'Californium', 'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium', 'Lawrencium', 'Rutherfordium', 'Dubnium', 'Seaborgium', 'Bohrium', 'Hassium', 'Meitnerium', 'Darmastadtium', 'Roentgenium', 'Copernicium', 'Nihonium', 'Flerovium', 'Moscovium', 'Livermorium', 'Tennessine', 'Oganesson']
    atomic number 23 is V Vanadium


## 空孔の作成

結晶構造に欠陥がある場合の挙動を調べることは産業上重要です。

ASEでは `pop` という関数を使うことにより、特定のindexの原子の削除ができます。
(または、`del`関数も使用可能です。)

Si結晶から、以下は13番目のSi原子を除いた欠陥構造を作る例です。


```python
atoms = bulk("Si") * (2, 2, 2)
# `del atoms[13]` also works
atoms.pop(13)
view_ngl(atoms, representations=["ball+stick"], w=400, h=300)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), value='All'), D…



以上、ASEに実装されているいろいろな機能をみて、原子をどのように取り扱うかを学びました。

ここから先は実践で、計算手法・事例とともに学習をしていきましょう。
