Copyright Preferred Computational Chemistry, Inc. and Preferred Networks, Inc. as contributors to Matlantis contrib project

# ランダム構造をMatlantis上で生成するexample
※原子数が増えるほど実行時間がかかります。


```python
!pip install ase pymatgen

# # 初回使用時のみ、ライブラリのインストールをお願いします。
```

## 1. ASEを使ったランダム構造生成


```python
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from ase.ga.utilities import closest_distances_generator, CellBounds
from ase.data import atomic_numbers, covalent_radii
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.visualize import view
```

### 1-1. 結晶


```python
blocks = ['Ti'] * 4 + ['O'] * 8 # 作りたい結晶の組成
box_volume = 12 * 12 # 10~12 * 原子数くらいが経験的に上手くいきます

blmin = closest_distances_generator(atom_numbers=Atoms(blocks).get_atomic_numbers(),
                                    ratio_of_covalent_radii=0.8) # 原子同士の近づいていい距離リスト
cellbounds = CellBounds(bounds={'phi': [30, 150], 'chi': [30, 150],
                                'psi': [30, 150], 'a': [3, 10],
                                'b': [3, 10], 'c': [3, 10]}) # cellの可動範囲

slab = Atoms('', pbc=True) # 原子を詰めるための雛形を用意します
sg = StartGenerator(slab, blocks, blmin, box_volume=box_volume,
                    cellbounds=cellbounds, 
                    number_of_variable_cell_vectors=3,
                    test_too_far=False)
```


```python
atoms = sg.get_new_candidate()
atoms
```




    Atoms(symbols='Ti4O8', pbc=True, cell=[[3.399468661931798, 0.0, 0.0], [-0.7876779702811255, 5.57401653937108, 0.0], [-2.4380262757553357, -1.7769862187552858, 7.59946811350224]], tags=...)




```python
v = view(atoms, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


    



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'Ti'), value='All…


### 1-2. 分子結晶


```python
blocks = ['H2O'] * 8
box_volume = 20 * 8 # どの程度の密度で詰めたいかに依存します

blmin = closest_distances_generator(atom_numbers=[1,8], # HとO
                                    ratio_of_covalent_radii=1.2) # 大きめにとると分子同士がくっつきにくくなります
cellbounds = CellBounds(bounds={'phi': [30, 150], 'chi': [30, 150],
                                'psi': [30, 150], 'a': [3, 10],
                                'b': [3, 10], 'c': [3, 10]})

slab = Atoms('', pbc=True)
sg = StartGenerator(slab, blocks, blmin, box_volume=box_volume,
                    cellbounds=cellbounds, 
                    number_of_variable_cell_vectors=3,
                    test_too_far=False)
```


```python
atoms = sg.get_new_candidate()
atoms
```




    Atoms(symbols='OH2OH2OH2OH2OH2OH2OH2OH2', pbc=True, cell=[[5.351229037091787, 0.0, 0.0], [0.4781589609694113, 5.361917945042938, 0.0], [-0.2782740527196683, -1.1820004956289019, 5.576301915370867]], tags=...)




```python
v = view(atoms, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H'), value='All'…


## 2. pyxtalを用いたランダム構造生成

pyxtalライブラリ( https://github.com/qzhu2017/PyXtal )はMITライセンスで公開されているオープンソースソフトです。
現在も開発が活発に行われているので使用するバージョンには注意してください。

pyxtal.XRDが `numba&gt;=0.50.1` に依存しており、`numba 0.54.1` が `numpy&lt;1.21,&gt;=1.17` を要求するため、Matlantisのデフォルトnumpyバージョンではインストールできません。


```python
!pip install numpy==1.20
!pip install pyxtal
```


```python
from pyxtal import pyxtal
from ase.visualize import view
from pymatgen.io.ase import AseAtomsAdaptor
```


```python
!pip install -U numpy
```

### 2-1. 結晶

空間群を番号で指定して、対称性を保ちながら構造を生成します。
例えばルチル型はP42/mnm (136)に属します。


```python
crystal = pyxtal()
crystal.from_random(3, 136, ['Ti','O'], [4,8]) # 次元, 空間群, 元素 ,組成
crystal
```




    
    ------Crystal from random------
    Dimension: 3
    Composition: O8Ti4
    Group: P42/mnm (136)
    tetragonal lattice:   5.9228   5.9228   4.1791  90.0000  90.0000  90.0000
    Wyckoff sites:
    	Ti @ [ 0.2499 -0.2499  0.0000], WP [4g] Site [m.m2]
    	 O @ [ 0.0000  0.5000  0.2500], WP [4d] Site [-4..]
    	 O @ [ 0.0000  0.0000  0.0000], WP [2a] Site [2/m.2/m2/m]
    	 O @ [ 0.0000  0.0000  0.5000], WP [2b] Site [2/m.2/m2/m]




```python
atoms = crystal.to_ase()
atoms.wrap()
v = view(atoms, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'Ti'), value='All…


### 2-2. 分子結晶


```python
mol_crystal = pyxtal(molecular=True)
mol_crystal.from_random(3, 36, ['H2O'], [8])
mol_crystal
```




    
    ------Crystal from random------
    Dimension: 3
    Composition: [H2O]8
    Group: Cmc21 (36)
    orthorhombic lattice:   7.2906   4.8405  11.3199  90.0000  90.0000  90.0000
    Wyckoff sites:
    	H2O1         @ [ 0.7033  0.3840  0.9518]  WP [8b] Site [1] Euler [-138.0    7.9   21.3]




```python
v = view(mol_crystal.to_ase(), viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H'), value='All'…


### 2-3. クラスター

次元を0にし、空間群ではなく点群で指定します。何も対称性を入れない場合はC1 (1)を用います。


```python
cluster = pyxtal()
cluster.from_random(0, 1, ['Pt'], [13]) # 次元, 点群, 元素 ,組成
cluster
```




    
    ------Crystal from random------
    Dimension: 0
    Composition: Pt13
    Group: C1 (1)
    spherical lattice:   3.6669   3.6669   3.6669  90.0000  90.0000  90.0000
    Wyckoff sites:
    	Pt @ [-0.7798  0.3045 -0.2505], WP [1a] Site [1]
    	Pt @ [-0.2504  0.0355 -0.0971], WP [1a] Site [1]
    	Pt @ [ 0.1569 -0.0117  0.5000], WP [1a] Site [1]
    	Pt @ [-0.2366 -0.3016 -0.2748], WP [1a] Site [1]
    	Pt @ [ 0.6771 -0.6440  0.1037], WP [1a] Site [1]
    	Pt @ [ 0.4976 -0.7905 -0.2079], WP [1a] Site [1]
    	Pt @ [ 0.3236  0.7859 -0.3149], WP [1a] Site [1]
    	Pt @ [ 0.7292  0.4299  0.4916], WP [1a] Site [1]
    	Pt @ [-0.4418 -0.1537 -0.5743], WP [1a] Site [1]
    	Pt @ [ 0.6301 -0.0632  0.5775], WP [1a] Site [1]
    	Pt @ [ 0.5740 -0.1216 -0.7349], WP [1a] Site [1]
    	Pt @ [ 0.6156  0.1639  0.2582], WP [1a] Site [1]
    	Pt @ [-0.1879  0.2728  0.4535], WP [1a] Site [1]




```python
v = view(cluster.to_ase(), viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```


    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pt'), value='All'), D…

