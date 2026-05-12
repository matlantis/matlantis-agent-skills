# 組成とスーパーセルを限定した置換構造探索

このnotebookでは、比較的大きなスーパーセルを用いて置換構造探索を行う例を示します。
また、置換構造の組成を指定することで、探索範囲を限定する方法も紹介します。


```python
%pip install -U mtcsp[matlantis]
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain

from datetime import datetime
import io

from ase.data import atomic_numbers
import ase.io
import ase.visualize
import numpy as np

from mtcsp.analysis import standardize_atoms
from mtcsp.convex_hull_search import SystemConfig
from mtcsp.derivative_structure import append_derivative_structures_with_random_sampling
from mtcsp.derivative_structure import append_derivative_structures_with_supercell
from mtcsp.derivative_structure import DerivativeStructureSearchConfig
from mtcsp.derivative_structure import initialize_derivative_structure_search
from mtcsp.derivative_structure import perform_derivative_structure_search
from mtcsp.derivative_structure import prepare_base_structure_from_atoms
from mtcsp.derivative_structure import VACANCY_NUMBER
from mtcsp_derst import enumerate_derivative_structures_with_sublattice
```


## 探索の準備

### 探索システム設定


```python
experiment_name = "derivative_structure_search_with_composition_" + datetime.now().strftime(
    "%Y%m%d%H%M%S"
)

system_config = SystemConfig.from_experiment_name(experiment_name)
```


### 母構造

この例では六方晶ペロブスカイト構造BaNiO3([mp-19138](https://next-gen.materialsproject.org/materials/mp-19138))を母構造とします。


```python
file_bno = io.StringIO("""# generated using pymatgen
data_BaNiO3
_symmetry_space_group_name_H-M   P6_3/mmc
_cell_length_a   5.67658435
_cell_length_b   5.67658435
_cell_length_c   4.78201573
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   120.00000000
_symmetry_Int_Tables_number   194
_chemical_formula_structural   BaNiO3
_chemical_formula_sum   'Ba2 Ni2 O6'
_cell_volume   133.44915350
_cell_formula_units_Z   2
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
  2  '-x, -y, -z'
  3  'x-y, x, z+1/2'
  4  '-x+y, -x, -z+1/2'
  5  '-y, x-y, z'
  6  'y, -x+y, -z'
  7  '-x, -y, z+1/2'
  8  'x, y, -z+1/2'
  9  '-x+y, -x, z'
  10  'x-y, x, -z'
  11  'y, -x+y, z+1/2'
  12  '-y, x-y, -z+1/2'
  13  '-y, -x, -z+1/2'
  14  'y, x, z+1/2'
  15  '-x, -x+y, -z'
  16  'x, x-y, z'
  17  '-x+y, y, -z+1/2'
  18  'x-y, -y, z+1/2'
  19  'y, x, -z'
  20  '-y, -x, z'
  21  'x, x-y, -z+1/2'
  22  '-x, -x+y, z+1/2'
  23  'x-y, -y, -z'
  24  '-x+y, y, z'
loop_
 _atom_type_symbol
 _atom_type_oxidation_number
  Ba2+  2.0
  Ni4+  4.0
  O2-  -2.0
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Ba2+  Ba0  2  0.33333333  0.66666667  0.75000000  1
  Ni4+  Ni1  2  0.00000000  0.00000000  0.00000000  1
  O2-  O2  6  0.14655756  0.29311512  0.25000000  1
""")
parent_atoms = ase.io.read(file_bno, format="cif")
assert isinstance(parent_atoms, ase.Atoms)
primitive_atoms = standardize_atoms(parent_atoms, primitive=True)

ase.visualize.view(primitive_atoms, viewer="x3d")
```


### 置換構造の設定

元素の置換方法を設定します。
この例では、OサイトをOまたは空孔(`VACANCY_NUMBER`)で置換することでO欠損構造を生成します。
置換を行わない元素も含めて母構造に存在するすべての元素について記述してください。


```python
feasible_number_mapping = {
    atomic_numbers["Ba"]: [
        atomic_numbers["Ba"],
    ],
    atomic_numbers["Ni"]: [
        atomic_numbers["Ni"],
    ],
    atomic_numbers["O"]: [
        atomic_numbers["O"],
        VACANCY_NUMBER,
    ],
}

search_config = DerivativeStructureSearchConfig.from_feasible_number_mapping(
    feasible_number_mapping=feasible_number_mapping,
    experiment_name=experiment_name,
)
```


ここでは`composition`で各元素ごとの組成比を指定することで、組成$\mathrm{Ba}\mathrm{Ni}\mathrm{O}_{5/2}$に限定して置換構造を考えます。
`composition`には整数で組成比を与えて下さい。


```python
composition = {
    atomic_numbers["Ba"]: 2,
    atomic_numbers["Ni"]: 2,
    # Substitute one of the six O sites in Ba2Ni2O6 with a vacancy
    atomic_numbers["O"]: 5,
    VACANCY_NUMBER: 1,
}

base_structure = prepare_base_structure_from_atoms(
    atoms=primitive_atoms,
    feasible_number_mapping=feasible_number_mapping,
    composition=composition,
)
```


## 置換構造探索

ここまで準備した入力をもとに置換構造探索を行います。
まず、探索の初期化を行います。


```python
experiment = initialize_derivative_structure_search(
    search_config=search_config,
    system_config=system_config,
)
```


### 全列挙によるスーパーセル指定での初期構造の追加

与えられた組成とスーパーセルの範囲で、対称性を考慮して非等価な置換構造を`append_derivative_structures_with_supercell`を用いて初期構造に追加することができます。
得られる置換構造の数は、`enumerate_derivative_structures_with_sublattice`で事前に確認することができます。


```python
supercell = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 2],
]

num_structures = len(enumerate_derivative_structures_with_sublattice(base_structure, supercell))
print(f"Number of derivative structures with the specified supercell: {num_structures}")
```


```python
append_derivative_structures_with_supercell(
    experiment=experiment,
    base_structure=base_structure,
    supercell=supercell,
)
```


### ランダムサンプリングによるスーパーセル指定での初期構造の追加

与えられた組成を満たすために大きなスーパーセルが必要な場合、`append_derivative_structures`で置換構造を全列挙すると数が膨大になることがあります。
ここでは、置換構造のスーパーセルの取り方(3x3行列)を指定して、その取り方の範囲で置換構造をランダムに生成する方法を紹介します。

スーパーセルは得られる置換構造のサイト数(`len(primitive_atoms) * np.linalg.det(supercell)`)がおよそ300サイト以下の範囲で指定して下さい。


```python
# Specify a supercell such that the number of atoms is within an appropriate range.
# The j-th axis of the resulting supercell is obtained by summing the i-th axes of the original unit cell
# multiplied by supercell[i][j].
supercell = [
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2],
]
n_crystal_structures = 4
print(
    f"Generate {n_crystal_structures} derivative structures with {len(primitive_atoms) * np.around(np.abs(np.linalg.det(supercell))).astype(int)} sites"
)
```


```python
# Generate `n_crystal_structures` derivative structures randomly with the specified `supercell`
append_derivative_structures_with_random_sampling(
    experiment=experiment,
    base_structure=base_structure,
    supercell=supercell,
    n_crystal_structures=n_crystal_structures,
    seed=0,
)
```

### Executing the Search


By adding to the initial structure, you can perform the search as usual.

### 探索の実行

初期構造に追加した後は、通常通り探索を行うことができます。


```python
perform_derivative_structure_search(
    search_config=search_config,
    system_config=system_config,
)
```
