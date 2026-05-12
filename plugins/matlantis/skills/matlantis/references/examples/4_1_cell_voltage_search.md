# 置換構造探索を用いた平衡電位計算

平衡電位計算のためにLiCoO2由来の置換構造探索を行います。

## インストール


```python
%pip install -U mtcsp[matlantis]
%pip install pandas
```


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain

from datetime import datetime
import io

from ase.data import atomic_numbers
import ase.io
import ase.visualize
import pandas as pd

from mtcsp.analysis import standardize_atoms
from mtcsp.convex_hull_search import SystemConfig
from mtcsp.derivative_structure import append_derivative_structures
from mtcsp.derivative_structure import DerivativeStructureSearchConfig
from mtcsp.derivative_structure import initialize_derivative_structure_search
from mtcsp.derivative_structure import list_num_derivative_structures
from mtcsp.derivative_structure import perform_derivative_structure_search
from mtcsp.derivative_structure import prepare_base_structure_from_atoms
from mtcsp.derivative_structure import VACANCY_NUMBER
```


## 探索の準備

### 探索システム設定

システム側の設定をします。

- experiment_name: 探索名
- db_file: メタデータを保存するデータベースファイル、拡張子は`.journal`としてください。
- atoms_store_dir: 原子構造の情報を保存するディレクトリ


```python
experiment_name = "cell_voltage_" + datetime.now().strftime("%Y%m%d%H%M%S")

system_config = SystemConfig.from_experiment_name(experiment_name)
```


### 母構造

この例では三方晶LiCoO2([mp-22526](https://next-gen.materialsproject.org/materials/mp-22526))を母構造として設定します。


```python
file_lco = io.StringIO("""# generated using pymatgen
data_LiCoO2
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   4.99271302
_cell_length_b   4.99271302
_cell_length_c   4.99271302
_cell_angle_alpha   33.08240238
_cell_angle_beta   33.08240238
_cell_angle_gamma   33.08240238
_symmetry_Int_Tables_number   1
_chemical_formula_structural   LiCoO2
_chemical_formula_sum   'Li1 Co1 O2'
_cell_volume   33.00303392
_cell_formula_units_Z   1
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Li  Li0  1  0.00000000  0.00000000  0.00000000  1
  Co  Co1  1  0.50000000  0.50000000  0.50000000  1
  O  O2  1  0.23958700  0.23958700  0.23958700  1
  O  O3  1  0.76041300  0.76041300  0.76041300  1
""")
parent_atoms = ase.io.read(file_lco, format="cif")
assert isinstance(parent_atoms, ase.Atoms)
ase.visualize.view(parent_atoms, viewer="x3d")
```


置換探索では母構造を primitive cell で取る必要があります。


```python
primitive_atoms = standardize_atoms(parent_atoms, primitive=True)
ase.visualize.view(primitive_atoms, viewer="x3d")
```


### 置換構造の設定

元素の置換方法を設定します。
この例では、LiサイトをLiまたは空孔(`VACANCY_NUMBER`)で置換することでLi欠損構造を生成します。
置換を行わない元素も含めて母構造に存在するすべての元素について記述してください。


```python
feasible_number_mapping = {
    atomic_numbers["Li"]: [
        atomic_numbers["Li"],
        VACANCY_NUMBER,
    ],
    atomic_numbers["Co"]: [
        atomic_numbers["Co"],
    ],
    atomic_numbers["O"]: [
        atomic_numbers["O"],
    ],
}

search_config = DerivativeStructureSearchConfig.from_feasible_number_mapping(
    feasible_number_mapping=feasible_number_mapping,
    experiment_name=experiment_name,
)
```


Primitive cell の何倍までの supercell を取るかを示す `max_index` 設定します。
`max_index` を大きくすると置換構造数が指数関数的に大きくなってしまうので設定には注意が必要です。
以下では index をいくつにすると合計の構造数 `total_num_structures` がいくつになるかを示しています。


```python
base_structure = prepare_base_structure_from_atoms(
    atoms=primitive_atoms,
    feasible_number_mapping=feasible_number_mapping,
)
data = list_num_derivative_structures(base_structure=base_structure)
pd.DataFrame(data).set_index("index")
```


`total_num_structures` をみて max_index を設定します。
この例では数分で探索が終了するように100構造程度を評価するように設定します。


```python
max_index = 5
```


## 置換構造探索

ここまで準備した入力をもとに置換構造探索を行います。
まず、探索の初期化を行います。
初期化中に置換構造を列挙し探索リストに追加していくので、構造数が多い場合には少し時間がかかります。


```python
experiment = initialize_derivative_structure_search(
    search_config=search_config,
    system_config=system_config,
)
append_derivative_structures(
    experiment=experiment,
    base_structure=base_structure,
    max_index=max_index,
)
```


実際の探索を行います。
探索はプロセス並列化して実行します。
`parallelism` は並列プロセス数を指定します。
通常 notebook では 6, 2x notebook では 12, 4x notebook では 24 程度に設定してください。


```python
perform_derivative_structure_search(
    search_config=search_config,
    system_config=system_config,
    parallelism=6,
)
```


以上で探索は終了です。
探索結果の結果の解析は`cell_voltage_analysis.ipynb`で行います。
