Copyright 2023 Matlantis Corp. (formerly Preferred Computational Chemistry Inc.) and Preferred Networks Inc. as contributors to Matlantis contrib project

# 樹脂-金属界面の応力-ひずみ曲線計算

本ノートブックでは、金属スラブに挟まれた樹脂の引っ張りシミュレーションを実行し、応力-ひずみ曲線を作成します。

&lt;figure&gt;
  &lt;img src="./input/image.png" alt="Example Image"/&gt;
  &lt;figcaption&gt;Fig.1 - 引っ張りシミュレーションで得られた最終構造 (OVITOを用いて可視化)&lt;/figcaption&gt;
&lt;/figure&gt;

## 事前セットアップ

pfp-api-client, matlantis-featuresをインストールします。


```python
# ! pip install -U pfp-api-client
# ! pip install -U matlantis-features
```

## 1. パッケージのインポート、クラス・関数の定義


```python
import os
import numpy as np
import pandas as pd
import random
import pathlib
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from sklearn.linear_model import LinearRegression

from ase import Atoms
from ase import units
from ase.io import read, write, Trajectory
from ase.build import bulk, surface, sort
from ase.units import Pascal
from ase.visualize import view
from ase.constraints import FixSymmetry
from ase import neighborlist
from ase.constraints import FixAtoms

from rdkit import Chem
from rdkit.Chem import AllChem

from matlantis_features.utils.atoms_util import convert_atoms_to_upper
from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.common.opt import LBFGSASEOptFeature, FireLBFGSASEOptFeature
from matlantis_features.features.md import ASEMDSystem, LangevinIntegrator, MDFeature, MDExtensionBase, NPTIntegrator
from matlantis_features.features.md.md_extensions import DeformScheduler
from matlantis_features.utils.calculators import get_calculator, pfp_estimator_fn

from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
estimator_fn = pfp_estimator_fn(model_version='v5.0.0', calc_mode='pbe_plus_d3')

from pfcc_extras.visualize.view import view_ngl
from pfcc_extras.structure.ase_rdkit_converter import smiles_to_atoms
from pfcc_extras.structure.surface import makesurface

def density(atoms):
    vol = atoms.get_cell().volume *(1e-8)**3   # as cm3
    wt = sum(atoms.get_masses())    # as g/mol
    dens = wt / units.mol / vol
    return(dens)

class PrintCellShape(MDExtensionBase):
    def __init__(self, cell_log=None):
        self.cell_log = cell_log
    def __call__(self, system, integrator) -&gt; None:
        cell_par = system.ase_atoms.get_cell_lengths_and_angles()
        istep = system.current_total_step
        print(f"Dyn step {istep:4d} a {cell_par[0]:3.2f} b {cell_par[1]:3.2f} c {cell_par[2]:3.2f} alpha {cell_par[3]:3.2f} beta {cell_par[4]:3.2f} gamma {cell_par[5]:3.2f} ")
        if self.cell_log is not None:
            self.cell_log.append(cell_par)
            
try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()
    
input_dir = dir_path / "input"
output_dir = dir_path / "output"
```


    


## 2. 構造作成

真空層で挟まれた金属スラブモデルを作成し、樹脂を充填します。その構造を圧着することで初期構造を生成します。

### 2-1 入出力ファイルパス


```python
# input file
# 樹脂構造
resin_path = input_dir / "PEO.xyz"

# output file
# 保存先ディレクトリ
results = "output_ja"
# 金属スラブ
metal_opt_path = output_dir / "model" / "slab_Fe.xyz"
# 最適後樹脂構造
resin_opt_path = output_dir / "model" / "resin_opt.xyz"
# 金属スラブ及び樹脂構造
slab_opt_path = output_dir / "model" / "slab.xyz"
# 高温圧着MDトラジェクトリ
npt1_path = output_dir / "md" / "npt1.traj"
# 室温平衡化MDトラジェクトリ
npt2_path = output_dir / "md" / "npt2.traj"
# 引っ張りMDトラジェクトリ
tensile_path = output_dir / "md" / "tensile.traj"
# 応力-ひずみ曲線図
ss_path = output_dir / "figure" / "ss.png"

# ディレクトリ作成
os.makedirs(output_dir / "model", exist_ok = True)
os.makedirs(output_dir / "md", exist_ok = True)
os.makedirs(output_dir / "figure", exist_ok = True)
```

### 2-2. 真空層を有する金属スラブ作成


```python
# バルク構造の緩和
print("バルク構造の緩和")
blk = bulk("Fe", "bcc")
blk.set_constraint([FixSymmetry(blk)])
opt = LBFGSASEOptFeature(n_run = 10000, fmax = 0.01, filter=True,
                         show_progress_bar=True, estimator_fn = estimator_fn)
result_opt = opt(blk)

# 真空層を有する表面モデルの作成
vacuum = 500
slab = result_opt.atoms.ase_atoms.copy()
slab.set_constraint([])
slab = makesurface(slab, miller_indices=(1,0,0), layers=10, rep=[8,8,1], vacuum=vacuum)
```

    バルク構造の緩和



      0%|          | 0/10001 [00:00&lt;?, ?it/s]



```python
# 層数の半分だけ下方にシフト
atoms = slab.copy()
# Sort atoms by their z-coordinate
z_coordinates = [atom.position[2] for atom in atoms]
unique_layers = np.unique(np.round(z_coordinates, decimals=3))
layer_z_value = unique_layers[int(len(unique_layers) / 2)]
atoms_in_layer = [atom for atom in slab if np.round(atom.position[2], decimals=3) == layer_z_value]
z_coor = atoms_in_layer[0].position[2]
slab_2 = slab.copy()
slab_2.positions += [0,0,-z_coor]
slab_2.wrap()

# 表面構造の緩和
print("表面構造の緩和")
indices = [atom.index for atom in slab_2 if atom.position[2] &lt;= 3]  # 3A以下を固定
slab_2.set_constraint(FixAtoms(indices))
result_opt_2 = opt(slab_2)
os.makedirs("output", exist_ok=True)
write(metal_opt_path, result_opt_2.atoms.ase_atoms, format='extxyz')

# 構造の可視化
v = view_ngl(result_opt_2.atoms.ase_atoms, representations=["ball+stick"])
v.view.control.spin([1,0,0], np.deg2rad(-90))
display(v)
```

    表面構造の緩和



      0%|          | 0/10001 [00:00&lt;?, ?it/s]



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Fe'), value='All'), D…


### 2-3. 樹脂の挿入


```python
atoms_organic = read(resin_path)

print("樹脂構造の緩和")
opt = LBFGSASEOptFeature(n_run = 10000, fmax = 0.01, filter=False,
                         show_progress_bar=True, estimator_fn = estimator_fn)
result_opt_organic = opt(atoms_organic)
write(resin_opt_path, result_opt_organic.atoms.ase_atoms, format='xyz')

view_ngl(result_opt_organic.atoms.ase_atoms)
```

    樹脂構造の緩和



      0%|          | 0/10001 [00:00&lt;?, ?it/s]





    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H', 'C'), value=…




```python
# 真空層と同サイズの真空層セル (`vacuum_cell`)を作成
empty_cell = slab_2.copy()
slab_height = empty_cell.cell.cellpar()[2] - vacuum
slabatoms = len(empty_cell)
cellsize = empty_cell.cell.cellpar()[:3]
cellsize[2] = vacuum

vacuum_cell = Atoms(symbols="", 
                  pbc=True, 
                  cell=cellsize)

NoAtoms0=len(atoms_organic)

# 分子を詰め込んだときの密度と原子数を確認
target_density = 0.10 # unit = g/cm3
sum_mass = atoms_organic.get_masses().sum()
sum_atoms = len(atoms_organic)
vac_vol = vacuum_cell.get_volume()*(1e-8)**3 # as cm3
d0 = sum_mass / units.mol / vac_vol # 真空層に分子を詰め込んだ時の推定密度
mlp = int(np.floor(target_density / d0)) # 有機分子を何個挿入するかを計算
t_atoms = NoAtoms0 * mlp + slabatoms
print ("liquid phase initial density =", np.round(d0*mlp, 3), "g/cm3")
print ("slab + resin atoms =", t_atoms)

# 真空層セル (`vacuum_cell`)に樹脂を充填
target_atoms = t_atoms    
vac_vol = empty_cell.get_volume()*(1e-8)**3 # as cm3

filled_cell = empty_cell.copy()
No_atoms, n = 0, 0
while No_atoms &lt; target_atoms:
    # 真空層における分子の挿入座標をランダムに決定
    ranx = random.uniform(0,1)
    rany = random.uniform(0,1)
    ranz = random.uniform(0.04,0.96)                           # slabのz軸範囲を除外。
    xyz_position=np.matmul([ranx,rany,ranz], filled_cell.cell)
    # 分子の回転角をランダムに決定
    phi, theta, psi = Rotation.random().as_euler('zxz', degrees=True)
    
    # 挿入分子を読み込み、並進移動+回転
    # molNo = n
    # m = mol_list[molNo]
    m = result_opt_organic.atoms.ase_atoms
    m1 = m.copy()
    m1.euler_rotate(phi=phi, theta=theta, psi=psi)
    m1.positions += xyz_position
    
    # 力が閾値以下であれば挿入
    Fcheck = filled_cell + m1
    Fcheck.calc = get_calculator(estimator_fn)
    if (np.sum(Fcheck.get_forces()**2, axis=1)**0.5).max() &lt;= 100:
        filled_cell = filled_cell + m1
        n+=1
        No_atoms = len(filled_cell)
        print("No_of_atoms =", No_atoms, ", including slab", len(empty_cell), "atoms")    
    else:
        n+=0

assert mlp == n
```

    liquid phase initial density = 0.093 g/cm3
    slab + resin atoms = 2336
    No_of_atoms = 852 , including slab 640 atoms
    No_of_atoms = 1064 , including slab 640 atoms
    No_of_atoms = 1276 , including slab 640 atoms
    No_of_atoms = 1488 , including slab 640 atoms
    No_of_atoms = 1700 , including slab 640 atoms
    No_of_atoms = 1912 , including slab 640 atoms
    No_of_atoms = 2124 , including slab 640 atoms
    No_of_atoms = 2336 , including slab 640 atoms



```python
# 充填構造の最適化
indices = [atom.index for atom in filled_cell if atom.position[2] &lt;= 3.0]
c = FixAtoms(indices)
filled_cell.set_constraint(c)

opt = LBFGSASEOptFeature(n_run = 10000, fmax = 0.1, filter = True,
                         show_progress_bar=True, estimator_fn = estimator_fn)
result_opt_filled = opt(filled_cell)

# 構造の保存
result_opt_filled.atoms.ase_atoms.set_initial_charges()
write(slab_opt_path, result_opt_filled.atoms.ase_atoms, format='extxyz')

# 構造の可視化
v = view_ngl(result_opt_filled.atoms.ase_atoms, representations=["ball+stick"])
v.view.control.spin([1,0,0], np.deg2rad(-90))
display(v)
```


      0%|          | 0/10001 [00:00&lt;?, ?it/s]


    /home/jovyan/.py38/lib/python3.8/site-packages/ase/io/extxyz.py:1000: UserWarning: write_xyz() overwriting array "forces" present in atoms.arrays with stored results from calculator
      warnings.warn('write_xyz() overwriting array "{0}" present '



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H', 'Fe', 'C'), …


### 2-4. 圧着シミュレーション


```python
ase_atoms =read(slab_opt_path)
indices = [atom.index for atom in ase_atoms if atom.position[2] &lt;= 3]
c = FixAtoms(indices)
ase_atoms.set_constraint(c)

atoms = MatlantisAtoms.from_ase_atoms(ase_atoms)
atoms.rotate_atoms_to_upper()

# NPT
t_step = 1     # as fs
press = 0.000101325 * 1000 # as GPa 大気圧の1000倍
temperature = 473.0

system = ASEMDSystem(atoms)
system.init_temperature(temperature)

integrator = NPTIntegrator(timestep=t_step,
                           temperature=temperature,
                           pressure=press * units.GPa,
                           ttime=20,
                           pfactor=2e6*units.GPa*(units.fs**2),
                           mask=np.array([[0,0,0], [0,0,0], [0,0,1]]))

info = PrintCellShape()

md = MDFeature(integrator, n_run=75000, traj_file_name=npt1_path,
               traj_freq=100, estimator_fn = estimator_fn)

md(system, extensions=[(info, 1000)])
```

    The MD trajectory will be saved at /home/jovyan/work/2023/sscurve/contrib/output_ja/md/npt1.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    WARNING: NPT: Setting the center-of-mass momentum to zero (was 13.6739 25.2356 -22.7151)
    /home/jovyan/.py38/lib/python3.8/site-packages/ase/utils/__init__.py:62: FutureWarning: Please use atoms.cell.cellpar() instead
      warnings.warn(warning)


    Dyn step    0 a 19.45 b 19.45 c 517.89 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 1000 a 19.45 b 19.45 c 505.81 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 2000 a 19.45 b 19.45 c 470.82 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 3000 a 19.45 b 19.45 c 419.98 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 4000 a 19.45 b 19.45 c 361.05 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 5000 a 19.45 b 19.45 c 300.88 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 6000 a 19.45 b 19.45 c 244.74 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 7000 a 19.45 b 19.45 c 195.28 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 8000 a 19.45 b 19.45 c 153.42 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 9000 a 19.45 b 19.45 c 119.34 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 10000 a 19.45 b 19.45 c 93.07 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 11000 a 19.45 b 19.45 c 74.56 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 12000 a 19.45 b 19.45 c 63.65 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 13000 a 19.45 b 19.45 c 59.61 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 14000 a 19.45 b 19.45 c 62.15 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 15000 a 19.45 b 19.45 c 69.11 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 16000 a 19.45 b 19.45 c 77.25 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 17000 a 19.45 b 19.45 c 85.22 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 18000 a 19.45 b 19.45 c 92.66 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 19000 a 19.45 b 19.45 c 99.27 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 20000 a 19.45 b 19.45 c 104.96 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 21000 a 19.45 b 19.45 c 109.54 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 22000 a 19.45 b 19.45 c 112.80 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 23000 a 19.45 b 19.45 c 114.81 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 24000 a 19.45 b 19.45 c 115.47 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 25000 a 19.45 b 19.45 c 115.03 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 26000 a 19.45 b 19.45 c 113.16 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 27000 a 19.45 b 19.45 c 110.14 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 28000 a 19.45 b 19.45 c 106.17 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 29000 a 19.45 b 19.45 c 101.53 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 30000 a 19.45 b 19.45 c 96.31 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 31000 a 19.45 b 19.45 c 90.92 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 32000 a 19.45 b 19.45 c 85.38 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 33000 a 19.45 b 19.45 c 80.30 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 34000 a 19.45 b 19.45 c 75.66 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 35000 a 19.45 b 19.45 c 71.43 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 36000 a 19.45 b 19.45 c 68.38 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 37000 a 19.45 b 19.45 c 66.95 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 38000 a 19.45 b 19.45 c 66.89 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 39000 a 19.45 b 19.45 c 68.20 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 40000 a 19.45 b 19.45 c 70.01 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 41000 a 19.45 b 19.45 c 71.98 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 42000 a 19.45 b 19.45 c 73.51 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 43000 a 19.45 b 19.45 c 74.26 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 44000 a 19.45 b 19.45 c 74.57 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 45000 a 19.45 b 19.45 c 74.44 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 46000 a 19.45 b 19.45 c 73.82 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 47000 a 19.45 b 19.45 c 72.86 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 48000 a 19.45 b 19.45 c 71.72 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 49000 a 19.45 b 19.45 c 70.53 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 50000 a 19.45 b 19.45 c 69.64 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 51000 a 19.45 b 19.45 c 69.03 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 52000 a 19.45 b 19.45 c 68.86 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 53000 a 19.45 b 19.45 c 68.97 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 54000 a 19.45 b 19.45 c 69.43 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 55000 a 19.45 b 19.45 c 69.98 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 56000 a 19.45 b 19.45 c 70.35 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 57000 a 19.45 b 19.45 c 70.57 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 58000 a 19.45 b 19.45 c 70.74 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 59000 a 19.45 b 19.45 c 70.90 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 60000 a 19.45 b 19.45 c 71.03 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 61000 a 19.45 b 19.45 c 71.09 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 62000 a 19.45 b 19.45 c 70.92 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 63000 a 19.45 b 19.45 c 70.52 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 64000 a 19.45 b 19.45 c 69.99 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 65000 a 19.45 b 19.45 c 69.66 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 66000 a 19.45 b 19.45 c 69.81 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 67000 a 19.45 b 19.45 c 70.20 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 68000 a 19.45 b 19.45 c 70.42 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 69000 a 19.45 b 19.45 c 70.34 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 70000 a 19.45 b 19.45 c 70.05 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 71000 a 19.45 b 19.45 c 69.81 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 72000 a 19.45 b 19.45 c 69.50 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 73000 a 19.45 b 19.45 c 69.25 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 74000 a 19.45 b 19.45 c 69.18 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 75000 a 19.45 b 19.45 c 69.10 alpha 90.00 beta 90.00 gamma 109.50 





    MDFeatureResult(traj_path=PosixPath('/home/jovyan/work/2023/sscurve/contrib/output_ja/md/npt1.traj'), checkpoint_path='/tmp/matlantis_j006kati/tmpa7tbr__t.chkp', temp_dir=&lt;TemporaryDirectory '/tmp/matlantis_j006kati'&gt;)




```python
# 圧着時のc軸長の推移と最終構造を確認
traj = read(npt1_path, "::1")

c_length = [atoms.cell.cellpar()[2] for atoms in traj]
plt.plot(c_length)
plt.ylabel("c-axis (Å)")
plt.xlabel("MD steps (×100)")
plt.show()

v = view_ngl(traj[-1], representations=["ball+stick"])
v.view.control.spin([1,0,0], np.deg2rad(-90))
display(v)
```


    
![png](output_17_0.png)
    



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H', 'Fe', 'C'), …


## 3. 室温NPTシミュレーション


```python
ase_atoms = read(npt1_path, "-1")
atoms = MatlantisAtoms.from_ase_atoms(ase_atoms)
atoms.rotate_atoms_to_upper()

# NPT
t_step = 1     # as fs
temperature = 300.0

system = ASEMDSystem(atoms)
system.init_temperature(temperature)

integrator = NPTIntegrator(timestep=1.0,
                           temperature=temperature,
                           pressure=units.bar,
                           ttime=20,
                           pfactor=2e6*units.GPa*(units.fs**2),
                           mask=np.array([[0,0,0], [0,0,0], [0,0,1]]))

info = PrintCellShape()

md = MDFeature(integrator, n_run=50000, traj_file_name=npt2_path,
               traj_freq=100, estimator_fn=estimator_fn)
md(system, extensions=[(info, 1000)])
```

    The MD trajectory will be saved at /home/jovyan/work/2023/sscurve/contrib/output_ja/md/npt2.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    WARNING: NPT: Setting the center-of-mass momentum to zero (was 17.0962 22.0152 -49.3797)
    /home/jovyan/.py38/lib/python3.8/site-packages/ase/utils/__init__.py:62: FutureWarning: Please use atoms.cell.cellpar() instead
      warnings.warn(warning)


    Dyn step    0 a 19.45 b 19.45 c 69.10 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 1000 a 19.45 b 19.45 c 69.17 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 2000 a 19.45 b 19.45 c 69.20 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 3000 a 19.45 b 19.45 c 69.43 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 4000 a 19.45 b 19.45 c 69.48 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 5000 a 19.45 b 19.45 c 68.48 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 6000 a 19.45 b 19.45 c 67.77 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 7000 a 19.45 b 19.45 c 68.35 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 8000 a 19.45 b 19.45 c 69.04 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 9000 a 19.45 b 19.45 c 69.14 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 10000 a 19.45 b 19.45 c 67.53 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 11000 a 19.45 b 19.45 c 66.73 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 12000 a 19.45 b 19.45 c 67.41 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 13000 a 19.45 b 19.45 c 68.57 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 14000 a 19.45 b 19.45 c 68.52 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 15000 a 19.45 b 19.45 c 68.02 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 16000 a 19.45 b 19.45 c 67.87 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 17000 a 19.45 b 19.45 c 67.66 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 18000 a 19.45 b 19.45 c 67.42 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 44000 a 19.45 b 19.45 c 67.58 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 45000 a 19.45 b 19.45 c 67.77 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 46000 a 19.45 b 19.45 c 68.05 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 47000 a 19.45 b 19.45 c 68.05 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 48000 a 19.45 b 19.45 c 67.48 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 49000 a 19.45 b 19.45 c 67.26 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 50000 a 19.45 b 19.45 c 67.56 alpha 90.00 beta 90.00 gamma 109.50 





    MDFeatureResult(traj_path=PosixPath('/home/jovyan/work/2023/sscurve/contrib/output_ja/md/npt2.traj'), checkpoint_path='/tmp/matlantis_lusvzudy/tmpf3dhkql4.chkp', temp_dir=&lt;TemporaryDirectory '/tmp/matlantis_lusvzudy'&gt;)



## 4. 引っ張りシミュレーション


```python
ase_atoms = read(npt2_path, "-1")
atoms = MatlantisAtoms.from_ase_atoms(ase_atoms)
atoms.rotate_atoms_to_upper()

system = ASEMDSystem(atoms)
system.init_temperature(temperature)

# z軸方向に20%引っ張ります。
latt = np.array(system.ase_atoms.cell)
latt[2,2] = latt[2,2]*3.0

# 引っ張り計算。
integrator = LangevinIntegrator(timestep=1.0, temperature=300.0)

info = PrintCellShape()
deform = DeformScheduler(latt, 100000)

md = MDFeature(integrator, n_run=100000, traj_file_name=tensile_path,
               traj_freq=100, estimator_fn = estimator_fn)
md(system, extensions=[(info, 2500), (deform, 1)])
```

    The MD trajectory will be saved at /home/jovyan/work/2023/sscurve/contrib/output_ja/md/tensile.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    /home/jovyan/.py38/lib/python3.8/site-packages/ase/utils/__init__.py:62: FutureWarning: Please use atoms.cell.cellpar() instead
      warnings.warn(warning)


    Dyn step    0 a 19.45 b 19.45 c 67.56 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 1000 a 19.45 b 19.45 c 68.23 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 2000 a 19.45 b 19.45 c 68.91 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 3000 a 19.45 b 19.45 c 69.58 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 4000 a 19.45 b 19.45 c 70.26 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 5000 a 19.45 b 19.45 c 70.93 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 6000 a 19.45 b 19.45 c 71.61 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 7000 a 19.45 b 19.45 c 72.29 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 8000 a 19.45 b 19.45 c 72.96 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 9000 a 19.45 b 19.45 c 73.64 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 10000 a 19.45 b 19.45 c 74.31 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 11000 a 19.45 b 19.45 c 74.99 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 12000 a 19.45 b 19.45 c 75.66 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 13000 a 19.45 b 19.45 c 76.34 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 14000 a 19.45 b 19.45 c 77.01 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 15000 a 19.45 b 19.45 c 77.69 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 16000 a 19.45 b 19.45 c 78.37 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 17000 a 19.45 b 19.45 c 79.04 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 18000 a 19.45 b 19.45 c 79.72 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 19000 a 19.45 b 19.45 c 80.39 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 20000 a 19.45 b 19.45 c 81.07 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 21000 a 19.45 b 19.45 c 81.74 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 22000 a 19.45 b 19.45 c 82.42 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 23000 a 19.45 b 19.45 c 83.09 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 24000 a 19.45 b 19.45 c 83.77 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 25000 a 19.45 b 19.45 c 84.45 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 26000 a 19.45 b 19.45 c 85.12 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 27000 a 19.45 b 19.45 c 85.80 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 28000 a 19.45 b 19.45 c 86.47 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 29000 a 19.45 b 19.45 c 87.15 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 30000 a 19.45 b 19.45 c 87.82 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 31000 a 19.45 b 19.45 c 88.50 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 32000 a 19.45 b 19.45 c 89.17 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 33000 a 19.45 b 19.45 c 89.85 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 34000 a 19.45 b 19.45 c 90.53 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 35000 a 19.45 b 19.45 c 91.20 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 36000 a 19.45 b 19.45 c 91.88 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 37000 a 19.45 b 19.45 c 92.55 alpha 90.00 beta 90.00 gamma 109.50 
    Dyn step 38000 a 19.45 b 19.45 c 93.23 alpha 90.00 beta 90.00 gamma 109.50 


## 5. 応力-ひずみ曲線のプロット


```python
traj = Trajectory(tensile_path)

stress = []
zlen = []
for i in traj:
    stress.append(i.get_stress()[2])
    zlen.append(i.cell.cellpar()[2])
zlen = np.array(zlen)
strain = (zlen-zlen[0])/zlen[0]
stress = np.array(stress)/Pascal/10e6

concatenated_columns = np.concatenate((strain.reshape(-1, 1), stress.reshape(-1, 1)), axis=1)
data = pd.DataFrame(concatenated_columns, columns=['strain', 'stress'])

os.makedirs("output/figure", exist_ok = True)

window_size = 20
data['moving_average'] = data['stress'].rolling(window=window_size).mean()

# Linear regression with no intercept in specified strain range
linear_range_data = data[(data['strain'] &gt;= 0.0) &amp; (data['strain'] &lt;= 0.025)]
X = linear_range_data[['strain']]
y = linear_range_data['stress']
linear_model_no_intercept = LinearRegression(fit_intercept=False).fit(X, y)
linear_fit = linear_model_no_intercept.predict(X)

# Youngs Modulus
youngs_modulus = linear_model_no_intercept.coef_[0]
print(f"Youngs modulus:{youngs_modulus / 1000:.4f}(GPa)")

# Plot
plt.scatter(data['strain'], data['stress'], label='Original Data', alpha=0.25, s = 10)
plt.plot(data['strain'], data['moving_average'], label='Moving Average (20 points)', color='green')
plt.plot(linear_range_data['strain'], linear_fit, label='Linear Fit through Origin', color='red', linewidth=2)
plt.title('Stress-Strain Curve with Moving Average and Linear Fit')
plt.xlabel('Strain')
plt.ylabel('Stress (MPa)')
plt.legend()
plt.grid(True)
plt.show()
plt.savefig(ss_path)
```

    Youngs modulus:4.8634(GPa)



    
![png](output_23_1.png)
    
