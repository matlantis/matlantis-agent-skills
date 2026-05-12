Copyright 2024 Matlantis Corp. (formerly Preferred Computational Chemistry Inc.) as contributors to Matlantis contrib project
# エポキシ樹脂（Bisphenol-A dimer diglycidyl ether）の熱分解シミュレーション
- エポキシ樹脂の熱分解反応は、プラスチックリサイクルのプロセス設計において重要です。通常、ポリマー分子の熱分解は高温で進行し、様々な小分子が生成されます。このような反応過程を扱う方法として、反応力場（ReaxFF）を用いた古典分子動力学（MD）シミュレーションが広く用いられています。これにより、分解反応を詳細に解析することが可能です。一方、ReaxFFは力場の作成が難しいという課題があります。そこで、学習済みの汎用ポテンシャルであるPFPを用いて同様の計算を行い、結果の検証を行いました（HPに公開済の事例から、PFPのバージョンと分子構造が変更されています）。  
https://matlantis.com/calculation/thermal-decomposition-of-epoxy-molecules  

- モデル化合物: Bisphenol-A dimer diglycidyl ether (formula: C&lt;sub&gt;39&lt;/sub&gt;H&lt;sub&gt;44&lt;/sub&gt;O&lt;sub&gt;7&lt;/sub&gt;)  
![image.png](7d3fcc9a-0fda-4df1-b9b3-4b1b3d36463d.png)
- reaxFFを使用したオリジナルの報告（Z. Diao et al., &lt;I&gt;Journal of Analytical and Applied Pyrolysis&lt;/I&gt; &lt;B&gt;2013&lt;/B&gt;, &lt;I&gt;104&lt;/I&gt;, 618.）  
https://www.sciencedirect.com/science/article/abs/pii/S016523701300096X?via%3Dihub

- 計算の流れ  
1) ライブラリの読み込み  
2) 系の作成  
    - エポキシ樹脂モデルの読み込み・構造最適化
    - エポキシ樹脂モデル15分子の凝集構造作成
3) 熱分解シミュレーション  
    - 昇温速度違い（100, 250, 500 K/ps で 300 K から 2300 K まで）
    - 到達温度違い（500 ps/K で 300 K から 2800, 3300, 4300 K まで）
4) トラジェクトリー解析  
    - 昇温速度違い（元論文 Fig. 4）
    - 到達温度違い（元論文 Fig. 5）


## 1. ライブラリの読み込み


```python
import os
```


```python
## PFP-API
from pfp_api_client import Estimator, ASECalculator

## Matlantis Features
from matlantis_features.utils.calculators import pfp_estimator_fn
from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.common.opt import FireLBFGSASEOptFeature
from matlantis_features.features.md import (ASEMDSystem, MDFeature, MDExtensionBase,
                                            NVTBerendsenIntegrator, NPTIntegrator, NPTBerendsenIntegrator)
from matlantis_features.features.md.md_extensions import TemperatureScheduler

## PFCC_extras
from pfcc_extras.liquidgenerator.liquid_generator import LiquidGenerator
from pfcc_extras.structure.molecule import get_mol_list
from pfcc_extras.visualize.view import view_ngl
from pfcc_extras.structure.connectivity import exists_colliding_atom

## ASE
from ase import units
from ase.io import read, write, Trajectory
from ase.optimize import LBFGS, FIRE
from ase.md import MDLogger
from ase.md.nptberendsen import NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

## Python General
from collections import Counter
from joblib import Parallel, delayed
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from time import perf_counter
```


    



```python
## PFP model_version &amp; calc_mode
model_version = "v6.0.0"
calc_mode = "PBE_PLUS_D3"
method = 'v6cU0d3'

## set PFP calculator
estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
calculator = ASECalculator(estimator)

## set pfp_estimator_fn for matlantis-feature
estimator_fn = pfp_estimator_fn(model_version=model_version, calc_mode=calc_mode)
```

## 2. 系の作成
- エポキシ樹脂のモデル:  
(&lt;I&gt;r&lt;/I&gt;)-1-(4-(2-(4-((&lt;I&gt;R&lt;/I&gt;)-oxiran-2-ylmethoxy)phenyl)propan-2-yl)phenoxy)-3-(4-(2-(4-((&lt;I&gt;S&lt;/I&gt;)-oxiran-2-ylmethoxy)phenyl)propan-2-yl)phenoxy)propan-2-ol  
https://pubchem.ncbi.nlm.nih.gov/compound/14367528 (立体化学の情報なし)
- モデル構造はAvogadro2を使用、立体化学を考慮して作成
- pfcc_extras 内の LiquidGenerator (Torch version)を使用して、凝集構造作成（&lt;B&gt;pfcc_extrasのインストールが必要&lt;/B&gt;）

### エポキシ樹脂モデルの読み込み・構造最適化


```python
## read the epoxy structure ('BisphenolA-dimer-diglycidyl-ether.xyz')
atoms_epoxy = read('./input/BisphenolA-dimer-diglycidyl-ether.xyz')
view_ngl(atoms_epoxy, representations=['ball+stick'])
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'C', 'H'), value=…




```python
## rough optimization
dirout = 'output'
fname = f'C39H44O7_opt_{method}'
os.makedirs(dirout, exist_ok=True)

atoms_epoxy.calc = calculator
opt = FIRE(atoms_epoxy, logfile=f'{dirout}/{fname}.log', trajectory=f'{dirout}/{fname}.traj')
opt.run(fmax=0.01)
```




    True



### エポキシ樹脂モデル15分子の凝集構造作成
- 元論文（31.7×22.0×22.6 Å&lt;sup&gt;3&lt;/sup&gt;, 1.0 g cm&lt;sup&gt;-3&lt;/sup&gt;）に近い緩和構造を作成
- x:y:z = 3:2:2 となるように初期構造セルを設定し, LiquidGeneratorで初期凝集構造作成
- セル長比を維持しつつNPT-MDで平衡化
- NVT-MDによる平衡化



```python
# Packing epoxy mols. into a periodic box with x:y:z = 3:2:2
## LiquidGenerator Torch ver.
n_mol = 15
a_unit = 14.5
params = {
    "cell": [[a_unit*3,0,0], [0,a_unit*2,0], [0,0,a_unit*2]],  ## cell format for Torch ver.
    "composition": [atoms_epoxy] * n_mol,
    "epochs": 100
}

## run LiquidGenerator (~5 min.)
generator = LiquidGenerator("torch", **params)
atoms_agg = generator.run()
```

    step  score  cell_x  cell_y  cell_z


    /home/jovyan/.py39/lib/python3.9/site-packages/torch/functional.py:507: UserWarning: torch.meshgrid: in an upcoming release, it will be required to pass the indexing argument. (Triggered internally at ../aten/src/ATen/native/TensorShape.cpp:3549.)
      return _VF.meshgrid(tensors, **kwargs)  # type: ignore[attr-defined]


       0   28215.98 43.50 29.00 29.00
       1   27910.40 43.50 29.00 29.00
       2   28146.26 43.50 29.00 29.00
       3   27022.24 43.50 29.00 29.00
       4   26758.96 43.50 29.00 29.00
       5   27441.61 43.50 29.00 29.00
       6   27487.13 43.50 29.00 29.00
       7   27210.45 43.50 29.00 29.00
       8   27503.83 43.50 29.00 29.00
       9   27405.76 43.50 29.00 29.00
      10   27085.61 43.50 29.00 29.00
      11   26512.78 43.50 29.00 29.00
      12   26004.75 43.50 29.00 29.00
      13   25651.72 43.50 29.00 29.00
      14   25611.56 43.50 29.00 29.00
      15   25528.54 43.50 29.00 29.00
      16   25514.75 43.50 29.00 29.00
      17   25229.21 43.50 29.00 29.00
      18   25369.65 43.50 29.00 29.00
      19   25092.83 43.50 29.00 29.00
      20   24996.01 43.50 29.00 29.00
      21   24948.09 43.50 29.00 29.00
      22   24833.09 43.50 29.00 29.00
      23   24844.09 43.50 29.00 29.00
      24   24777.93 43.50 29.00 29.00
      25   24740.56 43.50 29.00 29.00
      26   24757.24 43.50 29.00 29.00
      27   24748.17 43.50 29.00 29.00
      28   24869.59 43.50 29.00 29.00
      29   24711.13 43.50 29.00 29.00
      30   24663.55 43.50 29.00 29.00
      31   24639.65 43.50 29.00 29.00
      32   24632.05 43.50 29.00 29.00
      33   24627.64 43.50 29.00 29.00
      34   24622.01 43.50 29.00 29.00
      35   24600.74 43.50 29.00 29.00
      36   24589.52 43.50 29.00 29.00
      37   24574.86 43.50 29.00 29.00
      38   24571.55 43.50 29.00 29.00
      39   24569.40 43.50 29.00 29.00
      40   24565.23 43.50 29.00 29.00
      41   24563.93 43.50 29.00 29.00
      42   24566.68 43.50 29.00 29.00
      43   24565.03 43.50 29.00 29.00
      44   24557.81 43.50 29.00 29.00
      45   24549.81 43.50 29.00 29.00
      46   24543.21 43.50 29.00 29.00
      47   24543.67 43.50 29.00 29.00
      48   24543.33 43.50 29.00 29.00
      49   24540.64 43.50 29.00 29.00



```python
## check whether generated atoms collide each other or not 
exists_colliding_atom(atoms_agg)
```




    False




```python
## check generated system
dens = sum(atoms_agg.get_masses()) / units.mol / (atoms_agg.get_volume()*(1e-8)**3)
print(f'density = {dens:.3f} g cm-3')
view_ngl(atoms_agg, representations=['ball+stick'])
```

    density = 0.425 g cm-3





    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'C', 'H'), value=…




```python
## optimization of atomic positions
sysname = '15C39H44O7'
fname = f'{sysname}_opt_{method}'
os.makedirs(dirout, exist_ok=True)
atoms_agg.calc = calculator
fmax = 0.03  ## rough optimization

## run opt
opt = FIRE(atoms_agg)
opt = FIRE(atoms_agg, logfile=f'{dirout}/{fname}.log', trajectory=f'{dirout}/{fname}.traj')
opt.run(fmax=fmax)
```




    True




```python
## set MD-eq print&amp;logger
def print_dyn():  ## general
    line  = f'Dyn  step: {dyn.get_number_of_steps():8d}, '
    line += f'time[ps]: {dyn.get_time() / units.fs / 1000.:7.2f}, '    
    line += f'Etot[eV]: {atoms.get_total_energy():11.4f}, '
    line += f'Epot[eV]: {atoms.get_potential_energy():11.4f}, '
    line += f'Ekin[eV]: {atoms.get_kinetic_energy():9.4f}, '
    line += f'T[K]: {atoms.get_temperature():7.2f}, '
    line += f'density[g/cm3]: {atoms.get_masses().sum()/units.mol/(atoms.get_volume()*(1e-8)**3):5.3f}, '
    line += f'calc_time[min]: {(perf_counter() - t_start)/60.:8.2f}'
    print(line)

class AddMDLogger(MDLogger):  ## density&amp;cell parameters
    def __call__(self):
        line  = f'{dyn.get_number_of_steps():8d}, '
        line += f'{dyn.get_time() / units.fs / 1000.:7.2f}, '    
        line += f'{atoms.get_total_energy():11.4f}, '
        line += f'{atoms.get_temperature():7.2f}, '
        line += f'{atoms.get_masses().sum()/units.mol/(atoms.get_volume()*(1e-8)**3):5.3f}, '
        line += f'{self.atoms.cell.cellpar()[0:6]}, '
        line += f'{(perf_counter() - t_start)/60.:8.2f}\n'
        if dyn.get_number_of_steps() == 0:
            hdr = 'step,time[ps],Etot[eV],T[K],density[g/cm3],cell_params,calc_time[min]\n'
            self.logfile.write(hdr)
        self.logfile.write(line)
        self.logfile.flush()
```


```python
# NPT-MD (100 ps, Berendsen thermostat and barostat at 300 K, 1 bar)
## set parameters NPT-MD
dirout = 'output'
sysname = '15C39H44O7'
fname = f'{sysname}_npt-eq_{method}'
traj_path = f'{dirout}/{fname}.traj'
log_path = f'{dirout}/{fname}.log'
addlog_path = f'{dirout}/{fname}-add.log'
temp = 300.
dt = 1.0  ## 1 fs
steps = 100_000  ## 100 ps
atoms = atoms_agg.copy()
atoms.calc = calculator

## set NPTBerendsen dynamics
dyn = NPTBerendsen(
    atoms=atoms, 
    timestep=dt * units.fs, 
    temperature_K=temp,
    taut=30 * units.fs, 
    pressure_au=1.01325 * units.bar,
    taup=500.0 * units.fs, 
    compressibility_au=5e-5 / units.bar,
    trajectory=f'{traj_path}', 
    loginterval=100,
    )
dyn.attach(print_dyn, interval=100)
dyn.attach(MDLogger(
            dyn=dyn, atoms=atoms, logfile=f'{log_path}', 
            header=True, stress=True, peratom=False, mode="w"), 
            interval=100
          )
dyn.attach(AddMDLogger(
            dyn=dyn, atoms=atoms, logfile=f'{addlog_path}', 
            header=False, mode="w"),
            interval=100
          ) 

## run MD
MaxwellBoltzmannDistribution(
    atoms=atoms,
    temperature_K=temp,
    rng=np.random.RandomState(42)
    )
Stationary(atoms)
t_start = perf_counter()
dyn.run(steps)
```

    Dyn  step:        0, time[ps]:    0.00, Etot[eV]:  -6609.9879, Epot[eV]:  -6662.1609, Ekin[eV]:   52.1730, T[K]:  298.98, density[g/cm3]: 0.425, calc_time[min]:     0.00
    Dyn  step:      100, time[ps]:    0.10, Etot[eV]:  -6572.9789, Epot[eV]:  -6619.0284, Ekin[eV]:   46.0495, T[K]:  263.89, density[g/cm3]: 0.427, calc_time[min]:     0.80
    Dyn  step:      200, time[ps]:    0.20, Etot[eV]:  -6565.5684, Epot[eV]:  -6618.4116, Ekin[eV]:   52.8432, T[K]:  302.82, density[g/cm3]: 0.429, calc_time[min]:     1.60
    Dyn  step:      300, time[ps]:    0.30, Etot[eV]:  -6563.4217, Epot[eV]:  -6617.0061, Ekin[eV]:   53.5844, T[K]:  307.07, density[g/cm3]: 0.431, calc_time[min]:     2.39
    Dyn  step:      400, time[ps]:    0.40, Etot[eV]:  -6562.8542, Epot[eV]:  -6614.8160, Ekin[eV]:   51.9618, T[K]:  297.77, density[g/cm3]: 0.432, calc_time[min]:     3.17
    Dyn  step:      500, time[ps]:    0.50, Etot[eV]:  -6561.8918, Epot[eV]:  -6617.4888, Ekin[eV]:   55.5969, T[K]:  318.61, density[g/cm3]: 0.433, calc_time[min]:     3.97
    Dyn  step:      600, time[ps]:    0.60, Etot[eV]:  -6562.9333, Epot[eV]:  -6612.8713, Ekin[eV]:   49.9380, T[K]:  286.18, density[g/cm3]: 0.434, calc_time[min]:     4.78
    Dyn  step:      700, time[ps]:    0.70, Etot[eV]:  -6563.5809, Epot[eV]:  -6616.0952, Ekin[eV]:   52.5142, T[K]:  300.94, density[g/cm3]: 0.435, calc_time[min]:     5.59
    Dyn  step:      800, time[ps]:    0.80, Etot[eV]:  -6564.2729, Epot[eV]:  -6617.7952, Ekin[eV]:   53.5223, T[K]:  306.72, density[g/cm3]: 0.437, calc_time[min]:     6.41
    Dyn  step:      900, time[ps]:    0.90, Etot[eV]:  -6563.9382, Epot[eV]:  -6614.8632, Ekin[eV]:   50.9250, T[K]:  291.83, density[g/cm3]: 0.437, calc_time[min]:     7.23
    Dyn  step:     1000, time[ps]:    1.00, Etot[eV]:  -6564.1350, Epot[eV]:  -6618.1161, Ekin[eV]:   53.9810, T[K]:  309.35, density[g/cm3]: 0.439, calc_time[min]:     8.05
    Dyn  step:     1100, time[ps]:    1.10, Etot[eV]:  -6564.0912, Epot[eV]:  -6615.9716, Ekin[eV]:   51.8804, T[K]:  297.31, density[g/cm3]: 0.440, calc_time[min]:     8.87
    Dyn  step:     1200, time[ps]:    1.20, Etot[eV]:  -6564.1335, Epot[eV]:  -6617.1887, Ekin[eV]:   53.0552, T[K]:  304.04, density[g/cm3]: 0.442, calc_time[min]:     9.69
    Dyn  step:     1300, time[ps]:    1.30, Etot[eV]:  -6564.6420, Epot[eV]:  -6616.7207, Ekin[eV]:   52.0787, T[K]:  298.44, density[g/cm3]: 0.444, calc_time[min]:    10.51
    Dyn  step:     1400, time[ps]:    1.40, Etot[eV]:  -6564.3096, Epot[eV]:  -6617.1812, Ekin[eV]:   52.8716, T[K]:  302.99, density[g/cm3]: 0.445, calc_time[min]:    11.33
    Dyn  step:     1500, time[ps]:    1.50, Etot[eV]:  -6564.0768, Epot[eV]:  -6616.0240, Ekin[eV]:   51.9473, T[K]:  297.69, density[g/cm3]: 0.448, calc_time[min]:    12.16
    Dyn  step:     1600, time[ps]:    1.60, Etot[eV]:  -6564.9108, Epot[eV]:  -6616.8155, Ekin[eV]:   51.9047, T[K]:  297.45, density[g/cm3]: 0.451, calc_time[min]:    13.00
    Dyn  step:     1700, time[ps]:    1.70, Etot[eV]:  -6564.8741, Epot[eV]:  -6616.9026, Ekin[eV]:   52.0285, T[K]:  298.16, density[g/cm3]: 0.451, calc_time[min]:    13.83
    Dyn  step:     1800, time[ps]:    1.80, Etot[eV]:  -6565.7237, Epot[eV]:  -6617.6742, Ekin[eV]:   51.9504, T[K]:  297.71, density[g/cm3]: 0.453, calc_time[min]:    14.68
    Dyn  step:     1900, time[ps]:    1.90, Etot[eV]:  -6566.8379, Epot[eV]:  -6619.4983, Ekin[eV]:   52.6604, T[K]:  301.78, density[g/cm3]: 0.452, calc_time[min]:    15.50
    Dyn  step:     2000, time[ps]:    2.00, Etot[eV]:  -6565.5610, Epot[eV]:  -6618.2445, Ekin[eV]:   52.6835, T[K]:  301.91, density[g/cm3]: 0.453, calc_time[min]:    16.33
    Dyn  step:     2100, time[ps]:    2.10, Etot[eV]:  -6565.9991, Epot[eV]:  -6619.0050, Ekin[eV]:   53.0059, T[K]:  303.76, density[g/cm3]: 0.452, calc_time[min]:    17.17
    Dyn  step:     2200, time[ps]:    2.20, Etot[eV]:  -6565.6602, Epot[eV]:  -6617.8947, Ekin[eV]:   52.2345, T[K]:  299.34, density[g/cm3]: 0.453, calc_time[min]:    18.00
    Dyn  step:     2300, time[ps]:    2.30, Etot[eV]:  -6566.3223, Epot[eV]:  -6619.1171, Ekin[eV]:   52.7949, T[K]:  302.55, density[g/cm3]: 0.454, calc_time[min]:    18.84
    Dyn  step:     2400, time[ps]:    2.40, Etot[eV]:  -6566.5058, Epot[eV]:  -6618.6208, Ekin[eV]:   52.1150, T[K]:  298.65, density[g/cm3]: 0.457, calc_time[min]:    19.68
    Dyn  step:     2500, time[ps]:    2.50, Etot[eV]:  -6565.9573, Epot[eV]:  -6617.0337, Ekin[eV]:   51.0764, T[K]:  292.70, density[g/cm3]: 0.459, calc_time[min]:    20.53
    Dyn  step:     2600, time[ps]:    2.60, Etot[eV]:  -6566.3250, Epot[eV]:  -6618.2440, Ekin[eV]:   51.9190, T[K]:  297.53, density[g/cm3]: 0.461, calc_time[min]:    21.36
    Dyn  step:     2700, time[ps]:    2.70, Etot[eV]:  -6566.5375, Epot[eV]:  -6619.3732, Ekin[eV]:   52.8357, T[K]:  302.78, density[g/cm3]: 0.462, calc_time[min]:    22.20
    Dyn  step:     2800, time[ps]:    2.80, Etot[eV]:  -6566.1844, Epot[eV]:  -6619.3832, Ekin[eV]:   53.1988, T[K]:  304.86, density[g/cm3]: 0.460, calc_time[min]:    23.05
    Dyn  step:     2900, time[ps]:    2.90, Etot[eV]:  -6567.0611, Epot[eV]:  -6619.3399, Ekin[eV]:   52.2788, T[K]:  299.59, density[g/cm3]: 0.463, calc_time[min]:    23.88
    Dyn  step:     3000, time[ps]:    3.00, Etot[eV]:  -6567.4373, Epot[eV]:  -6619.6474, Ekin[eV]:   52.2101, T[K]:  299.20, density[g/cm3]: 0.465, calc_time[min]:    24.73
    Dyn  step:     3100, time[ps]:    3.10, Etot[eV]:  -6567.9748, Epot[eV]:  -6621.0781, Ekin[eV]:   53.1033, T[K]:  304.32, density[g/cm3]: 0.465, calc_time[min]:    25.58
    Dyn  step:     3200, time[ps]:    3.20, Etot[eV]:  -6567.3865, Epot[eV]:  -6619.6521, Ekin[eV]:   52.2656, T[K]:  299.51, density[g/cm3]: 0.467, calc_time[min]:    26.43
    Dyn  step:     3300, time[ps]:    3.30, Etot[eV]:  -6567.4557, Epot[eV]:  -6619.4930, Ekin[eV]:   52.0373, T[K]:  298.21, density[g/cm3]: 0.469, calc_time[min]:    27.29
    Dyn  step:     3400, time[ps]:    3.40, Etot[eV]:  -6567.2980, Epot[eV]:  -6619.0371, Ekin[eV]:   51.7391, T[K]:  296.50, density[g/cm3]: 0.467, calc_time[min]:    28.14
    Dyn  step:     3500, time[ps]:    3.50, Etot[eV]:  -6567.1615, Epot[eV]:  -6618.5584, Ekin[eV]:   51.3969, T[K]:  294.54, density[g/cm3]: 0.469, calc_time[min]:    29.01
    Dyn  step:     3600, time[ps]:    3.60, Etot[eV]:  -6566.3219, Epot[eV]:  -6618.0995, Ekin[eV]:   51.7776, T[K]:  296.72, density[g/cm3]: 0.471, calc_time[min]:    29.86
    Dyn  step:     3700, time[ps]:    3.70, Etot[eV]:  -6566.9590, Epot[eV]:  -6619.9028, Ekin[eV]:   52.9438, T[K]:  303.40, density[g/cm3]: 0.472, calc_time[min]:    30.73
    Dyn  step:     3800, time[ps]:    3.80, Etot[eV]:  -6567.7018, Epot[eV]:  -6621.5528, Ekin[eV]:   53.8509, T[K]:  308.60, density[g/cm3]: 0.474, calc_time[min]:    31.57
    Dyn  step:     3900, time[ps]:    3.90, Etot[eV]:  -6567.5685, Epot[eV]:  -6619.1311, Ekin[eV]:   51.5625, T[K]:  295.49, density[g/cm3]: 0.476, calc_time[min]:    32.44
    Dyn  step:     4000, time[ps]:    4.00, Etot[eV]:  -6568.2487, Epot[eV]:  -6620.5144, Ekin[eV]:   52.2657, T[K]:  299.52, density[g/cm3]: 0.477, calc_time[min]:    33.30
    Dyn  step:     4100, time[ps]:    4.10, Etot[eV]:  -6567.7580, Epot[eV]:  -6619.9685, Ekin[eV]:   52.2105, T[K]:  299.20, density[g/cm3]: 0.478, calc_time[min]:    34.17
    Dyn  step:     4200, time[ps]:    4.20, Etot[eV]:  -6567.8796, Epot[eV]:  -6620.0603, Ekin[eV]:   52.1807, T[K]:  299.03, density[g/cm3]: 0.478, calc_time[min]:    35.03
    Dyn  step:     4300, time[ps]:    4.30, Etot[eV]:  -6568.1498, Epot[eV]:  -6621.2364, Ekin[eV]:   53.0866, T[K]:  304.22, density[g/cm3]: 0.480, calc_time[min]:    35.87
    Dyn  step:     4400, time[ps]:    4.40, Etot[eV]:  -6567.8280, Epot[eV]:  -6618.7995, Ekin[eV]:   50.9715, T[K]:  292.10, density[g/cm3]: 0.480, calc_time[min]:    36.73
    Dyn  step:     4500, time[ps]:    4.50, Etot[eV]:  -6568.0156, Epot[eV]:  -6619.5551, Ekin[eV]:   51.5396, T[K]:  295.35, density[g/cm3]: 0.481, calc_time[min]:    37.60
    Dyn  step:     4600, time[ps]:    4.60, Etot[eV]:  -6568.4472, Epot[eV]:  -6621.1556, Ekin[eV]:   52.7085, T[K]:  302.05, density[g/cm3]: 0.482, calc_time[min]:    38.46
    Dyn  step:     4700, time[ps]:    4.70, Etot[eV]:  -6568.8602, Epot[eV]:  -6621.4293, Ekin[eV]:   52.5691, T[K]:  301.25, density[g/cm3]: 0.482, calc_time[min]:    39.32
    Dyn  step:     4800, time[ps]:    4.80, Etot[eV]:  -6569.1416, Epot[eV]:  -6620.7896, Ekin[eV]:   51.6480, T[K]:  295.98, density[g/cm3]: 0.486, calc_time[min]:    40.20
    Dyn  step:     4900, time[ps]:    4.90, Etot[eV]:  -6568.7907, Epot[eV]:  -6620.4847, Ekin[eV]:   51.6940, T[K]:  296.24, density[g/cm3]: 0.489, calc_time[min]:    41.07
    Dyn  step:     5000, time[ps]:    5.00, Etot[eV]:  -6568.9138, Epot[eV]:  -6620.8227, Ekin[eV]:   51.9089, T[K]:  297.47, density[g/cm3]: 0.492, calc_time[min]:    41.93
    Dyn  step:     5100, time[ps]:    5.10, Etot[eV]:  -6569.3716, Epot[eV]:  -6621.3622, Ekin[eV]:   51.9906, T[K]:  297.94, density[g/cm3]: 0.492, calc_time[min]:    42.80
    Dyn  step:     5200, time[ps]:    5.20, Etot[eV]:  -6569.0881, Epot[eV]:  -6621.1038, Ekin[eV]:   52.0158, T[K]:  298.08, density[g/cm3]: 0.494, calc_time[min]:    43.68
    Dyn  step:     5300, time[ps]:    5.30, Etot[eV]:  -6569.1938, Epot[eV]:  -6621.3758, Ekin[eV]:   52.1819, T[K]:  299.04, density[g/cm3]: 0.496, calc_time[min]:    44.54
    Dyn  step:     5400, time[ps]:    5.40, Etot[eV]:  -6569.6656, Epot[eV]:  -6622.1813, Ekin[eV]:   52.5157, T[K]:  300.95, density[g/cm3]: 0.497, calc_time[min]:    45.43
    Dyn  step:     5500, time[ps]:    5.50, Etot[eV]:  -6569.3238, Epot[eV]:  -6620.6948, Ekin[eV]:   51.3710, T[K]:  294.39, density[g/cm3]: 0.496, calc_time[min]:    46.30
    Dyn  step:     5600, time[ps]:    5.60, Etot[eV]:  -6569.9005, Epot[eV]:  -6621.5741, Ekin[eV]:   51.6736, T[K]:  296.12, density[g/cm3]: 0.496, calc_time[min]:    47.16
    Dyn  step:     5700, time[ps]:    5.70, Etot[eV]:  -6569.8025, Epot[eV]:  -6623.3258, Ekin[eV]:   53.5233, T[K]:  306.72, density[g/cm3]: 0.497, calc_time[min]:    48.03
    Dyn  step:     5800, time[ps]:    5.80, Etot[eV]:  -6570.6769, Epot[eV]:  -6622.7547, Ekin[eV]:   52.0778, T[K]:  298.44, density[g/cm3]: 0.499, calc_time[min]:    48.92
    Dyn  step:     5900, time[ps]:    5.90, Etot[eV]:  -6569.8068, Epot[eV]:  -6622.6000, Ekin[eV]:   52.7932, T[K]:  302.54, density[g/cm3]: 0.500, calc_time[min]:    49.80
    Dyn  step:     6000, time[ps]:    6.00, Etot[eV]:  -6570.2403, Epot[eV]:  -6623.6235, Ekin[eV]:   53.3832, T[K]:  305.92, density[g/cm3]: 0.501, calc_time[min]:    50.70
    Dyn  step:     6100, time[ps]:    6.10, Etot[eV]:  -6569.2184, Epot[eV]:  -6619.9070, Ekin[eV]:   50.6886, T[K]:  290.48, density[g/cm3]: 0.503, calc_time[min]:    51.59
    Dyn  step:     6200, time[ps]:    6.20, Etot[eV]:  -6569.5912, Epot[eV]:  -6621.6173, Ekin[eV]:   52.0261, T[K]:  298.14, density[g/cm3]: 0.505, calc_time[min]:    52.47
    Dyn  step:     6300, time[ps]:    6.30, Etot[eV]:  -6570.2906, Epot[eV]:  -6621.5586, Ekin[eV]:   51.2680, T[K]:  293.80, density[g/cm3]: 0.508, calc_time[min]:    53.35
    Dyn  step:     6400, time[ps]:    6.40, Etot[eV]:  -6570.8172, Epot[eV]:  -6623.9179, Ekin[eV]:   53.1007, T[K]:  304.30, density[g/cm3]: 0.508, calc_time[min]:    54.24
    Dyn  step:     6500, time[ps]:    6.50, Etot[eV]:  -6571.0814, Epot[eV]:  -6624.0492, Ekin[eV]:   52.9678, T[K]:  303.54, density[g/cm3]: 0.509, calc_time[min]:    55.13
    Dyn  step:     6600, time[ps]:    6.60, Etot[eV]:  -6571.5881, Epot[eV]:  -6623.9047, Ekin[eV]:   52.3166, T[K]:  299.81, density[g/cm3]: 0.509, calc_time[min]:    56.02
    Dyn  step:     6700, time[ps]:    6.70, Etot[eV]:  -6572.1613, Epot[eV]:  -6623.1391, Ekin[eV]:   50.9778, T[K]:  292.13, density[g/cm3]: 0.512, calc_time[min]:    56.92
    Dyn  step:     6800, time[ps]:    6.80, Etot[eV]:  -6572.0711, Epot[eV]:  -6624.4864, Ekin[eV]:   52.4153, T[K]:  300.37, density[g/cm3]: 0.513, calc_time[min]:    57.80
    Dyn  step:     6900, time[ps]:    6.90, Etot[eV]:  -6571.3576, Epot[eV]:  -6623.1938, Ekin[eV]:   51.8362, T[K]:  297.05, density[g/cm3]: 0.514, calc_time[min]:    58.69
    Dyn  step:     7000, time[ps]:    7.00, Etot[eV]:  -6571.9007, Epot[eV]:  -6624.1406, Ekin[eV]:   52.2399, T[K]:  299.37, density[g/cm3]: 0.513, calc_time[min]:    59.61
    Dyn  step:     7100, time[ps]:    7.10, Etot[eV]:  -6571.3979, Epot[eV]:  -6623.6988, Ekin[eV]:   52.3009, T[K]:  299.72, density[g/cm3]: 0.514, calc_time[min]:    60.51
    Dyn  step:     7200, time[ps]:    7.20, Etot[eV]:  -6570.7382, Epot[eV]:  -6623.2257, Ekin[eV]:   52.4875, T[K]:  300.79, density[g/cm3]: 0.515, calc_time[min]:    61.43
    Dyn  step:     7300, time[ps]:    7.30, Etot[eV]:  -6570.9599, Epot[eV]:  -6623.3217, Ekin[eV]:   52.3618, T[K]:  300.07, density[g/cm3]: 0.514, calc_time[min]:    62.32
    Dyn  step:     7400, time[ps]:    7.40, Etot[eV]:  -6570.7702, Epot[eV]:  -6622.8061, Ekin[eV]:   52.0358, T[K]:  298.20, density[g/cm3]: 0.517, calc_time[min]:    63.23
    Dyn  step:     7500, time[ps]:    7.50, Etot[eV]:  -6570.8077, Epot[eV]:  -6622.6438, Ekin[eV]:   51.8361, T[K]:  297.05, density[g/cm3]: 0.520, calc_time[min]:    64.12
    Dyn  step:     7600, time[ps]:    7.60, Etot[eV]:  -6570.6624, Epot[eV]:  -6623.2048, Ekin[eV]:   52.5424, T[K]:  301.10, density[g/cm3]: 0.520, calc_time[min]:    65.03
    Dyn  step:     7700, time[ps]:    7.70, Etot[eV]:  -6570.7481, Epot[eV]:  -6624.4169, Ekin[eV]:   53.6689, T[K]:  307.56, density[g/cm3]: 0.523, calc_time[min]:    65.93
    Dyn  step:     7800, time[ps]:    7.80, Etot[eV]:  -6570.7706, Epot[eV]:  -6622.9235, Ekin[eV]:   52.1529, T[K]:  298.87, density[g/cm3]: 0.526, calc_time[min]:    66.85
    Dyn  step:     7900, time[ps]:    7.90, Etot[eV]:  -6571.5336, Epot[eV]:  -6624.2330, Ekin[eV]:   52.6994, T[K]:  302.00, density[g/cm3]: 0.527, calc_time[min]:    67.74
    Dyn  step:     8000, time[ps]:    8.00, Etot[eV]:  -6571.2334, Epot[eV]:  -6623.8930, Ekin[eV]:   52.6596, T[K]:  301.77, density[g/cm3]: 0.529, calc_time[min]:    68.63
    Dyn  step:     8100, time[ps]:    8.10, Etot[eV]:  -6572.2192, Epot[eV]:  -6623.4319, Ekin[eV]:   51.2127, T[K]:  293.48, density[g/cm3]: 0.531, calc_time[min]:    69.54
    Dyn  step:     8200, time[ps]:    8.20, Etot[eV]:  -6572.1786, Epot[eV]:  -6625.2922, Ekin[eV]:   53.1136, T[K]:  304.37, density[g/cm3]: 0.534, calc_time[min]:    70.44
    Dyn  step:     8300, time[ps]:    8.30, Etot[eV]:  -6572.6281, Epot[eV]:  -6624.5839, Ekin[eV]:   51.9558, T[K]:  297.74, density[g/cm3]: 0.534, calc_time[min]:    71.36
    Dyn  step:     8400, time[ps]:    8.40, Etot[eV]:  -6572.5577, Epot[eV]:  -6626.0274, Ekin[eV]:   53.4697, T[K]:  306.42, density[g/cm3]: 0.533, calc_time[min]:    72.27
    Dyn  step:     8500, time[ps]:    8.50, Etot[eV]:  -6572.2935, Epot[eV]:  -6623.4990, Ekin[eV]:   51.2055, T[K]:  293.44, density[g/cm3]: 0.536, calc_time[min]:    73.16
    Dyn  step:     8600, time[ps]:    8.60, Etot[eV]:  -6572.0408, Epot[eV]:  -6624.8732, Ekin[eV]:   52.8323, T[K]:  302.76, density[g/cm3]: 0.537, calc_time[min]:    74.08
    Dyn  step:     8700, time[ps]:    8.70, Etot[eV]:  -6571.8651, Epot[eV]:  -6624.7193, Ekin[eV]:   52.8543, T[K]:  302.89, density[g/cm3]: 0.538, calc_time[min]:    74.99
    Dyn  step:     8800, time[ps]:    8.80, Etot[eV]:  -6571.5512, Epot[eV]:  -6624.1172, Ekin[eV]:   52.5660, T[K]:  301.24, density[g/cm3]: 0.541, calc_time[min]:    75.90
    Dyn  step:     8900, time[ps]:    8.90, Etot[eV]:  -6572.5844, Epot[eV]:  -6623.7929, Ekin[eV]:   51.2085, T[K]:  293.46, density[g/cm3]: 0.540, calc_time[min]:    76.81
    Dyn  step:     9000, time[ps]:    9.00, Etot[eV]:  -6571.6986, Epot[eV]:  -6625.9040, Ekin[eV]:   54.2054, T[K]:  310.63, density[g/cm3]: 0.542, calc_time[min]:    77.72
    Dyn  step:     9100, time[ps]:    9.10, Etot[eV]:  -6572.2389, Epot[eV]:  -6623.7617, Ekin[eV]:   51.5228, T[K]:  295.26, density[g/cm3]: 0.542, calc_time[min]:    78.63
    Dyn  step:     9200, time[ps]:    9.20, Etot[eV]:  -6572.1325, Epot[eV]:  -6624.6392, Ekin[eV]:   52.5066, T[K]:  300.90, density[g/cm3]: 0.545, calc_time[min]:    79.54
    Dyn  step:     9300, time[ps]:    9.30, Etot[eV]:  -6572.6542, Epot[eV]:  -6624.7297, Ekin[eV]:   52.0756, T[K]:  298.43, density[g/cm3]: 0.544, calc_time[min]:    80.46
    Dyn  step:     9400, time[ps]:    9.40, Etot[eV]:  -6573.1423, Epot[eV]:  -6624.4768, Ekin[eV]:   51.3346, T[K]:  294.18, density[g/cm3]: 0.546, calc_time[min]:    81.35
    Dyn  step:     9500, time[ps]:    9.50, Etot[eV]:  -6573.2527, Epot[eV]:  -6625.5040, Ekin[eV]:   52.2513, T[K]:  299.43, density[g/cm3]: 0.547, calc_time[min]:    82.26
    Dyn  step:     9600, time[ps]:    9.60, Etot[eV]:  -6573.6456, Epot[eV]:  -6625.5554, Ekin[eV]:   51.9098, T[K]:  297.48, density[g/cm3]: 0.548, calc_time[min]:    83.16
    Dyn  step:     9700, time[ps]:    9.70, Etot[eV]:  -6574.4050, Epot[eV]:  -6626.9554, Ekin[eV]:   52.5504, T[K]:  301.15, density[g/cm3]: 0.548, calc_time[min]:    84.08
    Dyn  step:     9800, time[ps]:    9.80, Etot[eV]:  -6573.4274, Epot[eV]:  -6624.8258, Ekin[eV]:   51.3984, T[K]:  294.55, density[g/cm3]: 0.550, calc_time[min]:    84.99
    Dyn  step:     9900, time[ps]:    9.90, Etot[eV]:  -6573.1040, Epot[eV]:  -6623.9247, Ekin[eV]:   50.8206, T[K]:  291.23, density[g/cm3]: 0.550, calc_time[min]:    85.90
    Dyn  step:    10000, time[ps]:   10.00, Etot[eV]:  -6573.0704, Epot[eV]:  -6624.6961, Ekin[eV]:   51.6257, T[K]:  295.85, density[g/cm3]: 0.550, calc_time[min]:    86.81
    Dyn  step:    10100, time[ps]:   10.10, Etot[eV]:  -6572.5037, Epot[eV]:  -6624.3445, Ekin[eV]:   51.8408, T[K]:  297.08, density[g/cm3]: 0.551, calc_time[min]:    87.72
    Dyn  step:    10200, time[ps]:   10.20, Etot[eV]:  -6572.3760, Epot[eV]:  -6624.3323, Ekin[eV]:   51.9563, T[K]:  297.74, density[g/cm3]: 0.552, calc_time[min]:    88.63
    Dyn  step:    10300, time[ps]:   10.30, Etot[eV]:  -6573.2289, Epot[eV]:  -6625.7849, Ekin[eV]:   52.5560, T[K]:  301.18, density[g/cm3]: 0.555, calc_time[min]:    89.52
    Dyn  step:    10400, time[ps]:   10.40, Etot[eV]:  -6573.1877, Epot[eV]:  -6626.1173, Ekin[eV]:   52.9296, T[K]:  303.32, density[g/cm3]: 0.556, calc_time[min]:    90.43
    Dyn  step:    10500, time[ps]:   10.50, Etot[eV]:  -6573.4965, Epot[eV]:  -6625.9309, Ekin[eV]:   52.4344, T[K]:  300.48, density[g/cm3]: 0.556, calc_time[min]:    91.35
    Dyn  step:    10600, time[ps]:   10.60, Etot[eV]:  -6572.2171, Epot[eV]:  -6624.0847, Ekin[eV]:   51.8676, T[K]:  297.23, density[g/cm3]: 0.558, calc_time[min]:    92.26
    Dyn  step:    10700, time[ps]:   10.70, Etot[eV]:  -6572.0825, Epot[eV]:  -6623.5727, Ekin[eV]:   51.4902, T[K]:  295.07, density[g/cm3]: 0.559, calc_time[min]:    93.18
    Dyn  step:    10800, time[ps]:   10.80, Etot[eV]:  -6572.7437, Epot[eV]:  -6625.7472, Ekin[eV]:   53.0035, T[K]:  303.74, density[g/cm3]: 0.559, calc_time[min]:    94.11
    Dyn  step:    10900, time[ps]:   10.90, Etot[eV]:  -6572.5059, Epot[eV]:  -6625.3794, Ekin[eV]:   52.8735, T[K]:  303.00, density[g/cm3]: 0.560, calc_time[min]:    95.03
    Dyn  step:    11000, time[ps]:   11.00, Etot[eV]:  -6573.0784, Epot[eV]:  -6624.8072, Ekin[eV]:   51.7288, T[K]:  296.44, density[g/cm3]: 0.561, calc_time[min]:    95.94
    Dyn  step:    11100, time[ps]:   11.10, Etot[eV]:  -6574.1239, Epot[eV]:  -6625.5868, Ekin[eV]:   51.4629, T[K]:  294.91, density[g/cm3]: 0.563, calc_time[min]:    96.87
    Dyn  step:    11200, time[ps]:   11.20, Etot[eV]:  -6574.0900, Epot[eV]:  -6627.0693, Ekin[eV]:   52.9793, T[K]:  303.60, density[g/cm3]: 0.565, calc_time[min]:    97.78
    Dyn  step:    11300, time[ps]:   11.30, Etot[eV]:  -6572.9262, Epot[eV]:  -6625.2470, Ekin[eV]:   52.3207, T[K]:  299.83, density[g/cm3]: 0.567, calc_time[min]:    98.70
    Dyn  step:    11400, time[ps]:   11.40, Etot[eV]:  -6572.8000, Epot[eV]:  -6624.4576, Ekin[eV]:   51.6576, T[K]:  296.03, density[g/cm3]: 0.569, calc_time[min]:    99.61
    Dyn  step:    11500, time[ps]:   11.50, Etot[eV]:  -6573.2383, Epot[eV]:  -6625.2157, Ekin[eV]:   51.9774, T[K]:  297.86, density[g/cm3]: 0.570, calc_time[min]:   100.52
    Dyn  step:    11600, time[ps]:   11.60, Etot[eV]:  -6573.6503, Epot[eV]:  -6625.4312, Ekin[eV]:   51.7809, T[K]:  296.74, density[g/cm3]: 0.572, calc_time[min]:   101.45
    Dyn  step:    11700, time[ps]:   11.70, Etot[eV]:  -6573.5331, Epot[eV]:  -6625.8460, Ekin[eV]:   52.3129, T[K]:  299.79, density[g/cm3]: 0.575, calc_time[min]:   102.35
    Dyn  step:    11800, time[ps]:   11.80, Etot[eV]:  -6573.8044, Epot[eV]:  -6627.7928, Ekin[eV]:   53.9883, T[K]:  309.39, density[g/cm3]: 0.577, calc_time[min]:   103.27
    Dyn  step:    11900, time[ps]:   11.90, Etot[eV]:  -6573.1002, Epot[eV]:  -6625.6560, Ekin[eV]:   52.5558, T[K]:  301.18, density[g/cm3]: 0.578, calc_time[min]:   104.19
    Dyn  step:    12000, time[ps]:   12.00, Etot[eV]:  -6575.0678, Epot[eV]:  -6627.0458, Ekin[eV]:   51.9780, T[K]:  297.87, density[g/cm3]: 0.579, calc_time[min]:   105.11
    Dyn  step:    12100, time[ps]:   12.10, Etot[eV]:  -6574.4231, Epot[eV]:  -6625.9248, Ekin[eV]:   51.5017, T[K]:  295.14, density[g/cm3]: 0.579, calc_time[min]:   106.03
    Dyn  step:    12200, time[ps]:   12.20, Etot[eV]:  -6574.2641, Epot[eV]:  -6626.0455, Ekin[eV]:   51.7814, T[K]:  296.74, density[g/cm3]: 0.582, calc_time[min]:   106.97
    Dyn  step:    12300, time[ps]:   12.30, Etot[eV]:  -6574.0818, Epot[eV]:  -6626.4198, Ekin[eV]:   52.3380, T[K]:  299.93, density[g/cm3]: 0.583, calc_time[min]:   107.88
    Dyn  step:    12400, time[ps]:   12.40, Etot[eV]:  -6574.5443, Epot[eV]:  -6626.9896, Ekin[eV]:   52.4453, T[K]:  300.54, density[g/cm3]: 0.584, calc_time[min]:   108.80
    Dyn  step:    12500, time[ps]:   12.50, Etot[eV]:  -6574.1214, Epot[eV]:  -6626.6039, Ekin[eV]:   52.4825, T[K]:  300.76, density[g/cm3]: 0.583, calc_time[min]:   109.71
    Dyn  step:    12600, time[ps]:   12.60, Etot[eV]:  -6575.1718, Epot[eV]:  -6628.8916, Ekin[eV]:   53.7197, T[K]:  307.85, density[g/cm3]: 0.584, calc_time[min]:   110.62
    Dyn  step:    12700, time[ps]:   12.70, Etot[eV]:  -6574.9120, Epot[eV]:  -6626.6801, Ekin[eV]:   51.7681, T[K]:  296.66, density[g/cm3]: 0.586, calc_time[min]:   111.52
    Dyn  step:    12800, time[ps]:   12.80, Etot[eV]:  -6573.8169, Epot[eV]:  -6625.3942, Ekin[eV]:   51.5772, T[K]:  295.57, density[g/cm3]: 0.587, calc_time[min]:   112.41
    Dyn  step:    12900, time[ps]:   12.90, Etot[eV]:  -6573.8048, Epot[eV]:  -6626.1908, Ekin[eV]:   52.3861, T[K]:  300.20, density[g/cm3]: 0.591, calc_time[min]:   113.32
    Dyn  step:    13000, time[ps]:   13.00, Etot[eV]:  -6574.6978, Epot[eV]:  -6627.2972, Ekin[eV]:   52.5994, T[K]:  301.43, density[g/cm3]: 0.594, calc_time[min]:   114.25
    Dyn  step:    13100, time[ps]:   13.10, Etot[eV]:  -6574.3640, Epot[eV]:  -6626.9383, Ekin[eV]:   52.5743, T[K]:  301.28, density[g/cm3]: 0.599, calc_time[min]:   115.18
    Dyn  step:    13200, time[ps]:   13.20, Etot[eV]:  -6574.6303, Epot[eV]:  -6626.7914, Ekin[eV]:   52.1611, T[K]:  298.92, density[g/cm3]: 0.603, calc_time[min]:   116.11
    Dyn  step:    13300, time[ps]:   13.30, Etot[eV]:  -6574.1603, Epot[eV]:  -6626.9137, Ekin[eV]:   52.7533, T[K]:  302.31, density[g/cm3]: 0.605, calc_time[min]:   117.05
    Dyn  step:    13400, time[ps]:   13.40, Etot[eV]:  -6574.2736, Epot[eV]:  -6626.7613, Ekin[eV]:   52.4877, T[K]:  300.79, density[g/cm3]: 0.607, calc_time[min]:   117.98
    Dyn  step:    13500, time[ps]:   13.50, Etot[eV]:  -6574.0132, Epot[eV]:  -6624.4997, Ekin[eV]:   50.4864, T[K]:  289.32, density[g/cm3]: 0.608, calc_time[min]:   118.91
    Dyn  step:    13600, time[ps]:   13.60, Etot[eV]:  -6574.0316, Epot[eV]:  -6625.9221, Ekin[eV]:   51.8905, T[K]:  297.37, density[g/cm3]: 0.612, calc_time[min]:   119.85
    Dyn  step:    13700, time[ps]:   13.70, Etot[eV]:  -6573.0959, Epot[eV]:  -6625.4997, Ekin[eV]:   52.4038, T[K]:  300.31, density[g/cm3]: 0.616, calc_time[min]:   120.78
    Dyn  step:    13800, time[ps]:   13.80, Etot[eV]:  -6573.0591, Epot[eV]:  -6625.2917, Ekin[eV]:   52.2325, T[K]:  299.33, density[g/cm3]: 0.621, calc_time[min]:   121.73
    Dyn  step:    13900, time[ps]:   13.90, Etot[eV]:  -6573.7441, Epot[eV]:  -6626.1290, Ekin[eV]:   52.3848, T[K]:  300.20, density[g/cm3]: 0.625, calc_time[min]:   122.68
    Dyn  step:    14000, time[ps]:   14.00, Etot[eV]:  -6573.6327, Epot[eV]:  -6625.7066, Ekin[eV]:   52.0739, T[K]:  298.42, density[g/cm3]: 0.627, calc_time[min]:   123.63
    Dyn  step:    14100, time[ps]:   14.10, Etot[eV]:  -6574.0529, Epot[eV]:  -6625.8470, Ekin[eV]:   51.7941, T[K]:  296.81, density[g/cm3]: 0.631, calc_time[min]:   124.56
    Dyn  step:    14200, time[ps]:   14.20, Etot[eV]:  -6575.1398, Epot[eV]:  -6628.0246, Ekin[eV]:   52.8848, T[K]:  303.06, density[g/cm3]: 0.633, calc_time[min]:   125.51
    Dyn  step:    14300, time[ps]:   14.30, Etot[eV]:  -6575.2519, Epot[eV]:  -6627.4182, Ekin[eV]:   52.1664, T[K]:  298.95, density[g/cm3]: 0.636, calc_time[min]:   126.45
    Dyn  step:    14400, time[ps]:   14.40, Etot[eV]:  -6575.6989, Epot[eV]:  -6627.2438, Ekin[eV]:   51.5449, T[K]:  295.38, density[g/cm3]: 0.635, calc_time[min]:   127.38
    Dyn  step:    14500, time[ps]:   14.50, Etot[eV]:  -6575.2357, Epot[eV]:  -6626.7253, Ekin[eV]:   51.4896, T[K]:  295.07, density[g/cm3]: 0.638, calc_time[min]:   128.34
    Dyn  step:    14600, time[ps]:   14.60, Etot[eV]:  -6575.1414, Epot[eV]:  -6628.8437, Ekin[eV]:   53.7023, T[K]:  307.75, density[g/cm3]: 0.640, calc_time[min]:   129.28
    Dyn  step:    14700, time[ps]:   14.70, Etot[eV]:  -6575.2079, Epot[eV]:  -6627.2341, Ekin[eV]:   52.0262, T[K]:  298.14, density[g/cm3]: 0.639, calc_time[min]:   130.22
    Dyn  step:    14800, time[ps]:   14.80, Etot[eV]:  -6575.3145, Epot[eV]:  -6626.9093, Ekin[eV]:   51.5948, T[K]:  295.67, density[g/cm3]: 0.641, calc_time[min]:   131.18
    Dyn  step:    14900, time[ps]:   14.90, Etot[eV]:  -6574.9082, Epot[eV]:  -6626.1888, Ekin[eV]:   51.2806, T[K]:  293.87, density[g/cm3]: 0.641, calc_time[min]:   132.14
    Dyn  step:    15000, time[ps]:   15.00, Etot[eV]:  -6575.3910, Epot[eV]:  -6628.7827, Ekin[eV]:   53.3917, T[K]:  305.97, density[g/cm3]: 0.645, calc_time[min]:   133.08
    Dyn  step:    15100, time[ps]:   15.10, Etot[eV]:  -6575.9493, Epot[eV]:  -6627.7619, Ekin[eV]:   51.8127, T[K]:  296.92, density[g/cm3]: 0.647, calc_time[min]:   134.02
    Dyn  step:    15200, time[ps]:   15.20, Etot[eV]:  -6576.2952, Epot[eV]:  -6627.9034, Ekin[eV]:   51.6082, T[K]:  295.75, density[g/cm3]: 0.648, calc_time[min]:   134.99
    Dyn  step:    15300, time[ps]:   15.30, Etot[eV]:  -6576.4444, Epot[eV]:  -6629.1257, Ekin[eV]:   52.6813, T[K]:  301.90, density[g/cm3]: 0.649, calc_time[min]:   135.96
    Dyn  step:    15400, time[ps]:   15.40, Etot[eV]:  -6576.5484, Epot[eV]:  -6627.7039, Ekin[eV]:   51.1555, T[K]:  293.15, density[g/cm3]: 0.651, calc_time[min]:   136.91
    Dyn  step:    15500, time[ps]:   15.50, Etot[eV]:  -6576.9837, Epot[eV]:  -6629.2906, Ekin[eV]:   52.3070, T[K]:  299.75, density[g/cm3]: 0.652, calc_time[min]:   137.86
    Dyn  step:    15600, time[ps]:   15.60, Etot[eV]:  -6576.1931, Epot[eV]:  -6628.5129, Ekin[eV]:   52.3198, T[K]:  299.83, density[g/cm3]: 0.654, calc_time[min]:   138.83
    Dyn  step:    15700, time[ps]:   15.70, Etot[eV]:  -6576.5078, Epot[eV]:  -6627.8306, Ekin[eV]:   51.3228, T[K]:  294.11, density[g/cm3]: 0.660, calc_time[min]:   139.79
    Dyn  step:    15800, time[ps]:   15.80, Etot[eV]:  -6576.3195, Epot[eV]:  -6629.2473, Ekin[eV]:   52.9277, T[K]:  303.31, density[g/cm3]: 0.663, calc_time[min]:   140.76
    Dyn  step:    15900, time[ps]:   15.90, Etot[eV]:  -6576.6839, Epot[eV]:  -6627.9796, Ekin[eV]:   51.2957, T[K]:  293.96, density[g/cm3]: 0.662, calc_time[min]:   141.71
    Dyn  step:    16000, time[ps]:   16.00, Etot[eV]:  -6576.1447, Epot[eV]:  -6629.1763, Ekin[eV]:   53.0316, T[K]:  303.90, density[g/cm3]: 0.666, calc_time[min]:   142.67
    Dyn  step:    16100, time[ps]:   16.10, Etot[eV]:  -6575.5922, Epot[eV]:  -6628.0962, Ekin[eV]:   52.5040, T[K]:  300.88, density[g/cm3]: 0.671, calc_time[min]:   143.65
    Dyn  step:    16200, time[ps]:   16.20, Etot[eV]:  -6576.0199, Epot[eV]:  -6627.5246, Ekin[eV]:   51.5047, T[K]:  295.15, density[g/cm3]: 0.674, calc_time[min]:   144.63
    Dyn  step:    16300, time[ps]:   16.30, Etot[eV]:  -6576.0769, Epot[eV]:  -6627.9201, Ekin[eV]:   51.8432, T[K]:  297.09, density[g/cm3]: 0.680, calc_time[min]:   145.59
    Dyn  step:    16400, time[ps]:   16.40, Etot[eV]:  -6576.4711, Epot[eV]:  -6628.2197, Ekin[eV]:   51.7486, T[K]:  296.55, density[g/cm3]: 0.684, calc_time[min]:   146.56
    Dyn  step:    16500, time[ps]:   16.50, Etot[eV]:  -6576.7499, Epot[eV]:  -6629.3836, Ekin[eV]:   52.6337, T[K]:  301.62, density[g/cm3]: 0.691, calc_time[min]:   147.52
    Dyn  step:    16600, time[ps]:   16.60, Etot[eV]:  -6576.3288, Epot[eV]:  -6628.5996, Ekin[eV]:   52.2709, T[K]:  299.54, density[g/cm3]: 0.693, calc_time[min]:   148.49
    Dyn  step:    16700, time[ps]:   16.70, Etot[eV]:  -6576.9663, Epot[eV]:  -6629.2316, Ekin[eV]:   52.2652, T[K]:  299.51, density[g/cm3]: 0.696, calc_time[min]:   149.46
    Dyn  step:    16800, time[ps]:   16.80, Etot[eV]:  -6575.5465, Epot[eV]:  -6625.6063, Ekin[eV]:   50.0598, T[K]:  286.87, density[g/cm3]: 0.695, calc_time[min]:   150.44
    Dyn  step:    16900, time[ps]:   16.90, Etot[eV]:  -6575.0748, Epot[eV]:  -6628.1101, Ekin[eV]:   53.0353, T[K]:  303.93, density[g/cm3]: 0.694, calc_time[min]:   151.43
    Dyn  step:    17000, time[ps]:   17.00, Etot[eV]:  -6576.5248, Epot[eV]:  -6630.5550, Ekin[eV]:   54.0302, T[K]:  309.63, density[g/cm3]: 0.697, calc_time[min]:   152.38
    Dyn  step:    17100, time[ps]:   17.10, Etot[eV]:  -6576.3318, Epot[eV]:  -6627.9143, Ekin[eV]:   51.5824, T[K]:  295.60, density[g/cm3]: 0.700, calc_time[min]:   153.37
    Dyn  step:    17200, time[ps]:   17.20, Etot[eV]:  -6575.8446, Epot[eV]:  -6628.3925, Ekin[eV]:   52.5479, T[K]:  301.13, density[g/cm3]: 0.703, calc_time[min]:   154.35
    Dyn  step:    17300, time[ps]:   17.30, Etot[eV]:  -6577.0708, Epot[eV]:  -6628.8986, Ekin[eV]:   51.8279, T[K]:  297.01, density[g/cm3]: 0.703, calc_time[min]:   155.35
    Dyn  step:    17400, time[ps]:   17.40, Etot[eV]:  -6576.6671, Epot[eV]:  -6628.8976, Ekin[eV]:   52.2305, T[K]:  299.31, density[g/cm3]: 0.709, calc_time[min]:   156.33
    Dyn  step:    17500, time[ps]:   17.50, Etot[eV]:  -6576.4842, Epot[eV]:  -6628.8230, Ekin[eV]:   52.3388, T[K]:  299.93, density[g/cm3]: 0.711, calc_time[min]:   157.33
    Dyn  step:    17600, time[ps]:   17.60, Etot[eV]:  -6576.7765, Epot[eV]:  -6628.0646, Ekin[eV]:   51.2881, T[K]:  293.91, density[g/cm3]: 0.718, calc_time[min]:   158.31
    Dyn  step:    17700, time[ps]:   17.70, Etot[eV]:  -6576.8980, Epot[eV]:  -6628.7925, Ekin[eV]:   51.8944, T[K]:  297.39, density[g/cm3]: 0.716, calc_time[min]:   159.31
    Dyn  step:    17800, time[ps]:   17.80, Etot[eV]:  -6576.6162, Epot[eV]:  -6628.3698, Ekin[eV]:   51.7536, T[K]:  296.58, density[g/cm3]: 0.720, calc_time[min]:   160.31
    Dyn  step:    17900, time[ps]:   17.90, Etot[eV]:  -6576.1432, Epot[eV]:  -6630.9984, Ekin[eV]:   54.8552, T[K]:  314.35, density[g/cm3]: 0.722, calc_time[min]:   161.30
    Dyn  step:    18000, time[ps]:   18.00, Etot[eV]:  -6577.0268, Epot[eV]:  -6628.4153, Ekin[eV]:   51.3885, T[K]:  294.49, density[g/cm3]: 0.722, calc_time[min]:   162.27
    Dyn  step:    18100, time[ps]:   18.10, Etot[eV]:  -6577.0977, Epot[eV]:  -6629.3192, Ekin[eV]:   52.2215, T[K]:  299.26, density[g/cm3]: 0.725, calc_time[min]:   163.26
    Dyn  step:    18200, time[ps]:   18.20, Etot[eV]:  -6577.6654, Epot[eV]:  -6629.6796, Ekin[eV]:   52.0143, T[K]:  298.07, density[g/cm3]: 0.725, calc_time[min]:   164.26
    Dyn  step:    18300, time[ps]:   18.30, Etot[eV]:  -6577.4562, Epot[eV]:  -6630.4972, Ekin[eV]:   53.0410, T[K]:  303.96, density[g/cm3]: 0.728, calc_time[min]:   165.25
    Dyn  step:    18400, time[ps]:   18.40, Etot[eV]:  -6578.1946, Epot[eV]:  -6629.0962, Ekin[eV]:   50.9016, T[K]:  291.70, density[g/cm3]: 0.729, calc_time[min]:   166.24
    Dyn  step:    18500, time[ps]:   18.50, Etot[eV]:  -6576.7943, Epot[eV]:  -6629.1115, Ekin[eV]:   52.3172, T[K]:  299.81, density[g/cm3]: 0.728, calc_time[min]:   167.24
    Dyn  step:    18600, time[ps]:   18.60, Etot[eV]:  -6578.5378, Epot[eV]:  -6630.9530, Ekin[eV]:   52.4152, T[K]:  300.37, density[g/cm3]: 0.730, calc_time[min]:   168.22
    Dyn  step:    18700, time[ps]:   18.70, Etot[eV]:  -6578.6269, Epot[eV]:  -6631.5964, Ekin[eV]:   52.9696, T[K]:  303.55, density[g/cm3]: 0.736, calc_time[min]:   169.21
    Dyn  step:    18800, time[ps]:   18.80, Etot[eV]:  -6578.7092, Epot[eV]:  -6630.9190, Ekin[eV]:   52.2097, T[K]:  299.19, density[g/cm3]: 0.737, calc_time[min]:   170.22
    Dyn  step:    18900, time[ps]:   18.90, Etot[eV]:  -6577.5178, Epot[eV]:  -6630.2388, Ekin[eV]:   52.7210, T[K]:  302.12, density[g/cm3]: 0.737, calc_time[min]:   171.23
    Dyn  step:    19000, time[ps]:   19.00, Etot[eV]:  -6577.8400, Epot[eV]:  -6631.0622, Ekin[eV]:   53.2222, T[K]:  305.00, density[g/cm3]: 0.738, calc_time[min]:   172.24
    Dyn  step:    19100, time[ps]:   19.10, Etot[eV]:  -6578.0805, Epot[eV]:  -6629.6589, Ekin[eV]:   51.5784, T[K]:  295.58, density[g/cm3]: 0.741, calc_time[min]:   173.24
    Dyn  step:    19200, time[ps]:   19.20, Etot[eV]:  -6577.8305, Epot[eV]:  -6631.5914, Ekin[eV]:   53.7609, T[K]:  308.08, density[g/cm3]: 0.746, calc_time[min]:   174.25
    Dyn  step:    19300, time[ps]:   19.30, Etot[eV]:  -6578.4881, Epot[eV]:  -6630.4884, Ekin[eV]:   52.0003, T[K]:  297.99, density[g/cm3]: 0.743, calc_time[min]:   175.24
    Dyn  step:    19400, time[ps]:   19.40, Etot[eV]:  -6577.8073, Epot[eV]:  -6630.2793, Ekin[eV]:   52.4720, T[K]:  300.70, density[g/cm3]: 0.741, calc_time[min]:   176.25
    Dyn  step:    19500, time[ps]:   19.50, Etot[eV]:  -6577.1227, Epot[eV]:  -6629.1006, Ekin[eV]:   51.9779, T[K]:  297.87, density[g/cm3]: 0.737, calc_time[min]:   177.24
    Dyn  step:    19600, time[ps]:   19.60, Etot[eV]:  -6577.5124, Epot[eV]:  -6630.0225, Ekin[eV]:   52.5101, T[K]:  300.92, density[g/cm3]: 0.739, calc_time[min]:   178.25
    Dyn  step:    19700, time[ps]:   19.70, Etot[eV]:  -6577.0900, Epot[eV]:  -6628.4613, Ekin[eV]:   51.3713, T[K]:  294.39, density[g/cm3]: 0.738, calc_time[min]:   179.24
    Dyn  step:    19800, time[ps]:   19.80, Etot[eV]:  -6577.9580, Epot[eV]:  -6631.5799, Ekin[eV]:   53.6219, T[K]:  307.29, density[g/cm3]: 0.742, calc_time[min]:   180.25
    Dyn  step:    19900, time[ps]:   19.90, Etot[eV]:  -6576.8476, Epot[eV]:  -6629.4633, Ekin[eV]:   52.6157, T[K]:  301.52, density[g/cm3]: 0.745, calc_time[min]:   181.26
    Dyn  step:    20000, time[ps]:   20.00, Etot[eV]:  -6577.0552, Epot[eV]:  -6631.1184, Ekin[eV]:   54.0632, T[K]:  309.82, density[g/cm3]: 0.749, calc_time[min]:   182.27
    Dyn  step:    20100, time[ps]:   20.10, Etot[eV]:  -6578.2678, Epot[eV]:  -6630.7770, Ekin[eV]:   52.5092, T[K]:  300.91, density[g/cm3]: 0.750, calc_time[min]:   183.28
    Dyn  step:    20200, time[ps]:   20.20, Etot[eV]:  -6578.2620, Epot[eV]:  -6629.3041, Ekin[eV]:   51.0421, T[K]:  292.50, density[g/cm3]: 0.752, calc_time[min]:   184.29
    Dyn  step:    20300, time[ps]:   20.30, Etot[eV]:  -6577.6416, Epot[eV]:  -6629.3203, Ekin[eV]:   51.6787, T[K]:  296.15, density[g/cm3]: 0.751, calc_time[min]:   185.31
    Dyn  step:    20400, time[ps]:   20.40, Etot[eV]:  -6578.4019, Epot[eV]:  -6629.6597, Ekin[eV]:   51.2578, T[K]:  293.74, density[g/cm3]: 0.753, calc_time[min]:   186.33
    Dyn  step:    20500, time[ps]:   20.50, Etot[eV]:  -6578.3308, Epot[eV]:  -6629.6251, Ekin[eV]:   51.2943, T[K]:  293.95, density[g/cm3]: 0.754, calc_time[min]:   187.35
    Dyn  step:    20600, time[ps]:   20.60, Etot[eV]:  -6578.7170, Epot[eV]:  -6631.5938, Ekin[eV]:   52.8768, T[K]:  303.02, density[g/cm3]: 0.751, calc_time[min]:   188.36
    Dyn  step:    20700, time[ps]:   20.70, Etot[eV]:  -6577.7104, Epot[eV]:  -6629.6688, Ekin[eV]:   51.9583, T[K]:  297.75, density[g/cm3]: 0.751, calc_time[min]:   189.38
    Dyn  step:    20800, time[ps]:   20.80, Etot[eV]:  -6578.3096, Epot[eV]:  -6629.9405, Ekin[eV]:   51.6309, T[K]:  295.88, density[g/cm3]: 0.750, calc_time[min]:   190.42
    Dyn  step:    20900, time[ps]:   20.90, Etot[eV]:  -6577.8518, Epot[eV]:  -6629.7947, Ekin[eV]:   51.9429, T[K]:  297.67, density[g/cm3]: 0.750, calc_time[min]:   191.44
    Dyn  step:    21000, time[ps]:   21.00, Etot[eV]:  -6578.1952, Epot[eV]:  -6630.0308, Ekin[eV]:   51.8356, T[K]:  297.05, density[g/cm3]: 0.751, calc_time[min]:   192.45
    Dyn  step:    21100, time[ps]:   21.10, Etot[eV]:  -6578.0875, Epot[eV]:  -6630.0215, Ekin[eV]:   51.9340, T[K]:  297.61, density[g/cm3]: 0.754, calc_time[min]:   193.46
    Dyn  step:    21200, time[ps]:   21.20, Etot[eV]:  -6578.4030, Epot[eV]:  -6631.1335, Ekin[eV]:   52.7305, T[K]:  302.18, density[g/cm3]: 0.754, calc_time[min]:   194.47
    Dyn  step:    21300, time[ps]:   21.30, Etot[eV]:  -6578.5826, Epot[eV]:  -6629.9496, Ekin[eV]:   51.3670, T[K]:  294.37, density[g/cm3]: 0.752, calc_time[min]:   195.46
    Dyn  step:    21400, time[ps]:   21.40, Etot[eV]:  -6578.4075, Epot[eV]:  -6631.4512, Ekin[eV]:   53.0438, T[K]:  303.97, density[g/cm3]: 0.753, calc_time[min]:   196.45
    Dyn  step:    21500, time[ps]:   21.50, Etot[eV]:  -6578.7582, Epot[eV]:  -6628.9690, Ekin[eV]:   50.2108, T[K]:  287.74, density[g/cm3]: 0.757, calc_time[min]:   197.45
    Dyn  step:    21600, time[ps]:   21.60, Etot[eV]:  -6578.5015, Epot[eV]:  -6629.8282, Ekin[eV]:   51.3268, T[K]:  294.13, density[g/cm3]: 0.756, calc_time[min]:   198.46
    Dyn  step:    21700, time[ps]:   21.70, Etot[eV]:  -6578.2891, Epot[eV]:  -6630.1222, Ekin[eV]:   51.8331, T[K]:  297.04, density[g/cm3]: 0.759, calc_time[min]:   199.48
    Dyn  step:    21800, time[ps]:   21.80, Etot[eV]:  -6577.8818, Epot[eV]:  -6629.7050, Ekin[eV]:   51.8232, T[K]:  296.98, density[g/cm3]: 0.762, calc_time[min]:   200.48
    Dyn  step:    21900, time[ps]:   21.90, Etot[eV]:  -6578.6738, Epot[eV]:  -6632.1546, Ekin[eV]:   53.4808, T[K]:  306.48, density[g/cm3]: 0.764, calc_time[min]:   201.50
    Dyn  step:    22000, time[ps]:   22.00, Etot[eV]:  -6578.4515, Epot[eV]:  -6629.7483, Ekin[eV]:   51.2968, T[K]:  293.96, density[g/cm3]: 0.763, calc_time[min]:   202.51
    Dyn  step:    22100, time[ps]:   22.10, Etot[eV]:  -6578.4766, Epot[eV]:  -6630.5073, Ekin[eV]:   52.0307, T[K]:  298.17, density[g/cm3]: 0.768, calc_time[min]:   203.50
    Dyn  step:    22200, time[ps]:   22.20, Etot[eV]:  -6578.3747, Epot[eV]:  -6632.0662, Ekin[eV]:   53.6915, T[K]:  307.69, density[g/cm3]: 0.770, calc_time[min]:   204.51
    Dyn  step:    22300, time[ps]:   22.30, Etot[eV]:  -6578.5227, Epot[eV]:  -6631.3602, Ekin[eV]:   52.8374, T[K]:  302.79, density[g/cm3]: 0.772, calc_time[min]:   205.51
    Dyn  step:    22400, time[ps]:   22.40, Etot[eV]:  -6579.0932, Epot[eV]:  -6631.2913, Ekin[eV]:   52.1981, T[K]:  299.13, density[g/cm3]: 0.773, calc_time[min]:   206.53
    Dyn  step:    22500, time[ps]:   22.50, Etot[eV]:  -6578.8152, Epot[eV]:  -6630.1410, Ekin[eV]:   51.3259, T[K]:  294.13, density[g/cm3]: 0.774, calc_time[min]:   207.55
    Dyn  step:    22600, time[ps]:   22.60, Etot[eV]:  -6579.6496, Epot[eV]:  -6632.9179, Ekin[eV]:   53.2683, T[K]:  305.26, density[g/cm3]: 0.770, calc_time[min]:   208.58
    Dyn  step:    22700, time[ps]:   22.70, Etot[eV]:  -6580.0798, Epot[eV]:  -6632.3955, Ekin[eV]:   52.3157, T[K]:  299.80, density[g/cm3]: 0.773, calc_time[min]:   209.60
    Dyn  step:    22800, time[ps]:   22.80, Etot[eV]:  -6579.4816, Epot[eV]:  -6632.0249, Ekin[eV]:   52.5434, T[K]:  301.11, density[g/cm3]: 0.771, calc_time[min]:   210.63
    Dyn  step:    22900, time[ps]:   22.90, Etot[eV]:  -6579.6483, Epot[eV]:  -6630.9047, Ekin[eV]:   51.2563, T[K]:  293.73, density[g/cm3]: 0.772, calc_time[min]:   211.64
    Dyn  step:    23000, time[ps]:   23.00, Etot[eV]:  -6579.4165, Epot[eV]:  -6633.1711, Ekin[eV]:   53.7546, T[K]:  308.05, density[g/cm3]: 0.772, calc_time[min]:   212.65
    Dyn  step:    23100, time[ps]:   23.10, Etot[eV]:  -6578.8994, Epot[eV]:  -6630.0870, Ekin[eV]:   51.1876, T[K]:  293.34, density[g/cm3]: 0.774, calc_time[min]:   213.65
    Dyn  step:    23200, time[ps]:   23.20, Etot[eV]:  -6578.7583, Epot[eV]:  -6630.6361, Ekin[eV]:   51.8777, T[K]:  297.29, density[g/cm3]: 0.777, calc_time[min]:   214.66
    Dyn  step:    23300, time[ps]:   23.30, Etot[eV]:  -6578.5421, Epot[eV]:  -6631.0726, Ekin[eV]:   52.5305, T[K]:  301.03, density[g/cm3]: 0.777, calc_time[min]:   215.67
    Dyn  step:    23400, time[ps]:   23.40, Etot[eV]:  -6578.9892, Epot[eV]:  -6631.3780, Ekin[eV]:   52.3888, T[K]:  300.22, density[g/cm3]: 0.777, calc_time[min]:   216.68
    Dyn  step:    23500, time[ps]:   23.50, Etot[eV]:  -6578.9810, Epot[eV]:  -6630.7171, Ekin[eV]:   51.7360, T[K]:  296.48, density[g/cm3]: 0.778, calc_time[min]:   217.68
    Dyn  step:    23600, time[ps]:   23.60, Etot[eV]:  -6579.4189, Epot[eV]:  -6630.3925, Ekin[eV]:   50.9736, T[K]:  292.11, density[g/cm3]: 0.782, calc_time[min]:   218.70
    Dyn  step:    23700, time[ps]:   23.70, Etot[eV]:  -6579.1720, Epot[eV]:  -6632.3009, Ekin[eV]:   53.1289, T[K]:  304.46, density[g/cm3]: 0.780, calc_time[min]:   219.73
    Dyn  step:    23800, time[ps]:   23.80, Etot[eV]:  -6579.8038, Epot[eV]:  -6630.5241, Ekin[eV]:   50.7204, T[K]:  290.66, density[g/cm3]: 0.782, calc_time[min]:   220.75
    Dyn  step:    23900, time[ps]:   23.90, Etot[eV]:  -6579.3671, Epot[eV]:  -6631.3315, Ekin[eV]:   51.9644, T[K]:  297.79, density[g/cm3]: 0.784, calc_time[min]:   221.77
    Dyn  step:    24000, time[ps]:   24.00, Etot[eV]:  -6579.1990, Epot[eV]:  -6630.8857, Ekin[eV]:   51.6867, T[K]:  296.20, density[g/cm3]: 0.786, calc_time[min]:   222.80
    Dyn  step:    24100, time[ps]:   24.10, Etot[eV]:  -6579.0862, Epot[eV]:  -6631.4238, Ekin[eV]:   52.3375, T[K]:  299.93, density[g/cm3]: 0.784, calc_time[min]:   223.84
    Dyn  step:    24200, time[ps]:   24.20, Etot[eV]:  -6578.6180, Epot[eV]:  -6631.3575, Ekin[eV]:   52.7395, T[K]:  302.23, density[g/cm3]: 0.786, calc_time[min]:   224.85
    Dyn  step:    24300, time[ps]:   24.30, Etot[eV]:  -6579.6637, Epot[eV]:  -6631.0412, Ekin[eV]:   51.3776, T[K]:  294.43, density[g/cm3]: 0.784, calc_time[min]:   225.87
    Dyn  step:    24400, time[ps]:   24.40, Etot[eV]:  -6578.8588, Epot[eV]:  -6630.8383, Ekin[eV]:   51.9794, T[K]:  297.87, density[g/cm3]: 0.785, calc_time[min]:   226.89
    Dyn  step:    24500, time[ps]:   24.50, Etot[eV]:  -6579.1119, Epot[eV]:  -6631.8745, Ekin[eV]:   52.7626, T[K]:  302.36, density[g/cm3]: 0.785, calc_time[min]:   227.91
    Dyn  step:    24600, time[ps]:   24.60, Etot[eV]:  -6578.9845, Epot[eV]:  -6629.9706, Ekin[eV]:   50.9862, T[K]:  292.18, density[g/cm3]: 0.784, calc_time[min]:   228.94
    Dyn  step:    24700, time[ps]:   24.70, Etot[eV]:  -6579.0993, Epot[eV]:  -6631.2200, Ekin[eV]:   52.1207, T[K]:  298.68, density[g/cm3]: 0.788, calc_time[min]:   229.97
    Dyn  step:    24800, time[ps]:   24.80, Etot[eV]:  -6579.7497, Epot[eV]:  -6631.3143, Ekin[eV]:   51.5646, T[K]:  295.50, density[g/cm3]: 0.790, calc_time[min]:   230.98
    Dyn  step:    24900, time[ps]:   24.90, Etot[eV]:  -6578.8562, Epot[eV]:  -6632.0004, Ekin[eV]:   53.1442, T[K]:  304.55, density[g/cm3]: 0.791, calc_time[min]:   231.98
    Dyn  step:    25000, time[ps]:   25.00, Etot[eV]:  -6579.3749, Epot[eV]:  -6631.2153, Ekin[eV]:   51.8404, T[K]:  297.08, density[g/cm3]: 0.792, calc_time[min]:   232.99
    Dyn  step:    25100, time[ps]:   25.10, Etot[eV]:  -6578.8719, Epot[eV]:  -6631.5679, Ekin[eV]:   52.6960, T[K]:  301.98, density[g/cm3]: 0.789, calc_time[min]:   234.02
    Dyn  step:    25200, time[ps]:   25.20, Etot[eV]:  -6578.4884, Epot[eV]:  -6631.9872, Ekin[eV]:   53.4989, T[K]:  306.58, density[g/cm3]: 0.789, calc_time[min]:   235.03
    Dyn  step:    25300, time[ps]:   25.30, Etot[eV]:  -6578.9497, Epot[eV]:  -6630.3347, Ekin[eV]:   51.3850, T[K]:  294.47, density[g/cm3]: 0.787, calc_time[min]:   236.07
    Dyn  step:    25400, time[ps]:   25.40, Etot[eV]:  -6579.0460, Epot[eV]:  -6631.7357, Ekin[eV]:   52.6898, T[K]:  301.95, density[g/cm3]: 0.790, calc_time[min]:   237.11
    Dyn  step:    25500, time[ps]:   25.50, Etot[eV]:  -6579.0350, Epot[eV]:  -6629.7292, Ekin[eV]:   50.6942, T[K]:  290.51, density[g/cm3]: 0.794, calc_time[min]:   238.14
    Dyn  step:    25600, time[ps]:   25.60, Etot[eV]:  -6578.1708, Epot[eV]:  -6629.9043, Ekin[eV]:   51.7335, T[K]:  296.47, density[g/cm3]: 0.793, calc_time[min]:   239.17
    Dyn  step:    25700, time[ps]:   25.70, Etot[eV]:  -6578.4068, Epot[eV]:  -6631.1588, Ekin[eV]:   52.7521, T[K]:  302.30, density[g/cm3]: 0.790, calc_time[min]:   240.21
    Dyn  step:    25800, time[ps]:   25.80, Etot[eV]:  -6578.2200, Epot[eV]:  -6630.8593, Ekin[eV]:   52.6392, T[K]:  301.66, density[g/cm3]: 0.786, calc_time[min]:   241.23
    Dyn  step:    25900, time[ps]:   25.90, Etot[eV]:  -6578.4206, Epot[eV]:  -6629.7091, Ekin[eV]:   51.2886, T[K]:  293.92, density[g/cm3]: 0.788, calc_time[min]:   242.25
    Dyn  step:    26000, time[ps]:   26.00, Etot[eV]:  -6578.3750, Epot[eV]:  -6632.0928, Ekin[eV]:   53.7178, T[K]:  307.84, density[g/cm3]: 0.793, calc_time[min]:   243.28
    Dyn  step:    26100, time[ps]:   26.10, Etot[eV]:  -6578.7345, Epot[eV]:  -6631.9082, Ekin[eV]:   53.1738, T[K]:  304.72, density[g/cm3]: 0.793, calc_time[min]:   244.31
    Dyn  step:    26200, time[ps]:   26.20, Etot[eV]:  -6579.0015, Epot[eV]:  -6632.6064, Ekin[eV]:   53.6049, T[K]:  307.19, density[g/cm3]: 0.795, calc_time[min]:   245.34
    Dyn  step:    26300, time[ps]:   26.30, Etot[eV]:  -6579.1136, Epot[eV]:  -6630.5196, Ekin[eV]:   51.4060, T[K]:  294.59, density[g/cm3]: 0.798, calc_time[min]:   246.36
    Dyn  step:    26400, time[ps]:   26.40, Etot[eV]:  -6579.7741, Epot[eV]:  -6633.3179, Ekin[eV]:   53.5439, T[K]:  306.84, density[g/cm3]: 0.799, calc_time[min]:   247.39
    Dyn  step:    26500, time[ps]:   26.50, Etot[eV]:  -6579.3315, Epot[eV]:  -6631.5914, Ekin[eV]:   52.2600, T[K]:  299.48, density[g/cm3]: 0.800, calc_time[min]:   248.40
    Dyn  step:    26600, time[ps]:   26.60, Etot[eV]:  -6579.3469, Epot[eV]:  -6630.4054, Ekin[eV]:   51.0585, T[K]:  292.60, density[g/cm3]: 0.801, calc_time[min]:   249.43
    Dyn  step:    26700, time[ps]:   26.70, Etot[eV]:  -6578.7550, Epot[eV]:  -6629.8295, Ekin[eV]:   51.0746, T[K]:  292.69, density[g/cm3]: 0.803, calc_time[min]:   250.45
    Dyn  step:    26800, time[ps]:   26.80, Etot[eV]:  -6578.3762, Epot[eV]:  -6631.6558, Ekin[eV]:   53.2797, T[K]:  305.33, density[g/cm3]: 0.802, calc_time[min]:   251.47
    Dyn  step:    26900, time[ps]:   26.90, Etot[eV]:  -6578.6636, Epot[eV]:  -6631.1070, Ekin[eV]:   52.4433, T[K]:  300.53, density[g/cm3]: 0.804, calc_time[min]:   252.50
    Dyn  step:    27000, time[ps]:   27.00, Etot[eV]:  -6578.5943, Epot[eV]:  -6631.1121, Ekin[eV]:   52.5177, T[K]:  300.96, density[g/cm3]: 0.801, calc_time[min]:   253.53
    Dyn  step:    27100, time[ps]:   27.10, Etot[eV]:  -6579.1658, Epot[eV]:  -6631.4923, Ekin[eV]:   52.3265, T[K]:  299.86, density[g/cm3]: 0.802, calc_time[min]:   254.55
    Dyn  step:    27200, time[ps]:   27.20, Etot[eV]:  -6578.8851, Epot[eV]:  -6630.5739, Ekin[eV]:   51.6887, T[K]:  296.21, density[g/cm3]: 0.800, calc_time[min]:   255.58
    Dyn  step:    27300, time[ps]:   27.30, Etot[eV]:  -6578.9326, Epot[eV]:  -6630.7360, Ekin[eV]:   51.8034, T[K]:  296.87, density[g/cm3]: 0.804, calc_time[min]:   256.60
    Dyn  step:    27400, time[ps]:   27.40, Etot[eV]:  -6579.0154, Epot[eV]:  -6631.3427, Ekin[eV]:   52.3274, T[K]:  299.87, density[g/cm3]: 0.802, calc_time[min]:   257.62
    Dyn  step:    27500, time[ps]:   27.50, Etot[eV]:  -6579.5073, Epot[eV]:  -6632.2214, Ekin[eV]:   52.7141, T[K]:  302.09, density[g/cm3]: 0.805, calc_time[min]:   258.65
    Dyn  step:    27600, time[ps]:   27.60, Etot[eV]:  -6579.4285, Epot[eV]:  -6629.7350, Ekin[eV]:   50.3065, T[K]:  288.29, density[g/cm3]: 0.806, calc_time[min]:   259.69
    Dyn  step:    27700, time[ps]:   27.70, Etot[eV]:  -6578.9942, Epot[eV]:  -6630.9431, Ekin[eV]:   51.9490, T[K]:  297.70, density[g/cm3]: 0.802, calc_time[min]:   260.72
    Dyn  step:    27800, time[ps]:   27.80, Etot[eV]:  -6579.7575, Epot[eV]:  -6632.6484, Ekin[eV]:   52.8909, T[K]:  303.10, density[g/cm3]: 0.806, calc_time[min]:   261.74
    Dyn  step:    27900, time[ps]:   27.90, Etot[eV]:  -6580.6035, Epot[eV]:  -6632.6947, Ekin[eV]:   52.0912, T[K]:  298.52, density[g/cm3]: 0.807, calc_time[min]:   262.78
    Dyn  step:    28000, time[ps]:   28.00, Etot[eV]:  -6580.0239, Epot[eV]:  -6631.6967, Ekin[eV]:   51.6729, T[K]:  296.12, density[g/cm3]: 0.805, calc_time[min]:   263.82
    Dyn  step:    28100, time[ps]:   28.10, Etot[eV]:  -6579.4175, Epot[eV]:  -6630.3925, Ekin[eV]:   50.9750, T[K]:  292.12, density[g/cm3]: 0.803, calc_time[min]:   264.85
    Dyn  step:    28200, time[ps]:   28.20, Etot[eV]:  -6579.0033, Epot[eV]:  -6631.5624, Ekin[eV]:   52.5592, T[K]:  301.20, density[g/cm3]: 0.801, calc_time[min]:   265.90
    Dyn  step:    28300, time[ps]:   28.30, Etot[eV]:  -6578.2320, Epot[eV]:  -6630.0876, Ekin[eV]:   51.8556, T[K]:  297.17, density[g/cm3]: 0.795, calc_time[min]:   266.90
    Dyn  step:    28400, time[ps]:   28.40, Etot[eV]:  -6578.5488, Epot[eV]:  -6630.1742, Ekin[eV]:   51.6254, T[K]:  295.85, density[g/cm3]: 0.795, calc_time[min]:   267.92
    Dyn  step:    28500, time[ps]:   28.50, Etot[eV]:  -6579.6550, Epot[eV]:  -6631.0300, Ekin[eV]:   51.3750, T[K]:  294.41, density[g/cm3]: 0.796, calc_time[min]:   268.94
    Dyn  step:    28600, time[ps]:   28.60, Etot[eV]:  -6578.5956, Epot[eV]:  -6631.4876, Ekin[eV]:   52.8920, T[K]:  303.10, density[g/cm3]: 0.796, calc_time[min]:   269.94
    Dyn  step:    28700, time[ps]:   28.70, Etot[eV]:  -6578.6490, Epot[eV]:  -6631.9070, Ekin[eV]:   53.2580, T[K]:  305.20, density[g/cm3]: 0.794, calc_time[min]:   270.96
    Dyn  step:    28800, time[ps]:   28.80, Etot[eV]:  -6578.8155, Epot[eV]:  -6630.7451, Ekin[eV]:   51.9297, T[K]:  297.59, density[g/cm3]: 0.793, calc_time[min]:   271.97
    Dyn  step:    28900, time[ps]:   28.90, Etot[eV]:  -6579.0037, Epot[eV]:  -6630.9655, Ekin[eV]:   51.9619, T[K]:  297.77, density[g/cm3]: 0.790, calc_time[min]:   272.99
    Dyn  step:    29000, time[ps]:   29.00, Etot[eV]:  -6578.9247, Epot[eV]:  -6631.4508, Ekin[eV]:   52.5261, T[K]:  301.01, density[g/cm3]: 0.792, calc_time[min]:   274.02
    Dyn  step:    29100, time[ps]:   29.10, Etot[eV]:  -6579.2790, Epot[eV]:  -6631.6822, Ekin[eV]:   52.4032, T[K]:  300.30, density[g/cm3]: 0.793, calc_time[min]:   275.04
    Dyn  step:    29200, time[ps]:   29.20, Etot[eV]:  -6578.7103, Epot[eV]:  -6631.0044, Ekin[eV]:   52.2941, T[K]:  299.68, density[g/cm3]: 0.791, calc_time[min]:   276.07
    Dyn  step:    29300, time[ps]:   29.30, Etot[eV]:  -6578.7760, Epot[eV]:  -6631.0756, Ekin[eV]:   52.2996, T[K]:  299.71, density[g/cm3]: 0.792, calc_time[min]:   277.10
    Dyn  step:    29400, time[ps]:   29.40, Etot[eV]:  -6579.4428, Epot[eV]:  -6631.4546, Ekin[eV]:   52.0118, T[K]:  298.06, density[g/cm3]: 0.792, calc_time[min]:   278.11
    Dyn  step:    29500, time[ps]:   29.50, Etot[eV]:  -6578.7404, Epot[eV]:  -6631.0044, Ekin[eV]:   52.2640, T[K]:  299.51, density[g/cm3]: 0.791, calc_time[min]:   279.13
    Dyn  step:    29600, time[ps]:   29.60, Etot[eV]:  -6579.4966, Epot[eV]:  -6632.6510, Ekin[eV]:   53.1544, T[K]:  304.61, density[g/cm3]: 0.792, calc_time[min]:   280.14
    Dyn  step:    29700, time[ps]:   29.70, Etot[eV]:  -6579.5428, Epot[eV]:  -6630.9088, Ekin[eV]:   51.3661, T[K]:  294.36, density[g/cm3]: 0.792, calc_time[min]:   281.16
    Dyn  step:    29800, time[ps]:   29.80, Etot[eV]:  -6579.2056, Epot[eV]:  -6631.7645, Ekin[eV]:   52.5589, T[K]:  301.20, density[g/cm3]: 0.792, calc_time[min]:   282.17
    Dyn  step:    29900, time[ps]:   29.90, Etot[eV]:  -6579.3892, Epot[eV]:  -6630.8616, Ekin[eV]:   51.4724, T[K]:  294.97, density[g/cm3]: 0.793, calc_time[min]:   283.21
    Dyn  step:    30000, time[ps]:   30.00, Etot[eV]:  -6579.3897, Epot[eV]:  -6631.6175, Ekin[eV]:   52.2277, T[K]:  299.30, density[g/cm3]: 0.790, calc_time[min]:   284.23
    Dyn  step:    30100, time[ps]:   30.10, Etot[eV]:  -6579.3644, Epot[eV]:  -6632.0935, Ekin[eV]:   52.7291, T[K]:  302.17, density[g/cm3]: 0.794, calc_time[min]:   285.25
    Dyn  step:    30200, time[ps]:   30.20, Etot[eV]:  -6579.2710, Epot[eV]:  -6631.0112, Ekin[eV]:   51.7402, T[K]:  296.50, density[g/cm3]: 0.797, calc_time[min]:   286.28
    Dyn  step:    30300, time[ps]:   30.30, Etot[eV]:  -6578.7376, Epot[eV]:  -6632.2407, Ekin[eV]:   53.5031, T[K]:  306.61, density[g/cm3]: 0.798, calc_time[min]:   287.31
    Dyn  step:    30400, time[ps]:   30.40, Etot[eV]:  -6580.2827, Epot[eV]:  -6632.9996, Ekin[eV]:   52.7169, T[K]:  302.10, density[g/cm3]: 0.801, calc_time[min]:   288.35
    Dyn  step:    30500, time[ps]:   30.50, Etot[eV]:  -6579.5868, Epot[eV]:  -6632.7295, Ekin[eV]:   53.1426, T[K]:  304.54, density[g/cm3]: 0.803, calc_time[min]:   289.38
    Dyn  step:    30600, time[ps]:   30.60, Etot[eV]:  -6580.0232, Epot[eV]:  -6632.1126, Ekin[eV]:   52.0894, T[K]:  298.50, density[g/cm3]: 0.806, calc_time[min]:   290.41
    Dyn  step:    30700, time[ps]:   30.70, Etot[eV]:  -6579.4549, Epot[eV]:  -6631.6621, Ekin[eV]:   52.2072, T[K]:  299.18, density[g/cm3]: 0.807, calc_time[min]:   291.45
    Dyn  step:    30800, time[ps]:   30.80, Etot[eV]:  -6579.8424, Epot[eV]:  -6632.7677, Ekin[eV]:   52.9253, T[K]:  303.30, density[g/cm3]: 0.809, calc_time[min]:   292.48
    Dyn  step:    30900, time[ps]:   30.90, Etot[eV]:  -6579.1920, Epot[eV]:  -6632.4346, Ekin[eV]:   53.2426, T[K]:  305.11, density[g/cm3]: 0.808, calc_time[min]:   293.50
    Dyn  step:    31000, time[ps]:   31.00, Etot[eV]:  -6579.1769, Epot[eV]:  -6631.9524, Ekin[eV]:   52.7755, T[K]:  302.44, density[g/cm3]: 0.810, calc_time[min]:   294.53
    Dyn  step:    31100, time[ps]:   31.10, Etot[eV]:  -6579.0452, Epot[eV]:  -6631.9937, Ekin[eV]:   52.9485, T[K]:  303.43, density[g/cm3]: 0.811, calc_time[min]:   295.54
    Dyn  step:    31200, time[ps]:   31.20, Etot[eV]:  -6579.2410, Epot[eV]:  -6632.3957, Ekin[eV]:   53.1547, T[K]:  304.61, density[g/cm3]: 0.809, calc_time[min]:   296.53
    Dyn  step:    31300, time[ps]:   31.30, Etot[eV]:  -6579.3939, Epot[eV]:  -6630.7753, Ekin[eV]:   51.3815, T[K]:  294.45, density[g/cm3]: 0.810, calc_time[min]:   297.53
    Dyn  step:    31400, time[ps]:   31.40, Etot[eV]:  -6579.8279, Epot[eV]:  -6630.3696, Ekin[eV]:   50.5416, T[K]:  289.64, density[g/cm3]: 0.809, calc_time[min]:   298.55
    Dyn  step:    31500, time[ps]:   31.50, Etot[eV]:  -6579.7889, Epot[eV]:  -6632.4218, Ekin[eV]:   52.6329, T[K]:  301.62, density[g/cm3]: 0.808, calc_time[min]:   299.55
    Dyn  step:    31600, time[ps]:   31.60, Etot[eV]:  -6580.1673, Epot[eV]:  -6632.3370, Ekin[eV]:   52.1697, T[K]:  298.97, density[g/cm3]: 0.809, calc_time[min]:   300.56
    Dyn  step:    31700, time[ps]:   31.70, Etot[eV]:  -6580.7940, Epot[eV]:  -6632.5876, Ekin[eV]:   51.7936, T[K]:  296.81, density[g/cm3]: 0.809, calc_time[min]:   301.59
    Dyn  step:    31800, time[ps]:   31.80, Etot[eV]:  -6580.2144, Epot[eV]:  -6631.9891, Ekin[eV]:   51.7747, T[K]:  296.70, density[g/cm3]: 0.808, calc_time[min]:   302.59
    Dyn  step:    31900, time[ps]:   31.90, Etot[eV]:  -6579.4536, Epot[eV]:  -6632.0382, Ekin[eV]:   52.5846, T[K]:  301.34, density[g/cm3]: 0.804, calc_time[min]:   303.60
    Dyn  step:    32000, time[ps]:   32.00, Etot[eV]:  -6579.7070, Epot[eV]:  -6632.9986, Ekin[eV]:   53.2916, T[K]:  305.39, density[g/cm3]: 0.807, calc_time[min]:   304.60
    Dyn  step:    32100, time[ps]:   32.10, Etot[eV]:  -6579.2951, Epot[eV]:  -6630.2946, Ekin[eV]:   50.9995, T[K]:  292.26, density[g/cm3]: 0.809, calc_time[min]:   305.61
    Dyn  step:    32200, time[ps]:   32.20, Etot[eV]:  -6578.7702, Epot[eV]:  -6629.9358, Ekin[eV]:   51.1656, T[K]:  293.21, density[g/cm3]: 0.810, calc_time[min]:   306.65
    Dyn  step:    32300, time[ps]:   32.30, Etot[eV]:  -6578.9641, Epot[eV]:  -6631.9861, Ekin[eV]:   53.0221, T[K]:  303.85, density[g/cm3]: 0.810, calc_time[min]:   307.68
    Dyn  step:    32400, time[ps]:   32.40, Etot[eV]:  -6579.1086, Epot[eV]:  -6631.7786, Ekin[eV]:   52.6700, T[K]:  301.83, density[g/cm3]: 0.809, calc_time[min]:   308.72
    Dyn  step:    32500, time[ps]:   32.50, Etot[eV]:  -6579.3721, Epot[eV]:  -6630.2587, Ekin[eV]:   50.8866, T[K]:  291.61, density[g/cm3]: 0.811, calc_time[min]:   309.74
    Dyn  step:    32600, time[ps]:   32.60, Etot[eV]:  -6579.1481, Epot[eV]:  -6632.0601, Ekin[eV]:   52.9121, T[K]:  303.22, density[g/cm3]: 0.813, calc_time[min]:   310.78
    Dyn  step:    32700, time[ps]:   32.70, Etot[eV]:  -6578.8909, Epot[eV]:  -6632.0475, Ekin[eV]:   53.1566, T[K]:  304.62, density[g/cm3]: 0.815, calc_time[min]:   311.78
    Dyn  step:    32800, time[ps]:   32.80, Etot[eV]:  -6579.4359, Epot[eV]:  -6632.9537, Ekin[eV]:   53.5178, T[K]:  306.69, density[g/cm3]: 0.816, calc_time[min]:   312.77
    Dyn  step:    32900, time[ps]:   32.90, Etot[eV]:  -6579.7413, Epot[eV]:  -6631.4935, Ekin[eV]:   51.7522, T[K]:  296.57, density[g/cm3]: 0.815, calc_time[min]:   313.78
    Dyn  step:    33000, time[ps]:   33.00, Etot[eV]:  -6578.5445, Epot[eV]:  -6630.7144, Ekin[eV]:   52.1698, T[K]:  298.97, density[g/cm3]: 0.812, calc_time[min]:   314.81
    Dyn  step:    33100, time[ps]:   33.10, Etot[eV]:  -6579.2285, Epot[eV]:  -6631.3874, Ekin[eV]:   52.1589, T[K]:  298.90, density[g/cm3]: 0.815, calc_time[min]:   315.85
    Dyn  step:    33200, time[ps]:   33.20, Etot[eV]:  -6579.6805, Epot[eV]:  -6632.2621, Ekin[eV]:   52.5816, T[K]:  301.33, density[g/cm3]: 0.813, calc_time[min]:   316.88
    Dyn  step:    33300, time[ps]:   33.30, Etot[eV]:  -6579.9281, Epot[eV]:  -6631.3832, Ekin[eV]:   51.4551, T[K]:  294.87, density[g/cm3]: 0.818, calc_time[min]:   317.89
    Dyn  step:    33400, time[ps]:   33.40, Etot[eV]:  -6579.4723, Epot[eV]:  -6631.7848, Ekin[eV]:   52.3125, T[K]:  299.78, density[g/cm3]: 0.820, calc_time[min]:   318.93
    Dyn  step:    33500, time[ps]:   33.50, Etot[eV]:  -6579.7901, Epot[eV]:  -6631.2182, Ekin[eV]:   51.4280, T[K]:  294.71, density[g/cm3]: 0.821, calc_time[min]:   319.96
    Dyn  step:    33600, time[ps]:   33.60, Etot[eV]:  -6579.8296, Epot[eV]:  -6632.0935, Ekin[eV]:   52.2639, T[K]:  299.50, density[g/cm3]: 0.824, calc_time[min]:   320.99
    Dyn  step:    33700, time[ps]:   33.70, Etot[eV]:  -6579.5209, Epot[eV]:  -6632.0394, Ekin[eV]:   52.5186, T[K]:  300.96, density[g/cm3]: 0.825, calc_time[min]:   322.03
    Dyn  step:    33800, time[ps]:   33.80, Etot[eV]:  -6579.7699, Epot[eV]:  -6633.2496, Ekin[eV]:   53.4797, T[K]:  306.47, density[g/cm3]: 0.827, calc_time[min]:   323.06
    Dyn  step:    33900, time[ps]:   33.90, Etot[eV]:  -6579.8399, Epot[eV]:  -6629.9285, Ekin[eV]:   50.0886, T[K]:  287.04, density[g/cm3]: 0.826, calc_time[min]:   324.11
    Dyn  step:    34000, time[ps]:   34.00, Etot[eV]:  -6580.1454, Epot[eV]:  -6633.4016, Ekin[eV]:   53.2561, T[K]:  305.19, density[g/cm3]: 0.825, calc_time[min]:   325.15
    Dyn  step:    34100, time[ps]:   34.10, Etot[eV]:  -6580.1788, Epot[eV]:  -6632.9342, Ekin[eV]:   52.7554, T[K]:  302.32, density[g/cm3]: 0.826, calc_time[min]:   326.18
    Dyn  step:    34200, time[ps]:   34.20, Etot[eV]:  -6580.3138, Epot[eV]:  -6631.3223, Ekin[eV]:   51.0085, T[K]:  292.31, density[g/cm3]: 0.828, calc_time[min]:   327.21
    Dyn  step:    34300, time[ps]:   34.30, Etot[eV]:  -6579.8968, Epot[eV]:  -6631.3449, Ekin[eV]:   51.4481, T[K]:  294.83, density[g/cm3]: 0.828, calc_time[min]:   328.23
    Dyn  step:    34400, time[ps]:   34.40, Etot[eV]:  -6579.7509, Epot[eV]:  -6630.0726, Ekin[eV]:   50.3217, T[K]:  288.37, density[g/cm3]: 0.833, calc_time[min]:   329.29
    Dyn  step:    34500, time[ps]:   34.50, Etot[eV]:  -6579.4369, Epot[eV]:  -6632.0246, Ekin[eV]:   52.5877, T[K]:  301.36, density[g/cm3]: 0.836, calc_time[min]:   330.34
    Dyn  step:    34600, time[ps]:   34.60, Etot[eV]:  -6580.1437, Epot[eV]:  -6633.0198, Ekin[eV]:   52.8761, T[K]:  303.01, density[g/cm3]: 0.839, calc_time[min]:   331.39
    Dyn  step:    34700, time[ps]:   34.70, Etot[eV]:  -6579.7666, Epot[eV]:  -6631.4305, Ekin[eV]:   51.6639, T[K]:  296.07, density[g/cm3]: 0.843, calc_time[min]:   332.42
    Dyn  step:    34800, time[ps]:   34.80, Etot[eV]:  -6579.8532, Epot[eV]:  -6633.4277, Ekin[eV]:   53.5745, T[K]:  307.02, density[g/cm3]: 0.843, calc_time[min]:   333.47
    Dyn  step:    34900, time[ps]:   34.90, Etot[eV]:  -6580.4377, Epot[eV]:  -6632.4794, Ekin[eV]:   52.0417, T[K]:  298.23, density[g/cm3]: 0.844, calc_time[min]:   334.51
    Dyn  step:    35000, time[ps]:   35.00, Etot[eV]:  -6579.3303, Epot[eV]:  -6631.2660, Ekin[eV]:   51.9357, T[K]:  297.62, density[g/cm3]: 0.845, calc_time[min]:   335.56
    Dyn  step:    35100, time[ps]:   35.10, Etot[eV]:  -6580.0041, Epot[eV]:  -6632.0708, Ekin[eV]:   52.0666, T[K]:  298.37, density[g/cm3]: 0.851, calc_time[min]:   336.61
    Dyn  step:    35200, time[ps]:   35.20, Etot[eV]:  -6579.5569, Epot[eV]:  -6631.7310, Ekin[eV]:   52.1741, T[K]:  298.99, density[g/cm3]: 0.851, calc_time[min]:   337.66
    Dyn  step:    35300, time[ps]:   35.30, Etot[eV]:  -6579.8311, Epot[eV]:  -6632.5740, Ekin[eV]:   52.7429, T[K]:  302.25, density[g/cm3]: 0.852, calc_time[min]:   338.68
    Dyn  step:    35400, time[ps]:   35.40, Etot[eV]:  -6579.7709, Epot[eV]:  -6632.8102, Ekin[eV]:   53.0393, T[K]:  303.95, density[g/cm3]: 0.857, calc_time[min]:   339.73
    Dyn  step:    35500, time[ps]:   35.50, Etot[eV]:  -6580.0144, Epot[eV]:  -6631.6645, Ekin[eV]:   51.6502, T[K]:  295.99, density[g/cm3]: 0.862, calc_time[min]:   340.77
    Dyn  step:    35600, time[ps]:   35.60, Etot[eV]:  -6579.7872, Epot[eV]:  -6632.2228, Ekin[eV]:   52.4356, T[K]:  300.49, density[g/cm3]: 0.864, calc_time[min]:   341.81
    Dyn  step:    35700, time[ps]:   35.70, Etot[eV]:  -6579.6957, Epot[eV]:  -6632.3777, Ekin[eV]:   52.6820, T[K]:  301.90, density[g/cm3]: 0.869, calc_time[min]:   342.86
    Dyn  step:    35800, time[ps]:   35.80, Etot[eV]:  -6580.2546, Epot[eV]:  -6633.4032, Ekin[eV]:   53.1487, T[K]:  304.58, density[g/cm3]: 0.872, calc_time[min]:   343.92
    Dyn  step:    35900, time[ps]:   35.90, Etot[eV]:  -6580.6635, Epot[eV]:  -6635.5589, Ekin[eV]:   54.8954, T[K]:  314.59, density[g/cm3]: 0.869, calc_time[min]:   344.98
    Dyn  step:    36000, time[ps]:   36.00, Etot[eV]:  -6580.5888, Epot[eV]:  -6632.9030, Ekin[eV]:   52.3143, T[K]:  299.79, density[g/cm3]: 0.871, calc_time[min]:   346.03
    Dyn  step:    36100, time[ps]:   36.10, Etot[eV]:  -6580.5286, Epot[eV]:  -6633.0878, Ekin[eV]:   52.5592, T[K]:  301.20, density[g/cm3]: 0.871, calc_time[min]:   347.09
    Dyn  step:    36200, time[ps]:   36.20, Etot[eV]:  -6580.5066, Epot[eV]:  -6631.5503, Ekin[eV]:   51.0437, T[K]:  292.51, density[g/cm3]: 0.872, calc_time[min]:   348.15
    Dyn  step:    36300, time[ps]:   36.30, Etot[eV]:  -6579.2683, Epot[eV]:  -6631.1582, Ekin[eV]:   51.8899, T[K]:  297.36, density[g/cm3]: 0.874, calc_time[min]:   349.21
    Dyn  step:    36400, time[ps]:   36.40, Etot[eV]:  -6580.1645, Epot[eV]:  -6632.6420, Ekin[eV]:   52.4775, T[K]:  300.73, density[g/cm3]: 0.879, calc_time[min]:   350.25
    Dyn  step:    36500, time[ps]:   36.50, Etot[eV]:  -6580.6069, Epot[eV]:  -6632.8180, Ekin[eV]:   52.2111, T[K]:  299.20, density[g/cm3]: 0.873, calc_time[min]:   351.32
    Dyn  step:    36600, time[ps]:   36.60, Etot[eV]:  -6580.4347, Epot[eV]:  -6631.2885, Ekin[eV]:   50.8539, T[K]:  291.42, density[g/cm3]: 0.872, calc_time[min]:   352.38
    Dyn  step:    36700, time[ps]:   36.70, Etot[eV]:  -6580.0739, Epot[eV]:  -6632.1766, Ekin[eV]:   52.1027, T[K]:  298.58, density[g/cm3]: 0.869, calc_time[min]:   353.44
    Dyn  step:    36800, time[ps]:   36.80, Etot[eV]:  -6580.7865, Epot[eV]:  -6634.2742, Ekin[eV]:   53.4877, T[K]:  306.52, density[g/cm3]: 0.868, calc_time[min]:   354.51
    Dyn  step:    36900, time[ps]:   36.90, Etot[eV]:  -6580.3202, Epot[eV]:  -6632.9507, Ekin[eV]:   52.6305, T[K]:  301.61, density[g/cm3]: 0.870, calc_time[min]:   355.57
    Dyn  step:    37000, time[ps]:   37.00, Etot[eV]:  -6580.4644, Epot[eV]:  -6633.1052, Ekin[eV]:   52.6408, T[K]:  301.66, density[g/cm3]: 0.872, calc_time[min]:   356.62
    Dyn  step:    37100, time[ps]:   37.10, Etot[eV]:  -6580.3997, Epot[eV]:  -6631.9751, Ekin[eV]:   51.5754, T[K]:  295.56, density[g/cm3]: 0.868, calc_time[min]:   357.64
    Dyn  step:    37200, time[ps]:   37.20, Etot[eV]:  -6580.4859, Epot[eV]:  -6632.3786, Ekin[eV]:   51.8928, T[K]:  297.38, density[g/cm3]: 0.870, calc_time[min]:   358.65
    Dyn  step:    37300, time[ps]:   37.30, Etot[eV]:  -6580.4305, Epot[eV]:  -6633.1957, Ekin[eV]:   52.7651, T[K]:  302.38, density[g/cm3]: 0.867, calc_time[min]:   359.67
    Dyn  step:    37400, time[ps]:   37.40, Etot[eV]:  -6580.7209, Epot[eV]:  -6632.8941, Ekin[eV]:   52.1731, T[K]:  298.98, density[g/cm3]: 0.869, calc_time[min]:   360.70
    Dyn  step:    37500, time[ps]:   37.50, Etot[eV]:  -6580.6519, Epot[eV]:  -6634.0682, Ekin[eV]:   53.4162, T[K]:  306.11, density[g/cm3]: 0.869, calc_time[min]:   361.73
    Dyn  step:    37600, time[ps]:   37.60, Etot[eV]:  -6580.8173, Epot[eV]:  -6632.3059, Ekin[eV]:   51.4886, T[K]:  295.06, density[g/cm3]: 0.870, calc_time[min]:   362.76
    Dyn  step:    37700, time[ps]:   37.70, Etot[eV]:  -6580.8747, Epot[eV]:  -6632.2077, Ekin[eV]:   51.3330, T[K]:  294.17, density[g/cm3]: 0.872, calc_time[min]:   363.81
    Dyn  step:    37800, time[ps]:   37.80, Etot[eV]:  -6580.9000, Epot[eV]:  -6633.3309, Ekin[eV]:   52.4309, T[K]:  300.46, density[g/cm3]: 0.878, calc_time[min]:   364.86
    Dyn  step:    37900, time[ps]:   37.90, Etot[eV]:  -6581.2071, Epot[eV]:  -6633.4561, Ekin[eV]:   52.2490, T[K]:  299.42, density[g/cm3]: 0.884, calc_time[min]:   365.90
    Dyn  step:    38000, time[ps]:   38.00, Etot[eV]:  -6580.3246, Epot[eV]:  -6632.5562, Ekin[eV]:   52.2317, T[K]:  299.32, density[g/cm3]: 0.878, calc_time[min]:   366.95
    Dyn  step:    38100, time[ps]:   38.10, Etot[eV]:  -6581.3480, Epot[eV]:  -6634.5018, Ekin[eV]:   53.1538, T[K]:  304.60, density[g/cm3]: 0.880, calc_time[min]:   367.97
    Dyn  step:    38200, time[ps]:   38.20, Etot[eV]:  -6580.9161, Epot[eV]:  -6632.4297, Ekin[eV]:   51.5136, T[K]:  295.21, density[g/cm3]: 0.882, calc_time[min]:   369.02
    Dyn  step:    38300, time[ps]:   38.30, Etot[eV]:  -6580.7888, Epot[eV]:  -6633.6480, Ekin[eV]:   52.8592, T[K]:  302.92, density[g/cm3]: 0.884, calc_time[min]:   370.08
    Dyn  step:    38400, time[ps]:   38.40, Etot[eV]:  -6581.6579, Epot[eV]:  -6634.6717, Ekin[eV]:   53.0138, T[K]:  303.80, density[g/cm3]: 0.886, calc_time[min]:   371.12
    Dyn  step:    38500, time[ps]:   38.50, Etot[eV]:  -6580.7862, Epot[eV]:  -6633.4577, Ekin[eV]:   52.6716, T[K]:  301.84, density[g/cm3]: 0.886, calc_time[min]:   372.17
    Dyn  step:    38600, time[ps]:   38.60, Etot[eV]:  -6580.5015, Epot[eV]:  -6632.2674, Ekin[eV]:   51.7659, T[K]:  296.65, density[g/cm3]: 0.881, calc_time[min]:   373.22
    Dyn  step:    38700, time[ps]:   38.70, Etot[eV]:  -6580.7378, Epot[eV]:  -6632.7809, Ekin[eV]:   52.0431, T[K]:  298.24, density[g/cm3]: 0.881, calc_time[min]:   374.27
    Dyn  step:    38800, time[ps]:   38.80, Etot[eV]:  -6581.3726, Epot[eV]:  -6634.2526, Ekin[eV]:   52.8800, T[K]:  303.04, density[g/cm3]: 0.885, calc_time[min]:   375.33
    Dyn  step:    38900, time[ps]:   38.90, Etot[eV]:  -6581.0356, Epot[eV]:  -6633.7844, Ekin[eV]:   52.7489, T[K]:  302.28, density[g/cm3]: 0.886, calc_time[min]:   376.39
    Dyn  step:    39000, time[ps]:   39.00, Etot[eV]:  -6581.7154, Epot[eV]:  -6633.9180, Ekin[eV]:   52.2026, T[K]:  299.15, density[g/cm3]: 0.885, calc_time[min]:   377.42
    Dyn  step:    39100, time[ps]:   39.10, Etot[eV]:  -6581.2797, Epot[eV]:  -6634.8209, Ekin[eV]:   53.5413, T[K]:  306.83, density[g/cm3]: 0.885, calc_time[min]:   378.50
    Dyn  step:    39200, time[ps]:   39.20, Etot[eV]:  -6580.7369, Epot[eV]:  -6631.7323, Ekin[eV]:   50.9954, T[K]:  292.24, density[g/cm3]: 0.883, calc_time[min]:   379.54
    Dyn  step:    39300, time[ps]:   39.30, Etot[eV]:  -6581.4415, Epot[eV]:  -6634.5466, Ekin[eV]:   53.1051, T[K]:  304.33, density[g/cm3]: 0.887, calc_time[min]:   380.58
    Dyn  step:    39400, time[ps]:   39.40, Etot[eV]:  -6581.4569, Epot[eV]:  -6633.0950, Ekin[eV]:   51.6381, T[K]:  295.92, density[g/cm3]: 0.885, calc_time[min]:   381.63
    Dyn  step:    39500, time[ps]:   39.50, Etot[eV]:  -6581.5142, Epot[eV]:  -6633.7798, Ekin[eV]:   52.2656, T[K]:  299.51, density[g/cm3]: 0.880, calc_time[min]:   382.67
    Dyn  step:    39600, time[ps]:   39.60, Etot[eV]:  -6581.3505, Epot[eV]:  -6634.4428, Ekin[eV]:   53.0924, T[K]:  304.25, density[g/cm3]: 0.883, calc_time[min]:   383.71
    Dyn  step:    39700, time[ps]:   39.70, Etot[eV]:  -6580.9276, Epot[eV]:  -6632.2717, Ekin[eV]:   51.3441, T[K]:  294.23, density[g/cm3]: 0.886, calc_time[min]:   384.78
    Dyn  step:    39800, time[ps]:   39.80, Etot[eV]:  -6581.5959, Epot[eV]:  -6634.7816, Ekin[eV]:   53.1857, T[K]:  304.79, density[g/cm3]: 0.885, calc_time[min]:   385.84
    Dyn  step:    39900, time[ps]:   39.90, Etot[eV]:  -6581.5403, Epot[eV]:  -6632.4930, Ekin[eV]:   50.9526, T[K]:  291.99, density[g/cm3]: 0.886, calc_time[min]:   386.86
    Dyn  step:    40000, time[ps]:   40.00, Etot[eV]:  -6581.4640, Epot[eV]:  -6634.4356, Ekin[eV]:   52.9716, T[K]:  303.56, density[g/cm3]: 0.889, calc_time[min]:   387.91
    Dyn  step:    40100, time[ps]:   40.10, Etot[eV]:  -6582.1149, Epot[eV]:  -6634.9488, Ekin[eV]:   52.8339, T[K]:  302.77, density[g/cm3]: 0.891, calc_time[min]:   388.93
    Dyn  step:    40200, time[ps]:   40.20, Etot[eV]:  -6581.5371, Epot[eV]:  -6633.1439, Ekin[eV]:   51.6068, T[K]:  295.74, density[g/cm3]: 0.896, calc_time[min]:   389.98
    Dyn  step:    40300, time[ps]:   40.30, Etot[eV]:  -6581.0987, Epot[eV]:  -6633.4850, Ekin[eV]:   52.3863, T[K]:  300.21, density[g/cm3]: 0.895, calc_time[min]:   391.04
    Dyn  step:    40400, time[ps]:   40.40, Etot[eV]:  -6581.9386, Epot[eV]:  -6635.4557, Ekin[eV]:   53.5171, T[K]:  306.69, density[g/cm3]: 0.897, calc_time[min]:   392.10
    Dyn  step:    40500, time[ps]:   40.50, Etot[eV]:  -6581.1306, Epot[eV]:  -6633.5893, Ekin[eV]:   52.4587, T[K]:  300.62, density[g/cm3]: 0.899, calc_time[min]:   393.17
    Dyn  step:    40600, time[ps]:   40.60, Etot[eV]:  -6581.1217, Epot[eV]:  -6632.5762, Ekin[eV]:   51.4545, T[K]:  294.87, density[g/cm3]: 0.899, calc_time[min]:   394.22
    Dyn  step:    40700, time[ps]:   40.70, Etot[eV]:  -6581.3739, Epot[eV]:  -6633.4041, Ekin[eV]:   52.0302, T[K]:  298.17, density[g/cm3]: 0.897, calc_time[min]:   395.28
    Dyn  step:    40800, time[ps]:   40.80, Etot[eV]:  -6580.3831, Epot[eV]:  -6632.7553, Ekin[eV]:   52.3723, T[K]:  300.13, density[g/cm3]: 0.898, calc_time[min]:   396.34
    Dyn  step:    40900, time[ps]:   40.90, Etot[eV]:  -6580.1751, Epot[eV]:  -6631.3967, Ekin[eV]:   51.2216, T[K]:  293.53, density[g/cm3]: 0.894, calc_time[min]:   397.43
    Dyn  step:    41000, time[ps]:   41.00, Etot[eV]:  -6581.6958, Epot[eV]:  -6635.2070, Ekin[eV]:   53.5112, T[K]:  306.65, density[g/cm3]: 0.901, calc_time[min]:   398.49
    Dyn  step:    41100, time[ps]:   41.10, Etot[eV]:  -6581.1963, Epot[eV]:  -6634.0412, Ekin[eV]:   52.8449, T[K]:  302.83, density[g/cm3]: 0.902, calc_time[min]:   399.56
    Dyn  step:    41200, time[ps]:   41.20, Etot[eV]:  -6581.4144, Epot[eV]:  -6632.9175, Ekin[eV]:   51.5031, T[K]:  295.14, density[g/cm3]: 0.901, calc_time[min]:   400.64
    Dyn  step:    41300, time[ps]:   41.30, Etot[eV]:  -6580.9172, Epot[eV]:  -6632.7323, Ekin[eV]:   51.8152, T[K]:  296.93, density[g/cm3]: 0.900, calc_time[min]:   401.71
    Dyn  step:    41400, time[ps]:   41.40, Etot[eV]:  -6581.0083, Epot[eV]:  -6635.0128, Ekin[eV]:   54.0044, T[K]:  309.48, density[g/cm3]: 0.900, calc_time[min]:   402.77
    Dyn  step:    41500, time[ps]:   41.50, Etot[eV]:  -6581.5717, Epot[eV]:  -6634.3693, Ekin[eV]:   52.7977, T[K]:  302.56, density[g/cm3]: 0.899, calc_time[min]:   403.84
    Dyn  step:    41600, time[ps]:   41.60, Etot[eV]:  -6581.0782, Epot[eV]:  -6633.7629, Ekin[eV]:   52.6846, T[K]:  301.92, density[g/cm3]: 0.893, calc_time[min]:   404.92
    Dyn  step:    41700, time[ps]:   41.70, Etot[eV]:  -6581.0860, Epot[eV]:  -6632.5736, Ekin[eV]:   51.4876, T[K]:  295.06, density[g/cm3]: 0.893, calc_time[min]:   406.00
    Dyn  step:    41800, time[ps]:   41.80, Etot[eV]:  -6581.1897, Epot[eV]:  -6632.4558, Ekin[eV]:   51.2661, T[K]:  293.79, density[g/cm3]: 0.890, calc_time[min]:   407.07
    Dyn  step:    41900, time[ps]:   41.90, Etot[eV]:  -6580.7291, Epot[eV]:  -6632.0458, Ekin[eV]:   51.3167, T[K]:  294.08, density[g/cm3]: 0.889, calc_time[min]:   408.13
    Dyn  step:    42000, time[ps]:   42.00, Etot[eV]:  -6580.5752, Epot[eV]:  -6632.5959, Ekin[eV]:   52.0207, T[K]:  298.11, density[g/cm3]: 0.893, calc_time[min]:   409.19
    Dyn  step:    42100, time[ps]:   42.10, Etot[eV]:  -6581.6699, Epot[eV]:  -6633.4706, Ekin[eV]:   51.8006, T[K]:  296.85, density[g/cm3]: 0.894, calc_time[min]:   410.24
    Dyn  step:    42200, time[ps]:   42.20, Etot[eV]:  -6580.2375, Epot[eV]:  -6632.8864, Ekin[eV]:   52.6489, T[K]:  301.71, density[g/cm3]: 0.890, calc_time[min]:   411.31
    Dyn  step:    42300, time[ps]:   42.30, Etot[eV]:  -6580.9225, Epot[eV]:  -6633.3149, Ekin[eV]:   52.3924, T[K]:  300.24, density[g/cm3]: 0.888, calc_time[min]:   412.38
    Dyn  step:    42400, time[ps]:   42.40, Etot[eV]:  -6581.4638, Epot[eV]:  -6634.1412, Ekin[eV]:   52.6774, T[K]:  301.87, density[g/cm3]: 0.890, calc_time[min]:   413.46
    Dyn  step:    42500, time[ps]:   42.50, Etot[eV]:  -6581.2087, Epot[eV]:  -6633.5797, Ekin[eV]:   52.3710, T[K]:  300.12, density[g/cm3]: 0.894, calc_time[min]:   414.53
    Dyn  step:    42600, time[ps]:   42.60, Etot[eV]:  -6581.1840, Epot[eV]:  -6633.6208, Ekin[eV]:   52.4368, T[K]:  300.50, density[g/cm3]: 0.890, calc_time[min]:   415.61
    Dyn  step:    42700, time[ps]:   42.70, Etot[eV]:  -6580.9888, Epot[eV]:  -6634.5675, Ekin[eV]:   53.5788, T[K]:  307.04, density[g/cm3]: 0.892, calc_time[min]:   416.68
    Dyn  step:    42800, time[ps]:   42.80, Etot[eV]:  -6581.2415, Epot[eV]:  -6633.0142, Ekin[eV]:   51.7726, T[K]:  296.69, density[g/cm3]: 0.897, calc_time[min]:   417.75
    Dyn  step:    42900, time[ps]:   42.90, Etot[eV]:  -6581.4312, Epot[eV]:  -6632.8889, Ekin[eV]:   51.4577, T[K]:  294.88, density[g/cm3]: 0.901, calc_time[min]:   418.84
    Dyn  step:    43000, time[ps]:   43.00, Etot[eV]:  -6580.7394, Epot[eV]:  -6633.6657, Ekin[eV]:   52.9264, T[K]:  303.30, density[g/cm3]: 0.902, calc_time[min]:   419.90
    Dyn  step:    43100, time[ps]:   43.10, Etot[eV]:  -6581.5115, Epot[eV]:  -6633.7443, Ekin[eV]:   52.2327, T[K]:  299.33, density[g/cm3]: 0.902, calc_time[min]:   420.95
    Dyn  step:    43200, time[ps]:   43.20, Etot[eV]:  -6581.0059, Epot[eV]:  -6632.6309, Ekin[eV]:   51.6251, T[K]:  295.84, density[g/cm3]: 0.901, calc_time[min]:   422.02
    Dyn  step:    43300, time[ps]:   43.30, Etot[eV]:  -6581.4941, Epot[eV]:  -6634.7562, Ekin[eV]:   53.2621, T[K]:  305.23, density[g/cm3]: 0.898, calc_time[min]:   423.06
    Dyn  step:    43400, time[ps]:   43.40, Etot[eV]:  -6582.6694, Epot[eV]:  -6635.7468, Ekin[eV]:   53.0774, T[K]:  304.17, density[g/cm3]: 0.900, calc_time[min]:   424.05
    Dyn  step:    43500, time[ps]:   43.50, Etot[eV]:  -6581.8523, Epot[eV]:  -6633.6211, Ekin[eV]:   51.7688, T[K]:  296.67, density[g/cm3]: 0.901, calc_time[min]:   425.07
    Dyn  step:    43600, time[ps]:   43.60, Etot[eV]:  -6581.2904, Epot[eV]:  -6632.8054, Ekin[eV]:   51.5150, T[K]:  295.21, density[g/cm3]: 0.904, calc_time[min]:   426.10
    Dyn  step:    43700, time[ps]:   43.70, Etot[eV]:  -6581.0401, Epot[eV]:  -6633.3395, Ekin[eV]:   52.2994, T[K]:  299.71, density[g/cm3]: 0.908, calc_time[min]:   427.14
    Dyn  step:    43800, time[ps]:   43.80, Etot[eV]:  -6581.3980, Epot[eV]:  -6634.3022, Ekin[eV]:   52.9042, T[K]:  303.17, density[g/cm3]: 0.912, calc_time[min]:   428.15
    Dyn  step:    43900, time[ps]:   43.90, Etot[eV]:  -6581.6097, Epot[eV]:  -6632.9862, Ekin[eV]:   51.3765, T[K]:  294.42, density[g/cm3]: 0.912, calc_time[min]:   429.14
    Dyn  step:    44000, time[ps]:   44.00, Etot[eV]:  -6580.5577, Epot[eV]:  -6632.7918, Ekin[eV]:   52.2341, T[K]:  299.33, density[g/cm3]: 0.910, calc_time[min]:   430.15
    Dyn  step:    44100, time[ps]:   44.10, Etot[eV]:  -6580.8765, Epot[eV]:  -6632.7939, Ekin[eV]:   51.9173, T[K]:  297.52, density[g/cm3]: 0.919, calc_time[min]:   431.15
    Dyn  step:    44200, time[ps]:   44.20, Etot[eV]:  -6581.6071, Epot[eV]:  -6634.1008, Ekin[eV]:   52.4937, T[K]:  300.82, density[g/cm3]: 0.920, calc_time[min]:   432.19
    Dyn  step:    44300, time[ps]:   44.30, Etot[eV]:  -6581.9272, Epot[eV]:  -6634.4280, Ekin[eV]:   52.5009, T[K]:  300.86, density[g/cm3]: 0.922, calc_time[min]:   433.23
    Dyn  step:    44400, time[ps]:   44.40, Etot[eV]:  -6581.3702, Epot[eV]:  -6634.5743, Ekin[eV]:   53.2040, T[K]:  304.89, density[g/cm3]: 0.922, calc_time[min]:   434.26
    Dyn  step:    44500, time[ps]:   44.50, Etot[eV]:  -6582.1155, Epot[eV]:  -6633.0538, Ekin[eV]:   50.9383, T[K]:  291.91, density[g/cm3]: 0.924, calc_time[min]:   435.30
    Dyn  step:    44600, time[ps]:   44.60, Etot[eV]:  -6581.3880, Epot[eV]:  -6633.9810, Ekin[eV]:   52.5930, T[K]:  301.39, density[g/cm3]: 0.925, calc_time[min]:   436.33
    Dyn  step:    44700, time[ps]:   44.70, Etot[eV]:  -6582.2516, Epot[eV]:  -6633.9432, Ekin[eV]:   51.6916, T[K]:  296.23, density[g/cm3]: 0.924, calc_time[min]:   437.37
    Dyn  step:    44800, time[ps]:   44.80, Etot[eV]:  -6581.3745, Epot[eV]:  -6632.8256, Ekin[eV]:   51.4511, T[K]:  294.85, density[g/cm3]: 0.925, calc_time[min]:   438.40
    Dyn  step:    44900, time[ps]:   44.90, Etot[eV]:  -6581.6552, Epot[eV]:  -6634.2102, Ekin[eV]:   52.5550, T[K]:  301.17, density[g/cm3]: 0.929, calc_time[min]:   439.42
    Dyn  step:    45000, time[ps]:   45.00, Etot[eV]:  -6581.0741, Epot[eV]:  -6633.6524, Ekin[eV]:   52.5783, T[K]:  301.31, density[g/cm3]: 0.929, calc_time[min]:   440.46
    Dyn  step:    45100, time[ps]:   45.10, Etot[eV]:  -6580.1114, Epot[eV]:  -6632.7657, Ekin[eV]:   52.6543, T[K]:  301.74, density[g/cm3]: 0.927, calc_time[min]:   441.47
    Dyn  step:    45200, time[ps]:   45.20, Etot[eV]:  -6580.7163, Epot[eV]:  -6632.2703, Ekin[eV]:   51.5540, T[K]:  295.44, density[g/cm3]: 0.925, calc_time[min]:   442.46
    Dyn  step:    45300, time[ps]:   45.30, Etot[eV]:  -6580.5896, Epot[eV]:  -6633.1258, Ekin[eV]:   52.5362, T[K]:  301.07, density[g/cm3]: 0.928, calc_time[min]:   443.47
    Dyn  step:    45400, time[ps]:   45.40, Etot[eV]:  -6580.6226, Epot[eV]:  -6633.3960, Ekin[eV]:   52.7734, T[K]:  302.42, density[g/cm3]: 0.919, calc_time[min]:   444.45
    Dyn  step:    45500, time[ps]:   45.50, Etot[eV]:  -6581.0265, Epot[eV]:  -6633.1541, Ekin[eV]:   52.1276, T[K]:  298.72, density[g/cm3]: 0.921, calc_time[min]:   445.49
    Dyn  step:    45600, time[ps]:   45.60, Etot[eV]:  -6580.8837, Epot[eV]:  -6632.4244, Ekin[eV]:   51.5407, T[K]:  295.36, density[g/cm3]: 0.922, calc_time[min]:   446.51
    Dyn  step:    45700, time[ps]:   45.70, Etot[eV]:  -6580.4555, Epot[eV]:  -6633.0926, Ekin[eV]:   52.6371, T[K]:  301.64, density[g/cm3]: 0.925, calc_time[min]:   447.50
    Dyn  step:    45800, time[ps]:   45.80, Etot[eV]:  -6581.3963, Epot[eV]:  -6632.1546, Ekin[eV]:   50.7583, T[K]:  290.88, density[g/cm3]: 0.933, calc_time[min]:   448.53
    Dyn  step:    45900, time[ps]:   45.90, Etot[eV]:  -6581.4874, Epot[eV]:  -6632.4827, Ekin[eV]:   50.9953, T[K]:  292.24, density[g/cm3]: 0.934, calc_time[min]:   449.59
    Dyn  step:    46000, time[ps]:   46.00, Etot[eV]:  -6581.6498, Epot[eV]:  -6634.1114, Ekin[eV]:   52.4615, T[K]:  300.64, density[g/cm3]: 0.934, calc_time[min]:   450.68
    Dyn  step:    46100, time[ps]:   46.10, Etot[eV]:  -6581.2074, Epot[eV]:  -6633.8725, Ekin[eV]:   52.6651, T[K]:  301.80, density[g/cm3]: 0.934, calc_time[min]:   451.75
    Dyn  step:    46200, time[ps]:   46.20, Etot[eV]:  -6581.1243, Epot[eV]:  -6634.3119, Ekin[eV]:   53.1876, T[K]:  304.80, density[g/cm3]: 0.934, calc_time[min]:   452.82
    Dyn  step:    46300, time[ps]:   46.30, Etot[eV]:  -6581.2146, Epot[eV]:  -6633.0534, Ekin[eV]:   51.8388, T[K]:  297.07, density[g/cm3]: 0.938, calc_time[min]:   453.92
    Dyn  step:    46400, time[ps]:   46.40, Etot[eV]:  -6580.6423, Epot[eV]:  -6633.1398, Ekin[eV]:   52.4975, T[K]:  300.84, density[g/cm3]: 0.942, calc_time[min]:   455.01
    Dyn  step:    46500, time[ps]:   46.50, Etot[eV]:  -6580.7778, Epot[eV]:  -6633.7249, Ekin[eV]:   52.9471, T[K]:  303.42, density[g/cm3]: 0.944, calc_time[min]:   456.08
    Dyn  step:    46600, time[ps]:   46.60, Etot[eV]:  -6580.9642, Epot[eV]:  -6632.8299, Ekin[eV]:   51.8657, T[K]:  297.22, density[g/cm3]: 0.942, calc_time[min]:   457.18
    Dyn  step:    46700, time[ps]:   46.70, Etot[eV]:  -6581.2966, Epot[eV]:  -6634.2040, Ekin[eV]:   52.9074, T[K]:  303.19, density[g/cm3]: 0.945, calc_time[min]:   458.26
    Dyn  step:    46800, time[ps]:   46.80, Etot[eV]:  -6581.7873, Epot[eV]:  -6634.0372, Ekin[eV]:   52.2499, T[K]:  299.42, density[g/cm3]: 0.948, calc_time[min]:   459.34
    Dyn  step:    46900, time[ps]:   46.90, Etot[eV]:  -6581.5195, Epot[eV]:  -6633.1290, Ekin[eV]:   51.6095, T[K]:  295.75, density[g/cm3]: 0.957, calc_time[min]:   460.43
    Dyn  step:    47000, time[ps]:   47.00, Etot[eV]:  -6582.4187, Epot[eV]:  -6635.2694, Ekin[eV]:   52.8506, T[K]:  302.87, density[g/cm3]: 0.963, calc_time[min]:   461.51
    Dyn  step:    47100, time[ps]:   47.10, Etot[eV]:  -6582.0342, Epot[eV]:  -6634.1169, Ekin[eV]:   52.0828, T[K]:  298.47, density[g/cm3]: 0.963, calc_time[min]:   462.61
    Dyn  step:    47200, time[ps]:   47.20, Etot[eV]:  -6581.7112, Epot[eV]:  -6634.0322, Ekin[eV]:   52.3209, T[K]:  299.83, density[g/cm3]: 0.962, calc_time[min]:   463.69
    Dyn  step:    47300, time[ps]:   47.30, Etot[eV]:  -6580.8244, Epot[eV]:  -6633.8294, Ekin[eV]:   53.0051, T[K]:  303.75, density[g/cm3]: 0.959, calc_time[min]:   464.79
    Dyn  step:    47400, time[ps]:   47.40, Etot[eV]:  -6581.4367, Epot[eV]:  -6633.9765, Ekin[eV]:   52.5398, T[K]:  301.09, density[g/cm3]: 0.961, calc_time[min]:   465.91
    Dyn  step:    47500, time[ps]:   47.50, Etot[eV]:  -6580.5036, Epot[eV]:  -6632.8224, Ekin[eV]:   52.3189, T[K]:  299.82, density[g/cm3]: 0.958, calc_time[min]:   467.03
    Dyn  step:    47600, time[ps]:   47.60, Etot[eV]:  -6581.8604, Epot[eV]:  -6634.4027, Ekin[eV]:   52.5423, T[K]:  301.10, density[g/cm3]: 0.962, calc_time[min]:   468.15
    Dyn  step:    47700, time[ps]:   47.70, Etot[eV]:  -6581.4693, Epot[eV]:  -6633.9715, Ekin[eV]:   52.5022, T[K]:  300.87, density[g/cm3]: 0.958, calc_time[min]:   469.26
    Dyn  step:    47800, time[ps]:   47.80, Etot[eV]:  -6581.9710, Epot[eV]:  -6633.0417, Ekin[eV]:   51.0707, T[K]:  292.67, density[g/cm3]: 0.955, calc_time[min]:   470.37
    Dyn  step:    47900, time[ps]:   47.90, Etot[eV]:  -6582.1589, Epot[eV]:  -6634.7802, Ekin[eV]:   52.6213, T[K]:  301.55, density[g/cm3]: 0.956, calc_time[min]:   471.47
    Dyn  step:    48000, time[ps]:   48.00, Etot[eV]:  -6581.8161, Epot[eV]:  -6634.9070, Ekin[eV]:   53.0909, T[K]:  304.24, density[g/cm3]: 0.954, calc_time[min]:   472.56
    Dyn  step:    48100, time[ps]:   48.10, Etot[eV]:  -6582.0418, Epot[eV]:  -6634.2818, Ekin[eV]:   52.2399, T[K]:  299.37, density[g/cm3]: 0.952, calc_time[min]:   473.65
    Dyn  step:    48200, time[ps]:   48.20, Etot[eV]:  -6581.3182, Epot[eV]:  -6634.0940, Ekin[eV]:   52.7758, T[K]:  302.44, density[g/cm3]: 0.956, calc_time[min]:   474.70
    Dyn  step:    48300, time[ps]:   48.30, Etot[eV]:  -6582.1302, Epot[eV]:  -6634.6556, Ekin[eV]:   52.5254, T[K]:  301.00, density[g/cm3]: 0.951, calc_time[min]:   475.74
    Dyn  step:    48400, time[ps]:   48.40, Etot[eV]:  -6581.5049, Epot[eV]:  -6632.7099, Ekin[eV]:   51.2050, T[K]:  293.44, density[g/cm3]: 0.957, calc_time[min]:   476.78
    Dyn  step:    48500, time[ps]:   48.50, Etot[eV]:  -6582.0879, Epot[eV]:  -6634.6608, Ekin[eV]:   52.5729, T[K]:  301.28, density[g/cm3]: 0.959, calc_time[min]:   477.84
    Dyn  step:    48600, time[ps]:   48.60, Etot[eV]:  -6582.1752, Epot[eV]:  -6633.2702, Ekin[eV]:   51.0950, T[K]:  292.81, density[g/cm3]: 0.958, calc_time[min]:   478.90
    Dyn  step:    48700, time[ps]:   48.70, Etot[eV]:  -6581.8476, Epot[eV]:  -6633.8511, Ekin[eV]:   52.0035, T[K]:  298.01, density[g/cm3]: 0.958, calc_time[min]:   479.97
    Dyn  step:    48800, time[ps]:   48.80, Etot[eV]:  -6581.8808, Epot[eV]:  -6633.7267, Ekin[eV]:   51.8459, T[K]:  297.11, density[g/cm3]: 0.958, calc_time[min]:   481.05
    Dyn  step:    48900, time[ps]:   48.90, Etot[eV]:  -6582.0793, Epot[eV]:  -6634.4062, Ekin[eV]:   52.3269, T[K]:  299.87, density[g/cm3]: 0.958, calc_time[min]:   482.11
    Dyn  step:    49000, time[ps]:   49.00, Etot[eV]:  -6582.2307, Epot[eV]:  -6635.0407, Ekin[eV]:   52.8100, T[K]:  302.63, density[g/cm3]: 0.962, calc_time[min]:   483.16
    Dyn  step:    49100, time[ps]:   49.10, Etot[eV]:  -6581.8022, Epot[eV]:  -6633.5455, Ekin[eV]:   51.7434, T[K]:  296.52, density[g/cm3]: 0.960, calc_time[min]:   484.21
    Dyn  step:    49200, time[ps]:   49.20, Etot[eV]:  -6581.7524, Epot[eV]:  -6632.9328, Ekin[eV]:   51.1804, T[K]:  293.30, density[g/cm3]: 0.956, calc_time[min]:   485.25
    Dyn  step:    49300, time[ps]:   49.30, Etot[eV]:  -6581.7497, Epot[eV]:  -6633.7279, Ekin[eV]:   51.9781, T[K]:  297.87, density[g/cm3]: 0.961, calc_time[min]:   486.34
    Dyn  step:    49400, time[ps]:   49.40, Etot[eV]:  -6581.9175, Epot[eV]:  -6633.9306, Ekin[eV]:   52.0130, T[K]:  298.07, density[g/cm3]: 0.961, calc_time[min]:   487.43
    Dyn  step:    49500, time[ps]:   49.50, Etot[eV]:  -6581.0362, Epot[eV]:  -6632.1538, Ekin[eV]:   51.1176, T[K]:  292.94, density[g/cm3]: 0.961, calc_time[min]:   488.53
    Dyn  step:    49600, time[ps]:   49.60, Etot[eV]:  -6581.4943, Epot[eV]:  -6633.7810, Ekin[eV]:   52.2867, T[K]:  299.64, density[g/cm3]: 0.962, calc_time[min]:   489.62
    Dyn  step:    49700, time[ps]:   49.70, Etot[eV]:  -6582.1855, Epot[eV]:  -6635.0868, Ekin[eV]:   52.9013, T[K]:  303.16, density[g/cm3]: 0.965, calc_time[min]:   490.70
    Dyn  step:    49800, time[ps]:   49.80, Etot[eV]:  -6582.3564, Epot[eV]:  -6633.7117, Ekin[eV]:   51.3553, T[K]:  294.30, density[g/cm3]: 0.966, calc_time[min]:   491.82
    Dyn  step:    49900, time[ps]:   49.90, Etot[eV]:  -6582.5723, Epot[eV]:  -6634.9142, Ekin[eV]:   52.3419, T[K]:  299.95, density[g/cm3]: 0.968, calc_time[min]:   492.92
    Dyn  step:    50000, time[ps]:   50.00, Etot[eV]:  -6581.6758, Epot[eV]:  -6633.8298, Ekin[eV]:   52.1541, T[K]:  298.88, density[g/cm3]: 0.974, calc_time[min]:   493.98
    Dyn  step:    50100, time[ps]:   50.10, Etot[eV]:  -6582.9692, Epot[eV]:  -6635.1116, Ekin[eV]:   52.1424, T[K]:  298.81, density[g/cm3]: 0.984, calc_time[min]:   495.05
    Dyn  step:    50200, time[ps]:   50.20, Etot[eV]:  -6582.6658, Epot[eV]:  -6633.6093, Ekin[eV]:   50.9435, T[K]:  291.94, density[g/cm3]: 0.985, calc_time[min]:   496.10
    Dyn  step:    50300, time[ps]:   50.30, Etot[eV]:  -6582.4595, Epot[eV]:  -6634.9986, Ekin[eV]:   52.5392, T[K]:  301.08, density[g/cm3]: 0.986, calc_time[min]:   497.18
    Dyn  step:    50400, time[ps]:   50.40, Etot[eV]:  -6582.3186, Epot[eV]:  -6634.9816, Ekin[eV]:   52.6630, T[K]:  301.79, density[g/cm3]: 0.992, calc_time[min]:   498.26
    Dyn  step:    50500, time[ps]:   50.50, Etot[eV]:  -6583.1201, Epot[eV]:  -6637.4982, Ekin[eV]:   54.3781, T[K]:  311.62, density[g/cm3]: 0.998, calc_time[min]:   499.37
    Dyn  step:    50600, time[ps]:   50.60, Etot[eV]:  -6582.6389, Epot[eV]:  -6633.5728, Ekin[eV]:   50.9339, T[K]:  291.88, density[g/cm3]: 1.005, calc_time[min]:   500.46
    Dyn  step:    50700, time[ps]:   50.70, Etot[eV]:  -6582.6216, Epot[eV]:  -6634.8302, Ekin[eV]:   52.2086, T[K]:  299.19, density[g/cm3]: 1.001, calc_time[min]:   501.55
    Dyn  step:    50800, time[ps]:   50.80, Etot[eV]:  -6583.0780, Epot[eV]:  -6636.4044, Ekin[eV]:   53.3265, T[K]:  305.59, density[g/cm3]: 1.005, calc_time[min]:   502.65
    Dyn  step:    50900, time[ps]:   50.90, Etot[eV]:  -6582.8661, Epot[eV]:  -6636.0614, Ekin[eV]:   53.1953, T[K]:  304.84, density[g/cm3]: 1.009, calc_time[min]:   503.77
    Dyn  step:    51000, time[ps]:   51.00, Etot[eV]:  -6582.0266, Epot[eV]:  -6633.8964, Ekin[eV]:   51.8698, T[K]:  297.25, density[g/cm3]: 1.010, calc_time[min]:   504.91
    Dyn  step:    51100, time[ps]:   51.10, Etot[eV]:  -6582.0259, Epot[eV]:  -6634.1341, Ekin[eV]:   52.1082, T[K]:  298.61, density[g/cm3]: 1.009, calc_time[min]:   506.03
    Dyn  step:    51200, time[ps]:   51.20, Etot[eV]:  -6582.5213, Epot[eV]:  -6635.0176, Ekin[eV]:   52.4963, T[K]:  300.84, density[g/cm3]: 1.011, calc_time[min]:   507.16
    Dyn  step:    51300, time[ps]:   51.30, Etot[eV]:  -6583.2066, Epot[eV]:  -6634.9896, Ekin[eV]:   51.7830, T[K]:  296.75, density[g/cm3]: 1.012, calc_time[min]:   508.26
    Dyn  step:    51400, time[ps]:   51.40, Etot[eV]:  -6582.7088, Epot[eV]:  -6633.9438, Ekin[eV]:   51.2350, T[K]:  293.61, density[g/cm3]: 1.007, calc_time[min]:   509.38
    Dyn  step:    51500, time[ps]:   51.50, Etot[eV]:  -6583.3250, Epot[eV]:  -6635.3097, Ekin[eV]:   51.9847, T[K]:  297.91, density[g/cm3]: 1.014, calc_time[min]:   510.50
    Dyn  step:    51600, time[ps]:   51.60, Etot[eV]:  -6583.4852, Epot[eV]:  -6635.5848, Ekin[eV]:   52.0995, T[K]:  298.56, density[g/cm3]: 1.016, calc_time[min]:   511.63
    Dyn  step:    51700, time[ps]:   51.70, Etot[eV]:  -6583.3294, Epot[eV]:  -6635.3193, Ekin[eV]:   51.9899, T[K]:  297.93, density[g/cm3]: 1.018, calc_time[min]:   512.73
    Dyn  step:    51800, time[ps]:   51.80, Etot[eV]:  -6582.5437, Epot[eV]:  -6633.2064, Ekin[eV]:   50.6627, T[K]:  290.33, density[g/cm3]: 1.028, calc_time[min]:   513.82
    Dyn  step:    51900, time[ps]:   51.90, Etot[eV]:  -6582.7031, Epot[eV]:  -6636.3781, Ekin[eV]:   53.6750, T[K]:  307.59, density[g/cm3]: 1.031, calc_time[min]:   514.93
    Dyn  step:    52000, time[ps]:   52.00, Etot[eV]:  -6583.3352, Epot[eV]:  -6636.0293, Ekin[eV]:   52.6941, T[K]:  301.97, density[g/cm3]: 1.040, calc_time[min]:   516.03
    Dyn  step:    52100, time[ps]:   52.10, Etot[eV]:  -6583.3131, Epot[eV]:  -6635.7067, Ekin[eV]:   52.3937, T[K]:  300.25, density[g/cm3]: 1.038, calc_time[min]:   517.14
    Dyn  step:    52200, time[ps]:   52.20, Etot[eV]:  -6583.0941, Epot[eV]:  -6634.5989, Ekin[eV]:   51.5048, T[K]:  295.16, density[g/cm3]: 1.034, calc_time[min]:   518.24
    Dyn  step:    52300, time[ps]:   52.30, Etot[eV]:  -6583.1691, Epot[eV]:  -6636.2217, Ekin[eV]:   53.0526, T[K]:  304.02, density[g/cm3]: 1.036, calc_time[min]:   519.31
    Dyn  step:    52400, time[ps]:   52.40, Etot[eV]:  -6582.7315, Epot[eV]:  -6634.7785, Ekin[eV]:   52.0470, T[K]:  298.26, density[g/cm3]: 1.038, calc_time[min]:   520.40
    Dyn  step:    52500, time[ps]:   52.50, Etot[eV]:  -6583.2255, Epot[eV]:  -6636.7994, Ekin[eV]:   53.5739, T[K]:  307.01, density[g/cm3]: 1.035, calc_time[min]:   521.50
    Dyn  step:    52600, time[ps]:   52.60, Etot[eV]:  -6583.2232, Epot[eV]:  -6634.4422, Ekin[eV]:   51.2190, T[K]:  293.52, density[g/cm3]: 1.038, calc_time[min]:   522.59
    Dyn  step:    52700, time[ps]:   52.70, Etot[eV]:  -6583.6792, Epot[eV]:  -6637.6631, Ekin[eV]:   53.9839, T[K]:  309.36, density[g/cm3]: 1.042, calc_time[min]:   523.69
    Dyn  step:    52800, time[ps]:   52.80, Etot[eV]:  -6583.4244, Epot[eV]:  -6635.8056, Ekin[eV]:   52.3812, T[K]:  300.18, density[g/cm3]: 1.040, calc_time[min]:   524.79
    Dyn  step:    52900, time[ps]:   52.90, Etot[eV]:  -6583.5578, Epot[eV]:  -6633.5140, Ekin[eV]:   49.9561, T[K]:  286.28, density[g/cm3]: 1.045, calc_time[min]:   525.90
    Dyn  step:    53000, time[ps]:   53.00, Etot[eV]:  -6583.2864, Epot[eV]:  -6636.3387, Ekin[eV]:   53.0524, T[K]:  304.02, density[g/cm3]: 1.041, calc_time[min]:   527.04
    Dyn  step:    53100, time[ps]:   53.10, Etot[eV]:  -6582.7907, Epot[eV]:  -6635.1783, Ekin[eV]:   52.3876, T[K]:  300.21, density[g/cm3]: 1.042, calc_time[min]:   528.16
    Dyn  step:    53200, time[ps]:   53.20, Etot[eV]:  -6582.8789, Epot[eV]:  -6634.6515, Ekin[eV]:   51.7725, T[K]:  296.69, density[g/cm3]: 1.042, calc_time[min]:   529.29
    Dyn  step:    53300, time[ps]:   53.30, Etot[eV]:  -6583.4790, Epot[eV]:  -6635.9265, Ekin[eV]:   52.4475, T[K]:  300.56, density[g/cm3]: 1.045, calc_time[min]:   530.41
    Dyn  step:    53400, time[ps]:   53.40, Etot[eV]:  -6583.1076, Epot[eV]:  -6635.5991, Ekin[eV]:   52.4914, T[K]:  300.81, density[g/cm3]: 1.049, calc_time[min]:   531.52
    Dyn  step:    53500, time[ps]:   53.50, Etot[eV]:  -6583.0647, Epot[eV]:  -6634.7861, Ekin[eV]:   51.7214, T[K]:  296.40, density[g/cm3]: 1.054, calc_time[min]:   532.62
    Dyn  step:    53600, time[ps]:   53.60, Etot[eV]:  -6582.6779, Epot[eV]:  -6634.7194, Ekin[eV]:   52.0415, T[K]:  298.23, density[g/cm3]: 1.055, calc_time[min]:   533.78
    Dyn  step:    53700, time[ps]:   53.70, Etot[eV]:  -6583.2318, Epot[eV]:  -6636.1790, Ekin[eV]:   52.9471, T[K]:  303.42, density[g/cm3]: 1.055, calc_time[min]:   534.91
    Dyn  step:    53800, time[ps]:   53.80, Etot[eV]:  -6584.0036, Epot[eV]:  -6637.2299, Ekin[eV]:   53.2263, T[K]:  305.02, density[g/cm3]: 1.058, calc_time[min]:   536.04
    Dyn  step:    53900, time[ps]:   53.90, Etot[eV]:  -6583.6509, Epot[eV]:  -6635.6738, Ekin[eV]:   52.0229, T[K]:  298.12, density[g/cm3]: 1.060, calc_time[min]:   537.17
    Dyn  step:    54000, time[ps]:   54.00, Etot[eV]:  -6583.3667, Epot[eV]:  -6634.5090, Ekin[eV]:   51.1423, T[K]:  293.08, density[g/cm3]: 1.062, calc_time[min]:   538.32
    Dyn  step:    54100, time[ps]:   54.10, Etot[eV]:  -6584.4932, Epot[eV]:  -6637.6068, Ekin[eV]:   53.1136, T[K]:  304.37, density[g/cm3]: 1.069, calc_time[min]:   539.46
    Dyn  step:    54200, time[ps]:   54.20, Etot[eV]:  -6584.1529, Epot[eV]:  -6636.3667, Ekin[eV]:   52.2138, T[K]:  299.22, density[g/cm3]: 1.066, calc_time[min]:   540.59
    Dyn  step:    54300, time[ps]:   54.30, Etot[eV]:  -6584.0490, Epot[eV]:  -6636.4082, Ekin[eV]:   52.3592, T[K]:  300.05, density[g/cm3]: 1.069, calc_time[min]:   541.74
    Dyn  step:    54400, time[ps]:   54.40, Etot[eV]:  -6584.6285, Epot[eV]:  -6636.4105, Ekin[eV]:   51.7820, T[K]:  296.74, density[g/cm3]: 1.067, calc_time[min]:   542.88
    Dyn  step:    54500, time[ps]:   54.50, Etot[eV]:  -6584.3469, Epot[eV]:  -6635.8631, Ekin[eV]:   51.5162, T[K]:  295.22, density[g/cm3]: 1.066, calc_time[min]:   544.02
    Dyn  step:    54600, time[ps]:   54.60, Etot[eV]:  -6583.6918, Epot[eV]:  -6636.5557, Ekin[eV]:   52.8638, T[K]:  302.94, density[g/cm3]: 1.063, calc_time[min]:   545.17
    Dyn  step:    54700, time[ps]:   54.70, Etot[eV]:  -6582.8140, Epot[eV]:  -6635.1436, Ekin[eV]:   52.3296, T[K]:  299.88, density[g/cm3]: 1.056, calc_time[min]:   546.30
    Dyn  step:    54800, time[ps]:   54.80, Etot[eV]:  -6583.3709, Epot[eV]:  -6636.0314, Ekin[eV]:   52.6605, T[K]:  301.78, density[g/cm3]: 1.061, calc_time[min]:   547.45
    Dyn  step:    54900, time[ps]:   54.90, Etot[eV]:  -6583.9073, Epot[eV]:  -6636.3880, Ekin[eV]:   52.4807, T[K]:  300.75, density[g/cm3]: 1.060, calc_time[min]:   548.66
    Dyn  step:    55000, time[ps]:   55.00, Etot[eV]:  -6583.7219, Epot[eV]:  -6635.1765, Ekin[eV]:   51.4546, T[K]:  294.87, density[g/cm3]: 1.067, calc_time[min]:   549.87
    Dyn  step:    55100, time[ps]:   55.10, Etot[eV]:  -6583.2855, Epot[eV]:  -6635.5588, Ekin[eV]:   52.2733, T[K]:  299.56, density[g/cm3]: 1.070, calc_time[min]:   551.06
    Dyn  step:    55200, time[ps]:   55.20, Etot[eV]:  -6583.6083, Epot[eV]:  -6635.3469, Ekin[eV]:   51.7385, T[K]:  296.49, density[g/cm3]: 1.073, calc_time[min]:   552.21
    Dyn  step:    55300, time[ps]:   55.30, Etot[eV]:  -6583.5466, Epot[eV]:  -6634.7455, Ekin[eV]:   51.1989, T[K]:  293.40, density[g/cm3]: 1.071, calc_time[min]:   553.34
    Dyn  step:    55400, time[ps]:   55.40, Etot[eV]:  -6584.1429, Epot[eV]:  -6636.4763, Ekin[eV]:   52.3334, T[K]:  299.90, density[g/cm3]: 1.074, calc_time[min]:   554.45
    Dyn  step:    55500, time[ps]:   55.50, Etot[eV]:  -6583.8084, Epot[eV]:  -6633.9840, Ekin[eV]:   50.1756, T[K]:  287.54, density[g/cm3]: 1.077, calc_time[min]:   555.54
    Dyn  step:    55600, time[ps]:   55.60, Etot[eV]:  -6583.4775, Epot[eV]:  -6635.5050, Ekin[eV]:   52.0276, T[K]:  298.15, density[g/cm3]: 1.077, calc_time[min]:   556.60
    Dyn  step:    55700, time[ps]:   55.70, Etot[eV]:  -6583.8703, Epot[eV]:  -6635.3906, Ekin[eV]:   51.5203, T[K]:  295.24, density[g/cm3]: 1.083, calc_time[min]:   557.67
    Dyn  step:    55800, time[ps]:   55.80, Etot[eV]:  -6583.7991, Epot[eV]:  -6635.4666, Ekin[eV]:   51.6675, T[K]:  296.09, density[g/cm3]: 1.091, calc_time[min]:   558.73
    Dyn  step:    55900, time[ps]:   55.90, Etot[eV]:  -6584.7649, Epot[eV]:  -6637.1991, Ekin[eV]:   52.4342, T[K]:  300.48, density[g/cm3]: 1.092, calc_time[min]:   559.75
    Dyn  step:    56000, time[ps]:   56.00, Etot[eV]:  -6584.1346, Epot[eV]:  -6634.8008, Ekin[eV]:   50.6662, T[K]:  290.35, density[g/cm3]: 1.092, calc_time[min]:   560.74
    Dyn  step:    56100, time[ps]:   56.10, Etot[eV]:  -6584.4020, Epot[eV]:  -6636.8133, Ekin[eV]:   52.4114, T[K]:  300.35, density[g/cm3]: 1.087, calc_time[min]:   561.82
    Dyn  step:    56200, time[ps]:   56.20, Etot[eV]:  -6584.1596, Epot[eV]:  -6636.4679, Ekin[eV]:   52.3083, T[K]:  299.76, density[g/cm3]: 1.093, calc_time[min]:   562.89
    Dyn  step:    56300, time[ps]:   56.30, Etot[eV]:  -6584.9947, Epot[eV]:  -6638.0298, Ekin[eV]:   53.0350, T[K]:  303.92, density[g/cm3]: 1.090, calc_time[min]:   563.96
    Dyn  step:    56400, time[ps]:   56.40, Etot[eV]:  -6584.8683, Epot[eV]:  -6636.9065, Ekin[eV]:   52.0382, T[K]:  298.21, density[g/cm3]: 1.097, calc_time[min]:   565.04
    Dyn  step:    56500, time[ps]:   56.50, Etot[eV]:  -6584.4614, Epot[eV]:  -6635.5932, Ekin[eV]:   51.1318, T[K]:  293.02, density[g/cm3]: 1.094, calc_time[min]:   566.12
    Dyn  step:    56600, time[ps]:   56.60, Etot[eV]:  -6584.9195, Epot[eV]:  -6636.8951, Ekin[eV]:   51.9755, T[K]:  297.85, density[g/cm3]: 1.095, calc_time[min]:   567.18
    Dyn  step:    56700, time[ps]:   56.70, Etot[eV]:  -6585.0922, Epot[eV]:  -6636.8476, Ekin[eV]:   51.7554, T[K]:  296.59, density[g/cm3]: 1.093, calc_time[min]:   568.28
    Dyn  step:    56800, time[ps]:   56.80, Etot[eV]:  -6584.3378, Epot[eV]:  -6636.7677, Ekin[eV]:   52.4298, T[K]:  300.46, density[g/cm3]: 1.095, calc_time[min]:   569.39
    Dyn  step:    56900, time[ps]:   56.90, Etot[eV]:  -6585.0942, Epot[eV]:  -6637.8392, Ekin[eV]:   52.7451, T[K]:  302.26, density[g/cm3]: 1.090, calc_time[min]:   570.48
    Dyn  step:    57000, time[ps]:   57.00, Etot[eV]:  -6585.4852, Epot[eV]:  -6637.7218, Ekin[eV]:   52.2366, T[K]:  299.35, density[g/cm3]: 1.086, calc_time[min]:   571.54
    Dyn  step:    57100, time[ps]:   57.10, Etot[eV]:  -6584.9676, Epot[eV]:  -6638.0402, Ekin[eV]:   53.0726, T[K]:  304.14, density[g/cm3]: 1.082, calc_time[min]:   572.63
    Dyn  step:    57200, time[ps]:   57.20, Etot[eV]:  -6584.6910, Epot[eV]:  -6636.6173, Ekin[eV]:   51.9263, T[K]:  297.57, density[g/cm3]: 1.082, calc_time[min]:   573.72
    Dyn  step:    57300, time[ps]:   57.30, Etot[eV]:  -6584.5374, Epot[eV]:  -6637.1606, Ekin[eV]:   52.6232, T[K]:  301.56, density[g/cm3]: 1.081, calc_time[min]:   574.84
    Dyn  step:    57400, time[ps]:   57.40, Etot[eV]:  -6583.9039, Epot[eV]:  -6635.6622, Ekin[eV]:   51.7583, T[K]:  296.61, density[g/cm3]: 1.084, calc_time[min]:   575.96
    Dyn  step:    57500, time[ps]:   57.50, Etot[eV]:  -6584.3872, Epot[eV]:  -6637.5139, Ekin[eV]:   53.1267, T[K]:  304.45, density[g/cm3]: 1.081, calc_time[min]:   577.07
    Dyn  step:    57600, time[ps]:   57.60, Etot[eV]:  -6584.4462, Epot[eV]:  -6637.0847, Ekin[eV]:   52.6386, T[K]:  301.65, density[g/cm3]: 1.084, calc_time[min]:   578.19
    Dyn  step:    57700, time[ps]:   57.70, Etot[eV]:  -6584.6076, Epot[eV]:  -6636.3719, Ekin[eV]:   51.7643, T[K]:  296.64, density[g/cm3]: 1.087, calc_time[min]:   579.30
    Dyn  step:    57800, time[ps]:   57.80, Etot[eV]:  -6584.7983, Epot[eV]:  -6638.4215, Ekin[eV]:   53.6233, T[K]:  307.30, density[g/cm3]: 1.092, calc_time[min]:   580.40
    Dyn  step:    57900, time[ps]:   57.90, Etot[eV]:  -6585.0037, Epot[eV]:  -6635.8572, Ekin[eV]:   50.8536, T[K]:  291.42, density[g/cm3]: 1.094, calc_time[min]:   581.54
    Dyn  step:    58000, time[ps]:   58.00, Etot[eV]:  -6583.5588, Epot[eV]:  -6635.6423, Ekin[eV]:   52.0834, T[K]:  298.47, density[g/cm3]: 1.088, calc_time[min]:   582.65
    Dyn  step:    58100, time[ps]:   58.10, Etot[eV]:  -6584.3316, Epot[eV]:  -6636.9102, Ekin[eV]:   52.5785, T[K]:  301.31, density[g/cm3]: 1.088, calc_time[min]:   583.75
    Dyn  step:    58200, time[ps]:   58.20, Etot[eV]:  -6585.2746, Epot[eV]:  -6638.3351, Ekin[eV]:   53.0605, T[K]:  304.07, density[g/cm3]: 1.087, calc_time[min]:   584.86
    Dyn  step:    58300, time[ps]:   58.30, Etot[eV]:  -6584.3862, Epot[eV]:  -6636.3263, Ekin[eV]:   51.9401, T[K]:  297.65, density[g/cm3]: 1.081, calc_time[min]:   585.94
    Dyn  step:    58400, time[ps]:   58.40, Etot[eV]:  -6585.2702, Epot[eV]:  -6637.4185, Ekin[eV]:   52.1483, T[K]:  298.84, density[g/cm3]: 1.084, calc_time[min]:   587.00
    Dyn  step:    58500, time[ps]:   58.50, Etot[eV]:  -6585.3694, Epot[eV]:  -6638.7735, Ekin[eV]:   53.4041, T[K]:  306.04, density[g/cm3]: 1.087, calc_time[min]:   588.05
    Dyn  step:    58600, time[ps]:   58.60, Etot[eV]:  -6584.4942, Epot[eV]:  -6636.2567, Ekin[eV]:   51.7625, T[K]:  296.63, density[g/cm3]: 1.087, calc_time[min]:   589.10
    Dyn  step:    58700, time[ps]:   58.70, Etot[eV]:  -6584.6848, Epot[eV]:  -6637.0951, Ekin[eV]:   52.4103, T[K]:  300.34, density[g/cm3]: 1.091, calc_time[min]:   590.17
    Dyn  step:    58800, time[ps]:   58.80, Etot[eV]:  -6584.1527, Epot[eV]:  -6636.6117, Ekin[eV]:   52.4589, T[K]:  300.62, density[g/cm3]: 1.090, calc_time[min]:   591.23
    Dyn  step:    58900, time[ps]:   58.90, Etot[eV]:  -6585.3600, Epot[eV]:  -6638.5468, Ekin[eV]:   53.1868, T[K]:  304.79, density[g/cm3]: 1.092, calc_time[min]:   592.29
    Dyn  step:    59000, time[ps]:   59.00, Etot[eV]:  -6584.5617, Epot[eV]:  -6635.3957, Ekin[eV]:   50.8340, T[K]:  291.31, density[g/cm3]: 1.097, calc_time[min]:   593.37
    Dyn  step:    59100, time[ps]:   59.10, Etot[eV]:  -6584.0137, Epot[eV]:  -6636.3353, Ekin[eV]:   52.3216, T[K]:  299.84, density[g/cm3]: 1.098, calc_time[min]:   594.43
    Dyn  step:    59200, time[ps]:   59.20, Etot[eV]:  -6584.4624, Epot[eV]:  -6637.2979, Ekin[eV]:   52.8355, T[K]:  302.78, density[g/cm3]: 1.093, calc_time[min]:   595.49
    Dyn  step:    59300, time[ps]:   59.30, Etot[eV]:  -6585.0123, Epot[eV]:  -6637.2168, Ekin[eV]:   52.2044, T[K]:  299.16, density[g/cm3]: 1.102, calc_time[min]:   596.58
    Dyn  step:    59400, time[ps]:   59.40, Etot[eV]:  -6584.8231, Epot[eV]:  -6636.4300, Ekin[eV]:   51.6069, T[K]:  295.74, density[g/cm3]: 1.099, calc_time[min]:   597.66
    Dyn  step:    59500, time[ps]:   59.50, Etot[eV]:  -6584.5902, Epot[eV]:  -6637.0373, Ekin[eV]:   52.4472, T[K]:  300.56, density[g/cm3]: 1.089, calc_time[min]:   598.74
    Dyn  step:    59600, time[ps]:   59.60, Etot[eV]:  -6584.4407, Epot[eV]:  -6637.7076, Ekin[eV]:   53.2668, T[K]:  305.25, density[g/cm3]: 1.088, calc_time[min]:   599.84
    Dyn  step:    59700, time[ps]:   59.70, Etot[eV]:  -6584.9487, Epot[eV]:  -6637.1886, Ekin[eV]:   52.2399, T[K]:  299.37, density[g/cm3]: 1.090, calc_time[min]:   600.93
    Dyn  step:    59800, time[ps]:   59.80, Etot[eV]:  -6584.1885, Epot[eV]:  -6637.6802, Ekin[eV]:   53.4917, T[K]:  306.54, density[g/cm3]: 1.093, calc_time[min]:   602.03
    Dyn  step:    59900, time[ps]:   59.90, Etot[eV]:  -6584.5636, Epot[eV]:  -6636.0209, Ekin[eV]:   51.4573, T[K]:  294.88, density[g/cm3]: 1.097, calc_time[min]:   603.10
    Dyn  step:    60000, time[ps]:   60.00, Etot[eV]:  -6584.8776, Epot[eV]:  -6636.6514, Ekin[eV]:   51.7738, T[K]:  296.70, density[g/cm3]: 1.092, calc_time[min]:   604.20
    Dyn  step:    60100, time[ps]:   60.10, Etot[eV]:  -6584.8724, Epot[eV]:  -6638.1963, Ekin[eV]:   53.3239, T[K]:  305.58, density[g/cm3]: 1.092, calc_time[min]:   605.32
    Dyn  step:    60200, time[ps]:   60.20, Etot[eV]:  -6585.0776, Epot[eV]:  -6636.7940, Ekin[eV]:   51.7164, T[K]:  296.37, density[g/cm3]: 1.090, calc_time[min]:   606.44
    Dyn  step:    60300, time[ps]:   60.30, Etot[eV]:  -6584.3976, Epot[eV]:  -6636.8947, Ekin[eV]:   52.4971, T[K]:  300.84, density[g/cm3]: 1.081, calc_time[min]:   607.55
    Dyn  step:    60400, time[ps]:   60.40, Etot[eV]:  -6585.0497, Epot[eV]:  -6637.5168, Ekin[eV]:   52.4671, T[K]:  300.67, density[g/cm3]: 1.084, calc_time[min]:   608.63
    Dyn  step:    60500, time[ps]:   60.50, Etot[eV]:  -6584.9148, Epot[eV]:  -6636.6088, Ekin[eV]:   51.6940, T[K]:  296.24, density[g/cm3]: 1.090, calc_time[min]:   609.69
    Dyn  step:    60600, time[ps]:   60.60, Etot[eV]:  -6584.9618, Epot[eV]:  -6636.7524, Ekin[eV]:   51.7906, T[K]:  296.79, density[g/cm3]: 1.090, calc_time[min]:   610.77
    Dyn  step:    60700, time[ps]:   60.70, Etot[eV]:  -6584.8535, Epot[eV]:  -6637.4769, Ekin[eV]:   52.6234, T[K]:  301.57, density[g/cm3]: 1.092, calc_time[min]:   611.83
    Dyn  step:    60800, time[ps]:   60.80, Etot[eV]:  -6584.7885, Epot[eV]:  -6638.0678, Ekin[eV]:   53.2793, T[K]:  305.32, density[g/cm3]: 1.089, calc_time[min]:   612.93
    Dyn  step:    60900, time[ps]:   60.90, Etot[eV]:  -6585.3933, Epot[eV]:  -6637.6628, Ekin[eV]:   52.2694, T[K]:  299.54, density[g/cm3]: 1.090, calc_time[min]:   614.00
    Dyn  step:    61000, time[ps]:   61.00, Etot[eV]:  -6585.4481, Epot[eV]:  -6637.7580, Ekin[eV]:   52.3099, T[K]:  299.77, density[g/cm3]: 1.084, calc_time[min]:   615.05
    Dyn  step:    61100, time[ps]:   61.10, Etot[eV]:  -6585.1562, Epot[eV]:  -6636.9639, Ekin[eV]:   51.8077, T[K]:  296.89, density[g/cm3]: 1.079, calc_time[min]:   616.11
    Dyn  step:    61200, time[ps]:   61.20, Etot[eV]:  -6584.7400, Epot[eV]:  -6636.5588, Ekin[eV]:   51.8189, T[K]:  296.95, density[g/cm3]: 1.077, calc_time[min]:   617.13
    Dyn  step:    61300, time[ps]:   61.30, Etot[eV]:  -6584.8724, Epot[eV]:  -6637.5143, Ekin[eV]:   52.6419, T[K]:  301.67, density[g/cm3]: 1.081, calc_time[min]:   618.14
    Dyn  step:    61400, time[ps]:   61.40, Etot[eV]:  -6585.0507, Epot[eV]:  -6638.1578, Ekin[eV]:   53.1071, T[K]:  304.34, density[g/cm3]: 1.080, calc_time[min]:   619.17
    Dyn  step:    61500, time[ps]:   61.50, Etot[eV]:  -6585.6881, Epot[eV]:  -6636.3653, Ekin[eV]:   50.6772, T[K]:  290.41, density[g/cm3]: 1.083, calc_time[min]:   620.21
    Dyn  step:    61600, time[ps]:   61.60, Etot[eV]:  -6585.1591, Epot[eV]:  -6637.8366, Ekin[eV]:   52.6775, T[K]:  301.87, density[g/cm3]: 1.082, calc_time[min]:   621.22
    Dyn  step:    61700, time[ps]:   61.70, Etot[eV]:  -6583.8116, Epot[eV]:  -6636.0201, Ekin[eV]:   52.2085, T[K]:  299.19, density[g/cm3]: 1.079, calc_time[min]:   622.22
    Dyn  step:    61800, time[ps]:   61.80, Etot[eV]:  -6584.7975, Epot[eV]:  -6636.4294, Ekin[eV]:   51.6320, T[K]:  295.88, density[g/cm3]: 1.079, calc_time[min]:   623.26
    Dyn  step:    61900, time[ps]:   61.90, Etot[eV]:  -6584.9220, Epot[eV]:  -6637.7764, Ekin[eV]:   52.8544, T[K]:  302.89, density[g/cm3]: 1.083, calc_time[min]:   624.28
    Dyn  step:    62000, time[ps]:   62.00, Etot[eV]:  -6584.5959, Epot[eV]:  -6636.9973, Ekin[eV]:   52.4014, T[K]:  300.29, density[g/cm3]: 1.083, calc_time[min]:   625.31
    Dyn  step:    62100, time[ps]:   62.10, Etot[eV]:  -6584.5224, Epot[eV]:  -6636.3332, Ekin[eV]:   51.8108, T[K]:  296.91, density[g/cm3]: 1.075, calc_time[min]:   626.35
    Dyn  step:    62200, time[ps]:   62.20, Etot[eV]:  -6585.5567, Epot[eV]:  -6637.1217, Ekin[eV]:   51.5650, T[K]:  295.50, density[g/cm3]: 1.080, calc_time[min]:   627.35
    Dyn  step:    62300, time[ps]:   62.30, Etot[eV]:  -6585.5731, Epot[eV]:  -6637.6811, Ekin[eV]:   52.1080, T[K]:  298.61, density[g/cm3]: 1.078, calc_time[min]:   628.37
    Dyn  step:    62400, time[ps]:   62.40, Etot[eV]:  -6585.4756, Epot[eV]:  -6637.6287, Ekin[eV]:   52.1530, T[K]:  298.87, density[g/cm3]: 1.082, calc_time[min]:   629.41
    Dyn  step:    62500, time[ps]:   62.50, Etot[eV]:  -6585.4466, Epot[eV]:  -6637.3628, Ekin[eV]:   51.9162, T[K]:  297.51, density[g/cm3]: 1.082, calc_time[min]:   630.41
    Dyn  step:    62600, time[ps]:   62.60, Etot[eV]:  -6584.5136, Epot[eV]:  -6635.3882, Ekin[eV]:   50.8746, T[K]:  291.54, density[g/cm3]: 1.082, calc_time[min]:   631.45
    Dyn  step:    62700, time[ps]:   62.70, Etot[eV]:  -6585.2480, Epot[eV]:  -6637.1392, Ekin[eV]:   51.8912, T[K]:  297.37, density[g/cm3]: 1.079, calc_time[min]:   632.48
    Dyn  step:    62800, time[ps]:   62.80, Etot[eV]:  -6585.8211, Epot[eV]:  -6638.7837, Ekin[eV]:   52.9626, T[K]:  303.51, density[g/cm3]: 1.083, calc_time[min]:   633.50
    Dyn  step:    62900, time[ps]:   62.90, Etot[eV]:  -6585.3219, Epot[eV]:  -6637.0108, Ekin[eV]:   51.6889, T[K]:  296.21, density[g/cm3]: 1.078, calc_time[min]:   634.55
    Dyn  step:    63000, time[ps]:   63.00, Etot[eV]:  -6584.9938, Epot[eV]:  -6636.6127, Ekin[eV]:   51.6188, T[K]:  295.81, density[g/cm3]: 1.075, calc_time[min]:   635.60
    Dyn  step:    63100, time[ps]:   63.10, Etot[eV]:  -6584.7018, Epot[eV]:  -6636.7323, Ekin[eV]:   52.0305, T[K]:  298.17, density[g/cm3]: 1.077, calc_time[min]:   636.62
    Dyn  step:    63200, time[ps]:   63.20, Etot[eV]:  -6584.9090, Epot[eV]:  -6637.7393, Ekin[eV]:   52.8303, T[K]:  302.75, density[g/cm3]: 1.081, calc_time[min]:   637.67
    Dyn  step:    63300, time[ps]:   63.30, Etot[eV]:  -6585.2745, Epot[eV]:  -6637.3075, Ekin[eV]:   52.0330, T[K]:  298.18, density[g/cm3]: 1.080, calc_time[min]:   638.71
    Dyn  step:    63400, time[ps]:   63.40, Etot[eV]:  -6585.5390, Epot[eV]:  -6638.7625, Ekin[eV]:   53.2236, T[K]:  305.00, density[g/cm3]: 1.081, calc_time[min]:   639.77
    Dyn  step:    63500, time[ps]:   63.50, Etot[eV]:  -6585.4450, Epot[eV]:  -6636.8718, Ekin[eV]:   51.4267, T[K]:  294.71, density[g/cm3]: 1.088, calc_time[min]:   640.85
    Dyn  step:    63600, time[ps]:   63.60, Etot[eV]:  -6585.8147, Epot[eV]:  -6639.4435, Ekin[eV]:   53.6288, T[K]:  307.33, density[g/cm3]: 1.087, calc_time[min]:   641.89
    Dyn  step:    63700, time[ps]:   63.70, Etot[eV]:  -6586.0600, Epot[eV]:  -6636.2272, Ekin[eV]:   50.1672, T[K]:  287.49, density[g/cm3]: 1.088, calc_time[min]:   642.95
    Dyn  step:    63800, time[ps]:   63.80, Etot[eV]:  -6585.4477, Epot[eV]:  -6637.8511, Ekin[eV]:   52.4034, T[K]:  300.30, density[g/cm3]: 1.094, calc_time[min]:   644.05
    Dyn  step:    63900, time[ps]:   63.90, Etot[eV]:  -6584.6984, Epot[eV]:  -6636.7420, Ekin[eV]:   52.0437, T[K]:  298.24, density[g/cm3]: 1.093, calc_time[min]:   645.16
    Dyn  step:    64000, time[ps]:   64.00, Etot[eV]:  -6584.7876, Epot[eV]:  -6637.4686, Ekin[eV]:   52.6810, T[K]:  301.90, density[g/cm3]: 1.095, calc_time[min]:   646.24
    Dyn  step:    64100, time[ps]:   64.10, Etot[eV]:  -6584.2018, Epot[eV]:  -6636.1884, Ekin[eV]:   51.9866, T[K]:  297.92, density[g/cm3]: 1.095, calc_time[min]:   647.36
    Dyn  step:    64200, time[ps]:   64.20, Etot[eV]:  -6585.3117, Epot[eV]:  -6638.6633, Ekin[eV]:   53.3516, T[K]:  305.74, density[g/cm3]: 1.092, calc_time[min]:   648.45
    Dyn  step:    64300, time[ps]:   64.30, Etot[eV]:  -6585.1085, Epot[eV]:  -6636.6575, Ekin[eV]:   51.5490, T[K]:  295.41, density[g/cm3]: 1.100, calc_time[min]:   649.55
    Dyn  step:    64400, time[ps]:   64.40, Etot[eV]:  -6585.0576, Epot[eV]:  -6637.7376, Ekin[eV]:   52.6800, T[K]:  301.89, density[g/cm3]: 1.102, calc_time[min]:   650.63
    Dyn  step:    64500, time[ps]:   64.50, Etot[eV]:  -6584.5893, Epot[eV]:  -6636.1141, Ekin[eV]:   51.5249, T[K]:  295.27, density[g/cm3]: 1.105, calc_time[min]:   651.73
    Dyn  step:    64600, time[ps]:   64.60, Etot[eV]:  -6584.9785, Epot[eV]:  -6637.7908, Ekin[eV]:   52.8124, T[K]:  302.65, density[g/cm3]: 1.104, calc_time[min]:   652.81
    Dyn  step:    64700, time[ps]:   64.70, Etot[eV]:  -6585.1723, Epot[eV]:  -6636.0077, Ekin[eV]:   50.8354, T[K]:  291.32, density[g/cm3]: 1.107, calc_time[min]:   653.89
    Dyn  step:    64800, time[ps]:   64.80, Etot[eV]:  -6585.5916, Epot[eV]:  -6638.6841, Ekin[eV]:   53.0925, T[K]:  304.25, density[g/cm3]: 1.107, calc_time[min]:   655.06
    Dyn  step:    64900, time[ps]:   64.90, Etot[eV]:  -6585.5432, Epot[eV]:  -6637.7864, Ekin[eV]:   52.2432, T[K]:  299.39, density[g/cm3]: 1.107, calc_time[min]:   656.20
    Dyn  step:    65000, time[ps]:   65.00, Etot[eV]:  -6584.6434, Epot[eV]:  -6636.3578, Ekin[eV]:   51.7144, T[K]:  296.36, density[g/cm3]: 1.102, calc_time[min]:   657.30
    Dyn  step:    65100, time[ps]:   65.10, Etot[eV]:  -6585.4254, Epot[eV]:  -6638.6576, Ekin[eV]:   53.2322, T[K]:  305.05, density[g/cm3]: 1.108, calc_time[min]:   658.39
    Dyn  step:    65200, time[ps]:   65.20, Etot[eV]:  -6585.1973, Epot[eV]:  -6637.2764, Ekin[eV]:   52.0791, T[K]:  298.45, density[g/cm3]: 1.106, calc_time[min]:   659.50
    Dyn  step:    65300, time[ps]:   65.30, Etot[eV]:  -6585.0672, Epot[eV]:  -6637.8223, Ekin[eV]:   52.7551, T[K]:  302.32, density[g/cm3]: 1.105, calc_time[min]:   660.58
    Dyn  step:    65400, time[ps]:   65.40, Etot[eV]:  -6584.6086, Epot[eV]:  -6636.8258, Ekin[eV]:   52.2173, T[K]:  299.24, density[g/cm3]: 1.106, calc_time[min]:   661.64
    Dyn  step:    65500, time[ps]:   65.50, Etot[eV]:  -6585.0071, Epot[eV]:  -6636.8546, Ekin[eV]:   51.8475, T[K]:  297.12, density[g/cm3]: 1.110, calc_time[min]:   662.71
    Dyn  step:    65600, time[ps]:   65.60, Etot[eV]:  -6585.7664, Epot[eV]:  -6638.0403, Ekin[eV]:   52.2739, T[K]:  299.56, density[g/cm3]: 1.109, calc_time[min]:   663.77
    Dyn  step:    65700, time[ps]:   65.70, Etot[eV]:  -6584.7149, Epot[eV]:  -6636.7485, Ekin[eV]:   52.0337, T[K]:  298.19, density[g/cm3]: 1.108, calc_time[min]:   664.81
    Dyn  step:    65800, time[ps]:   65.80, Etot[eV]:  -6584.3284, Epot[eV]:  -6637.6583, Ekin[eV]:   53.3299, T[K]:  305.61, density[g/cm3]: 1.102, calc_time[min]:   665.87
    Dyn  step:    65900, time[ps]:   65.90, Etot[eV]:  -6584.8153, Epot[eV]:  -6636.3231, Ekin[eV]:   51.5078, T[K]:  295.17, density[g/cm3]: 1.102, calc_time[min]:   666.86
    Dyn  step:    66000, time[ps]:   66.00, Etot[eV]:  -6584.8326, Epot[eV]:  -6636.8587, Ekin[eV]:   52.0261, T[K]:  298.14, density[g/cm3]: 1.102, calc_time[min]:   667.89
    Dyn  step:    66100, time[ps]:   66.10, Etot[eV]:  -6585.2210, Epot[eV]:  -6637.6353, Ekin[eV]:   52.4143, T[K]:  300.37, density[g/cm3]: 1.101, calc_time[min]:   668.92
    Dyn  step:    66200, time[ps]:   66.20, Etot[eV]:  -6585.1196, Epot[eV]:  -6638.0707, Ekin[eV]:   52.9512, T[K]:  303.44, density[g/cm3]: 1.100, calc_time[min]:   669.93
    Dyn  step:    66300, time[ps]:   66.30, Etot[eV]:  -6584.5692, Epot[eV]:  -6635.3062, Ekin[eV]:   50.7371, T[K]:  290.76, density[g/cm3]: 1.102, calc_time[min]:   670.95
    Dyn  step:    66400, time[ps]:   66.40, Etot[eV]:  -6585.1576, Epot[eV]:  -6638.2028, Ekin[eV]:   53.0451, T[K]:  303.98, density[g/cm3]: 1.100, calc_time[min]:   671.95
    Dyn  step:    66500, time[ps]:   66.50, Etot[eV]:  -6584.1461, Epot[eV]:  -6635.7351, Ekin[eV]:   51.5890, T[K]:  295.64, density[g/cm3]: 1.093, calc_time[min]:   672.97
    Dyn  step:    66600, time[ps]:   66.60, Etot[eV]:  -6585.0344, Epot[eV]:  -6637.3376, Ekin[eV]:   52.3032, T[K]:  299.73, density[g/cm3]: 1.093, calc_time[min]:   673.98
    Dyn  step:    66700, time[ps]:   66.70, Etot[eV]:  -6584.2542, Epot[eV]:  -6636.8838, Ekin[eV]:   52.6296, T[K]:  301.60, density[g/cm3]: 1.091, calc_time[min]:   674.97
    Dyn  step:    66800, time[ps]:   66.80, Etot[eV]:  -6584.2773, Epot[eV]:  -6636.3640, Ekin[eV]:   52.0867, T[K]:  298.49, density[g/cm3]: 1.094, calc_time[min]:   675.96
    Dyn  step:    66900, time[ps]:   66.90, Etot[eV]:  -6584.7166, Epot[eV]:  -6637.1904, Ekin[eV]:   52.4738, T[K]:  300.71, density[g/cm3]: 1.093, calc_time[min]:   676.97
    Dyn  step:    67000, time[ps]:   67.00, Etot[eV]:  -6585.7914, Epot[eV]:  -6638.7295, Ekin[eV]:   52.9381, T[K]:  303.37, density[g/cm3]: 1.093, calc_time[min]:   677.99
    Dyn  step:    67100, time[ps]:   67.10, Etot[eV]:  -6585.6550, Epot[eV]:  -6637.0639, Ekin[eV]:   51.4089, T[K]:  294.61, density[g/cm3]: 1.094, calc_time[min]:   678.97
    Dyn  step:    67200, time[ps]:   67.20, Etot[eV]:  -6585.9849, Epot[eV]:  -6638.8383, Ekin[eV]:   52.8534, T[K]:  302.88, density[g/cm3]: 1.097, calc_time[min]:   679.96
    Dyn  step:    67300, time[ps]:   67.30, Etot[eV]:  -6585.3206, Epot[eV]:  -6638.1663, Ekin[eV]:   52.8456, T[K]:  302.84, density[g/cm3]: 1.101, calc_time[min]:   680.98
    Dyn  step:    67400, time[ps]:   67.40, Etot[eV]:  -6585.1579, Epot[eV]:  -6636.8736, Ekin[eV]:   51.7157, T[K]:  296.36, density[g/cm3]: 1.103, calc_time[min]:   682.00
    Dyn  step:    67500, time[ps]:   67.50, Etot[eV]:  -6585.1082, Epot[eV]:  -6637.8184, Ekin[eV]:   52.7102, T[K]:  302.06, density[g/cm3]: 1.103, calc_time[min]:   683.02
    Dyn  step:    67600, time[ps]:   67.60, Etot[eV]:  -6585.1303, Epot[eV]:  -6636.0646, Ekin[eV]:   50.9343, T[K]:  291.89, density[g/cm3]: 1.100, calc_time[min]:   684.05
    Dyn  step:    67700, time[ps]:   67.70, Etot[eV]:  -6585.5161, Epot[eV]:  -6637.8040, Ekin[eV]:   52.2878, T[K]:  299.64, density[g/cm3]: 1.095, calc_time[min]:   685.07
    Dyn  step:    67800, time[ps]:   67.80, Etot[eV]:  -6586.3079, Epot[eV]:  -6639.0576, Ekin[eV]:   52.7497, T[K]:  302.29, density[g/cm3]: 1.097, calc_time[min]:   686.07
    Dyn  step:    67900, time[ps]:   67.90, Etot[eV]:  -6585.2662, Epot[eV]:  -6638.1896, Ekin[eV]:   52.9234, T[K]:  303.28, density[g/cm3]: 1.098, calc_time[min]:   687.08
    Dyn  step:    68000, time[ps]:   68.00, Etot[eV]:  -6585.2484, Epot[eV]:  -6636.8985, Ekin[eV]:   51.6501, T[K]:  295.99, density[g/cm3]: 1.094, calc_time[min]:   688.13
    Dyn  step:    68100, time[ps]:   68.10, Etot[eV]:  -6586.2316, Epot[eV]:  -6638.7689, Ekin[eV]:   52.5374, T[K]:  301.07, density[g/cm3]: 1.093, calc_time[min]:   689.19
    Dyn  step:    68200, time[ps]:   68.20, Etot[eV]:  -6585.8114, Epot[eV]:  -6638.1592, Ekin[eV]:   52.3477, T[K]:  299.99, density[g/cm3]: 1.099, calc_time[min]:   690.23
    Dyn  step:    68300, time[ps]:   68.30, Etot[eV]:  -6584.8657, Epot[eV]:  -6636.5372, Ekin[eV]:   51.6715, T[K]:  296.11, density[g/cm3]: 1.100, calc_time[min]:   691.30
    Dyn  step:    68400, time[ps]:   68.40, Etot[eV]:  -6585.8088, Epot[eV]:  -6637.8950, Ekin[eV]:   52.0861, T[K]:  298.49, density[g/cm3]: 1.094, calc_time[min]:   692.39
    Dyn  step:    68500, time[ps]:   68.50, Etot[eV]:  -6585.8170, Epot[eV]:  -6637.9858, Ekin[eV]:   52.1688, T[K]:  298.96, density[g/cm3]: 1.101, calc_time[min]:   693.47
    Dyn  step:    68600, time[ps]:   68.60, Etot[eV]:  -6586.0079, Epot[eV]:  -6638.4603, Ekin[eV]:   52.4524, T[K]:  300.59, density[g/cm3]: 1.093, calc_time[min]:   694.56
    Dyn  step:    68700, time[ps]:   68.70, Etot[eV]:  -6586.1405, Epot[eV]:  -6637.7255, Ekin[eV]:   51.5850, T[K]:  295.61, density[g/cm3]: 1.097, calc_time[min]:   695.59
    Dyn  step:    68800, time[ps]:   68.80, Etot[eV]:  -6585.9384, Epot[eV]:  -6639.0596, Ekin[eV]:   53.1212, T[K]:  304.42, density[g/cm3]: 1.097, calc_time[min]:   696.62
    Dyn  step:    68900, time[ps]:   68.90, Etot[eV]:  -6584.8229, Epot[eV]:  -6636.5832, Ekin[eV]:   51.7603, T[K]:  296.62, density[g/cm3]: 1.098, calc_time[min]:   697.65
    Dyn  step:    69000, time[ps]:   69.00, Etot[eV]:  -6585.0216, Epot[eV]:  -6637.6480, Ekin[eV]:   52.6265, T[K]:  301.58, density[g/cm3]: 1.097, calc_time[min]:   698.67
    Dyn  step:    69100, time[ps]:   69.10, Etot[eV]:  -6585.5915, Epot[eV]:  -6637.8929, Ekin[eV]:   52.3014, T[K]:  299.72, density[g/cm3]: 1.097, calc_time[min]:   699.67
    Dyn  step:    69200, time[ps]:   69.20, Etot[eV]:  -6585.6123, Epot[eV]:  -6637.2766, Ekin[eV]:   51.6643, T[K]:  296.07, density[g/cm3]: 1.094, calc_time[min]:   700.70
    Dyn  step:    69300, time[ps]:   69.30, Etot[eV]:  -6585.2580, Epot[eV]:  -6637.1526, Ekin[eV]:   51.8946, T[K]:  297.39, density[g/cm3]: 1.091, calc_time[min]:   701.68
    Dyn  step:    69400, time[ps]:   69.40, Etot[eV]:  -6585.6358, Epot[eV]:  -6637.4366, Ekin[eV]:   51.8008, T[K]:  296.85, density[g/cm3]: 1.097, calc_time[min]:   702.67
    Dyn  step:    69500, time[ps]:   69.50, Etot[eV]:  -6584.6042, Epot[eV]:  -6637.1853, Ekin[eV]:   52.5811, T[K]:  301.32, density[g/cm3]: 1.095, calc_time[min]:   703.63
    Dyn  step:    69600, time[ps]:   69.60, Etot[eV]:  -6585.1895, Epot[eV]:  -6638.0284, Ekin[eV]:   52.8389, T[K]:  302.80, density[g/cm3]: 1.095, calc_time[min]:   704.57
    Dyn  step:    69700, time[ps]:   69.70, Etot[eV]:  -6585.3372, Epot[eV]:  -6637.0028, Ekin[eV]:   51.6656, T[K]:  296.08, density[g/cm3]: 1.090, calc_time[min]:   705.53
    Dyn  step:    69800, time[ps]:   69.80, Etot[eV]:  -6586.3267, Epot[eV]:  -6638.8040, Ekin[eV]:   52.4773, T[K]:  300.73, density[g/cm3]: 1.088, calc_time[min]:   706.46
    Dyn  step:    69900, time[ps]:   69.90, Etot[eV]:  -6586.4073, Epot[eV]:  -6638.8256, Ekin[eV]:   52.4183, T[K]:  300.39, density[g/cm3]: 1.092, calc_time[min]:   707.40
    Dyn  step:    70000, time[ps]:   70.00, Etot[eV]:  -6585.9713, Epot[eV]:  -6638.2165, Ekin[eV]:   52.2452, T[K]:  299.40, density[g/cm3]: 1.093, calc_time[min]:   708.33
    Dyn  step:    70100, time[ps]:   70.10, Etot[eV]:  -6585.6521, Epot[eV]:  -6637.4187, Ekin[eV]:   51.7666, T[K]:  296.66, density[g/cm3]: 1.097, calc_time[min]:   709.25
    Dyn  step:    70200, time[ps]:   70.20, Etot[eV]:  -6585.5804, Epot[eV]:  -6638.0085, Ekin[eV]:   52.4281, T[K]:  300.45, density[g/cm3]: 1.100, calc_time[min]:   710.18
    Dyn  step:    70300, time[ps]:   70.30, Etot[eV]:  -6586.2673, Epot[eV]:  -6638.0636, Ekin[eV]:   51.7963, T[K]:  296.83, density[g/cm3]: 1.097, calc_time[min]:   711.12
    Dyn  step:    70400, time[ps]:   70.40, Etot[eV]:  -6585.3139, Epot[eV]:  -6637.9871, Ekin[eV]:   52.6732, T[K]:  301.85, density[g/cm3]: 1.093, calc_time[min]:   712.07
    Dyn  step:    70500, time[ps]:   70.50, Etot[eV]:  -6585.7169, Epot[eV]:  -6638.2169, Ekin[eV]:   52.5000, T[K]:  300.86, density[g/cm3]: 1.096, calc_time[min]:   713.01
    Dyn  step:    70600, time[ps]:   70.60, Etot[eV]:  -6585.4345, Epot[eV]:  -6637.4737, Ekin[eV]:   52.0393, T[K]:  298.22, density[g/cm3]: 1.098, calc_time[min]:   713.91
    Dyn  step:    70700, time[ps]:   70.70, Etot[eV]:  -6585.1345, Epot[eV]:  -6637.7002, Ekin[eV]:   52.5657, T[K]:  301.23, density[g/cm3]: 1.094, calc_time[min]:   714.83
    Dyn  step:    70800, time[ps]:   70.80, Etot[eV]:  -6584.7935, Epot[eV]:  -6637.6281, Ekin[eV]:   52.8346, T[K]:  302.78, density[g/cm3]: 1.098, calc_time[min]:   715.76
    Dyn  step:    70900, time[ps]:   70.90, Etot[eV]:  -6585.8007, Epot[eV]:  -6637.0601, Ekin[eV]:   51.2594, T[K]:  293.75, density[g/cm3]: 1.100, calc_time[min]:   716.68
    Dyn  step:    71000, time[ps]:   71.00, Etot[eV]:  -6585.7630, Epot[eV]:  -6637.8991, Ekin[eV]:   52.1361, T[K]:  298.77, density[g/cm3]: 1.101, calc_time[min]:   717.59
    Dyn  step:    71100, time[ps]:   71.10, Etot[eV]:  -6585.6303, Epot[eV]:  -6637.3881, Ekin[eV]:   51.7578, T[K]:  296.60, density[g/cm3]: 1.102, calc_time[min]:   718.51
    Dyn  step:    71200, time[ps]:   71.20, Etot[eV]:  -6586.4873, Epot[eV]:  -6639.4999, Ekin[eV]:   53.0126, T[K]:  303.80, density[g/cm3]: 1.105, calc_time[min]:   719.42
    Dyn  step:    71300, time[ps]:   71.30, Etot[eV]:  -6585.1423, Epot[eV]:  -6636.2332, Ekin[eV]:   51.0909, T[K]:  292.78, density[g/cm3]: 1.104, calc_time[min]:   720.40
    Dyn  step:    71400, time[ps]:   71.40, Etot[eV]:  -6586.0800, Epot[eV]:  -6639.3515, Ekin[eV]:   53.2715, T[K]:  305.28, density[g/cm3]: 1.103, calc_time[min]:   721.31
    Dyn  step:    71500, time[ps]:   71.50, Etot[eV]:  -6586.4891, Epot[eV]:  -6637.8122, Ekin[eV]:   51.3231, T[K]:  294.11, density[g/cm3]: 1.110, calc_time[min]:   722.22
    Dyn  step:    71600, time[ps]:   71.60, Etot[eV]:  -6585.0035, Epot[eV]:  -6635.9387, Ekin[eV]:   50.9353, T[K]:  291.89, density[g/cm3]: 1.105, calc_time[min]:   723.15
    Dyn  step:    71700, time[ps]:   71.70, Etot[eV]:  -6584.9327, Epot[eV]:  -6637.5452, Ekin[eV]:   52.6125, T[K]:  301.50, density[g/cm3]: 1.102, calc_time[min]:   724.06
    Dyn  step:    71800, time[ps]:   71.80, Etot[eV]:  -6585.7456, Epot[eV]:  -6638.5575, Ekin[eV]:   52.8119, T[K]:  302.65, density[g/cm3]: 1.103, calc_time[min]:   724.99
    Dyn  step:    71900, time[ps]:   71.90, Etot[eV]:  -6585.3486, Epot[eV]:  -6638.1961, Ekin[eV]:   52.8476, T[K]:  302.85, density[g/cm3]: 1.095, calc_time[min]:   725.93
    Dyn  step:    72000, time[ps]:   72.00, Etot[eV]:  -6584.5535, Epot[eV]:  -6637.3616, Ekin[eV]:   52.8082, T[K]:  302.62, density[g/cm3]: 1.091, calc_time[min]:   726.89
    Dyn  step:    72100, time[ps]:   72.10, Etot[eV]:  -6584.8985, Epot[eV]:  -6636.4521, Ekin[eV]:   51.5536, T[K]:  295.43, density[g/cm3]: 1.087, calc_time[min]:   727.81
    Dyn  step:    72200, time[ps]:   72.20, Etot[eV]:  -6584.8593, Epot[eV]:  -6636.0958, Ekin[eV]:   51.2366, T[K]:  293.62, density[g/cm3]: 1.087, calc_time[min]:   728.74
    Dyn  step:    72300, time[ps]:   72.30, Etot[eV]:  -6584.6875, Epot[eV]:  -6638.8135, Ekin[eV]:   54.1260, T[K]:  310.18, density[g/cm3]: 1.085, calc_time[min]:   729.66
    Dyn  step:    72400, time[ps]:   72.40, Etot[eV]:  -6585.0795, Epot[eV]:  -6637.6593, Ekin[eV]:   52.5797, T[K]:  301.31, density[g/cm3]: 1.089, calc_time[min]:   730.55
    Dyn  step:    72500, time[ps]:   72.50, Etot[eV]:  -6585.2835, Epot[eV]:  -6637.1687, Ekin[eV]:   51.8852, T[K]:  297.33, density[g/cm3]: 1.098, calc_time[min]:   731.50
    Dyn  step:    72600, time[ps]:   72.60, Etot[eV]:  -6585.0759, Epot[eV]:  -6638.1310, Ekin[eV]:   53.0551, T[K]:  304.04, density[g/cm3]: 1.097, calc_time[min]:   732.43
    Dyn  step:    72700, time[ps]:   72.70, Etot[eV]:  -6585.5935, Epot[eV]:  -6637.5851, Ekin[eV]:   51.9915, T[K]:  297.94, density[g/cm3]: 1.099, calc_time[min]:   733.39
    Dyn  step:    72800, time[ps]:   72.80, Etot[eV]:  -6585.6538, Epot[eV]:  -6638.1596, Ekin[eV]:   52.5058, T[K]:  300.89, density[g/cm3]: 1.094, calc_time[min]:   734.30
    Dyn  step:    72900, time[ps]:   72.90, Etot[eV]:  -6586.0916, Epot[eV]:  -6637.4220, Ekin[eV]:   51.3305, T[K]:  294.16, density[g/cm3]: 1.094, calc_time[min]:   735.25
    Dyn  step:    73000, time[ps]:   73.00, Etot[eV]:  -6585.5411, Epot[eV]:  -6636.8563, Ekin[eV]:   51.3153, T[K]:  294.07, density[g/cm3]: 1.093, calc_time[min]:   736.19
    Dyn  step:    73100, time[ps]:   73.10, Etot[eV]:  -6585.8782, Epot[eV]:  -6637.0180, Ekin[eV]:   51.1397, T[K]:  293.06, density[g/cm3]: 1.094, calc_time[min]:   737.13
    Dyn  step:    73200, time[ps]:   73.20, Etot[eV]:  -6586.0706, Epot[eV]:  -6639.7399, Ekin[eV]:   53.6693, T[K]:  307.56, density[g/cm3]: 1.093, calc_time[min]:   738.08
    Dyn  step:    73300, time[ps]:   73.30, Etot[eV]:  -6586.9956, Epot[eV]:  -6639.3365, Ekin[eV]:   52.3409, T[K]:  299.95, density[g/cm3]: 1.104, calc_time[min]:   739.03
    Dyn  step:    73400, time[ps]:   73.40, Etot[eV]:  -6586.8406, Epot[eV]:  -6640.1218, Ekin[eV]:   53.2812, T[K]:  305.33, density[g/cm3]: 1.108, calc_time[min]:   739.99
    Dyn  step:    73500, time[ps]:   73.50, Etot[eV]:  -6585.9270, Epot[eV]:  -6637.4580, Ekin[eV]:   51.5311, T[K]:  295.31, density[g/cm3]: 1.107, calc_time[min]:   740.93
    Dyn  step:    73600, time[ps]:   73.60, Etot[eV]:  -6585.8579, Epot[eV]:  -6636.7366, Ekin[eV]:   50.8787, T[K]:  291.57, density[g/cm3]: 1.101, calc_time[min]:   741.90
    Dyn  step:    73700, time[ps]:   73.70, Etot[eV]:  -6585.8647, Epot[eV]:  -6638.1041, Ekin[eV]:   52.2394, T[K]:  299.36, density[g/cm3]: 1.097, calc_time[min]:   742.85
    Dyn  step:    73800, time[ps]:   73.80, Etot[eV]:  -6586.3658, Epot[eV]:  -6637.7914, Ekin[eV]:   51.4256, T[K]:  294.70, density[g/cm3]: 1.105, calc_time[min]:   743.79
    Dyn  step:    73900, time[ps]:   73.90, Etot[eV]:  -6585.5458, Epot[eV]:  -6637.8694, Ekin[eV]:   52.3236, T[K]:  299.85, density[g/cm3]: 1.101, calc_time[min]:   744.75
    Dyn  step:    74000, time[ps]:   74.00, Etot[eV]:  -6586.2492, Epot[eV]:  -6638.4079, Ekin[eV]:   52.1587, T[K]:  298.90, density[g/cm3]: 1.104, calc_time[min]:   745.69
    Dyn  step:    74100, time[ps]:   74.10, Etot[eV]:  -6586.0960, Epot[eV]:  -6637.7118, Ekin[eV]:   51.6158, T[K]:  295.79, density[g/cm3]: 1.099, calc_time[min]:   746.68
    Dyn  step:    74200, time[ps]:   74.20, Etot[eV]:  -6585.2770, Epot[eV]:  -6636.8599, Ekin[eV]:   51.5830, T[K]:  295.60, density[g/cm3]: 1.100, calc_time[min]:   747.66
    Dyn  step:    74300, time[ps]:   74.30, Etot[eV]:  -6585.5032, Epot[eV]:  -6638.4219, Ekin[eV]:   52.9187, T[K]:  303.26, density[g/cm3]: 1.099, calc_time[min]:   748.62
    Dyn  step:    74400, time[ps]:   74.40, Etot[eV]:  -6585.4837, Epot[eV]:  -6638.5282, Ekin[eV]:   53.0445, T[K]:  303.98, density[g/cm3]: 1.097, calc_time[min]:   749.59
    Dyn  step:    74500, time[ps]:   74.50, Etot[eV]:  -6585.5072, Epot[eV]:  -6637.7094, Ekin[eV]:   52.2021, T[K]:  299.15, density[g/cm3]: 1.096, calc_time[min]:   750.56
    Dyn  step:    74600, time[ps]:   74.60, Etot[eV]:  -6586.5973, Epot[eV]:  -6639.9259, Ekin[eV]:   53.3286, T[K]:  305.61, density[g/cm3]: 1.098, calc_time[min]:   751.54
    Dyn  step:    74700, time[ps]:   74.70, Etot[eV]:  -6585.1829, Epot[eV]:  -6636.9854, Ekin[eV]:   51.8025, T[K]:  296.86, density[g/cm3]: 1.097, calc_time[min]:   752.52
    Dyn  step:    74800, time[ps]:   74.80, Etot[eV]:  -6585.7451, Epot[eV]:  -6636.5154, Ekin[eV]:   50.7703, T[K]:  290.95, density[g/cm3]: 1.101, calc_time[min]:   753.51
    Dyn  step:    74900, time[ps]:   74.90, Etot[eV]:  -6584.2119, Epot[eV]:  -6636.7567, Ekin[eV]:   52.5448, T[K]:  301.11, density[g/cm3]: 1.098, calc_time[min]:   754.49
    Dyn  step:    75000, time[ps]:   75.00, Etot[eV]:  -6585.1859, Epot[eV]:  -6636.5542, Ekin[eV]:   51.3683, T[K]:  294.37, density[g/cm3]: 1.098, calc_time[min]:   755.48
    Dyn  step:    75100, time[ps]:   75.10, Etot[eV]:  -6585.6203, Epot[eV]:  -6637.5256, Ekin[eV]:   51.9053, T[K]:  297.45, density[g/cm3]: 1.098, calc_time[min]:   756.46
    Dyn  step:    75200, time[ps]:   75.20, Etot[eV]:  -6585.9641, Epot[eV]:  -6637.3492, Ekin[eV]:   51.3851, T[K]:  294.47, density[g/cm3]: 1.102, calc_time[min]:   757.47
    Dyn  step:    75300, time[ps]:   75.30, Etot[eV]:  -6585.5251, Epot[eV]:  -6637.9780, Ekin[eV]:   52.4529, T[K]:  300.59, density[g/cm3]: 1.099, calc_time[min]:   758.44
    Dyn  step:    75400, time[ps]:   75.40, Etot[eV]:  -6585.3664, Epot[eV]:  -6637.8317, Ekin[eV]:   52.4653, T[K]:  300.66, density[g/cm3]: 1.097, calc_time[min]:   759.43
    Dyn  step:    75500, time[ps]:   75.50, Etot[eV]:  -6585.4965, Epot[eV]:  -6638.0834, Ekin[eV]:   52.5869, T[K]:  301.36, density[g/cm3]: 1.094, calc_time[min]:   760.42
    Dyn  step:    75600, time[ps]:   75.60, Etot[eV]:  -6585.5977, Epot[eV]:  -6637.3360, Ekin[eV]:   51.7383, T[K]:  296.49, density[g/cm3]: 1.090, calc_time[min]:   761.40
    Dyn  step:    75700, time[ps]:   75.70, Etot[eV]:  -6584.7496, Epot[eV]:  -6636.7579, Ekin[eV]:   52.0083, T[K]:  298.04, density[g/cm3]: 1.085, calc_time[min]:   762.41
    Dyn  step:    75800, time[ps]:   75.80, Etot[eV]:  -6584.7536, Epot[eV]:  -6635.5763, Ekin[eV]:   50.8227, T[K]:  291.25, density[g/cm3]: 1.087, calc_time[min]:   763.42
    Dyn  step:    75900, time[ps]:   75.90, Etot[eV]:  -6586.2004, Epot[eV]:  -6639.5005, Ekin[eV]:   53.3001, T[K]:  305.44, density[g/cm3]: 1.081, calc_time[min]:   764.39
    Dyn  step:    76000, time[ps]:   76.00, Etot[eV]:  -6586.4509, Epot[eV]:  -6638.3354, Ekin[eV]:   51.8845, T[K]:  297.33, density[g/cm3]: 1.090, calc_time[min]:   765.39
    Dyn  step:    76100, time[ps]:   76.10, Etot[eV]:  -6586.4856, Epot[eV]:  -6638.6876, Ekin[eV]:   52.2020, T[K]:  299.15, density[g/cm3]: 1.098, calc_time[min]:   766.40
    Dyn  step:    76200, time[ps]:   76.20, Etot[eV]:  -6585.5012, Epot[eV]:  -6637.4813, Ekin[eV]:   51.9801, T[K]:  297.88, density[g/cm3]: 1.094, calc_time[min]:   767.41
    Dyn  step:    76300, time[ps]:   76.30, Etot[eV]:  -6585.9847, Epot[eV]:  -6638.3920, Ekin[eV]:   52.4073, T[K]:  300.33, density[g/cm3]: 1.105, calc_time[min]:   768.44
    Dyn  step:    76400, time[ps]:   76.40, Etot[eV]:  -6586.1138, Epot[eV]:  -6638.2169, Ekin[eV]:   52.1030, T[K]:  298.58, density[g/cm3]: 1.106, calc_time[min]:   769.46
    Dyn  step:    76500, time[ps]:   76.50, Etot[eV]:  -6586.0583, Epot[eV]:  -6637.2407, Ekin[eV]:   51.1824, T[K]:  293.31, density[g/cm3]: 1.103, calc_time[min]:   770.48
    Dyn  step:    76600, time[ps]:   76.60, Etot[eV]:  -6586.0981, Epot[eV]:  -6638.3834, Ekin[eV]:   52.2853, T[K]:  299.63, density[g/cm3]: 1.108, calc_time[min]:   771.49
    Dyn  step:    76700, time[ps]:   76.70, Etot[eV]:  -6586.2021, Epot[eV]:  -6638.0585, Ekin[eV]:   51.8564, T[K]:  297.17, density[g/cm3]: 1.107, calc_time[min]:   772.48
    Dyn  step:    76800, time[ps]:   76.80, Etot[eV]:  -6586.0397, Epot[eV]:  -6639.5735, Ekin[eV]:   53.5337, T[K]:  306.78, density[g/cm3]: 1.110, calc_time[min]:   773.49
    Dyn  step:    76900, time[ps]:   76.90, Etot[eV]:  -6585.4877, Epot[eV]:  -6636.4667, Ekin[eV]:   50.9790, T[K]:  292.14, density[g/cm3]: 1.108, calc_time[min]:   774.50
    Dyn  step:    77000, time[ps]:   77.00, Etot[eV]:  -6585.5730, Epot[eV]:  -6639.0741, Ekin[eV]:   53.5011, T[K]:  306.59, density[g/cm3]: 1.108, calc_time[min]:   775.49
    Dyn  step:    77100, time[ps]:   77.10, Etot[eV]:  -6586.3581, Epot[eV]:  -6639.2052, Ekin[eV]:   52.8471, T[K]:  302.85, density[g/cm3]: 1.111, calc_time[min]:   776.45
    Dyn  step:    77200, time[ps]:   77.20, Etot[eV]:  -6585.3209, Epot[eV]:  -6639.6525, Ekin[eV]:   54.3316, T[K]:  311.35, density[g/cm3]: 1.101, calc_time[min]:   777.47
    Dyn  step:    77300, time[ps]:   77.30, Etot[eV]:  -6585.7913, Epot[eV]:  -6638.2458, Ekin[eV]:   52.4544, T[K]:  300.60, density[g/cm3]: 1.099, calc_time[min]:   778.48
    Dyn  step:    77400, time[ps]:   77.40, Etot[eV]:  -6586.1784, Epot[eV]:  -6639.3099, Ekin[eV]:   53.1315, T[K]:  304.48, density[g/cm3]: 1.103, calc_time[min]:   779.43
    Dyn  step:    77500, time[ps]:   77.50, Etot[eV]:  -6586.3977, Epot[eV]:  -6639.1331, Ekin[eV]:   52.7354, T[K]:  302.21, density[g/cm3]: 1.103, calc_time[min]:   780.43
    Dyn  step:    77600, time[ps]:   77.60, Etot[eV]:  -6585.7909, Epot[eV]:  -6639.0550, Ekin[eV]:   53.2641, T[K]:  305.24, density[g/cm3]: 1.108, calc_time[min]:   781.43
    Dyn  step:    77700, time[ps]:   77.70, Etot[eV]:  -6586.0377, Epot[eV]:  -6638.2925, Ekin[eV]:   52.2548, T[K]:  299.45, density[g/cm3]: 1.109, calc_time[min]:   782.40
    Dyn  step:    77800, time[ps]:   77.80, Etot[eV]:  -6585.6589, Epot[eV]:  -6638.9706, Ekin[eV]:   53.3117, T[K]:  305.51, density[g/cm3]: 1.107, calc_time[min]:   783.38
    Dyn  step:    77900, time[ps]:   77.90, Etot[eV]:  -6585.3304, Epot[eV]:  -6638.2400, Ekin[eV]:   52.9095, T[K]:  303.20, density[g/cm3]: 1.104, calc_time[min]:   784.37
    Dyn  step:    78000, time[ps]:   78.00, Etot[eV]:  -6585.2156, Epot[eV]:  -6635.6444, Ekin[eV]:   50.4289, T[K]:  288.99, density[g/cm3]: 1.105, calc_time[min]:   785.36
    Dyn  step:    78100, time[ps]:   78.10, Etot[eV]:  -6585.7136, Epot[eV]:  -6639.2503, Ekin[eV]:   53.5368, T[K]:  306.80, density[g/cm3]: 1.101, calc_time[min]:   786.35
    Dyn  step:    78200, time[ps]:   78.20, Etot[eV]:  -6586.0342, Epot[eV]:  -6636.6992, Ekin[eV]:   50.6650, T[K]:  290.34, density[g/cm3]: 1.108, calc_time[min]:   787.35
    Dyn  step:    78300, time[ps]:   78.30, Etot[eV]:  -6585.3621, Epot[eV]:  -6637.1577, Ekin[eV]:   51.7956, T[K]:  296.82, density[g/cm3]: 1.108, calc_time[min]:   788.32
    Dyn  step:    78400, time[ps]:   78.40, Etot[eV]:  -6585.3362, Epot[eV]:  -6637.8904, Ekin[eV]:   52.5541, T[K]:  301.17, density[g/cm3]: 1.106, calc_time[min]:   789.26
    Dyn  step:    78500, time[ps]:   78.50, Etot[eV]:  -6585.8900, Epot[eV]:  -6637.6365, Ekin[eV]:   51.7464, T[K]:  296.54, density[g/cm3]: 1.109, calc_time[min]:   790.24
    Dyn  step:    78600, time[ps]:   78.60, Etot[eV]:  -6586.1792, Epot[eV]:  -6639.0811, Ekin[eV]:   52.9019, T[K]:  303.16, density[g/cm3]: 1.111, calc_time[min]:   791.19
    Dyn  step:    78700, time[ps]:   78.70, Etot[eV]:  -6586.5878, Epot[eV]:  -6639.9473, Ekin[eV]:   53.3595, T[K]:  305.78, density[g/cm3]: 1.121, calc_time[min]:   792.19
    Dyn  step:    78800, time[ps]:   78.80, Etot[eV]:  -6586.3925, Epot[eV]:  -6639.3141, Ekin[eV]:   52.9216, T[K]:  303.27, density[g/cm3]: 1.116, calc_time[min]:   793.20
    Dyn  step:    78900, time[ps]:   78.90, Etot[eV]:  -6586.1681, Epot[eV]:  -6638.5638, Ekin[eV]:   52.3957, T[K]:  300.26, density[g/cm3]: 1.110, calc_time[min]:   794.20
    Dyn  step:    79000, time[ps]:   79.00, Etot[eV]:  -6585.9227, Epot[eV]:  -6637.7305, Ekin[eV]:   51.8078, T[K]:  296.89, density[g/cm3]: 1.108, calc_time[min]:   795.21
    Dyn  step:    79100, time[ps]:   79.10, Etot[eV]:  -6585.9514, Epot[eV]:  -6638.3658, Ekin[eV]:   52.4144, T[K]:  300.37, density[g/cm3]: 1.106, calc_time[min]:   796.21
    Dyn  step:    79200, time[ps]:   79.20, Etot[eV]:  -6585.7235, Epot[eV]:  -6636.9833, Ekin[eV]:   51.2597, T[K]:  293.75, density[g/cm3]: 1.105, calc_time[min]:   797.18
    Dyn  step:    79300, time[ps]:   79.30, Etot[eV]:  -6585.7017, Epot[eV]:  -6638.8492, Ekin[eV]:   53.1475, T[K]:  304.57, density[g/cm3]: 1.111, calc_time[min]:   798.16
    Dyn  step:    79400, time[ps]:   79.40, Etot[eV]:  -6586.4337, Epot[eV]:  -6639.2334, Ekin[eV]:   52.7997, T[K]:  302.58, density[g/cm3]: 1.113, calc_time[min]:   799.14
    Dyn  step:    79500, time[ps]:   79.50, Etot[eV]:  -6586.3961, Epot[eV]:  -6638.4060, Ekin[eV]:   52.0099, T[K]:  298.05, density[g/cm3]: 1.106, calc_time[min]:   800.13
    Dyn  step:    79600, time[ps]:   79.60, Etot[eV]:  -6586.4722, Epot[eV]:  -6639.1647, Ekin[eV]:   52.6924, T[K]:  301.96, density[g/cm3]: 1.107, calc_time[min]:   801.13
    Dyn  step:    79700, time[ps]:   79.70, Etot[eV]:  -6586.1360, Epot[eV]:  -6638.5701, Ekin[eV]:   52.4341, T[K]:  300.48, density[g/cm3]: 1.107, calc_time[min]:   802.15
    Dyn  step:    79800, time[ps]:   79.80, Etot[eV]:  -6586.3952, Epot[eV]:  -6638.3168, Ekin[eV]:   51.9215, T[K]:  297.54, density[g/cm3]: 1.110, calc_time[min]:   803.15
    Dyn  step:    79900, time[ps]:   79.90, Etot[eV]:  -6586.2924, Epot[eV]:  -6638.4950, Ekin[eV]:   52.2027, T[K]:  299.15, density[g/cm3]: 1.105, calc_time[min]:   804.18
    Dyn  step:    80000, time[ps]:   80.00, Etot[eV]:  -6586.8986, Epot[eV]:  -6640.2624, Ekin[eV]:   53.3638, T[K]:  305.81, density[g/cm3]: 1.106, calc_time[min]:   805.22
    Dyn  step:    80100, time[ps]:   80.10, Etot[eV]:  -6586.2272, Epot[eV]:  -6639.2201, Ekin[eV]:   52.9929, T[K]:  303.68, density[g/cm3]: 1.108, calc_time[min]:   806.22
    Dyn  step:    80200, time[ps]:   80.20, Etot[eV]:  -6586.0783, Epot[eV]:  -6637.5656, Ekin[eV]:   51.4874, T[K]:  295.05, density[g/cm3]: 1.105, calc_time[min]:   807.24
    Dyn  step:    80300, time[ps]:   80.30, Etot[eV]:  -6585.9616, Epot[eV]:  -6637.3425, Ekin[eV]:   51.3809, T[K]:  294.45, density[g/cm3]: 1.105, calc_time[min]:   808.25
    Dyn  step:    80400, time[ps]:   80.40, Etot[eV]:  -6585.9968, Epot[eV]:  -6638.9725, Ekin[eV]:   52.9757, T[K]:  303.58, density[g/cm3]: 1.101, calc_time[min]:   809.26
    Dyn  step:    80500, time[ps]:   80.50, Etot[eV]:  -6586.6728, Epot[eV]:  -6639.6224, Ekin[eV]:   52.9496, T[K]:  303.43, density[g/cm3]: 1.101, calc_time[min]:   810.30
    Dyn  step:    80600, time[ps]:   80.60, Etot[eV]:  -6585.9185, Epot[eV]:  -6638.4553, Ekin[eV]:   52.5368, T[K]:  301.07, density[g/cm3]: 1.104, calc_time[min]:   811.30
    Dyn  step:    80700, time[ps]:   80.70, Etot[eV]:  -6585.5474, Epot[eV]:  -6637.2264, Ekin[eV]:   51.6791, T[K]:  296.15, density[g/cm3]: 1.104, calc_time[min]:   812.33
    Dyn  step:    80800, time[ps]:   80.80, Etot[eV]:  -6585.8151, Epot[eV]:  -6638.5049, Ekin[eV]:   52.6897, T[K]:  301.95, density[g/cm3]: 1.106, calc_time[min]:   813.33
    Dyn  step:    80900, time[ps]:   80.90, Etot[eV]:  -6586.6588, Epot[eV]:  -6638.2506, Ekin[eV]:   51.5918, T[K]:  295.65, density[g/cm3]: 1.107, calc_time[min]:   814.32
    Dyn  step:    81000, time[ps]:   81.00, Etot[eV]:  -6585.9159, Epot[eV]:  -6638.4538, Ekin[eV]:   52.5379, T[K]:  301.07, density[g/cm3]: 1.108, calc_time[min]:   815.34
    Dyn  step:    81100, time[ps]:   81.10, Etot[eV]:  -6585.7165, Epot[eV]:  -6637.2364, Ekin[eV]:   51.5199, T[K]:  295.24, density[g/cm3]: 1.107, calc_time[min]:   816.33
    Dyn  step:    81200, time[ps]:   81.20, Etot[eV]:  -6586.0489, Epot[eV]:  -6639.0976, Ekin[eV]:   53.0487, T[K]:  304.00, density[g/cm3]: 1.106, calc_time[min]:   817.34
    Dyn  step:    81300, time[ps]:   81.30, Etot[eV]:  -6586.0682, Epot[eV]:  -6637.1433, Ekin[eV]:   51.0751, T[K]:  292.69, density[g/cm3]: 1.103, calc_time[min]:   818.38
    Dyn  step:    81400, time[ps]:   81.40, Etot[eV]:  -6586.5673, Epot[eV]:  -6637.9563, Ekin[eV]:   51.3891, T[K]:  294.49, density[g/cm3]: 1.100, calc_time[min]:   819.41
    Dyn  step:    81500, time[ps]:   81.50, Etot[eV]:  -6586.3397, Epot[eV]:  -6638.1776, Ekin[eV]:   51.8379, T[K]:  297.06, density[g/cm3]: 1.111, calc_time[min]:   820.46
    Dyn  step:    81600, time[ps]:   81.60, Etot[eV]:  -6585.7692, Epot[eV]:  -6639.2137, Ekin[eV]:   53.4445, T[K]:  306.27, density[g/cm3]: 1.118, calc_time[min]:   821.45
    Dyn  step:    81700, time[ps]:   81.70, Etot[eV]:  -6586.2726, Epot[eV]:  -6638.8817, Ekin[eV]:   52.6091, T[K]:  301.48, density[g/cm3]: 1.123, calc_time[min]:   822.47
    Dyn  step:    81800, time[ps]:   81.80, Etot[eV]:  -6586.1771, Epot[eV]:  -6640.7452, Ekin[eV]:   54.5681, T[K]:  312.71, density[g/cm3]: 1.125, calc_time[min]:   823.51
    Dyn  step:    81900, time[ps]:   81.90, Etot[eV]:  -6586.7137, Epot[eV]:  -6639.2088, Ekin[eV]:   52.4952, T[K]:  300.83, density[g/cm3]: 1.127, calc_time[min]:   824.57
    Dyn  step:    82000, time[ps]:   82.00, Etot[eV]:  -6586.9462, Epot[eV]:  -6639.6496, Ekin[eV]:   52.7035, T[K]:  302.02, density[g/cm3]: 1.128, calc_time[min]:   825.91
    Dyn  step:    82100, time[ps]:   82.10, Etot[eV]:  -6586.4519, Epot[eV]:  -6638.6634, Ekin[eV]:   52.2115, T[K]:  299.20, density[g/cm3]: 1.122, calc_time[min]:   827.46
    Dyn  step:    82200, time[ps]:   82.20, Etot[eV]:  -6585.9214, Epot[eV]:  -6637.7213, Ekin[eV]:   51.7999, T[K]:  296.85, density[g/cm3]: 1.118, calc_time[min]:   829.19
    Dyn  step:    82300, time[ps]:   82.30, Etot[eV]:  -6585.3159, Epot[eV]:  -6639.0707, Ekin[eV]:   53.7548, T[K]:  308.05, density[g/cm3]: 1.119, calc_time[min]:   830.90
    Dyn  step:    82400, time[ps]:   82.40, Etot[eV]:  -6585.8347, Epot[eV]:  -6638.3711, Ekin[eV]:   52.5364, T[K]:  301.07, density[g/cm3]: 1.115, calc_time[min]:   832.18
    Dyn  step:    82500, time[ps]:   82.50, Etot[eV]:  -6586.4024, Epot[eV]:  -6639.4265, Ekin[eV]:   53.0241, T[K]:  303.86, density[g/cm3]: 1.117, calc_time[min]:   833.17
    Dyn  step:    82600, time[ps]:   82.60, Etot[eV]:  -6586.0778, Epot[eV]:  -6638.9199, Ekin[eV]:   52.8421, T[K]:  302.82, density[g/cm3]: 1.120, calc_time[min]:   834.12
    Dyn  step:    82700, time[ps]:   82.70, Etot[eV]:  -6585.8010, Epot[eV]:  -6639.1148, Ekin[eV]:   53.3138, T[K]:  305.52, density[g/cm3]: 1.120, calc_time[min]:   835.03
    Dyn  step:    82800, time[ps]:   82.80, Etot[eV]:  -6585.9721, Epot[eV]:  -6637.1195, Ekin[eV]:   51.1474, T[K]:  293.11, density[g/cm3]: 1.119, calc_time[min]:   835.99
    Dyn  step:    82900, time[ps]:   82.90, Etot[eV]:  -6585.7260, Epot[eV]:  -6637.7193, Ekin[eV]:   51.9932, T[K]:  297.95, density[g/cm3]: 1.120, calc_time[min]:   836.94
    Dyn  step:    83000, time[ps]:   83.00, Etot[eV]:  -6586.1344, Epot[eV]:  -6638.7413, Ekin[eV]:   52.6069, T[K]:  301.47, density[g/cm3]: 1.113, calc_time[min]:   837.82
    Dyn  step:    83100, time[ps]:   83.10, Etot[eV]:  -6586.0055, Epot[eV]:  -6638.2357, Ekin[eV]:   52.2302, T[K]:  299.31, density[g/cm3]: 1.119, calc_time[min]:   838.73
    Dyn  step:    83200, time[ps]:   83.20, Etot[eV]:  -6586.1194, Epot[eV]:  -6638.1863, Ekin[eV]:   52.0669, T[K]:  298.38, density[g/cm3]: 1.120, calc_time[min]:   839.61
    Dyn  step:    83300, time[ps]:   83.30, Etot[eV]:  -6586.0117, Epot[eV]:  -6638.2271, Ekin[eV]:   52.2154, T[K]:  299.23, density[g/cm3]: 1.118, calc_time[min]:   840.47
    Dyn  step:    83400, time[ps]:   83.40, Etot[eV]:  -6586.6558, Epot[eV]:  -6638.7364, Ekin[eV]:   52.0806, T[K]:  298.45, density[g/cm3]: 1.115, calc_time[min]:   841.34
    Dyn  step:    83500, time[ps]:   83.50, Etot[eV]:  -6586.5811, Epot[eV]:  -6639.3659, Ekin[eV]:   52.7848, T[K]:  302.49, density[g/cm3]: 1.117, calc_time[min]:   842.20
    Dyn  step:    83600, time[ps]:   83.60, Etot[eV]:  -6586.0085, Epot[eV]:  -6638.7356, Ekin[eV]:   52.7271, T[K]:  302.16, density[g/cm3]: 1.116, calc_time[min]:   843.06
    Dyn  step:    83700, time[ps]:   83.70, Etot[eV]:  -6586.1363, Epot[eV]:  -6639.0298, Ekin[eV]:   52.8936, T[K]:  303.11, density[g/cm3]: 1.113, calc_time[min]:   843.93
    Dyn  step:    83800, time[ps]:   83.80, Etot[eV]:  -6586.5479, Epot[eV]:  -6638.1235, Ekin[eV]:   51.5756, T[K]:  295.56, density[g/cm3]: 1.109, calc_time[min]:   844.78
    Dyn  step:    83900, time[ps]:   83.90, Etot[eV]:  -6585.7838, Epot[eV]:  -6637.5998, Ekin[eV]:   51.8160, T[K]:  296.94, density[g/cm3]: 1.107, calc_time[min]:   845.64
    Dyn  step:    84000, time[ps]:   84.00, Etot[eV]:  -6585.5321, Epot[eV]:  -6638.4787, Ekin[eV]:   52.9466, T[K]:  303.42, density[g/cm3]: 1.107, calc_time[min]:   846.51
    Dyn  step:    84100, time[ps]:   84.10, Etot[eV]:  -6585.4479, Epot[eV]:  -6637.8197, Ekin[eV]:   52.3719, T[K]:  300.12, density[g/cm3]: 1.108, calc_time[min]:   847.33
    Dyn  step:    84200, time[ps]:   84.20, Etot[eV]:  -6585.7087, Epot[eV]:  -6640.4204, Ekin[eV]:   54.7117, T[K]:  313.53, density[g/cm3]: 1.104, calc_time[min]:   848.18
    Dyn  step:    84300, time[ps]:   84.30, Etot[eV]:  -6585.2809, Epot[eV]:  -6637.1540, Ekin[eV]:   51.8731, T[K]:  297.27, density[g/cm3]: 1.105, calc_time[min]:   849.04
    Dyn  step:    84400, time[ps]:   84.40, Etot[eV]:  -6586.0924, Epot[eV]:  -6638.6029, Ekin[eV]:   52.5105, T[K]:  300.92, density[g/cm3]: 1.106, calc_time[min]:   849.89
    Dyn  step:    84500, time[ps]:   84.50, Etot[eV]:  -6586.2563, Epot[eV]:  -6638.0395, Ekin[eV]:   51.7832, T[K]:  296.75, density[g/cm3]: 1.108, calc_time[min]:   850.76
    Dyn  step:    84600, time[ps]:   84.60, Etot[eV]:  -6586.3705, Epot[eV]:  -6638.9216, Ekin[eV]:   52.5512, T[K]:  301.15, density[g/cm3]: 1.116, calc_time[min]:   851.64
    Dyn  step:    84700, time[ps]:   84.70, Etot[eV]:  -6586.2313, Epot[eV]:  -6638.2880, Ekin[eV]:   52.0567, T[K]:  298.32, density[g/cm3]: 1.116, calc_time[min]:   852.49
    Dyn  step:    84800, time[ps]:   84.80, Etot[eV]:  -6585.9522, Epot[eV]:  -6638.7132, Ekin[eV]:   52.7610, T[K]:  302.35, density[g/cm3]: 1.115, calc_time[min]:   853.36
    Dyn  step:    84900, time[ps]:   84.90, Etot[eV]:  -6585.9270, Epot[eV]:  -6638.7199, Ekin[eV]:   52.7929, T[K]:  302.54, density[g/cm3]: 1.126, calc_time[min]:   854.20
    Dyn  step:    85000, time[ps]:   85.00, Etot[eV]:  -6585.8619, Epot[eV]:  -6638.4895, Ekin[eV]:   52.6276, T[K]:  301.59, density[g/cm3]: 1.129, calc_time[min]:   855.04
    Dyn  step:    85100, time[ps]:   85.10, Etot[eV]:  -6585.9502, Epot[eV]:  -6638.3333, Ekin[eV]:   52.3831, T[K]:  300.19, density[g/cm3]: 1.132, calc_time[min]:   855.86
    Dyn  step:    85200, time[ps]:   85.20, Etot[eV]:  -6586.0133, Epot[eV]:  -6638.4700, Ekin[eV]:   52.4567, T[K]:  300.61, density[g/cm3]: 1.126, calc_time[min]:   856.72
    Dyn  step:    85300, time[ps]:   85.30, Etot[eV]:  -6586.0763, Epot[eV]:  -6637.9550, Ekin[eV]:   51.8787, T[K]:  297.30, density[g/cm3]: 1.129, calc_time[min]:   857.58
    Dyn  step:    85400, time[ps]:   85.40, Etot[eV]:  -6586.4640, Epot[eV]:  -6637.8814, Ekin[eV]:   51.4174, T[K]:  294.65, density[g/cm3]: 1.119, calc_time[min]:   858.44
    Dyn  step:    85500, time[ps]:   85.50, Etot[eV]:  -6585.9945, Epot[eV]:  -6637.9982, Ekin[eV]:   52.0037, T[K]:  298.01, density[g/cm3]: 1.119, calc_time[min]:   859.32
    Dyn  step:    85600, time[ps]:   85.60, Etot[eV]:  -6585.7459, Epot[eV]:  -6637.0192, Ekin[eV]:   51.2733, T[K]:  293.83, density[g/cm3]: 1.120, calc_time[min]:   860.17
    Dyn  step:    85700, time[ps]:   85.70, Etot[eV]:  -6586.1514, Epot[eV]:  -6637.8373, Ekin[eV]:   51.6859, T[K]:  296.19, density[g/cm3]: 1.116, calc_time[min]:   861.00
    Dyn  step:    85800, time[ps]:   85.80, Etot[eV]:  -6585.5377, Epot[eV]:  -6637.7903, Ekin[eV]:   52.2525, T[K]:  299.44, density[g/cm3]: 1.121, calc_time[min]:   861.86
    Dyn  step:    85900, time[ps]:   85.90, Etot[eV]:  -6586.1757, Epot[eV]:  -6637.6630, Ekin[eV]:   51.4873, T[K]:  295.05, density[g/cm3]: 1.127, calc_time[min]:   862.75
    Dyn  step:    86000, time[ps]:   86.00, Etot[eV]:  -6585.8749, Epot[eV]:  -6638.6763, Ekin[eV]:   52.8014, T[K]:  302.59, density[g/cm3]: 1.122, calc_time[min]:   863.60
    Dyn  step:    86100, time[ps]:   86.10, Etot[eV]:  -6586.2279, Epot[eV]:  -6637.9783, Ekin[eV]:   51.7504, T[K]:  296.56, density[g/cm3]: 1.120, calc_time[min]:   864.46
    Dyn  step:    86200, time[ps]:   86.20, Etot[eV]:  -6586.2719, Epot[eV]:  -6639.5164, Ekin[eV]:   53.2445, T[K]:  305.12, density[g/cm3]: 1.120, calc_time[min]:   865.32
    Dyn  step:    86300, time[ps]:   86.30, Etot[eV]:  -6586.6582, Epot[eV]:  -6639.2374, Ekin[eV]:   52.5792, T[K]:  301.31, density[g/cm3]: 1.121, calc_time[min]:   866.17
    Dyn  step:    86400, time[ps]:   86.40, Etot[eV]:  -6586.1029, Epot[eV]:  -6638.0763, Ekin[eV]:   51.9734, T[K]:  297.84, density[g/cm3]: 1.117, calc_time[min]:   867.01
    Dyn  step:    86500, time[ps]:   86.50, Etot[eV]:  -6585.7135, Epot[eV]:  -6637.5890, Ekin[eV]:   51.8755, T[K]:  297.28, density[g/cm3]: 1.110, calc_time[min]:   867.86
    Dyn  step:    86600, time[ps]:   86.60, Etot[eV]:  -6586.1844, Epot[eV]:  -6638.7924, Ekin[eV]:   52.6080, T[K]:  301.48, density[g/cm3]: 1.114, calc_time[min]:   868.71
    Dyn  step:    86700, time[ps]:   86.70, Etot[eV]:  -6586.6915, Epot[eV]:  -6639.8164, Ekin[eV]:   53.1249, T[K]:  304.44, density[g/cm3]: 1.120, calc_time[min]:   869.57
    Dyn  step:    86800, time[ps]:   86.80, Etot[eV]:  -6585.8640, Epot[eV]:  -6637.6499, Ekin[eV]:   51.7859, T[K]:  296.77, density[g/cm3]: 1.121, calc_time[min]:   870.45
    Dyn  step:    86900, time[ps]:   86.90, Etot[eV]:  -6586.6499, Epot[eV]:  -6638.5340, Ekin[eV]:   51.8841, T[K]:  297.33, density[g/cm3]: 1.121, calc_time[min]:   871.31
    Dyn  step:    87000, time[ps]:   87.00, Etot[eV]:  -6586.3104, Epot[eV]:  -6638.4906, Ekin[eV]:   52.1802, T[K]:  299.03, density[g/cm3]: 1.120, calc_time[min]:   872.15
    Dyn  step:    87100, time[ps]:   87.10, Etot[eV]:  -6586.2859, Epot[eV]:  -6638.8445, Ekin[eV]:   52.5586, T[K]:  301.19, density[g/cm3]: 1.110, calc_time[min]:   873.01
    Dyn  step:    87200, time[ps]:   87.20, Etot[eV]:  -6586.2058, Epot[eV]:  -6637.9288, Ekin[eV]:   51.7230, T[K]:  296.41, density[g/cm3]: 1.109, calc_time[min]:   873.84
    Dyn  step:    87300, time[ps]:   87.30, Etot[eV]:  -6586.1209, Epot[eV]:  -6638.2476, Ekin[eV]:   52.1267, T[K]:  298.72, density[g/cm3]: 1.108, calc_time[min]:   874.71
    Dyn  step:    87400, time[ps]:   87.40, Etot[eV]:  -6586.2559, Epot[eV]:  -6639.0869, Ekin[eV]:   52.8310, T[K]:  302.75, density[g/cm3]: 1.111, calc_time[min]:   875.58
    Dyn  step:    87500, time[ps]:   87.50, Etot[eV]:  -6586.1077, Epot[eV]:  -6637.7326, Ekin[eV]:   51.6249, T[K]:  295.84, density[g/cm3]: 1.111, calc_time[min]:   876.44
    Dyn  step:    87600, time[ps]:   87.60, Etot[eV]:  -6586.1918, Epot[eV]:  -6639.2868, Ekin[eV]:   53.0949, T[K]:  304.27, density[g/cm3]: 1.110, calc_time[min]:   877.27
    Dyn  step:    87700, time[ps]:   87.70, Etot[eV]:  -6586.4367, Epot[eV]:  -6638.7137, Ekin[eV]:   52.2770, T[K]:  299.58, density[g/cm3]: 1.115, calc_time[min]:   878.13
    Dyn  step:    87800, time[ps]:   87.80, Etot[eV]:  -6585.7218, Epot[eV]:  -6637.1428, Ekin[eV]:   51.4210, T[K]:  294.67, density[g/cm3]: 1.113, calc_time[min]:   879.05
    Dyn  step:    87900, time[ps]:   87.90, Etot[eV]:  -6585.8030, Epot[eV]:  -6637.0337, Ekin[eV]:   51.2307, T[K]:  293.58, density[g/cm3]: 1.116, calc_time[min]:   879.96
    Dyn  step:    88000, time[ps]:   88.00, Etot[eV]:  -6585.8404, Epot[eV]:  -6638.2954, Ekin[eV]:   52.4550, T[K]:  300.60, density[g/cm3]: 1.118, calc_time[min]:   880.86
    Dyn  step:    88100, time[ps]:   88.10, Etot[eV]:  -6585.8687, Epot[eV]:  -6637.0220, Ekin[eV]:   51.1532, T[K]:  293.14, density[g/cm3]: 1.121, calc_time[min]:   881.78
    Dyn  step:    88200, time[ps]:   88.20, Etot[eV]:  -6585.6041, Epot[eV]:  -6638.3818, Ekin[eV]:   52.7777, T[K]:  302.45, density[g/cm3]: 1.122, calc_time[min]:   882.68
    Dyn  step:    88300, time[ps]:   88.30, Etot[eV]:  -6586.2451, Epot[eV]:  -6638.9511, Ekin[eV]:   52.7059, T[K]:  302.04, density[g/cm3]: 1.118, calc_time[min]:   883.57
    Dyn  step:    88400, time[ps]:   88.40, Etot[eV]:  -6585.6197, Epot[eV]:  -6637.5086, Ekin[eV]:   51.8889, T[K]:  297.36, density[g/cm3]: 1.119, calc_time[min]:   884.48
    Dyn  step:    88500, time[ps]:   88.50, Etot[eV]:  -6585.9575, Epot[eV]:  -6639.2033, Ekin[eV]:   53.2458, T[K]:  305.13, density[g/cm3]: 1.112, calc_time[min]:   885.39
    Dyn  step:    88600, time[ps]:   88.60, Etot[eV]:  -6586.5379, Epot[eV]:  -6638.5680, Ekin[eV]:   52.0301, T[K]:  298.17, density[g/cm3]: 1.111, calc_time[min]:   886.29
    Dyn  step:    88700, time[ps]:   88.70, Etot[eV]:  -6585.8830, Epot[eV]:  -6637.9053, Ekin[eV]:   52.0224, T[K]:  298.12, density[g/cm3]: 1.113, calc_time[min]:   887.20
    Dyn  step:    88800, time[ps]:   88.80, Etot[eV]:  -6586.2573, Epot[eV]:  -6638.9974, Ekin[eV]:   52.7401, T[K]:  302.23, density[g/cm3]: 1.117, calc_time[min]:   888.12
    Dyn  step:    88900, time[ps]:   88.90, Etot[eV]:  -6586.5100, Epot[eV]:  -6638.2465, Ekin[eV]:   51.7364, T[K]:  296.48, density[g/cm3]: 1.117, calc_time[min]:   889.03
    Dyn  step:    89000, time[ps]:   89.00, Etot[eV]:  -6585.7712, Epot[eV]:  -6638.8952, Ekin[eV]:   53.1240, T[K]:  304.43, density[g/cm3]: 1.111, calc_time[min]:   889.92
    Dyn  step:    89100, time[ps]:   89.10, Etot[eV]:  -6586.3904, Epot[eV]:  -6638.4324, Ekin[eV]:   52.0420, T[K]:  298.23, density[g/cm3]: 1.113, calc_time[min]:   890.85
    Dyn  step:    89200, time[ps]:   89.20, Etot[eV]:  -6586.6911, Epot[eV]:  -6638.5365, Ekin[eV]:   51.8454, T[K]:  297.11, density[g/cm3]: 1.112, calc_time[min]:   891.76
    Dyn  step:    89300, time[ps]:   89.30, Etot[eV]:  -6586.3659, Epot[eV]:  -6639.4653, Ekin[eV]:   53.0994, T[K]:  304.29, density[g/cm3]: 1.108, calc_time[min]:   892.64
    Dyn  step:    89400, time[ps]:   89.40, Etot[eV]:  -6586.1378, Epot[eV]:  -6636.8924, Ekin[eV]:   50.7547, T[K]:  290.86, density[g/cm3]: 1.109, calc_time[min]:   893.57
    Dyn  step:    89500, time[ps]:   89.50, Etot[eV]:  -6585.8280, Epot[eV]:  -6637.1929, Ekin[eV]:   51.3649, T[K]:  294.35, density[g/cm3]: 1.109, calc_time[min]:   894.47
    Dyn  step:    89600, time[ps]:   89.60, Etot[eV]:  -6585.1656, Epot[eV]:  -6636.6673, Ekin[eV]:   51.5017, T[K]:  295.14, density[g/cm3]: 1.104, calc_time[min]:   895.39
    Dyn  step:    89700, time[ps]:   89.70, Etot[eV]:  -6586.0761, Epot[eV]:  -6638.0225, Ekin[eV]:   51.9464, T[K]:  297.69, density[g/cm3]: 1.108, calc_time[min]:   896.31
    Dyn  step:    89800, time[ps]:   89.80, Etot[eV]:  -6586.5290, Epot[eV]:  -6639.5753, Ekin[eV]:   53.0463, T[K]:  303.99, density[g/cm3]: 1.107, calc_time[min]:   897.23
    Dyn  step:    89900, time[ps]:   89.90, Etot[eV]:  -6586.3063, Epot[eV]:  -6639.4469, Ekin[eV]:   53.1406, T[K]:  304.53, density[g/cm3]: 1.108, calc_time[min]:   898.15
    Dyn  step:    90000, time[ps]:   90.00, Etot[eV]:  -6585.5125, Epot[eV]:  -6637.6302, Ekin[eV]:   52.1177, T[K]:  298.67, density[g/cm3]: 1.111, calc_time[min]:   899.08
    Dyn  step:    90100, time[ps]:   90.10, Etot[eV]:  -6585.1878, Epot[eV]:  -6637.7644, Ekin[eV]:   52.5766, T[K]:  301.30, density[g/cm3]: 1.105, calc_time[min]:   900.01
    Dyn  step:    90200, time[ps]:   90.20, Etot[eV]:  -6586.4788, Epot[eV]:  -6639.4302, Ekin[eV]:   52.9513, T[K]:  303.44, density[g/cm3]: 1.105, calc_time[min]:   900.91
    Dyn  step:    90300, time[ps]:   90.30, Etot[eV]:  -6585.6750, Epot[eV]:  -6637.6303, Ekin[eV]:   51.9552, T[K]:  297.74, density[g/cm3]: 1.111, calc_time[min]:   901.85
    Dyn  step:    90400, time[ps]:   90.40, Etot[eV]:  -6585.2181, Epot[eV]:  -6638.7986, Ekin[eV]:   53.5804, T[K]:  307.05, density[g/cm3]: 1.112, calc_time[min]:   902.75
    Dyn  step:    90500, time[ps]:   90.50, Etot[eV]:  -6585.9362, Epot[eV]:  -6637.9477, Ekin[eV]:   52.0116, T[K]:  298.06, density[g/cm3]: 1.114, calc_time[min]:   903.71
    Dyn  step:    90600, time[ps]:   90.60, Etot[eV]:  -6586.1214, Epot[eV]:  -6638.5921, Ekin[eV]:   52.4708, T[K]:  300.69, density[g/cm3]: 1.115, calc_time[min]:   904.63
    Dyn  step:    90700, time[ps]:   90.70, Etot[eV]:  -6585.9805, Epot[eV]:  -6638.7352, Ekin[eV]:   52.7547, T[K]:  302.32, density[g/cm3]: 1.118, calc_time[min]:   905.57
    Dyn  step:    90800, time[ps]:   90.80, Etot[eV]:  -6585.8989, Epot[eV]:  -6639.4589, Ekin[eV]:   53.5600, T[K]:  306.93, density[g/cm3]: 1.119, calc_time[min]:   906.49
    Dyn  step:    90900, time[ps]:   90.90, Etot[eV]:  -6586.1228, Epot[eV]:  -6639.4551, Ekin[eV]:   53.3323, T[K]:  305.63, density[g/cm3]: 1.110, calc_time[min]:   907.43
    Dyn  step:    91000, time[ps]:   91.00, Etot[eV]:  -6586.6302, Epot[eV]:  -6638.4979, Ekin[eV]:   51.8678, T[K]:  297.23, density[g/cm3]: 1.102, calc_time[min]:   908.34
    Dyn  step:    91100, time[ps]:   91.10, Etot[eV]:  -6585.7099, Epot[eV]:  -6637.6698, Ekin[eV]:   51.9599, T[K]:  297.76, density[g/cm3]: 1.096, calc_time[min]:   909.25
    Dyn  step:    91200, time[ps]:   91.20, Etot[eV]:  -6585.3275, Epot[eV]:  -6638.6518, Ekin[eV]:   53.3244, T[K]:  305.58, density[g/cm3]: 1.094, calc_time[min]:   910.16
    Dyn  step:    91300, time[ps]:   91.30, Etot[eV]:  -6586.2293, Epot[eV]:  -6638.2953, Ekin[eV]:   52.0659, T[K]:  298.37, density[g/cm3]: 1.096, calc_time[min]:   911.07
    Dyn  step:    91400, time[ps]:   91.40, Etot[eV]:  -6586.3299, Epot[eV]:  -6638.5747, Ekin[eV]:   52.2448, T[K]:  299.40, density[g/cm3]: 1.101, calc_time[min]:   911.96
    Dyn  step:    91500, time[ps]:   91.50, Etot[eV]:  -6585.4812, Epot[eV]:  -6638.4175, Ekin[eV]:   52.9363, T[K]:  303.36, density[g/cm3]: 1.100, calc_time[min]:   912.87
    Dyn  step:    91600, time[ps]:   91.60, Etot[eV]:  -6585.5071, Epot[eV]:  -6638.2453, Ekin[eV]:   52.7382, T[K]:  302.22, density[g/cm3]: 1.096, calc_time[min]:   913.79
    Dyn  step:    91700, time[ps]:   91.70, Etot[eV]:  -6585.7547, Epot[eV]:  -6638.3510, Ekin[eV]:   52.5964, T[K]:  301.41, density[g/cm3]: 1.098, calc_time[min]:   914.73
    Dyn  step:    91800, time[ps]:   91.80, Etot[eV]:  -6585.1403, Epot[eV]:  -6636.7893, Ekin[eV]:   51.6490, T[K]:  295.98, density[g/cm3]: 1.100, calc_time[min]:   915.70
    Dyn  step:    91900, time[ps]:   91.90, Etot[eV]:  -6586.1289, Epot[eV]:  -6638.3462, Ekin[eV]:   52.2174, T[K]:  299.24, density[g/cm3]: 1.099, calc_time[min]:   916.63
    Dyn  step:    92000, time[ps]:   92.00, Etot[eV]:  -6585.4774, Epot[eV]:  -6638.5435, Ekin[eV]:   53.0661, T[K]:  304.10, density[g/cm3]: 1.098, calc_time[min]:   917.57
    Dyn  step:    92100, time[ps]:   92.10, Etot[eV]:  -6585.8902, Epot[eV]:  -6636.6457, Ekin[eV]:   50.7555, T[K]:  290.86, density[g/cm3]: 1.100, calc_time[min]:   918.52
    Dyn  step:    92200, time[ps]:   92.20, Etot[eV]:  -6585.8079, Epot[eV]:  -6639.4680, Ekin[eV]:   53.6601, T[K]:  307.51, density[g/cm3]: 1.103, calc_time[min]:   919.47
    Dyn  step:    92300, time[ps]:   92.30, Etot[eV]:  -6586.0462, Epot[eV]:  -6638.5565, Ekin[eV]:   52.5103, T[K]:  300.92, density[g/cm3]: 1.101, calc_time[min]:   920.42
    Dyn  step:    92400, time[ps]:   92.40, Etot[eV]:  -6585.8538, Epot[eV]:  -6639.0364, Ekin[eV]:   53.1826, T[K]:  304.77, density[g/cm3]: 1.103, calc_time[min]:   921.40
    Dyn  step:    92500, time[ps]:   92.50, Etot[eV]:  -6586.1916, Epot[eV]:  -6639.3322, Ekin[eV]:   53.1407, T[K]:  304.53, density[g/cm3]: 1.106, calc_time[min]:   922.33
    Dyn  step:    92600, time[ps]:   92.60, Etot[eV]:  -6586.9767, Epot[eV]:  -6638.1326, Ekin[eV]:   51.1559, T[K]:  293.16, density[g/cm3]: 1.108, calc_time[min]:   923.29
    Dyn  step:    92700, time[ps]:   92.70, Etot[eV]:  -6586.7427, Epot[eV]:  -6639.6590, Ekin[eV]:   52.9163, T[K]:  303.24, density[g/cm3]: 1.117, calc_time[min]:   924.20
    Dyn  step:    92800, time[ps]:   92.80, Etot[eV]:  -6586.4063, Epot[eV]:  -6640.2834, Ekin[eV]:   53.8771, T[K]:  308.75, density[g/cm3]: 1.115, calc_time[min]:   925.15
    Dyn  step:    92900, time[ps]:   92.90, Etot[eV]:  -6586.4753, Epot[eV]:  -6637.9591, Ekin[eV]:   51.4837, T[K]:  295.03, density[g/cm3]: 1.114, calc_time[min]:   926.07
    Dyn  step:    93000, time[ps]:   93.00, Etot[eV]:  -6585.9163, Epot[eV]:  -6639.3086, Ekin[eV]:   53.3923, T[K]:  305.97, density[g/cm3]: 1.110, calc_time[min]:   927.01
    Dyn  step:    93100, time[ps]:   93.10, Etot[eV]:  -6586.1912, Epot[eV]:  -6638.0241, Ekin[eV]:   51.8329, T[K]:  297.04, density[g/cm3]: 1.110, calc_time[min]:   927.97
    Dyn  step:    93200, time[ps]:   93.20, Etot[eV]:  -6586.4742, Epot[eV]:  -6638.0084, Ekin[eV]:   51.5342, T[K]:  295.32, density[g/cm3]: 1.117, calc_time[min]:   928.96
    Dyn  step:    93300, time[ps]:   93.30, Etot[eV]:  -6586.2729, Epot[eV]:  -6640.3430, Ekin[eV]:   54.0702, T[K]:  309.86, density[g/cm3]: 1.116, calc_time[min]:   929.95
    Dyn  step:    93400, time[ps]:   93.40, Etot[eV]:  -6586.1956, Epot[eV]:  -6636.8702, Ekin[eV]:   50.6747, T[K]:  290.40, density[g/cm3]: 1.115, calc_time[min]:   930.91
    Dyn  step:    93500, time[ps]:   93.50, Etot[eV]:  -6586.1163, Epot[eV]:  -6638.4202, Ekin[eV]:   52.3039, T[K]:  299.73, density[g/cm3]: 1.120, calc_time[min]:   931.88
    Dyn  step:    93600, time[ps]:   93.60, Etot[eV]:  -6586.2324, Epot[eV]:  -6638.9014, Ekin[eV]:   52.6690, T[K]:  301.83, density[g/cm3]: 1.121, calc_time[min]:   932.85
    Dyn  step:    93700, time[ps]:   93.70, Etot[eV]:  -6586.4548, Epot[eV]:  -6639.4547, Ekin[eV]:   52.9999, T[K]:  303.72, density[g/cm3]: 1.119, calc_time[min]:   933.82
    Dyn  step:    93800, time[ps]:   93.80, Etot[eV]:  -6585.8657, Epot[eV]:  -6639.0260, Ekin[eV]:   53.1603, T[K]:  304.64, density[g/cm3]: 1.116, calc_time[min]:   934.78
    Dyn  step:    93900, time[ps]:   93.90, Etot[eV]:  -6586.4396, Epot[eV]:  -6637.4051, Ekin[eV]:   50.9656, T[K]:  292.06, density[g/cm3]: 1.120, calc_time[min]:   935.79
    Dyn  step:    94000, time[ps]:   94.00, Etot[eV]:  -6585.7618, Epot[eV]:  -6639.0623, Ekin[eV]:   53.3004, T[K]:  305.44, density[g/cm3]: 1.112, calc_time[min]:   936.77
    Dyn  step:    94100, time[ps]:   94.10, Etot[eV]:  -6586.1003, Epot[eV]:  -6638.7918, Ekin[eV]:   52.6915, T[K]:  301.96, density[g/cm3]: 1.111, calc_time[min]:   937.76
    Dyn  step:    94200, time[ps]:   94.20, Etot[eV]:  -6586.3387, Epot[eV]:  -6640.1027, Ekin[eV]:   53.7640, T[K]:  308.10, density[g/cm3]: 1.107, calc_time[min]:   938.72
    Dyn  step:    94300, time[ps]:   94.30, Etot[eV]:  -6586.3380, Epot[eV]:  -6638.6724, Ekin[eV]:   52.3344, T[K]:  299.91, density[g/cm3]: 1.111, calc_time[min]:   939.67
    Dyn  step:    94400, time[ps]:   94.40, Etot[eV]:  -6586.1258, Epot[eV]:  -6638.9861, Ekin[eV]:   52.8603, T[K]:  302.92, density[g/cm3]: 1.111, calc_time[min]:   940.66
    Dyn  step:    94500, time[ps]:   94.50, Etot[eV]:  -6586.0305, Epot[eV]:  -6638.7805, Ekin[eV]:   52.7500, T[K]:  302.29, density[g/cm3]: 1.116, calc_time[min]:   941.66
    Dyn  step:    94600, time[ps]:   94.60, Etot[eV]:  -6585.9507, Epot[eV]:  -6637.7258, Ekin[eV]:   51.7751, T[K]:  296.70, density[g/cm3]: 1.110, calc_time[min]:   942.66
    Dyn  step:    94700, time[ps]:   94.70, Etot[eV]:  -6585.7435, Epot[eV]:  -6637.6979, Ekin[eV]:   51.9543, T[K]:  297.73, density[g/cm3]: 1.112, calc_time[min]:   943.67
    Dyn  step:    94800, time[ps]:   94.80, Etot[eV]:  -6586.0128, Epot[eV]:  -6640.0442, Ekin[eV]:   54.0315, T[K]:  309.63, density[g/cm3]: 1.113, calc_time[min]:   944.67
    Dyn  step:    94900, time[ps]:   94.90, Etot[eV]:  -6585.4607, Epot[eV]:  -6637.4632, Ekin[eV]:   52.0024, T[K]:  298.01, density[g/cm3]: 1.114, calc_time[min]:   945.69
    Dyn  step:    95000, time[ps]:   95.00, Etot[eV]:  -6586.2652, Epot[eV]:  -6638.1260, Ekin[eV]:   51.8608, T[K]:  297.19, density[g/cm3]: 1.113, calc_time[min]:   946.70
    Dyn  step:    95100, time[ps]:   95.10, Etot[eV]:  -6585.8609, Epot[eV]:  -6638.9804, Ekin[eV]:   53.1195, T[K]:  304.41, density[g/cm3]: 1.118, calc_time[min]:   947.67
    Dyn  step:    95200, time[ps]:   95.20, Etot[eV]:  -6586.2103, Epot[eV]:  -6638.9179, Ekin[eV]:   52.7075, T[K]:  302.05, density[g/cm3]: 1.116, calc_time[min]:   948.69
    Dyn  step:    95300, time[ps]:   95.30, Etot[eV]:  -6586.0991, Epot[eV]:  -6637.9780, Ekin[eV]:   51.8789, T[K]:  297.30, density[g/cm3]: 1.113, calc_time[min]:   949.71
    Dyn  step:    95400, time[ps]:   95.40, Etot[eV]:  -6586.2581, Epot[eV]:  -6638.3211, Ekin[eV]:   52.0630, T[K]:  298.35, density[g/cm3]: 1.114, calc_time[min]:   950.71
    Dyn  step:    95500, time[ps]:   95.50, Etot[eV]:  -6586.3438, Epot[eV]:  -6637.6856, Ekin[eV]:   51.3418, T[K]:  294.22, density[g/cm3]: 1.110, calc_time[min]:   951.71
    Dyn  step:    95600, time[ps]:   95.60, Etot[eV]:  -6585.9405, Epot[eV]:  -6637.9111, Ekin[eV]:   51.9706, T[K]:  297.82, density[g/cm3]: 1.111, calc_time[min]:   952.70
    Dyn  step:    95700, time[ps]:   95.70, Etot[eV]:  -6585.7329, Epot[eV]:  -6638.8931, Ekin[eV]:   53.1602, T[K]:  304.64, density[g/cm3]: 1.113, calc_time[min]:   953.72
    Dyn  step:    95800, time[ps]:   95.80, Etot[eV]:  -6586.2252, Epot[eV]:  -6638.8111, Ekin[eV]:   52.5859, T[K]:  301.35, density[g/cm3]: 1.112, calc_time[min]:   954.73
    Dyn  step:    95900, time[ps]:   95.90, Etot[eV]:  -6586.2969, Epot[eV]:  -6638.7537, Ekin[eV]:   52.4568, T[K]:  300.61, density[g/cm3]: 1.116, calc_time[min]:   955.71
    Dyn  step:    96000, time[ps]:   96.00, Etot[eV]:  -6585.7661, Epot[eV]:  -6638.3435, Ekin[eV]:   52.5774, T[K]:  301.30, density[g/cm3]: 1.116, calc_time[min]:   956.69
    Dyn  step:    96100, time[ps]:   96.10, Etot[eV]:  -6586.6692, Epot[eV]:  -6637.6443, Ekin[eV]:   50.9750, T[K]:  292.12, density[g/cm3]: 1.124, calc_time[min]:   957.70
    Dyn  step:    96200, time[ps]:   96.20, Etot[eV]:  -6585.6119, Epot[eV]:  -6637.5430, Ekin[eV]:   51.9311, T[K]:  297.60, density[g/cm3]: 1.121, calc_time[min]:   958.67
    Dyn  step:    96300, time[ps]:   96.30, Etot[eV]:  -6586.3477, Epot[eV]:  -6640.0617, Ekin[eV]:   53.7140, T[K]:  307.81, density[g/cm3]: 1.121, calc_time[min]:   959.65
    Dyn  step:    96400, time[ps]:   96.40, Etot[eV]:  -6585.6267, Epot[eV]:  -6636.8511, Ekin[eV]:   51.2244, T[K]:  293.55, density[g/cm3]: 1.125, calc_time[min]:   960.62
    Dyn  step:    96500, time[ps]:   96.50, Etot[eV]:  -6586.1341, Epot[eV]:  -6637.9278, Ekin[eV]:   51.7937, T[K]:  296.81, density[g/cm3]: 1.122, calc_time[min]:   961.57
    Dyn  step:    96600, time[ps]:   96.60, Etot[eV]:  -6586.0298, Epot[eV]:  -6639.2569, Ekin[eV]:   53.2271, T[K]:  305.02, density[g/cm3]: 1.122, calc_time[min]:   962.52
    Dyn  step:    96700, time[ps]:   96.70, Etot[eV]:  -6586.0865, Epot[eV]:  -6639.1024, Ekin[eV]:   53.0159, T[K]:  303.81, density[g/cm3]: 1.119, calc_time[min]:   963.45
    Dyn  step:    96800, time[ps]:   96.80, Etot[eV]:  -6585.8002, Epot[eV]:  -6637.1015, Ekin[eV]:   51.3013, T[K]:  293.99, density[g/cm3]: 1.117, calc_time[min]:   964.41
    Dyn  step:    96900, time[ps]:   96.90, Etot[eV]:  -6586.2055, Epot[eV]:  -6639.0398, Ekin[eV]:   52.8343, T[K]:  302.77, density[g/cm3]: 1.118, calc_time[min]:   965.36
    Dyn  step:    97000, time[ps]:   97.00, Etot[eV]:  -6586.3315, Epot[eV]:  -6639.2587, Ekin[eV]:   52.9272, T[K]:  303.31, density[g/cm3]: 1.122, calc_time[min]:   966.29
    Dyn  step:    97100, time[ps]:   97.10, Etot[eV]:  -6586.0434, Epot[eV]:  -6637.5312, Ekin[eV]:   51.4878, T[K]:  295.06, density[g/cm3]: 1.112, calc_time[min]:   967.20
    Dyn  step:    97200, time[ps]:   97.20, Etot[eV]:  -6585.5045, Epot[eV]:  -6637.4524, Ekin[eV]:   51.9479, T[K]:  297.69, density[g/cm3]: 1.116, calc_time[min]:   968.11
    Dyn  step:    97300, time[ps]:   97.30, Etot[eV]:  -6585.4131, Epot[eV]:  -6637.7195, Ekin[eV]:   52.3064, T[K]:  299.75, density[g/cm3]: 1.119, calc_time[min]:   969.04
    Dyn  step:    97400, time[ps]:   97.40, Etot[eV]:  -6586.2161, Epot[eV]:  -6638.3140, Ekin[eV]:   52.0979, T[K]:  298.55, density[g/cm3]: 1.114, calc_time[min]:   969.90
    Dyn  step:    97500, time[ps]:   97.50, Etot[eV]:  -6586.2590, Epot[eV]:  -6640.2513, Ekin[eV]:   53.9923, T[K]:  309.41, density[g/cm3]: 1.114, calc_time[min]:   970.75
    Dyn  step:    97600, time[ps]:   97.60, Etot[eV]:  -6585.7558, Epot[eV]:  -6637.6643, Ekin[eV]:   51.9084, T[K]:  297.47, density[g/cm3]: 1.113, calc_time[min]:   971.61
    Dyn  step:    97700, time[ps]:   97.70, Etot[eV]:  -6586.2260, Epot[eV]:  -6638.3713, Ekin[eV]:   52.1453, T[K]:  298.83, density[g/cm3]: 1.111, calc_time[min]:   972.45
    Dyn  step:    97800, time[ps]:   97.80, Etot[eV]:  -6585.6401, Epot[eV]:  -6637.8523, Ekin[eV]:   52.2121, T[K]:  299.21, density[g/cm3]: 1.107, calc_time[min]:   973.28
    Dyn  step:    97900, time[ps]:   97.90, Etot[eV]:  -6585.9192, Epot[eV]:  -6639.8597, Ekin[eV]:   53.9405, T[K]:  309.11, density[g/cm3]: 1.108, calc_time[min]:   974.15
    Dyn  step:    98000, time[ps]:   98.00, Etot[eV]:  -6586.2920, Epot[eV]:  -6638.1778, Ekin[eV]:   51.8858, T[K]:  297.34, density[g/cm3]: 1.108, calc_time[min]:   974.97
    Dyn  step:    98100, time[ps]:   98.10, Etot[eV]:  -6586.0234, Epot[eV]:  -6638.3431, Ekin[eV]:   52.3197, T[K]:  299.82, density[g/cm3]: 1.109, calc_time[min]:   975.85
    Dyn  step:    98200, time[ps]:   98.20, Etot[eV]:  -6585.6634, Epot[eV]:  -6638.6019, Ekin[eV]:   52.9385, T[K]:  303.37, density[g/cm3]: 1.107, calc_time[min]:   976.70
    Dyn  step:    98300, time[ps]:   98.30, Etot[eV]:  -6586.3782, Epot[eV]:  -6639.6006, Ekin[eV]:   53.2223, T[K]:  305.00, density[g/cm3]: 1.107, calc_time[min]:   977.54
    Dyn  step:    98400, time[ps]:   98.40, Etot[eV]:  -6585.8270, Epot[eV]:  -6638.1271, Ekin[eV]:   52.3001, T[K]:  299.71, density[g/cm3]: 1.106, calc_time[min]:   978.39
    Dyn  step:    98500, time[ps]:   98.50, Etot[eV]:  -6585.9618, Epot[eV]:  -6638.9339, Ekin[eV]:   52.9721, T[K]:  303.56, density[g/cm3]: 1.109, calc_time[min]:   979.25
    Dyn  step:    98600, time[ps]:   98.60, Etot[eV]:  -6585.9476, Epot[eV]:  -6636.7745, Ekin[eV]:   50.8269, T[K]:  291.27, density[g/cm3]: 1.106, calc_time[min]:   980.13
    Dyn  step:    98700, time[ps]:   98.70, Etot[eV]:  -6585.8946, Epot[eV]:  -6638.2476, Ekin[eV]:   52.3530, T[K]:  300.02, density[g/cm3]: 1.104, calc_time[min]:   981.03
    Dyn  step:    98800, time[ps]:   98.80, Etot[eV]:  -6586.3193, Epot[eV]:  -6639.3454, Ekin[eV]:   53.0261, T[K]:  303.87, density[g/cm3]: 1.108, calc_time[min]:   981.92
    Dyn  step:    98900, time[ps]:   98.90, Etot[eV]:  -6585.9353, Epot[eV]:  -6639.2746, Ekin[eV]:   53.3393, T[K]:  305.67, density[g/cm3]: 1.109, calc_time[min]:   982.82
    Dyn  step:    99000, time[ps]:   99.00, Etot[eV]:  -6585.9497, Epot[eV]:  -6638.9083, Ekin[eV]:   52.9586, T[K]:  303.49, density[g/cm3]: 1.107, calc_time[min]:   983.70
    Dyn  step:    99100, time[ps]:   99.10, Etot[eV]:  -6585.7957, Epot[eV]:  -6637.6487, Ekin[eV]:   51.8530, T[K]:  297.15, density[g/cm3]: 1.110, calc_time[min]:   984.60
    Dyn  step:    99200, time[ps]:   99.20, Etot[eV]:  -6585.8588, Epot[eV]:  -6638.7902, Ekin[eV]:   52.9313, T[K]:  303.33, density[g/cm3]: 1.110, calc_time[min]:   985.50
    Dyn  step:    99300, time[ps]:   99.30, Etot[eV]:  -6585.7681, Epot[eV]:  -6638.7667, Ekin[eV]:   52.9986, T[K]:  303.72, density[g/cm3]: 1.108, calc_time[min]:   986.41
    Dyn  step:    99400, time[ps]:   99.40, Etot[eV]:  -6586.7030, Epot[eV]:  -6638.7093, Ekin[eV]:   52.0063, T[K]:  298.03, density[g/cm3]: 1.111, calc_time[min]:   987.34
    Dyn  step:    99500, time[ps]:   99.50, Etot[eV]:  -6586.6668, Epot[eV]:  -6639.6658, Ekin[eV]:   52.9989, T[K]:  303.72, density[g/cm3]: 1.114, calc_time[min]:   988.25
    Dyn  step:    99600, time[ps]:   99.60, Etot[eV]:  -6586.3641, Epot[eV]:  -6639.8583, Ekin[eV]:   53.4942, T[K]:  306.56, density[g/cm3]: 1.104, calc_time[min]:   989.16
    Dyn  step:    99700, time[ps]:   99.70, Etot[eV]:  -6585.5362, Epot[eV]:  -6638.3249, Ekin[eV]:   52.7887, T[K]:  302.51, density[g/cm3]: 1.099, calc_time[min]:   990.08
    Dyn  step:    99800, time[ps]:   99.80, Etot[eV]:  -6585.6430, Epot[eV]:  -6637.6717, Ekin[eV]:   52.0287, T[K]:  298.16, density[g/cm3]: 1.099, calc_time[min]:   990.96
    Dyn  step:    99900, time[ps]:   99.90, Etot[eV]:  -6586.5814, Epot[eV]:  -6639.9463, Ekin[eV]:   53.3649, T[K]:  305.81, density[g/cm3]: 1.102, calc_time[min]:   991.85
    Dyn  step:   100000, time[ps]:  100.00, Etot[eV]:  -6586.3856, Epot[eV]:  -6638.0475, Ekin[eV]:   51.6619, T[K]:  296.06, density[g/cm3]: 1.101, calc_time[min]:   992.77





    True




```python
## check density
df = pd.read_csv('output/15C39H44O7_npt-eq_v6cU0d3-add.log', sep=',')

## set figure
fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(1, 1, 1)

## set plots
ax1.plot(df['time[ps]'], df['density[g/cm3]'], c='blue', linewidth = 2, label='density')
ax2 = ax1.twinx()
ax2.plot(df['time[ps]'], df['T[K]'], c='green', linewidth = 1, label='T')

## set plot format
ax1.set_xlabel("Time /ps", fontsize=14)
ax1.set_ylabel("Density /g cm$^{-3}$", fontsize=14)
ax1.set_xlim(0,100)
ax1.set_ylim(0.4,1.2)
ax1.tick_params(labelsize=14)
ax2.set_ylabel('T /K', fontsize=14)
ax2.set_ylim(270,330)
ax2.tick_params(labelsize=14)
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()    
ax1.legend(h1+h2, l1+l2, loc='upper left', fontsize=12, frameon=False)

print(f'mean density(80-100 ps) = {df["density[g/cm3]"][-200:].mean():.3f}')
```

    mean density(80-100 ps) = 1.112



    
![png](output_16_1.png)
    


類似化合物の実験値と近い密度~1.1に収束している  
- Bisphenol-A diglycidyl ether  
![image.png](da6b8f0a-f652-4f59-bd6c-9a2dcdf0195d.png)  
density = 1.16 g/mL at 25 °C (lit.)  
https://www.sigmaaldrich.com/JP/ja/product/sigma/d3415  



```python
# NVT-MD (20 ps, Berendsen thermostat at 300 K)
## set parameters NVT-MD
dirout = 'output'
sysname = '15C39H44O7'
fname = f'{sysname}_nvt-eq_{method}'
traj_path = f'{dirout}/{fname}.traj'
log_path = f'{dirout}/{fname}.log'
addlog_path = f'{dirout}/{fname}-add.log'
temp = 300.
dt = 1.0  ## 1 fs
steps = 20000  ## 20 ps
atoms = read('output/15C39H44O7_npt-eq_v6cU0d3.traj', index=-1)
atoms.calc = calculator

## set NVTBerendsen dynamics
dyn = NVTBerendsen(
    atoms=atoms, 
    timestep=dt * units.fs, 
    temperature_K=temp,
    taut=30 * units.fs,
    fixcm=True,
    trajectory=f'{traj_path}', 
    loginterval=100,
    )
dyn.attach(print_dyn, interval=100)
dyn.attach(MDLogger(
            dyn=dyn, atoms=atoms, logfile=f'{log_path}', 
            header=True, stress=True, peratom=False, mode="w"), 
            interval=100
          )
dyn.attach(AddMDLogger(
            dyn=dyn, atoms=atoms, logfile=f'{addlog_path}', 
            header=False, mode="w"),
            interval=100
          ) 

## run MD
Stationary(atoms)
t_start = perf_counter()
dyn.run(steps)
```

    Dyn  step:        0, time[ps]:    0.00, Etot[eV]:  -6586.3856, Epot[eV]:  -6638.0475, Ekin[eV]:   51.6619, T[K]:  296.06, density[g/cm3]: 1.101, calc_time[min]:     0.00
    Dyn  step:      100, time[ps]:    0.10, Etot[eV]:  -6586.2645, Epot[eV]:  -6637.8719, Ekin[eV]:   51.6074, T[K]:  295.74, density[g/cm3]: 1.101, calc_time[min]:     0.45
    Dyn  step:      200, time[ps]:    0.20, Etot[eV]:  -6586.1473, Epot[eV]:  -6638.9645, Ekin[eV]:   52.8173, T[K]:  302.68, density[g/cm3]: 1.101, calc_time[min]:     0.91
    Dyn  step:      300, time[ps]:    0.30, Etot[eV]:  -6586.1753, Epot[eV]:  -6637.5190, Ekin[eV]:   51.3437, T[K]:  294.23, density[g/cm3]: 1.101, calc_time[min]:     1.36
    Dyn  step:      400, time[ps]:    0.40, Etot[eV]:  -6585.8027, Epot[eV]:  -6638.2641, Ekin[eV]:   52.4615, T[K]:  300.64, density[g/cm3]: 1.101, calc_time[min]:     1.83
    Dyn  step:      500, time[ps]:    0.50, Etot[eV]:  -6586.6401, Epot[eV]:  -6639.0046, Ekin[eV]:   52.3645, T[K]:  300.08, density[g/cm3]: 1.101, calc_time[min]:     2.30
    Dyn  step:      600, time[ps]:    0.60, Etot[eV]:  -6586.1295, Epot[eV]:  -6638.2527, Ekin[eV]:   52.1231, T[K]:  298.70, density[g/cm3]: 1.101, calc_time[min]:     2.76
    Dyn  step:      700, time[ps]:    0.70, Etot[eV]:  -6586.1018, Epot[eV]:  -6639.1127, Ekin[eV]:   53.0109, T[K]:  303.79, density[g/cm3]: 1.101, calc_time[min]:     3.24
    Dyn  step:      800, time[ps]:    0.80, Etot[eV]:  -6585.2448, Epot[eV]:  -6637.5022, Ekin[eV]:   52.2574, T[K]:  299.47, density[g/cm3]: 1.101, calc_time[min]:     3.71
    Dyn  step:      900, time[ps]:    0.90, Etot[eV]:  -6585.6470, Epot[eV]:  -6639.2234, Ekin[eV]:   53.5764, T[K]:  307.03, density[g/cm3]: 1.101, calc_time[min]:     4.16
    Dyn  step:     1000, time[ps]:    1.00, Etot[eV]:  -6586.5304, Epot[eV]:  -6638.6721, Ekin[eV]:   52.1417, T[K]:  298.80, density[g/cm3]: 1.101, calc_time[min]:     4.63
    Dyn  step:     1100, time[ps]:    1.10, Etot[eV]:  -6586.1750, Epot[eV]:  -6639.8432, Ekin[eV]:   53.6681, T[K]:  307.55, density[g/cm3]: 1.101, calc_time[min]:     5.09
    Dyn  step:     1200, time[ps]:    1.20, Etot[eV]:  -6585.8639, Epot[eV]:  -6638.1360, Ekin[eV]:   52.2721, T[K]:  299.55, density[g/cm3]: 1.101, calc_time[min]:     5.52
    Dyn  step:     1300, time[ps]:    1.30, Etot[eV]:  -6585.7087, Epot[eV]:  -6637.0735, Ekin[eV]:   51.3648, T[K]:  294.35, density[g/cm3]: 1.101, calc_time[min]:     5.99
    Dyn  step:     1400, time[ps]:    1.40, Etot[eV]:  -6586.1461, Epot[eV]:  -6638.6246, Ekin[eV]:   52.4784, T[K]:  300.73, density[g/cm3]: 1.101, calc_time[min]:     6.43
    Dyn  step:     1500, time[ps]:    1.50, Etot[eV]:  -6585.8395, Epot[eV]:  -6636.7581, Ekin[eV]:   50.9186, T[K]:  291.80, density[g/cm3]: 1.101, calc_time[min]:     6.89
    Dyn  step:     1600, time[ps]:    1.60, Etot[eV]:  -6586.1825, Epot[eV]:  -6639.2821, Ekin[eV]:   53.0996, T[K]:  304.29, density[g/cm3]: 1.101, calc_time[min]:     7.35
    Dyn  step:     1700, time[ps]:    1.70, Etot[eV]:  -6585.4596, Epot[eV]:  -6637.5034, Ekin[eV]:   52.0438, T[K]:  298.24, density[g/cm3]: 1.101, calc_time[min]:     7.83
    Dyn  step:     1800, time[ps]:    1.80, Etot[eV]:  -6585.2793, Epot[eV]:  -6638.1663, Ekin[eV]:   52.8870, T[K]:  303.08, density[g/cm3]: 1.101, calc_time[min]:     8.32
    Dyn  step:     1900, time[ps]:    1.90, Etot[eV]:  -6585.9725, Epot[eV]:  -6639.0177, Ekin[eV]:   53.0452, T[K]:  303.98, density[g/cm3]: 1.101, calc_time[min]:     8.79
    Dyn  step:     2000, time[ps]:    2.00, Etot[eV]:  -6585.8306, Epot[eV]:  -6638.3801, Ekin[eV]:   52.5495, T[K]:  301.14, density[g/cm3]: 1.101, calc_time[min]:     9.23
    Dyn  step:     2100, time[ps]:    2.10, Etot[eV]:  -6586.0114, Epot[eV]:  -6639.1210, Ekin[eV]:   53.1096, T[K]:  304.35, density[g/cm3]: 1.101, calc_time[min]:     9.68
    Dyn  step:     2200, time[ps]:    2.20, Etot[eV]:  -6586.4787, Epot[eV]:  -6638.4717, Ekin[eV]:   51.9929, T[K]:  297.95, density[g/cm3]: 1.101, calc_time[min]:    10.12
    Dyn  step:     2300, time[ps]:    2.30, Etot[eV]:  -6585.9489, Epot[eV]:  -6638.0851, Ekin[eV]:   52.1362, T[K]:  298.77, density[g/cm3]: 1.101, calc_time[min]:    10.57
    Dyn  step:     2400, time[ps]:    2.40, Etot[eV]:  -6585.6530, Epot[eV]:  -6637.8377, Ekin[eV]:   52.1847, T[K]:  299.05, density[g/cm3]: 1.101, calc_time[min]:    11.03
    Dyn  step:     2500, time[ps]:    2.50, Etot[eV]:  -6585.3636, Epot[eV]:  -6639.1684, Ekin[eV]:   53.8048, T[K]:  308.34, density[g/cm3]: 1.101, calc_time[min]:    11.48
    Dyn  step:     2600, time[ps]:    2.60, Etot[eV]:  -6585.9775, Epot[eV]:  -6638.4866, Ekin[eV]:   52.5091, T[K]:  300.91, density[g/cm3]: 1.101, calc_time[min]:    11.95
    Dyn  step:     2700, time[ps]:    2.70, Etot[eV]:  -6586.7201, Epot[eV]:  -6638.4715, Ekin[eV]:   51.7514, T[K]:  296.57, density[g/cm3]: 1.101, calc_time[min]:    12.40
    Dyn  step:     2800, time[ps]:    2.80, Etot[eV]:  -6586.5155, Epot[eV]:  -6640.4794, Ekin[eV]:   53.9639, T[K]:  309.25, density[g/cm3]: 1.101, calc_time[min]:    12.87
    Dyn  step:     2900, time[ps]:    2.90, Etot[eV]:  -6586.4758, Epot[eV]:  -6639.6815, Ekin[eV]:   53.2057, T[K]:  304.90, density[g/cm3]: 1.101, calc_time[min]:    13.32
    Dyn  step:     3000, time[ps]:    3.00, Etot[eV]:  -6586.2395, Epot[eV]:  -6638.1466, Ekin[eV]:   51.9071, T[K]:  297.46, density[g/cm3]: 1.101, calc_time[min]:    13.80
    Dyn  step:     3100, time[ps]:    3.10, Etot[eV]:  -6585.7456, Epot[eV]:  -6637.7837, Ekin[eV]:   52.0381, T[K]:  298.21, density[g/cm3]: 1.101, calc_time[min]:    14.29
    Dyn  step:     3200, time[ps]:    3.20, Etot[eV]:  -6585.7971, Epot[eV]:  -6638.2491, Ekin[eV]:   52.4520, T[K]:  300.58, density[g/cm3]: 1.101, calc_time[min]:    14.74
    Dyn  step:     3300, time[ps]:    3.30, Etot[eV]:  -6586.2582, Epot[eV]:  -6638.0659, Ekin[eV]:   51.8076, T[K]:  296.89, density[g/cm3]: 1.101, calc_time[min]:    15.22
    Dyn  step:     3400, time[ps]:    3.40, Etot[eV]:  -6586.0045, Epot[eV]:  -6637.4281, Ekin[eV]:   51.4236, T[K]:  294.69, density[g/cm3]: 1.101, calc_time[min]:    15.69
    Dyn  step:     3500, time[ps]:    3.50, Etot[eV]:  -6585.8235, Epot[eV]:  -6638.5210, Ekin[eV]:   52.6975, T[K]:  301.99, density[g/cm3]: 1.101, calc_time[min]:    16.15
    Dyn  step:     3600, time[ps]:    3.60, Etot[eV]:  -6586.6051, Epot[eV]:  -6638.3191, Ekin[eV]:   51.7140, T[K]:  296.35, density[g/cm3]: 1.101, calc_time[min]:    16.62
    Dyn  step:     3700, time[ps]:    3.70, Etot[eV]:  -6586.0710, Epot[eV]:  -6637.6797, Ekin[eV]:   51.6088, T[K]:  295.75, density[g/cm3]: 1.101, calc_time[min]:    17.09
    Dyn  step:     3800, time[ps]:    3.80, Etot[eV]:  -6586.2924, Epot[eV]:  -6636.3764, Ekin[eV]:   50.0840, T[K]:  287.01, density[g/cm3]: 1.101, calc_time[min]:    17.56
    Dyn  step:     3900, time[ps]:    3.90, Etot[eV]:  -6585.7353, Epot[eV]:  -6638.9811, Ekin[eV]:   53.2458, T[K]:  305.13, density[g/cm3]: 1.101, calc_time[min]:    18.00
    Dyn  step:     4000, time[ps]:    4.00, Etot[eV]:  -6585.7054, Epot[eV]:  -6639.0848, Ekin[eV]:   53.3794, T[K]:  305.90, density[g/cm3]: 1.101, calc_time[min]:    18.47
    Dyn  step:     4100, time[ps]:    4.10, Etot[eV]:  -6585.8167, Epot[eV]:  -6639.2256, Ekin[eV]:   53.4089, T[K]:  306.07, density[g/cm3]: 1.101, calc_time[min]:    18.93
    Dyn  step:     4200, time[ps]:    4.20, Etot[eV]:  -6586.3322, Epot[eV]:  -6639.0149, Ekin[eV]:   52.6826, T[K]:  301.90, density[g/cm3]: 1.101, calc_time[min]:    19.36
    Dyn  step:     4300, time[ps]:    4.30, Etot[eV]:  -6586.2522, Epot[eV]:  -6638.3595, Ekin[eV]:   52.1073, T[K]:  298.61, density[g/cm3]: 1.101, calc_time[min]:    19.85
    Dyn  step:     4400, time[ps]:    4.40, Etot[eV]:  -6586.0800, Epot[eV]:  -6638.0997, Ekin[eV]:   52.0197, T[K]:  298.11, density[g/cm3]: 1.101, calc_time[min]:    20.33
    Dyn  step:     4500, time[ps]:    4.50, Etot[eV]:  -6585.8790, Epot[eV]:  -6637.5675, Ekin[eV]:   51.6885, T[K]:  296.21, density[g/cm3]: 1.101, calc_time[min]:    20.81
    Dyn  step:     4600, time[ps]:    4.60, Etot[eV]:  -6585.8419, Epot[eV]:  -6638.1123, Ekin[eV]:   52.2704, T[K]:  299.54, density[g/cm3]: 1.101, calc_time[min]:    21.28
    Dyn  step:     4700, time[ps]:    4.70, Etot[eV]:  -6586.4800, Epot[eV]:  -6639.1194, Ekin[eV]:   52.6394, T[K]:  301.66, density[g/cm3]: 1.101, calc_time[min]:    21.74
    Dyn  step:     4800, time[ps]:    4.80, Etot[eV]:  -6586.4228, Epot[eV]:  -6639.1832, Ekin[eV]:   52.7604, T[K]:  302.35, density[g/cm3]: 1.101, calc_time[min]:    22.19
    Dyn  step:     4900, time[ps]:    4.90, Etot[eV]:  -6586.0577, Epot[eV]:  -6637.8590, Ekin[eV]:   51.8013, T[K]:  296.85, density[g/cm3]: 1.101, calc_time[min]:    22.65
    Dyn  step:     5000, time[ps]:    5.00, Etot[eV]:  -6586.0159, Epot[eV]:  -6637.8212, Ekin[eV]:   51.8053, T[K]:  296.88, density[g/cm3]: 1.101, calc_time[min]:    23.11
    Dyn  step:     5100, time[ps]:    5.10, Etot[eV]:  -6586.3270, Epot[eV]:  -6639.9297, Ekin[eV]:   53.6027, T[K]:  307.18, density[g/cm3]: 1.101, calc_time[min]:    23.56
    Dyn  step:     5200, time[ps]:    5.20, Etot[eV]:  -6586.9414, Epot[eV]:  -6640.0813, Ekin[eV]:   53.1399, T[K]:  304.52, density[g/cm3]: 1.101, calc_time[min]:    24.03
    Dyn  step:     5300, time[ps]:    5.30, Etot[eV]:  -6586.1936, Epot[eV]:  -6640.3683, Ekin[eV]:   54.1747, T[K]:  310.46, density[g/cm3]: 1.101, calc_time[min]:    24.52
    Dyn  step:     5400, time[ps]:    5.40, Etot[eV]:  -6586.0815, Epot[eV]:  -6638.0180, Ekin[eV]:   51.9365, T[K]:  297.63, density[g/cm3]: 1.101, calc_time[min]:    25.02
    Dyn  step:     5500, time[ps]:    5.50, Etot[eV]:  -6586.0136, Epot[eV]:  -6638.7813, Ekin[eV]:   52.7677, T[K]:  302.39, density[g/cm3]: 1.101, calc_time[min]:    25.50
    Dyn  step:     5600, time[ps]:    5.60, Etot[eV]:  -6586.2763, Epot[eV]:  -6639.5168, Ekin[eV]:   53.2405, T[K]:  305.10, density[g/cm3]: 1.101, calc_time[min]:    25.98
    Dyn  step:     5700, time[ps]:    5.70, Etot[eV]:  -6585.8903, Epot[eV]:  -6637.1929, Ekin[eV]:   51.3026, T[K]:  294.00, density[g/cm3]: 1.101, calc_time[min]:    26.47
    Dyn  step:     5800, time[ps]:    5.80, Etot[eV]:  -6585.8085, Epot[eV]:  -6638.3174, Ekin[eV]:   52.5089, T[K]:  300.91, density[g/cm3]: 1.101, calc_time[min]:    26.93
    Dyn  step:     5900, time[ps]:    5.90, Etot[eV]:  -6586.1563, Epot[eV]:  -6637.7252, Ekin[eV]:   51.5689, T[K]:  295.52, density[g/cm3]: 1.101, calc_time[min]:    27.41
    Dyn  step:     6000, time[ps]:    6.00, Etot[eV]:  -6587.1807, Epot[eV]:  -6639.5188, Ekin[eV]:   52.3382, T[K]:  299.93, density[g/cm3]: 1.101, calc_time[min]:    27.90
    Dyn  step:     6100, time[ps]:    6.10, Etot[eV]:  -6586.2093, Epot[eV]:  -6637.7839, Ekin[eV]:   51.5746, T[K]:  295.55, density[g/cm3]: 1.101, calc_time[min]:    28.37
    Dyn  step:     6200, time[ps]:    6.20, Etot[eV]:  -6585.4387, Epot[eV]:  -6638.2429, Ekin[eV]:   52.8043, T[K]:  302.60, density[g/cm3]: 1.101, calc_time[min]:    28.85
    Dyn  step:     6300, time[ps]:    6.30, Etot[eV]:  -6586.1718, Epot[eV]:  -6639.7075, Ekin[eV]:   53.5356, T[K]:  306.79, density[g/cm3]: 1.101, calc_time[min]:    29.32
    Dyn  step:     6400, time[ps]:    6.40, Etot[eV]:  -6586.4849, Epot[eV]:  -6638.9896, Ekin[eV]:   52.5048, T[K]:  300.89, density[g/cm3]: 1.101, calc_time[min]:    29.79
    Dyn  step:     6500, time[ps]:    6.50, Etot[eV]:  -6585.0100, Epot[eV]:  -6635.6247, Ekin[eV]:   50.6148, T[K]:  290.05, density[g/cm3]: 1.101, calc_time[min]:    30.27
    Dyn  step:     6600, time[ps]:    6.60, Etot[eV]:  -6585.1212, Epot[eV]:  -6638.1750, Ekin[eV]:   53.0538, T[K]:  304.03, density[g/cm3]: 1.101, calc_time[min]:    30.73
    Dyn  step:     6700, time[ps]:    6.70, Etot[eV]:  -6585.3556, Epot[eV]:  -6637.2520, Ekin[eV]:   51.8964, T[K]:  297.40, density[g/cm3]: 1.101, calc_time[min]:    31.21
    Dyn  step:     6800, time[ps]:    6.80, Etot[eV]:  -6585.8787, Epot[eV]:  -6637.9876, Ekin[eV]:   52.1089, T[K]:  298.62, density[g/cm3]: 1.101, calc_time[min]:    31.68
    Dyn  step:     6900, time[ps]:    6.90, Etot[eV]:  -6585.5010, Epot[eV]:  -6638.3175, Ekin[eV]:   52.8165, T[K]:  302.67, density[g/cm3]: 1.101, calc_time[min]:    32.17
    Dyn  step:     7000, time[ps]:    7.00, Etot[eV]:  -6586.5902, Epot[eV]:  -6638.4866, Ekin[eV]:   51.8963, T[K]:  297.40, density[g/cm3]: 1.101, calc_time[min]:    32.63
    Dyn  step:     7100, time[ps]:    7.10, Etot[eV]:  -6586.4498, Epot[eV]:  -6638.7605, Ekin[eV]:   52.3107, T[K]:  299.77, density[g/cm3]: 1.101, calc_time[min]:    33.10
    Dyn  step:     7200, time[ps]:    7.20, Etot[eV]:  -6586.2560, Epot[eV]:  -6639.3517, Ekin[eV]:   53.0957, T[K]:  304.27, density[g/cm3]: 1.101, calc_time[min]:    33.60
    Dyn  step:     7300, time[ps]:    7.30, Etot[eV]:  -6585.7083, Epot[eV]:  -6637.7752, Ekin[eV]:   52.0669, T[K]:  298.38, density[g/cm3]: 1.101, calc_time[min]:    34.08
    Dyn  step:     7400, time[ps]:    7.40, Etot[eV]:  -6585.5635, Epot[eV]:  -6637.8628, Ekin[eV]:   52.2993, T[K]:  299.71, density[g/cm3]: 1.101, calc_time[min]:    34.56
    Dyn  step:     7500, time[ps]:    7.50, Etot[eV]:  -6586.0407, Epot[eV]:  -6637.3134, Ekin[eV]:   51.2727, T[K]:  293.82, density[g/cm3]: 1.101, calc_time[min]:    35.04
    Dyn  step:     7600, time[ps]:    7.60, Etot[eV]:  -6585.9724, Epot[eV]:  -6637.8420, Ekin[eV]:   51.8696, T[K]:  297.25, density[g/cm3]: 1.101, calc_time[min]:    35.52
    Dyn  step:     7700, time[ps]:    7.70, Etot[eV]:  -6586.3491, Epot[eV]:  -6639.2336, Ekin[eV]:   52.8845, T[K]:  303.06, density[g/cm3]: 1.101, calc_time[min]:    36.02
    Dyn  step:     7800, time[ps]:    7.80, Etot[eV]:  -6585.9281, Epot[eV]:  -6638.8394, Ekin[eV]:   52.9113, T[K]:  303.21, density[g/cm3]: 1.101, calc_time[min]:    36.51
    Dyn  step:     7900, time[ps]:    7.90, Etot[eV]:  -6585.6702, Epot[eV]:  -6636.8695, Ekin[eV]:   51.1992, T[K]:  293.40, density[g/cm3]: 1.101, calc_time[min]:    36.99
    Dyn  step:     8000, time[ps]:    8.00, Etot[eV]:  -6585.8638, Epot[eV]:  -6637.5219, Ekin[eV]:   51.6581, T[K]:  296.03, density[g/cm3]: 1.101, calc_time[min]:    37.47
    Dyn  step:     8100, time[ps]:    8.10, Etot[eV]:  -6586.4908, Epot[eV]:  -6639.8460, Ekin[eV]:   53.3552, T[K]:  305.76, density[g/cm3]: 1.101, calc_time[min]:    37.97
    Dyn  step:     8200, time[ps]:    8.20, Etot[eV]:  -6586.3422, Epot[eV]:  -6638.5625, Ekin[eV]:   52.2203, T[K]:  299.26, density[g/cm3]: 1.101, calc_time[min]:    38.45
    Dyn  step:     8300, time[ps]:    8.30, Etot[eV]:  -6585.8123, Epot[eV]:  -6638.1122, Ekin[eV]:   52.2999, T[K]:  299.71, density[g/cm3]: 1.101, calc_time[min]:    38.94
    Dyn  step:     8400, time[ps]:    8.40, Etot[eV]:  -6585.8062, Epot[eV]:  -6639.6102, Ekin[eV]:   53.8040, T[K]:  308.33, density[g/cm3]: 1.101, calc_time[min]:    39.43
    Dyn  step:     8500, time[ps]:    8.50, Etot[eV]:  -6585.9978, Epot[eV]:  -6637.3263, Ekin[eV]:   51.3285, T[K]:  294.14, density[g/cm3]: 1.101, calc_time[min]:    39.91
    Dyn  step:     8600, time[ps]:    8.60, Etot[eV]:  -6586.1362, Epot[eV]:  -6639.7888, Ekin[eV]:   53.6526, T[K]:  307.46, density[g/cm3]: 1.101, calc_time[min]:    40.39
    Dyn  step:     8700, time[ps]:    8.70, Etot[eV]:  -6586.1565, Epot[eV]:  -6638.5776, Ekin[eV]:   52.4211, T[K]:  300.41, density[g/cm3]: 1.101, calc_time[min]:    40.90
    Dyn  step:     8800, time[ps]:    8.80, Etot[eV]:  -6586.2225, Epot[eV]:  -6637.9945, Ekin[eV]:   51.7719, T[K]:  296.69, density[g/cm3]: 1.101, calc_time[min]:    41.38
    Dyn  step:     8900, time[ps]:    8.90, Etot[eV]:  -6586.2047, Epot[eV]:  -6640.0394, Ekin[eV]:   53.8347, T[K]:  308.51, density[g/cm3]: 1.101, calc_time[min]:    41.88
    Dyn  step:     9000, time[ps]:    9.00, Etot[eV]:  -6585.9230, Epot[eV]:  -6638.3809, Ekin[eV]:   52.4579, T[K]:  300.62, density[g/cm3]: 1.101, calc_time[min]:    42.37
    Dyn  step:     9100, time[ps]:    9.10, Etot[eV]:  -6585.9162, Epot[eV]:  -6637.6415, Ekin[eV]:   51.7253, T[K]:  296.42, density[g/cm3]: 1.101, calc_time[min]:    42.84
    Dyn  step:     9200, time[ps]:    9.20, Etot[eV]:  -6585.7621, Epot[eV]:  -6637.1456, Ekin[eV]:   51.3834, T[K]:  294.46, density[g/cm3]: 1.101, calc_time[min]:    43.32
    Dyn  step:     9300, time[ps]:    9.30, Etot[eV]:  -6585.5407, Epot[eV]:  -6637.8534, Ekin[eV]:   52.3127, T[K]:  299.78, density[g/cm3]: 1.101, calc_time[min]:    43.82
    Dyn  step:     9400, time[ps]:    9.40, Etot[eV]:  -6586.2102, Epot[eV]:  -6638.4312, Ekin[eV]:   52.2210, T[K]:  299.26, density[g/cm3]: 1.101, calc_time[min]:    44.32
    Dyn  step:     9500, time[ps]:    9.50, Etot[eV]:  -6586.7277, Epot[eV]:  -6639.2772, Ekin[eV]:   52.5495, T[K]:  301.14, density[g/cm3]: 1.101, calc_time[min]:    44.79
    Dyn  step:     9600, time[ps]:    9.60, Etot[eV]:  -6585.7353, Epot[eV]:  -6637.6167, Ekin[eV]:   51.8815, T[K]:  297.31, density[g/cm3]: 1.101, calc_time[min]:    45.30
    Dyn  step:     9700, time[ps]:    9.70, Etot[eV]:  -6585.4412, Epot[eV]:  -6638.2298, Ekin[eV]:   52.7886, T[K]:  302.51, density[g/cm3]: 1.101, calc_time[min]:    45.80
    Dyn  step:     9800, time[ps]:    9.80, Etot[eV]:  -6585.2422, Epot[eV]:  -6638.0777, Ekin[eV]:   52.8355, T[K]:  302.78, density[g/cm3]: 1.101, calc_time[min]:    46.29
    Dyn  step:     9900, time[ps]:    9.90, Etot[eV]:  -6586.0376, Epot[eV]:  -6638.8719, Ekin[eV]:   52.8343, T[K]:  302.77, density[g/cm3]: 1.101, calc_time[min]:    46.78
    Dyn  step:    10000, time[ps]:   10.00, Etot[eV]:  -6585.4889, Epot[eV]:  -6638.2583, Ekin[eV]:   52.7694, T[K]:  302.40, density[g/cm3]: 1.101, calc_time[min]:    47.26
    Dyn  step:    10100, time[ps]:   10.10, Etot[eV]:  -6586.3958, Epot[eV]:  -6638.2954, Ekin[eV]:   51.8996, T[K]:  297.42, density[g/cm3]: 1.101, calc_time[min]:    47.76
    Dyn  step:    10200, time[ps]:   10.20, Etot[eV]:  -6586.0648, Epot[eV]:  -6637.8572, Ekin[eV]:   51.7924, T[K]:  296.80, density[g/cm3]: 1.101, calc_time[min]:    48.26
    Dyn  step:    10300, time[ps]:   10.30, Etot[eV]:  -6585.9278, Epot[eV]:  -6638.6296, Ekin[eV]:   52.7019, T[K]:  302.01, density[g/cm3]: 1.101, calc_time[min]:    48.76
    Dyn  step:    10400, time[ps]:   10.40, Etot[eV]:  -6585.6596, Epot[eV]:  -6636.5571, Ekin[eV]:   50.8975, T[K]:  291.67, density[g/cm3]: 1.101, calc_time[min]:    49.26
    Dyn  step:    10500, time[ps]:   10.50, Etot[eV]:  -6586.4514, Epot[eV]:  -6637.3697, Ekin[eV]:   50.9183, T[K]:  291.79, density[g/cm3]: 1.101, calc_time[min]:    49.75
    Dyn  step:    10600, time[ps]:   10.60, Etot[eV]:  -6586.2384, Epot[eV]:  -6638.2920, Ekin[eV]:   52.0536, T[K]:  298.30, density[g/cm3]: 1.101, calc_time[min]:    50.25
    Dyn  step:    10700, time[ps]:   10.70, Etot[eV]:  -6586.0292, Epot[eV]:  -6638.8838, Ekin[eV]:   52.8545, T[K]:  302.89, density[g/cm3]: 1.101, calc_time[min]:    50.76
    Dyn  step:    10800, time[ps]:   10.80, Etot[eV]:  -6586.4110, Epot[eV]:  -6639.4874, Ekin[eV]:   53.0764, T[K]:  304.16, density[g/cm3]: 1.101, calc_time[min]:    51.26
    Dyn  step:    10900, time[ps]:   10.90, Etot[eV]:  -6585.8295, Epot[eV]:  -6637.4006, Ekin[eV]:   51.5711, T[K]:  295.53, density[g/cm3]: 1.101, calc_time[min]:    51.77
    Dyn  step:    11000, time[ps]:   11.00, Etot[eV]:  -6586.0140, Epot[eV]:  -6637.3921, Ekin[eV]:   51.3781, T[K]:  294.43, density[g/cm3]: 1.101, calc_time[min]:    52.26
    Dyn  step:    11100, time[ps]:   11.10, Etot[eV]:  -6585.8867, Epot[eV]:  -6638.1102, Ekin[eV]:   52.2235, T[K]:  299.27, density[g/cm3]: 1.101, calc_time[min]:    52.75
    Dyn  step:    11200, time[ps]:   11.20, Etot[eV]:  -6586.0729, Epot[eV]:  -6639.3840, Ekin[eV]:   53.3111, T[K]:  305.51, density[g/cm3]: 1.101, calc_time[min]:    53.25
    Dyn  step:    11300, time[ps]:   11.30, Etot[eV]:  -6586.1677, Epot[eV]:  -6637.1184, Ekin[eV]:   50.9508, T[K]:  291.98, density[g/cm3]: 1.101, calc_time[min]:    53.76
    Dyn  step:    11400, time[ps]:   11.40, Etot[eV]:  -6585.8074, Epot[eV]:  -6637.9786, Ekin[eV]:   52.1713, T[K]:  298.97, density[g/cm3]: 1.101, calc_time[min]:    54.26
    Dyn  step:    11500, time[ps]:   11.50, Etot[eV]:  -6586.1588, Epot[eV]:  -6638.2617, Ekin[eV]:   52.1030, T[K]:  298.58, density[g/cm3]: 1.101, calc_time[min]:    54.78
    Dyn  step:    11600, time[ps]:   11.60, Etot[eV]:  -6586.6738, Epot[eV]:  -6639.9285, Ekin[eV]:   53.2547, T[K]:  305.18, density[g/cm3]: 1.101, calc_time[min]:    55.27
    Dyn  step:    11700, time[ps]:   11.70, Etot[eV]:  -6586.5999, Epot[eV]:  -6637.9668, Ekin[eV]:   51.3669, T[K]:  294.36, density[g/cm3]: 1.101, calc_time[min]:    55.77
    Dyn  step:    11800, time[ps]:   11.80, Etot[eV]:  -6586.1622, Epot[eV]:  -6638.5331, Ekin[eV]:   52.3709, T[K]:  300.12, density[g/cm3]: 1.101, calc_time[min]:    56.29
    Dyn  step:    11900, time[ps]:   11.90, Etot[eV]:  -6586.0967, Epot[eV]:  -6638.9184, Ekin[eV]:   52.8217, T[K]:  302.70, density[g/cm3]: 1.101, calc_time[min]:    56.76
    Dyn  step:    12000, time[ps]:   12.00, Etot[eV]:  -6586.4301, Epot[eV]:  -6638.7344, Ekin[eV]:   52.3042, T[K]:  299.74, density[g/cm3]: 1.101, calc_time[min]:    57.28
    Dyn  step:    12100, time[ps]:   12.10, Etot[eV]:  -6586.3433, Epot[eV]:  -6638.5921, Ekin[eV]:   52.2487, T[K]:  299.42, density[g/cm3]: 1.101, calc_time[min]:    57.77
    Dyn  step:    12200, time[ps]:   12.20, Etot[eV]:  -6585.7425, Epot[eV]:  -6638.0667, Ekin[eV]:   52.3242, T[K]:  299.85, density[g/cm3]: 1.101, calc_time[min]:    58.27
    Dyn  step:    12300, time[ps]:   12.30, Etot[eV]:  -6586.5099, Epot[eV]:  -6640.8116, Ekin[eV]:   54.3017, T[K]:  311.18, density[g/cm3]: 1.101, calc_time[min]:    58.76
    Dyn  step:    12400, time[ps]:   12.40, Etot[eV]:  -6586.6869, Epot[eV]:  -6639.8530, Ekin[eV]:   53.1662, T[K]:  304.68, density[g/cm3]: 1.101, calc_time[min]:    59.26
    Dyn  step:    12500, time[ps]:   12.50, Etot[eV]:  -6586.8447, Epot[eV]:  -6638.4773, Ekin[eV]:   51.6326, T[K]:  295.89, density[g/cm3]: 1.101, calc_time[min]:    59.76
    Dyn  step:    12600, time[ps]:   12.60, Etot[eV]:  -6586.0466, Epot[eV]:  -6638.0316, Ekin[eV]:   51.9850, T[K]:  297.91, density[g/cm3]: 1.101, calc_time[min]:    60.27
    Dyn  step:    12700, time[ps]:   12.70, Etot[eV]:  -6585.9100, Epot[eV]:  -6638.7875, Ekin[eV]:   52.8775, T[K]:  303.02, density[g/cm3]: 1.101, calc_time[min]:    60.76
    Dyn  step:    12800, time[ps]:   12.80, Etot[eV]:  -6585.5472, Epot[eV]:  -6638.7638, Ekin[eV]:   53.2166, T[K]:  304.96, density[g/cm3]: 1.101, calc_time[min]:    61.24
    Dyn  step:    12900, time[ps]:   12.90, Etot[eV]:  -6585.5061, Epot[eV]:  -6639.5088, Ekin[eV]:   54.0027, T[K]:  309.47, density[g/cm3]: 1.101, calc_time[min]:    61.75
    Dyn  step:    13000, time[ps]:   13.00, Etot[eV]:  -6585.6812, Epot[eV]:  -6637.8207, Ekin[eV]:   52.1395, T[K]:  298.79, density[g/cm3]: 1.101, calc_time[min]:    62.22
    Dyn  step:    13100, time[ps]:   13.10, Etot[eV]:  -6586.3098, Epot[eV]:  -6638.4416, Ekin[eV]:   52.1319, T[K]:  298.75, density[g/cm3]: 1.101, calc_time[min]:    62.70
    Dyn  step:    13200, time[ps]:   13.20, Etot[eV]:  -6586.3669, Epot[eV]:  -6639.5977, Ekin[eV]:   53.2308, T[K]:  305.05, density[g/cm3]: 1.101, calc_time[min]:    63.20
    Dyn  step:    13300, time[ps]:   13.30, Etot[eV]:  -6585.4106, Epot[eV]:  -6637.6303, Ekin[eV]:   52.2197, T[K]:  299.25, density[g/cm3]: 1.101, calc_time[min]:    63.70
    Dyn  step:    13400, time[ps]:   13.40, Etot[eV]:  -6585.3412, Epot[eV]:  -6638.2835, Ekin[eV]:   52.9423, T[K]:  303.39, density[g/cm3]: 1.101, calc_time[min]:    64.19
    Dyn  step:    13500, time[ps]:   13.50, Etot[eV]:  -6585.4705, Epot[eV]:  -6637.3949, Ekin[eV]:   51.9244, T[K]:  297.56, density[g/cm3]: 1.101, calc_time[min]:    64.68
    Dyn  step:    13600, time[ps]:   13.60, Etot[eV]:  -6586.2471, Epot[eV]:  -6639.3456, Ekin[eV]:   53.0985, T[K]:  304.29, density[g/cm3]: 1.101, calc_time[min]:    65.20
    Dyn  step:    13700, time[ps]:   13.70, Etot[eV]:  -6584.9598, Epot[eV]:  -6637.4533, Ekin[eV]:   52.4936, T[K]:  300.82, density[g/cm3]: 1.101, calc_time[min]:    65.70
    Dyn  step:    13800, time[ps]:   13.80, Etot[eV]:  -6585.0439, Epot[eV]:  -6637.7187, Ekin[eV]:   52.6748, T[K]:  301.86, density[g/cm3]: 1.101, calc_time[min]:    66.18
    Dyn  step:    13900, time[ps]:   13.90, Etot[eV]:  -6585.8356, Epot[eV]:  -6637.9206, Ekin[eV]:   52.0849, T[K]:  298.48, density[g/cm3]: 1.101, calc_time[min]:    66.67
    Dyn  step:    14000, time[ps]:   14.00, Etot[eV]:  -6586.0905, Epot[eV]:  -6638.7828, Ekin[eV]:   52.6924, T[K]:  301.96, density[g/cm3]: 1.101, calc_time[min]:    67.16
    Dyn  step:    14100, time[ps]:   14.10, Etot[eV]:  -6586.0686, Epot[eV]:  -6638.1524, Ekin[eV]:   52.0839, T[K]:  298.47, density[g/cm3]: 1.101, calc_time[min]:    67.67
    Dyn  step:    14200, time[ps]:   14.20, Etot[eV]:  -6585.1956, Epot[eV]:  -6636.9399, Ekin[eV]:   51.7443, T[K]:  296.53, density[g/cm3]: 1.101, calc_time[min]:    68.15
    Dyn  step:    14300, time[ps]:   14.30, Etot[eV]:  -6585.6622, Epot[eV]:  -6639.1764, Ekin[eV]:   53.5142, T[K]:  306.67, density[g/cm3]: 1.101, calc_time[min]:    68.64
    Dyn  step:    14400, time[ps]:   14.40, Etot[eV]:  -6586.3243, Epot[eV]:  -6639.3655, Ekin[eV]:   53.0411, T[K]:  303.96, density[g/cm3]: 1.101, calc_time[min]:    69.13
    Dyn  step:    14500, time[ps]:   14.50, Etot[eV]:  -6585.9189, Epot[eV]:  -6636.9237, Ekin[eV]:   51.0048, T[K]:  292.29, density[g/cm3]: 1.101, calc_time[min]:    69.62
    Dyn  step:    14600, time[ps]:   14.60, Etot[eV]:  -6586.1604, Epot[eV]:  -6639.0811, Ekin[eV]:   52.9207, T[K]:  303.27, density[g/cm3]: 1.101, calc_time[min]:    70.09
    Dyn  step:    14700, time[ps]:   14.70, Etot[eV]:  -6585.8228, Epot[eV]:  -6636.9169, Ekin[eV]:   51.0940, T[K]:  292.80, density[g/cm3]: 1.101, calc_time[min]:    70.59
    Dyn  step:    14800, time[ps]:   14.80, Etot[eV]:  -6586.0807, Epot[eV]:  -6638.3407, Ekin[eV]:   52.2600, T[K]:  299.48, density[g/cm3]: 1.101, calc_time[min]:    71.08
    Dyn  step:    14900, time[ps]:   14.90, Etot[eV]:  -6586.1593, Epot[eV]:  -6638.8350, Ekin[eV]:   52.6757, T[K]:  301.86, density[g/cm3]: 1.101, calc_time[min]:    71.58
    Dyn  step:    15000, time[ps]:   15.00, Etot[eV]:  -6586.0977, Epot[eV]:  -6637.6520, Ekin[eV]:   51.5543, T[K]:  295.44, density[g/cm3]: 1.101, calc_time[min]:    72.08
    Dyn  step:    15100, time[ps]:   15.10, Etot[eV]:  -6585.7560, Epot[eV]:  -6639.1303, Ekin[eV]:   53.3742, T[K]:  305.87, density[g/cm3]: 1.101, calc_time[min]:    72.55
    Dyn  step:    15200, time[ps]:   15.20, Etot[eV]:  -6584.9438, Epot[eV]:  -6636.0386, Ekin[eV]:   51.0949, T[K]:  292.81, density[g/cm3]: 1.101, calc_time[min]:    73.06
    Dyn  step:    15300, time[ps]:   15.30, Etot[eV]:  -6586.1203, Epot[eV]:  -6638.6410, Ekin[eV]:   52.5207, T[K]:  300.98, density[g/cm3]: 1.101, calc_time[min]:    73.55
    Dyn  step:    15400, time[ps]:   15.40, Etot[eV]:  -6585.5236, Epot[eV]:  -6637.9692, Ekin[eV]:   52.4455, T[K]:  300.55, density[g/cm3]: 1.101, calc_time[min]:    74.04
    Dyn  step:    15500, time[ps]:   15.50, Etot[eV]:  -6585.6375, Epot[eV]:  -6639.9007, Ekin[eV]:   54.2632, T[K]:  310.96, density[g/cm3]: 1.101, calc_time[min]:    74.49
    Dyn  step:    15600, time[ps]:   15.60, Etot[eV]:  -6586.0544, Epot[eV]:  -6637.5147, Ekin[eV]:   51.4603, T[K]:  294.90, density[g/cm3]: 1.101, calc_time[min]:    74.98
    Dyn  step:    15700, time[ps]:   15.70, Etot[eV]:  -6585.3995, Epot[eV]:  -6636.4651, Ekin[eV]:   51.0656, T[K]:  292.64, density[g/cm3]: 1.101, calc_time[min]:    75.47
    Dyn  step:    15800, time[ps]:   15.80, Etot[eV]:  -6585.4800, Epot[eV]:  -6638.9157, Ekin[eV]:   53.4357, T[K]:  306.22, density[g/cm3]: 1.101, calc_time[min]:    75.94
    Dyn  step:    15900, time[ps]:   15.90, Etot[eV]:  -6585.8291, Epot[eV]:  -6637.6484, Ekin[eV]:   51.8192, T[K]:  296.96, density[g/cm3]: 1.101, calc_time[min]:    76.45
    Dyn  step:    16000, time[ps]:   16.00, Etot[eV]:  -6585.5719, Epot[eV]:  -6638.0097, Ekin[eV]:   52.4377, T[K]:  300.50, density[g/cm3]: 1.101, calc_time[min]:    76.94
    Dyn  step:    16100, time[ps]:   16.10, Etot[eV]:  -6585.0627, Epot[eV]:  -6637.9862, Ekin[eV]:   52.9234, T[K]:  303.28, density[g/cm3]: 1.101, calc_time[min]:    77.45
    Dyn  step:    16200, time[ps]:   16.20, Etot[eV]:  -6585.7768, Epot[eV]:  -6638.2947, Ekin[eV]:   52.5179, T[K]:  300.96, density[g/cm3]: 1.101, calc_time[min]:    77.95
    Dyn  step:    16300, time[ps]:   16.30, Etot[eV]:  -6586.2467, Epot[eV]:  -6638.3225, Ekin[eV]:   52.0758, T[K]:  298.43, density[g/cm3]: 1.101, calc_time[min]:    78.44
    Dyn  step:    16400, time[ps]:   16.40, Etot[eV]:  -6586.2541, Epot[eV]:  -6637.9452, Ekin[eV]:   51.6911, T[K]:  296.22, density[g/cm3]: 1.101, calc_time[min]:    78.94
    Dyn  step:    16500, time[ps]:   16.50, Etot[eV]:  -6585.2888, Epot[eV]:  -6637.7332, Ekin[eV]:   52.4444, T[K]:  300.54, density[g/cm3]: 1.101, calc_time[min]:    79.43
    Dyn  step:    16600, time[ps]:   16.60, Etot[eV]:  -6586.1841, Epot[eV]:  -6638.6214, Ekin[eV]:   52.4373, T[K]:  300.50, density[g/cm3]: 1.101, calc_time[min]:    79.93
    Dyn  step:    16700, time[ps]:   16.70, Etot[eV]:  -6586.0648, Epot[eV]:  -6639.3425, Ekin[eV]:   53.2777, T[K]:  305.31, density[g/cm3]: 1.101, calc_time[min]:    80.42
    Dyn  step:    16800, time[ps]:   16.80, Etot[eV]:  -6585.0810, Epot[eV]:  -6636.5636, Ekin[eV]:   51.4826, T[K]:  295.03, density[g/cm3]: 1.101, calc_time[min]:    80.93
    Dyn  step:    16900, time[ps]:   16.90, Etot[eV]:  -6586.2870, Epot[eV]:  -6638.4506, Ekin[eV]:   52.1636, T[K]:  298.93, density[g/cm3]: 1.101, calc_time[min]:    81.43
    Dyn  step:    17000, time[ps]:   17.00, Etot[eV]:  -6586.0653, Epot[eV]:  -6638.8880, Ekin[eV]:   52.8227, T[K]:  302.71, density[g/cm3]: 1.101, calc_time[min]:    81.91
    Dyn  step:    17100, time[ps]:   17.10, Etot[eV]:  -6586.1440, Epot[eV]:  -6639.0599, Ekin[eV]:   52.9160, T[K]:  303.24, density[g/cm3]: 1.101, calc_time[min]:    82.43
    Dyn  step:    17200, time[ps]:   17.20, Etot[eV]:  -6585.5145, Epot[eV]:  -6638.2910, Ekin[eV]:   52.7765, T[K]:  302.44, density[g/cm3]: 1.101, calc_time[min]:    82.93
    Dyn  step:    17300, time[ps]:   17.30, Etot[eV]:  -6585.1398, Epot[eV]:  -6636.9217, Ekin[eV]:   51.7820, T[K]:  296.74, density[g/cm3]: 1.101, calc_time[min]:    83.43
    Dyn  step:    17400, time[ps]:   17.40, Etot[eV]:  -6585.7331, Epot[eV]:  -6638.1834, Ekin[eV]:   52.4503, T[K]:  300.57, density[g/cm3]: 1.101, calc_time[min]:    83.93
    Dyn  step:    17500, time[ps]:   17.50, Etot[eV]:  -6585.6562, Epot[eV]:  -6639.1262, Ekin[eV]:   53.4700, T[K]:  306.42, density[g/cm3]: 1.101, calc_time[min]:    84.44
    Dyn  step:    17600, time[ps]:   17.60, Etot[eV]:  -6585.8679, Epot[eV]:  -6637.3837, Ekin[eV]:   51.5158, T[K]:  295.22, density[g/cm3]: 1.101, calc_time[min]:    84.93
    Dyn  step:    17700, time[ps]:   17.70, Etot[eV]:  -6585.4451, Epot[eV]:  -6638.5368, Ekin[eV]:   53.0917, T[K]:  304.25, density[g/cm3]: 1.101, calc_time[min]:    85.41
    Dyn  step:    17800, time[ps]:   17.80, Etot[eV]:  -6586.4025, Epot[eV]:  -6638.3919, Ekin[eV]:   51.9894, T[K]:  297.93, density[g/cm3]: 1.101, calc_time[min]:    85.91
    Dyn  step:    17900, time[ps]:   17.90, Etot[eV]:  -6586.8623, Epot[eV]:  -6640.2422, Ekin[eV]:   53.3799, T[K]:  305.90, density[g/cm3]: 1.101, calc_time[min]:    86.42
    Dyn  step:    18000, time[ps]:   18.00, Etot[eV]:  -6585.8592, Epot[eV]:  -6639.2220, Ekin[eV]:   53.3628, T[K]:  305.80, density[g/cm3]: 1.101, calc_time[min]:    86.90
    Dyn  step:    18100, time[ps]:   18.10, Etot[eV]:  -6586.2592, Epot[eV]:  -6639.5702, Ekin[eV]:   53.3110, T[K]:  305.51, density[g/cm3]: 1.101, calc_time[min]:    87.40
    Dyn  step:    18200, time[ps]:   18.20, Etot[eV]:  -6585.7948, Epot[eV]:  -6638.1944, Ekin[eV]:   52.3996, T[K]:  300.28, density[g/cm3]: 1.101, calc_time[min]:    87.89
    Dyn  step:    18300, time[ps]:   18.30, Etot[eV]:  -6585.7396, Epot[eV]:  -6637.7519, Ekin[eV]:   52.0123, T[K]:  298.06, density[g/cm3]: 1.101, calc_time[min]:    88.37
    Dyn  step:    18400, time[ps]:   18.40, Etot[eV]:  -6585.9808, Epot[eV]:  -6639.0840, Ekin[eV]:   53.1033, T[K]:  304.32, density[g/cm3]: 1.101, calc_time[min]:    88.87
    Dyn  step:    18500, time[ps]:   18.50, Etot[eV]:  -6585.8824, Epot[eV]:  -6637.3841, Ekin[eV]:   51.5017, T[K]:  295.14, density[g/cm3]: 1.101, calc_time[min]:    89.35
    Dyn  step:    18600, time[ps]:   18.60, Etot[eV]:  -6585.9610, Epot[eV]:  -6638.8770, Ekin[eV]:   52.9160, T[K]:  303.24, density[g/cm3]: 1.101, calc_time[min]:    89.84
    Dyn  step:    18700, time[ps]:   18.70, Etot[eV]:  -6585.6438, Epot[eV]:  -6637.7917, Ekin[eV]:   52.1480, T[K]:  298.84, density[g/cm3]: 1.101, calc_time[min]:    90.34
    Dyn  step:    18800, time[ps]:   18.80, Etot[eV]:  -6585.9248, Epot[eV]:  -6637.5145, Ekin[eV]:   51.5896, T[K]:  295.64, density[g/cm3]: 1.101, calc_time[min]:    90.85
    Dyn  step:    18900, time[ps]:   18.90, Etot[eV]:  -6585.7954, Epot[eV]:  -6638.6196, Ekin[eV]:   52.8242, T[K]:  302.72, density[g/cm3]: 1.101, calc_time[min]:    91.35
    Dyn  step:    19000, time[ps]:   19.00, Etot[eV]:  -6585.7872, Epot[eV]:  -6638.3115, Ekin[eV]:   52.5243, T[K]:  301.00, density[g/cm3]: 1.101, calc_time[min]:    91.85
    Dyn  step:    19100, time[ps]:   19.10, Etot[eV]:  -6585.2522, Epot[eV]:  -6638.4181, Ekin[eV]:   53.1659, T[K]:  304.67, density[g/cm3]: 1.101, calc_time[min]:    92.35
    Dyn  step:    19200, time[ps]:   19.20, Etot[eV]:  -6586.1333, Epot[eV]:  -6638.6923, Ekin[eV]:   52.5589, T[K]:  301.20, density[g/cm3]: 1.101, calc_time[min]:    92.85
    Dyn  step:    19300, time[ps]:   19.30, Etot[eV]:  -6586.0925, Epot[eV]:  -6638.6557, Ekin[eV]:   52.5632, T[K]:  301.22, density[g/cm3]: 1.101, calc_time[min]:    93.35
    Dyn  step:    19400, time[ps]:   19.40, Etot[eV]:  -6586.2240, Epot[eV]:  -6638.4355, Ekin[eV]:   52.2114, T[K]:  299.20, density[g/cm3]: 1.101, calc_time[min]:    93.85
    Dyn  step:    19500, time[ps]:   19.50, Etot[eV]:  -6586.2346, Epot[eV]:  -6639.5452, Ekin[eV]:   53.3106, T[K]:  305.50, density[g/cm3]: 1.101, calc_time[min]:    94.33
    Dyn  step:    19600, time[ps]:   19.60, Etot[eV]:  -6586.1229, Epot[eV]:  -6638.3111, Ekin[eV]:   52.1882, T[K]:  299.07, density[g/cm3]: 1.101, calc_time[min]:    94.81
    Dyn  step:    19700, time[ps]:   19.70, Etot[eV]:  -6585.6112, Epot[eV]:  -6637.5615, Ekin[eV]:   51.9503, T[K]:  297.71, density[g/cm3]: 1.101, calc_time[min]:    95.30
    Dyn  step:    19800, time[ps]:   19.80, Etot[eV]:  -6585.9581, Epot[eV]:  -6637.4503, Ekin[eV]:   51.4923, T[K]:  295.08, density[g/cm3]: 1.101, calc_time[min]:    95.80
    Dyn  step:    19900, time[ps]:   19.90, Etot[eV]:  -6587.0721, Epot[eV]:  -6638.9866, Ekin[eV]:   51.9145, T[K]:  297.50, density[g/cm3]: 1.101, calc_time[min]:    96.31
    Dyn  step:    20000, time[ps]:   20.00, Etot[eV]:  -6586.6361, Epot[eV]:  -6638.5211, Ekin[eV]:   51.8850, T[K]:  297.33, density[g/cm3]: 1.101, calc_time[min]:    96.83





    True




```python
## check structure
atoms_tmp = read('output/15C39H44O7_nvt-eq_v6cU0d3.traj',index=-1)
get_mol_list(atoms_tmp)[0]
#view_ngl(get_mol_list(atoms_tmp)[0], representations=['ball+stick'], replace_structure=True)
```




    [Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
     Atoms(symbols='C39H44O7', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...)]



平衡化過程で意図しない化学反応は起こっていない  
No undesired chemical reactions occurred in the equilibration process.

## 3. 熱分解シミュレーション
- 昇温速度違い（100, 250, 500 K/ps で 300 K から 2300 K まで）
- 到達温度違い（500 K/ps で 300 K から 2800, 3300, 4300 K まで）


### 昇温速度違い（100, 250, 500 K/ps で 300 K から 2300 K まで）
1) NVT-MD@昇温速度 100 K/ps, 300 K から 2300 K まで, 0.1 fs/step, 25 ps
2) NVT-MD@昇温速度 250 K/ps, 300 K から 2300 K まで, 0.1 fs/step, 25 ps
3) NVT-MD@昇温速度 500 K/ps, 300 K から 2300 K まで, 0.1 fs/step, 25 ps



```python
## set NVT-MD print&amp;logger
class PrintWriteLog(MDExtensionBase):
    def __init__(self, fname, dirout='.', stdout=False):
        self.fname   = fname
        self.dirout  = dirout
        self.t_start = perf_counter()
        self.stdout  = stdout
        
    def __call__(self, system, integrator):
        n_step    = system.current_total_step
        sim_time  = system.current_total_time /1000  ## ps
        E_tot     = system.ase_atoms.get_total_energy()
        E_pot     = system.ase_atoms.get_potential_energy()
        E_kin     = system.ase_atoms.get_kinetic_energy()
        temp      = system.ase_atoms.get_temperature()
        density   = system.ase_atoms.get_masses().sum()/units.mol / (system.ase_atoms.cell.volume*(1e-8**3))
        calc_time = (perf_counter() - self.t_start) / 60.  ## min.

        ## write header
        if n_step == 0:
            hdr  = 'step,time[ps],E_tot[eV],E_pot[eV],E_kin[eV],'
            hdr += 'T[K],density[g/cm3],calc_time[min]'
            with open(f'{self.dirout}/{self.fname}.log', 'w') as f_log:
                f_log.write(f'{hdr}\n')
                
        ## write&amp;print result
        line  = f'{n_step:8d},{sim_time:7.2f},'
        line += f'{E_tot:11.4f},{E_pot:11.4f},{E_kin:9.4f},'
        line += f'{temp:8.2f},{density:7.3f},{calc_time:8.2f}'
        with open(f'{self.dirout}/{self.fname}.log', 'a') as f_log:
            f_log.write(f'{line}\n')
        if self.stdout:
            print(line)
```


```python
## MD with TempSchedular
def pyrolysis_md(md_params):
    ## set pfp_estimator_fn
    estimator_fn = pfp_estimator_fn(
        model_version=md_params['model_version'], 
        calc_mode=md_params['calc_mode'], 
        )

    ## set system
    system = ASEMDSystem(
        atoms=md_params['atoms'],
        step=0,
        time=0.0,
        )

    ## set integrator
    integrator = NVTBerendsenIntegrator(
        timestep=md_params['timestep'],
        temperature=md_params['temp_start'],
        taut=100.,  ## 100 fs (same as the reference)
        fixcm=True,
        )
    
    ## set MD feature
    md = MDFeature(
        integrator=integrator,
        n_run=md_params['n_run'],
        show_progress_bar=False,
        show_logger=False,
        logger_interval=md_params['logger_interval'],
        estimator_fn=estimator_fn,
        traj_file_name=f"{md_params['dirout']}/{md_params['fname']}.traj",
        traj_freq=md_params['traj_freq'],  ## Trajectory saving frequency
        )

    ## set chedular for temperature rising (as md_extensions)
    step_ext = int( (md_params['temp_target']-md_params['temp_start'])/md_params['kelvin_per_ps']*1000/md_params['timestep'] )  ## temporary
    TempSchedular = TemperatureScheduler(
        start_value=md_params['temp_start'], 
        end_value=md_params['temp_target'],
        num_total_steps=step_ext,  ## Number of steps for the linear temperature change
        )

    ## perform MD simulations with previous settings
    md_results = md(
        system=system,
        extensions=[(TempSchedular, 1), 
                    (PrintWriteLog(fname=md_params['fname'],dirout=md_params['dirout']), md_params['logger_interval'])
                   ]
        )
```


```python
## set MD parameter list
dirout = 'output'
sysname = '15C39H44O7'
atoms = read(f'{dirout}/{sysname}_nvt-eq_v6cU0d3.traj', index=-1)
temp_rate_target = [
    [100, 2300],  ## 100 K/ps, 2300 K
    [250, 2300],  ## 250 K/ps, 2300 K
    [500, 2300],  ## 500 K/ps, 2300 K
]
jobparam_list = []
for kelvin_per_ps, temp_target in temp_rate_target:
    md_params = {
        'model_version': 'v6.0.0', 
        'calc_mode': 'PBE_PLUS_D3', 
        "atoms": atoms.copy(),
        "timestep": 0.1,  ## fs
        "temp_start": 300,
        "temp_target": temp_target,
        "n_run": 250000,
        "logger_interval": 500,
        "traj_freq": 500,
        "kelvin_per_ps": kelvin_per_ps,
        "dirout": dirout
        }
    tot_ps = int(md_params['n_run'] * md_params['timestep'] / 1000)
    md_params["fname"] = f'{sysname}_nvt-prod_{temp_target}K-{kelvin_per_ps}Kps-{tot_ps}ps'
    jobparam_list.append(md_params)
jobparam_list
```




    [{'model_version': 'v6.0.0',
      'calc_mode': 'PBE_PLUS_D3',
      'atoms': Atoms(symbols='C585H660O105', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
      'timestep': 0.1,
      'temp_start': 300,
      'temp_target': 2300,
      'n_run': 250000,
      'logger_interval': 500,
      'traj_freq': 500,
      'kelvin_per_ps': 100,
      'dirout': 'output',
      'fname': '15C39H44O7_nvt-prod_2300K-100Kps-25ps'},
     {'model_version': 'v6.0.0',
      'calc_mode': 'PBE_PLUS_D3',
      'atoms': Atoms(symbols='C585H660O105', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
      'timestep': 0.1,
      'temp_start': 300,
      'temp_target': 2300,
      'n_run': 250000,
      'logger_interval': 500,
      'traj_freq': 500,
      'kelvin_per_ps': 250,
      'dirout': 'output',
      'fname': '15C39H44O7_nvt-prod_2300K-250Kps-25ps'},
     {'model_version': 'v6.0.0',
      'calc_mode': 'PBE_PLUS_D3',
      'atoms': Atoms(symbols='C585H660O105', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
      'timestep': 0.1,
      'temp_start': 300,
      'temp_target': 2300,
      'n_run': 250000,
      'logger_interval': 500,
      'traj_freq': 500,
      'kelvin_per_ps': 500,
      'dirout': 'output',
      'fname': '15C39H44O7_nvt-prod_2300K-500Kps-25ps'}]




```python
## run MD (-15 hours each)
Parallel(n_jobs=3, backend="threading", verbose=1)(
    delayed(pyrolysis_md)(md_params) for md_params in jobparam_list
    )
```

    [Parallel(n_jobs=3)]: Using backend ThreadingBackend with 3 concurrent workers.
    The MD trajectory will be saved at /home/jovyan/work/epoxy_pyrolysis/output/15C39H44O7_nvt-prod_2300K-100Kps-25ps.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    The MD trajectory will be saved at /home/jovyan/work/epoxy_pyrolysis/output/15C39H44O7_nvt-prod_2300K-250Kps-25ps.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    The MD trajectory will be saved at /home/jovyan/work/epoxy_pyrolysis/output/15C39H44O7_nvt-prod_2300K-500Kps-25ps.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    [Parallel(n_jobs=3)]: Done   3 out of   3 | elapsed: 1404.8min finished





    [None, None, None]



### 到達温度違い（500 K/ps で 300 K から 2800, 3300, 4300 K まで）
1) NVT-MD@昇温速度 500 K/ps, 300 K から 2800 K まで, 0.1 fs/step, 25 ps
2) NVT-MD@昇温速度 500 K/ps, 300 K から 3300 K まで, 0.1 fs/step, 25 ps
3) NVT-MD@昇温速度 500 K/ps, 300 K から 4300 K まで, 0.1 fs/step, 25 ps



```python
## set MD parameter list
dirout = 'output'
sysname = '15C39H44O7'
atoms = read(f'{dirout}/{sysname}_nvt-eq_v6cU0d3.traj', index=-1)
temp_rate_target = [
    [500, 2800],  ## 500 K/ps, 2800 K
    [500, 3300],  ## 500 K/ps, 3300 K
    [500, 4300],  ## 500 K/ps, 4300 K
]
jobparam_list = []
for kelvin_per_ps, temp_target in temp_rate_target:
    md_params = {
        'model_version': 'v6.0.0', 
        'calc_mode': 'PBE_PLUS_D3', 
        "atoms": atoms.copy(),
        "timestep": 0.1,  ## fs
        "temp_start": 300,
        "temp_target": temp_target,
        "n_run": 250000,
        "logger_interval": 500,
        "traj_freq": 500,
        "kelvin_per_ps": kelvin_per_ps,
        "dirout": dirout
        }
    tot_ps = int(md_params['n_run'] * md_params['timestep'] / 1000)
    md_params["fname"] = f'{sysname}_nvt-prod_{temp_target}K-{kelvin_per_ps}Kps-{tot_ps}ps'
    jobparam_list.append(md_params)
jobparam_list
```




    [{'model_version': 'v6.0.0',
      'calc_mode': 'PBE_PLUS_D3',
      'atoms': Atoms(symbols='C585H660O105', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
      'timestep': 0.1,
      'temp_start': 300,
      'temp_target': 2800,
      'n_run': 250000,
      'logger_interval': 500,
      'traj_freq': 500,
      'kelvin_per_ps': 500,
      'dirout': 'output',
      'fname': '15C39H44O7_nvt-prod_2800K-500Kps-25ps'},
     {'model_version': 'v6.0.0',
      'calc_mode': 'PBE_PLUS_D3',
      'atoms': Atoms(symbols='C585H660O105', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
      'timestep': 0.1,
      'temp_start': 300,
      'temp_target': 3300,
      'n_run': 250000,
      'logger_interval': 500,
      'traj_freq': 500,
      'kelvin_per_ps': 500,
      'dirout': 'output',
      'fname': '15C39H44O7_nvt-prod_3300K-500Kps-25ps'},
     {'model_version': 'v6.0.0',
      'calc_mode': 'PBE_PLUS_D3',
      'atoms': Atoms(symbols='C585H660O105', pbc=True, cell=[31.686357020733265, 21.124238013822207, 21.124238013822207], momenta=..., tags=...),
      'timestep': 0.1,
      'temp_start': 300,
      'temp_target': 4300,
      'n_run': 250000,
      'logger_interval': 500,
      'traj_freq': 500,
      'kelvin_per_ps': 500,
      'dirout': 'output',
      'fname': '15C39H44O7_nvt-prod_4300K-500Kps-25ps'}]




```python
## run MD (-15 hours each)
Parallel(n_jobs=3, backend="threading", verbose=1)(
    delayed(pyrolysis_md)(md_params) for md_params in jobparam_list
    )
```

    [Parallel(n_jobs=3)]: Using backend ThreadingBackend with 3 concurrent workers.
    The MD trajectory will be saved at /home/jovyan/work/epoxy_pyrolysis/output/15C39H44O7_nvt-prod_3300K-500Kps-25ps.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    The MD trajectory will be saved at /home/jovyan/work/epoxy_pyrolysis/output/15C39H44O7_nvt-prod_2800K-500Kps-25ps.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    The MD trajectory will be saved at /home/jovyan/work/epoxy_pyrolysis/output/15C39H44O7_nvt-prod_4300K-500Kps-25ps.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    [Parallel(n_jobs=3)]: Done   3 out of   3 | elapsed: 1422.2min finished





    [None, None, None]



## 4. トラジェクトリー解析
- 熱分解シミュレーションで生成した小分子のフラグメントをカウント("get_mol_list")  
1) 昇温速度違い（元論文 Fig. 4）  
2) 到達温度違い（元論文 Fig. 5）


### 昇温速度違い（元論文 Fig. 4）


```python
## count fragments &amp; return resulting dataframe (get_mol_list)
def count_fragms(l_atoms, ref_fragm):
    for idx, atoms in enumerate(l_atoms):
        d_formula = dict(zip(ref_fragm, np.zeros(len(ref_fragm), dtype=int)))
        d_formula['time[ps]'] = idx * 0.1 * 500 / 1000  ## ps
        d_formula['T[K]'] = atoms.get_temperature()

        l_fragm_atoms = get_mol_list(atoms)[0]
        d_count = Counter([x.get_chemical_formula() for x in l_fragm_atoms])
        for fragm in ref_fragm:
            d_formula[fragm] = d_count[fragm]
            
        if idx == 0:
            df = pd.DataFrame(d_formula, index=[idx])
        else:
            df = pd.concat([df, pd.DataFrame(d_formula, index=[idx])])
    return df
```


```python
## set traj files
l_ftraj_path = [
    'output/15C39H44O7_nvt-prod_2300K-100Kps-25ps.traj',
    'output/15C39H44O7_nvt-prod_2300K-250Kps-25ps.traj',
    'output/15C39H44O7_nvt-prod_2300K-500Kps-25ps.traj',
]
atoms_ini = read('output/15C39H44O7_nvt-eq_v6cU0d3.traj',index=-1)

## set chemical formulas of fragments listed up in the original paper
ref_fragm = ["C39H44O7",  ## original epoxy
           "CH2O",
           "CO",
           "H2O",
           "H2",
           "CH4"
          ]

## get results (~15 min.)
l_rate_df = []
for ftraj in l_ftraj_path:
    l_atoms = [atoms_ini] + read(ftraj, index=':')
    df = count_fragms(l_atoms, ref_fragm)
    l_rate_df.append(df)
```


```python
## plot result
fig = plt.figure(figsize=(9,18))

## set fig
for idx, df in enumerate(l_rate_df,1):
    ax1 = fig.add_subplot(len(l_rate_df),1,idx)
    ax2 = ax1.twinx()

    ## plot
    ax1.plot(df['time[ps]'], df['C39H44O7'], label = "C39H44O7", c='#C65D4D')
    ax1.plot(df['time[ps]'], df['CH2O'], label = "CH2O", c='#E66D0F')
    ax1.plot(df['time[ps]'], df['CO'], label = "CO", c='#729238')
    ax1.plot(df['time[ps]'], df['H2O'], label = "H2O", c='#7230A4')
    ax1.plot(df['time[ps]'], df['H2'], label = "H2", c='#0AB450')
    ax1.plot(df['time[ps]'], df['CH4'], label = "CH4", c='#1E4E7E')
    ax2.plot(df['time[ps]'], df['T[K]'], label = "T[K]", linewidth = 2, c='#8BB3E5')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()    

    ## format
    ax1.set_xlabel("Time /ps", fontsize=16)
    ax1.set_ylabel("Number of fragments", fontsize=16)
    ax2.set_ylabel("T /K", fontsize=16)
    ax1.set_xlim(0,25)
    ax1.set_ylim(0,25)
    ax2.set_ylim(0,2500)
    ax1.legend(h1+h2, l1+l2, loc='upper left', fontsize=14, bbox_to_anchor=(1.14, 1.02))
    ax1.tick_params(labelsize=16)
    ax2.tick_params(labelsize=16)
    
fig.savefig('output/heatingrate_dependency.png', bbox_inches='tight',)
```


    
![png](output_34_0.png)
    


### 到達温度違い（元論文 Fig. 5）


```python
## read structures
l_traj_path = [
    'output/15C39H44O7_nvt-prod_2800K-500Kps-25ps.traj',
    'output/15C39H44O7_nvt-prod_3300K-500Kps-25ps.traj',
    'output/15C39H44O7_nvt-prod_4300K-500Kps-25ps.traj',
]
atoms_ini = read('output/15C39H44O7_nvt-eq_v6cU0d3.traj',index=-1)

## set chemical formulas of fragments listed up in the original paper
ref_fragm = ["C39H44O7",  ## original epoxy
           "CH2O",
           "CO",
           "H2O",
           "H2",
           "CH4"
          ]

## get results (dataframe)
l_temp_df = []
for traj in l_traj_path:
    l_atoms = [atoms_ini] + read(traj, index=':')
    df = count_fragms(l_atoms, ref_fragm)
    l_temp_df.append(df)
```


```python
## plot result
l_ylim = [
    [25, 3000],
    [40, 3500],
    [150, 4500],
]

fig = plt.figure(figsize=(9,18))
idx = 0

## set fig
for ylim, df in zip(l_ylim, l_temp_df):
    idx += 1
    ax1 = fig.add_subplot(len(l_ylim),1,idx)
    ax2 = ax1.twinx()

    ## plot
    ax1.plot(df['time[ps]'], df['C39H44O7'], label = "C39H44O7", c='#C65D4D')
    ax1.plot(df['time[ps]'], df['CH2O'], label = "CH2O", c='#E66D0F')
    ax1.plot(df['time[ps]'], df['CO'], label = "CO", c='#729238')
    ax1.plot(df['time[ps]'], df['H2O'], label = "H2O", c='#7230A4')
    ax1.plot(df['time[ps]'], df['H2'], label = "H2", c='#0AB450')
    ax1.plot(df['time[ps]'], df['CH4'], label = "CH4", c='#1E4E7E')
    ax2.plot(df['time[ps]'], df['T[K]'], label = "T[K]", linewidth = 2, c='#8BB3E5')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()    

    ## format
    ax1.set_xlabel("Time /ps", fontsize=16)
    ax1.set_ylabel("Number of fragments", fontsize=16)
    ax2.set_ylabel("T /K", fontsize=16)
    ax1.set_xlim(0,25)
    ax1.set_ylim(0,ylim[0])
    ax1.set_yticks(np.arange(0, ylim[0]+1, step=ylim[0]/5))
    ax2.set_ylim(0,ylim[1])
    ax2.set_yticks(np.arange(0, ylim[1]+1, step=ylim[1]/5))
    ax1.legend(h1+h2, l1+l2, loc='upper left', fontsize=14, bbox_to_anchor=(1.14, 1.02))
    ax1.tick_params(labelsize=16)
    ax2.tick_params(labelsize=16)

fig.savefig('output/temperature_dependency.png', bbox_inches='tight',)
```


    
![png](output_37_0.png)
    
