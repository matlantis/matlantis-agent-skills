Copyright 2024 ENEOS Corp. and Matlantis Corp. (formerly Preferred Computational Chemistry Inc.) as contributors to Matlantis contrib project

# Fe摺動面におけるリン酸のトライボケミカル反応

以下の系に対して摺動シミュレーションを実施します。
- 原子数：Fe 288原子(スラブ上下)、Trimethyl phosphite 16原子

&lt;img src="input/image.png" alt="代替テキスト" width="500" height="500"&gt;


```python
# 必要なライブラリのインストール。インストール済みであればSkip
!pip install -U pfp-api-client
```

## ライブラリの読み込み


```python
import os, sys, csv, glob, shutil, re, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ase import Atoms
from ase.calculators.dftd3 import DFTD3
from ase.md.langevin import Langevin
from ase.optimize import BFGS, FIRE, LBFGS
from ase import units
from ase.md.npt import NPT
from ase.md.verlet import VelocityVerlet
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.io.dmol import write_dmol_arc
from ase.io.dmol import write_dmol_car
from ase.build import surface, bulk, add_adsorbate
from ase.constraints import FixAtoms, StrainFilter, ExpCellFilter, Filter
from ase.visualize import view

import nglview as nv
from nglview.datafiles import PDB, XTC
from IPython.display import Image, display_png

# PFP
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

estimator = Estimator(model_version="v5.0.0", calc_mode="PBE_PLUS_D3")
calculator = ASECalculator(estimator)

os.makedirs("output", exist_ok = True)
```

## 各種メソッドの定義


```python
def surview (a):
    v = nv.show_ase (a, gui = False)
    v.control.spin ([1, 0, 0], np.pi * 1.4)
    v.add_representation (repr_type = "unitcell",)
    if len (a) &gt; 400:
        v.add_representation (repr_type = "spacefill",)
    v._remote_call ("setSize", args = ["250px", "400px"])
    v.background = '#161616' 
    return v
        
def add_dis (v, a, a1, a2):
    v.shape.add ("cylinder", [a[a1].x, a[a1].y, a[a1].z], [a[a2].x, a[a2].y, a[a2].z], [0.9, 0.1, 0.1], .05)
    x = 0.5 * (a[a1].x + a[a2].x)
    y = 0.5 * (a[a1].y + a[a2].y)
    z = 0.5 * (a[a1].z + a[a2].z)
    d = np.round (a.get_distance (a1, a2), 2)
    v.shape.add ("text", [x, y, z], [0.9, 0.1, 0.1], 3, " " + str (d) + "A")

def add_no (v, a, vec = []):
    if vec == []:
        vec = list (range (len (a)))                 
    for i, atom in enumerate (a):
        if i in vec:
            v.shape.add ("text", [atom.x, atom.y, atom.z], [0.1, 0.5, 0.1], 2, " " + str (i))
    
def add_charge (v, a):
    for i, atom in enumerate (a):
        v.shape.add ("text", [atom.x, atom.y, atom.z], [0.9, 0.1, 0.1], 1, "(" + str (i) + ")" + str (np.round (atom.charge, 3)))
    
def MakeInterface (a, b, h):
    pos1 = a.get_positions ()
    sbl1 = a.get_chemical_symbols ()
    pos2 = b.get_positions () + [0.0, 0.0, h]
    sbl2 = b.get_chemical_symbols ()
    sbl = np.hstack ([sbl1, sbl2])
    pos = np.vstack ([pos1, pos2])
    c = Atoms (symbols = sbl, positions = pos, cell = a.get_cell (), pbc = True)
    return c

def PrintLog ():
    step  = dyn.get_number_of_steps ()
    temp  = atoms.get_temperature ()
    e_pot = atoms.get_potential_energy ()
    e_kin = atoms.get_kinetic_energy ()
    e_tot = atoms.get_total_energy ()
    print ("STEP: %7d  T: %10.2f  PE: %12.5e  KE: %12.5e  TE: %12.5e" % (step, temp, e_pot, e_kin, e_tot))

def WriteFrictionalForces ():
    frc = c3.get_forces ()
    with open ("TMP_fric.dat", "a") as fp:
        np.savetxt (fp, frc)
        fp.close ()
        
class ApplyFriction:
    def __init__ (self, indices = None, mask = None):
        #self.removed_dof = 0
        if indices is None and mask is None:
            raise ValueError ('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValueError ('Use only one of "indices" and "mask".')
        if mask is not None:
            indices = np.arange (len (mask))[np.asarray(mask, bool)]
        else:
            srt = np.sort (indices)
            if (np.diff (srt) == 0).any ():
                raise ValueError (
                    'ApplyFriction: The indices array contained duplicates. '
                    'Perhaps you wanted to specify a mask instead, but '
                    'forgot the mask= keyword.')
                
        self.index = np.asarray (indices, int)

        if self.index.ndim != 1:
            raise ValueError ('Wrong argument to ApplyFriction class!')
            
        self.removed_dof = 3 * len (self.index)
        #pass

    def adjust_positions (self, atoms, newpositions):
        newpositions[self.index, 1] = atoms.positions[self.index, 1] + slid_step

    def adjust_momenta (self, atoms, momenta):
        momenta[self.index, 0] = 0.0

    def adjust_forces (self, atoms, forces):
        forces[self.index, 2] = forces[self.index, 2] + load_step

    def get_removed_dof(self, atoms):
        return 3 * len(self.index)
```

## 吸着分子の構造最適化


```python
mol = read ("./input/TMP.gjf")
mol.set_calculator (calculator)
opt = BFGS (mol)
opt.run (fmax = 0.05)
mol_E = mol.get_potential_energy ()
print ("atom#: {}   PE: {:.4f} eV".format (len (mol), mol_E))
v = surview (mol)
display (v)
```

## bcc-Feのスラブモデルを作成


```python
unitcell = bulk ("Fe", "bcc", a = 2.8664, cubic = True)
unitcell.set_calculator (calculator)
slab = surface (unitcell, (1, 1, 0), layers = 3, vacuum = 20.0)
slab.set_positions (slab.get_positions () - [0.0, 0.0, np.min (slab.get_positions ()[:, 2] - 2.0)])
slab = slab.repeat ([4, 6, 1])
slab.set_pbc (True)
lx,ly,lz = slab.cell.cellpar()[:3]
print ("Number of atoms in slab: {}".format (len (slab)))
print ("X: {:.3f}   Y: {:.3f}   Z: {:.3f}".format (lx, ly, lz))
v = surview (slab)
display (v)
```


```python
# slabの最適化。面方向はバルクの最適化された構造で決まっているので原子座標のみ最適化する。
slab.set_calculator (calculator)
opt = BFGS (slab)
opt.run (fmax = 0.05)
slab_E = slab.get_potential_energy ()
print ("atom#: {}   PE: {:.4f} eV".format (len (slab), slab_E))
```


```python
# 目視による確認
v = surview (slab)
display (v)
```

## 吸着分子をスラブ上に配置


```python
# 吸着分子の配置。
mol_i = mol.copy ()
mol_i.rotate (0.0, "x")
mol_i.rotate (0.0, "y")
b = mol_i.positions[:, 2].argmin ()
ads_i = slab.copy ()
plane_center = [0.5 * lx, 0.5 * ly]
print (plane_center)
add_adsorbate (ads_i, mol_i, 2.5, position = plane_center, mol_index = b)
surview (ads_i)
```


```python
# 配置後に一度構造最適化。極端にEnergyの高い状態などが起こらないようにするため。初期状態次第で時間がかかる可能性あり。
t1=time.time()
c = FixAtoms (mask = [atom.index for atom in ads_i if atom.z &lt;= 7])
ads_i.set_constraint (c)
ads_i.set_calculator (calculator)
opt = BFGS (ads_i)
opt.run (fmax = 0.01)
PotE = ads_i.get_potential_energy ()
print ("atom#: {}   PE: {:.4f} eV".format (len (ads_i), PotE))
write("tmp_on_fe_opt.POSCAR",ads_i)  # 一度ファイルに書き出しておく。必須ではない。
t2=time.time()
print(f"Elapsed time : {t2-t1} sec")
```


```python
ads_i = read("tmp_on_fe_opt.POSCAR")
print ("atom#: {}   AdsE: {:.4f} eV".format (len (ads_i), PotE - (mol_E + slab_E)))
surview (ads_i)
```

## 上部スラブの設定


```python
atoms = MakeInterface (ads_i, slab, 11.0)  # slab下部、上部、slabの相対位置（Å）
v = surview (atoms)
display (v)

from ase.io import write

write("output/slab.xyz", atoms)
```

## MDシミュレーションの準備

設定：slab下部は固定、上部は可動としてコントロールします。


```python
# 各パーツのindexでコントロール
btm = [atom.index for atom in atoms if atom.z &lt;=  2.5]
top = [atom.index for atom in atoms if atom.z &gt;= 15.5]
print ("Fixed atoms:   {:5d}".format (len (btm)))
print ("Sliding atoms: {:5d}".format (len (top)))

#下部は固定
c1 = FixAtoms (btm)
c2 = ApplyFriction (top)
atoms.set_constraint ([c1, c2])
atoms.set_calculator (calculator)
```


```python
# MDの計算パラメータ設定
# 以下、Scriptのテストのため短時間、超高温状態に設定されています

log_interval  = 10     # logの出力頻度
duration_pres = 20000  # 加圧熱平衡状態計算の長さ (MD steps)
duration_slid = 50000  # 摩擦のMD計算の長さ(MD steps)
rescale_intvl = 100    # Interval毎のMD steps
delta_t       = 1.0    # MDのstep size (fs)
temp          = 300.0  # Temperature (K)
slid_vlu      = 10.0   # 水平方向速度 （m/s）
load_vlu      = -2.0   # 垂直方向圧力 （GPa）
```

## 加圧状態での熱平衡状態

鉛直方向に加圧します。


```python
print ("&lt;&lt;&lt; Equilibrate system under pressure &gt;&gt;&gt;")

area = atoms.cell.cellpar()[0] * atoms.cell.cellpar()[1]  # 加圧面の面積の定義
load_step = (load_vlu * 1.0e9) * (area * 1.0e-20) / (1.0 / (units.kJ / 1000.0)) / 1.0e10 / float (len (top))  
slid_step = 0.0

if (os.path.exists ("output/TMP_press.traj") == True): os.remove ("output/TMP_press.traj")
if (os.path.exists ("output/TMP_press.log") == True): os.remove ("output/TMP_press.log")
dyn = VelocityVerlet (atoms, timestep = delta_t * units.fs, logfile = "output/TMP_press.log", loginterval = log_interval)
traj = Trajectory ("output/TMP_press.traj", 'a', atoms)
dyn.attach (PrintLog, interval = log_interval)
dyn.attach (traj.write, interval = log_interval)

# 定期的にMaxwell-Boltzmann速度分布を用いて温度制御、そしてMDの実行
n = 0
while n &lt; int (duration_pres / rescale_intvl):
    MaxwellBoltzmannDistribution (atoms, temp * units.kB)
    dyn.run (rescale_intvl)
    n += 1
```


```python
traj = Trajectory ("output/TMP_press.traj")
v = view (traj, viewer = "ngl")
v.view.add_representation ("ball+stick")
display (v)
```

## 摩擦シミュレーション

鉛直方向に加圧した状態で、水平方向の速度を片側のFeスラブに設定します。


```python
print ("&lt;&lt;&lt; Friction simulation &gt;&gt;&gt;")
area = atoms.cell.cellpar()[0] * atoms.cell.cellpar()[1]
load_step = ((load_vlu - 0.5) * 1.0e9) * (area * 1.0e-20) / (1.0 / (units.kJ / 1000.0)) / 1.0e10 / float (len (top))
slid_step = (slid_vlu * 1.0e10 * 1.0e-15) * (delta_t)

if (os.path.exists ("output/TMP_slid.traj") == True): os.remove ("output/TMP_slid.traj")
if (os.path.exists ("output/TMP_slid.log") == True): os.remove ("output/TMP_slid.log")
dyn = VelocityVerlet (atoms, timestep = delta_t * units.fs, logfile = "TMP_slid.log", loginterval = log_interval)
traj = Trajectory ("output/TMP_slid.traj", 'a', atoms)
dyn.attach (PrintLog, interval = log_interval)
dyn.attach (traj.write, interval = log_interval)

# 定期的にMaxwell-Boltzmann速度分布を用いて温度制御、そしてMDの実行
start_time = time.time ()
n = 0
while n &lt; int (duration_slid / rescale_intvl):
    MaxwellBoltzmannDistribution (atoms, temp * units.kB)
    dyn.run (rescale_intvl)
    n += 1
print ("Elapsed time: %12.6f min" % ((time.time () - start_time) / 60.0))
```


```python
traj = Trajectory ("output/TMP_slid.traj")
v = view (traj, viewer = "ngl")
v.view.add_representation ("ball+stick")
display (v)
```
