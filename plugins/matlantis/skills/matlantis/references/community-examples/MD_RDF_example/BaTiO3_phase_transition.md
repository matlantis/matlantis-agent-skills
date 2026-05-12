Copyright Preferred Computational Chemistry, Inc. as contributors to Matlantis contrib project

# BaTiO3の相転移解析（Tetragonal → Cubic）

ペロブスカイト構造で知られるBaTiO3は約4ooK において結晶対称性がtetragonalからcubicに相転移することが知られています。  
2つの構造の違いを格子定数と動径分布関数（RDF）から比較し、実際にNPTアンサンブルMDで相転移が生じることを確認します。  
Keyword : 相転移、NPT、RDF


```python
!pip install pfp_api_client
!pip install pandas joblib

# # 初回使用時のみ、ライブラリのインストールをお願いします。
```


```python
import os, sys,csv,glob,shutil,re,time
from collections import OrderedDict

# PFP
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator
estimator = Estimator()
calculator = ASECalculator(estimator)

# ASE
import ase
from ase.io import read,write,Trajectory
from ase import Atoms, Atom
from ase.constraints import FixAtoms, FixedPlane, FixBondLengths 
from ase.build import bulk
from ase.visualize import view
from ase.build import surface
from ase.constraints import StrainFilter,ExpCellFilter
from ase.optimize import BFGS,FIRE,MDMin
from ase.geometry.analysis import Analysis
from ase.md.nptberendsen import Inhomogeneous_NPTBerendsen
from ase.md import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase import units


# Other external packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import itertools
import json
from IPython.display import Image, display_png
import nglview as nv
import ipywidgets as widgets
from joblib import Parallel, delayed


try:
    dir_path = Path(__file__).parent
except:
    dir_path = Path("").resolve()

# definition of input and output directries
inp="./input/"
out="./output/"
```


    



```python
# 少しずつ変位幅を減らしながら格子定数最適化を複数回繰り返して構造が真に収束したか確認します
def opt_cell_size(m,sn = 10, iter_count = False): # m:Atomsオブジェクト
    m.set_constraint() # clear constraint
    m.set_calculator(calculator)
    maxf = np.sqrt(((m.get_forces())**2).sum(axis=1).max()) # √(fx^2 + fy^2 + fz^2)の一番大きいものを取得
    ucf = ExpCellFilter(m)
    print("ini   pot:{:.4f},maxforce:{:.4f}".format(m.get_potential_energy(),maxf))
    de = -1 
    s = 1
    ita = 100
    while ( de  &lt; -0.01 or de &gt; 0.01 ) and s &lt;= sn :
        opt = BFGS(ucf,maxstep=0.04*(0.9**s),logfile=None)
        old  =  m.get_potential_energy() 
        opt.run(fmax=0.005,steps =ita)
        maxf = np.sqrt(((m.get_forces())**2).sum(axis=1).max())
        de =  m.get_potential_energy()  - old
        print("{} pot:{:.4f},maxforce:{:.4f},delta:{:.4f}".format(s*ita,m.get_potential_energy(),maxf,de))
        s += 1
    if iter_count == True:
        return m, s*ita
    else:
        return m
```

## 初期構造（Materials Project）の確認

Tetragonal, Cubic 各構造は次の通りMaterials Projectより取得しました。  
Tetragonal : [mp-5986](https://materialsproject.org/materials/mp-5986/)   
Cubic : [mp-2998](https://materialsproject.org/materials/mp-2998/)   

Input cif files are from  
A. Jain*, S.P. Ong*, G. Hautier, W. Chen, W.D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, K.A. Persson (*=equal contributions)  
The Materials Project: A materials genome approach to accelerating materials innovation  
APL Materials, 2013, 1(1), 011002.  
[doi:10.1063/1.4812323](http://dx.doi.org/10.1063/1.4812323)  
[[bibtex]](https://materialsproject.org/static/docs/jain_ong2013.349ca3156250.bib)  
Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)  


```python
# Tetragonal structure
bto_tetra = read(inp+"BaTiO3_mp-5986_computed.cif")
bto_tetra = opt_cell_size(bto_tetra)
v = view(bto_tetra, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```

    ini   pot:-31.1870,maxforce:0.3099
    100 pot:-31.1990,maxforce:0.0006,delta:-0.0120
    200 pot:-31.1990,maxforce:0.0006,delta:-0.0000



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Ti', 'O', 'Ba'), valu…



```python
# Cubic structure
bto_cubic = read(inp+"BaTiO3_mp-2998_computed.cif")
bto_cubic = opt_cell_size(bto_cubic)
v = view(bto_cubic, viewer='ngl')
v.view.add_representation("ball+stick")
display(v)
```

    ini   pot:-31.1468,maxforce:0.0000
    100 pot:-31.1470,maxforce:0.0000,delta:-0.0002



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Ti', 'O', 'Ba'), valu…


### 格子定数と動径分布関数の取得

構造の比較のため、格子定数とRDFを取得します。  
RDFはase.geometry.analysisのAnalysisクラスから取得が可能です。  
使用する [get_rdf 関数](https://wiki.fysik.dtu.dk/ase/ase/geometry.html?highlight=rdf#ase.geometry.analysis.Analysis.get_rdf)の変数は次の通りです。
- rmax: RDFプロット範囲最大値（Å）。この値の2倍がセルサイズを超えてはいけない。
- nbins: ヒストグラムのビン数。出力はこの数の一次元array形式となる。
- elements: 元素種の指定。デフォルトでは全元素を指定。例えば["C", "O"]などのリスト表記などで特定の元素を抽出できる。


```python
bto_list = [bto_cubic, bto_tetra]
tags = ["Cubic", "Tetragonal"]
rmax = 5
nbins = 100
delta = [(rmax/nbins)*i for i in range(nbins)]
for bto, name in zip(bto_list, tags):
    atoms = bto.copy()
    print(name, ":",atoms.cell.get_bravais_lattice())
    atoms *= (4,4,4)
    ana = Analysis(atoms)
    rdf = ana.get_rdf(rmax,nbins,elements=["Ti","O"])[0]
    plt.plot(delta, rdf, label=name)
plt.xlabel("d / Å")
plt.ylabel("RDF")
plt.legend()
plt.show()
```

    Cubic : CUB(a=4.03842)
    Tetragonal : TET(a=3.9859, c=4.28734)



    
![png](output_12_1.png)
    


## MDによる解析

MD計算においては、結晶対称性が変化するためにInhomogeneous NPTアンサンブルを用います。  
基本的な使い方はNPT Berensenと同様です。  
また、制御温度を単調に上昇させるクラスを導入しdynに加えています。  
これにより300K → 600K にシミュレーションを通して徐々に温度が上昇していきます。


```python
def convert_atoms_to_upper(atoms: Atoms) -&gt; Atoms:
    atoms2 = atoms.copy()
    
    # cell "c" -&gt; z-axis
    atoms2.rotate(atoms2.cell[2], (0, 0, 1), rotate_cell=True)

    # cell "b" -&gt; yz-plane
    bx, by, bz = atoms2.cell[1, :]
    angle = 90.0 - np.rad2deg(np.arctan2(by, bx))
    atoms2.rotate(angle, 'z', rotate_cell=True)
    # [Note] cell "a" can be arbitrary.

    # supress numerical precision, lower triangular values must be 0.0.
    atoms2.cell = np.where(np.abs(atoms2.cell) &lt; 1e-10, 0.0, atoms2.cell)

    # Check that cell is uppper triangular form.
    m = atoms2.cell
    assert m[1, 0] == m[2, 0] == m[2, 1] == 0.0, f"cell {m} is not upper triangular!"
    return atoms2
```

### MD計算の実行


```python
#NPT ensamble MD with constant Temperature rising
atoms = bto_tetra.copy()
atoms *=(4,4,4)
atoms = convert_atoms_to_upper(atoms)
atoms.set_calculator(calculator)
os.makedirs(out, exist_ok=True)
tag = out+"BaTiO3_tetragonal_NPT"


mass = 137.33 + 47.867 + 15.999*3
Navo = 6.022
N_atom = 64

dt = 1.0
steps = 100000
temp = 300.0
temp_end = 600.0

# Set the momenta corresponding to T(K).
MaxwellBoltzmannDistribution(atoms, temp * units.kB)
# Sets the center-of-mass momentum to zero.
Stationary(atoms)

# Run MD by NPT
dyn = Inhomogeneous_NPTBerendsen(atoms, timestep=dt * units.fs, temperature_K=temp,
                   taut=100 * units.fs, pressure_au=1.01325 * units.bar,
                   taup=1000 * units.fs, compressibility_au=4.57e-5 / units.bar,
                                logfile=tag + '.log', trajectory=tag + ".traj")


#vol = atoms.get_volume() * 1e-3
rho_volinv = 10 * N_atom * mass / (Navo)
def print_dyn():
    print(f"Dyn  step: {dyn.get_number_of_steps(): &gt;5}, energy: {atoms.get_total_energy():.4f}, density: {rho_volinv/(atoms.get_volume()*1e-3):5.3f}")


dyn.attach(print_dyn, interval=100)

dyn.attach(MDLogger(dyn, atoms, tag+'_md.log', header=True, stress=True,
           peratom=True, mode="w"), interval=100)
class MyMDLogger(MDLogger): 
    def __call__(self):

        dat = f"Dyn  step: {self.dyn.get_number_of_steps(): &gt;5}, "
        dat += f"energy: {self.atoms.get_total_energy():.4f}, "
        dat += f"density: {rho_volinv/(self.atoms.get_volume()*1e-3):5.3f}, "
        dat += f"cell: {self.atoms.cell.cellpar()[0:6]}\n"

        self.logfile.write(dat)
        self.logfile.flush()

dyn.attach(MyMDLogger(dyn, atoms, tag+'_mymdlogger.log', header=False, stress=True,
           peratom=True, mode="w"), interval=100)

class Change_temp:
    def __init__(self, dyn, temp_ini, temp_end, steps):
        self.dyn = dyn
        self.delta = (temp_end - temp_ini)/steps
    def __call__(self):
        self.dyn.temperature += self.delta
    
dyn.attach(Change_temp(dyn, temp, temp_end, steps))

dyn.run(steps)
```

    /home/jovyan/.local/lib/python3.7/site-packages/ase/md/md.py:48: FutureWarning: Specify the temperature in K using the 'temperature_K' argument
      warnings.warn(FutureWarning(w))


    Dyn  step:     0, energy: -1984.5833, density: 5685.069
    Dyn  step:   100, energy: -1979.8720, density: 5611.623
    Dyn  step:   200, energy: -1976.8357, density: 5621.512
    Dyn  step:   300, energy: -1975.1508, density: 5611.820
    Dyn  step:   400, energy: -1974.0454, density: 5629.211
    Dyn  step:   500, energy: -1973.5791, density: 5611.292
    Dyn  step:   600, energy: -1973.0947, density: 5634.116
    Dyn  step:   700, energy: -1972.8615, density: 5627.771
    Dyn  step:   800, energy: -1972.6682, density: 5623.832
    Dyn  step:   900, energy: -1972.5162, density: 5633.815
    Dyn  step:  1000, energy: -1972.4103, density: 5641.118
    Dyn  step:  1100, energy: -1972.4471, density: 5635.657
    Dyn  step:  1200, energy: -1972.2633, density: 5630.431
    Dyn  step:  1300, energy: -1972.3258, density: 5633.563
    Dyn  step:  1400, energy: -1972.3179, density: 5633.698
    Dyn  step:  1500, energy: -1972.2197, density: 5626.144
    Dyn  step:  1600, energy: -1972.1274, density: 5629.216
    Dyn  step:  1700, energy: -1972.3060, density: 5641.174
    Dyn  step:  1800, energy: -1972.1148, density: 5637.335
    Dyn  step:  1900, energy: -1972.2047, density: 5628.483
    Dyn  step:  2000, energy: -1972.1843, density: 5639.580
    Dyn  step:  2100, energy: -1972.0648, density: 5636.023
    Dyn  step:  2200, energy: -1972.0823, density: 5634.102
    Dyn  step:  2300, energy: -1972.0360, density: 5647.178
    Dyn  step:  2400, energy: -1971.9075, density: 5630.785
    Dyn  step:  2500, energy: -1972.0365, density: 5628.827
    Dyn  step:  2600, energy: -1971.9299, density: 5630.841
    Dyn  step:  2700, energy: -1971.9972, density: 5626.643
    Dyn  step:  2800, energy: -1971.8408, density: 5636.213
    Dyn  step:  2900, energy: -1971.8475, density: 5635.641
    Dyn  step:  3000, energy: -1971.7948, density: 5641.946
    Dyn  step:  3100, energy: -1971.7697, density: 5623.109
    Dyn  step:  3200, energy: -1971.7452, density: 5631.597
    Dyn  step:  3300, energy: -1971.6448, density: 5640.789
    Dyn  step:  3400, energy: -1971.7728, density: 5630.121
    Dyn  step:  3500, energy: -1971.7627, density: 5630.397
    Dyn  step:  3600, energy: -1971.8682, density: 5624.891
    Dyn  step:  3700, energy: -1971.7987, density: 5629.595
    Dyn  step:  3800, energy: -1971.7408, density: 5623.187
    Dyn  step:  3900, energy: -1971.7430, density: 5623.645
    Dyn  step:  4000, energy: -1971.6689, density: 5626.383
    Dyn  step:  4100, energy: -1971.7055, density: 5624.457
    Dyn  step:  4200, energy: -1971.7055, density: 5631.332
    Dyn  step:  4300, energy: -1971.5842, density: 5622.702
    Dyn  step:  4400, energy: -1971.5292, density: 5624.467
    Dyn  step:  4500, energy: -1971.4860, density: 5630.226
    Dyn  step:  4600, energy: -1971.3831, density: 5626.817
    Dyn  step:  4700, energy: -1971.4711, density: 5634.701
    Dyn  step:  4800, energy: -1971.4166, density: 5629.169
    Dyn  step:  4900, energy: -1971.3495, density: 5637.161
    Dyn  step:  5000, energy: -1971.2796, density: 5638.844
    Dyn  step:  5100, energy: -1971.3116, density: 5629.042
    Dyn  step:  5200, energy: -1971.1166, density: 5612.650
    Dyn  step:  5300, energy: -1971.2709, density: 5637.264
    Dyn  step:  5400, energy: -1971.2439, density: 5640.497
    Dyn  step:  5500, energy: -1971.2599, density: 5632.473
    Dyn  step:  5600, energy: -1971.2003, density: 5619.332
    Dyn  step:  5700, energy: -1971.2907, density: 5634.236
    Dyn  step:  5800, energy: -1971.0282, density: 5628.495
    Dyn  step:  5900, energy: -1971.2391, density: 5624.361
    Dyn  step:  6000, energy: -1971.0550, density: 5631.540
    Dyn  step:  6100, energy: -1970.9743, density: 5627.834
    Dyn  step:  6200, energy: -1971.0203, density: 5627.817
    Dyn  step:  6300, energy: -1970.9832, density: 5629.626
    Dyn  step:  6400, energy: -1970.9774, density: 5621.822
    Dyn  step:  6500, energy: -1971.0219, density: 5631.742
    Dyn  step:  6600, energy: -1970.9631, density: 5627.681
    Dyn  step:  6700, energy: -1970.9685, density: 5631.456
    Dyn  step:  6800, energy: -1970.9669, density: 5637.363
    Dyn  step:  6900, energy: -1970.9208, density: 5619.277
    Dyn  step:  7000, energy: -1970.8577, density: 5631.396
    Dyn  step:  7100, energy: -1970.6876, density: 5624.114
    Dyn  step:  7200, energy: -1970.7489, density: 5622.149
    Dyn  step:  7300, energy: -1970.7906, density: 5619.585
    Dyn  step:  7400, energy: -1970.7844, density: 5619.740
    Dyn  step:  7500, energy: -1970.7170, density: 5634.527
    Dyn  step:  7600, energy: -1970.6352, density: 5629.484
    Dyn  step:  7700, energy: -1970.6899, density: 5629.939
    Dyn  step:  7800, energy: -1970.6599, density: 5628.248
    Dyn  step:  7900, energy: -1970.7058, density: 5625.427
    Dyn  step:  8000, energy: -1970.5999, density: 5628.478
    Dyn  step:  8100, energy: -1970.6331, density: 5624.969
    Dyn  step:  8200, energy: -1970.5038, density: 5631.063
    Dyn  step:  8300, energy: -1970.5188, density: 5627.904
    Dyn  step:  8400, energy: -1970.4648, density: 5626.835
    Dyn  step:  8500, energy: -1970.4589, density: 5646.239
    Dyn  step:  8600, energy: -1970.4248, density: 5609.977
    Dyn  step:  8700, energy: -1970.3774, density: 5639.122
    Dyn  step:  8800, energy: -1970.4194, density: 5636.914
    Dyn  step:  8900, energy: -1970.4098, density: 5630.888
    Dyn  step:  9000, energy: -1970.4125, density: 5627.427
    Dyn  step:  9100, energy: -1970.2791, density: 5636.618
    Dyn  step:  9200, energy: -1970.4302, density: 5623.578
    Dyn  step:  9300, energy: -1970.1607, density: 5623.231
    Dyn  step:  9400, energy: -1970.2530, density: 5643.373
    Dyn  step:  9500, energy: -1970.0123, density: 5616.182
    Dyn  step:  9600, energy: -1970.0911, density: 5627.264
    Dyn  step:  9700, energy: -1970.1319, density: 5625.186
    Dyn  step:  9800, energy: -1970.1156, density: 5633.361
    Dyn  step:  9900, energy: -1970.1578, density: 5629.877
    Dyn  step: 10000, energy: -1969.9877, density: 5615.370
    Dyn  step: 10100, energy: -1970.1360, density: 5633.292
    Dyn  step: 10200, energy: -1970.0586, density: 5632.272
    Dyn  step: 10300, energy: -1969.9661, density: 5637.749
    Dyn  step: 10400, energy: -1970.0653, density: 5630.162
    Dyn  step: 10500, energy: -1969.9380, density: 5634.058
    Dyn  step: 10600, energy: -1969.9783, density: 5627.967
    Dyn  step: 10700, energy: -1969.9392, density: 5630.495
    Dyn  step: 10800, energy: -1970.0678, density: 5641.976
    Dyn  step: 10900, energy: -1969.8107, density: 5614.123
    Dyn  step: 11000, energy: -1969.8177, density: 5624.168
    Dyn  step: 11100, energy: -1969.7818, density: 5642.931
    Dyn  step: 11200, energy: -1969.8454, density: 5629.640
    Dyn  step: 11300, energy: -1969.7774, density: 5616.619
    Dyn  step: 11400, energy: -1969.7623, density: 5626.879
    Dyn  step: 11500, energy: -1969.7664, density: 5638.839
    Dyn  step: 11600, energy: -1969.5785, density: 5631.289
    Dyn  step: 11700, energy: -1969.7430, density: 5618.598
    Dyn  step: 11800, energy: -1969.6480, density: 5631.592
    Dyn  step: 11900, energy: -1969.5564, density: 5628.749
    Dyn  step: 12000, energy: -1969.6067, density: 5627.174
    Dyn  step: 12100, energy: -1969.5585, density: 5632.684
    Dyn  step: 12200, energy: -1969.4808, density: 5634.062
    Dyn  step: 12300, energy: -1969.4804, density: 5642.279
    Dyn  step: 12400, energy: -1969.4825, density: 5621.507
    Dyn  step: 12500, energy: -1969.4955, density: 5652.399
    Dyn  step: 12600, energy: -1969.3834, density: 5625.390
    Dyn  step: 12700, energy: -1969.4503, density: 5636.670
    Dyn  step: 12800, energy: -1969.3210, density: 5628.541
    Dyn  step: 12900, energy: -1969.3576, density: 5626.314
    Dyn  step: 13000, energy: -1969.3985, density: 5633.259
    Dyn  step: 13100, energy: -1969.1655, density: 5628.694
    Dyn  step: 13200, energy: -1969.2825, density: 5643.519
    Dyn  step: 13300, energy: -1969.1895, density: 5624.480
    Dyn  step: 13400, energy: -1969.2056, density: 5622.189
    Dyn  step: 13500, energy: -1969.2432, density: 5638.464
    Dyn  step: 13600, energy: -1969.2598, density: 5634.771
    Dyn  step: 13700, energy: -1969.1519, density: 5625.585
    Dyn  step: 13800, energy: -1969.2628, density: 5641.234
    Dyn  step: 13900, energy: -1969.0877, density: 5640.294
    Dyn  step: 14000, energy: -1969.1235, density: 5630.289
    Dyn  step: 14100, energy: -1968.8664, density: 5635.936
    Dyn  step: 14200, energy: -1968.9984, density: 5640.966
    Dyn  step: 14300, energy: -1968.9368, density: 5650.054
    Dyn  step: 14400, energy: -1968.9143, density: 5643.156
    Dyn  step: 14500, energy: -1968.8621, density: 5632.215
    Dyn  step: 14600, energy: -1968.7753, density: 5649.217
    Dyn  step: 14700, energy: -1968.7042, density: 5629.989
    Dyn  step: 14800, energy: -1968.8556, density: 5655.410
    Dyn  step: 14900, energy: -1968.6271, density: 5647.561
    Dyn  step: 15000, energy: -1968.6346, density: 5649.217
    Dyn  step: 15100, energy: -1968.7577, density: 5641.867
    Dyn  step: 15200, energy: -1968.5578, density: 5628.360
    Dyn  step: 15300, energy: -1968.7294, density: 5636.337
    Dyn  step: 15400, energy: -1968.5680, density: 5636.964
    Dyn  step: 15500, energy: -1968.3941, density: 5622.121
    Dyn  step: 15600, energy: -1968.5311, density: 5653.347
    Dyn  step: 15700, energy: -1968.5128, density: 5642.632
    Dyn  step: 15800, energy: -1968.4279, density: 5644.505
    Dyn  step: 15900, energy: -1968.4977, density: 5639.297
    Dyn  step: 16000, energy: -1968.5960, density: 5652.614
    Dyn  step: 16100, energy: -1968.6002, density: 5633.919
    Dyn  step: 16200, energy: -1968.5165, density: 5639.656
    Dyn  step: 16300, energy: -1968.4466, density: 5641.627
    Dyn  step: 16400, energy: -1968.3709, density: 5617.382
    Dyn  step: 16500, energy: -1968.2890, density: 5629.061
    Dyn  step: 16600, energy: -1968.3719, density: 5628.590
    Dyn  step: 16700, energy: -1968.3799, density: 5640.457
    Dyn  step: 16800, energy: -1968.2605, density: 5636.114
    Dyn  step: 16900, energy: -1968.2299, density: 5636.003
    Dyn  step: 17000, energy: -1968.2657, density: 5622.960
    Dyn  step: 17100, energy: -1968.1132, density: 5638.145
    Dyn  step: 17200, energy: -1968.1210, density: 5624.833
    Dyn  step: 17300, energy: -1968.2128, density: 5628.226
    Dyn  step: 17400, energy: -1968.0624, density: 5623.431
    Dyn  step: 17500, energy: -1968.0851, density: 5622.863
    Dyn  step: 17600, energy: -1968.0407, density: 5618.135
    Dyn  step: 17700, energy: -1968.1086, density: 5611.256
    Dyn  step: 17800, energy: -1968.0285, density: 5626.740
    Dyn  step: 17900, energy: -1968.0474, density: 5621.443
    Dyn  step: 18000, energy: -1968.0682, density: 5624.706
    Dyn  step: 18100, energy: -1968.0100, density: 5627.849
    Dyn  step: 18200, energy: -1967.7433, density: 5620.892
    Dyn  step: 18300, energy: -1967.9722, density: 5638.477
    Dyn  step: 18400, energy: -1967.8683, density: 5636.491
    Dyn  step: 18500, energy: -1968.1195, density: 5616.898
    Dyn  step: 18600, energy: -1967.7503, density: 5615.729
    Dyn  step: 18700, energy: -1967.9921, density: 5631.510
    Dyn  step: 18800, energy: -1967.8726, density: 5624.537
    Dyn  step: 18900, energy: -1967.8340, density: 5630.603
    Dyn  step: 19000, energy: -1967.8199, density: 5638.192
    Dyn  step: 19100, energy: -1967.5976, density: 5611.096
    Dyn  step: 19200, energy: -1967.6938, density: 5631.188
    Dyn  step: 19300, energy: -1967.7795, density: 5622.374
    Dyn  step: 19400, energy: -1967.7319, density: 5618.065
    Dyn  step: 19500, energy: -1967.5919, density: 5615.065
    Dyn  step: 19600, energy: -1967.5074, density: 5614.725
    Dyn  step: 19700, energy: -1967.7141, density: 5626.620
    Dyn  step: 19800, energy: -1967.4281, density: 5620.062
    Dyn  step: 19900, energy: -1967.7589, density: 5649.414
    Dyn  step: 20000, energy: -1967.5376, density: 5627.598
    Dyn  step: 20100, energy: -1967.6252, density: 5638.147
    Dyn  step: 20200, energy: -1967.3447, density: 5634.994
    Dyn  step: 20300, energy: -1967.3506, density: 5641.779
    Dyn  step: 20400, energy: -1967.1423, density: 5628.569
    Dyn  step: 20500, energy: -1967.2650, density: 5632.889
    Dyn  step: 20600, energy: -1967.2562, density: 5618.480
    Dyn  step: 20700, energy: -1967.2139, density: 5635.566
    Dyn  step: 20800, energy: -1967.3030, density: 5625.738
    Dyn  step: 20900, energy: -1967.2975, density: 5629.654
    Dyn  step: 21000, energy: -1967.1922, density: 5624.409
    Dyn  step: 21100, energy: -1967.3140, density: 5641.026
    Dyn  step: 21200, energy: -1967.1137, density: 5630.055
    Dyn  step: 21300, energy: -1967.1987, density: 5619.910
    Dyn  step: 21400, energy: -1967.0614, density: 5634.414
    Dyn  step: 21500, energy: -1967.1144, density: 5620.659
    Dyn  step: 21600, energy: -1966.9120, density: 5623.246
    Dyn  step: 21700, energy: -1966.9948, density: 5631.184
    Dyn  step: 21800, energy: -1967.1490, density: 5630.278
    Dyn  step: 21900, energy: -1967.0392, density: 5620.611
    Dyn  step: 22000, energy: -1966.9888, density: 5630.798
    Dyn  step: 22100, energy: -1966.8716, density: 5624.796
    Dyn  step: 22200, energy: -1966.9779, density: 5628.861
    Dyn  step: 22300, energy: -1966.8254, density: 5634.355
    Dyn  step: 22400, energy: -1966.9035, density: 5631.017
    Dyn  step: 22500, energy: -1966.9817, density: 5619.024
    Dyn  step: 22600, energy: -1966.8582, density: 5635.089
    Dyn  step: 22700, energy: -1966.8784, density: 5633.971
    Dyn  step: 22800, energy: -1966.7457, density: 5619.288
    Dyn  step: 22900, energy: -1966.7961, density: 5623.806
    Dyn  step: 23000, energy: -1966.8479, density: 5630.474
    Dyn  step: 23100, energy: -1966.7135, density: 5615.735
    Dyn  step: 23200, energy: -1966.7741, density: 5625.372
    Dyn  step: 23300, energy: -1966.8131, density: 5634.331
    Dyn  step: 23400, energy: -1966.7094, density: 5627.695
    Dyn  step: 23500, energy: -1966.5532, density: 5629.976
    Dyn  step: 23600, energy: -1966.6040, density: 5615.431
    Dyn  step: 23700, energy: -1966.4440, density: 5629.144
    Dyn  step: 23800, energy: -1966.6859, density: 5622.409
    Dyn  step: 23900, energy: -1966.4248, density: 5623.950
    Dyn  step: 24000, energy: -1966.6463, density: 5642.406
    Dyn  step: 24100, energy: -1966.3168, density: 5609.534
    Dyn  step: 24200, energy: -1966.5593, density: 5618.583
    Dyn  step: 24300, energy: -1966.5057, density: 5620.905
    Dyn  step: 24400, energy: -1966.4611, density: 5637.791
    Dyn  step: 24500, energy: -1966.4304, density: 5612.807
    Dyn  step: 24600, energy: -1966.4045, density: 5616.785
    Dyn  step: 24700, energy: -1966.3524, density: 5628.436
    Dyn  step: 24800, energy: -1966.2604, density: 5625.687
    Dyn  step: 24900, energy: -1966.3439, density: 5621.772
    Dyn  step: 25000, energy: -1966.1529, density: 5646.839
    Dyn  step: 25100, energy: -1966.3553, density: 5638.054
    Dyn  step: 25200, energy: -1966.2497, density: 5624.951
    Dyn  step: 25300, energy: -1966.1157, density: 5617.948
    Dyn  step: 25400, energy: -1966.2409, density: 5630.914
    Dyn  step: 25500, energy: -1966.1334, density: 5614.228
    Dyn  step: 25600, energy: -1966.2834, density: 5644.170
    Dyn  step: 25700, energy: -1965.9069, density: 5621.000
    Dyn  step: 25800, energy: -1966.0623, density: 5620.283
    Dyn  step: 25900, energy: -1966.1134, density: 5619.071
    Dyn  step: 26000, energy: -1966.1798, density: 5643.965
    Dyn  step: 26100, energy: -1965.9919, density: 5631.229
    Dyn  step: 26200, energy: -1965.9542, density: 5627.043
    Dyn  step: 26300, energy: -1965.9051, density: 5615.990
    Dyn  step: 26400, energy: -1965.9579, density: 5643.389
    Dyn  step: 26500, energy: -1965.8723, density: 5639.377
    Dyn  step: 26600, energy: -1965.9123, density: 5612.724
    Dyn  step: 26700, energy: -1965.8907, density: 5626.690
    Dyn  step: 26800, energy: -1965.8107, density: 5605.055
    Dyn  step: 26900, energy: -1965.7977, density: 5640.859
    Dyn  step: 27000, energy: -1965.5967, density: 5635.194
    Dyn  step: 27100, energy: -1965.7154, density: 5642.295
    Dyn  step: 27200, energy: -1965.4480, density: 5618.222
    Dyn  step: 27300, energy: -1965.6220, density: 5651.799
    Dyn  step: 27400, energy: -1965.6047, density: 5629.241
    Dyn  step: 27500, energy: -1965.6614, density: 5631.555
    Dyn  step: 27600, energy: -1965.4924, density: 5632.963
    Dyn  step: 27700, energy: -1965.7449, density: 5637.820
    Dyn  step: 27800, energy: -1965.4245, density: 5611.980
    Dyn  step: 27900, energy: -1965.7137, density: 5639.703
    Dyn  step: 28000, energy: -1965.2614, density: 5634.239
    Dyn  step: 28100, energy: -1965.5459, density: 5635.103
    Dyn  step: 28200, energy: -1965.3321, density: 5636.389
    Dyn  step: 28300, energy: -1965.4013, density: 5637.571
    Dyn  step: 28400, energy: -1965.0822, density: 5628.150
    Dyn  step: 28500, energy: -1965.1571, density: 5635.178
    Dyn  step: 28600, energy: -1965.2142, density: 5634.111
    Dyn  step: 28700, energy: -1965.1592, density: 5632.995
    Dyn  step: 28800, energy: -1965.1570, density: 5637.163
    Dyn  step: 28900, energy: -1965.3599, density: 5651.683
    Dyn  step: 29000, energy: -1965.0111, density: 5615.746
    Dyn  step: 29100, energy: -1965.0180, density: 5661.429
    Dyn  step: 29200, energy: -1965.1726, density: 5616.037
    Dyn  step: 29300, energy: -1965.0332, density: 5631.290
    Dyn  step: 29400, energy: -1965.1162, density: 5628.957
    Dyn  step: 29500, energy: -1964.9569, density: 5650.244
    Dyn  step: 29600, energy: -1965.1117, density: 5635.549
    Dyn  step: 29700, energy: -1964.9421, density: 5630.111
    Dyn  step: 29800, energy: -1964.8815, density: 5627.638
    Dyn  step: 29900, energy: -1964.8272, density: 5630.867
    Dyn  step: 30000, energy: -1964.7557, density: 5615.481
    Dyn  step: 30100, energy: -1964.7775, density: 5622.574
    Dyn  step: 30200, energy: -1964.8580, density: 5643.222
    Dyn  step: 30300, energy: -1964.6482, density: 5640.437
    Dyn  step: 30400, energy: -1964.6759, density: 5620.934
    Dyn  step: 30500, energy: -1964.6922, density: 5612.838
    Dyn  step: 30600, energy: -1964.6168, density: 5615.321
    Dyn  step: 30700, energy: -1964.6338, density: 5642.293
    Dyn  step: 30800, energy: -1964.6543, density: 5612.522
    Dyn  step: 30900, energy: -1964.6080, density: 5620.362
    Dyn  step: 31000, energy: -1964.6418, density: 5611.993
    Dyn  step: 31100, energy: -1964.5224, density: 5613.130
    Dyn  step: 31200, energy: -1964.6889, density: 5606.151
    Dyn  step: 31300, energy: -1964.7187, density: 5633.896
    Dyn  step: 31400, energy: -1964.5108, density: 5601.741
    Dyn  step: 31500, energy: -1964.6779, density: 5617.964
    Dyn  step: 31600, energy: -1964.8910, density: 5621.770
    Dyn  step: 31700, energy: -1964.5642, density: 5616.150
    Dyn  step: 31800, energy: -1964.6093, density: 5628.973
    Dyn  step: 31900, energy: -1964.5894, density: 5630.796
    Dyn  step: 32000, energy: -1964.4410, density: 5610.995
    Dyn  step: 32100, energy: -1964.3824, density: 5621.181
    Dyn  step: 32200, energy: -1964.5380, density: 5613.917
    Dyn  step: 32300, energy: -1964.2948, density: 5603.685
    Dyn  step: 32400, energy: -1964.4612, density: 5616.645
    Dyn  step: 32500, energy: -1964.3971, density: 5608.129
    Dyn  step: 32600, energy: -1964.2351, density: 5618.184
    Dyn  step: 32700, energy: -1964.2670, density: 5615.669
    Dyn  step: 32800, energy: -1964.1045, density: 5612.272
    Dyn  step: 32900, energy: -1964.0851, density: 5609.806
    Dyn  step: 33000, energy: -1964.1014, density: 5606.705
    Dyn  step: 33100, energy: -1964.1884, density: 5611.762
    Dyn  step: 33200, energy: -1964.0870, density: 5615.607
    Dyn  step: 33300, energy: -1964.1891, density: 5593.268
    Dyn  step: 33400, energy: -1964.1897, density: 5617.146
    Dyn  step: 33500, energy: -1964.2686, density: 5621.315
    Dyn  step: 33600, energy: -1964.0528, density: 5614.968
    Dyn  step: 33700, energy: -1964.2310, density: 5619.242
    Dyn  step: 33800, energy: -1964.0090, density: 5614.071
    Dyn  step: 33900, energy: -1963.9870, density: 5606.472
    Dyn  step: 34000, energy: -1964.0132, density: 5626.536
    Dyn  step: 34100, energy: -1964.1547, density: 5624.291
    Dyn  step: 34200, energy: -1963.9267, density: 5627.268
    Dyn  step: 34300, energy: -1963.9583, density: 5625.553
    Dyn  step: 34400, energy: -1963.8633, density: 5631.442
    Dyn  step: 34500, energy: -1963.8594, density: 5623.343
    Dyn  step: 34600, energy: -1963.8009, density: 5616.684
    Dyn  step: 34700, energy: -1963.8878, density: 5647.633
    Dyn  step: 34800, energy: -1963.7934, density: 5618.892
    Dyn  step: 34900, energy: -1963.6890, density: 5629.047
    Dyn  step: 35000, energy: -1963.6517, density: 5615.166
    Dyn  step: 35100, energy: -1963.5878, density: 5610.958
    Dyn  step: 35200, energy: -1963.8337, density: 5636.038
    Dyn  step: 35300, energy: -1963.6471, density: 5629.932
    Dyn  step: 35400, energy: -1963.6787, density: 5612.386
    Dyn  step: 35500, energy: -1963.6286, density: 5636.300
    Dyn  step: 35600, energy: -1963.6029, density: 5602.254
    Dyn  step: 35700, energy: -1963.4903, density: 5622.214
    Dyn  step: 35800, energy: -1963.6457, density: 5632.556
    Dyn  step: 35900, energy: -1963.5445, density: 5626.605
    Dyn  step: 36000, energy: -1963.3420, density: 5606.037
    Dyn  step: 36100, energy: -1963.5389, density: 5617.948
    Dyn  step: 36200, energy: -1963.2370, density: 5598.886
    Dyn  step: 36300, energy: -1963.6157, density: 5619.192
    Dyn  step: 36400, energy: -1963.5164, density: 5595.794
    Dyn  step: 36500, energy: -1963.3709, density: 5614.058
    Dyn  step: 36600, energy: -1963.4589, density: 5619.882
    Dyn  step: 36700, energy: -1963.2272, density: 5625.710
    Dyn  step: 36800, energy: -1963.2959, density: 5630.085
    Dyn  step: 36900, energy: -1963.1703, density: 5620.585
    Dyn  step: 37000, energy: -1963.2168, density: 5636.656
    Dyn  step: 37100, energy: -1963.0744, density: 5630.840
    Dyn  step: 37200, energy: -1963.0802, density: 5615.204
    Dyn  step: 37300, energy: -1963.0565, density: 5621.030
    Dyn  step: 37400, energy: -1962.9069, density: 5614.933
    Dyn  step: 37500, energy: -1962.6357, density: 5618.061
    Dyn  step: 37600, energy: -1962.8086, density: 5622.085
    Dyn  step: 37700, energy: -1962.6131, density: 5614.234
    Dyn  step: 37800, energy: -1962.8778, density: 5610.096
    Dyn  step: 37900, energy: -1962.4588, density: 5601.901
    Dyn  step: 38000, energy: -1962.8169, density: 5612.565
    Dyn  step: 38100, energy: -1962.8476, density: 5629.400
    Dyn  step: 38200, energy: -1962.6501, density: 5613.536
    Dyn  step: 38300, energy: -1962.8478, density: 5631.724
    Dyn  step: 38400, energy: -1962.6837, density: 5614.094
    Dyn  step: 38500, energy: -1962.5548, density: 5616.257
    Dyn  step: 38600, energy: -1962.6749, density: 5620.258
    Dyn  step: 38700, energy: -1962.4952, density: 5627.117
    Dyn  step: 38800, energy: -1962.5968, density: 5632.006
    Dyn  step: 38900, energy: -1962.4797, density: 5618.620
    Dyn  step: 39000, energy: -1962.8018, density: 5618.045
    Dyn  step: 39100, energy: -1962.3989, density: 5626.769
    Dyn  step: 39200, energy: -1962.5123, density: 5641.807
    Dyn  step: 39300, energy: -1962.2834, density: 5625.049
    Dyn  step: 39400, energy: -1962.3046, density: 5605.991
    Dyn  step: 39500, energy: -1962.4729, density: 5616.952
    Dyn  step: 39600, energy: -1962.2827, density: 5621.157
    Dyn  step: 39700, energy: -1962.3607, density: 5621.971
    Dyn  step: 39800, energy: -1962.3175, density: 5620.047
    Dyn  step: 39900, energy: -1962.4812, density: 5614.541
    Dyn  step: 40000, energy: -1962.2240, density: 5613.663
    Dyn  step: 40100, energy: -1962.4526, density: 5619.986
    Dyn  step: 40200, energy: -1962.1968, density: 5621.962
    Dyn  step: 40300, energy: -1962.1265, density: 5605.430
    Dyn  step: 40400, energy: -1962.1847, density: 5621.837
    Dyn  step: 40500, energy: -1962.2319, density: 5620.301
    Dyn  step: 40600, energy: -1962.0250, density: 5604.552
    Dyn  step: 40700, energy: -1962.3107, density: 5622.756
    Dyn  step: 40800, energy: -1961.9321, density: 5613.572
    Dyn  step: 40900, energy: -1962.2001, density: 5617.740
    Dyn  step: 41000, energy: -1962.1599, density: 5616.169
    Dyn  step: 41100, energy: -1962.0309, density: 5595.479
    Dyn  step: 41200, energy: -1962.0611, density: 5624.648
    Dyn  step: 41300, energy: -1962.1532, density: 5598.162
    Dyn  step: 41400, energy: -1961.8971, density: 5595.198
    Dyn  step: 41500, energy: -1962.0501, density: 5627.643
    Dyn  step: 41600, energy: -1961.9172, density: 5625.334
    Dyn  step: 41700, energy: -1962.0299, density: 5624.639
    Dyn  step: 41800, energy: -1961.9487, density: 5628.925
    Dyn  step: 41900, energy: -1961.9637, density: 5615.360
    Dyn  step: 42000, energy: -1961.7729, density: 5634.956
    Dyn  step: 42100, energy: -1961.8098, density: 5607.825
    Dyn  step: 42200, energy: -1961.7695, density: 5622.461
    Dyn  step: 42300, energy: -1961.8521, density: 5622.575
    Dyn  step: 42400, energy: -1961.7465, density: 5602.705
    Dyn  step: 42500, energy: -1961.8578, density: 5616.996
    Dyn  step: 42600, energy: -1961.7745, density: 5606.048
    Dyn  step: 42700, energy: -1961.9223, density: 5626.332
    Dyn  step: 42800, energy: -1961.8209, density: 5607.095
    Dyn  step: 42900, energy: -1961.6852, density: 5625.598
    Dyn  step: 43000, energy: -1961.8833, density: 5609.850
    Dyn  step: 43100, energy: -1961.6136, density: 5608.725
    Dyn  step: 43200, energy: -1961.6359, density: 5622.234
    Dyn  step: 43300, energy: -1961.6857, density: 5615.877
    Dyn  step: 43400, energy: -1961.6523, density: 5628.196
    Dyn  step: 43500, energy: -1961.6595, density: 5624.658
    Dyn  step: 43600, energy: -1961.5603, density: 5631.371
    Dyn  step: 43700, energy: -1961.5459, density: 5617.709
    Dyn  step: 43800, energy: -1961.3841, density: 5627.746
    Dyn  step: 43900, energy: -1961.4370, density: 5618.401
    Dyn  step: 44000, energy: -1961.2173, density: 5624.980
    Dyn  step: 44100, energy: -1961.2968, density: 5612.703
    Dyn  step: 44200, energy: -1961.2415, density: 5621.658
    Dyn  step: 44300, energy: -1961.3365, density: 5604.520
    Dyn  step: 44400, energy: -1960.9445, density: 5627.456
    Dyn  step: 44500, energy: -1961.3003, density: 5635.252
    Dyn  step: 44600, energy: -1960.9231, density: 5618.059
    Dyn  step: 44700, energy: -1961.2038, density: 5634.806
    Dyn  step: 44800, energy: -1961.0026, density: 5614.530
    Dyn  step: 44900, energy: -1961.1951, density: 5620.146
    Dyn  step: 45000, energy: -1961.2409, density: 5612.353
    Dyn  step: 45100, energy: -1961.2828, density: 5617.437
    Dyn  step: 45200, energy: -1961.2018, density: 5616.654
    Dyn  step: 45300, energy: -1960.9075, density: 5608.932
    Dyn  step: 45400, energy: -1961.3446, density: 5615.089
    Dyn  step: 45500, energy: -1961.0148, density: 5602.804
    Dyn  step: 45600, energy: -1960.8994, density: 5604.214
    Dyn  step: 45700, energy: -1961.1190, density: 5606.551
    Dyn  step: 45800, energy: -1961.1015, density: 5611.319
    Dyn  step: 45900, energy: -1960.8592, density: 5616.733
    Dyn  step: 46000, energy: -1961.0001, density: 5602.444
    Dyn  step: 46100, energy: -1960.7582, density: 5590.191
    Dyn  step: 46200, energy: -1960.8834, density: 5626.423
    Dyn  step: 46300, energy: -1961.0056, density: 5611.598
    Dyn  step: 46400, energy: -1960.9772, density: 5600.753
    Dyn  step: 46500, energy: -1960.9027, density: 5605.826
    Dyn  step: 46600, energy: -1960.9059, density: 5615.524
    Dyn  step: 46700, energy: -1960.8314, density: 5612.332
    Dyn  step: 46800, energy: -1960.5262, density: 5612.676
    Dyn  step: 46900, energy: -1960.7539, density: 5592.564
    Dyn  step: 47000, energy: -1960.6839, density: 5601.521
    Dyn  step: 47100, energy: -1960.5319, density: 5611.344
    Dyn  step: 47200, energy: -1960.6340, density: 5608.244
    Dyn  step: 47300, energy: -1960.5258, density: 5610.218
    Dyn  step: 47400, energy: -1960.5049, density: 5598.689
    Dyn  step: 47500, energy: -1960.5862, density: 5609.212
    Dyn  step: 47600, energy: -1960.2418, density: 5633.804
    Dyn  step: 47700, energy: -1960.7071, density: 5603.404
    Dyn  step: 47800, energy: -1960.3929, density: 5616.446
    Dyn  step: 47900, energy: -1960.5669, density: 5625.282
    Dyn  step: 48000, energy: -1960.4730, density: 5601.697
    Dyn  step: 48100, energy: -1960.4610, density: 5632.036
    Dyn  step: 48200, energy: -1960.2386, density: 5624.630
    Dyn  step: 48300, energy: -1960.3167, density: 5617.761
    Dyn  step: 48400, energy: -1960.1191, density: 5610.469
    Dyn  step: 48500, energy: -1960.0424, density: 5612.340
    Dyn  step: 48600, energy: -1960.0158, density: 5639.142
    Dyn  step: 48700, energy: -1960.2104, density: 5638.973
    Dyn  step: 48800, energy: -1959.9868, density: 5605.610
    Dyn  step: 48900, energy: -1960.1035, density: 5627.728
    Dyn  step: 49000, energy: -1960.1334, density: 5631.472
    Dyn  step: 49100, energy: -1959.9806, density: 5623.772
    Dyn  step: 49200, energy: -1960.0661, density: 5627.702
    Dyn  step: 49300, energy: -1959.9786, density: 5622.825
    Dyn  step: 49400, energy: -1960.0156, density: 5624.639
    Dyn  step: 49500, energy: -1960.0235, density: 5630.648
    Dyn  step: 49600, energy: -1959.8480, density: 5616.721
    Dyn  step: 49700, energy: -1959.9218, density: 5614.813
    Dyn  step: 49800, energy: -1959.8460, density: 5605.500
    Dyn  step: 49900, energy: -1959.9274, density: 5607.531
    Dyn  step: 50000, energy: -1959.7758, density: 5624.164
    Dyn  step: 50100, energy: -1959.9134, density: 5626.250
    Dyn  step: 50200, energy: -1959.7874, density: 5616.733
    Dyn  step: 50300, energy: -1959.9640, density: 5621.463
    Dyn  step: 50400, energy: -1959.7149, density: 5627.571
    Dyn  step: 50500, energy: -1959.8502, density: 5616.221
    Dyn  step: 50600, energy: -1959.7733, density: 5615.774
    Dyn  step: 50700, energy: -1959.5581, density: 5608.587
    Dyn  step: 50800, energy: -1959.6870, density: 5640.991
    Dyn  step: 50900, energy: -1959.5831, density: 5613.528
    Dyn  step: 51000, energy: -1959.5567, density: 5628.227
    Dyn  step: 51100, energy: -1959.4544, density: 5604.211
    Dyn  step: 51200, energy: -1959.4860, density: 5621.214
    Dyn  step: 51300, energy: -1959.4236, density: 5611.328
    Dyn  step: 51400, energy: -1959.2830, density: 5610.273
    Dyn  step: 51500, energy: -1959.4903, density: 5628.349
    Dyn  step: 51600, energy: -1959.2519, density: 5611.719
    Dyn  step: 51700, energy: -1959.5654, density: 5610.685
    Dyn  step: 51800, energy: -1959.2732, density: 5609.185
    Dyn  step: 51900, energy: -1959.4863, density: 5623.766
    Dyn  step: 52000, energy: -1959.2261, density: 5619.628
    Dyn  step: 52100, energy: -1959.4124, density: 5601.104
    Dyn  step: 52200, energy: -1959.1766, density: 5620.119
    Dyn  step: 52300, energy: -1959.4183, density: 5616.941
    Dyn  step: 52400, energy: -1959.2965, density: 5606.590
    Dyn  step: 52500, energy: -1959.2442, density: 5586.679
    Dyn  step: 52600, energy: -1959.2768, density: 5616.473
    Dyn  step: 52700, energy: -1959.3580, density: 5597.652
    Dyn  step: 52800, energy: -1959.0881, density: 5613.493
    Dyn  step: 52900, energy: -1959.3129, density: 5598.306
    Dyn  step: 53000, energy: -1959.2466, density: 5597.973
    Dyn  step: 53100, energy: -1959.0819, density: 5588.435
    Dyn  step: 53200, energy: -1959.1254, density: 5611.945
    Dyn  step: 53300, energy: -1959.1153, density: 5585.482
    Dyn  step: 53400, energy: -1959.0368, density: 5595.027
    Dyn  step: 53500, energy: -1959.2861, density: 5591.627
    Dyn  step: 53600, energy: -1958.9538, density: 5581.821
    Dyn  step: 53700, energy: -1959.2127, density: 5600.491
    Dyn  step: 53800, energy: -1958.9299, density: 5569.974
    Dyn  step: 53900, energy: -1959.1174, density: 5588.469
    Dyn  step: 54000, energy: -1958.9751, density: 5584.900
    Dyn  step: 54100, energy: -1958.7452, density: 5578.252
    Dyn  step: 54200, energy: -1958.8975, density: 5600.843
    Dyn  step: 54300, energy: -1958.8750, density: 5589.772
    Dyn  step: 54400, energy: -1958.8896, density: 5589.864
    Dyn  step: 54500, energy: -1958.7303, density: 5603.900
    Dyn  step: 54600, energy: -1959.0816, density: 5605.034
    Dyn  step: 54700, energy: -1958.6482, density: 5606.678
    Dyn  step: 54800, energy: -1958.6820, density: 5592.540
    Dyn  step: 54900, energy: -1958.7775, density: 5613.822
    Dyn  step: 55000, energy: -1958.7205, density: 5587.616
    Dyn  step: 55100, energy: -1958.6170, density: 5590.861
    Dyn  step: 55200, energy: -1958.6222, density: 5615.036
    Dyn  step: 55300, energy: -1958.5040, density: 5596.780
    Dyn  step: 55400, energy: -1958.1659, density: 5599.690
    Dyn  step: 55500, energy: -1958.6619, density: 5616.158
    Dyn  step: 55600, energy: -1958.5994, density: 5605.435
    Dyn  step: 55700, energy: -1958.3642, density: 5613.589
    Dyn  step: 55800, energy: -1958.5531, density: 5612.059
    Dyn  step: 55900, energy: -1958.3179, density: 5616.565
    Dyn  step: 56000, energy: -1958.3143, density: 5590.476
    Dyn  step: 56100, energy: -1958.6559, density: 5597.585
    Dyn  step: 56200, energy: -1958.4229, density: 5594.146
    Dyn  step: 56300, energy: -1958.4693, density: 5613.340
    Dyn  step: 56400, energy: -1958.4546, density: 5610.464
    Dyn  step: 56500, energy: -1958.3698, density: 5605.945
    Dyn  step: 56600, energy: -1958.2768, density: 5608.440
    Dyn  step: 56700, energy: -1958.3056, density: 5616.899
    Dyn  step: 56800, energy: -1958.4152, density: 5614.125
    Dyn  step: 56900, energy: -1958.5078, density: 5624.201
    Dyn  step: 57000, energy: -1958.1665, density: 5613.387
    Dyn  step: 57100, energy: -1958.2388, density: 5600.115
    Dyn  step: 57200, energy: -1958.0288, density: 5595.387
    Dyn  step: 57300, energy: -1958.1871, density: 5585.220
    Dyn  step: 57400, energy: -1958.2455, density: 5610.184
    Dyn  step: 57500, energy: -1958.1654, density: 5594.443
    Dyn  step: 57600, energy: -1958.2095, density: 5598.409
    Dyn  step: 57700, energy: -1958.0847, density: 5600.650
    Dyn  step: 57800, energy: -1958.2141, density: 5599.803
    Dyn  step: 57900, energy: -1957.9380, density: 5606.276
    Dyn  step: 58000, energy: -1957.9588, density: 5582.310
    Dyn  step: 58100, energy: -1957.6902, density: 5587.607
    Dyn  step: 58200, energy: -1957.8136, density: 5579.706
    Dyn  step: 58300, energy: -1957.8562, density: 5590.747
    Dyn  step: 58400, energy: -1957.6409, density: 5575.617
    Dyn  step: 58500, energy: -1957.5865, density: 5577.199
    Dyn  step: 58600, energy: -1957.5638, density: 5605.929
    Dyn  step: 58700, energy: -1957.5595, density: 5599.369
    Dyn  step: 58800, energy: -1957.4937, density: 5599.722
    Dyn  step: 58900, energy: -1957.4938, density: 5593.199
    Dyn  step: 59000, energy: -1957.4266, density: 5599.506
    Dyn  step: 59100, energy: -1957.6496, density: 5619.589
    Dyn  step: 59200, energy: -1957.6061, density: 5623.246
    Dyn  step: 59300, energy: -1957.5175, density: 5594.936
    Dyn  step: 59400, energy: -1957.3450, density: 5593.588
    Dyn  step: 59500, energy: -1957.3843, density: 5610.295
    Dyn  step: 59600, energy: -1957.4628, density: 5610.979
    Dyn  step: 59700, energy: -1957.2795, density: 5622.804
    Dyn  step: 59800, energy: -1957.1339, density: 5603.189
    Dyn  step: 59900, energy: -1957.1516, density: 5612.850
    Dyn  step: 60000, energy: -1957.2908, density: 5612.673
    Dyn  step: 60100, energy: -1957.3437, density: 5624.193
    Dyn  step: 60200, energy: -1957.2244, density: 5596.199
    Dyn  step: 60300, energy: -1957.3391, density: 5624.245
    Dyn  step: 60400, energy: -1957.1703, density: 5596.473
    Dyn  step: 60500, energy: -1957.1161, density: 5614.795
    Dyn  step: 60600, energy: -1956.9722, density: 5619.824
    Dyn  step: 60700, energy: -1957.2022, density: 5632.494
    Dyn  step: 60800, energy: -1956.9791, density: 5620.551
    Dyn  step: 60900, energy: -1957.0932, density: 5616.070
    Dyn  step: 61000, energy: -1956.9957, density: 5610.799
    Dyn  step: 61100, energy: -1957.0870, density: 5570.964
    Dyn  step: 61200, energy: -1956.7155, density: 5592.847
    Dyn  step: 61300, energy: -1957.0754, density: 5594.168
    Dyn  step: 61400, energy: -1956.7745, density: 5573.856
    Dyn  step: 61500, energy: -1957.0008, density: 5588.352
    Dyn  step: 61600, energy: -1957.0052, density: 5580.402
    Dyn  step: 61700, energy: -1957.1388, density: 5613.314
    Dyn  step: 61800, energy: -1956.9723, density: 5590.553
    Dyn  step: 61900, energy: -1957.2155, density: 5591.911
    Dyn  step: 62000, energy: -1956.8547, density: 5586.402
    Dyn  step: 62100, energy: -1957.0732, density: 5591.093
    Dyn  step: 62200, energy: -1956.9196, density: 5591.135
    Dyn  step: 62300, energy: -1956.9027, density: 5594.512
    Dyn  step: 62400, energy: -1956.8162, density: 5592.965
    Dyn  step: 62500, energy: -1956.6743, density: 5589.171
    Dyn  step: 62600, energy: -1956.9279, density: 5606.231
    Dyn  step: 62700, energy: -1956.3587, density: 5599.144
    Dyn  step: 62800, energy: -1956.9371, density: 5585.843
    Dyn  step: 62900, energy: -1956.5812, density: 5587.686
    Dyn  step: 63000, energy: -1956.7251, density: 5593.334
    Dyn  step: 63100, energy: -1956.6423, density: 5586.167
    Dyn  step: 63200, energy: -1956.7331, density: 5604.047
    Dyn  step: 63300, energy: -1956.3336, density: 5588.126
    Dyn  step: 63400, energy: -1956.5518, density: 5602.686
    Dyn  step: 63500, energy: -1956.4835, density: 5594.444
    Dyn  step: 63600, energy: -1956.1955, density: 5578.602
    Dyn  step: 63700, energy: -1956.5090, density: 5602.578
    Dyn  step: 63800, energy: -1956.2963, density: 5583.162
    Dyn  step: 63900, energy: -1956.3614, density: 5598.318
    Dyn  step: 64000, energy: -1956.1024, density: 5581.613
    Dyn  step: 64100, energy: -1956.5036, density: 5604.777
    Dyn  step: 64200, energy: -1956.1400, density: 5587.951
    Dyn  step: 64300, energy: -1956.1806, density: 5594.336
    Dyn  step: 64400, energy: -1956.2761, density: 5602.834
    Dyn  step: 64500, energy: -1955.9885, density: 5573.619
    Dyn  step: 64600, energy: -1956.0839, density: 5592.611
    Dyn  step: 64700, energy: -1955.9633, density: 5602.966
    Dyn  step: 64800, energy: -1956.3961, density: 5607.744
    Dyn  step: 64900, energy: -1955.9438, density: 5578.558
    Dyn  step: 65000, energy: -1956.3430, density: 5616.836
    Dyn  step: 65100, energy: -1956.2180, density: 5579.600
    Dyn  step: 65200, energy: -1956.2270, density: 5613.568
    Dyn  step: 65300, energy: -1955.9812, density: 5583.704
    Dyn  step: 65400, energy: -1955.9131, density: 5596.596
    Dyn  step: 65500, energy: -1955.7710, density: 5598.157
    Dyn  step: 65600, energy: -1955.8270, density: 5604.711
    Dyn  step: 65700, energy: -1955.7949, density: 5580.568
    Dyn  step: 65800, energy: -1956.0379, density: 5585.562
    Dyn  step: 65900, energy: -1955.8226, density: 5587.274
    Dyn  step: 66000, energy: -1955.7630, density: 5606.488
    Dyn  step: 66100, energy: -1955.6091, density: 5603.546
    Dyn  step: 66200, energy: -1955.7639, density: 5610.982
    Dyn  step: 66300, energy: -1955.7765, density: 5600.487
    Dyn  step: 66400, energy: -1955.7290, density: 5590.557
    Dyn  step: 66500, energy: -1955.6920, density: 5594.132
    Dyn  step: 66600, energy: -1955.7709, density: 5604.894
    Dyn  step: 66700, energy: -1955.4329, density: 5594.894
    Dyn  step: 66800, energy: -1955.5999, density: 5593.211
    Dyn  step: 66900, energy: -1955.4603, density: 5573.618
    Dyn  step: 67000, energy: -1955.6767, density: 5599.976
    Dyn  step: 67100, energy: -1955.4671, density: 5590.567
    Dyn  step: 67200, energy: -1955.8938, density: 5613.440
    Dyn  step: 67300, energy: -1955.6825, density: 5603.377
    Dyn  step: 67400, energy: -1956.0855, density: 5606.721
    Dyn  step: 67500, energy: -1955.6345, density: 5583.056
    Dyn  step: 67600, energy: -1955.6951, density: 5600.475
    Dyn  step: 67700, energy: -1955.5326, density: 5599.371
    Dyn  step: 67800, energy: -1955.6227, density: 5613.622
    Dyn  step: 67900, energy: -1955.3926, density: 5609.713
    Dyn  step: 68000, energy: -1955.6219, density: 5602.111
    Dyn  step: 68100, energy: -1955.2658, density: 5596.160
    Dyn  step: 68200, energy: -1955.3679, density: 5616.925
    Dyn  step: 68300, energy: -1954.8461, density: 5595.521
    Dyn  step: 68400, energy: -1955.3476, density: 5594.669
    Dyn  step: 68500, energy: -1955.3541, density: 5591.550
    Dyn  step: 68600, energy: -1955.2697, density: 5611.408
    Dyn  step: 68700, energy: -1955.0252, density: 5618.617
    Dyn  step: 68800, energy: -1955.1701, density: 5590.923
    Dyn  step: 68900, energy: -1955.2554, density: 5595.919
    Dyn  step: 69000, energy: -1955.1019, density: 5586.668
    Dyn  step: 69100, energy: -1955.0881, density: 5594.001
    Dyn  step: 69200, energy: -1954.8298, density: 5598.089
    Dyn  step: 69300, energy: -1954.9630, density: 5595.149
    Dyn  step: 69400, energy: -1954.8875, density: 5603.134
    Dyn  step: 69500, energy: -1954.7627, density: 5589.086
    Dyn  step: 69600, energy: -1954.9935, density: 5619.555
    Dyn  step: 69700, energy: -1954.8294, density: 5587.110
    Dyn  step: 69800, energy: -1954.8781, density: 5605.314
    Dyn  step: 69900, energy: -1955.0338, density: 5584.290
    Dyn  step: 70000, energy: -1954.7487, density: 5594.347
    Dyn  step: 70100, energy: -1954.6705, density: 5605.534
    Dyn  step: 70200, energy: -1954.6326, density: 5598.383
    Dyn  step: 70300, energy: -1954.5321, density: 5599.936
    Dyn  step: 70400, energy: -1954.7661, density: 5600.961
    Dyn  step: 70500, energy: -1954.6282, density: 5593.711
    Dyn  step: 70600, energy: -1954.6012, density: 5591.419
    Dyn  step: 70700, energy: -1954.4677, density: 5604.208
    Dyn  step: 70800, energy: -1954.6148, density: 5590.536
    Dyn  step: 70900, energy: -1954.9242, density: 5614.281
    Dyn  step: 71000, energy: -1954.5462, density: 5589.961
    Dyn  step: 71100, energy: -1954.6938, density: 5592.335
    Dyn  step: 71200, energy: -1954.2387, density: 5582.336
    Dyn  step: 71300, energy: -1954.3439, density: 5584.045
    Dyn  step: 71400, energy: -1954.3086, density: 5596.492
    Dyn  step: 71500, energy: -1954.3477, density: 5589.956
    Dyn  step: 71600, energy: -1954.7459, density: 5609.827
    Dyn  step: 71700, energy: -1954.4702, density: 5592.848
    Dyn  step: 71800, energy: -1954.5654, density: 5603.643
    Dyn  step: 71900, energy: -1954.4505, density: 5612.575
    Dyn  step: 72000, energy: -1954.3565, density: 5578.530
    Dyn  step: 72100, energy: -1954.2630, density: 5592.609
    Dyn  step: 72200, energy: -1954.0607, density: 5585.596
    Dyn  step: 72300, energy: -1954.2838, density: 5582.028
    Dyn  step: 72400, energy: -1953.9690, density: 5589.455
    Dyn  step: 72500, energy: -1954.2908, density: 5595.292
    Dyn  step: 72600, energy: -1953.8741, density: 5600.859
    Dyn  step: 72700, energy: -1954.0279, density: 5586.172
    Dyn  step: 72800, energy: -1953.9993, density: 5596.431
    Dyn  step: 72900, energy: -1954.1743, density: 5592.900
    Dyn  step: 73000, energy: -1954.0731, density: 5606.472
    Dyn  step: 73100, energy: -1953.9863, density: 5595.256
    Dyn  step: 73200, energy: -1954.0125, density: 5587.987
    Dyn  step: 73300, energy: -1953.5682, density: 5580.944
    Dyn  step: 73400, energy: -1953.9688, density: 5601.443
    Dyn  step: 73500, energy: -1954.0286, density: 5581.659
    Dyn  step: 73600, energy: -1953.9900, density: 5595.487
    Dyn  step: 73700, energy: -1953.8984, density: 5594.979
    Dyn  step: 73800, energy: -1953.9078, density: 5619.149
    Dyn  step: 73900, energy: -1953.8589, density: 5595.869
    Dyn  step: 74000, energy: -1953.8631, density: 5593.039
    Dyn  step: 74100, energy: -1953.9491, density: 5582.328
    Dyn  step: 74200, energy: -1954.1052, density: 5588.525
    Dyn  step: 74300, energy: -1953.6441, density: 5585.554
    Dyn  step: 74400, energy: -1953.8719, density: 5567.370
    Dyn  step: 74500, energy: -1953.6442, density: 5586.246
    Dyn  step: 74600, energy: -1954.0215, density: 5583.021
    Dyn  step: 74700, energy: -1953.7519, density: 5597.352
    Dyn  step: 74800, energy: -1953.8318, density: 5586.060
    Dyn  step: 74900, energy: -1953.6713, density: 5578.060
    Dyn  step: 75000, energy: -1953.6903, density: 5568.666
    Dyn  step: 75100, energy: -1953.6453, density: 5580.557
    Dyn  step: 75200, energy: -1953.7127, density: 5588.739
    Dyn  step: 75300, energy: -1953.5850, density: 5584.751
    Dyn  step: 75400, energy: -1953.4635, density: 5601.790
    Dyn  step: 75500, energy: -1953.4277, density: 5602.077
    Dyn  step: 75600, energy: -1953.3458, density: 5586.104
    Dyn  step: 75700, energy: -1953.5813, density: 5582.647
    Dyn  step: 75800, energy: -1953.1676, density: 5601.640
    Dyn  step: 75900, energy: -1953.3312, density: 5602.981
    Dyn  step: 76000, energy: -1953.2855, density: 5586.939
    Dyn  step: 76100, energy: -1953.1821, density: 5596.201
    Dyn  step: 76200, energy: -1952.8274, density: 5591.383
    Dyn  step: 76300, energy: -1953.0680, density: 5598.015
    Dyn  step: 76400, energy: -1953.0125, density: 5583.124
    Dyn  step: 76500, energy: -1952.9821, density: 5601.292
    Dyn  step: 76600, energy: -1953.0372, density: 5597.517
    Dyn  step: 76700, energy: -1952.7702, density: 5562.932
    Dyn  step: 76800, energy: -1953.1422, density: 5602.262
    Dyn  step: 76900, energy: -1952.9120, density: 5584.877
    Dyn  step: 77000, energy: -1953.1607, density: 5587.109
    Dyn  step: 77100, energy: -1952.9234, density: 5592.730
    Dyn  step: 77200, energy: -1953.2222, density: 5586.784
    Dyn  step: 77300, energy: -1952.9846, density: 5587.141
    Dyn  step: 77400, energy: -1953.2229, density: 5603.363
    Dyn  step: 77500, energy: -1952.9532, density: 5578.864
    Dyn  step: 77600, energy: -1952.9951, density: 5592.557
    Dyn  step: 77700, energy: -1952.4864, density: 5589.117
    Dyn  step: 77800, energy: -1953.2513, density: 5595.929
    Dyn  step: 77900, energy: -1952.7902, density: 5595.143
    Dyn  step: 78000, energy: -1952.9775, density: 5594.111
    Dyn  step: 78100, energy: -1952.7461, density: 5604.847
    Dyn  step: 78200, energy: -1952.8700, density: 5579.442
    Dyn  step: 78300, energy: -1952.8699, density: 5573.066
    Dyn  step: 78400, energy: -1952.6576, density: 5590.442
    Dyn  step: 78500, energy: -1952.6228, density: 5573.559
    Dyn  step: 78600, energy: -1952.7548, density: 5578.279
    Dyn  step: 78700, energy: -1952.8009, density: 5580.696
    Dyn  step: 78800, energy: -1952.7020, density: 5549.131
    Dyn  step: 78900, energy: -1952.5624, density: 5550.057
    Dyn  step: 79000, energy: -1952.7014, density: 5577.276
    Dyn  step: 79100, energy: -1952.8265, density: 5594.223
    Dyn  step: 79200, energy: -1952.5589, density: 5588.920
    Dyn  step: 79300, energy: -1952.4337, density: 5568.959
    Dyn  step: 79400, energy: -1952.6061, density: 5569.172
    Dyn  step: 79500, energy: -1952.5686, density: 5582.177
    Dyn  step: 79600, energy: -1952.4002, density: 5566.256
    Dyn  step: 79700, energy: -1952.6721, density: 5573.865
    Dyn  step: 79800, energy: -1952.4526, density: 5572.733
    Dyn  step: 79900, energy: -1952.5646, density: 5570.393
    Dyn  step: 80000, energy: -1952.3857, density: 5580.442
    Dyn  step: 80100, energy: -1952.1662, density: 5579.828
    Dyn  step: 80200, energy: -1952.3304, density: 5568.008
    Dyn  step: 80300, energy: -1952.1455, density: 5589.540
    Dyn  step: 80400, energy: -1952.2993, density: 5569.795
    Dyn  step: 80500, energy: -1952.3954, density: 5580.692
    Dyn  step: 80600, energy: -1952.3602, density: 5572.837
    Dyn  step: 80700, energy: -1952.4242, density: 5569.356
    Dyn  step: 80800, energy: -1952.1961, density: 5575.605
    Dyn  step: 80900, energy: -1952.4573, density: 5576.938
    Dyn  step: 81000, energy: -1952.3411, density: 5576.217
    Dyn  step: 81100, energy: -1952.3828, density: 5595.221
    Dyn  step: 81200, energy: -1952.1836, density: 5572.757
    Dyn  step: 81300, energy: -1952.3286, density: 5578.247
    Dyn  step: 81400, energy: -1952.2768, density: 5594.166
    Dyn  step: 81500, energy: -1952.3009, density: 5559.543
    Dyn  step: 81600, energy: -1951.8514, density: 5564.776
    Dyn  step: 81700, energy: -1952.2882, density: 5577.790
    Dyn  step: 81800, energy: -1951.8948, density: 5601.579
    Dyn  step: 81900, energy: -1951.9753, density: 5562.647
    Dyn  step: 82000, energy: -1951.6850, density: 5557.558
    Dyn  step: 82100, energy: -1951.7971, density: 5566.299
    Dyn  step: 82200, energy: -1951.7544, density: 5561.835
    Dyn  step: 82300, energy: -1952.0016, density: 5583.803
    Dyn  step: 82400, energy: -1951.6483, density: 5565.191
    Dyn  step: 82500, energy: -1951.8680, density: 5588.729
    Dyn  step: 82600, energy: -1951.8206, density: 5592.352
    Dyn  step: 82700, energy: -1951.7525, density: 5580.995
    Dyn  step: 82800, energy: -1951.8413, density: 5576.097
    Dyn  step: 82900, energy: -1951.5087, density: 5569.149
    Dyn  step: 83000, energy: -1951.6793, density: 5599.911
    Dyn  step: 83100, energy: -1951.5518, density: 5580.961
    Dyn  step: 83200, energy: -1951.6000, density: 5578.459
    Dyn  step: 83300, energy: -1951.3654, density: 5587.611
    Dyn  step: 83400, energy: -1951.5738, density: 5601.127
    Dyn  step: 83500, energy: -1951.2251, density: 5589.532
    Dyn  step: 83600, energy: -1951.5956, density: 5593.838
    Dyn  step: 83700, energy: -1951.2360, density: 5597.423
    Dyn  step: 83800, energy: -1951.3719, density: 5614.673
    Dyn  step: 83900, energy: -1951.2838, density: 5599.519
    Dyn  step: 84000, energy: -1951.2803, density: 5611.340
    Dyn  step: 84100, energy: -1951.2114, density: 5599.671
    Dyn  step: 84200, energy: -1951.1345, density: 5601.091
    Dyn  step: 84300, energy: -1951.0118, density: 5581.078
    Dyn  step: 84400, energy: -1951.0261, density: 5589.732
    Dyn  step: 84500, energy: -1950.9941, density: 5609.054
    Dyn  step: 84600, energy: -1950.8105, density: 5581.923
    Dyn  step: 84700, energy: -1951.0252, density: 5583.471
    Dyn  step: 84800, energy: -1951.2199, density: 5589.532
    Dyn  step: 84900, energy: -1951.1159, density: 5578.512
    Dyn  step: 85000, energy: -1951.4534, density: 5595.867
    Dyn  step: 85100, energy: -1950.9510, density: 5573.930
    Dyn  step: 85200, energy: -1951.0460, density: 5578.910
    Dyn  step: 85300, energy: -1950.8255, density: 5604.792
    Dyn  step: 85400, energy: -1950.9537, density: 5590.190
    Dyn  step: 85500, energy: -1950.7137, density: 5591.773
    Dyn  step: 85600, energy: -1950.8465, density: 5591.753
    Dyn  step: 85700, energy: -1950.9985, density: 5588.457
    Dyn  step: 85800, energy: -1950.5949, density: 5572.520
    Dyn  step: 85900, energy: -1950.7914, density: 5596.814
    Dyn  step: 86000, energy: -1950.7732, density: 5580.093
    Dyn  step: 86100, energy: -1951.0157, density: 5610.594
    Dyn  step: 86200, energy: -1950.6364, density: 5579.659
    Dyn  step: 86300, energy: -1950.9174, density: 5580.947
    Dyn  step: 86400, energy: -1950.8877, density: 5602.576
    Dyn  step: 86500, energy: -1950.6814, density: 5562.561
    Dyn  step: 86600, energy: -1950.9691, density: 5592.331
    Dyn  step: 86700, energy: -1950.5628, density: 5579.484
    Dyn  step: 86800, energy: -1950.8400, density: 5579.206
    Dyn  step: 86900, energy: -1950.7325, density: 5589.579
    Dyn  step: 87000, energy: -1950.9007, density: 5589.906
    Dyn  step: 87100, energy: -1950.8340, density: 5608.913
    Dyn  step: 87200, energy: -1950.6064, density: 5601.648
    Dyn  step: 87300, energy: -1950.5506, density: 5577.151
    Dyn  step: 87400, energy: -1950.5075, density: 5594.023
    Dyn  step: 87500, energy: -1950.6619, density: 5591.606
    Dyn  step: 87600, energy: -1950.3194, density: 5579.725
    Dyn  step: 87700, energy: -1950.4465, density: 5567.418
    Dyn  step: 87800, energy: -1950.2478, density: 5589.112
    Dyn  step: 87900, energy: -1950.2994, density: 5561.739
    Dyn  step: 88000, energy: -1950.2121, density: 5574.156
    Dyn  step: 88100, energy: -1950.2224, density: 5567.255
    Dyn  step: 88200, energy: -1949.9564, density: 5569.532
    Dyn  step: 88300, energy: -1949.8935, density: 5581.269
    Dyn  step: 88400, energy: -1949.9627, density: 5556.505
    Dyn  step: 88500, energy: -1950.0514, density: 5576.195
    Dyn  step: 88600, energy: -1950.0580, density: 5586.443
    Dyn  step: 88700, energy: -1949.9882, density: 5564.291
    Dyn  step: 88800, energy: -1950.0786, density: 5582.431
    Dyn  step: 88900, energy: -1950.0181, density: 5588.979
    Dyn  step: 89000, energy: -1950.1278, density: 5566.329
    Dyn  step: 89100, energy: -1949.8901, density: 5590.345
    Dyn  step: 89200, energy: -1950.1802, density: 5571.361
    Dyn  step: 89300, energy: -1949.7517, density: 5564.441
    Dyn  step: 89400, energy: -1950.1843, density: 5583.623
    Dyn  step: 89500, energy: -1950.0662, density: 5556.395
    Dyn  step: 89600, energy: -1949.9111, density: 5569.918
    Dyn  step: 89700, energy: -1950.2224, density: 5553.612
    Dyn  step: 89800, energy: -1949.7912, density: 5565.264
    Dyn  step: 89900, energy: -1950.1659, density: 5583.774
    Dyn  step: 90000, energy: -1949.8724, density: 5562.923
    Dyn  step: 90100, energy: -1949.8896, density: 5574.388
    Dyn  step: 90200, energy: -1949.8226, density: 5589.772
    Dyn  step: 90300, energy: -1949.8860, density: 5586.709
    Dyn  step: 90400, energy: -1949.6986, density: 5574.818
    Dyn  step: 90500, energy: -1949.6811, density: 5577.220
    Dyn  step: 90600, energy: -1949.5247, density: 5578.194
    Dyn  step: 90700, energy: -1949.2766, density: 5574.882
    Dyn  step: 90800, energy: -1949.4940, density: 5584.033
    Dyn  step: 90900, energy: -1949.6857, density: 5592.757
    Dyn  step: 91000, energy: -1949.2359, density: 5587.076
    Dyn  step: 91100, energy: -1949.6526, density: 5571.746
    Dyn  step: 91200, energy: -1949.4079, density: 5551.467
    Dyn  step: 91300, energy: -1949.5630, density: 5572.494
    Dyn  step: 91400, energy: -1949.1155, density: 5574.437
    Dyn  step: 91500, energy: -1949.4786, density: 5584.170
    Dyn  step: 91600, energy: -1949.3255, density: 5568.608
    Dyn  step: 91700, energy: -1949.3882, density: 5580.475
    Dyn  step: 91800, energy: -1949.1936, density: 5576.033
    Dyn  step: 91900, energy: -1949.1315, density: 5568.321
    Dyn  step: 92000, energy: -1949.0484, density: 5567.760
    Dyn  step: 92100, energy: -1949.4007, density: 5624.814
    Dyn  step: 92200, energy: -1949.3201, density: 5581.510
    Dyn  step: 92300, energy: -1949.2247, density: 5572.795
    Dyn  step: 92400, energy: -1949.2827, density: 5589.198
    Dyn  step: 92500, energy: -1949.0837, density: 5590.716
    Dyn  step: 92600, energy: -1949.5512, density: 5563.808
    Dyn  step: 92700, energy: -1948.9534, density: 5567.310
    Dyn  step: 92800, energy: -1949.3249, density: 5564.754
    Dyn  step: 92900, energy: -1948.9468, density: 5553.508
    Dyn  step: 93000, energy: -1949.2227, density: 5575.898
    Dyn  step: 93100, energy: -1949.0166, density: 5557.179
    Dyn  step: 93200, energy: -1948.8215, density: 5554.171
    Dyn  step: 93300, energy: -1949.0624, density: 5575.852
    Dyn  step: 93400, energy: -1948.8433, density: 5583.947
    Dyn  step: 93500, energy: -1948.9845, density: 5570.556
    Dyn  step: 93600, energy: -1948.8730, density: 5563.210
    Dyn  step: 93700, energy: -1948.9402, density: 5579.697
    Dyn  step: 93800, energy: -1948.6810, density: 5579.968
    Dyn  step: 93900, energy: -1948.6728, density: 5586.081
    Dyn  step: 94000, energy: -1948.8570, density: 5574.925
    Dyn  step: 94100, energy: -1948.9034, density: 5562.266
    Dyn  step: 94200, energy: -1948.8480, density: 5552.657
    Dyn  step: 94300, energy: -1948.8249, density: 5561.882
    Dyn  step: 94400, energy: -1948.4016, density: 5552.247
    Dyn  step: 94500, energy: -1948.7273, density: 5554.850
    Dyn  step: 94600, energy: -1948.4014, density: 5565.994
    Dyn  step: 94700, energy: -1948.5118, density: 5561.064
    Dyn  step: 94800, energy: -1948.7610, density: 5575.463
    Dyn  step: 94900, energy: -1948.3757, density: 5564.079
    Dyn  step: 95000, energy: -1948.7697, density: 5573.713
    Dyn  step: 95100, energy: -1948.7621, density: 5573.350
    Dyn  step: 95200, energy: -1948.5219, density: 5581.187
    Dyn  step: 95300, energy: -1948.3478, density: 5559.405
    Dyn  step: 95400, energy: -1948.3799, density: 5578.499
    Dyn  step: 95500, energy: -1948.3157, density: 5551.566
    Dyn  step: 95600, energy: -1948.3044, density: 5568.872
    Dyn  step: 95700, energy: -1948.2829, density: 5559.549
    Dyn  step: 95800, energy: -1948.4026, density: 5580.400
    Dyn  step: 95900, energy: -1948.5914, density: 5595.882
    Dyn  step: 96000, energy: -1948.3561, density: 5560.729
    Dyn  step: 96100, energy: -1948.3983, density: 5590.427
    Dyn  step: 96200, energy: -1948.1941, density: 5577.259
    Dyn  step: 96300, energy: -1948.1726, density: 5555.641
    Dyn  step: 96400, energy: -1948.1405, density: 5576.667
    Dyn  step: 96500, energy: -1948.2485, density: 5567.652
    Dyn  step: 96600, energy: -1948.1411, density: 5561.599
    Dyn  step: 96700, energy: -1948.0311, density: 5574.032
    Dyn  step: 96800, energy: -1947.9159, density: 5560.517
    Dyn  step: 96900, energy: -1947.8023, density: 5573.325
    Dyn  step: 97000, energy: -1948.0153, density: 5558.244
    Dyn  step: 97100, energy: -1947.9308, density: 5560.638
    Dyn  step: 97200, energy: -1947.8382, density: 5572.483
    Dyn  step: 97300, energy: -1948.0369, density: 5550.152
    Dyn  step: 97400, energy: -1947.9298, density: 5558.990
    Dyn  step: 97500, energy: -1947.5024, density: 5546.247
    Dyn  step: 97600, energy: -1948.0358, density: 5574.032
    Dyn  step: 97700, energy: -1947.7092, density: 5538.887
    Dyn  step: 97800, energy: -1947.8526, density: 5543.995
    Dyn  step: 97900, energy: -1947.5540, density: 5561.197
    Dyn  step: 98000, energy: -1947.9925, density: 5551.678
    Dyn  step: 98100, energy: -1947.4339, density: 5529.232
    Dyn  step: 98200, energy: -1947.4383, density: 5536.276
    Dyn  step: 98300, energy: -1947.6271, density: 5540.688
    Dyn  step: 98400, energy: -1947.4332, density: 5546.672
    Dyn  step: 98500, energy: -1947.6573, density: 5515.290
    Dyn  step: 98600, energy: -1947.7385, density: 5556.119
    Dyn  step: 98700, energy: -1947.5045, density: 5536.630
    Dyn  step: 98800, energy: -1947.4416, density: 5548.676
    Dyn  step: 98900, energy: -1947.5769, density: 5568.480
    Dyn  step: 99000, energy: -1947.5297, density: 5553.034
    Dyn  step: 99100, energy: -1947.2910, density: 5564.180
    Dyn  step: 99200, energy: -1947.6288, density: 5554.287
    Dyn  step: 99300, energy: -1947.3849, density: 5574.066
    Dyn  step: 99400, energy: -1947.6595, density: 5581.671
    Dyn  step: 99500, energy: -1947.3564, density: 5558.499
    Dyn  step: 99600, energy: -1947.3424, density: 5562.978
    Dyn  step: 99700, energy: -1947.4017, density: 5579.387
    Dyn  step: 99800, energy: -1947.4773, density: 5570.695
    Dyn  step: 99900, energy: -1947.4484, density: 5575.208
    Dyn  step: 100000, energy: -1947.5257, density: 5573.151





    True



### 結果のプロット

※バックグラウンド実行後に集計する際はコメントアウトを外してライブラリをインポートして下さい


```python
# import matplotlib.pyplot as plt
# from ase.io import Trajectory
# tag = "./output/BaTiO3_tetragonal_NPT"

traj = Trajectory(tag+".traj")
x = [ss.cell.cellpar()[0] / 4 for ss in traj] 
y = [ss.cell.cellpar()[1] / 4 for ss in traj] 
z = [ss.cell.cellpar()[2] / 4 for ss in traj] 
t = [1.0 * i for i in range(len(x))]
temp = [ss.get_temperature() for ss in traj]
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax2 = ax1.twinx()

ax1.plot(t,temp,label="T", color="red")

ax2.plot(t,x,label="Cell X")
ax2.plot(t,y,label="Cell Y")
ax2.plot(t,z,label="Cell Z")


ax1.set_xlabel("Time / fs")
ax2.set_ylabel("Cell / Å")
ax1.set_ylabel("T / K")
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()

ax1.legend(h1+h2, l1+l2, loc='lower right')
fig.show()
```


    
![png](output_20_0.png)
    



```python
# from ase.geometry.analysis import Analysis

atoms1 = traj[-10]
atoms2 = traj[1000]    # MD計算が安定した頃のタイムステップから構造を抽出しています
tags = ["Cubic", "Tetragonal"]
rmax = 5
nbins = 100
delta = [(rmax/nbins)*i for i in range(nbins)]
for bto, name in zip([atoms1, atoms2], tags):
    atoms = bto.copy()
    print(name, ":",atoms.cell.cellpar())
    #atoms *= (4,4,4)
    ana = Analysis(atoms)
    rdf = ana.get_rdf(rmax,nbins,elements=["Ti","O"])[0]
    plt.plot(delta, rdf, label=name)
plt.xlabel("d / Å")
plt.ylabel("RDF")
plt.legend()
plt.show()
```

    Cubic : [16.49943864 16.38873524 16.42598681 90.00004028 89.99995033 90.00000356]
    Tetragonal : [16.22180295 16.177742   16.74072992 90.00004028 89.99995033 90.00000356]



    
![png](output_21_1.png)
    


Tetragonalを初期構造としてMD計算を行うことで
- 400K 周辺で格子定数が同じ値に近づく（＝cubicへの相転移が確認できる）
- tetragonal構造はRDFの2Å付近で3か所のピークに分かれる
- cubic構造もMD計算でRDFピークは広がりを見せる

ことを確認して本事例の説明を終えます。
