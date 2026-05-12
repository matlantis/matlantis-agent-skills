Copyright Preferred Computational Chemistry, Inc. as contributors to Matlantis contrib project

# MD: LGPS中のLi拡散

本計算事例では、Matlantisを用いたLi拡散のMD(分子動力学法)による解析を扱います。&lt;br/&gt;
詳細は以下で記載されています。

 - [硫化物固体電解質中のLi拡散](https://matlantis.com/ja/calculation/li-diffusion-in-li10gep2s12-sulfide-solid-electrolyte)

## Setup

必要なライブラリのインストール、モジュールのimport、関数定義を行います。


```python
# Please install these libraries only for first time of execution
!pip install pfp_api_client
!pip install pandas tqdm matplotlib seaborn optuna sklearn ase
```




```python
import pathlib
EXAMPLE_DIR = pathlib.Path("__file__").resolve().parent
INPUT_DIR = EXAMPLE_DIR / "input"
OUTPUT_DIR = EXAMPLE_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)
```


```python
# 汎用モジュール
import numpy as np
import pandas as pd
from tqdm import tqdm_notebook as tqdm
from IPython.display import Image, display_png
import ipywidgets as widgets
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
from matplotlib.animation import PillowWriter
import seaborn as sns
import math
import optuna
import nglview as nv
import os,sys,csv,glob,shutil,re,time
from time import perf_counter
from joblib import Parallel, delayed

# sklearn
from sklearn.metrics import mean_absolute_error

import ase
from ase.visualize import view
from ase.optimize import BFGS
from ase.constraints import FixAtoms, FixedPlane, FixBondLength, ExpCellFilter

from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase.md import MDLogger
from ase import Atoms
from ase.io import read, write
from ase.io import Trajectory
from ase import units

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

estimator = Estimator(model_version="v2.0.0", calc_mode=EstimatorCalcMode.PBE_U)
calculator = ASECalculator(estimator)
```


```python
def myopt(m,sn = 10,constraintatoms=[],cbonds=[]):
    fa = FixAtoms(indices=constraintatoms)
    fb = FixBondLengths(cbonds,tolerance=1e-5,)
    m.set_constraint([fa,fb])
    m.set_calculator(calculator)
    maxf = np.sqrt(((m.get_forces())**2).sum(axis=1).max())
    print("ini   pot:{:.4f},maxforce:{:.4f}".format(m.get_potential_energy(),maxf))
    de = -1
    s = 1
    ita = 50
    while ( de  &lt; -0.001 or de &gt; 0.001 ) and s &lt;= sn :
        opt = BFGS(m,maxstep=0.04*(0.9**s),logfile=None)
        old  =  m.get_potential_energy()
        opt.run(fmax=0.0005,steps =ita)
        maxf = np.sqrt(((m.get_forces())**2).sum(axis=1).max())
        de =  m.get_potential_energy()  - old
        print("{} pot:{:.4f},maxforce:{:.4f},delta:{:.4f}".format(s*ita,m.get_potential_energy(),maxf,de))
        s += 1
    return m

def opt_cell_size(m,sn = 10, iter_count = False): # m:Atomsオブジェクト
    m.set_constraint() # clear constraint
    m.set_calculator(calculator)
    maxf = np.sqrt(((m.get_forces())**2).sum(axis=1).max()) # √(fx^2 + fy^2 + fz^2)の一番大きいものを取得
    ucf = ExpCellFilter(m)
    print("ini   pot:{:.4f},maxforce:{:.4f}".format(m.get_potential_energy(),maxf))
    de = -1
    s = 1
    ita = 50
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

## LGPS構造の準備

LGPS構造はMaterials Projectから用意しました。

ASE io moduleを用いて、構造ファイルを読み込んだ後、構造緩和を行い安定なセルサイズなどを得ます。

 - https://wiki.fysik.dtu.dk/ase/ase/io/io.html
 - https://wiki.fysik.dtu.dk/ase/ase/optimize.html

Input cif file is from
A. Jain*, S.P. Ong*, G. Hautier, W. Chen, W.D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, K.A. Persson (*=equal contributions)
The Materials Project: A materials genome approach to accelerating materials innovation
APL Materials, 2013, 1(1), 011002.
[doi:10.1063/1.4812323](http://dx.doi.org/10.1063/1.4812323)
[[bibtex]](https://materialsproject.org/static/docs/jain_ong2013.349ca3156250.bib)
Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)


```python
bulk = read(INPUT_DIR / "Li10Ge(PS6)2_mp-696138_computed.cif")
bulk.calc = calculator

print("# atoms =", len(bulk))
print("Initial lattice constant =", bulk.cell.cellpar())

opt_cell_size(bulk)
print ("Optimized lattice constant =", bulk.cell.cellpar())
```

    # atoms = 50
    Initial lattice constant = [ 8.589951    8.87954092 12.97496188 91.98415256 90.64261126 90.24910835]
    ini   pot:-174.9558,maxforce:0.2321
    50 pot:-175.0080,maxforce:0.0270,delta:-0.0522
    100 pot:-175.0092,maxforce:0.0048,delta:-0.0012
    Optimized lattice constant = [ 8.70048578  8.78800588 13.07175459 91.32661125 90.54564294 90.13344282]



```python
# Remove comment out below if you want to run MD with bigger systems.
#bulk = bulk.repeat([2,2,1])
```

構造緩和されたLGPSの構造の可視化


```python
v = view(bulk, viewer='ngl')
#v.view.add_representation("ball+stick")
display(v)
```


    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'P', 'Li', 'Ge', 'S'),…



```python
os.makedirs(OUTPUT_DIR / "structure/", exist_ok=True)
write(OUTPUT_DIR / "structure/opt_structure.xyz", bulk)
```

Check number of Li in this structure


```python
Li_index = [i for i, x in enumerate(bulk.get_chemical_symbols()) if x == 'Li']
print(len(Li_index))
```

    20


## 様々な温度でのMD シミュレーションの実行

今回は NVT のMD シミュレーションを行うため、ASEの `Langevin` クラスを利用してMD シミュレーションを実行します。

 - https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin

以下では、様々な温度域でのMDシミュレーションが実行されています。&lt;br/&gt;
次の解析スクリプトで行うように、Li拡散の温度依存性を調べることで、 Arrhenius plot を描画することができます。&lt;br/&gt;
MDシミュレーションは、`joblib`を用いて並列化されています。


```python
temp_list = [423, 523, 623, 723, 823, 923, 973, 1023]
```


```python
os.makedirs(OUTPUT_DIR / "traj_and_log/", exist_ok=True)

def run_md(i):
    s_time = perf_counter()
    estimator = Estimator(model_version="v2.0.0", calc_mode=EstimatorCalcMode.PBE_U)
    calculator = ASECalculator(estimator)

    t_step = 0.5     # as fs
    temp = i       # as K
    itrvl = 100

    structure = read(f"{OUTPUT_DIR.name}/" + "structure/opt_structure.xyz")
    structure.calc = calculator

    MaxwellBoltzmannDistribution(structure, temperature_K=temp)

    dyn = Langevin(
        structure,
        t_step * units.fs,
        temperature_K=temp,
        friction=0.02,
        trajectory=f"{OUTPUT_DIR.name}/" + "traj_and_log/MD_"+str(i).zfill(4)+".traj",
        loginterval=itrvl,
        append_trajectory=False,
    )
    dyn.attach(MDLogger(dyn, structure, f"{OUTPUT_DIR.name}/" + "traj_and_log/MD_"+str(i).zfill(4)+".log", header=False, stress=False,
               peratom=True, mode="w"), interval=itrvl)

    # dyn.run(500000)
    dyn.run(2_000_000)
    proctime = perf_counter() - s_time

    return([i, proctime/3600])
```


```python
results = Parallel(n_jobs=len(temp_list), verbose=1)(delayed(run_md)(i) for i in temp_list)
```

    [Parallel(n_jobs=8)]: Using backend LokyBackend with 8 concurrent workers.
    /home/jovyan/.local/lib/python3.7/site-packages/joblib/externals/loky/process_executor.py:705: UserWarning: A worker stopped while some jobs were given to the executor. This can be caused by a too short worker timeout or by a memory leak.
      "timeout or by a memory leak.", UserWarning


To check the calculation time (hours) for each job:


```python
# 計算時間
print(results)
```

    [[423, 69.65517497624501], [523, 69.4307892338875], [623, 69.65359850568028], [723, 69.5858971493625], [823, 69.63537499042417], [923, 69.65377913352611], [973, 69.64277460395472], [1023, 69.59344196447528]]
