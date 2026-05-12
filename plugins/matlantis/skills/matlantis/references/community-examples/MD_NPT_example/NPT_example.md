Copyright Preferred Computational Chemistry, Inc. and Preferred Networks, Inc. as contributors to Matlantis contrib project

# NPT Example

## Setup

SiO2のCIFファイルはMaterials Projectからとってきました。&lt;br/&gt;
Launcher menuから"Import material from Material Project" をご参照ください。

 - https://materialsproject.org/materials/mp-6945/
 - https://materialsproject.org/materials/mp-546794/

Input cif files are from  
A. Jain*, S.P. Ong*, G. Hautier, W. Chen, W.D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, K.A. Persson (*=equal contributions)  
The Materials Project: A materials genome approach to accelerating materials innovation  
APL Materials, 2013, 1(1), 011002.  
[doi:10.1063/1.4812323](http://dx.doi.org/10.1063/1.4812323)  
[[bibtex]](https://materialsproject.org/static/docs/jain_ong2013.349ca3156250.bib)  
Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)  


```python
!pip install pfp_api_client
!pip install matplotlib ase

# # 初回使用時のみ、ライブラリのインストールをお願いします。
```


```python
import ase
from ase.build import molecule
from ase.visualize import view

ase.__version__
```




    '3.22.0'




```python
from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

estimator = Estimator()
calculator = ASECalculator(estimator)
```


```python
import os

os.makedirs("output", exist_ok=True)
```

今回計算対象とするSiO2 構造


```python
import ase.io
from ase.visualize import view

atoms = ase.io.read("input/SiO2_mp-6945_computed.cif")
atoms *= (4, 3, 3)

v = view(atoms, viewer='ngl')
v.view.add_ball_and_stick()
display(v)
```


    



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Si', 'O'), value='All…


## NPTでの等方制御
### 方法1: mask指定 (Cell角度固定)

NPT クラスの mask パラメータに(3, 3) で指定することで、Cell の変化する方向の指定ができます。&lt;br/&gt;
`np.eye(3)` を指定すると、角度固定でMDができます。

 - https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.npt



```python
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase.md.npt import NPT
from ase.md.nptberendsen import NPTBerendsen
from ase.io import Trajectory
from ase import units

from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

from time import perf_counter
import numpy as np

atoms = ase.io.read("input/SiO2_mp-6945_computed.cif")
atoms *= (4, 3, 3)
atoms.calc = calculator


# ase version&gt;=3.21.0では temperature_K 指定が推奨されている
MaxwellBoltzmannDistribution(
    atoms,
    temperature_K=300.0
)
Stationary(atoms)
dyn = NPT(atoms,
          1.0 * units.fs,
          temperature_K=300,
          externalstress=101325 * units.Pascal,
          mask=np.eye(3),  # ここに np.eye(3) を指定する
          ttime=10 * units.fs,
          pfactor=100 * units.fs,
          logfile="output/sio2.log",
          trajectory="output/sio2.traj")


start_time = perf_counter()

def print_dyn():
    print(f"Dyn  step: {dyn.get_number_of_steps(): &gt;3}, energy: {atoms.get_total_energy():.3f}, elapsed_time: {perf_counter() - start_time:.2f} sec")

dyn.attach(print_dyn, interval=10)
dyn.run(100)
```

    Dyn  step:  10, energy: -2720.284, elapsed_time: 3.49 sec
    Dyn  step:  20, energy: -2710.637, elapsed_time: 6.69 sec
    Dyn  step:  30, energy: -2694.522, elapsed_time: 9.88 sec
    Dyn  step:  40, energy: -2691.940, elapsed_time: 13.04 sec
    Dyn  step:  50, energy: -2703.926, elapsed_time: 16.20 sec
    Dyn  step:  60, energy: -2708.250, elapsed_time: 19.36 sec
    Dyn  step:  70, energy: -2707.254, elapsed_time: 22.53 sec
    Dyn  step:  80, energy: -2704.370, elapsed_time: 25.72 sec
    Dyn  step:  90, energy: -2703.853, elapsed_time: 28.93 sec
    Dyn  step: 100, energy: -2705.025, elapsed_time: 32.13 sec


**Cell sizeの確認**

a,b,c 軸の長さは変わりつつも、角度は変化していません。


```python
traj = Trajectory("output/sio2.traj")
cell_size_history = [atoms.cell.cellpar() for atoms in traj]

# ｐｒｉｎｔ量を減らすため、10step毎に確認
cell_size_history[::10]
```




    [array([20.338804, 15.254103, 21.295737, 90.      , 90.      , 90.      ]),
     array([20.34709056, 15.25997189, 21.30078079, 90.        , 90.        ,
            90.        ]),
     array([20.36662827, 15.27483021, 21.31356124, 90.        , 90.        ,
            90.        ]),
     array([20.38826633, 15.29445575, 21.33151445, 90.        , 90.        ,
            90.        ]),
     array([20.4053848 , 15.31091376, 21.34388992, 90.        , 90.        ,
            90.        ]),
     array([20.40853203, 15.31629593, 21.33536484, 90.        , 90.        ,
            90.        ]),
     array([20.38617267, 15.30688091, 21.3034135 , 90.        , 90.        ,
            90.        ]),
     array([20.34711366, 15.28739519, 21.25958074, 90.        , 90.        ,
            90.        ]),
     array([20.31370621, 15.27012572, 21.22078718, 90.        , 90.        ,
            90.        ]),
     array([20.29884377, 15.26646011, 21.20456456, 90.        , 90.        ,
            90.        ])]




```python
v = view(traj, viewer='ngl')
v.view.add_ball_and_stick()
display(v)
```


    HBox(children=(NGLWidget(max_frame=99), VBox(children=(Dropdown(description='Show', options=('All', 'Si', 'O')…


### 方法2: Filterを活用 (セルの長さ比率を保つ)

Stressを後処理で同じ大きさになるようにすることにより、セルの長さの比率を保つことができます。方法を以下に示します。


```python
import numpy as np
from ase.constraints import Filter



class IsotropicFilter(Filter):
    def __init__(self, atoms, indices=None, mask=None):
        if indices is None and mask is None:
            # Apply filter for all atoms
            indices = np.arange(atoms.get_global_number_of_atoms())
        super(IsotropicFilter, self).__init__(atoms, indices, mask)

    def get_stress(self, **kwargs):
        cell = self.atoms.get_cell()
        stresses = self.atoms.get_stress(**kwargs)
        if stresses.shape != (6,):
            raise ValueError(f"Error: stresses.shape {stresses.shape}")
        stresses[3:] = 0.0
        pressure = np.mean(stresses[:3])
        stresses[:3] = pressure
        return stresses

    def get_temperature(self):
        return self.atoms.get_temperature()

    def get_total_energy(self):
        return self.atoms.get_total_energy()

    def get_kinetic_energy(self):
        return self.atoms.get_kinetic_energy()

    def get_global_number_of_atoms(self):
        return self.atoms.get_global_number_of_atoms()

    def set_cell(self, cell, scale_atoms=False, apply_constraint=True):
        return self.atoms.set_cell(
            cell, scale_atoms=scale_atoms, apply_constraint=apply_constraint
        )
```


```python
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase.md.npt import NPT
from ase.md.nptberendsen import NPTBerendsen
from ase.io import Trajectory
from ase.build import bulk
from ase import units

from pfp_api_client.pfp.estimator import Estimator
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

from time import perf_counter
import numpy as np

atoms = bulk("Si", cubic = True)
atoms *= (4, 3, 3)
atoms.calc = calculator

# Filterをatomsに適用
iso_atoms = IsotropicFilter(atoms)
MaxwellBoltzmannDistribution(iso_atoms, temperature_K=300.0)
Stationary(iso_atoms)
dyn = NPT(iso_atoms,
          1.0 * units.fs,
          temperature_K=300,
          externalstress=101325 * units.Pascal,
          mask=np.eye(3),
          ttime=10 * units.fs,
          pfactor=100 * units.fs,
          logfile="output/si_iso.log",
          trajectory="output/si_iso.traj")


start_time = perf_counter()
def print_dyn():
    print(f"Dyn  step: {dyn.get_number_of_steps(): &gt;3}, energy: {atoms.get_total_energy():.3f}, elapsed_time: {perf_counter() - start_time:.2f} sec")
dyn.attach(print_dyn, interval=10)
dyn.run(100)
```

    Dyn  step:  10, energy: -1298.112, elapsed_time: 2.18 sec
    Dyn  step:  20, energy: -1294.476, elapsed_time: 4.16 sec
    Dyn  step:  30, energy: -1279.642, elapsed_time: 6.09 sec
    Dyn  step:  40, energy: -1281.333, elapsed_time: 7.98 sec
    Dyn  step:  50, energy: -1286.203, elapsed_time: 9.81 sec
    Dyn  step:  60, energy: -1282.718, elapsed_time: 11.63 sec
    Dyn  step:  70, energy: -1273.227, elapsed_time: 13.45 sec
    Dyn  step:  80, energy: -1280.578, elapsed_time: 15.30 sec
    Dyn  step:  90, energy: -1282.381, elapsed_time: 17.12 sec
    Dyn  step: 100, energy: -1284.916, elapsed_time: 18.98 sec


Cell 長さ `a:b:c` の比率を確認します。比率が変わらずに保たれていることがわかります。


```python
traj = Trajectory("output/si_iso.traj")
cell_size_history = [atoms.cell.cellpar() for atoms in traj]

cell_length_history = np.array([atoms.get_cell_lengths_and_angles()[:3] for atoms in traj])

# ｐｒｉｎｔ量を減らすため、10step毎に確認
cell_length_ratio_history = cell_length_history / cell_length_history[:, 0:1]
print("a:b:c = \n", cell_length_ratio_history[::10])
```

    a:b:c = 
     [[1.   0.75 0.75]
     [1.   0.75 0.75]
     [1.   0.75 0.75]
     [1.   0.75 0.75]
     [1.   0.75 0.75]
     [1.   0.75 0.75]
     [1.   0.75 0.75]
     [1.   0.75 0.75]
     [1.   0.75 0.75]
     [1.   0.75 0.75]]


    /home/jovyan/.local/lib/python3.7/site-packages/ase/utils/__init__.py:62: FutureWarning: Please use atoms.cell.cellpar() instead
      warnings.warn(warning)


 実際のCell lengthは変化しています。


```python
print("a,b,c = \n", cell_length_history[::10])
```

    a,b,c = 
     [[21.72       16.29       16.29      ]
     [21.7356418  16.30173135 16.30173135]
     [21.77517368 16.33138026 16.33138026]
     [21.83240146 16.37430109 16.37430109]
     [21.89665807 16.42249356 16.42249356]
     [21.95295893 16.4647192  16.4647192 ]
     [21.98924709 16.49193532 16.49193532]
     [21.99966685 16.49975014 16.49975014]
     [21.98592141 16.48944105 16.48944105]
     [21.95324602 16.46493451 16.46493451]]


## 密度・体積を後処理で得る方法

Trajectory instanceから、各Frame時のAtoms情報が得られます。&lt;br/&gt;
Ase atoms に `get_volume`, `get_global_number_of_atoms` 関数があるので、これが使用できます。


```python
# cell_list = [atoms.cell for atoms in traj]

volume_list = [atoms.get_volume() for atoms in traj]
density_list = [atoms.get_global_number_of_atoms() / atoms.get_volume() for atoms in traj]
```


```python
print(volume_list)
print(density_list)
```

    [5763.708251999998, 5764.131707178368, 5764.553834457096, 5765.395785246378, 5766.235719316216, 5767.48975024803, 5768.742057145161, 5770.400739141297, 5772.058871091771, 5774.113555362282, 5776.169542055843, 5778.609989119666, 5781.053978471361, 5783.867838795086, 5786.687509375583, 5789.859725415075, 5793.039666921781, 5796.551906389886, 5800.0731450987, 5803.903457723219, 5807.743265920976, 5811.866075304465, 5815.998187935012, 5820.384763044434, 5824.780054012687, 5829.399356050736, 5834.026786906732, 5838.846429752368, 5843.673863376523, 5848.660450369243, 5853.654057532671, 5858.771470356455, 5863.893017166693, 5869.099187710506, 5874.30307612228, 5879.548151579635, 5884.781362135775, 5890.009561282161, 5895.215120088264, 5900.369365913755, 5905.49133912333, 5910.517666844609, 5915.504455254107, 5920.354094271415, 5925.159373090376, 5929.788761936238, 5934.370890390756, 5938.740762714563, 5943.061813055538, 5947.136302720512, 5951.161299392146, 5954.907390782732, 5958.60401212543, 5961.991460893361, 5965.3301964442935, 5968.33191155933, 5971.28665623559, 5973.879492119523, 5976.428553055096, 5978.594478313054, 5980.721850087601, 5982.449259471617, 5984.145847803037, 5985.43056227638, 5986.694841459074, 5987.540204958565, 5988.377580194568, 5988.793009284016, 5989.213615428513, 5989.21185008912, 5989.22790757915, 5988.822903005305, 5988.44732378247, 5987.65380274905, 5986.900584732638, 5985.734892207617, 5984.620098306944, 5983.100878268245, 5981.64275942176, 5979.790520462952, 5978.008665616815, 5975.844704818999, 5973.758972443026, 5971.304402240134, 5968.934168668893, 5966.2094088043505, 5963.573317154438, 5960.598091307133, 5957.714288408501, 5954.508058362449, 5951.3947803054625, 5947.977367668708, 5944.653567226575, 5941.04566665038, 5937.5313149813255, 5933.754546868684, 5930.070469849946, 5926.146955468612, 5922.314177233043, 5918.265860239761]
    [0.04996783102268656, 0.04996416019455955, 0.04996050141443839, 0.049953205422078865, 0.0499459290287481, 0.04993506923659719, 0.04992422908618064, 0.0499098785369381, 0.04989554099012949, 0.04987778595600034, 0.04986003231087563, 0.049838975210693345, 0.04981790536336656, 0.04979366887815976, 0.049769405991490434, 0.049742137747448324, 0.04971483306846286, 0.049684709919102135, 0.049654546209192525, 0.0496217764643828, 0.04958896886677888, 0.04955379154790875, 0.049518584891831835, 0.049481264851871674, 0.049443927037484785, 0.04940474694036272, 0.04936556010444733, 0.0493248115813545, 0.04928406456851636, 0.04924204481417924, 0.04920003764646671, 0.04915706329512759, 0.04911412932617849, 0.0490705627540002, 0.04902709245129949, 0.04898335596122709, 0.04893979610747604, 0.048896355261146135, 0.04885317908393613, 0.04881050357012678, 0.04876816905852131, 0.04872669641367501, 0.048685619659064, 0.04864573899028628, 0.048606287504767705, 0.04856834055349387, 0.04853083929525616, 0.04849512910348976, 0.04845986951832309, 0.04842666879322317, 0.048393915995759755, 0.048363472527848055, 0.048333468613443, 0.04830600679136922, 0.04827897040329232, 0.04825468895960831, 0.04823081131086755, 0.04820987774860822, 0.048189315314875976, 0.048171857289317825, 0.04815472232599843, 0.04814081783377078, 0.04812716924433478, 0.04811683921540103, 0.04810667782923253, 0.04809988578640251, 0.0480931598155243, 0.04808982370129228, 0.04808644648407558, 0.048086460657709836, 0.048086331734938065, 0.048089583656827806, 0.0480925997054761, 0.0480989732351883, 0.04810502461564784, 0.04811439283335548, 0.048123355412564205, 0.04813557482309391, 0.048147308621259204, 0.04816222224080572, 0.048176577872237655, 0.04819402347717521, 0.04821085037554162, 0.04823066797464836, 0.048249820128980526, 0.04827185575735871, 0.04829319347371104, 0.04831729896703082, 0.048340686722816, 0.04836671597001801, 0.048392017439854336, 0.04841982109170007, 0.048446893791720784, 0.04847631480375024, 0.0485050073375328, 0.04853588022982534, 0.04856603331516354, 0.04859818734907263, 0.0486296389183723, 0.04866290342494558]


可視化して変化履歴を確認


```python
energy_list = [atoms.get_potential_energy() for atoms in traj]
temperature_list = [atoms.get_temperature() for atoms in traj]
```


```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(4, 1, sharex=True, figsize=(10, 10))

axes[0].plot(volume_list)
axes[0].set_title("Volume history")

axes[1].plot(density_list)
axes[1].set_title("Density history")

axes[2].plot(energy_list)
axes[2].set_title("Energy history")

axes[3].plot(temperature_list)
axes[3].set_title("Temperature history")
axes[3].set_xlabel("step")

fig.savefig("output/md_stats_history.png")
```


    
![png](output_25_0.png)
    


※ 100step程度だと温度は安定しません。

## 三斜晶系でのNPT MD

aseのNPT moduleは `atoms` のcell が下記のように上三角行列である必要があるようです。
```
[[X X X]
 [0 X X]
 [0 0 X]]
```

`convert_atoms_to_upper`関数は `atoms`を回転させてCellを上三角行列にします。


```python
from ase import Atoms
import numpy as np


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


```python
atoms = ase.io.read("input/SiO2_mp-546794_computed.cif")
atoms.cell  # NOT upper tiangular
```




    Cell([[5.13584237, 0.0, 0.0], [0.15785281555277758, 5.133415952181283, 0.0], [-2.646847592991795, -2.566707976299526, 3.575344016085885]])




```python
atoms = convert_atoms_to_upper(atoms)
atoms.cell  # This is upper tiangular
```




    Cell([[4.170109258320234, -1.407573601142272, -2.646847592991795], [0.0, 4.401258305185092, -2.646847592991795], [0.0, 0.0, 5.135842369999999]])




```python
atoms = ase.io.read("input/SiO2_mp-546794_computed.cif") * (2, 2, 3)
# ↓をコメントアウトするとNPT Module がエラーとなります。
# NotImplementedError: Can (so far) only operate on lists of atoms where the computational box is an upper triangular matrix.
atoms = convert_atoms_to_upper(atoms)
atoms.calc = calculator
MaxwellBoltzmannDistribution(atoms, temperature_K=300.0)
Stationary(atoms)
dyn = NPT(atoms,
          1.0 * units.fs,
          temperature_K=300,
          externalstress=101325 * units.Pascal,
          mask=np.eye(3),
          ttime=10 * units.fs,
          pfactor=100 * units.fs,
          logfile="output/si_tri.log",
          trajectory="output/si_tri.traj")


start_time = perf_counter()
def print_dyn():
    print(f"Dyn  step: {dyn.get_number_of_steps(): &gt;3}, energy: {atoms.get_total_energy():.3f}, elapsed_time: {perf_counter() - start_time:.2f} sec")
dyn.attach(print_dyn, interval=10)
dyn.run(100)
```

    Dyn  step:  10, energy: -453.486, elapsed_time: 1.29 sec
    Dyn  step:  20, energy: -451.625, elapsed_time: 2.46 sec
    Dyn  step:  30, energy: -448.398, elapsed_time: 3.67 sec
    Dyn  step:  40, energy: -448.321, elapsed_time: 4.85 sec
    Dyn  step:  50, energy: -451.135, elapsed_time: 6.00 sec
    Dyn  step:  60, energy: -452.141, elapsed_time: 7.19 sec
    Dyn  step:  70, energy: -451.945, elapsed_time: 8.34 sec
    Dyn  step:  80, energy: -450.731, elapsed_time: 9.50 sec
    Dyn  step:  90, energy: -450.171, elapsed_time: 10.66 sec
    Dyn  step: 100, energy: -450.561, elapsed_time: 11.84 sec


三斜晶系でもMDを行うことができました。

## 可視化方法

### nglview

nglviewer を用いるとTrajectory をInteractiveに可視化できます。

※: nglviewerでは、はじめのフレームの状態の結合・Cellの状態を保持し続け、また原子の追加・削除の描画もされないようです。


```python
traj = Trajectory("output/si_iso.traj")
v = view(traj[::10], viewer='ngl')
v.view.add_ball_and_stick()
display(v)
```


    HBox(children=(NGLWidget(max_frame=9), VBox(children=(Dropdown(description='Show', options=('All', 'Si'), valu…


### gifファイルとして保存

nglviewではなく、ase デフォルトの機能をもちいてgifの作成も可能です。&lt;br/&gt;
可視化の負荷を下げるために10 stepごとの描画をしています。


```python
%%time
traj = Trajectory("output/si_iso.traj")
traj_small = [atoms for atoms in traj[::10]]
ase.io.write("output/si_iso.gif", traj_small, rotation="-90x,0y,0z")
```

    MovieWriter ffmpeg unavailable; using Pillow instead.


    CPU times: user 7.77 s, sys: 276 ms, total: 8.05 s
    Wall time: 8.03 s



    
![png](output_37_2.png)
    



```python
from IPython.display import Image
Image(url='output/si_iso.gif')  
```




&lt;img src="output/si_iso.gif"/&gt;



### 別の可視化ツールを使う

traj ファイルを別形式で保存します。&lt;br/&gt;
ここでは、xyz, pdbファイルへの変換を行なっています。保存したファイルをダウンロードすることで、[VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)や[OVITO](https://www.ovito.org/)などの別ソフトで可視化できます。


```python
from ase.io import write

traj = Trajectory("output/sio2.traj")
write("output/sio2.xyz", traj)
write("output/sio2.pdb", traj)
```
