Copyright Preferred Computational Chemistry, Inc. as contributors to Matlantis contrib project

# MD: LGPS中のLi拡散 - 解析スクリプト

本スクリプトでは、LGPS中のLi拡散 MD計算結果の解析を行います。

## Setup


```python
# Please install these libraries only for first time of execution
!pip install pandas matplotlib scipy ase
```




```python
import pathlib
EXAMPLE_DIR = pathlib.Path("__file__").resolve().parent
INPUT_DIR = EXAMPLE_DIR / "input"
OUTPUT_DIR = EXAMPLE_DIR / "output"
```


```python
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
import glob, os

import ase
from ase.visualize import view
from ase.io import read, write
from ase.io import Trajectory
```

## 特定温度の解析

以下では、MD計算結果として得られたtrajectory ファイルを20 step間隔で読み込んでいます。&lt;br/&gt;
Trajectory fileは 0.05 ps 間隔で保存されていたため、 `trj` は 1ps 間隔となります。


```python
temp = 523
#temp = 423

trj = read(OUTPUT_DIR / f"traj_and_log/MD_{temp:04}.traj", index="::20")
view(trj, viewer = "ngl")
```






    HBox(children=(NGLWidget(max_frame=1000), VBox(children=(Dropdown(description='Show', options=('All', 'P', 'Ge…



```python
Li_index = [i for i, x in enumerate(trj[0].get_chemical_symbols()) if x == 'Li']
print(len(Li_index))
```

    20


拡散係数の計算のために、MSD (mean square displacement) を計算します。


```python
# t0 = len(trj) // 2
t0 = 0

positions_all = np.array([trj[i].get_positions() for i in range(t0, len(trj))])

# shape is (n_traj, n_atoms, 3 (xyz))
print("positions_all.shape: ", positions_all.shape)

# position of Li
positions = positions_all[:, Li_index]
positions_x = positions[:, :, 0]
positions_y = positions[:, :, 1]
positions_z = positions[:, :, 2]

print("positions.shape    : ", positions.shape)
print("positions_x.shape  : ", positions_x.shape)
```

    positions_all.shape:  (1001, 50, 3)
    positions.shape    :  (1001, 20, 3)
    positions_x.shape  :  (1001, 20)



```python
# msd for each x,y,z axis
msd_x = np.mean((positions_x-positions_x[0])**2, axis=1)
msd_y = np.mean((positions_y-positions_y[0])**2, axis=1)
msd_z = np.mean((positions_z-positions_z[0])**2, axis=1)

# total msd. sum along xyz axis &amp; mean along Li atoms axis.
msd = np.mean(np.sum((positions-positions[0])**2, axis=2), axis=1)
```

まずx, y, z 軸方向それぞれの拡散をプロットした結果として、z軸方向の拡散が x軸・y軸方向と比べて大きいという、既知の結果と整合する事実が確認できます。


```python
plt.plot(range(len(msd_x)), msd_x, label="x")
plt.plot(range(len(msd_y)), msd_y, label="y")
plt.plot(range(len(msd_z)), msd_z, label="z")
plt.grid(True)
#plt.xlim(0,100)
#plt.ylim(0,10)
plt.xlabel("time (psec)")
plt.ylabel("MSD (A^2)")
plt.title(f"xyz MSD at {temp}K")
plt.legend()
plt.show()
```



![png](output_12_0.png)



次にtotal msd を線形Fittingすることで、その傾きから拡散係数を得ることができます。


```python
plt.plot(range(len(msd)), msd)
plt.grid(True)
#plt.xlim(0,100)
#plt.ylim(0,10)
plt.xlabel("time (psec)")
plt.ylabel("MSD (A^2)")
plt.title(f"MSD at {temp}K")
plt.show()
```



![png](output_14_0.png)



線形Fittingには `scipy.stats.linregress` を用いました。
 - https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html


ここで 6 という係数で割っているのは、xyz軸方向の自由度を考慮しています。


```python
slope, intercept, r_value, _, _ = stats.linregress(range(len(msd)), msd)
D = slope / 6
print(slope, intercept, r_value)
print(f"Diffusion coefficient {D:.2f} A^2/psec")
```

    2.733539684523921 -32.93407675583876 0.9884154195683691
    Diffusion coefficient 0.46 A^2/psec



```python
t = np.arange(len(msd))
#plt.scatter(t, msd, label="MSD")
plt.plot(t, msd, label="MSD")
plt.plot(t, t * slope + intercept, label="fitted line")
plt.grid(True)
plt.legend()
#plt.xlim(0,100)
#plt.ylim(0,10)
plt.xlabel("time (psec)")
plt.ylabel("MSD (A^2)")
plt.title(f"MSD at {temp}K")
plt.show()
```



![png](output_17_0.png)



Convert diffusion coefficient from A^2/ps to cm^2/sec unit, take log.


```python
np.log10(D*1e-16*1e12)
```




    -5.193181079582682



## 温度依存性の解析

拡散係数の温度依存性解析を行うことで、アレニウスの式から活性化エネルギー $E_A$を算出することができます。

$$D = D_0 \exp \left(- \frac{E_A}{RT} \right)$$

ここで、 $D$ は温度$T$ における拡散係数です。

 - [アレニウスの式](https://ja.wikipedia.org/wiki/%E3%82%A2%E3%83%AC%E3%83%8B%E3%82%A6%E3%82%B9%E3%81%AE%E5%BC%8F)

両辺のlog を取ると、

$$\log D = \log D_0 - \frac{E_A}{RT}$$

となるため、$\log D$ を y軸に、 $1/T$ をx軸にとってプロットすることで、その傾きから活性化エネルギー $E_A$ を計算できます。

まずは前節同様に、それぞれの温度での拡散係数の計算を行います。


```python
trj_list = sorted(glob.glob(f"{OUTPUT_DIR}/traj_and_log/*.traj"))
```


```python
trj_list
```




    ['/home/jovyan/MD_Li-diffusion_in_LGPS/output/traj_and_log/MD_0423.traj',
     '/home/jovyan/MD_Li-diffusion_in_LGPS/output/traj_and_log/MD_0523.traj',
     '/home/jovyan/MD_Li-diffusion_in_LGPS/output/traj_and_log/MD_0623.traj',
     '/home/jovyan/MD_Li-diffusion_in_LGPS/output/traj_and_log/MD_0723.traj',
     '/home/jovyan/MD_Li-diffusion_in_LGPS/output/traj_and_log/MD_0823.traj',
     '/home/jovyan/MD_Li-diffusion_in_LGPS/output/traj_and_log/MD_0923.traj',
     '/home/jovyan/MD_Li-diffusion_in_LGPS/output/traj_and_log/MD_0973.traj',
     '/home/jovyan/MD_Li-diffusion_in_LGPS/output/traj_and_log/MD_1023.traj']




```python
t0 = 0

os.makedirs(OUTPUT_DIR / "msd/", exist_ok=True)
D_list = []
for path in trj_list:
    trj = read(path, index="::20")
    Li_index = [Li_i for Li_i, x in enumerate(trj[0].get_chemical_symbols()) if x == 'Li']

    # msd for each x,y,z axis
    positions_all = np.array([trj[i].get_positions() for i in range(t0, len(trj))])
    positions = positions_all[:, Li_index]
    positions_x = positions[:, :, 0]
    positions_y = positions[:, :, 1]
    positions_z = positions[:, :, 2]
    # msd for each x,y,z axis
    msd_x = np.mean((positions_x-positions_x[0])**2, axis=1)
    msd_y = np.mean((positions_y-positions_y[0])**2, axis=1)
    msd_z = np.mean((positions_z-positions_z[0])**2, axis=1)

    # total msd. sum along xyz axis &amp; mean along Li atoms axis.
    msd = np.mean(np.sum((positions-positions[0])**2, axis=2), axis=1)

    slope, intercept, r_value, _, _ = stats.linregress(range(len(msd)), msd)
    logD = np.log10(slope*1e-16*1e12/6)
    T = int(os.path.basename(path.split(".")[0].replace("MD_","").replace("traj_and_log/","")))
    D_list.append([T, 1000/T, logD])

    fig = plt.figure(figsize=(10,4), facecolor='w')
    ax1 = fig.add_subplot(1,2,1)
    ax1.plot(range(len(msd_x)), msd_x, label="x")
    ax1.plot(range(len(msd_y)), msd_y, label="y")
    ax1.plot(range(len(msd_z)), msd_z, label="z")
    ax1.set_xlabel("time (psec)")
    ax1.set_ylabel("MSD (A^2)")
    ax1.legend()
    ax1.set_title(f"xyz MSD at {T} K")
    ax1.grid(True)

    #fig = plt.figure()
    ax2 = fig.add_subplot(1,2,2)
    ax2.plot(range(len(msd)), msd)
    ax2.set_xlabel("time (psec)")
    ax2.set_ylabel("MSD (A^2)")
    ax2.set_title(f"MSD at {T} K")
    ax2.grid(True)

    plt.show(fig)
    fig.savefig(path.replace("traj_and_log/", "msd/").replace("traj", "png"))
    plt.close(fig)
```



![png](output_24_0.png)





![png](output_24_1.png)





![png](output_24_2.png)





![png](output_24_3.png)





![png](output_24_4.png)





![png](output_24_5.png)





![png](output_24_6.png)





![png](output_24_7.png)




```python
df = pd.DataFrame(D_list, columns=["T", "1000/T", "logD"])
```


```python
df
```




&lt;div&gt;
&lt;style scoped&gt;
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
&lt;/style&gt;
&lt;table border="1" class="dataframe"&gt;
  &lt;thead&gt;
    &lt;tr style="text-align: right;"&gt;
      &lt;th&gt;&lt;/th&gt;
      &lt;th&gt;T&lt;/th&gt;
      &lt;th&gt;1000/T&lt;/th&gt;
      &lt;th&gt;logD&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;423&lt;/td&gt;
      &lt;td&gt;2.364066&lt;/td&gt;
      &lt;td&gt;-5.802397&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;523&lt;/td&gt;
      &lt;td&gt;1.912046&lt;/td&gt;
      &lt;td&gt;-5.193181&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;623&lt;/td&gt;
      &lt;td&gt;1.605136&lt;/td&gt;
      &lt;td&gt;-4.756214&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;723&lt;/td&gt;
      &lt;td&gt;1.383126&lt;/td&gt;
      &lt;td&gt;-4.583190&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;823&lt;/td&gt;
      &lt;td&gt;1.215067&lt;/td&gt;
      &lt;td&gt;-4.332002&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;5&lt;/th&gt;
      &lt;td&gt;923&lt;/td&gt;
      &lt;td&gt;1.083424&lt;/td&gt;
      &lt;td&gt;-4.343488&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;6&lt;/th&gt;
      &lt;td&gt;973&lt;/td&gt;
      &lt;td&gt;1.027749&lt;/td&gt;
      &lt;td&gt;-4.290990&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;7&lt;/th&gt;
      &lt;td&gt;1023&lt;/td&gt;
      &lt;td&gt;0.977517&lt;/td&gt;
      &lt;td&gt;-4.341426&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;




```python
sl, ic, rv, _, _ = stats.linregress(df["1000/T"], df["logD"])
print(sl, ic, rv)
```

    -1.0856951904952106 -3.1354278416677452 -0.9834568664784172


アレニウスプロットは、以下の参考文献に合わせたスケール軸をとっています。
 - [First Principles Study of the Li10GeP2S12 Lithium Super Ionic Conductor Material](https://pubs.acs.org/doi/10.1021/cm203303y)


```python
fig = plt.figure()
plt.scatter(df["1000/T"], df["logD"])
plt.plot(df["1000/T"], df["1000/T"]*sl+ic)

plt.grid(True)
plt.xlabel("1000/T (1/K)")
plt.ylabel("log(D(cm^2/sec))")
plt.title("Arrhenius plot")
fig.savefig(OUTPUT_DIR / "arrhenius_plot.png")
```



![png](output_29_0.png)



結果として、以下のように活性化エネルギーが算出できます。

以下の式で、 `1000 * np.log(10)` の項は、上記アレニウスプロットの x軸、y軸のスケール補正を行う項です。


```python
from ase.units import J, mol
R = 8.31446261815324  # J/(K・mol)


E_act = -sl * 1000 * np.log(10) * R * (J / mol)
print(f"Activation energy: {E_act* 1000:.1f} meV")
```

    Activation energy: 215.4 meV

