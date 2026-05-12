# Appendix 2: Opt algorithm fmax experiment

本節では構造最適化のパラメーターであるfmaxについて考察します。


```python
!pip install more_itertools
```




```python
from ase import Atoms
from ase.constraints import ExpCellFilter, FixAtoms
from ase.build import bulk, molecule, add_adsorbate, surface
from ase.optimize import FIRE
from ase.units import kB
import pandas as pd
from tqdm.auto import tqdm
import numpy as np
from more_itertools import windowed
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfcc_extras.structure.ase_rdkit_converter import smiles_to_atoms, atoms_to_smiles
from pfcc_extras.visualize.view import view_ngl
from matlantis_features.features.common.fire_lbfgs import FIRELBFGS
```





    /home/jovyan/.py311/lib/python3.11/site-packages/matlantis_features/features/common/fire_lbfgs.py:10: FutureWarning: matlantis_features.features.common.fire_lbfgs is deprecated; use ase_ext.optimize
      warnings.warn(



```python
def max_distance(atoms1: Atoms, atoms2: Atoms) -&gt; float:
    return float(np.max(np.linalg.norm(atoms1.positions - atoms2.positions, axis=1)))
```


```python
calc_mol = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.WB97XD, model_version="v8.0.0"))
calc = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE, model_version="v8.0.0"))
calc_d3 = ASECalculator(Estimator(calc_mode=EstimatorCalcMode.PBE_PLUS_D3, model_version="v8.0.0"))
```

## fmax=0.05eV/Aで十分な例

ASEにおいてfmaxのデフォルト値は0.05eV/Aですが、これは多くの場合において十分な精度を持っています。
分子と固体の例でそれぞれ確認してみましょう。

### 分子の例


```python
fmax_list = [0.05, 0.01, 0.001]
```


```python
atoms = molecule("C2H6")
atoms.rattle(stdev=0.1)
atoms.calc = calc_mol
images = []
for fmax in tqdm(fmax_list):
    with FIRE(atoms, logfile=None) as opt:
        opt.run(fmax=fmax)
        images.append(atoms.copy())
```


      0%|          | 0/3 [00:00&lt;?, ?it/s]



```python
for atoms in images:
    atoms.calc = calc_mol

results = {
    "fmax": [],
    "energy (eV)": [],
    "distance (A)": [],
}
for (fmax_old, old), (fmax_new, new) in windowed(zip(fmax_list, images), 2):
    de = new.get_potential_energy() - old.get_potential_energy()
    dx = max_distance(new, old)
    results["energy (eV)"].append(de)
    results["distance (A)"].append(dx)
    results["fmax"].append(f"fmax {fmax_old}-&gt;{fmax_new}")
pd.DataFrame(results)
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
      &lt;th&gt;fmax&lt;/th&gt;
      &lt;th&gt;energy (eV)&lt;/th&gt;
      &lt;th&gt;distance (A)&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;fmax 0.05-&amp;gt;0.01&lt;/td&gt;
      &lt;td&gt;-0.002578&lt;/td&gt;
      &lt;td&gt;0.031803&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;fmax 0.01-&amp;gt;0.001&lt;/td&gt;
      &lt;td&gt;-0.000340&lt;/td&gt;
      &lt;td&gt;0.014891&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;



### 固体の例


```python
fmax_list = [0.05, 0.01, 0.001]
```


```python
atoms = bulk("Pt") * (4, 4, 4)
atoms.rattle(stdev=0.1)
atoms.calc = calc_mol
images = []
for fmax in tqdm(fmax_list):
    with FIRE(atoms, logfile=None) as opt:
        opt.run(fmax=fmax)
        images.append(atoms.copy())
```


      0%|          | 0/3 [00:00&lt;?, ?it/s]


    [WARNING 2025-07-18 05:36:23,397] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:23,442] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:23,484] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:23,560] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:23,628] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:23,704] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:23,782] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:23,874] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:23,948] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,059] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,135] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,192] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,257] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,330] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,448] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,520] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,601] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,716] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,783] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,869] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:24,924] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,026] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,096] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,144] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,212] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,261] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,333] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,390] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,444] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,487] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,562] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,621] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,696] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,746] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,802] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,885] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:25,964] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,050] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,126] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,211] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,279] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,332] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,383] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,469] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,523] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,571] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,624] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,692] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,767] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:26,984] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,044] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,106] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,166] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,210] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,297] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,379] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,486] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,571] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,639] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,740] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,818] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,920] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:27,982] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,056] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,118] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,167] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,235] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,317] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,440] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,528] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,602] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,662] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,723] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,781] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,878] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:28,962] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:29,019] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:29,101] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:29,174] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:29,267] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.



```python
for atoms in images:
    atoms.calc = calc_mol

results = {
    "fmax": [],
    "energy (eV)": [],
    "distance (A)": [],
}
for (fmax_old, old), (fmax_new, new) in windowed(zip(fmax_list, images), 2):
    de = new.get_potential_energy() - old.get_potential_energy()
    dx = max_distance(new, old)
    results["energy (eV)"].append(de)
    results["distance (A)"].append(dx)
    results["fmax"].append(f"fmax {fmax_old}-&gt;{fmax_new}")
pd.DataFrame(results)
```

    [WARNING 2025-07-18 05:36:29,328] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:29,403] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:29,468] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:36:29,598] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.





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
      &lt;th&gt;fmax&lt;/th&gt;
      &lt;th&gt;energy (eV)&lt;/th&gt;
      &lt;th&gt;distance (A)&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;fmax 0.05-&amp;gt;0.01&lt;/td&gt;
      &lt;td&gt;-0.005576&lt;/td&gt;
      &lt;td&gt;0.018937&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;fmax 0.01-&amp;gt;0.001&lt;/td&gt;
      &lt;td&gt;-0.000062&lt;/td&gt;
      &lt;td&gt;0.001370&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;



いずれの場合も、最も動いた原子でも0.01A程度しか動かず、エネルギーも0.01eV程度しか変わりません。
0.01eVというと、300Kで存在比が1.5倍になる程度のエネルギー差です。常温で扱う材料を探索したい場合、この値は数値計算で欲しい精度を超えています。
また、実験化学の精度に対応するとされるchemical accuracyは0.04eV程度ですが、多くの数値計算手法はその精度に達していないのが現状であり、その意味でもこの精度は必要ないことが多いでしょう。

0.01eVのエネルギー差がある2つの構造の300Kにおける存在比はボルツマン因子から以下のように求めます。


```python
np.exp(0.01 / (300 * kB))
```




    1.4722876299701388



## fmax=0.05eVで不十分な例

### エネルギーに関して精度不足な例

少し意地悪な例ですが、2,3-Dimethyl-2-buteneがPt表面にゆるやかに吸着しているような例を考えてみましょう。
このような例ではゆるやかな力が働くため、fmax=0.05eV/Aでは構造最適化が途中で止まってしまい、さらに最適化した場合と比べて大きなエネルギー差が生じてしまいます。
このような状況が発生する可能性があるならば、なるべく小さなfmaxを使ったほうがより信頼のおける構造最適化が出来るでしょう。


```python
bulk1111 = bulk("Pt")
bulk1111.calc = calc_d3
with FIRELBFGS(ExpCellFilter(bulk1111), logfile=None) as opt:
    opt.run(0.0001)
```

    /tmp/ipykernel_11479/3077852052.py:3: FutureWarning: Import ExpCellFilter from ase.filters
      with FIRELBFGS(ExpCellFilter(bulk1111), logfile=None) as opt:



```python
atoms = surface(bulk1111, (1, 1, 1), 4, vacuum=20.0) * (4, 4, 1)
c = atoms.cell[2, 2] / 2
atoms.constraints = [FixAtoms(mask=[atom.position[2] &lt; c for atom in atoms])]
atoms.rattle(stdev=0.1)
ads = smiles_to_atoms("CC(=C(C)C)C")
v = (ads.positions[0] - ads.positions[1]) / 2
p = (atoms.cell[0, :2] + atoms.cell[1, :2]) / 2
add_adsorbate(atoms, ads, 4.0, position=tuple(p))
atoms.positions[64:] = atoms.positions[64:] + v
atoms.calc = calc_d3
```


```python
view_ngl(atoms, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pt', 'C', 'H'), value…




```python
fmax_list = [0.05, 0.01, 0.001]
```


```python
images = []
for fmax in tqdm(fmax_list):
    with FIRELBFGS(atoms, logfile=None) as opt:
        opt.run(fmax=fmax)
        images.append(atoms.copy())
```


      0%|          | 0/3 [00:00&lt;?, ?it/s]



```python
for atoms in images:
    atoms.calc = calc_mol

results = {
    "fmax": [],
    "energy (eV)": [],
    "distance (A)": [],
}
for (fmax_old, old), (fmax_new, new) in windowed(zip(fmax_list, images), 2):
    de = new.get_potential_energy() - old.get_potential_energy()
    dx = max_distance(new, old)
    results["energy (eV)"].append(de)
    results["distance (A)"].append(dx)
    results["fmax"].append(f"fmax {fmax_old}-&gt;{fmax_new}")
pd.DataFrame(results)
```

    [WARNING 2025-07-18 05:37:04,739] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:37:04,793] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:37:04,922] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:37:04,982] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.





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
      &lt;th&gt;fmax&lt;/th&gt;
      &lt;th&gt;energy (eV)&lt;/th&gt;
      &lt;th&gt;distance (A)&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;fmax 0.05-&amp;gt;0.01&lt;/td&gt;
      &lt;td&gt;-0.279559&lt;/td&gt;
      &lt;td&gt;2.917367&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;fmax 0.01-&amp;gt;0.001&lt;/td&gt;
      &lt;td&gt;-0.004070&lt;/td&gt;
      &lt;td&gt;0.044530&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;




```python
view_ngl(images, representations=["ball+stick"], replace_structure=True)
```




    HBox(children=(NGLWidget(max_frame=2), VBox(children=(Dropdown(description='Show', options=('All', 'Pt', 'C', …



### 座標に関して精度不足な例（１）

Pt(111)面にエチレンが吸着したような例では、fmax=0.05eVでは中途半端な構造で構造最適化が止まってしまいます。


```python
bulk1111 = bulk("Pt")
bulk1111.calc = calc
with FIRELBFGS(ExpCellFilter(bulk1111), logfile=None) as opt:
    opt.run(0.0001)
```

    /tmp/ipykernel_11479/2764059397.py:3: FutureWarning: Import ExpCellFilter from ase.filters
      with FIRELBFGS(ExpCellFilter(bulk1111), logfile=None) as opt:



```python
atoms = surface(bulk1111, (1, 1, 1), 4, vacuum=20.0) * (4, 4, 1)
c = atoms.cell[2, 2] / 2
atoms.constraints = [FixAtoms(mask=[atom.position[2] &lt; c for atom in atoms])]
atoms.calc = calc
atoms.rattle(stdev=0.1)
add_adsorbate(atoms, smiles_to_atoms("C=C"), 2.0, position=tuple((atoms.cell[0, :2] + atoms.cell[1, :2]) / 2))
```


```python
view_ngl(atoms, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pt', 'C', 'H'), value…




```python
images = []
for fmax in tqdm(fmax_list):
    with FIRELBFGS(atoms, logfile=None) as opt:
        opt.run(fmax=fmax)
        images.append(atoms.copy())
```


      0%|          | 0/3 [00:00&lt;?, ?it/s]



```python
for atoms in images:
    atoms.calc = calc_mol

results = {
    "fmax": [],
    "energy (eV)": [],
    "distance (A)": [],
}
for (fmax_old, old), (fmax_new, new) in windowed(zip(fmax_list, images), 2):
    de = new.get_potential_energy() - old.get_potential_energy()
    dx = max_distance(new, old)
    results["energy (eV)"].append(de)
    results["distance (A)"].append(dx)
    results["fmax"].append(f"fmax {fmax_old}-&gt;{fmax_new}")
pd.DataFrame(results)
```

    [WARNING 2025-07-18 05:37:21,743] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:37:21,810] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:37:21,924] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.
    [WARNING 2025-07-18 05:37:21,993] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.





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
      &lt;th&gt;fmax&lt;/th&gt;
      &lt;th&gt;energy (eV)&lt;/th&gt;
      &lt;th&gt;distance (A)&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;fmax 0.05-&amp;gt;0.01&lt;/td&gt;
      &lt;td&gt;-0.054356&lt;/td&gt;
      &lt;td&gt;0.434370&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;fmax 0.01-&amp;gt;0.001&lt;/td&gt;
      &lt;td&gt;-0.000365&lt;/td&gt;
      &lt;td&gt;0.016433&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;




```python
view_ngl(images, representations=["ball+stick"], replace_structure=True)
```




    HBox(children=(NGLWidget(max_frame=2), VBox(children=(Dropdown(description='Show', options=('All', 'Pt', 'C', …



中途半端で止まってしまった構造と最適化しきった構造を見比べると違った印象をうけると思います。
とはいえこの例ではエネルギー差は大きくないため、場合によってはこの違いを無視するという判断を行っても良いかもしれません。

### 座標について精度不足な例（２）

構造最適化が難しい系として良く知られているトルエンのメチル基の回転を試してみます。


```python
atoms = smiles_to_atoms("Cc1ccccc1")
tmp = atoms[7:10]
tmp.rotate([1.0, 0.0, 0.0], 15.0)
atoms.positions[7:10] = tmp.positions
atoms.calc = calc_mol
```


```python
view_ngl(atoms, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'C', 'H'), value='All'…




```python
images = []
for fmax in tqdm(fmax_list):
    with FIRELBFGS(atoms, logfile=None) as opt:
        opt.run(fmax=fmax)
        images.append(atoms.copy())
```


      0%|          | 0/3 [00:00&lt;?, ?it/s]



```python
for atoms in images:
    atoms.calc = calc_mol

results = {
    "fmax": [],
    "energy (eV)": [],
    "distance (A)": [],
}
for (fmax_old, old), (fmax_new, new) in windowed(zip(fmax_list, images), 2):
    de = new.get_potential_energy() - old.get_potential_energy()
    dx = max_distance(new, old)
    results["energy (eV)"].append(de)
    results["distance (A)"].append(dx)
    results["fmax"].append(f"fmax {fmax_old}-&gt;{fmax_new}")
pd.DataFrame(results)
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
      &lt;th&gt;fmax&lt;/th&gt;
      &lt;th&gt;energy (eV)&lt;/th&gt;
      &lt;th&gt;distance (A)&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;fmax 0.05-&amp;gt;0.01&lt;/td&gt;
      &lt;td&gt;-0.001871&lt;/td&gt;
      &lt;td&gt;0.030084&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;fmax 0.01-&amp;gt;0.001&lt;/td&gt;
      &lt;td&gt;-0.000052&lt;/td&gt;
      &lt;td&gt;0.004847&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;/div&gt;




```python
view_ngl(images, representations=["ball+stick"], replace_structure=True)
```




    HBox(children=(NGLWidget(max_frame=2), VBox(children=(Dropdown(description='Show', options=('All', 'C', 'H'), …



トルエンのメチル基の回転は角度が変わってもエネルギー差があまり変わらないことが判ると思います。
このような例では小さなfmaxを用いないと最安定構造の予測が出来ません。
0.002eVのエネルギー差は300Kで存在比が1.08倍になる程度の影響しかないので、ほとんど影響がないと言っても構わないでしょう。
従ってこのような小さなエネルギー差は場合によっては無視してかまいません。
ただし、振動解析のようなエネルギーが極小であることを前提とした解析を行いたい場合は時として解析に悪影響を及ぼすかもしれません。


```python
np.exp(0.002 / (300 * kB))
```




    1.080434721876578


