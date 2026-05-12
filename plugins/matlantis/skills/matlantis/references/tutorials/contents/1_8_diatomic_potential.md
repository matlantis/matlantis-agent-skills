# Diatomic Potential

2体ポテンシャルエネルギーとは、真空中に2つの原子だけを置いたときの相互作用のエネルギーです。原子間距離によって値が変わります。

## 水素分子の2体Potential energyを計算してみる

ここまで習ったことを利用して、水素分子の２体ポテンシャルエネルギーを計算してみましょう。

水素原子間が様々な距離になるように水素分子を生成しそのエネルギーを計算、プロットしてみます。


```python
from ase import Atoms
from ase.visualize import view

import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode


estimator = Estimator(model_version="v8.0.0")
calculator = ASECalculator(estimator)
```


```python
import numpy as np
from tqdm.auto import tqdm

distances = np.linspace(0.3, 6.5, 100)
energy_list = []
for d in tqdm(distances):
    atoms = Atoms("H2", [[0, 0, 0], [0, 0, d]])
    atoms.calc = calculator
    E_pot = atoms.get_potential_energy()
    energy_list.append(E_pot)

energies = np.array(energy_list)
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]



```python
import matplotlib.pyplot as plt


plt.plot(distances, energies)
# H2 bond length = 0.74A
plt.vlines(0.74, np.min(energies), np.max(energies), color="red")
plt.xlabel("bond distance: [A]")
plt.ylabel("Potential energy: [eV]")
plt.title("H2 potential curve")
plt.show()
```


    
![png](output_3_0.png)
    



```python
distances[np.argmin(energies)]
```




    0.7383838383838384



同じ水素原子2つからなる水素分子でも原子間距離が変わるだけでそのエネルギーは大きく異なることがわかります。&lt;br/&gt;
実験的には水素分子は結合長さ 0.74A で最安定であるということが知られていますが、実際にMatlantisで計算を行ってみると近い位置でエネルギーが一番低くなっていることが確認できました。

このように、様々な構造に対してエネルギーが計算できることで、自然界で物質がどのような構造を取っているのかを解析することができます。
