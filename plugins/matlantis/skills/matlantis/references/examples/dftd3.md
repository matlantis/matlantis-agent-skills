# DFT-D3 法を使った分散力補正

DFT-D3法を用いて分散力を補正する例を示します。

## 準備
DFT-D3法を使うには、`pfp-api-client&gt;= 1.1.0`をインストールする必要があります。


```python
# If you have already installed pfp-api-client&gt;=1.1.0, you can skip this.
! pip install 'pfp-api-client&gt;=1.1.0'
```




```python
import pfp_api_client
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator()
print(f"pfp_api_client: {pfp_api_client.__version__}")
print(f"current pfp model version: {estimator.model_version}")
```

    pfp_api_client: 2.0.1
    current pfp model version: v8.0.0



```python
import numpy as np

import ase
from ase import Atoms, units
from ase.visualize import view
from ase.optimize import BFGS

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

estimator_d3 = Estimator(calc_mode=EstimatorCalcMode.PBE_U_PLUS_D3)
calculator_d3 = ASECalculator(estimator_d3)

estimator = Estimator()
calculator = ASECalculator(estimator)
```


## 構造の準備
今回は、分散力補正の有無の違いをみるためベンゼン二量体を使用します。


```python
from ase.build import molecule
benzene1, benzene2 = molecule("C6H6"), molecule("C6H6")
benzene2.positions[:, 2] += 3.5
benzene2.positions[:, 1] += 1.395248
atoms = benzene1 + benzene2
```


```python
from ase.visualize import view
v = view(atoms, viewer='ngl')
display(v)
```






    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'H', 'C'), value='All'…



## ポテンシャルの計算
* 分散力補正**なし**のポテンシャル


```python
atoms1 = atoms.copy()
atoms1.calc = calculator
print(f"The potential energy w/o D3 is {atoms1.get_potential_energy()} eV")
```

    The potential energy w/o D3 is -122.1310962487184 eV



* 分散力補正**あり**のポテンシャル


```python
atoms2 = atoms.copy()
atoms2.calc = calculator_d3
print(f"The potential energy with D3 is {atoms2.get_potential_energy()} eV")
```

    The potential energy with D3 is -122.90645574608556 eV



## ふたつのベンゼン分子の層間距離の計算

* 分散力補正**なし**の層間距離


```python
from ase.optimize import BFGS

bfgs = BFGS(atoms1, logfile=None)
bfgs.run(fmax=0.01)

print(f"The interlayer distance w/o D3 is {atoms1.positions[12,2]-atoms1.positions[0,2]} Angstrom.")
```

    The interlayer distance w/o D3 is 4.648246688926918 Angstrom.



* 分散力補正**あり**の層間距離


```python
bfgs = BFGS(atoms2, logfile=None)
bfgs.run(fmax=0.01)

print(f"The interlayer distance with D3 is {atoms2.positions[12,2]-atoms2.positions[0,2]} Angstrom.")
```

    The interlayer distance with D3 is 3.508276358959358 Angstrom.



## まとめ
DFT-D3法により分散力補正を計算することができました。ベンゼン二量体においては、DFT-D3法は、
a) ポテンシャルを深くし、
b) 層間距離を小さくしました。DF-SCS-LMP2法を用いると層間距離は3.54 Å と計算されます(DF-SCS-LMP2). [1])

[1] Phys. Chem. Chem. Phys., 2006, 8, 4072–4078
