Copyright Matlantis Corp. as contributors to Matlantis contrib project


```python
import glob
import os
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

from ase import Atoms
from ase.io import read, write
from ase.visualize import view
from ase.data import atomic_numbers
from ase.optimize import LBFGS, FIRE
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from src.dataset import HDF5Dataset
from pfcc_extras import view_ngl

%load_ext autoreload
%autoreload 2
```


```python
# for PFP
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

calc_mode = "r2scan"
model_version = "v8.0.0"
estimator = Estimator(calc_mode=calc_mode, model_version=model_version)
calculator = ASECalculator(estimator)
```

## 2. Get PFP descriptors
From the molecule structure files, PFP descriptors were calculated.  
It takes around 2.5 h for 130,000 data.


```python
path_list = glob.glob(os.path.join("../qm9/molecule_xyz", "*.xyz"))
len(path_list)
```


```python
output_dir = f"descriptors_{calc_mode}" 
os.makedirs(output_dir, exist_ok=True)
```


```python
conv_dict = {v:k for k, v in atomic_numbers.items()}

def get_df_descript(calculator, mol_id, atoms):
    descriptors = calculator.get_descriptors(atoms)
    df_desc = pd.DataFrame(descriptors["scalar"], columns=[f"desc-{x}" for x in range(descriptors["scalar"].shape[1])])
    df_desc["ELEMENT"] = [conv_dict[x] for x in atoms.get_atomic_numbers()]
    df_desc["MOLECULE_IDX"] = mol_id 
    df_desc["ATOM_IDX"] = df_desc.index
    return df_desc
```


```python
%%time
for idx, path in enumerate(path_list):
    if idx % 1000 == 0:
        print(idx)

    mol_id = path.split("/")[-1][:-4]

    atoms = read(path)
    atoms.calc = None
    df_desc = get_df_descript(calculator, mol_id, atoms)

    df_desc.to_csv(f"{output_dir}/{mol_id}.csv", index=False)

    del atoms
    del df_desc

```


```python
csv_list = glob.glob(f"{output_dir}/*.csv")
len(csv_list), csv_list[:3]
```


```python
m_keys = [x.split("/")[-1][:-4] for x in csv_list]
m_keys[:3]
```


```python
os.makedirs("dataset", exist_ok=True)
f_name = f"dataset/qm9_vapor_{calc_mode}.h5"

with h5py.File(f_name, "w") as f:
    for m_k in tqdm(m_keys):
        desc_path = f"{output_dir}/{m_k}.csv"
        value_path = f"../qm9/molecule_csv/{m_k}.csv"
        df_desc = pd.read_csv(desc_path, index_col=-1)
        df_v = pd.read_csv(value_path, index_col=0)

        if df_desc.shape[0] != df_v.shape[0]:  # At least check the atom number
            print(f"Something wrong {m_k}")
            continue
        
        m_g = f.create_group(m_k)
        
        for a_k in df_desc.index:
            if df_desc.at[a_k, "ELEMENT"] == df_v.at[a_k, "ELEMENT"]:
                a_k_str = str(a_k)
                a_g = m_g.create_group(a_k_str)
                a_g.create_dataset("value", data=df_v.at[a_k, "VALUE"])
                a_g.create_dataset("element", data=df_v.at[a_k, "ELEMENT"])
                a_g.create_dataset("scalar", data=df_desc.iloc[a_k, :-2].values.astype("float"))
            else:
                print(f"Something wrong {m_k}, {a_k}")

        del df_desc, df_v

```


```python
with h5py.File(f_name, "r") as f:
    print(len(f.keys()))
```
