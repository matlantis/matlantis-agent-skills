Copyright Matlantis Corp. as contributors to Matlantis contrib project


```python
import os
import glob
import pandas as pd
import re
```


```python
from ase.io import read, write
from pfcc_extras import view_ngl
```


    


## 1. Data extraction from raw data
The raw data was obtained from [QM9NMR dataset](https://moldis-group.github.io/qm9nmr/).


```python
base_path = "."
path_list = glob.glob(f"{base_path}/raw_data/71Zaisc1SPORF1wnneyUtg/*/*") + glob.glob(f"{base_path}/raw_data/gHm0Zu65Q--8x91ap90IMg/*/*") + \
            glob.glob(f"{base_path}/raw_data/ry7E8IsyRcGlvOyWZkxfTQ/*/*") + glob.glob(f"{base_path}/raw_data/V1FvXN_2ScasSB36-oDzeQ/*/*")
len(path_list)
```




    130831




```python
# Following scripts can not be executed in normal instance (6.5G Memory) with full size dataset
# Smaller dataset is better for testing and demonstration
# To make smaller dataset, please execute below 2 lines
# import random
# path_list = random.sample(path_list, 10000)
```


```python
os.makedirs("molecule_csv", exist_ok=True)
os.makedirs("molecule_xyz", exist_ok=True)
```

### Extract structure and value from the dataset


```python
for idx, path in enumerate(path_list):
    if idx % 1000 == 0:
        print(idx)
    
    molecule_id = "m" + path.split("_")[-1][:-4]

    # Extract value
    with open(path) as f:
        cnt = f.read()

    str_list = cnt.split("\n")
    extracted_list = [s for s in str_list if re.match('.*Anisotropy =.*', s)]
    # not Anisotoropy, but isotoropy was used with correction
    info_list = [[line[9], float(line[26:37].strip())] for line in extracted_list]

    df = pd.DataFrame(info_list, columns=["ELEMENT", "VALUE"]).reset_index().rename(columns={"index": "ATOM_IDX"})
    df["MOLECULE_IDX"] = molecule_id
    df.to_csv(f"molecule_csv/{molecule_id}.csv", index=False)
    
    #Extract strucutre
    atoms = read(path)
    atoms.calc = None  # To escape the error, calculator must be None
    write(f"molecule_xyz/{molecule_id}.xyz", atoms)

```

    0
    1000
    2000
    3000
    4000
    5000
    6000
    7000
    8000
    9000


### Confirm the extracted data


```python
atoms = read(f"molecule_xyz/{molecule_id}.xyz")
print(atoms, len(atoms))
view_ngl(atoms, representations=["ball+stick"])
```

    Atoms(symbols='C8OH14', pbc=False) 23





    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H', 'C'), value=…


