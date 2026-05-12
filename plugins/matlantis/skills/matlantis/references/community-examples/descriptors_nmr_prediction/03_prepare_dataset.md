Copyright Matlantis Corp. as contributors to Matlantis contrib project


```python
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import pickle
import h5py
import numpy as np
import plotly.express as px
from sklearn.model_selection import train_test_split

from tqdm import tqdm
from src.dataset import HDF5Dataset

%load_ext autoreload
%autoreload 2
```

## 3. Prepare 1H dataset
dataset for 1H nmr prediction was prepared.  
The outlier data was removed. (23 data from 1208486 data, which magneticshielding value is under 15)  


```python
#input_file = "dataset/qm9_vapor_pbev8.h5"
input_file = "dataset/qm9_vapor_r2scan.h5"
```


```python
X = []
y = []
skipped = 0
skipped_other_element = 0

print("Preparing datasets...")
with h5py.File(input_file, "r") as f:
    m_k = list(f.keys())
    for k_id in tqdm(m_k):
        for a_k in f[k_id].keys():
            atoms = f[k_id][a_k]

            if atoms["element"][()] == b"H":
                if atoms["value"][()] &lt; 15:  # remove here and use not filtered data
                    skipped += 1
                    continue
                X.append(atoms["scalar"][()])
                y.append([k_id, atoms["value"][()]])
            else:
                skipped_other_element += 1
                continue

print(f"{skipped} data were skipped in {len(X) + skipped} data")
```

    Preparing datasets...


    100%|██████████| 130831/130831 [30:00&lt;00:00, 72.68it/s] 


    23 data were skipped in 1208486 data



```python
df_x = pd.DataFrame(X, columns=[f"desc-{x}" for x in range(256)])
df_y = pd.DataFrame(y, columns=["molecule_idx", "value"])

df_all = pd.concat([df_x, df_y], axis=1)
del df_x
del df_y
df_all.to_csv(input_file[:-3]+"_H.csv",index=False)
```


```python
df_all.shape
```




    (1208463, 258)




```python
df_all.head()
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
      &lt;th&gt;desc-0&lt;/th&gt;
      &lt;th&gt;desc-1&lt;/th&gt;
      &lt;th&gt;desc-2&lt;/th&gt;
      &lt;th&gt;desc-3&lt;/th&gt;
      &lt;th&gt;desc-4&lt;/th&gt;
      &lt;th&gt;desc-5&lt;/th&gt;
      &lt;th&gt;desc-6&lt;/th&gt;
      &lt;th&gt;desc-7&lt;/th&gt;
      &lt;th&gt;desc-8&lt;/th&gt;
      &lt;th&gt;desc-9&lt;/th&gt;
      &lt;th&gt;...&lt;/th&gt;
      &lt;th&gt;desc-248&lt;/th&gt;
      &lt;th&gt;desc-249&lt;/th&gt;
      &lt;th&gt;desc-250&lt;/th&gt;
      &lt;th&gt;desc-251&lt;/th&gt;
      &lt;th&gt;desc-252&lt;/th&gt;
      &lt;th&gt;desc-253&lt;/th&gt;
      &lt;th&gt;desc-254&lt;/th&gt;
      &lt;th&gt;desc-255&lt;/th&gt;
      &lt;th&gt;molecule_idx&lt;/th&gt;
      &lt;th&gt;value&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;0.725259&lt;/td&gt;
      &lt;td&gt;0.257291&lt;/td&gt;
      &lt;td&gt;0.238289&lt;/td&gt;
      &lt;td&gt;-0.157121&lt;/td&gt;
      &lt;td&gt;0.562056&lt;/td&gt;
      &lt;td&gt;0.041166&lt;/td&gt;
      &lt;td&gt;1.364311&lt;/td&gt;
      &lt;td&gt;0.292885&lt;/td&gt;
      &lt;td&gt;0.311306&lt;/td&gt;
      &lt;td&gt;0.320144&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.090678&lt;/td&gt;
      &lt;td&gt;-0.079473&lt;/td&gt;
      &lt;td&gt;0.432806&lt;/td&gt;
      &lt;td&gt;1.070554&lt;/td&gt;
      &lt;td&gt;0.012076&lt;/td&gt;
      &lt;td&gt;-0.087990&lt;/td&gt;
      &lt;td&gt;0.165963&lt;/td&gt;
      &lt;td&gt;-0.098839&lt;/td&gt;
      &lt;td&gt;m000001&lt;/td&gt;
      &lt;td&gt;31.5743&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;0.725260&lt;/td&gt;
      &lt;td&gt;0.257290&lt;/td&gt;
      &lt;td&gt;0.238290&lt;/td&gt;
      &lt;td&gt;-0.157120&lt;/td&gt;
      &lt;td&gt;0.562059&lt;/td&gt;
      &lt;td&gt;0.041163&lt;/td&gt;
      &lt;td&gt;1.364309&lt;/td&gt;
      &lt;td&gt;0.292885&lt;/td&gt;
      &lt;td&gt;0.311308&lt;/td&gt;
      &lt;td&gt;0.320145&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.090679&lt;/td&gt;
      &lt;td&gt;-0.079473&lt;/td&gt;
      &lt;td&gt;0.432805&lt;/td&gt;
      &lt;td&gt;1.070554&lt;/td&gt;
      &lt;td&gt;0.012075&lt;/td&gt;
      &lt;td&gt;-0.087990&lt;/td&gt;
      &lt;td&gt;0.165964&lt;/td&gt;
      &lt;td&gt;-0.098838&lt;/td&gt;
      &lt;td&gt;m000001&lt;/td&gt;
      &lt;td&gt;31.5743&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;0.725275&lt;/td&gt;
      &lt;td&gt;0.257281&lt;/td&gt;
      &lt;td&gt;0.238300&lt;/td&gt;
      &lt;td&gt;-0.157112&lt;/td&gt;
      &lt;td&gt;0.562067&lt;/td&gt;
      &lt;td&gt;0.041146&lt;/td&gt;
      &lt;td&gt;1.364298&lt;/td&gt;
      &lt;td&gt;0.292880&lt;/td&gt;
      &lt;td&gt;0.311308&lt;/td&gt;
      &lt;td&gt;0.320150&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.090683&lt;/td&gt;
      &lt;td&gt;-0.079477&lt;/td&gt;
      &lt;td&gt;0.432794&lt;/td&gt;
      &lt;td&gt;1.070552&lt;/td&gt;
      &lt;td&gt;0.012067&lt;/td&gt;
      &lt;td&gt;-0.087986&lt;/td&gt;
      &lt;td&gt;0.165965&lt;/td&gt;
      &lt;td&gt;-0.098826&lt;/td&gt;
      &lt;td&gt;m000001&lt;/td&gt;
      &lt;td&gt;31.5744&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;0.725273&lt;/td&gt;
      &lt;td&gt;0.257281&lt;/td&gt;
      &lt;td&gt;0.238300&lt;/td&gt;
      &lt;td&gt;-0.157114&lt;/td&gt;
      &lt;td&gt;0.562064&lt;/td&gt;
      &lt;td&gt;0.041147&lt;/td&gt;
      &lt;td&gt;1.364298&lt;/td&gt;
      &lt;td&gt;0.292879&lt;/td&gt;
      &lt;td&gt;0.311307&lt;/td&gt;
      &lt;td&gt;0.320149&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.090682&lt;/td&gt;
      &lt;td&gt;-0.079476&lt;/td&gt;
      &lt;td&gt;0.432795&lt;/td&gt;
      &lt;td&gt;1.070552&lt;/td&gt;
      &lt;td&gt;0.012067&lt;/td&gt;
      &lt;td&gt;-0.087986&lt;/td&gt;
      &lt;td&gt;0.165965&lt;/td&gt;
      &lt;td&gt;-0.098829&lt;/td&gt;
      &lt;td&gt;m000001&lt;/td&gt;
      &lt;td&gt;31.5744&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;0.702751&lt;/td&gt;
      &lt;td&gt;0.412856&lt;/td&gt;
      &lt;td&gt;0.481192&lt;/td&gt;
      &lt;td&gt;0.238448&lt;/td&gt;
      &lt;td&gt;0.789768&lt;/td&gt;
      &lt;td&gt;0.111927&lt;/td&gt;
      &lt;td&gt;1.303263&lt;/td&gt;
      &lt;td&gt;-0.103865&lt;/td&gt;
      &lt;td&gt;-0.169121&lt;/td&gt;
      &lt;td&gt;0.390848&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;-0.118581&lt;/td&gt;
      &lt;td&gt;-0.326023&lt;/td&gt;
      &lt;td&gt;0.321385&lt;/td&gt;
      &lt;td&gt;0.600912&lt;/td&gt;
      &lt;td&gt;0.007953&lt;/td&gt;
      &lt;td&gt;-0.082237&lt;/td&gt;
      &lt;td&gt;0.067479&lt;/td&gt;
      &lt;td&gt;0.129453&lt;/td&gt;
      &lt;td&gt;m000002&lt;/td&gt;
      &lt;td&gt;31.9746&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;p&gt;5 rows × 258 columns&lt;/p&gt;
&lt;/div&gt;


