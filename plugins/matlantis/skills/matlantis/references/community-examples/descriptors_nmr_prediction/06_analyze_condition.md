Copyright Matlantis Corp. as contributors to Matlantis contrib project


```python
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
```

## 6. Analyze the condition screening results
The condition screening results were confirmed.  
As a result the condition was choosed as follows:  
- Descriptors: pbev8 / r2scan
- Preprocess: MinMaxScaling + PCA with whiting
- model: lightGBM


```python
with open("sample_output/05_res_data.txt", "r") as f:
    cnt = f.read().split("\n")
```


```python
## temporary functions to extract info from log txt
def split_info(model_str):
    res = model_str.split("_")
    res[-2] = int(res[-2][1:])
    res[-1] = int(res[-1][1:])
    res[-1] = 256 if res[-1] &lt; 0 else res[-1]
    return res

def extract_value(acc_str):
    return float(acc_str.split(":")[-1])

def summarize_info(cnt):
    x = cnt.split(" ")
    return split_info(x[0]) + [extract_value(x[1])] + [float(x[3][:-1])] + [extract_value(x[4])] + [float(x[6])]

info_list = [summarize_info(x) for x in cnt[:-1]]
```


```python
df = pd.DataFrame(info_list, columns=["desc", "scaler", "white", "model", "data_size", "n_feature", "train_mae", "train_std", "val_mae", "val_std"])
df.head()
df["log_dsize"] = [np.log10(x) for x in df["data_size"]]
```


```python
df.sort_values("val_mae")
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
      &lt;th&gt;desc&lt;/th&gt;
      &lt;th&gt;scaler&lt;/th&gt;
      &lt;th&gt;white&lt;/th&gt;
      &lt;th&gt;model&lt;/th&gt;
      &lt;th&gt;data_size&lt;/th&gt;
      &lt;th&gt;n_feature&lt;/th&gt;
      &lt;th&gt;train_mae&lt;/th&gt;
      &lt;th&gt;train_std&lt;/th&gt;
      &lt;th&gt;val_mae&lt;/th&gt;
      &lt;th&gt;val_std&lt;/th&gt;
      &lt;th&gt;log_dsize&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;154&lt;/th&gt;
      &lt;td&gt;pbev8&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202908&lt;/td&gt;
      &lt;td&gt;256&lt;/td&gt;
      &lt;td&gt;0.079&lt;/td&gt;
      &lt;td&gt;0.005&lt;/td&gt;
      &lt;td&gt;0.124&lt;/td&gt;
      &lt;td&gt;0.002&lt;/td&gt;
      &lt;td&gt;6.080232&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;466&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202914&lt;/td&gt;
      &lt;td&gt;256&lt;/td&gt;
      &lt;td&gt;0.081&lt;/td&gt;
      &lt;td&gt;0.004&lt;/td&gt;
      &lt;td&gt;0.125&lt;/td&gt;
      &lt;td&gt;0.002&lt;/td&gt;
      &lt;td&gt;6.080235&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;310&lt;/th&gt;
      &lt;td&gt;pbev8&lt;/td&gt;
      &lt;td&gt;STN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202863&lt;/td&gt;
      &lt;td&gt;256&lt;/td&gt;
      &lt;td&gt;0.085&lt;/td&gt;
      &lt;td&gt;0.006&lt;/td&gt;
      &lt;td&gt;0.128&lt;/td&gt;
      &lt;td&gt;0.003&lt;/td&gt;
      &lt;td&gt;6.080216&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;620&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;STN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202829&lt;/td&gt;
      &lt;td&gt;128&lt;/td&gt;
      &lt;td&gt;0.088&lt;/td&gt;
      &lt;td&gt;0.004&lt;/td&gt;
      &lt;td&gt;0.131&lt;/td&gt;
      &lt;td&gt;0.002&lt;/td&gt;
      &lt;td&gt;6.080204&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;152&lt;/th&gt;
      &lt;td&gt;pbev8&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202908&lt;/td&gt;
      &lt;td&gt;128&lt;/td&gt;
      &lt;td&gt;0.094&lt;/td&gt;
      &lt;td&gt;0.010&lt;/td&gt;
      &lt;td&gt;0.132&lt;/td&gt;
      &lt;td&gt;0.005&lt;/td&gt;
      &lt;td&gt;6.080232&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;...&lt;/th&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;479&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;STN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;xgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;256&lt;/td&gt;
      &lt;td&gt;0.191&lt;/td&gt;
      &lt;td&gt;0.129&lt;/td&gt;
      &lt;td&gt;1.217&lt;/td&gt;
      &lt;td&gt;0.191&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;164&lt;/th&gt;
      &lt;td&gt;pbev8&lt;/td&gt;
      &lt;td&gt;STN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;128&lt;/td&gt;
      &lt;td&gt;0.241&lt;/td&gt;
      &lt;td&gt;0.152&lt;/td&gt;
      &lt;td&gt;1.243&lt;/td&gt;
      &lt;td&gt;0.299&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;477&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;STN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;xgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;128&lt;/td&gt;
      &lt;td&gt;0.064&lt;/td&gt;
      &lt;td&gt;0.035&lt;/td&gt;
      &lt;td&gt;1.263&lt;/td&gt;
      &lt;td&gt;0.217&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;478&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;STN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;256&lt;/td&gt;
      &lt;td&gt;0.427&lt;/td&gt;
      &lt;td&gt;0.069&lt;/td&gt;
      &lt;td&gt;1.293&lt;/td&gt;
      &lt;td&gt;0.290&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;166&lt;/th&gt;
      &lt;td&gt;pbev8&lt;/td&gt;
      &lt;td&gt;STN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;256&lt;/td&gt;
      &lt;td&gt;0.332&lt;/td&gt;
      &lt;td&gt;0.074&lt;/td&gt;
      &lt;td&gt;1.315&lt;/td&gt;
      &lt;td&gt;0.359&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;p&gt;623 rows × 11 columns&lt;/p&gt;
&lt;/div&gt;



## Appendix
The effect of datasize and num_features were confirmed with lightGBM and XGBoost models.



```python
## data_size
cond_list = [
    ["MMN+PCA", "r2scan", "xgb", "256"],
    ["MMN+PCA", "r2scan", "lgb", "256"],
]

color_list = ["red", "blue"]

for idx, cond in enumerate(cond_list):
    df_tmp = df.query(f"scaler == '{cond[0]}' &amp; desc == '{cond[1]}' &amp; model == '{cond[2]}' &amp; n_feature == {cond[3]}")
    plt.plot(df_tmp["log_dsize"], df_tmp["train_mae"], label="_".join(cond) + "_train", 
             color=color_list[idx], marker="o", linestyle=":")
    plt.plot(df_tmp["log_dsize"], df_tmp["val_mae"], label="_".join(cond) + "_val",
             color=color_list[idx], marker="x")
plt.xlabel("log10(datasize)")
plt.ylabel("MAE / ppm")
plt.title("Effect of datasize")
plt.legend()
```




    &lt;matplotlib.legend.Legend at 0x7f01b2a7a7d0&gt;




    
![png](output_8_1.png)
    



```python
df.query(f"scaler == '{cond[0]}' &amp; desc == '{cond[1]}' &amp; model == '{cond[2]}' ")
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
      &lt;th&gt;desc&lt;/th&gt;
      &lt;th&gt;scaler&lt;/th&gt;
      &lt;th&gt;white&lt;/th&gt;
      &lt;th&gt;model&lt;/th&gt;
      &lt;th&gt;data_size&lt;/th&gt;
      &lt;th&gt;n_feature&lt;/th&gt;
      &lt;th&gt;train_mae&lt;/th&gt;
      &lt;th&gt;train_std&lt;/th&gt;
      &lt;th&gt;val_mae&lt;/th&gt;
      &lt;th&gt;val_std&lt;/th&gt;
      &lt;th&gt;log_dsize&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;312&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;8&lt;/td&gt;
      &lt;td&gt;0.515&lt;/td&gt;
      &lt;td&gt;0.110&lt;/td&gt;
      &lt;td&gt;0.887&lt;/td&gt;
      &lt;td&gt;0.150&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;314&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;16&lt;/td&gt;
      &lt;td&gt;0.381&lt;/td&gt;
      &lt;td&gt;0.078&lt;/td&gt;
      &lt;td&gt;0.786&lt;/td&gt;
      &lt;td&gt;0.080&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;316&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;32&lt;/td&gt;
      &lt;td&gt;0.381&lt;/td&gt;
      &lt;td&gt;0.128&lt;/td&gt;
      &lt;td&gt;0.911&lt;/td&gt;
      &lt;td&gt;0.171&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;318&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;64&lt;/td&gt;
      &lt;td&gt;0.489&lt;/td&gt;
      &lt;td&gt;0.044&lt;/td&gt;
      &lt;td&gt;0.998&lt;/td&gt;
      &lt;td&gt;0.203&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;320&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;100&lt;/td&gt;
      &lt;td&gt;128&lt;/td&gt;
      &lt;td&gt;0.427&lt;/td&gt;
      &lt;td&gt;0.043&lt;/td&gt;
      &lt;td&gt;0.997&lt;/td&gt;
      &lt;td&gt;0.181&lt;/td&gt;
      &lt;td&gt;2.000000&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;...&lt;/th&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;458&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202914&lt;/td&gt;
      &lt;td&gt;16&lt;/td&gt;
      &lt;td&gt;0.182&lt;/td&gt;
      &lt;td&gt;0.004&lt;/td&gt;
      &lt;td&gt;0.216&lt;/td&gt;
      &lt;td&gt;0.002&lt;/td&gt;
      &lt;td&gt;6.080235&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;460&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202914&lt;/td&gt;
      &lt;td&gt;32&lt;/td&gt;
      &lt;td&gt;0.123&lt;/td&gt;
      &lt;td&gt;0.003&lt;/td&gt;
      &lt;td&gt;0.165&lt;/td&gt;
      &lt;td&gt;0.001&lt;/td&gt;
      &lt;td&gt;6.080235&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;462&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202914&lt;/td&gt;
      &lt;td&gt;64&lt;/td&gt;
      &lt;td&gt;0.111&lt;/td&gt;
      &lt;td&gt;0.013&lt;/td&gt;
      &lt;td&gt;0.147&lt;/td&gt;
      &lt;td&gt;0.007&lt;/td&gt;
      &lt;td&gt;6.080235&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;464&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202914&lt;/td&gt;
      &lt;td&gt;128&lt;/td&gt;
      &lt;td&gt;0.092&lt;/td&gt;
      &lt;td&gt;0.006&lt;/td&gt;
      &lt;td&gt;0.132&lt;/td&gt;
      &lt;td&gt;0.003&lt;/td&gt;
      &lt;td&gt;6.080235&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;466&lt;/th&gt;
      &lt;td&gt;r2scan&lt;/td&gt;
      &lt;td&gt;MMN+PCA&lt;/td&gt;
      &lt;td&gt;white&lt;/td&gt;
      &lt;td&gt;lgb&lt;/td&gt;
      &lt;td&gt;1202914&lt;/td&gt;
      &lt;td&gt;256&lt;/td&gt;
      &lt;td&gt;0.081&lt;/td&gt;
      &lt;td&gt;0.004&lt;/td&gt;
      &lt;td&gt;0.125&lt;/td&gt;
      &lt;td&gt;0.002&lt;/td&gt;
      &lt;td&gt;6.080235&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;p&gt;78 rows × 11 columns&lt;/p&gt;
&lt;/div&gt;




```python
## n_feature
cond_list = [
    ["MMN+PCA", "r2scan", "xgb", "1202914"],
    ["MMN+PCA", "r2scan", "lgb", "1202914"],
]

color_list = ["red", "blue"]

for idx, cond in enumerate(cond_list):
    df_tmp = df.query(f"scaler == '{cond[0]}' &amp; desc == '{cond[1]}' &amp; model == '{cond[2]}' &amp; data_size == {cond[3]}")
    plt.plot(df_tmp["n_feature"], df_tmp["train_mae"], label="_".join(cond[:3]) + "_train", 
             color=color_list[idx], marker="o", linestyle=":")
    plt.plot(df_tmp["n_feature"], df_tmp["val_mae"], label="_".join(cond[:3]) + "_val", 
             color=color_list[idx], marker="x")

plt.xlabel("num_features")
plt.ylabel("MAE / ppm")
plt.title("Effect of num_features")

plt.legend()
```




    &lt;matplotlib.legend.Legend at 0x7f01b294a710&gt;




    
![png](output_10_1.png)
    

