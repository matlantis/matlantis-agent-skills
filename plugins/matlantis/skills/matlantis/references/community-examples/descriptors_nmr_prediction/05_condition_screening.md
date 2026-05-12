Copyright Matlantis Corp. as contributors to Matlantis contrib project


```python
import random
import pickle
import glob
import pandas as pd
import numpy as np
import xgboost as xgb
import lightgbm as lgb
from matplotlib import pyplot as plt
from sklearn import svm
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.decomposition import PCA
%load_ext autoreload
%autoreload 2
```

## 5. Condition screening
The second screening was conducted to confirm the effect of descriptor type, preprocess, datasize and num_feature.  
The results are summarized in the next notebook.


```python
lgb_params = {
    "objective": "regression",
    "metric": "mae",
    "random_state": 12345,
    "boosting_type": "gbdt",
    "n_estimators": 10000,
    "n_jobs": -1,
    "verbose": -1,
}

random.seed(12345)
```


```python
xgb_params = {
    #booster="gbtree",
    "objective": "reg:absoluteerror",
    "eval_metric": "mae",
    "n_estimators": 5000,
    "early_stopping_rounds": 20,
    "seed": 12345,
    "max_depth": 8,
    "n_jobs": -1,
    "verbosity": 0,
}
```


```python
## conditions
data_source_name_list = ["pbev8", "r2scan"] 
data_source_list = ["dataset/qm9_vapor_pbev8_H.csv", "dataset/qm9_vapor_r2scan_H.csv"]
preprocess_list = ["MMN+PCA_white", "STN+PCA_white"]  #None, Normarize
scaler_name_list = ["minmax", "standard"]
model_name_list =["lgb", "xgb"]  # rf, svm
num_feature_list = [8, 16, 32, 64, 128, -1] # 5, 10, 25, 50, 100, 250, 8, 16, 32, 64, 128, 
data_size_list = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 2000000] # 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 
```


```python
for data_source, data_source_name in zip(data_source_list, data_source_name_list):    
    for preprocess in preprocess_list:
        df_h = pd.read_csv(data_source)
        ex_cols = [x for x in df_h.columns if x.startswith("desc")]
        rest_col = [x for x in df_h.columns if x not in ex_cols]
        t_col = "value"
        group_col = "molecule_idx"
        X = df_h[ex_cols]

        if preprocess is not None:
            if preprocess.startswith("MMN"):
                scaler = MinMaxScaler()
            elif preprocess.startswith("STN"):
                scaler = StandardScaler()
            
            X_norm = X.copy()
            X_norm.iloc[:,:] = scaler.fit_transform(X)
            
            if "PCA" in preprocess:
                white_flag = True if preprocess.endswith("white") else False
                pca = PCA(whiten=white_flag)
                X_norm.iloc[:,:] = pca.fit_transform(X_norm)
            df_norm = pd.concat([X_norm, df_h[rest_col]], axis=1)
        else:
            scaler_name = "non"
            df_norm = df_h
            X_norm = X

        for data_size in data_size_list:
            mol_list = random.sample(list(df_norm["molecule_idx"].unique()), 600)
    
            df_early = (
                df_norm[df_norm["molecule_idx"].isin(mol_list)]
                .reset_index()
                .drop(columns="index")
            )
            df_cv = (
                df_norm[~df_norm["molecule_idx"].isin(mol_list)]
                .reset_index()
                .drop(columns="index")
            )
            data_size = data_size if len(df_cv.index) &gt; data_size else len(df_cv.index)
            input_id_list = random.sample(range(len(df_cv.index)), data_size)
            df_cv = df_cv.loc[input_id_list, :].reset_index().drop(columns="index")
            
            #del df_norm
            
            for num_feat in num_feature_list:
                if num_feat != -1:
                    if preprocess is None:
                        continue
                    selected = [True] * num_feat + [False] * (X_norm.shape[1] - num_feat)
                else:
                    selected = [True] * X_norm.shape[1]
    
                X = df_cv.copy()[ex_cols]
                X = X[X.columns[selected]]
                y = df_cv.copy()[t_col]
                group = df_cv[group_col]
    
                X_ea = df_early.copy()[ex_cols]
                X_ea = X_ea[X_ea.columns[selected]]
                y_ea = df_early.copy()[t_col]
                
                for model_name in model_name_list:                
                    kf = GroupKFold(n_splits=5)
                    res_dict = {}
    
                    for fold, (t_i, v_i) in enumerate(kf.split(X, y, group)):
                        X_train, X_valid = X.loc[t_i, :], X.loc[v_i, :]
                        y_train, y_valid = y[t_i], y[v_i]
    
                        if model_name == "rf":
                            model = RandomForestRegressor(n_jobs=-1, max_depth=10)
                        elif model_name == "lgb":
                            model = lgb.LGBMRegressor(**lgb_params)
                        elif model_name == "xgb":
                            model = xgb.XGBRegressor(**xgb_params)
                        
                        if isinstance(model, type(lgb.LGBMRegressor())):
                            model.fit(
                                X_train,
                                y_train,
                                eval_set=[(X_ea, y_ea)],
                                callbacks=[
                                    lgb.early_stopping(stopping_rounds=20, verbose=True),
                                    lgb.log_evaluation(100),
                                ],
                            )
                        elif isinstance(model, type(xgb.XGBRegressor())):
                            model.fit(
                                X_train,
                                y_train,
                                eval_set=[(X_ea, y_ea)],
                            )
                        else:
                            model.fit(X_train, y_train)
    
                        res_dict.update({fold: [model, X_train, X_valid, y_train, y_valid]})
                    
                    print(f"{data_source_name}_{preprocess}_{model_name}_{num_feat} completed")
    
                    # evaluate
                    train_res_list = []
                    val_res_list = []
    
                    for i in range(5):
                        model, X_train, X_valid, y_train, y_valid = res_dict[i]
                        train_res_list.append(mean_absolute_error(model.predict(X_train), y_train))
                        val_res_list.append(mean_absolute_error(model.predict(X_valid), y_valid))
    
                    with open("res_data_size_feature.txt", "a") as f:
                        cnt = f.write(f"{data_source_name}_{preprocess}_{model_name}_d{data_size}_f{num_feat} Train:{np.mean(train_res_list):.3f} +- {np.std(train_res_list):.3f}, Val:{np.mean(val_res_list):.3f} +- {np.std(val_res_list):.3f}\n")
            del df_cv, df_early
        del X_norm, X, X_ea, y, y_ea, df_norm
        del res_dict
print("completed")
```
