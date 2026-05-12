Copyright Matlantis Corp. as contributors to Matlantis contrib project


```python
# If the libraries are not enough to run the following script, please run the install command below
#!pip install lightgbm xgboost pycaret
```


```python
import warnings
import os
warnings.filterwarnings("ignore")
from pycaret.regression import *
import pandas as pd
import pickle
```

## 4. Screening with pycaret
The first screening was conducted with pycaret.  
The dataset was divided by groupkfold with moledule id to compensate for the accuracy against unknown molecules.  


```python
dataset_file = "dataset/qm9_vapor_r2scan_H.csv"

model_type = dataset_file.split("_")[-2]
seed = 12345
model_type
```




    'r2scan'




```python
df_all = pd.read_csv(dataset_file)
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




```python
reg = setup(df_all, target ="value", fold = 5, 
            ignore_features = ["molecule_idx"],
            fold_strategy="groupkfold", fold_groups="molecule_idx",
            )
```


&lt;style type="text/css"&gt;
#T_61891_row9_col1 {
  background-color: lightgreen;
}
&lt;/style&gt;
&lt;table id="T_61891"&gt;
  &lt;thead&gt;
    &lt;tr&gt;
      &lt;th class="blank level0" &gt;&amp;nbsp;&lt;/th&gt;
      &lt;th id="T_61891_level0_col0" class="col_heading level0 col0" &gt;Description&lt;/th&gt;
      &lt;th id="T_61891_level0_col1" class="col_heading level0 col1" &gt;Value&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row0" class="row_heading level0 row0" &gt;0&lt;/th&gt;
      &lt;td id="T_61891_row0_col0" class="data row0 col0" &gt;Session id&lt;/td&gt;
      &lt;td id="T_61891_row0_col1" class="data row0 col1" &gt;3218&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row1" class="row_heading level0 row1" &gt;1&lt;/th&gt;
      &lt;td id="T_61891_row1_col0" class="data row1 col0" &gt;Target&lt;/td&gt;
      &lt;td id="T_61891_row1_col1" class="data row1 col1" &gt;value&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row2" class="row_heading level0 row2" &gt;2&lt;/th&gt;
      &lt;td id="T_61891_row2_col0" class="data row2 col0" &gt;Target type&lt;/td&gt;
      &lt;td id="T_61891_row2_col1" class="data row2 col1" &gt;Regression&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row3" class="row_heading level0 row3" &gt;3&lt;/th&gt;
      &lt;td id="T_61891_row3_col0" class="data row3 col0" &gt;Original data shape&lt;/td&gt;
      &lt;td id="T_61891_row3_col1" class="data row3 col1" &gt;(1208463, 258)&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row4" class="row_heading level0 row4" &gt;4&lt;/th&gt;
      &lt;td id="T_61891_row4_col0" class="data row4 col0" &gt;Transformed data shape&lt;/td&gt;
      &lt;td id="T_61891_row4_col1" class="data row4 col1" &gt;(1208463, 257)&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row5" class="row_heading level0 row5" &gt;5&lt;/th&gt;
      &lt;td id="T_61891_row5_col0" class="data row5 col0" &gt;Transformed train set shape&lt;/td&gt;
      &lt;td id="T_61891_row5_col1" class="data row5 col1" &gt;(845924, 257)&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row6" class="row_heading level0 row6" &gt;6&lt;/th&gt;
      &lt;td id="T_61891_row6_col0" class="data row6 col0" &gt;Transformed test set shape&lt;/td&gt;
      &lt;td id="T_61891_row6_col1" class="data row6 col1" &gt;(362539, 257)&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row7" class="row_heading level0 row7" &gt;7&lt;/th&gt;
      &lt;td id="T_61891_row7_col0" class="data row7 col0" &gt;Ignore features&lt;/td&gt;
      &lt;td id="T_61891_row7_col1" class="data row7 col1" &gt;1&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row8" class="row_heading level0 row8" &gt;8&lt;/th&gt;
      &lt;td id="T_61891_row8_col0" class="data row8 col0" &gt;Numeric features&lt;/td&gt;
      &lt;td id="T_61891_row8_col1" class="data row8 col1" &gt;256&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row9" class="row_heading level0 row9" &gt;9&lt;/th&gt;
      &lt;td id="T_61891_row9_col0" class="data row9 col0" &gt;Preprocess&lt;/td&gt;
      &lt;td id="T_61891_row9_col1" class="data row9 col1" &gt;True&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row10" class="row_heading level0 row10" &gt;10&lt;/th&gt;
      &lt;td id="T_61891_row10_col0" class="data row10 col0" &gt;Imputation type&lt;/td&gt;
      &lt;td id="T_61891_row10_col1" class="data row10 col1" &gt;simple&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row11" class="row_heading level0 row11" &gt;11&lt;/th&gt;
      &lt;td id="T_61891_row11_col0" class="data row11 col0" &gt;Numeric imputation&lt;/td&gt;
      &lt;td id="T_61891_row11_col1" class="data row11 col1" &gt;mean&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row12" class="row_heading level0 row12" &gt;12&lt;/th&gt;
      &lt;td id="T_61891_row12_col0" class="data row12 col0" &gt;Categorical imputation&lt;/td&gt;
      &lt;td id="T_61891_row12_col1" class="data row12 col1" &gt;mode&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row13" class="row_heading level0 row13" &gt;13&lt;/th&gt;
      &lt;td id="T_61891_row13_col0" class="data row13 col0" &gt;Fold Generator&lt;/td&gt;
      &lt;td id="T_61891_row13_col1" class="data row13 col1" &gt;GroupKFold&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row14" class="row_heading level0 row14" &gt;14&lt;/th&gt;
      &lt;td id="T_61891_row14_col0" class="data row14 col0" &gt;Fold Number&lt;/td&gt;
      &lt;td id="T_61891_row14_col1" class="data row14 col1" &gt;5&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row15" class="row_heading level0 row15" &gt;15&lt;/th&gt;
      &lt;td id="T_61891_row15_col0" class="data row15 col0" &gt;CPU Jobs&lt;/td&gt;
      &lt;td id="T_61891_row15_col1" class="data row15 col1" &gt;-1&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row16" class="row_heading level0 row16" &gt;16&lt;/th&gt;
      &lt;td id="T_61891_row16_col0" class="data row16 col0" &gt;Use GPU&lt;/td&gt;
      &lt;td id="T_61891_row16_col1" class="data row16 col1" &gt;False&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row17" class="row_heading level0 row17" &gt;17&lt;/th&gt;
      &lt;td id="T_61891_row17_col0" class="data row17 col0" &gt;Log Experiment&lt;/td&gt;
      &lt;td id="T_61891_row17_col1" class="data row17 col1" &gt;False&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row18" class="row_heading level0 row18" &gt;18&lt;/th&gt;
      &lt;td id="T_61891_row18_col0" class="data row18 col0" &gt;Experiment Name&lt;/td&gt;
      &lt;td id="T_61891_row18_col1" class="data row18 col1" &gt;reg-default-name&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_61891_level0_row19" class="row_heading level0 row19" &gt;19&lt;/th&gt;
      &lt;td id="T_61891_row19_col0" class="data row19 col0" &gt;USI&lt;/td&gt;
      &lt;td id="T_61891_row19_col1" class="data row19 col1" &gt;23d4&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;




```python
best_model = compare_models(sort="MAE", n_select=3, 
                            include=["knn", "xgboost", "lightgbm", "lr", "ridge", "lasso", "mlp"])
best_model
```






&lt;style type="text/css"&gt;
#T_ee7dd th {
  text-align: left;
}
#T_ee7dd_row0_col0, #T_ee7dd_row0_col2, #T_ee7dd_row0_col3, #T_ee7dd_row0_col4, #T_ee7dd_row0_col5, #T_ee7dd_row1_col0, #T_ee7dd_row1_col1, #T_ee7dd_row1_col6, #T_ee7dd_row2_col0, #T_ee7dd_row2_col1, #T_ee7dd_row2_col2, #T_ee7dd_row2_col3, #T_ee7dd_row2_col4, #T_ee7dd_row2_col5, #T_ee7dd_row2_col6, #T_ee7dd_row3_col0, #T_ee7dd_row3_col1, #T_ee7dd_row3_col2, #T_ee7dd_row3_col3, #T_ee7dd_row3_col4, #T_ee7dd_row3_col5, #T_ee7dd_row3_col6, #T_ee7dd_row4_col0, #T_ee7dd_row4_col1, #T_ee7dd_row4_col2, #T_ee7dd_row4_col3, #T_ee7dd_row4_col4, #T_ee7dd_row4_col5, #T_ee7dd_row4_col6, #T_ee7dd_row5_col0, #T_ee7dd_row5_col1, #T_ee7dd_row5_col2, #T_ee7dd_row5_col3, #T_ee7dd_row5_col4, #T_ee7dd_row5_col5, #T_ee7dd_row5_col6, #T_ee7dd_row6_col0, #T_ee7dd_row6_col1, #T_ee7dd_row6_col2, #T_ee7dd_row6_col3, #T_ee7dd_row6_col4, #T_ee7dd_row6_col5, #T_ee7dd_row6_col6 {
  text-align: left;
}
#T_ee7dd_row0_col1, #T_ee7dd_row0_col6, #T_ee7dd_row1_col2, #T_ee7dd_row1_col3, #T_ee7dd_row1_col4, #T_ee7dd_row1_col5 {
  text-align: left;
  background-color: yellow;
}
#T_ee7dd_row0_col7, #T_ee7dd_row1_col7, #T_ee7dd_row2_col7, #T_ee7dd_row3_col7, #T_ee7dd_row4_col7, #T_ee7dd_row5_col7 {
  text-align: left;
  background-color: lightgrey;
}
#T_ee7dd_row6_col7 {
  text-align: left;
  background-color: yellow;
  background-color: lightgrey;
}
&lt;/style&gt;
&lt;table id="T_ee7dd"&gt;
  &lt;thead&gt;
    &lt;tr&gt;
      &lt;th class="blank level0" &gt;&amp;nbsp;&lt;/th&gt;
      &lt;th id="T_ee7dd_level0_col0" class="col_heading level0 col0" &gt;Model&lt;/th&gt;
      &lt;th id="T_ee7dd_level0_col1" class="col_heading level0 col1" &gt;MAE&lt;/th&gt;
      &lt;th id="T_ee7dd_level0_col2" class="col_heading level0 col2" &gt;MSE&lt;/th&gt;
      &lt;th id="T_ee7dd_level0_col3" class="col_heading level0 col3" &gt;RMSE&lt;/th&gt;
      &lt;th id="T_ee7dd_level0_col4" class="col_heading level0 col4" &gt;R2&lt;/th&gt;
      &lt;th id="T_ee7dd_level0_col5" class="col_heading level0 col5" &gt;RMSLE&lt;/th&gt;
      &lt;th id="T_ee7dd_level0_col6" class="col_heading level0 col6" &gt;MAPE&lt;/th&gt;
      &lt;th id="T_ee7dd_level0_col7" class="col_heading level0 col7" &gt;TT (Sec)&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th id="T_ee7dd_level0_row0" class="row_heading level0 row0" &gt;knn&lt;/th&gt;
      &lt;td id="T_ee7dd_row0_col0" class="data row0 col0" &gt;K Neighbors Regressor&lt;/td&gt;
      &lt;td id="T_ee7dd_row0_col1" class="data row0 col1" &gt;0.1263&lt;/td&gt;
      &lt;td id="T_ee7dd_row0_col2" class="data row0 col2" &gt;0.0381&lt;/td&gt;
      &lt;td id="T_ee7dd_row0_col3" class="data row0 col3" &gt;0.1953&lt;/td&gt;
      &lt;td id="T_ee7dd_row0_col4" class="data row0 col4" &gt;0.9903&lt;/td&gt;
      &lt;td id="T_ee7dd_row0_col5" class="data row0 col5" &gt;0.0067&lt;/td&gt;
      &lt;td id="T_ee7dd_row0_col6" class="data row0 col6" &gt;0.0044&lt;/td&gt;
      &lt;td id="T_ee7dd_row0_col7" class="data row0 col7" &gt;463.1520&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_ee7dd_level0_row1" class="row_heading level0 row1" &gt;mlp&lt;/th&gt;
      &lt;td id="T_ee7dd_row1_col0" class="data row1 col0" &gt;MLP Regressor&lt;/td&gt;
      &lt;td id="T_ee7dd_row1_col1" class="data row1 col1" &gt;0.1435&lt;/td&gt;
      &lt;td id="T_ee7dd_row1_col2" class="data row1 col2" &gt;0.0374&lt;/td&gt;
      &lt;td id="T_ee7dd_row1_col3" class="data row1 col3" &gt;0.1934&lt;/td&gt;
      &lt;td id="T_ee7dd_row1_col4" class="data row1 col4" &gt;0.9905&lt;/td&gt;
      &lt;td id="T_ee7dd_row1_col5" class="data row1 col5" &gt;0.0066&lt;/td&gt;
      &lt;td id="T_ee7dd_row1_col6" class="data row1 col6" &gt;0.0050&lt;/td&gt;
      &lt;td id="T_ee7dd_row1_col7" class="data row1 col7" &gt;293.9660&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_ee7dd_level0_row2" class="row_heading level0 row2" &gt;xgboost&lt;/th&gt;
      &lt;td id="T_ee7dd_row2_col0" class="data row2 col0" &gt;Extreme Gradient Boosting&lt;/td&gt;
      &lt;td id="T_ee7dd_row2_col1" class="data row2 col1" &gt;0.2409&lt;/td&gt;
      &lt;td id="T_ee7dd_row2_col2" class="data row2 col2" &gt;0.1070&lt;/td&gt;
      &lt;td id="T_ee7dd_row2_col3" class="data row2 col3" &gt;0.3270&lt;/td&gt;
      &lt;td id="T_ee7dd_row2_col4" class="data row2 col4" &gt;0.9728&lt;/td&gt;
      &lt;td id="T_ee7dd_row2_col5" class="data row2 col5" &gt;0.0112&lt;/td&gt;
      &lt;td id="T_ee7dd_row2_col6" class="data row2 col6" &gt;0.0084&lt;/td&gt;
      &lt;td id="T_ee7dd_row2_col7" class="data row2 col7" &gt;34.8340&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_ee7dd_level0_row3" class="row_heading level0 row3" &gt;lr&lt;/th&gt;
      &lt;td id="T_ee7dd_row3_col0" class="data row3 col0" &gt;Linear Regression&lt;/td&gt;
      &lt;td id="T_ee7dd_row3_col1" class="data row3 col1" &gt;0.2739&lt;/td&gt;
      &lt;td id="T_ee7dd_row3_col2" class="data row3 col2" &gt;0.1364&lt;/td&gt;
      &lt;td id="T_ee7dd_row3_col3" class="data row3 col3" &gt;0.3693&lt;/td&gt;
      &lt;td id="T_ee7dd_row3_col4" class="data row3 col4" &gt;0.9654&lt;/td&gt;
      &lt;td id="T_ee7dd_row3_col5" class="data row3 col5" &gt;0.0128&lt;/td&gt;
      &lt;td id="T_ee7dd_row3_col6" class="data row3 col6" &gt;0.0096&lt;/td&gt;
      &lt;td id="T_ee7dd_row3_col7" class="data row3 col7" &gt;6.3460&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_ee7dd_level0_row4" class="row_heading level0 row4" &gt;ridge&lt;/th&gt;
      &lt;td id="T_ee7dd_row4_col0" class="data row4 col0" &gt;Ridge Regression&lt;/td&gt;
      &lt;td id="T_ee7dd_row4_col1" class="data row4 col1" &gt;0.2739&lt;/td&gt;
      &lt;td id="T_ee7dd_row4_col2" class="data row4 col2" &gt;0.1364&lt;/td&gt;
      &lt;td id="T_ee7dd_row4_col3" class="data row4 col3" &gt;0.3693&lt;/td&gt;
      &lt;td id="T_ee7dd_row4_col4" class="data row4 col4" &gt;0.9654&lt;/td&gt;
      &lt;td id="T_ee7dd_row4_col5" class="data row4 col5" &gt;0.0128&lt;/td&gt;
      &lt;td id="T_ee7dd_row4_col6" class="data row4 col6" &gt;0.0096&lt;/td&gt;
      &lt;td id="T_ee7dd_row4_col7" class="data row4 col7" &gt;3.7460&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_ee7dd_level0_row5" class="row_heading level0 row5" &gt;lightgbm&lt;/th&gt;
      &lt;td id="T_ee7dd_row5_col0" class="data row5 col0" &gt;Light Gradient Boosting Machine&lt;/td&gt;
      &lt;td id="T_ee7dd_row5_col1" class="data row5 col1" &gt;0.2801&lt;/td&gt;
      &lt;td id="T_ee7dd_row5_col2" class="data row5 col2" &gt;0.1442&lt;/td&gt;
      &lt;td id="T_ee7dd_row5_col3" class="data row5 col3" &gt;0.3797&lt;/td&gt;
      &lt;td id="T_ee7dd_row5_col4" class="data row5 col4" &gt;0.9634&lt;/td&gt;
      &lt;td id="T_ee7dd_row5_col5" class="data row5 col5" &gt;0.0131&lt;/td&gt;
      &lt;td id="T_ee7dd_row5_col6" class="data row5 col6" &gt;0.0098&lt;/td&gt;
      &lt;td id="T_ee7dd_row5_col7" class="data row5 col7" &gt;73.4020&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th id="T_ee7dd_level0_row6" class="row_heading level0 row6" &gt;lasso&lt;/th&gt;
      &lt;td id="T_ee7dd_row6_col0" class="data row6 col0" &gt;Lasso Regression&lt;/td&gt;
      &lt;td id="T_ee7dd_row6_col1" class="data row6 col1" &gt;1.5081&lt;/td&gt;
      &lt;td id="T_ee7dd_row6_col2" class="data row6 col2" &gt;3.9387&lt;/td&gt;
      &lt;td id="T_ee7dd_row6_col3" class="data row6 col3" &gt;1.9846&lt;/td&gt;
      &lt;td id="T_ee7dd_row6_col4" class="data row6 col4" &gt;-0.0000&lt;/td&gt;
      &lt;td id="T_ee7dd_row6_col5" class="data row6 col5" &gt;0.0700&lt;/td&gt;
      &lt;td id="T_ee7dd_row6_col6" class="data row6 col6" &gt;0.0546&lt;/td&gt;
      &lt;td id="T_ee7dd_row6_col7" class="data row6 col7" &gt;3.1180&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;










    [KNeighborsRegressor(n_jobs=-1),
     MLPRegressor(max_iter=500, random_state=3218),
     XGBRegressor(base_score=None, booster='gbtree', callbacks=None,
                  colsample_bylevel=None, colsample_bynode=None,
                  colsample_bytree=None, device='cpu', early_stopping_rounds=None,
                  enable_categorical=False, eval_metric=None, feature_types=None,
                  feature_weights=None, gamma=None, grow_policy=None,
                  importance_type=None, interaction_constraints=None,
                  learning_rate=None, max_bin=None, max_cat_threshold=None,
                  max_cat_to_onehot=None, max_delta_step=None, max_depth=None,
                  max_leaves=None, min_child_weight=None, missing=nan,
                  monotone_constraints=None, multi_strategy=None, n_estimators=None,
                  n_jobs=-1, num_parallel_tree=None, ...)]




```python
os.makedirs("models", exist_ok=True)
with open(f'models/pycaret_model_compared_result_{model_type}.pkl', 'wb') as f:
    pickle.dump(best_model, f)
```
