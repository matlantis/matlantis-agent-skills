# 誘電特性予測のためのPFP記述子の使用

このチュートリアルでは、PFP記述子を使用してシステムレベルの特性を予測する方法を、勾配ブースティング決定木 (gbdt) を用いて説明します。原子レベルの特徴は集約され、トレーニングおよび推論の両方で異なる数の原子を持つ個別のシステムを統合し、システムレベルの特性予測を可能にします。



## 依存関係のインストール

モデルを構築してトレーニングするために、`xgboost`と Optuna を使用します。これらのパッケージはMatlantisノートブック環境に事前インストールされていないため、まずそれらをインストールします。

`pfp-api-client &gt;= 2.0.2` をご利用ください。


```python
#%pip install -U monty optuna "pymatgen==2025.4.10"  scikit-learn xgboost
#%pip install -U pfp-api-client
#%pip install "nglview==3.1.2"  # included already in Matlantis environment, install if this is not the case.
```



## データ準備

この例では、Matbenchの誘電タスクからランダムに100個のデータポイントをサンプリングします。
参考: Petousis, I., Mrdjenovich, D., Ballouz, E., Liu, M., Winston, D., Chen, W., Graf, T., Schladt, T. D., Persson, K. A. &amp; Prinz, F. B. High-throughput screening of inorganic compounds for the discovery of novel dielectric and optical materials. Sci. Data 4, 160134 (2017).
データはJSONファイルに事前保存されており、各データポイントは以下の形式で格納されています。

- data:
    - id: `&lt;unique_identifier_1&gt;`
        - structure: `&lt;pymatgen.Structure_1&gt;`
        - target: `&lt;dielectric_value_1&gt;`
  
    - ...

    - id: `&lt;unique_identifier_100&gt;`
        - structure: `&lt;pymatgen.Structure_100&gt;`
        - target: `&lt;dielectric_value_100&gt;`


```python
import json
from monty.json import MontyDecoder

with open("assets/dielectric_prediction_with_descriptors/matbench_dielectric_samples.json", "r") as f:
    data = json.load(f, cls=MontyDecoder)
```



PFPでASECalculatorを初期化し、3D構造からスカラーPFP descriptorsを簡単に生成するための`get_pfp_descriptors`関数を定義します。
現在、PFP descriptorsの生成には最大4000個の原子数を持つ構造体が可能です。


```python
import pfp_api_client
print(pfp_api_client.__version__)
```

    2.2.0



```python
from pymatgen.io.ase import AseAtomsAdaptor
from pfp_api_client import ASECalculator, Estimator
import numpy as np
from tqdm import tqdm

# Model versions v7.0.0 and up are supported for descriptors.
estimator = Estimator(calc_mode='PBE')
calculator = ASECalculator(estimator)

def get_pfp_descriptors(struct) -&gt; np.ndarray:
    struct = AseAtomsAdaptor.get_atoms(struct)
    descriptor = calculator.get_descriptors(struct)
    return descriptor

for id_ in tqdm(data.keys(), desc='Generating PFP descriptors'):
    struct = data[id_]['structure']
    descriptor = get_pfp_descriptors(struct)
    data[id_]['descriptors'] = descriptor['scalar']
```

    Generating PFP descriptors: 100%|██████████| 100/100 [00:06&lt;00:00, 14.54it/s]




この例では、100個のデータポイントを単純に60個のトレーニングデータポイントと20個のバリデーションデータポイントと20個のテストデータポイントに分割します。


```python
ids = list(data.keys())
train_ids = ids[:60]
val_ids = ids[60:80]
test_ids = ids[80:]
train_targets = [data[id_]['target'] for id_ in train_ids]
val_targets = [data[id_]['target'] for id_ in val_ids]
test_targets = [data[id_]['target'] for id_ in test_ids]
```



## モデル訓練
このセクションでは、Optunaと組み合わせてXGBoostを使用し、シンプルな勾配ブースティング決定木モデルを訓練します。このモデルは記述子を入力として受け取り、ターゲット特性の予測を行います。


```python
def readout(pfp_descriptors):
    """
    Aggregates atom-level features into material representation by summing along the atom dimension.
    Args:
        pfp_descriptors (np.ndarray): A 2D array of shape (N, feature_dim) containing 
                                      the features of the atoms, where N is the number of atoms.
    Returns:
        np.ndarray: A 1D array of shape (feature_dim,) containing the summed features.
    """
    return np.sum(pfp_descriptors, axis=0)
```


```python
import optuna
import xgboost as xgb

from sklearn.metrics import mean_absolute_error


train_desc = [data[id_]['descriptors'] for id_ in train_ids]
val_desc = [data[id_]['descriptors'] for id_ in val_ids]
test_desc = [data[id_]['descriptors'] for id_ in test_ids]

# Aggregate descriptors with np.sum() to make the shape work for xgboost
train_desc_sum = [readout(data[id_]['descriptors']) for id_ in train_ids]
val_desc_sum = [readout(data[id_]['descriptors']) for id_ in val_ids]
test_desc_sum = [readout(data[id_]['descriptors']) for id_ in test_ids]


def objective(trial):
    params = {
            'n_estimators': trial.suggest_int('n_estimators', 100, 500),
            'max_depth': trial.suggest_int('max_depth', 3, 15),
            'learning_rate': trial.suggest_float('learning_rate', 1e-5, 0.1, log=True),
            #'subsample': trial.suggest_float('subsample', 0.5, 1.0),
            #'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
            #'gamma': trial.suggest_float('gamma', 0, 5),
            #'reg_alpha': trial.suggest_float('reg_alpha', 0, 5),
            #'reg_lambda': trial.suggest_float('reg_lambda', 0, 5),
            'random_state': 42,
            'eval_metric': mean_absolute_error,
    }
    
    model = xgb.XGBRegressor(**params)
    model.fit(train_desc_sum, train_targets, eval_set=[(val_desc_sum, val_targets)], verbose=0)
    
    val_pred = model.predict(val_desc_sum)
    return mean_absolute_error(val_targets, val_pred)


study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=10, timeout=600)
```


```python
best_params = study.best_params
best_params['random_state'] = 42
model = xgb.XGBRegressor(**best_params)
model.fit(train_desc_sum, train_targets)

test_pred = model.predict(test_desc_sum)
test_score = mean_absolute_error(test_targets, test_pred)
print(f"Test score: {test_score}")
```

    Test score: 0.33085280656814575



```python
import matplotlib.pyplot as plt

test_targets = test_targets
test_pred = test_pred.flatten()

min_v = min(min(test_targets), min(test_pred))
max_v = max(max(test_targets), max(test_pred))

plt.figure(figsize=(5, 5))
plt.scatter(test_targets, test_pred, label='Pred vs True', color='blue')

plt.plot([min_v, max_v], [min_v, max_v], 'r--', label='y = x')

plt.xlim([min_v, max_v])  # Set x-axis limits
plt.ylim([min_v, max_v])  # Set y-axis limits

plt.xlabel('True Values')
plt.ylabel('Predicted Values')
plt.title('Parity Plot')
plt.legend()

plt.show()
```


    
![png](output_14_0.png)
    




## SHAP値とPFP記述子を用いた予測の説明
* いくつかの原子を系統的に取り除き、残りの原子のPFP記述子の合計に基づいて、学習済みモデルが予測値を計算します。

* これらの予測値にカーネルSHAPを適用して、各原子の貢献度を決定します。
    * 赤色の原子：目的特性を増加させるのに貢献します。
    * 青色の原子：目的特性を減少させるのに貢献します。

* Shapley valuesを計算するアルゴリズムは、カーネル近似のためのサンプル数を決定した後に、いくつかのステップから構成されています。
  1. 特別な連合の設定: 全てオン と 全てオフ
  2. 重みが大きい連合からのサンプリング と　Shapley Kernelの重み計算
  3. 加重線形回帰


```python
import numpy as np
from scipy.special import binom
import itertools


def batched_value_function(Z_prime, atom_features, model, readout_func, batch_size=None):
    """
    Evaluates the model predictions for a batch of binary feature masks.

    Args:
        Z_prime (np.ndarray): A 2D array of shape (num_samples, M) containing binary masks (0s and 1s).
        atom_features (np.ndarray): A 2D array of shape (M, feature_dim) containing the original atom features.
        model (object): The machine learning model, which must have a `predict` method.
        readout_func (callable): A function that takes a masked 2D array of shape (M, feature_dim) 
                                 and aggregates it into a 1D array of shape (feature_dim,).
        batch_size (int, optional): The number of samples to process at once to prevent memory overflow. 
                                    If None, processes all samples in a single batch. Defaults to None.

    Returns:
        np.ndarray: A 1D array of shape (num_samples,) containing the model predictions.
    """
    num_samples = Z_prime.shape[0]
    
    # If batch_size is None, process all samples at once
    if batch_size is None:
        batch_size = num_samples

    all_values = []
    
    # Process the data in chunks of the specified batch_size
    for start_idx in range(0, num_samples, batch_size):
        end_idx = min(start_idx + batch_size, num_samples)
        Z_chunk = Z_prime[start_idx:end_idx]
        
        # 1. Mask processing for each chunk (using broadcasting)
        masked_chunk = Z_chunk[:, :, np.newaxis] * atom_features[np.newaxis, :, :]
        
        # 2. Aggregate features using the readout_func
        chunk_inputs = np.array([readout_func(feat) for feat in masked_chunk])
        
        # 3. Model inference for each chunk
        chunk_values = model.predict(chunk_inputs)
        all_values.append(chunk_values.flatten())
        
    # Concatenate the chunked inference results into a single 1D array and return
    return np.concatenate(all_values)


def calculate_weights(s_array, M):
    """
    Calculates the Kernel SHAP weights for a batch of coalitions.

    Args:
        s_array (np.ndarray): A 1D array containing the number of 'ON' features (s) for each coalition.
        M (int): The total number of unmasked features (atoms).

    Returns:
        np.ndarray: A 1D array containing the Shapley kernel weights.
    """
    weights = np.zeros_like(s_array, dtype=float)

    # Special cases: s == 0 or s == M get an infinitely large weight (approximated as 1e9)
    special_mask = (s_array == 0) | (s_array == M)
    weights[special_mask] = 1e9

    # Normal cases
    normal_mask = ~special_mask
    s_normal = s_array[normal_mask]

    # Calculate the Shapley kernel weight
    combinations_term = binom(M, s_normal)
    denominator = combinations_term * s_normal * (M - s_normal)    
    weights[normal_mask] = (M - 1) / denominator
    
    return weights


def calc_kernel_shap(atom_features, model, readout_func, n_kernel_samples=None, batch_size=None):
    """
    Calculates the Shapley values using the Kernel SHAP method with deterministic symmetric sampling.

    Args:
        atom_features (np.ndarray): A 2D array of shape (M, feature_dim) representing the atoms.
        model (object): The trained model to explain.
        readout_func (callable): The aggregation function to apply after masking features.
        n_kernel_samples (int, optional): The maximum number of coalitions to sample. 
                                          If None, a default heuristic based on M is used. Defaults to None.
        batch_size (int, optional): The number of samples to process at a time. Defaults to None.

    Returns:
        tuple:
            - shapley_values (np.ndarray): A 1D array of length M containing the Shapley values for each atom.
            - phi_0 (float): The expected base value (intercept) of the model prediction.
    """
    M = atom_features.shape[0] 
    if M == 0:
        raise ValueError("The input 'atom_features' contains no atoms.")

    max_samples = 2**M
    if n_kernel_samples is None:
        n_kernel_samples = min(max_samples, 2*M + 2048)
    else:
        n_kernel_samples = min(max_samples, n_kernel_samples)
    
    Z_prime_list = []

    # 1. Generate sampling patterns systematically
    # Add the coalitions with all features ON and all features OFF
    Z_prime_list.append(np.ones(M))
    Z_prime_list.append(np.zeros(M))

    max_s = M // 2
    for s in range(1, max_s + 1):
        if len(Z_prime_list) &gt;= n_kernel_samples:
            break

        # Generate all combinations of size s
        for combo in itertools.combinations(range(M), s):
                 
            # Create the ON pattern (set selected indices to 1)
            z_s = np.zeros(M)
            z_s[list(combo)] = 1
            Z_prime_list.append(z_s)

            # Create the corresponding OFF pattern (invert 0s and 1s) to maintain symmetry
            if s != M - s:
                Z_prime_list.append(1 - z_s)

            # Check if the desired number of samples is reached
            if len(Z_prime_list) &gt;= n_kernel_samples:
                break

    # Optional: Truncate exactly to n_kernel_samples if it overshot by 1 due to the pair appending
    Z_prime_list = Z_prime_list[:n_kernel_samples]
    
    # 2. Convert the list to a NumPy array in preparation for batch processing
    Z_prime = np.array(Z_prime_list)
    s_array = np.sum(Z_prime, axis=1)

    # 3. Calculate weights and prediction values (passing batch_size here)
    kernel_weights = calculate_weights(s_array, M)
    y_target = batched_value_function(Z_prime, atom_features, model, readout_func, batch_size)

    # 4. Perform weighted linear regression
    n_actual_samples = len(Z_prime)
    X_matrix = np.hstack((np.ones((n_actual_samples, 1)), Z_prime))
    
    sqrt_weights = np.sqrt(kernel_weights)
    X_weighted = X_matrix * sqrt_weights[:, np.newaxis] 
    y_weighted = y_target * sqrt_weights

    # Solve the least squares problem
    coeffs, _, _, _ = np.linalg.lstsq(X_weighted, y_weighted, rcond=None)

    # The first coefficient is the intercept (phi_0), the rest are the Shapley values
    shapley_values = coeffs[1:] 
    phi_0 = coeffs[0]
    
    return shapley_values, phi_0
```


```python
from ase.build import molecule
import nglview as nv

i = 7 # Change this to select a different test structure
struct = data[test_ids[i]]['structure']
atoms = AseAtomsAdaptor.get_atoms(struct)
print(atoms.symbols)

view_ase = nv.show_ase(atoms)
view_ase
```


    


    Tl4Sb4Se8



    NGLWidget()



```python
shap_values, phi_0 = calc_kernel_shap(test_desc[i], model, readout, n_kernel_samples=10000)
y_target = model.predict(test_desc_sum[i].reshape(1,-1))

print(f"predict : {y_target[0]:.9f}")
print(f"shapley : {(sum(shap_values) + phi_0):.9f}")

min_val = np.min(shap_values)
max_val = np.max(shap_values)
limit = max(abs(min_val), abs(max_val))

scaled_shap_values = shap_values / limit
atoms.set_array('bfactor', scaled_shap_values)

view_ase = nv.show_ase(atoms)
view_ase.clear_representations()


view_ase.add_representation(
    repr_type='ball+stick',
    selection='all',
    color_scheme='bfactor',
    color_scale='rwb', # Red-White-Blue
    color_domain=[1, 0, -1]
)

view_ase.center()
view_ase
```

    predict : 2.930350780
    shapley : 2.930350780



    NGLWidget()

