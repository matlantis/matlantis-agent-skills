# PFP descriptorsを用いたTrajectoryの次元削減と可視化

このチュートリアルでは、PFP descriptorsを使用してTrajectoryを次元削減し、可視化する方法を紹介します。  
機械学習モデルを用いて、物性を予測する際にはPFP記述子をsystem-levelに集約する必要がありましたが、シミュレーションのTrajectoryでは原子数や原子種が変わらない場合が多いのでPFP descriptorsをそのまま使用できます。  
今回は、分子構造を最適化した際のTrajectoryの全ての構造をPFP descriptorsを用いてベクトル化し、次元削減して可視化していきます。



## 依存関係のインストール

今回、次元削減には[LocalMAP](https://arxiv.org/pdf/2412.15426) (Dimension Reduction with Locally Adjusted Graphs, Wang et al. 2025から)という手法を用います。この手法は局所的な情報だけでなく、大域的な情報も保存して次元削減できると言われています。  
Matlantisノートブック環境に事前インストールされていないため、まずインストールします。  
また、pfcc-extrasは`pfcc-extras &gt;= 0.12.3` をご利用ください。


```python
%pip install pacmap "pfcc-extras&gt;=0.12.3" "pfp-api-client&gt;=2.0.2"
```


```python
import datetime
import numpy as np
import os
import pandas as pd
import random

from ase.optimize import FIRE
from ase.io import Trajectory
from pacmap import LocalMAP
import plotly.express as px
from rdkit import Chem
from tqdm import tqdm

from pfcc_extras.structure.ase_rdkit_converter import smiles_to_atoms
from pfcc_extras.visualize import view_ngl
from pfp_api_client import ASECalculator, Estimator, EstimatorCalcMode
```


    


    /home/jovyan/.py313/lib/python3.13/site-packages/nglview/__init__.py:12: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools&lt;81.
      import pkg_resources




## .trajファイルの準備

分子構造を最適化し、出力されたtrajectoryファイルを用いて解析を行います。
今回は、`L-Alanyl-L-alanine`を対象にします。これは約10分かかりますが代わりにassetsフォルダーにある`L-Alanyl-L-alanine_opt.traj`をご利用でも結構です。


```python
# L-Alanyl-L-alanine
smiles ="C[C@@H](C(=O)N[C@@H](C)C(=O)O)N"
Chem.MolFromSmiles(smiles)
```




    
![png](output_5_0.png)
    




```python
#%%time
#def timestamp():
#    time_cur = datetime.datetime.now()
#    print("datetime : ", time_cur.strftime("%y/%m/%d/%H:%M:%S"))
#    stamp = time_cur.strftime("%y%m%d%H%M%S")
#    return stamp
#
#output_file = f"L-Alanyl-L-alanine_opt_{timestamp()}"
#
## opt structure
#atoms = smiles_to_atoms(smiles, randomSeed = SEED_VALUE, maxAttempts=10_000)
#atoms.calc = calculator
#opt = FIRE(atoms, trajectory=f"{output_file}.traj", logfile=f"{output_file}.log")
#opt.run(fmax=0.0001, steps=10_000)
```


```python
# modify this to your own trajectory if you ran the previous cell
traj = Trajectory("assets/trajectory_visualization_with_descriptors/L-Alanyl-L-alanine_opt.traj")
view_ngl(traj, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(max_frame=9329), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'C'…





## dataframeの準備

Trajectoryファイルが準備できたので、Trajectoryファイル内の各構造に対して、PFP descriptorsを計算し、dataframeにします。  
このとき、potential energyも可視化の際に使用するため、dataframeに追加します。
その後、`fit_transform()` を使用して LocalMAP 次元削減アルゴリズムを実行し、結果の埋め込み空間を視覚化します。


```python
# seed settings
SEED_VALUE = 0

# PFP settings
model_version = "v8.0.0"
mode = "R2SCAN"

random.seed(SEED_VALUE)
np.random.seed(SEED_VALUE)

estimator = Estimator(model_version=model_version, calc_mode=mode)
calculator = ASECalculator(estimator)

# create pandas dataframe containing the trajectory using PFP descriptors and potential energy
data_list = []

# Obtain descriptors from PFP
for i in tqdm(range(len(traj))):

    energy = traj[i].get_potential_energy()
    temp = calculator.get_descriptors(traj[i])
    descs = temp["scalar"].reshape(1, -1)

    energy_array = np.array([[energy]])
    combined_data = np.hstack([energy_array, descs])
    
    data_list.append(combined_data[0])

num_descriptors = descs.shape[1]
column_names = ['potential_energy'] + [f'desc_{j}' for j in range(num_descriptors)]

df = pd.DataFrame(data_list, columns=column_names)
df.head()  # 5 rows × 5889 columns = n_atoms x 256 descriptor features per atom
```

    100%|██████████| 9330/9330 [08:27&lt;00:00, 18.37it/s]





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
      &lt;th&gt;potential_energy&lt;/th&gt;
      &lt;th&gt;desc_0&lt;/th&gt;
      &lt;th&gt;desc_1&lt;/th&gt;
      &lt;th&gt;desc_2&lt;/th&gt;
      &lt;th&gt;desc_3&lt;/th&gt;
      &lt;th&gt;desc_4&lt;/th&gt;
      &lt;th&gt;desc_5&lt;/th&gt;
      &lt;th&gt;desc_6&lt;/th&gt;
      &lt;th&gt;desc_7&lt;/th&gt;
      &lt;th&gt;desc_8&lt;/th&gt;
      &lt;th&gt;...&lt;/th&gt;
      &lt;th&gt;desc_5878&lt;/th&gt;
      &lt;th&gt;desc_5879&lt;/th&gt;
      &lt;th&gt;desc_5880&lt;/th&gt;
      &lt;th&gt;desc_5881&lt;/th&gt;
      &lt;th&gt;desc_5882&lt;/th&gt;
      &lt;th&gt;desc_5883&lt;/th&gt;
      &lt;th&gt;desc_5884&lt;/th&gt;
      &lt;th&gt;desc_5885&lt;/th&gt;
      &lt;th&gt;desc_5886&lt;/th&gt;
      &lt;th&gt;desc_5887&lt;/th&gt;
    &lt;/tr&gt;
  &lt;/thead&gt;
  &lt;tbody&gt;
    &lt;tr&gt;
      &lt;th&gt;0&lt;/th&gt;
      &lt;td&gt;-98.176785&lt;/td&gt;
      &lt;td&gt;-0.855506&lt;/td&gt;
      &lt;td&gt;0.045976&lt;/td&gt;
      &lt;td&gt;-0.174566&lt;/td&gt;
      &lt;td&gt;-0.114590&lt;/td&gt;
      &lt;td&gt;1.155891&lt;/td&gt;
      &lt;td&gt;-0.061934&lt;/td&gt;
      &lt;td&gt;0.135960&lt;/td&gt;
      &lt;td&gt;0.111798&lt;/td&gt;
      &lt;td&gt;0.175189&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.015548&lt;/td&gt;
      &lt;td&gt;0.080175&lt;/td&gt;
      &lt;td&gt;-0.262020&lt;/td&gt;
      &lt;td&gt;-0.292123&lt;/td&gt;
      &lt;td&gt;-0.029150&lt;/td&gt;
      &lt;td&gt;0.048947&lt;/td&gt;
      &lt;td&gt;-0.130301&lt;/td&gt;
      &lt;td&gt;-0.099153&lt;/td&gt;
      &lt;td&gt;0.209052&lt;/td&gt;
      &lt;td&gt;-0.294981&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;1&lt;/th&gt;
      &lt;td&gt;-98.750714&lt;/td&gt;
      &lt;td&gt;-0.870115&lt;/td&gt;
      &lt;td&gt;0.056059&lt;/td&gt;
      &lt;td&gt;-0.170057&lt;/td&gt;
      &lt;td&gt;-0.108194&lt;/td&gt;
      &lt;td&gt;1.187468&lt;/td&gt;
      &lt;td&gt;-0.076991&lt;/td&gt;
      &lt;td&gt;0.125620&lt;/td&gt;
      &lt;td&gt;0.128484&lt;/td&gt;
      &lt;td&gt;0.159751&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.012423&lt;/td&gt;
      &lt;td&gt;0.052388&lt;/td&gt;
      &lt;td&gt;-0.256430&lt;/td&gt;
      &lt;td&gt;-0.291630&lt;/td&gt;
      &lt;td&gt;-0.038580&lt;/td&gt;
      &lt;td&gt;0.042853&lt;/td&gt;
      &lt;td&gt;-0.128245&lt;/td&gt;
      &lt;td&gt;-0.109500&lt;/td&gt;
      &lt;td&gt;0.205922&lt;/td&gt;
      &lt;td&gt;-0.252962&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;2&lt;/th&gt;
      &lt;td&gt;-98.719696&lt;/td&gt;
      &lt;td&gt;-0.876572&lt;/td&gt;
      &lt;td&gt;0.069110&lt;/td&gt;
      &lt;td&gt;-0.166180&lt;/td&gt;
      &lt;td&gt;-0.073010&lt;/td&gt;
      &lt;td&gt;1.213117&lt;/td&gt;
      &lt;td&gt;-0.071084&lt;/td&gt;
      &lt;td&gt;0.133842&lt;/td&gt;
      &lt;td&gt;0.153512&lt;/td&gt;
      &lt;td&gt;0.146715&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.011038&lt;/td&gt;
      &lt;td&gt;0.024367&lt;/td&gt;
      &lt;td&gt;-0.251365&lt;/td&gt;
      &lt;td&gt;-0.289564&lt;/td&gt;
      &lt;td&gt;-0.041478&lt;/td&gt;
      &lt;td&gt;0.040277&lt;/td&gt;
      &lt;td&gt;-0.116585&lt;/td&gt;
      &lt;td&gt;-0.123747&lt;/td&gt;
      &lt;td&gt;0.186281&lt;/td&gt;
      &lt;td&gt;-0.201915&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;3&lt;/th&gt;
      &lt;td&gt;-98.886978&lt;/td&gt;
      &lt;td&gt;-0.875966&lt;/td&gt;
      &lt;td&gt;0.068305&lt;/td&gt;
      &lt;td&gt;-0.169243&lt;/td&gt;
      &lt;td&gt;-0.083279&lt;/td&gt;
      &lt;td&gt;1.214657&lt;/td&gt;
      &lt;td&gt;-0.073376&lt;/td&gt;
      &lt;td&gt;0.139540&lt;/td&gt;
      &lt;td&gt;0.144769&lt;/td&gt;
      &lt;td&gt;0.149552&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.012039&lt;/td&gt;
      &lt;td&gt;0.030077&lt;/td&gt;
      &lt;td&gt;-0.252459&lt;/td&gt;
      &lt;td&gt;-0.289815&lt;/td&gt;
      &lt;td&gt;-0.039420&lt;/td&gt;
      &lt;td&gt;0.046315&lt;/td&gt;
      &lt;td&gt;-0.118689&lt;/td&gt;
      &lt;td&gt;-0.121864&lt;/td&gt;
      &lt;td&gt;0.189454&lt;/td&gt;
      &lt;td&gt;-0.209639&lt;/td&gt;
    &lt;/tr&gt;
    &lt;tr&gt;
      &lt;th&gt;4&lt;/th&gt;
      &lt;td&gt;-98.994454&lt;/td&gt;
      &lt;td&gt;-0.873594&lt;/td&gt;
      &lt;td&gt;0.066165&lt;/td&gt;
      &lt;td&gt;-0.173797&lt;/td&gt;
      &lt;td&gt;-0.099921&lt;/td&gt;
      &lt;td&gt;1.215143&lt;/td&gt;
      &lt;td&gt;-0.076155&lt;/td&gt;
      &lt;td&gt;0.147933&lt;/td&gt;
      &lt;td&gt;0.130804&lt;/td&gt;
      &lt;td&gt;0.153724&lt;/td&gt;
      &lt;td&gt;...&lt;/td&gt;
      &lt;td&gt;0.013637&lt;/td&gt;
      &lt;td&gt;0.039514&lt;/td&gt;
      &lt;td&gt;-0.253985&lt;/td&gt;
      &lt;td&gt;-0.290382&lt;/td&gt;
      &lt;td&gt;-0.035626&lt;/td&gt;
      &lt;td&gt;0.055569&lt;/td&gt;
      &lt;td&gt;-0.121586&lt;/td&gt;
      &lt;td&gt;-0.118354&lt;/td&gt;
      &lt;td&gt;0.194261&lt;/td&gt;
      &lt;td&gt;-0.223282&lt;/td&gt;
    &lt;/tr&gt;
  &lt;/tbody&gt;
&lt;/table&gt;
&lt;p&gt;5 rows × 5889 columns&lt;/p&gt;
&lt;/div&gt;




```python
descs_array = df.iloc[:,1:].values
energy_array = df.iloc[:,0].values

# Run dimension reduction algorithm
embedding = LocalMAP(n_components=2, n_neighbors=10, MN_ratio=0.5, FP_ratio=2.0, random_state=SEED_VALUE) 
lmap = embedding.fit_transform(descs_array, init="pca")
```

    Warning: random state is set to 0.



```python
df["LocalMAP_1"] = lmap[:,0]
df["LocalMAP_2"] = lmap[:,1]

# Calculate the interquartile range (IQR)

q1 = np.percentile(energy_array, 25)
q3 = np.percentile(energy_array, 75)
iqr = q3 - q1

# Calculate the upper and lower bound of the box plot
upper_bound = q3 + 1.5 * iqr
lower_bound = q1 - 1.5 * iqr


# Visualize the dimension-reduced embedding space.
fig = px.scatter(df,
                 x="LocalMAP_1",
                 y="LocalMAP_2",
                 color="potential_energy",
                 log_x=False,
                 hover_name=df.index,
                 width=800,  
                 height=600,
                 color_continuous_scale='jet',
                 range_color=[lower_bound, upper_bound]
                 )

fig.show()
```


    
![png](output_11_0.png)
    




## 結果

同じクラスターにはポテンシャルエネルギーの近い構造が集まっており、局所的な情報が保存されていることがわかります。加えて、隣接するクラスター同士もエネルギー的に近いことから、大域的な情報も保持されていると考えられます。

今回は一分子の構造最適化を例としましたが、この解析手法をMD計算のトラジェクトリなどに用いることで、シミュレーション中に系がどのような状態を遷移したのかを可視化できるでしょう。

なお、観察された不連続性は主に埋め込みアルゴリズム（LocalMAP）に起因するものであり、基盤となる構造の進化によるものではないことにご留意ください。連続的なTrajectoryが望ましい場合には、PCAなどの連続的な埋め込み手法がより適しています。Trajectoryの連続性が不要な場合、グラフベースの局所構造的な視点が得られる点でLocalMAPが有利となる可能性があります。
