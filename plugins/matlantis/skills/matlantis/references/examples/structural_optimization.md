# 構造最適化

原子シミュレーションにおける構造最適化とは、エネルギーを最小化するように原子配置を最適化する数値的手法のことを指します。

構造最適化の計算では、系のポテンシャルエネルギーを最小化するように逐次的に原子配置を更新していきます。通常、最適化に必要なエネルギーや力の導出には、古典力場や第一原理計算によって得られた電子状態が利用されます。このexampleでは、PFPを用いて水素分子、水分子、シリコン結晶についてエネルギー的に最も安定な構造を探索します。

## 初期設定

* 必要なライブラリ等の準備を行います。


```python
# If you have already installed matlantis-features, you can skip this.
!pip install 'matlantis-features&gt;=0.19.3'
```






```python
import matlantis_features
import pfp_api_client
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

estimator = Estimator(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
calc = ASECalculator(estimator)

print(f"matlantis_features: {matlantis_features.__version__}")
print(f"pfp_api_client: {pfp_api_client.__version__}")
print(f"current pfp model version: {estimator.model_version}")
```

    matlantis_features: 1.0.1
    pfp_api_client: 2.0.1
    current pfp model version: v8.0.0




## 水素分子
水素分子は2つの水素原子からなる二原子分子です。H2分子におけるH-H結合の長さは実験により0.741Åと測定されています。
このexampleでは二原子分子H-Hの系のポテンシャルエネルギーを計算することでH-H結合に関する理解を深めます。その後、BFGS、LBFGS、FIRE、FIRE-LBFGS (matlantis-featuresで提供されている独自手法) それぞれの手法でH2分子構造を最適化し、安定構造を求めます。

### 水素分子のポテンシャルエネルギー
H2分子のポテンシャルエネルギーを原子間の距離の関数として描画します。ここでは0.5Åから5.0Åまで0.01Å間隔でH2構造を変化させることにします。ポテンシャルエネルギーはPFP calculatorを使うことで簡単に計算することができます。


```python
import numpy as np
from tqdm.auto import tqdm
from ase import Atoms

distance = np.arange(0.5, 5.0, 0.01)

energy = []
for d in tqdm(distance):
    atoms = Atoms(
        positions=np.array([[0.0, 0.0, 0.0],[d, 0.0, 0.0]]),
        numbers=np.array([1,1]),
        pbc=False
    )
    atoms.calc = calc
    energy.append(atoms.get_potential_energy())

```


      0%|          | 0/450 [00:00&lt;?, ?it/s]




ポテンシャルエネルギーは無限遠点、つまり原子間の相互作用が無くなる極限において0に近付きます。2原子が近づいていくとポテンシャルエネルギーは小さくなっていき、平衡結合長においてエネルギーが最小になります。そこからさらに原子同士が近づくと急速にエネルギーは大きくなります。平衡結合長はこのグラフでエネルギーが最小になる点を見ることで大まかに見積もることができます。


```python
import matplotlib.pyplot as plt

plt.plot(distance, energy)
plt.xlabel("H-H distance (Å)")
plt.ylabel("potential energy (eV)")
```




    Text(0, 0.5, 'potential energy (eV)')





![png](output_6_1.png)




```python
print(f"The H-H bond energy is minimized at {distance[np.argmin(energy)]:.2f} Å")
```

    The H-H bond energy is minimized at 0.75 Å




### BFGS optimizer

BFGSは、Broyden–Fletcher–Goldfarb–Shanno（BFGS）アルゴリズムの略で、広く使用されている準ニュートン法の一つです。 最適化中のトラジェクトリからヘシアンを近似し、近似ヘシアンを用いてNewton法を実行します。

[`matlantis-features`](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html)では、BFGS法を使用して構造最適化を実行する`BFGSASEOptFeature`が提供されています。 以下の例では、各原子に作用する力に関する収束基準`fmax`を0.001eV/オングストロームに設定します。 また、最適化ステップの最大数`n_run`を1000に設定します。収束基準を満たしていない場合でも、`n_run`ステップ後にシミュレーションが終了します。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from matlantis_features.features.common.opt import BFGSASEOptFeature

atoms = Atoms(
    positions=np.array([[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]]),
    numbers=np.array([1,1]),
    pbc=False
)

opt_feature = BFGSASEOptFeature(
    fmax=0.00001, n_run=1000, show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE),
    trajectory="H2_bfgs.traj"
)
opt_result = opt_feature(atoms)
```


      0%|          | 0/1001 [00:00&lt;?, ?it/s]




最適化された構造は `opt_result.atoms` に保存されます。得られた構造からH-Hの距離を確認します。


```python
print("The H-H bond length in H2 molecule is: ")
print(f"    calc: {opt_result.atoms.ase_atoms.get_distance(0,1):.3} Å")
print("    expt: 0.741 Å")
```

    The H-H bond length in H2 molecule is:
        calc: 0.75 Å
        expt: 0.741 Å



```python
from ase.io import Trajectory
from ase.visualize import view

traj = Trajectory("H2_bfgs.traj")
view(traj, viewer="ngl")
```








    HBox(children=(NGLWidget(max_frame=8), VBox(children=(Dropdown(description='Show', options=('All', 'H'), value…





#### JSONファイルへの結果とパラメータの保存と読み込み

シミュレーション結果と呼び出し引数はJSONファイルに保存して、後から読み込むことができます。
なお、いくつかの `Feature` の初期化パラメータは読み込みに対応していません。
具体的には、`estimator_fn` はサポートされていません。
これを `Feature` 初期化時に設定したとしても、 `saved_result.feature` に入っている `Feature` オブジェクトのパラメーターはデフォルト値となっています。以下の例のような保存/再読み込みは、各OptFeaturesの他の結果に対しても同様に行うことができます。


```python
# Save the optimization result.
save_success = opt_result.save("opt_result.json")
print(save_success)
```

    True



```python
# Reload the saved result.
# UserWarning will be shown because estimator_fn does not support saving and reloading.
from matlantis_features.features.common.opt import OptFeatureResult

saved_result = OptFeatureResult.load("opt_result.json")
```

    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/common/opt.py:154: UserWarning: Model version is not specified, default 'v8.0.0' will be used.Please specify it by providing an estimator_fn or setting the 'MATLANTIS_PFP_MODEL_VERSION' environment variable.
      self.check_estimator_fn(estimator_fn=estimator_fn)
    /home/jovyan/.py313/lib/python3.13/site-packages/matlantis_features/features/common/opt.py:154: UserWarning: Calculation mode is not specified, default 'EstimatorCalcMode.PBE' will be used.Please specify it by providing an estimator_fn or setting the 'MATLANTIS_PFP_CALC_MODE' environment variable.
      self.check_estimator_fn(estimator_fn=estimator_fn)



```python
# Inspect loaded results
print("Atoms:", saved_result.atoms)
print("Is converged:", saved_result.converged)

# PFP parameters (model version, calc mode, etc.)
print("Calculator info:", saved_result.calculator_info)

# Package version information
print("Version info:", saved_result.version_info)

# Inspect the feature init parameters (as dict)
print("Feature init params:", saved_result.feature.to_dict())
```

    Atoms: MatlantisAtoms()
    Is converged: True
    Calculator info: {'model_version': 'v8.0.0', 'calc_mode': 'pbe', 'method_type': 'pfvm_d3_pfvm'}
    Version info: {'matlantis-feature-result-version': '0.0.1', 'ase': '3.26.0', 'pfp-api-client': '2.0.1', 'matlantis-features': '1.0.1', 'python': '3.13'}
    Feature init params: {'cls_path': 'matlantis_features.features.common.opt.BFGSASEOptFeature', 'n_run': 1000, 'fmax': 1e-05, 'show_progress_bar': True, 'tqdm_options': None, 'show_logger': False, 'logger_interval': 10, 'filter': None, 'trajectory': 'H2_bfgs.traj', 'maxstep': 0.2}



```python
# Rerunning the simulation
# result_rerun = saved_result.feature(**saved_result.call_params)
```


### LBFGS optimzer
Limited-memory BFGS(LBFGS)法は、BFGS法の計算手順に近似を取り込みメモリ使用量を抑えた手法です。
LBFGSはBFGSとほとんど同じ結果を返しますが、系が大きいほどLBFGSのほうが計算時間が短くなります。
LBFGSのメモリ使用量はO(N)なのに対して、BFGSのメモリ使用量はO(N^2)です。また、ASEの実装ではLBFGSの計算時間はO(N)ですが、BFGSはO(N^3)かかります。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from matlantis_features.features.common.opt import LBFGSASEOptFeature

atoms = Atoms(
    positions=np.array([[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]]),
    numbers=np.array([1,1]),
    pbc=False
)

opt_feature = LBFGSASEOptFeature(
    fmax=0.00001, n_run=1000, show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE),
    trajectory="H2_lbfgs.traj"
)
opt_result = opt_feature(atoms)

print("The H-H bond length in H2 molecule is: ")
print(f"    calc: {opt_result.atoms.ase_atoms.get_distance(0,1):.3} Å")
print("    expt: 0.741 Å")
```


      0%|          | 0/1001 [00:00&lt;?, ?it/s]


    The H-H bond length in H2 molecule is:
        calc: 0.75 Å
        expt: 0.741 Å




### FIRE optimizer

FIREは Fast Inertial Relaxation Engineの略で、分子動力学法(MD)のような方法です。 FIREでは基本的にMDを実行しますが、高速で収束するようにさまざまな工夫が施されています。

[`matlantis-features`](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html)では、FIRE法を使用して構造最適化を実行するために`FireASEOptFeature`が提供されています。 収束基準と最適化ステップの最大数の設定は、上の例と同じです。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from matlantis_features.features.common.opt import FireASEOptFeature

atoms = Atoms(
    positions=np.array([[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]]),
    numbers=np.array([1,1]),
    pbc=False
)

opt_feature = FireASEOptFeature(
    fmax=0.00001, n_run=1000, show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE),
    trajectory="H2_fire.traj"
)
opt_result = opt_feature(atoms)

print("The H-H bond length in H2 molecule is: ")
print(f"    calc: {opt_result.atoms.ase_atoms.get_distance(0,1):.3} Å")
print("    expt: 0.741 Å")
```


      0%|          | 0/1001 [00:00&lt;?, ?it/s]


    The H-H bond length in H2 molecule is:
        calc: 0.75 Å
        expt: 0.741 Å




### FIRE-LBFGS optimizer

上で述べたように、初期構造が安定点に近く二次関数近似が有効な場合、BFGSなどの2次ニュートン法を使用することで良好な収束を実現できます。
二次関数近似が適切ではない不安定な構造の場合、FIRE法の方が頑健にエネルギーを削減できます。

そこで、[`matlantis-features`](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html)では両者の利点を組み合わせた新しい手法を提供しています。この手法は`FIRE-LBFGS`といい、基本となる考え方は次のとおりです。
  - 構造最適化の初期段階でエネルギーとその勾配を減らすためにFIREを使用します。
  - ある程度安定点に達したとき（力が小さい値になったときなど）、LBFGS法を使用して高速に収束状態に向かおうとします。

[`matlantis-features`](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html)ではこのアルゴリズムを`FireLBFGSASEOptFeature`にて提供しています。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from matlantis_features.features.common.opt import FireLBFGSASEOptFeature

atoms = Atoms(
    positions=np.array([[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]]),
    numbers=np.array([1,1]),
    pbc=False
)

opt_feature = FireLBFGSASEOptFeature(
    fmax=0.00001, n_run=1000, show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE),
    trajectory="H2_lbfgsfire.traj"
)
opt_result = opt_feature(atoms)

print("The H-H bond length in H2 molecule is: ")
print(f"    calc: {opt_result.atoms.ase_atoms.get_distance(0,1):.3} Å")
print("    expt: 0.741 Å")
```


      0%|          | 0/1001 [00:00&lt;?, ?it/s]


    The H-H bond length in H2 molecule is:
        calc: 0.75 Å
        expt: 0.741 Å




## 水分子
水分子は折れ線型の構造を持っており、酸素原子を中心に2つの水素原子が酸素原子と結合することで折れ線構造を形成しています。2つのH-O結合のなす角は104.5度であり、H-O結合の結合長は0.958Åです。

水素分子の例と同様に、分子のポテンシャルエネルギーのグラフを描画した後に構造最適化を行います。ここでは簡単のためBFGS法のみを使用します。もし興味があれば他の方法での最適化を試してみてください。


```python
import matlantis_features
import pfp_api_client
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

estimator = Estimator(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
calc = ASECalculator(estimator)

print(f"matlantis_features: {matlantis_features.__version__}")
print(f"pfp_api_client: {pfp_api_client.__version__}")
print(f"current pfp model version: {estimator.model_version}")
```

    matlantis_features: 1.0.1
    pfp_api_client: 2.0.1
    current pfp model version: v8.0.0




### 水分子のポテンシャルエネルギー
水分子のポテンシャルエネルギーを可視化するにあたり、ポテンシャルエネルギーを酸素原子と水素原子の距離と2つのO-H結合のなす角の関数としてグラフを描画します。ここではO-H結合長を0.9Åから1.05Åまで0.01Å間隔で、2つの結合のなす角を90°から120°まで1°間隔で変化させることにします。対称性を考慮し、2つの結合の長さは等しいものとします。分子のポテンシャルエネルギーはこれまでと同様にPFP calculatorで計算します。


```python
import numpy as np
from tqdm.auto import tqdm
from ase import Atoms

def get_water_molecule(bond_length, bond_angle):
    angle = bond_angle / 180 * np.pi
    positions = np.array([
        [0.0, 0.0, 0.0],
        [bond_length, 0.0, 0.0],
        [bond_length * np.cos(angle), bond_length*np.sin(angle), 0.0]
    ])
    return Atoms(positions=positions, numbers=[8,1,1], pbc=False)


bond_lengths = np.round(np.arange(0.90, 1.05, 0.01), 3)
bond_angles = np.arange(90, 120, 1)
energy = np.empty([len(bond_lengths), len(bond_angles)])
for i, r in enumerate(tqdm(bond_lengths)):
    for j, a, in enumerate(tqdm(bond_angles, leave=False)):
        atoms = get_water_molecule(r, a)
        atoms.calc = calc
        energy[i,j] = atoms.get_potential_energy()

```


      0%|          | 0/16 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



      0%|          | 0/30 [00:00&lt;?, ?it/s]



```python
import matplotlib.pyplot as plt

im = plt.imshow(energy, cmap="jet")
plt.colorbar(im, fraction=0.025)
plt.xticks(ticks=np.arange(0, len(bond_angles), 10), labels=bond_angles[::10])
plt.yticks(ticks=np.arange(0, len(bond_lengths), 5), labels=bond_lengths[::5])
plt.xlabel("H-O-H angle")
plt.ylabel("H-O length (Å)")
```




    Text(0, 0.5, 'H-O length (Å)')





![png](output_28_1.png)




```python
ind = np.unravel_index(energy.argmin(), energy.shape)
print(f"The potential energy of H-O-H is minimized when H-O length is {bond_lengths[ind[0]]:.3f} Å and H-O-H angle is {bond_angles[ind[1]]:.1f}")
```

    The potential energy of H-O-H is minimized when H-O length is 0.970 Å and H-O-H angle is 105.0




### BFGS optimizer

BFGS法を用いて水分子の構造最適化を行います。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from matlantis_features.features.common.opt import BFGSASEOptFeature

atoms = get_water_molecule(bond_length=1.5, bond_angle=150.0)

opt_feature = BFGSASEOptFeature(
    fmax=0.00001, n_run=1000, show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE),
    trajectory="H2O_bfgs.traj"
)
opt_result = opt_feature(atoms)
```


      0%|          | 0/1001 [00:00&lt;?, ?it/s]




最適化によって得られた構造のH-O結合長と2つの結合のなす角を確認します。


```python
print("The H-O bond length in H2O molecule is: ")
print(f"    calc: {opt_result.atoms.ase_atoms.get_distance(0, 1):.3} Å")
print("    expt: 0.958 Å")
```

    The H-O bond length in H2O molecule is:
        calc: 0.972 Å
        expt: 0.958 Å



```python
print("The H-O-H bond angle in H2O molecule is: ")
print(f"    calc: {opt_result.atoms.ase_atoms.get_angle(1, 0, 2):.1f}")
print("    expt: 104.5")
```

    The H-O-H bond angle in H2O molecule is:
        calc: 104.5
        expt: 104.5



```python
from ase.io import Trajectory
from ase.visualize import view

traj = Trajectory("H2O_bfgs.traj")
view(traj, viewer="ngl")
```




    HBox(children=(NGLWidget(max_frame=11), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H'),…





## シリコン結晶
結晶のポテンシャルエネルギーは分子と違って、格子内の原子の相対的位置のみならず、格子の大きさや形状によっても決められます。以下の例では、Si結晶の格子定数によるポテンシャルエネルギーの変化を計算し、[`matlantis-features`](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html) を用いたセル形状の最適化の方法を示します。


```python
import matlantis_features
import pfp_api_client
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator

estimator = Estimator(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
calc = ASECalculator(estimator)

print(f"matlantis_features: {matlantis_features.__version__}")
print(f"pfp_api_client: {pfp_api_client.__version__}")
print(f"current pfp model version: {estimator.model_version}")
```

    matlantis_features: 1.0.1
    pfp_api_client: 2.0.1
    current pfp model version: v8.0.0




### 格子定数とポテンシャルエネルギーの関係

格子定数とポテンシャルエネルギーの関係性を可視化するために、primitive cell内の2個のSi原子の分率座標(scaled positions, fractional coordinates)を固定し、格子定数のみ変化させます。構造のポテンシャルエネルギーは、PFP calculator によって計算されます。


```python
import numpy as np
from ase.build import bulk
from tqdm.auto import tqdm
```


```python
atoms = bulk("Si")

lattice_constant = np.arange(5.4, 5.6, 0.01)

energy = []
for a in tqdm(lattice_constant):
    cell = [[0.0, a/2, a/2], [a/2, 0.0, a/2], [a/2, a/2, 0.0]]
    atoms.set_cell(cell, scale_atoms=True)
    atoms.calc = calc
    energy.append(atoms.get_potential_energy())

```


      0%|          | 0/20 [00:00&lt;?, ?it/s]



```python
import matplotlib.pyplot as plt

plt.plot(lattice_constant, energy)
plt.xlabel("lattice constant (Å)")
plt.ylabel("potential energy (eV)")
```




    Text(0, 0.5, 'potential energy (eV)')





![png](output_41_1.png)




```python
print(f"The energy of Si is minimized with the lattice constant of {lattice_constant[np.argmin(energy)]:.2f} Å")
```

    The energy of Si is minimized with the lattice constant of 5.47 Å




### フィルター付きのセル形状の最適化

`UnitCellASEFilter`を使用してフィルター付き最適化を実行できます。 フィルタを使用することにより、原子にかかる力とユニットセルにかかる応力を同時に扱い、最小化できます。
この例では、シリコン結晶を使用します。 最適化を実行する前に、格子定数を1.1倍に伸ばしてみましょう。


```python
atoms = bulk("Si")
cell = atoms.get_cell()[:] * 1.1
atoms.set_cell(cell, scale_atoms=True)
```


```python
from matlantis_features.filters import UnitCellASEFilter
from matlantis_features.features.common import BFGSASEOptFeature

cell_filter = UnitCellASEFilter()
opt_feature = BFGSASEOptFeature(
    fmax=0.01,
    filter=cell_filter,
    show_logger=True,
    show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
)
opt_result = opt_feature(atoms)
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]



```python
print("Cell shape before optimization: ")
print(atoms.get_cell()[:])

cell = opt_result.atoms.ase_atoms.get_cell()[:]
print("Cell shape after optimization: ")
print(np.round(cell, 3))
print("Lattice constant of crystal silicon: ")
print(f"    calc: {cell[0,1]*2:.3f} Å")
print("    expt: 5.431 Å")

```

    Cell shape before optimization:
    [[0.     2.9865 2.9865]
     [2.9865 0.     2.9865]
     [2.9865 2.9865 0.    ]]
    Cell shape after optimization:
    [[-0.     2.733  2.733]
     [ 2.733 -0.     2.733]
     [ 2.733  2.733 -0.   ]]
    Lattice constant of crystal silicon:
        calc: 5.465 Å
        expt: 5.431 Å




また、フィルターオブジェクトを初期化せず、手軽にセル形状の最適化を実行する方法も提供しています。Optimizerのオプションで`filter=True`と設定するだけで、内部で`UnitCellASEFilter`が呼び出されます。


```python
opt_feature = BFGSASEOptFeature(
    fmax=0.01,
    filter=True,
    show_logger=True,
    show_progress_bar=True,
    estimator_fn=pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
)
opt_result = opt_feature(atoms)
```


      0%|          | 0/201 [00:00&lt;?, ?it/s]




最適化前後のセル形状を比較すると、想定通り格子が縮小していることがわかります。


```python
print("Cell shape before optimization: ")
print(atoms.get_cell()[:])

cell = opt_result.atoms.ase_atoms.get_cell()[:]
print("Cell shape after optimization: ")
print(np.round(cell, 3))
```

    Cell shape before optimization:
    [[0.     2.9865 2.9865]
     [2.9865 0.     2.9865]
     [2.9865 2.9865 0.    ]]
    Cell shape after optimization:
    [[-0.     2.733  2.733]
     [ 2.733 -0.     2.733]
     [ 2.733  2.733 -0.   ]]

