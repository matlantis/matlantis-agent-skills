# Matlantis-features: MDの発展的な使い方
**Contents**

* [0. 基本的な使い方](#advanced_0)
* [1. 乱数シードの固定方法](#advanced_1)
* [2. MD extensionsの使い方](#advanced_2)
* [3. 独自に定義したMD extensionの使い方](#advanced_3)
* [4. NPT integratorの使い方](#advanced_4)
* [5. check point fileからの計算再開方法](#advanced_5)
* [6. MD integratorを構築する方法](#advanced_6)


&lt;a id='advanced_0'&gt;&lt;/a&gt;
## 0. 基本的な使い方
### 初期化
* 必要なライブラリを準備します。


```python
# If you have already installed matlantis-features, you can skip this.
# Version 0.8.1 or later is required to run this notebook.
!pip install 'matlantis-features&gt;=0.8.1'
```




```python
import matlantis_features
import pfp_api_client
from pfp_api_client.pfp.estimator import Estimator

estimator = Estimator()
print(f"matlantis_features: {matlantis_features.__version__}")
print(f"pfp_api_client: {pfp_api_client.__version__}")
print(f"current pfp model version: {estimator.model_version}")
```

    matlantis_features: 1.0.1
    pfp_api_client: 2.0.1
    current pfp model version: v8.0.0



### `estimator_fn` による計算モードとモデルバージョンの指定

* Feature に `estimator_fn` 引数を使うことで、matlantis-features で使われる PFP の計算モードとモデルバージョンを指定することができます。

* `estimator_fn` は、Estimator オブジェクトを作る factory method です。計算モードとモデルバージョンのみを指定する場合は、以下のように `pfp_estimator_fn` を利用できます。より細かい指定が必要な場合は、自身で factory method を定義できます。詳細に関しては [Matlantis Guidebook](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html#estimator-fnpfp) をご参照ください。

* `estimator_fn` が指定されない場合は、環境変数で指定された値が使用されます。`estimator_fn` も環境変数も指定されない場合は、デフォルトのモデルバージョンと計算モード（`PBE`）が使用されます。

* 環境変数と `estimator_fn` が同時に定義されている場合、 `estimator_fn` の設定が優先して使用されます。


```python
from matlantis_features.utils.calculators import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode

estimator_fn = pfp_estimator_fn(model_version='v8.0.0', calc_mode=EstimatorCalcMode.PBE)
```


### 基本的な構成要素
* MD計算を始めるには、次のクラスを用意します。
    * MD System
    * Integrator
    * MD Feature
    * MD Extension (optional)

* 以下、それぞれのクラスについて紹介します。


### MD System

* MD systemは、原子位置や時刻などシミュレーション対象の基本的な情報を保持します。
* `ASEMDSystem`は、ASEの`Atoms`と`MatlantisAtoms`の両方から初期化できます。


```python
from ase.build import bulk
from matlantis_features.features.md import ASEMDSystem

system = ASEMDSystem(bulk("Si", cubic=True), step = 0, time = 0.0)
```


* 原子の座標は`ASEMDSystem.ase_atoms`に格納されています。


```python
print(system.ase_atoms)
```

    Atoms(symbols='Si8', pbc=True, cell=[5.43, 5.43, 5.43])



* MD systemは現在のステップ数と時刻(fs)も記録しています。


```python
print(f"current MD step: {system.current_step}")
print(f"current MD time: {system.current_time} fs")
```

    current MD step: 0
    current MD time: 0.0 fs



* 現在のMD systemのすべての情報を`system.state`から確認することができます。
* `system.state`をcheckpointファイルに保存し、計算再開時に使用することができます。


```python
print(system.state)
```

    {'positions': array([[0.    , 0.    , 0.    ],
           [1.3575, 1.3575, 1.3575],
           [0.    , 2.715 , 2.715 ],
           [1.3575, 4.0725, 4.0725],
           [2.715 , 0.    , 2.715 ],
           [4.0725, 1.3575, 4.0725],
           [2.715 , 2.715 , 0.    ],
           [4.0725, 4.0725, 1.3575]]), 'velocities': array([[0., 0., 0.],
           [0., 0., 0.],
           [0., 0., 0.],
           [0., 0., 0.],
           [0., 0., 0.],
           [0., 0., 0.],
           [0., 0., 0.],
           [0., 0., 0.]]), 'cell': Cell([5.43, 5.43, 5.43]), 'pbc': array([ True,  True,  True]), 'total_time': 0.0, 'total_step': 0}



* MD systemの初期化時には各原子の速度は0になっています。 `init_temperature()`メソッドにより、ボルツマン分布に従って各原子の速度を初期化することができます。


```python
system.init_temperature(temperature=300, stationary=True, zero_rotation=True)
print(system.state["velocities"])
```

    [[ 0.00389069 -0.03350972  0.06374191]
     [ 0.01898614 -0.01518247 -0.06587098]
     [ 0.02174528 -0.01364948  0.01273067]
     [ 0.02710742 -0.03749348 -0.03745951]
     [-0.07987663  0.03000984 -0.05202031]
     [-0.00100243  0.01243713  0.04579406]
     [-0.03516471  0.01140989  0.00391086]
     [ 0.04431423  0.0459783   0.02917329]]



### Integrator
* integratorにはMDのtrajectoryを発展させる際に使用するアルゴリズムが格納されます(アンサンブル、熱浴、圧力浴など)。
* 現在、Matlantis-featuresでは以下のintegratorが提供されています。
    * VelocityVerletIntegrator
    * LangevinIntegrator
    * AndersenIntegrator
    * NVTBerendsenIntegrator
    * NPTBerendsenIntegrator
    * NPTIntegrator


* 例として、以下では`NPTBerendsenIntegrator`を使用します。


```python
from ase import units
from matlantis_features.features.md import NPTBerendsenIntegrator

integrator = NPTBerendsenIntegrator(
    timestep = 1.0,
    temperature = 300.0,
    pressure = 101325 * units.Pascal,
)
```


* `NPTBerendsenIntegrator`はASEの`NPTBerendsen`クラスのラッパーです。 `create_ase_dynamics`メソッドによりASEのクラスに変換することもできます。


```python
dyn = integrator.create_ase_dynamics(system.ase_atoms)
print(dyn)
```

    &lt;ase.md.nptberendsen.NPTBerendsen object at 0x7f9d68bf4590&gt;



### MD Feature

* MD systemとintegratorの準備ができたら、`MDFeature`クラスによりMDを実行できます.
* 以下の例では、MDを100ステップ行い、10ステップごとにsystemのスナップショットをエネルギーや力の値とともにtrajectryファイルに保存しています。
* MDの進捗状況をプログレスバーとログにより可視化することができます。


```python
import logging
from matlantis_features.features.md import MDFeature

logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

md = MDFeature(
    integrator=integrator,
    n_run=100,
    show_progress_bar=True,
    show_logger=True,
    logger_interval=10,
    traj_file_name="md.traj",
    traj_freq=10,
    traj_props=["energy", "forces"],
    estimator_fn=estimator_fn,
)
md_results = md(system)
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/md.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.49 eV/atom  temperature: 448.76 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:    10  energy：-4.53 eV/atom  total energy: -4.49 eV/atom  temperature: 273.28 K  volume:   161 Ang^3  density: 2.313 g/cm^3
    steps:    20  energy：-4.51 eV/atom  total energy: -4.49 eV/atom  temperature: 122.45 K  volume:   162 Ang^3  density: 2.300 g/cm^3
    steps:    30  energy：-4.53 eV/atom  total energy: -4.49 eV/atom  temperature: 257.67 K  volume:   163 Ang^3  density: 2.294 g/cm^3
    steps:    40  energy：-4.53 eV/atom  total energy: -4.49 eV/atom  temperature: 301.99 K  volume:   163 Ang^3  density: 2.293 g/cm^3
    steps:    50  energy：-4.51 eV/atom  total energy: -4.49 eV/atom  temperature: 111.54 K  volume:   163 Ang^3  density: 2.293 g/cm^3
    steps:    60  energy：-4.50 eV/atom  total energy: -4.49 eV/atom  temperature: 52.85 K  volume:   163 Ang^3  density: 2.293 g/cm^3
    steps:    70  energy：-4.52 eV/atom  total energy: -4.49 eV/atom  temperature: 211.90 K  volume:   163 Ang^3  density: 2.293 g/cm^3
    steps:    80  energy：-4.53 eV/atom  total energy: -4.49 eV/atom  temperature: 274.75 K  volume:   163 Ang^3  density: 2.293 g/cm^3
    steps:    90  energy：-4.52 eV/atom  total energy: -4.49 eV/atom  temperature: 213.95 K  volume:   163 Ang^3  density: 2.291 g/cm^3
    steps:   100  energy：-4.52 eV/atom  total energy: -4.49 eV/atom  temperature: 240.81 K  volume:   163 Ang^3  density: 2.287 g/cm^3


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-4.55 eV/atom  total energy: -4.50 eV/atom  temperature: 361.07 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:    10  energy：-4.53 eV/atom  total energy: -4.51 eV/atom  temperature: 210.93 K  volume:   161 Ang^3  density: 2.313 g/cm^3


    steps:    20  energy：-4.52 eV/atom  total energy: -4.51 eV/atom  temperature: 91.56 K  volume:   162 Ang^3  density: 2.301 g/cm^3


    steps:    30  energy：-4.53 eV/atom  total energy: -4.51 eV/atom  temperature: 218.89 K  volume:   163 Ang^3  density: 2.295 g/cm^3


    steps:    40  energy：-4.54 eV/atom  total energy: -4.51 eV/atom  temperature: 243.02 K  volume:   163 Ang^3  density: 2.293 g/cm^3


    steps:    50  energy：-4.52 eV/atom  total energy: -4.50 eV/atom  temperature: 102.41 K  volume:   163 Ang^3  density: 2.292 g/cm^3


    steps:    60  energy：-4.51 eV/atom  total energy: -4.50 eV/atom  temperature: 58.63 K  volume:   163 Ang^3  density: 2.291 g/cm^3


    steps:    70  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 177.77 K  volume:   163 Ang^3  density: 2.290 g/cm^3


    steps:    80  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 221.66 K  volume:   163 Ang^3  density: 2.290 g/cm^3


    steps:    90  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 180.67 K  volume:   163 Ang^3  density: 2.289 g/cm^3


    steps:   100  energy：-4.53 eV/atom  total energy: -4.50 eV/atom  temperature: 205.47 K  volume:   163 Ang^3  density: 2.286 g/cm^3



* MD完了後、MD systemの状態を出力してみます。原子位置が変化し、`total_time`が100fs、`total_step`が100になっていることがわかります。


```python
system.state
```




    {'positions': array([[ 0.03331952, -0.04873971,  0.07135978],
            [ 1.36739978,  1.35698283,  1.27844893],
            [ 0.05149021,  2.70382806,  2.73958557],
            [ 1.37447682,  4.06840978,  4.05131959],
            [ 2.6548874 ,  0.04099131,  2.66489255],
            [ 4.08106547,  1.36761821,  4.1648284 ],
            [ 2.71325241,  2.76013193, -0.00463353],
            [ 4.11740286,  4.14407205,  1.42749319]]),
     'velocities': array([[ 0.04585949,  0.03543682, -0.01010219],
            [-0.01653031, -0.00562304,  0.03595251],
            [ 0.02325738,  0.02577229, -0.02920272],
            [-0.02558901,  0.03021675,  0.03482884],
            [ 0.0327851 , -0.04055633,  0.01076809],
            [-0.00189992, -0.00470169, -0.03007233],
            [-0.01732267, -0.01911698,  0.01867509],
            [-0.04056006, -0.02142782, -0.03084729]]),
     'cell': Cell([5.464431486118609, 5.464431486118609, 5.464431486118609]),
     'pbc': array([ True,  True,  True]),
     'total_time': 100.00000000000001,
     'total_step': 100}




* ASEのtrajectoryと同様に、保存されたtrajectoryファイルからスナップショットを読み込むことができます。


```python
from ase.io import Trajectory

traj = Trajectory("md.traj")
atoms = traj[5]
```


### MD Extension
* MD extensionを用いることで、MD実行中に独自の関数を呼び出すことができます(たとえば、$n$ステップ毎に内部情報を標準出力させたり、変更したりするなど)。
* MD extensionの詳細については[2. MD extensionsの使い方](#advanced_2)、[3. 独自に定義したMD extensionの使い方](#advanced_3)にて紹介しています。


&lt;a id='advanced_1'&gt;&lt;/a&gt;
## 1. 乱数シードの固定方法
* 原子の初期速度は乱数を用いて生成されます。そのため、MDのtrajectoryは乱数シードを固定しない限り再現性がありません。
* 以下では乱数シードを固定する手順を紹介します。


```python
import numpy as np
from ase.build import bulk
from matlantis_features.features.md import MDFeature, ASEMDSystem, VelocityVerletIntegrator
```


* MDシミュレーションを100ステップ、温度1000.0Kで実行してみます。
* `MDSystem`の`init_temperature`メソッドのパラメータ`rng`を設定することで乱数シードを固定することができます。


```python
integrator = VelocityVerletIntegrator(timestep = 1.0)

system = ASEMDSystem(bulk("Si", cubic=True))
# system.init_temperature(1000.0, stationary=True)
system.init_temperature(1000.0, stationary=True, rng=np.random.RandomState(999))

md = MDFeature(integrator=integrator, n_run=100, show_progress_bar=True, estimator_fn=estimator_fn)
md_results = md(system)
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /tmp/matlantis_wsc07e_b/tmpy8jvkjkh.traj.
    Note: The max disk size of / is about 30G.
    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.


    Note: The max disk size of / is about 30G.


    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.



* 乱数シードを固定したので、異なるMD試行においても初期速度は同じになります。これにより、100ステップ後の原子配置は毎回同じような値になります。ただし、GPUの数値誤差により完全に一致するとは限りません。
* 正しく乱数シードを固定できていれば、100MDステップ後の原子座標は以下の値に近くなっているはずです。 &lt;br&gt;
\[-0.02566377  0.08580367  0.00859627] &lt;br&gt;
 [ 1.29153049  1.38985881  1.29309867] &lt;br&gt;
 [ 0.1078383   2.59242697  2.82837274] &lt;br&gt;
 [ 1.41601203  4.07919393  4.10945311] &lt;br&gt;
 [ 2.62866541  0.01489882  2.65615868] &lt;br&gt;
 [ 4.1075896   1.38416378  4.16920656] &lt;br&gt;
 [ 2.73973303  2.67174474 -0.03008335] &lt;br&gt;
 [ 4.02429492  4.07190928  1.25519733] &lt;br&gt;


```python
print(system.ase_atoms.get_positions())
```

    [[-0.02840698  0.08652461  0.01064221]
     [ 1.29081265  1.39367452  1.29087199]
     [ 0.11163018  2.59317641  2.83485908]
     [ 1.41408443  4.08539018  4.11547125]
     [ 2.627069    0.01505194  2.6467917 ]
     [ 4.10434162  1.37895897  4.16472741]
     [ 2.74285423  2.67525255 -0.03242459]
     [ 4.02761487  4.06197082  1.25906096]]



* 上記のセルを複数回実行し、毎回結果が同じような値になることが確認してみてください。
* 逆に、速度を初期化する際に`rng`を指定しないと、MDの結果毎回異なる原子配置が得られることも確認してみてください。
* 数値誤差はMDの進行に伴い累積していくことに注意してください。MDステップ数が非常に大きい場合、乱数シードを固定したとしても最終的な原子配置が同じになることは保証できません。


&lt;a id='advanced_2'&gt;&lt;/a&gt;
## 2. MD extensionsの使い方

* MD extensionsを用いて、MDの途中で系の状態を変更する方法を紹介します。
* ここでは、MatlantisFeaturesで提供されている`TemperatureScheduler`を例として用います。
* `TemperatureScheduler`は系を徐々に加熱したり、冷却したりすることができます。液体急冷法などのシミュレーションを行いたいときに有効です。


```python
import logging
import numpy as np
from ase import units
from ase.build import bulk

from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.md import NVTBerendsenIntegrator, ASEMDSystem, MDFeature
from matlantis_features.features.md.md_extensions import TemperatureScheduler
```


* `TemperatureScheduler`を初期化します。 次の3つのパラメータを指定する必要があります:
    * `start_value`: 初期温度
    * `end_value`: 最終温度
    * `num_total_steps`: 温度変化を起こす合計ステップ数
* `extensions`パラメータに、`TemperatureScheduler`および「extensionを何回実行するか」を指定する整数値を渡します。
* ログにて系の温度が徐々に減少していることが確認できます。


```python
logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

integrator = NVTBerendsenIntegrator(
    timestep=1.0,
    temperature=2000,
    taut=10.0
    )

system = ASEMDSystem(bulk("Si", cubic=True))
system.init_temperature(2000.0, stationary=True, rng=np.random.RandomState(seed=999))

temp_schedular = TemperatureScheduler(start_value=2000.0, end_value=1000.0, num_total_steps=1000)

md = MDFeature(integrator=integrator, n_run=1000, show_progress_bar=True, show_logger=True, logger_interval=100, estimator_fn=estimator_fn)
md_results = md(system, extensions=[(temp_schedular, 1)])
```


      0%|          | 0/1000 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /tmp/matlantis_fo_f7zik/tmpdh0i7igs.traj.
    Note: The max disk size of / is about 30G.
    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.30 eV/atom  temperature: 1964.30 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   100  energy：-4.10 eV/atom  total energy: -3.86 eV/atom  temperature: 1832.74 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   200  energy：-4.32 eV/atom  total energy: -4.13 eV/atom  temperature: 1485.25 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   300  energy：-4.27 eV/atom  total energy: -4.04 eV/atom  temperature: 1784.77 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   400  energy：-4.28 eV/atom  total energy: -4.09 eV/atom  temperature: 1427.28 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   500  energy：-4.31 eV/atom  total energy: -4.11 eV/atom  temperature: 1527.61 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   600  energy：-4.32 eV/atom  total energy: -4.16 eV/atom  temperature: 1233.33 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   700  energy：-4.38 eV/atom  total energy: -4.20 eV/atom  temperature: 1423.19 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   800  energy：-4.36 eV/atom  total energy: -4.20 eV/atom  temperature: 1298.44 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:   900  energy：-4.41 eV/atom  total energy: -4.25 eV/atom  temperature: 1188.66 K  volume:   160 Ang^3  density: 2.330 g/cm^3
    steps:  1000  energy：-4.42 eV/atom  total energy: -4.28 eV/atom  temperature: 1079.48 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    Note: The max disk size of / is about 30G.


    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.


    steps:     0  energy：-4.55 eV/atom  total energy: -4.30 eV/atom  temperature: 1964.30 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   100  energy：-4.10 eV/atom  total energy: -3.86 eV/atom  temperature: 1832.74 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   200  energy：-4.32 eV/atom  total energy: -4.13 eV/atom  temperature: 1485.25 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   300  energy：-4.27 eV/atom  total energy: -4.04 eV/atom  temperature: 1784.77 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   400  energy：-4.28 eV/atom  total energy: -4.09 eV/atom  temperature: 1427.28 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   500  energy：-4.31 eV/atom  total energy: -4.11 eV/atom  temperature: 1527.61 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   600  energy：-4.32 eV/atom  total energy: -4.16 eV/atom  temperature: 1233.33 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   700  energy：-4.38 eV/atom  total energy: -4.20 eV/atom  temperature: 1423.20 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   800  energy：-4.36 eV/atom  total energy: -4.20 eV/atom  temperature: 1298.44 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:   900  energy：-4.41 eV/atom  total energy: -4.25 eV/atom  temperature: 1188.66 K  volume:   160 Ang^3  density: 2.330 g/cm^3


    steps:  1000  energy：-4.42 eV/atom  total energy: -4.28 eV/atom  temperature: 1079.48 K  volume:   160 Ang^3  density: 2.330 g/cm^3



&lt;a id='advanced_3'&gt;&lt;/a&gt;
## 3. 独自に定義したMD extensionの使い方
* 内部状態を抽出や変更するため、独自のextensionを定義したくなることがあるかもしれません。ここでは、独自定義されたMD extensionの使い方を紹介します。


```python
import logging
import numpy as np
from ase import units
from ase.build import bulk

from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.md import NVTBerendsenIntegrator, ASEMDSystem, MDFeature, MDExtensionBase
```


* ここでは、内部状態(ステップ数、温度、全エネルギー、ポテンシャルエネルギー、密度)を出力する関数を定義してみます。
* 新しいextensionを定義するには、MDExtensionBaseを継承している必要があります。
* このextensionに`__call__`メソッドを定義し、入力として`MDSystem`と`MDIntegrator`を受け取るようにします。
* `__call__`メソッド内にてMD systemからASE atomsを受け取り、様々な内部状態を計算します。


```python
class PrintInfo(MDExtensionBase):
    def __call__(self, system, integrator) -&gt; None:
        atoms = system.ase_atoms
        istep = system.current_total_step
        temperature = atoms.get_temperature()
        total_energy = atoms.get_total_energy()
        potential_energy = atoms.get_potential_energy()
        density = atoms.get_masses().sum() / units.kg / atoms.get_volume() * 1e27
        print(f"Dyn step {istep:4d} temperature {temperature:3.2f} total energy {total_energy:3.2f} potential energy {potential_energy:3.2f} density {density:3.2f}")
```


* 上の`PrintInfo` extension を使ってみましょう。
* 10MDステップごとに情報が出力されていることが確認できます。


```python
logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

integrator = NVTBerendsenIntegrator(
    timestep=1.0,
    temperature=2000,
    taut=10.0
    )

system = ASEMDSystem(bulk("Si", cubic=True))
system.init_temperature(2000.0, stationary=True, rng=np.random.RandomState(seed=999))

md = MDFeature(integrator=integrator, n_run=100, show_progress_bar=True, show_logger=False, logger_interval=100, estimator_fn=estimator_fn)
md_results = md(system, extensions=[(PrintInfo(), 10)])
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /tmp/matlantis_ftmrhpnt/tmp8qzx7eb_.traj.
    Note: The max disk size of / is about 30G.
    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.


    Dyn step    0 temperature 1964.30 total energy -34.37 potential energy -36.40 density 2.33
    Dyn step   10 temperature 1253.39 total energy -34.08 potential energy -35.38 density 2.33
    Dyn step   20 temperature 862.98 total energy -33.02 potential energy -33.91 density 2.33
    Dyn step   30 temperature 1673.68 total energy -32.18 potential energy -33.91 density 2.33
    Dyn step   40 temperature 1867.15 total energy -32.02 potential energy -33.95 density 2.33
    Dyn step   50 temperature 1595.56 total energy -31.76 potential energy -33.41 density 2.33
    Dyn step   60 temperature 1590.52 total energy -31.27 potential energy -32.92 density 2.33
    Dyn step   70 temperature 2355.80 total energy -31.21 potential energy -33.64 density 2.33
    Dyn step   80 temperature 2108.33 total energy -31.57 potential energy -33.75 density 2.33
    Dyn step   90 temperature 1323.62 total energy -31.29 potential energy -32.66 density 2.33
    Dyn step  100 temperature 1843.60 total energy -30.69 potential energy -32.59 density 2.33


    Dyn step   10 temperature 1253.39 total energy -34.08 potential energy -35.38 density 2.33


    Dyn step   20 temperature 862.98 total energy -33.02 potential energy -33.91 density 2.33


    Dyn step   30 temperature 1673.68 total energy -32.18 potential energy -33.91 density 2.33


    Dyn step   40 temperature 1867.15 total energy -32.02 potential energy -33.95 density 2.33


    Dyn step   50 temperature 1595.56 total energy -31.76 potential energy -33.41 density 2.33


    Dyn step   60 temperature 1590.52 total energy -31.27 potential energy -32.92 density 2.33


    Dyn step   70 temperature 2355.80 total energy -31.21 potential energy -33.64 density 2.33


    Dyn step   80 temperature 2108.33 total energy -31.57 potential energy -33.75 density 2.33


    Dyn step   90 temperature 1323.62 total energy -31.29 potential energy -32.66 density 2.33


    Dyn step  100 temperature 1843.60 total energy -30.69 potential energy -32.59 density 2.33



&lt;a id='advanced_4'&gt;&lt;/a&gt;
## 4. NPT integratorの使い方
* この章では、NPT integratorの適切な使い方を紹介します。
* NPT integratorは`ase.md.NPT`クラスのラッパーです。 Nose-Hoover法およびParrinello-Rahman法による熱浴、圧力浴を提供しています。
* 次の事例を紹介します:
    * 4.1 NPT-MD中に格子を直交に保つ方法
    * 4.2 NPT-MD中に格子定数の比a:b:cを保つ方法
    * 4.3 三斜晶系においてセル行列を上三角行列に変換する方法


```python
import logging
import pathlib
import numpy as np
import ase.io
from ase import units
from ase.build import bulk

from matlantis_features.atoms import MatlantisAtoms
from matlantis_features.features.md import NPTIntegrator, ASEMDSystem, MDFeature, MDExtensionBase
try:
    dir_path = pathlib.Path(__file__).parent
except:
    dir_path = pathlib.Path("").resolve()

```


### 4.1 NPT-MD中に格子を直交に保つ方法
* NPT-MDシミュレーション中に格子を直交に保ちたい場合、NPT integratorの`mask`パラメータを設定します。



* 以下の例では、ベンチマークとしてクリストバライト型の二酸化ケイ素を用います。SiO2 (mp-6945) 構造は [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) の下でMaterials Projectからダウンロード可能です。 (https://materialsproject.org/materials/mp-6945/)
* クリストバライト型の二酸化ケイ素の空間群は、正方晶系である$P4_12_12$です。


```python
atoms = ase.io.read(str(dir_path/"assets/md_adv_usage/SiO2_mp-6945_computed.cif")) * (2,2,2)
```


* まず、格子定数(基本格子ベクトル$\vec a$, $\vec b$, $\vec c$の長さおよび軸どうしの角度)を出力するextensionを定義します。


```python
class PrintCellShape(MDExtensionBase):
    def __init__(self, cell_log=None):
        self.cell_log = cell_log
    def __call__(self, system, integrator) -&gt; None:
        cell_par = system.ase_atoms.cell.cellpar()
        istep = system.current_total_step
        print(f"Dyn step {istep:4d} a {cell_par[0]:3.2f} b {cell_par[1]:3.2f} c {cell_par[2]:3.2f} alpha {cell_par[3]:3.2f} beta {cell_par[4]:3.2f} gamma {cell_par[5]:3.2f} ")
        if self.cell_log is not None:
            self.cell_log.append(cell_par)
```


* 次に、NPT integratorの`mask`パラメータを`np.eye(3)`に設定します。これは、stress tensorの非対角項を無視することに相当します。したがって、基本格子ベクトルの長さのみが可変となり、軸どうしの角度は$90^{\circ}$に固定されます。
* 出力から角度が$90^{\circ}$度に固定されていることがわかります。


```python
integrator = NPTIntegrator(
    timestep=1.0,
    temperature=300,
    pressure=101325 * units.Pascal,
    ttime=10,
    pfactor=100 * units.fs,
    mask=np.eye(3)
)

system = ASEMDSystem(atoms.copy())
system.init_temperature(300.0, stationary=True, rng=np.random.RandomState(seed=999))

md = MDFeature(integrator=integrator, n_run=200, show_progress_bar=True, show_logger=True, logger_interval=100, estimator_fn=estimator_fn)
md_results = md(system, extensions=[(PrintCellShape(), 10)])
```


      0%|          | 0/200 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /tmp/matlantis_vbv6nnkx/tmpruaqyg4i.traj.
    Note: The max disk size of / is about 30G.
    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.
    steps:     0  energy：-6.35 eV/atom  total energy: -6.31 eV/atom  temperature: 278.01 K  volume:  1468 Ang^3  density: 2.174 g/cm^3
    steps:     0  energy：-6.35 eV/atom  total energy: -6.31 eV/atom  temperature: 278.01 K  volume:  1468 Ang^3  density: 2.174 g/cm^3


    Dyn step    0 a 10.17 b 10.17 c 14.20 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step    0 a 10.17 b 10.17 c 14.20 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   10 a 10.17 b 10.17 c 14.20 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   20 a 10.17 b 10.17 c 14.21 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   30 a 10.18 b 10.18 c 14.22 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   40 a 10.19 b 10.19 c 14.24 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   50 a 10.19 b 10.19 c 14.26 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   60 a 10.19 b 10.18 c 14.26 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   70 a 10.17 b 10.17 c 14.25 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   80 a 10.16 b 10.16 c 14.24 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   90 a 10.15 b 10.15 c 14.25 alpha 90.00 beta 90.00 gamma 90.00


    steps:   100  energy：-6.31 eV/atom  total energy: -6.27 eV/atom  temperature: 267.72 K  volume:  1473 Ang^3  density: 2.168 g/cm^3


    Dyn step  100 a 10.16 b 10.16 c 14.26 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  110 a 10.18 b 10.19 c 14.30 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  120 a 10.22 b 10.22 c 14.36 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  130 a 10.25 b 10.26 c 14.41 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  140 a 10.29 b 10.29 c 14.45 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  150 a 10.32 b 10.31 c 14.49 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  160 a 10.35 b 10.32 c 14.51 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  170 a 10.36 b 10.33 c 14.52 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  180 a 10.38 b 10.33 c 14.54 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step  190 a 10.40 b 10.32 c 14.55 alpha 90.00 beta 90.00 gamma 90.00


    steps:   200  energy：-6.32 eV/atom  total energy: -6.28 eV/atom  temperature: 264.49 K  volume:  1567 Ang^3  density: 2.037 g/cm^3


    Dyn step  200 a 10.42 b 10.33 c 14.57 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   90 a 10.15 b 10.15 c 14.25 alpha 90.00 beta 90.00 gamma 90.00


    steps:   100  energy：-6.31 eV/atom  total energy: -6.27 eV/atom  temperature: 267.72 K  volume:  1473 Ang^3  density: 2.168 g/cm^3


    Dyn step  100 a 10.16 b 10.16 c 14.26 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  110 a 10.18 b 10.19 c 14.30 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  120 a 10.22 b 10.22 c 14.36 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  130 a 10.25 b 10.26 c 14.41 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  140 a 10.29 b 10.29 c 14.45 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  150 a 10.32 b 10.31 c 14.49 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  160 a 10.35 b 10.32 c 14.51 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  170 a 10.36 b 10.33 c 14.52 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  180 a 10.38 b 10.33 c 14.54 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step  190 a 10.40 b 10.32 c 14.55 alpha 90.00 beta 90.00 gamma 90.00


    steps:   200  energy：-6.32 eV/atom  total energy: -6.28 eV/atom  temperature: 264.49 K  volume:  1567 Ang^3  density: 2.037 g/cm^3


    Dyn step  200 a 10.42 b 10.33 c 14.57 alpha 90.00 beta 90.00 gamma 90.00



### 4.2 NPT-MD中に格子定数の比a:b:cを保つ方法
* 上の例では、格子を直交に保つ方法を示しました。しかし、MDの途中で$\left | \vec a \right |:\left | \vec b \right |:\left | \vec c \right |$が変更さていることがわかります。格子の拡張が等方的ではないということになります。
* 格子が等方的であるということが知られている物質、たとえばSi結晶などでは、`hydrostatic_strain` パラメータを設定することで$\left | \vec a \right |:\left | \vec b \right |:\left | \vec c \right |$を固定する制約を加えることができます。


* まず、Si結晶の構造を用意します。
* 単位格子を2x2x3倍することでスーパーセルを作ります。$\left | \vec a \right |:\left | \vec b \right |:\left | \vec c \right |=2:2:3$となります。


```python
atoms = bulk("Si", cubic=True) * (2,2,3)
print(atoms.cell.cellpar())
```

    [10.86 10.86 16.29 90.   90.   90.  ]



* 次に、NPT integratorを`hydrostatic_strain=True`として初期化します。 この場合、stress tensorは`np.diag([np.trace(stress_tensor)]*3)`に置き換えられることになります。つまり、3つの対角項の値が同じになり、非対角項の値はゼロになります。結果として、$\vec a$, $\vec b$, $\vec c$の長さは等方的に変化することになります。


```python
integrator = NPTIntegrator(
    timestep=1.0,
    temperature=300,
    pressure=101325 * units.Pascal,
    ttime=10,
    pfactor=100 * units.fs,
    hydrostatic_strain=True
)

system = ASEMDSystem(atoms.copy())
system.init_temperature(300.0, stationary=True, rng=np.random.RandomState(seed=999))

cell_shape_log = []
md = MDFeature(integrator=integrator, n_run=100, show_progress_bar=True, show_logger=True, logger_interval=100, estimator_fn=estimator_fn)
md_results = md(system, extensions=[(PrintCellShape(cell_shape_log), 10)])
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /tmp/matlantis_1flndlmo/tmpl_k7vp2v.traj.
    Note: The max disk size of / is about 30G.
    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.51 eV/atom  temperature: 278.01 K  volume:  1921 Ang^3  density: 2.330 g/cm^3
    steps:     0  energy：-4.55 eV/atom  total energy: -4.51 eV/atom  temperature: 278.01 K  volume:  1921 Ang^3  density: 2.330 g/cm^3


    Dyn step    0 a 10.86 b 10.86 c 16.29 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step    0 a 10.86 b 10.86 c 16.29 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   10 a 10.87 b 10.87 c 16.30 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   20 a 10.89 b 10.89 c 16.33 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   30 a 10.92 b 10.92 c 16.37 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   40 a 10.95 b 10.95 c 16.42 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   50 a 10.98 b 10.98 c 16.47 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   60 a 11.00 b 11.00 c 16.50 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   70 a 11.00 b 11.00 c 16.51 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   80 a 11.00 b 11.00 c 16.50 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   90 a 10.98 b 10.98 c 16.47 alpha 90.00 beta 90.00 gamma 90.00


    steps:   100  energy：-4.51 eV/atom  total energy: -4.46 eV/atom  temperature: 383.36 K  volume:  1974 Ang^3  density: 2.268 g/cm^3


    Dyn step  100 a 10.96 b 10.96 c 16.44 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   50 a 10.98 b 10.98 c 16.47 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   60 a 11.00 b 11.00 c 16.50 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   70 a 11.00 b 11.00 c 16.51 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   80 a 11.00 b 11.00 c 16.50 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   90 a 10.98 b 10.98 c 16.47 alpha 90.00 beta 90.00 gamma 90.00


    steps:   100  energy：-4.51 eV/atom  total energy: -4.46 eV/atom  temperature: 383.36 K  volume:  1974 Ang^3  density: 2.268 g/cm^3


    Dyn step  100 a 10.96 b 10.96 c 16.44 alpha 90.00 beta 90.00 gamma 90.00



* シミュレーション中、格子定数の比$\left | \vec a \right |:\left | \vec b \right |:\left | \vec c \right |$ が2:2:3に固定されていることを確認できます。


```python
cell_shape = np.array(cell_shape_log)
print(cell_shape[:,:3] / cell_shape[:,:1])
```

    [[1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]
     [1.  1.  1.5]]



### 4.3 三斜晶系においてセル行列を上三角行列に変換する方法
* NPT integratorを三斜晶系(triclinic lattice)に適用する場合、セル行列が上三角行列になるように変換する必要があります。上三角行列でない場合、エラーが発生します。
* `MatlantisAtoms`にはこのような回転処理を行う`rotate_atoms_to_upper()`メソッドが実装されています。


* まず、三斜晶系SiO2の構造を用意します。MaterialsProjectからダウンロード可能です。(https://materialsproject.org/materials/mp-546794/#).
* 元のセルが下三角行列になっていることを確認します。


```python
ase_atoms = ase.io.read(str(dir_path/"assets/md_adv_usage/SiO2_mp-546794_computed.cif")) * (2,2,2)
print(ase_atoms.get_cell()[:])
print(ase_atoms.cell.cellpar())
atoms = MatlantisAtoms(ase_atoms.copy())
```

    [[10.27168474  0.          0.        ]
     [ 0.31570563 10.2668319   0.        ]
     [-5.29369519 -5.13341595  7.15068803]]
    [ 10.27168474  10.27168474  10.27168474 121.02203994 121.02203994
      88.23870671]



* 次に、`rotate_atoms_to_upper()`を用いて構造を回転させます。
* 基本格子ベクトルの長さや角度を保ったまま、セル行列が上三角行列に変換されていることが確認できます。
* `matlantis_feature.utils.atoms_utils`の`convert_atoms_to_upper`関数を使用して、ASEの`Atoms`に対して同じ回転操作を実行することもできます。


```python
atoms.rotate_atoms_to_upper()
print(atoms.ase_atoms.get_cell()[:])
print(atoms.ase_atoms.cell.cellpar())

# from matlantis_features.utils.atoms_util import convert_atoms_to_upper
# ase_atoms_upper = convert_atoms_to_upper(ase_atoms)
```

    [[ 8.34021852 -2.8151472  -5.29369519]
     [ 0.          8.80251661 -5.29369519]
     [ 0.          0.         10.27168474]]
    [ 10.27168474  10.27168474  10.27168474 121.02203994 121.02203994
      88.23870671]



* その後、NPT MDシミュレーションを通常通り実行することができます。
* 上のNotebookセルをコメントアウトして一連のコードを実行し直してみると、`ValueError`が発生します。


```python
integrator = NPTIntegrator(
    timestep=1.0,
    temperature=300,
    pressure=101325 * units.Pascal,
    ttime=10,
    pfactor=100 * units.fs,
)

system = ASEMDSystem(atoms)
system.init_temperature(300.0, stationary=True, rng=np.random.RandomState(seed=999))

md = MDFeature(integrator=integrator, n_run=100, show_progress_bar=True, show_logger=True, logger_interval=100, estimator_fn=estimator_fn)
md_results = md(system)
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /tmp/matlantis_hxjnpexv/tmptpfxun7g.traj.
    Note: The max disk size of / is about 30G.
    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.
    steps:     0  energy：-6.35 eV/atom  total energy: -6.31 eV/atom  temperature: 304.48 K  volume:   754 Ang^3  density: 2.117 g/cm^3
    steps:     0  energy：-6.35 eV/atom  total energy: -6.31 eV/atom  temperature: 304.48 K  volume:   754 Ang^3  density: 2.117 g/cm^3
    steps:   100  energy：-6.31 eV/atom  total energy: -6.28 eV/atom  temperature: 288.35 K  volume:   776 Ang^3  density: 2.058 g/cm^3


    Note: The max disk size of / is about 30G.


    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.


    steps:     0  energy：-6.35 eV/atom  total energy: -6.31 eV/atom  temperature: 304.48 K  volume:   754 Ang^3  density: 2.117 g/cm^3


    steps:     0  energy：-6.35 eV/atom  total energy: -6.31 eV/atom  temperature: 304.48 K  volume:   754 Ang^3  density: 2.117 g/cm^3


    steps:   100  energy：-6.31 eV/atom  total energy: -6.28 eV/atom  temperature: 288.35 K  volume:   776 Ang^3  density: 2.058 g/cm^3



&lt;a id='advanced_5'&gt;&lt;/a&gt;
## 5. check point fileからの計算再開方法

* まれに、trajectoryが生成される前にMDが異常終了することがあります。Matlantis-featuresではMDをcheckpoint fileから再開するメソッドを提供しています。以下ではこの機能の使い方を紹介します。


```python
import logging
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from matlantis_features.features.md import MDFeature, ASEMDSystem, NVTBerendsenIntegrator, MDExtensionBase
from matlantis_features.features.md.ase_md_system import ASEMDState
```


* 毎ステップごとにポテンシャルエネルギーを保存するextensionを作っておきます。


```python
class EnergyLog(MDExtensionBase):
    def __init__(self, energy_log):
        self.energy_log = energy_log
    def __call__(self, system, integrator):
        self.energy_log.append((system.current_total_step, system.ase_atoms.get_potential_energy()))
```


### 通常実行
* まず、MDシミュレーションを普通に実行し、各ステップのポテンシャルエネルギーの値をリスト`energy_log_1`に格納します。のちほど、MDが途中で再開されても同じtrajectoryが生成されるか確認するために用います。


```python
logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

atoms = bulk("Si", cubic=True) * (2,2,2)

nvtbi = NVTBerendsenIntegrator(timestep=1.0, temperature=2000, taut=10.0)
mdsys = ASEMDSystem(atoms)
mdsys.init_temperature(2000.0, stationary=True, rng=np.random.RandomState(seed=99))

energy_log_1 = []
extension = EnergyLog(energy_log_1)

md = MDFeature(
    integrator=nvtbi, n_run=100, show_progress_bar=True, show_logger=True, logger_interval=10, estimator_fn=estimator_fn
)
md_results = md(mdsys, extensions=[(extension, 1)])
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /tmp/matlantis_tmc06z_6/tmpmhvb_c22.traj.
    Note: The max disk size of / is about 30G.
    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.29 eV/atom  temperature: 1980.63 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    10  energy：-4.45 eV/atom  total energy: -4.27 eV/atom  temperature: 1439.68 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    20  energy：-4.33 eV/atom  total energy: -4.17 eV/atom  temperature: 1233.29 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    30  energy：-4.32 eV/atom  total energy: -4.10 eV/atom  temperature: 1712.70 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    40  energy：-4.32 eV/atom  total energy: -4.08 eV/atom  temperature: 1861.41 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    50  energy：-4.28 eV/atom  total energy: -4.05 eV/atom  temperature: 1727.32 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    60  energy：-4.20 eV/atom  total energy: -4.00 eV/atom  temperature: 1517.45 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    70  energy：-4.14 eV/atom  total energy: -3.94 eV/atom  temperature: 1522.62 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    80  energy：-4.12 eV/atom  total energy: -3.89 eV/atom  temperature: 1770.07 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    90  energy：-4.10 eV/atom  total energy: -3.87 eV/atom  temperature: 1781.70 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:   100  energy：-4.03 eV/atom  total energy: -3.83 eV/atom  temperature: 1589.89 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    Note: The max disk size of / is about 30G.


    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.


    steps:     0  energy：-4.55 eV/atom  total energy: -4.29 eV/atom  temperature: 1980.63 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    10  energy：-4.45 eV/atom  total energy: -4.27 eV/atom  temperature: 1439.68 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    20  energy：-4.33 eV/atom  total energy: -4.17 eV/atom  temperature: 1233.29 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    30  energy：-4.32 eV/atom  total energy: -4.10 eV/atom  temperature: 1712.70 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    40  energy：-4.32 eV/atom  total energy: -4.08 eV/atom  temperature: 1861.41 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    50  energy：-4.28 eV/atom  total energy: -4.05 eV/atom  temperature: 1727.32 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    60  energy：-4.20 eV/atom  total energy: -4.00 eV/atom  temperature: 1517.45 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    70  energy：-4.14 eV/atom  total energy: -3.94 eV/atom  temperature: 1522.62 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    80  energy：-4.12 eV/atom  total energy: -3.89 eV/atom  temperature: 1770.07 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    90  energy：-4.10 eV/atom  total energy: -3.87 eV/atom  temperature: 1781.69 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:   100  energy：-4.03 eV/atom  total energy: -3.83 eV/atom  temperature: 1589.89 K  volume:  1281 Ang^3  density: 2.330 g/cm^3



* 各ステップでのポテンシャルエネルギーをプロットしてみます。


```python
energy1 = np.array(energy_log_1)

plt.scatter(energy1[:,0], energy1[:,1], s=1)
plt.xlabel("MD steps")
plt.ylabel("Potential energy (eV)")
```




    Text(0, 0.5, 'Potential energy (eV)')





![png](output_76_1.png)




### 停止と再開
* 43回目のMDステップで例外を発生させるextensionを定義します。


```python
class RaiseException(MDExtensionBase):
    def __call__(self, system, integrator):
        if system.current_total_step == 43:
            raise ValueError
```


* 上の例と同じMDを実行します。MDを再開したい場合、 `MDFeature`のパラメータ`checkpoint_file_name`と`checkpoint_freq`を設定します。 MD featureは、MD systemとintegratorの内部状態を`checkpoint_freq`ステップごとにcheckpointファイルに保存します。
* また、パラメータ`traj_file_name`と`traj_freq`を設定しておき、MDのtrajectoryも保存しておきます。


```python
logger = logging.getLogger("matlantis_features")
logger.setLevel(logging.INFO)

atoms = bulk("Si", cubic=True) * (2,2,2)

nvtbi = NVTBerendsenIntegrator(
    timestep=1.0,
    temperature=2000,
    taut=10.0
    )
mdsys = ASEMDSystem(atoms)
mdsys.init_temperature(2000.0, stationary=True, rng=np.random.RandomState(seed=99))

energy_log_2 = []
extension = EnergyLog(energy_log_2)
exception = RaiseException()

md = MDFeature(
    integrator=nvtbi,
    n_run=100,
    checkpoint_file_name="Si_2000K.ckpt",
    checkpoint_freq=10,
    traj_file_name="Si_2000K.traj",
    traj_freq=10,
    show_progress_bar=True,
    show_logger=True,
    logger_interval=10,
    estimator_fn=estimator_fn,
)
try:
    md_results = md(mdsys, extensions=[(extension, 1), (exception, 1)])
except ValueError as e:
    print("Exception raised ", e)
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/Si_2000K.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.29 eV/atom  temperature: 1980.63 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    10  energy：-4.45 eV/atom  total energy: -4.27 eV/atom  temperature: 1439.68 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    20  energy：-4.33 eV/atom  total energy: -4.17 eV/atom  temperature: 1233.29 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    30  energy：-4.32 eV/atom  total energy: -4.10 eV/atom  temperature: 1712.70 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    40  energy：-4.32 eV/atom  total energy: -4.08 eV/atom  temperature: 1861.41 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    Exception raised


    steps:    40  energy：-4.32 eV/atom  total energy: -4.08 eV/atom  temperature: 1861.41 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    Exception raised



* 予め設定しておいたとおり、43ステップ目にMDが停止しました。checkpointファイル`Si_2000K.ckpt`とtrajectoryファイル`Si_2000K.traj`が左のファイル一覧から確認できます。
* `ASEMDState`を用いてcheckpointファイルからMD systemとintegratorの状態を復元します。
* checkpointファイルは10ステップごとに保存されるため、最後に保存されたのは40ステップ目であることがわかります。


```python
atoms = bulk("Si", cubic=True) * (2,2,2)

nvtbi = NVTBerendsenIntegrator(
    timestep=1.0,
    temperature=2000,
    taut=10.0
    )
mdsys = ASEMDSystem(atoms)

state = ASEMDState.from_file("Si_2000K.ckpt")
state.restore_to(mdsys, nvtbi)

state.system_state["total_step"]
```




    40




* MD featureを先の例と同様に初期化します。注意すべき点として、過去のtrajectoryファイルが上書きされないように`traj_append`を`True`に設定してください。
* さらに、`n_run`を100から60に変更します。MDを40ステップ目から再開するためです。


```python
md = MDFeature(
    integrator=nvtbi,
    n_run=60,
    checkpoint_file_name="Si_2000K.ckpt",
    checkpoint_freq=10,
    traj_file_name="Si_2000K.traj",
    traj_freq=10,
    traj_append=True,
    show_progress_bar=True,
    show_logger=True,
    logger_interval=10,
    estimator_fn=estimator_fn,
)
md_results = md(mdsys, extensions=[(extension, 1)])
```


      0%|          | 0/60 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /home/jovyan/matlantis-features/benchmarks/Si_2000K.traj.
    Note: The max disk size of /home/jovyan is about 99G.
    steps:     0  energy：-4.32 eV/atom  total energy: -4.08 eV/atom  temperature: 1861.41 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    10  energy：-4.28 eV/atom  total energy: -4.05 eV/atom  temperature: 1727.32 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    20  energy：-4.20 eV/atom  total energy: -4.00 eV/atom  temperature: 1517.45 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    30  energy：-4.14 eV/atom  total energy: -3.94 eV/atom  temperature: 1522.62 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    40  energy：-4.12 eV/atom  total energy: -3.89 eV/atom  temperature: 1770.07 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    50  energy：-4.10 eV/atom  total energy: -3.87 eV/atom  temperature: 1781.70 K  volume:  1281 Ang^3  density: 2.330 g/cm^3
    steps:    60  energy：-4.03 eV/atom  total energy: -3.83 eV/atom  temperature: 1589.89 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    Note: The max disk size of /home/jovyan is about 99G.


    steps:     0  energy：-4.32 eV/atom  total energy: -4.08 eV/atom  temperature: 1861.41 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    10  energy：-4.28 eV/atom  total energy: -4.05 eV/atom  temperature: 1727.32 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    20  energy：-4.20 eV/atom  total energy: -4.00 eV/atom  temperature: 1517.45 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    30  energy：-4.14 eV/atom  total energy: -3.94 eV/atom  temperature: 1522.62 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    40  energy：-4.12 eV/atom  total energy: -3.89 eV/atom  temperature: 1770.07 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    50  energy：-4.10 eV/atom  total energy: -3.87 eV/atom  temperature: 1781.69 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    steps:    60  energy：-4.03 eV/atom  total energy: -3.83 eV/atom  temperature: 1589.89 K  volume:  1281 Ang^3  density: 2.330 g/cm^3



* MDシミュレーションが40ステップ目に開始し100ステップ目に終了します。
* 計算再開前と後のポテンシャルエネルギーはどちらも `energy_log_2`に格納されています。これをプロットし、通常通り実行したMDのポテンシャルエネルギー`energy_log_1`と比較してみます。
* 両者のポテンシャルエネルギーが一致していることがわかります。このように、MDが途中で意図せず停止しても、上記の手順で計算を再開させることで同じ結果を得ることが可能です。


```python
energy2 = np.array(energy_log_2)

plt.scatter(energy1[:,0], energy1[:,1], s=1)
plt.scatter(energy2[:,0], energy2[:,1], s=1)

plt.xlabel("MD steps")
plt.ylabel("Potential energy (eV)")
```




    Text(0, 0.5, 'Potential energy (eV)')





![png](output_86_1.png)




&lt;a id='advanced_6'&gt;&lt;/a&gt;
## 6. MD integratorを構築する方法

matlantis-featuresの`integrators`にはすでにいくつかのASE Intergratorのラッパーが実装されていますが、 次の例では、必要なインテグレータを自由に定義できるよう、ASEのMDメソッドからmatlantis-featuresで`intergretor`クラスを作成する方法を示します。

ここでは、例としてASEの`Inhomogeneous_NPTBerendsen`を使用します。類似の手法としてNPTBerendsenが存在しますが、NPTBerendsenでは格子の形状（格子長の比率 a：b：c と3つの軸の角度 $\alpha$, $\beta$, $\gamma$）は変更されないのに対し、Inhomogeneous NPTBerendsenでは格子長の比率 a:b:cが変更されます。


```python
from typing import List

import numpy as np
from ase import Atoms, units
from ase.build import bulk
from ase.optimize.optimize import Dynamics
from ase.md.nptberendsen import Inhomogeneous_NPTBerendsen

from matlantis_features.features.md.ase_integrators import ASEIntegrator
from matlantis_features.features.md import NPTIntegrator, ASEMDSystem, MDFeature, MDExtensionBase, NPTBerendsenIntegrator
```



`InhomogeneousNPTBerendsenIntegrator`を作成します。
* まず、新しいintegratorを作成する際には`ASEIntegrator`を継承する必要があります。
* 重要な点として、`InhomogeneousNPTBerendsenIntegrator`の初期化関数`__init__`で`_get_dyn`関数を必ず定義してください。 `_get_dyn`関数は、入力パラメーターとしてASEの`Atoms`オブジェクトをとり、ASEの`Inhomogeneous_NPTBerendsen`オブジェクトを出力します。


```python
class InhomogeneousNPTBerendsenIntegrator(ASEIntegrator):
    def __init__(
        self,
        timestep: float,
        temperature: float,
        pressure: float,
        taut: float = 500,
        taup: float = 1000,
        compressibility: float = 67.2,
        fixcm: bool = True,
        mask: List[bool] = [True, True, True]
    ):
        def _get_dyn(atoms: Atoms) -&gt; Dynamics:
            return Inhomogeneous_NPTBerendsen(
                atoms,
                timestep=self.timestep * units.fs,
                temperature_K=self.temperature,
                pressure_au=self.pressure,
                taut=self.taut * units.fs,
                taup=self.taup * units.fs,
                compressibility_au=self.compressibility,
                fixcm=self.fixcm,
                mask=self.mask
            )

        super().__init__(get_dynamics_fn=_get_dyn, timestep=timestep)
        self.temperature = temperature
        self.pressure = pressure
        self.taut = taut
        self.taup = taup
        self.compressibility = compressibility
        self.fixcm = fixcm
        self.mask = mask
```



`InhomogeneousNPTBerendsenIntegrator`をテストしてみます。 まず、MDの実行中に格子形状の情報を出力するための`extension`を定義します。


```python
class PrintCellShape(MDExtensionBase):
    def __init__(self, cell_log=None):
        self.cell_log = cell_log
    def __call__(self, system, integrator) -&gt; None:
        cell_par = system.ase_atoms.cell.cellpar()
        istep = system.current_total_step
        print(f"Dyn step {istep:4d} a {cell_par[0]:3.5f} b {cell_par[1]:3.5f} c {cell_par[2]:3.5f} alpha {cell_par[3]:3.2f} beta {cell_par[4]:3.2f} gamma {cell_par[5]:3.2f} ")
        if self.cell_log is not None:
            self.cell_log.append(cell_par)
```



MDの実行を開始します。 InhomogeneousNPTBerendsenIntegratorを使用すると、格子の3つの軸の長さの変化が不均一であることがわかります。

興味があれば、NPTBerendsenIntegratorを使用して同じタスクを再実行し結果を比較してみましょう。 軸 $a$ $b$ $c$ の変化は一様であることが確認できます。


```python
integrator = InhomogeneousNPTBerendsenIntegrator(
    timestep=1.0,
    temperature=300,
    pressure=101325 * units.Pascal,
)
# integrator = NPTBerendsenIntegrator(
#     timestep=1.0,
#     temperature=300,
#     pressure=101325 * units.Pascal,
# )

system = ASEMDSystem(bulk("Si", cubic=True)*(2,2,2))
system.init_temperature(300.0, stationary=True, rng=np.random.RandomState(seed=999))

md = MDFeature(integrator=integrator, n_run=100, show_progress_bar=True, show_logger=True, logger_interval=100, estimator_fn=estimator_fn)
md_results = md(system, extensions=[(PrintCellShape(), 10)])
```


      0%|          | 0/100 [00:00&lt;?, ?it/s]


    The MD trajectory will be saved at /tmp/matlantis_g_oz_jt3/tmpotl7q_8h.traj.
    Note: The max disk size of / is about 30G.
    Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.
    steps:     0  energy：-4.55 eV/atom  total energy: -4.51 eV/atom  temperature: 278.85 K  volume:  1281 Ang^3  density: 2.330 g/cm^3


    Dyn step    0 a 10.86000 b 10.86000 c 10.86000 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   10 a 10.88653 b 10.88611 c 10.88563 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   20 a 10.90434 b 10.90424 c 10.90396 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   30 a 10.91466 b 10.91486 c 10.91531 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   40 a 10.91972 b 10.91982 c 10.92068 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   50 a 10.92332 b 10.92374 c 10.92380 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   60 a 10.92617 b 10.92778 c 10.92599 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   70 a 10.92716 b 10.92989 c 10.92668 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   80 a 10.92800 b 10.93063 c 10.92761 alpha 90.00 beta 90.00 gamma 90.00
    Dyn step   90 a 10.93055 b 10.93241 c 10.93037 alpha 90.00 beta 90.00 gamma 90.00


    steps:   100  energy：-4.54 eV/atom  total energy: -4.51 eV/atom  temperature: 175.11 K  volume:  1307 Ang^3  density: 2.284 g/cm^3


    Dyn step  100 a 10.93276 b 10.93451 c 10.93271 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   40 a 10.91972 b 10.91982 c 10.92068 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   50 a 10.92332 b 10.92374 c 10.92380 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   60 a 10.92617 b 10.92778 c 10.92599 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   70 a 10.92716 b 10.92989 c 10.92668 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   80 a 10.92800 b 10.93063 c 10.92761 alpha 90.00 beta 90.00 gamma 90.00


    Dyn step   90 a 10.93055 b 10.93241 c 10.93037 alpha 90.00 beta 90.00 gamma 90.00


    steps:   100  energy：-4.54 eV/atom  total energy: -4.51 eV/atom  temperature: 175.11 K  volume:  1307 Ang^3  density: 2.284 g/cm^3


    Dyn step  100 a 10.93276 b 10.93451 c 10.93271 alpha 90.00 beta 90.00 gamma 90.00

