# How to Use LightPFP

LightPFPは、従来のPFPモデルと比較して、計算速度が速く、扱えるスケールが大きいモデルです。このドキュメントでは、LightPFPモデルの訓練方法と、原子シミュレーションでの使用方法について、基本的なワークフローを説明します。

このプロセスは、3つの主要な部分に分けられます。

**パート1: 訓練データセットの生成**

このパートでは、LightPFPモデルの訓練に直接使用できるHDF5形式の訓練データセットを生成する方法を学びます。以下の2つの方法を紹介します。

1. **提供されるCLIツールの使用**: 訓練データセットを生成するためのコマンドラインインターフェース(CLI)ツールを提供しています。

2. **ユーザーデータ(MDトラジェクトリなど)を訓練データセットに変換**: 既存のデータ(分子動力学(MD)トラジェクトリなど)がある場合、それらを必要な訓練データセット形式に変換する方法を示します。

データセットの生成後、データセットの検証と検査方法についても説明し、データセットの品質を確認します。

**パート2: LightPFPモデルの訓練**

このパートでは、LightPFPモデルの訓練タスクを送信する方法を学びます。訓練データセットのアップロード、訓練パラメータの設定、ジョブの送信プロセスについて説明します。さらに、訓練の進行状況を監視する方法と、訓練が完了したらモデルライブラリからモデルにアクセスする方法についても示します。

**パート3: LightPFPモデルを使った原子シミュレーション**

このパートでは、訓練したLightPFPモデルを使って分子動力学(MD)シミュレーションを実行する方法を学びます。このガイドに従えば、シミュレーションパラメータの設定、システムの初期化、パート2で訓練したLightPFPモデルを使ったシミュレーションの実行方法がわかります。

**パート4：LightPFPモデルの評価**

最後のパートでは、ツールキットを使用してLightPFPの性能を評価する方法を示します。これは、実際の物理的特性計算においてLightPFPをPFPと比較することで行います。

このチュートリアルでは、ベンチマーク材料としてNi3Alを使用します。

このガイドに従えば、Ni3Alに対して、効率的で大規模な原子シミュレーションを行うためのLightPFPモデルの訓練と活用ができるようになります



LightPFPのトレーニングと推論に必要なパッケージをインストールしましょう。

* light-pfp-data：トレーニングデータセットコレクションのデータパッケージです。
* light-pfp-client：LightPFP推論のためのクライアントパッケージです。
* light-pfp-evaluate: LightPFPの性能を評価するためのパッケージ。


```
! pip install -U light-pfp-data light-pfp-client light-pfp-evaluate
```


## 1. 訓練データセットの生成

### 1.1 コマンドラインインターフェース（CLI）を使用する方法

以下のコマンドを使用して、`light-pfp-data` パッケージをインストールした後に `dataset-generation` CLI を使用することができます。

```
pip install -U light-pfp-data
```

`dataset-generation CLI` は、1つまたは複数の初期構造から訓練構造を生成するためのいくつかの方法を提供しています。これらの方法には、以下が含まれます：
* MDシミュレーションを使用した構造のサンプリング
* 入力構造の均等な圧縮と伸長
* 初期構造に歪みを加える
* 原子の変位を適用する
* 空孔を持つ構造の作成
* 表面構造の作成
* 元素の置換による構造の作成

これらの方法の詳細な説明は、ガイドブックに記載されています。

`dataset-generation` CLI を実行するためには、初期構造とJSON形式の制御ファイルを準備する必要があります。この例では、初期構造としてNi3Al結晶(mp-2593)の構造を使用します（[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) の下で &lt;a href="https://next-gen.materialsproject.org/materials/mp-2593"&gt;The Materials Project&lt;/a&gt;からダウンロード）。./assets/how_to_use_lightpfp/dataset.jsonには、制御ファイルの例が提供されています。制御ファイルの引数に関する情報については、ガイドブックを参照してください。

以下のコマンドでデータセット生成プロセスを開始できます。PFPを使用して訓練構造を生成するため、このステップは比較的時間がかかります（約1時間）。終了すると、AlNi3.h5という名前のHDF5データセットが生成され、次のステップでLightPFPモデルの訓練に使用できます。

（**処理結果のデータセットをあらかじめ assets/how_to_use_lightpfp/AlNi3.h5 に保存しているので、このコマンドを実行する必要はありません**）

```
PFP_DEFAULT_PRIORITY=1 dataset-generation -c assets/how_to_use_lightpfp/dataset.json --log-level info --progress-bar "True" -j 8
```

データ生成によるPFPトークン消費を抑えるために、上記のように PFP_DEFAULT_PRIORITY を設定することを推奨します。



## 1.2 ユーザーのデータの変換

ユーザーは、分子動力学（MD）トラジェクトリなどの自分のデータを、LightPFPモデルの訓練に使用するためのHDF5データセットに変換することもできます。

**ASE MDトラジェクトリファイルのCLIを使用する方法**
ASE MDトラジェクトリファイルの場合、変換タスクを実行するためのCLIツールを提供しています：
```
traj-to-dataset --traj assets/how_to_use_lightpfp/ase_traj.traj --dataset user_custom.h5 --calc_mode PBE --model_version v8.0.0
```

注意：PFPモデルバージョンと計算モードの情報は、既知の場合、提供することを推奨します。
これらの情報はデータセットの検証と適切な事前学習モデルの選択に役立ちます。
ただし、この情報が不明な場合や、データセットがDFT計算などの他のソースから生成された場合は、これらのフラグを省略しても問題ありません。

**DFT出力用のCLIを使用する**

ユーザーは`traj-to-dataset` CLIを使用して、DFTなどの外部ソースから分子動力学（MD）のトラジェクトリを学習データセット形式に変換することができます。
入力ファイルのフォーマットは`--format`引数で指定します。
例えば、Quantum ESPRESSOの出力を変換する場合は以下のようになります：

```
traj-to-dataset --traj assets/how_to_use_lightpfp/si.md8.out --dataset user_qe.h5 --format espresso-out
```

`traj-to-dataset` CLIツール以外にも、ユーザーは以下の方法でPythonの関数を使用して訓練データセットを作成することもできます。

**単一のASEアトムを保存する**
* PFP Calculatorを用いてASE atomsのエネルギー・力・応力を計算します。 次に、結果をデータセットに保存します。



```
from ase.io import read

import light_pfp_data
from light_pfp_data.sample.custom import traj_to_dataset, atoms_to_dataset
from light_pfp_data.utils.dataset import H5DatasetWriter

from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator as PFPASECalculator
from pfp_api_client.pfp.estimator import Estimator as PFPEstimator


print(f"light_pfp_data: {light_pfp_data.__version__}")

atoms = read("assets/how_to_use_lightpfp/AlNi3_mp-2593_conventional_standard.cif")
calc = PFPASECalculator(PFPEstimator(model_version="v8.0.0", calc_mode="PBE"))

atoms.calc = calc
E = atoms.get_potential_energy()
F = atoms.get_forces()
v = atoms.get_stress()
print(f"Energy {E} eV")

with H5DatasetWriter("AlNi3_manual.h5", mode="append") as dataset:
    atoms_to_dataset(
        atoms,
        dataset=dataset,
        model_version="v8.0.0",  # can be omitted if the model version is unknown
        calc_mode="PBE",  # can be omitted if the calculation mode is unknown
    )
```

    light_pfp_data: 1.3.1
    Energy -19.687408888675506 eV




**MD trajectoryを保存します**

* まずはASEを用いてtrajectoryを作成します。


```
from ase import units
from ase.io import Trajectory
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

# Set the momenta corresponding to the given temperature
MaxwellBoltzmannDistribution(atoms, temperature_K=300.0)
Stationary(atoms)  # Set zero total momentum to avoid drifting

traj = Trajectory("md.traj", "w", atoms=atoms)
md = Langevin(atoms, timestep=1.0*units.fs, temperature_K=300.0, friction=0.1, logfile="-", loginterval=100)
md.attach(traj, interval=100)
md.run(steps=1000)
traj.close()

with H5DatasetWriter("AlNi3_md.h5", mode="append") as dataset:
    traj_to_dataset(
        "md.traj", 
        dataset=dataset,
        model_version="v8.0.0",  # can be omitted if the model version is unknown
        calc_mode="PBE",  # can be omitted if the calculation mode is unknown
    )
```

    Time[ps]      Etot[eV]     Epot[eV]     Ekin[eV]    T[K]
    0.0000         -19.5643     -19.6874       0.1232   238.2


    0.1000         -19.2630     -19.5263       0.2633   509.2


    0.2000         -19.3071     -19.6002       0.2931   566.9


    0.3000         -19.4223     -19.5580       0.1357   262.4


    0.4000         -19.4952     -19.5536       0.0585   113.1


    0.5000         -19.5063     -19.5225       0.0162    31.4


    0.6000         -19.5529     -19.6276       0.0747   144.5


    0.7000         -19.5429     -19.6060       0.0630   121.9


    0.8000         -19.5178     -19.6190       0.1012   195.8


    0.9000         -19.5116     -19.5803       0.0687   132.9


    1.0000         -19.5204     -19.5767       0.0562   108.8



### 1.3 データセットの検証

ユーザーは、提供されたコマンドを使用するか、カスタムのPythonコードを記述することで、訓練データセットを検証および確認することができます。LightPFPモデルの訓練を進める前に、データセットの品質と一貫性を確認するために、このステップは重要です。

```
validate-dataset --datasets assets/how_to_use_lightpfp/AlNi3.h5 AlNi3_manual.h5 AlNi3_md.h5 --max-energy -3.0 --max-forces 30.0
```
このステップでは、以下のチェックが行われます：

* PFPモデルバージョンと計算モードの一貫性：これらの値がデータセットに保存されている場合、PFPモデルバージョンと計算モードがすべてのデータポイントで一貫していることを確認します。モデルバージョンと計算モードがデータセットに保存されていない場合、この一貫性チェックはスキップされます。
* 原子当たりの最大エネルギー：最大エネルギー（原子当たり）が指定された --max-energy 値（この例では -3.0 eV/atom）を超えていないことを検証します。
* 力の最大値：力の最大値が指定された --max-forces 値（この例では30.0 eV/Å）を超えていないことをチェックします。

これらの条件のいずれかを満たさない場合、警告メッセージが表示されます。問題のあるデータポイントを削除したい場合は、上記のコマンドを `--delete-invalid-keys` フラグを指定して再実行することができます。詳細については、ガイドブックを参照してください。

ユーザーはまた、以下の例に示されるようにPythonコードを使用してデータセットを調査することもできます：


```
import h5py

f = h5py.File('AlNi3_md.h5', 'r')

# Print number of groups in the dataset.
print(f"Number of groups: {len(f.keys())}")

# Inspect first group in the dataset.
dset = f[list(f.keys())[0]]
for k in dset.keys():
    print(f"{k}: {dset[k][()]}")
```

    Number of groups: 11
    cell: [[3.561084 0.       0.      ]
     [0.       3.561084 0.      ]
     [0.       0.       3.561084]]
    forces: [[-4.5945012e-07 -8.4279242e-07  6.3958328e-07]
     [ 1.1831097e-06  3.0871541e-07 -3.6671031e-07]
     [ 1.9860716e-08  3.3368090e-07 -4.7203193e-07]
     [-7.4352033e-07  2.0039606e-07  1.9915898e-07]]
    numbers: [13 28 28 28]
    pbc: [ True  True  True]
    positions: [[0.       0.       0.      ]
     [1.780542 1.780542 0.      ]
     [1.780542 0.       1.780542]
     [0.       1.780542 1.780542]]
    potential_energy: -19.687408888675506
    stress: [-5.252235e-03 -5.252240e-03 -5.252191e-03 -8.149200e-09  2.812724e-10
      8.557565e-09]



## 2. LightPFPモデルの訓練

### 2.1 LightPFP訓練タスクの送信

訓練データセットを取得したので、訓練データセットとJSON形式の制御ファイルを使用して、LightPFPモデルの訓練ジョブを開始することができます。この例では、事前に訓練されたモデルである "ALL_ELEMENTS_LARGE_6" から始めて、LightPFPモデルを500エポックだけ訓練します。JSON制御ファイルの詳細な説明については、ガイドブックを参照してください。

**訓練ジョブの制御ファイル**
```
{
    "common_config": {
        "total_epoch": 500
    },
    "mtp_config": {
        "pretrained_model": "ALL_ELEMENTS_SMALL_6"
    }
}
```

LightPFP訓練の送信は、以下のコマンドでトリガーできます：

```
submit-light-pfp-training -d assets/how_to_use_lightpfp/AlNi3.h5 AlNi3_manual.h5 AlNi3_md.h5 -c assets/how_to_use_lightpfp/job.json --model-name AlNi3_example
```

訓練プロセス中は、進捗状況を監視し、訓練ログを確認することができます。訓練が完了すると、訓練済みのLightPFPモデルがモデルライブラリで利用可能になり、生成されたIDを使用してアクセスすることができます。

訓練プロセスには、データセットのサイズや利用可能な計算リソースによって、かなりの時間がかかる場合があることに注意してください。この例では、データセットが比較的小さいため、エポック数を大きくしています。もしあなたのデータセットが大きい場合は、トレーニング時間とコストを短縮するために、エポック数を減らすことを検討してください。


```
!submit-light-pfp-training -d assets/how_to_use_lightpfp/AlNi3.h5 AlNi3_manual.h5 AlNi3_md.h5 --train-config-path assets/how_to_use_lightpfp/job.json --model-name AlNi3_example
```

    2025/09/19 00:49:59 Compressing file assets/how_to_use_lightpfp/AlNi3.h5...


    2025/09/19 00:49:59 Compressing file AlNi3_manual.h5...
    2025/09/19 00:49:59 Compressing file AlNi3_md.h5...


    2025/09/19 00:50:04 Number of atoms per dataset: [105450 4 44]
    2025/09/19 00:50:04 Completed creating the input zip file for job AlNi3_example


    2025/09/19 00:50:05 Complete! Training job ID: tiww3wo4v4njbm6n will be started.




### 2.2 訓練の進捗状況を確認する

訓練ジョブを送信した後、ユーザーはModelライブラリの `Launcher` -&gt; `LightPFP Models` で進捗状況を確認することができます。

&lt;img src="assets/how_to_use_lightpfp/model_library_1.png"&gt;

訓練中およびキュー中のすべての訓練タスクは、「In Progress」パネルで確認することができます。モデルアイテムをダブルクリックすると、完了したエポック数や現在のモデルのMAEを確認することができます。

&lt;img src="assets/how_to_use_lightpfp/model_library_2.png"&gt;

タスクが完了した後は、「Trained」パネルで結果を表示することができます。詳細を確認するには、モデルをダブルクリックします。

&lt;img src="assets/how_to_use_lightpfp/model_library_3.png"&gt;



# 学習済みLightPFPモデルを用いてMDシミュレーションを行う


### LightPFPのクライアントモジュールをインポートする


```
import light_pfp_client
from light_pfp_client import Estimator, ASECalculator

print(f"light_pfp_client: {light_pfp_client.__version__}")
```

    light_pfp_client: 1.1.4



### LightPFPのEstimatorを用いて推論を行う


```
from ase.io import read

# Replace with your trained model id.
model_id = "examplesv1alni3s"

# Get initial structure
atoms = read("assets/how_to_use_lightpfp/AlNi3_mp-2593_conventional_standard.cif")
print(f"Number of atoms: {len(atoms)}")

# Load MTP calculator
estimator = Estimator(model_id=model_id)
calc = ASECalculator(estimator)
calc.reset()
atoms.calc = calc

print(f"{atoms.get_potential_energy()=}")
print(f"{atoms.get_forces()=}")
```

    Number of atoms: 4


    atoms.get_potential_energy()=-19.67413330078125
    atoms.get_forces()=array([[ 4.47034836e-08, -7.45058060e-09,  1.49011612e-08],
           [-1.30385160e-08,  1.49011612e-08,  3.72529030e-08],
           [-5.58793545e-09, -1.86264515e-08,  3.72529030e-09],
           [-3.16649675e-08,  2.23517418e-08, -4.47034836e-08]])



### 結果をPFPのEstimatorと比較する


```
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator as PFPASECalculator
from pfp_api_client.pfp.estimator import Estimator as PFPEstimator

print(f"pfp_api_client: {pfp_api_client.__version__}")

pfp_calculator = PFPASECalculator(PFPEstimator(model_version="v8.0.0"))
atoms.calc = pfp_calculator
pfp_calculator.reset()
print(f"{atoms.get_potential_energy()=}")
print(f"{atoms.get_forces()=}")
```

    pfp_api_client: 1.23.1
    atoms.get_potential_energy()=-19.68740839419648
    atoms.get_forces()=array([[-2.77872020e-07,  5.22507128e-08, -3.01703441e-07],
           [ 5.62786813e-08, -2.66802763e-07,  1.16276027e-06],
           [ 7.60982725e-08, -3.65076872e-07, -5.78614338e-07],
           [ 1.45495066e-07,  5.79628922e-07, -2.82442495e-07]])



### LightPFPのEstimatorを用いて、MDシミュレーションを実行する


```
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

# Enable bookkeeping for estimator.
estimator = Estimator(model_id=model_id, book_keeping=True)
calc = ASECalculator(estimator)
calc.reset()

# Set calculator back to MTP calculator.
atoms = read("assets/how_to_use_lightpfp/AlNi3_mp-2593_conventional_standard.cif") * (4,4,4)
atoms.calc = calc

# MD settings
temperature = 300.0
steps = 10000

# Set the momenta corresponding to the given temperature
MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
Stationary(atoms)  # Set zero total momentum to avoid drifting

# Run MD
md = Langevin(
    atoms,
    timestep=units.fs,
    temperature_K=temperature,
    friction=0.2,
    logfile="-",
    loginterval = 100,
)
md.run(steps=steps)
atoms.write("final.xyz")
```



MDシミュレーションが終了したら、book_keepingキャッシュをクリアすることをお勧めします。そうする主な理由は2つあります：

* ブックキーピング・キャッシュをクリアすることで、隣接リストの情報が確実に削除され、GPUメモリが他のタスクや推論用に解放されます。
* 後続のタスクでブックキーピング・メソッドを意図せず使用し、パフォーマンスの問題につながることを防ぎます。


```
estimator.clean_up()
```



## 4. LightPFPモデルの評価



`light-pfp-evaluate`パッケージは、実世界の原子シミュレーションにおけるLightPFPモデルの信頼性を評価します。モデルの精度は訓練中にテストデータセットで測定できますが、実際のシミュレーションでの性能はより不確実です。このパッケージを使用することで、ユーザーは簡単な物理的特性の計算を行うことができ、応用シナリオでのモデルの信頼性を実践的にテストすることができます。


```
import os

from IPython import display

from ase.io import read
from matlantis_features.utils.calculators.pfp_api_calculator import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode
from light_pfp_client.estimator_fn import light_pfp_estimator_fn
from light_pfp_evaluate import (
    evaluate_elastic,
    evaluate_eos,
    evaluate_md,
    evaluate_phonon,
    evaluate_lattice,
    evaluate_vacancy,
    evaluate_surface
)

# Define estimator_fn for PFP and LightPFP
estimator_fn_pfp = pfp_estimator_fn(model_version="v8.0.0", calc_mode=EstimatorCalcMode.PBE)
estimator_fn_light_pfp = light_pfp_estimator_fn(model_id=model_id)
# Read input structure for validation
atoms = read("assets/how_to_use_lightpfp/AlNi3_mp-2593_conventional_standard.cif")
# Make working directory to save the evaluation results
os.makedirs("eval_results", exist_ok=True)
```



### Lattice

このテストでは、LightPFPとPFPの両方を使用してNi3Al超格子の形状を最適化します。その後、各モデルから得られたセルの長さと角度を比較します。


```
length_error, angle_error = evaluate_lattice(atoms, estimator_fn_pfp, estimator_fn_light_pfp, "eval_results/lattice.dat")
print(f"{length_error=}, {angle_error=}")
```

             LightPFP       PFP          err      rel err
    -----  ----------  --------  -----------  -----------
    a         3.56479   3.56563  0.000834876  0.000234146
    b         3.56479   3.56563  0.000834861  0.000234141
    c         3.56479   3.56563  0.000834856  0.00023414
    alpha    90        90        3.10915e-07  3.45461e-09
    beta     90        90        2.94716e-08  3.27462e-10
    gamma    90        90        5.81261e-08  6.45845e-10
    length_error=0.000234142435251996, angle_error=1.328374281683864e-07




### 状態方程式（EOS）

状態方程式（EOS）評価は、エネルギー-体積関係の予測におけるLightPFPの精度を検証します。このテストでは、LightPFPとPFPの両方を使用してNi3AlのEOSを計算します。体積弾性率（B）と平衡体積（V0）は、体積-エネルギー曲線のフィッティングによって導出され、これらの結果を二つのモデル間で比較します。LightPFPとPFPから得られた体積-エネルギー曲線の視覚的な比較は、`eval_results/eos.png`で提供されています。


```
evaluate_eos(atoms, estimator_fn_pfp, estimator_fn_light_pfp, save_fig="eval_results/eos.png")
```


      0%|          | 0/11 [00:00&lt;?, ?it/s]



      0%|          | 0/11 [00:00&lt;?, ?it/s]


          LightPFP       PFP         err      rel err
    --  ----------  --------  ----------  -----------
    B     182.814   181.012   1.80223     0.00995639
    V0     45.3276   45.3247  0.00288866  6.37326e-05



```
display.Image("eval_results/eos.png")
```




    
![png](output_31_0.png)
    





### 弾性（Elastic）

弾性評価は、材料の機械的特性の予測におけるLightPFPの精度を評価します。このテストでは、LightPFPとPFPの両方を使用してNi3Alの弾性テンソルと様々な弾性特性を計算します。弾性テンソルの成分と、それから導出される弾性特性（体積弾性率、せん断弾性率、ヤング率、ポアソン比を含む）を、二つのモデル間で比較します。


```
elastic_tensor_err, elasticity_err = evaluate_elastic(
    atoms, 
    estimator_fn_pfp, 
    estimator_fn_light_pfp, 
    save_result="eval_results/elastic.dat"
)
print(f"{elastic_tensor_err=}, {elasticity_err=}")
```

                     LightPFP         PFP        err    rel err
    -------------  ----------  ----------  ---------  ---------
    C11            252.468     242.548     9.92056    0.0409015
    C12            157.965     148.527     9.43772    0.0635422
    C13            157.964     148.528     9.43676    0.0635354
    C22            252.468     242.549     9.91895    0.0408946
    C23            157.965     148.529     9.43593    0.0635294
    C33            252.468     242.551     9.91691    0.0408859
    C44            124.497     128.87      4.373      0.0339334
    C55            124.497     128.87      4.37357    0.0339377
    C66            124.497     128.87      4.37337    0.0339362
    bulk_modulus   189.466     179.868     9.59747    0.0533584
    shear_modulus   84.4367     86.044     1.60731    0.0186801
    young_modulus  220.547     222.632     2.08444    0.0093627
    poisson_ratio    0.305992    0.293708  0.0122834  0.0418218
    elastic_tensor_err=0.046121801466313936, elasticity_err=0.03080572055593341




### フォノン

フォノン評価は、LightPFPが物質の振動特性と関連する熱力学的特性を予測する能力を調べるものです。このテストでは、LightPFPとPFPの両方を使用して、Ni3Alのフォノン分散とフォノン状態密度（DOS）を計算します。さらに、エントロピー（S）、比熱（Cv）、内部エネルギー（U）、およびヘルムホルツ自由エネルギー（F）など、いくつかの熱物性も計算されます。両モデルの結果は、図「eval_results/phonon.png」で可視化され、比較されます。


```
frequency_error = evaluate_phonon(
    atoms, 
    0.1,
    [5,5,5], 
    [10, 10, 10],
    estimator_fn_pfp, 
    estimator_fn_light_pfp, 
    save_fig="eval_results/phonon.png"
)
print(f"{frequency_error=}")
```


      0%|          | 0/25 [00:00&lt;?, ?it/s]



      0%|          | 0/25 [00:00&lt;?, ?it/s]


    frequency_error=0.31667089679183064



```
# Display phonon dispersion, DOS and thermochemistry properties
display.Image('eval_results/phonon.png')
```




    
![png](output_36_0.png)
    





### 表面

このテストでは、Ni3Alの低いMiller指数（最大指数≤2）の表面のポテンシャルエネルギーを、LightPFPとPFPの両方を使用して計算します。両モデルの結果を比較することにより、LightPFPが表面のエネルギーを予測する際の正確さを把握することができます。比較結果は、`eval_results/surface.png`にも可視化されています。


```
evaluate_surface(
    atoms, 
    2, 
    estimator_fn_pfp, 
    estimator_fn_light_pfp, 
    save_result="eval_results/surface.dat", 
    save_fig="eval_results/surface.png"
)
```


      0%|          | 0/6 [00:00&lt;?, ?it/s]



      0%|          | 0/6 [00:00&lt;?, ?it/s]


               LightPFP       PFP        err     rel err
    -------  ----------  --------  ---------  ----------
    (1 1 1)     5.16246   5.1087   0.0537582  0.0105229
    (1 1 0)     5.06755   5.03862  0.0289271  0.00574107
    (1 0 0)     3.4492    3.42169  0.0275142  0.00804111
    (2 2 1)    10.2329   10.094    0.138834   0.0137541
    (2 1 1)     9.30586   9.24766  0.0582024  0.00629375
    (2 1 0)     8.41818   8.48375  0.0655632  0.0077281





    0.008680161767822795




```
# Display potential energy of surface structures
display.Image('eval_results/surface.png')
```




    
![png](output_39_0.png)
    





### 空孔

この評価では、初期のNi3Al構造から各原子を系統的に取り除くことで空孔構造が生成され、LightPFPとPFPの両方を使用してこれらの欠陥構造のポテンシャルエネルギーが計算されます。両モデルの結果を比較することで、LightPFPが空孔形成エネルギーを予測する際の正確性を把握することができます。この比較結果は、`eval_results/vacancy.png`に視覚的に表現されています。


```
evaluate_vacancy(
    atoms, 
    [3,3,3], 
    estimator_fn_pfp, 
    estimator_fn_light_pfp, 
    save_result="eval_results/vacancy.dat", 
    save_fig="eval_results/vacancy.png"
)
```


      0%|          | 0/4 [00:00&lt;?, ?it/s]



      0%|          | 0/4 [00:00&lt;?, ?it/s]


           LightPFP      PFP        err    rel err
    ---  ----------  -------  ---------  ---------
    Al0     7.07991  7.16092  0.0810115  0.011313
    Ni1     6.39162  6.47411  0.0824961  0.0127425
    Ni2     6.39162  6.47407  0.0824554  0.0127362
    Ni3     6.39162  6.47413  0.0825117  0.0127448





    0.012384134779259258




```
# Display potential energy of vacancy structures
display.Image('eval_results/vacancy.png')
```




    
![png](output_42_0.png)
    





### MD

このテストでは、LightPFPとPFPの両方を使用して、Ni3Alの500.0 Kおよび1000.0 KにおけるMDシミュレーションが行われます。得られたMDトラジェクトリから、両モデル間でいくつかの基本的な特性を測定し、比較します。これらの特性には以下が含まれます：

* ポテンシャルエネルギーとスーパーセルの体積は、`eval_results/md.png`で可視化されます。
* 平均二乗変位は、`eval_results/msd.png`にプロットされます。
* 動径分布関数は、`eval_results/rdf.png`に表示されます。


```
evaluate_md(
    atoms, 
    supercell=[5, 5, 5], 
    temperatures=[500.0, 1000.0], 
    steps=[2000, 2000],
    pressures=None,
    save_dir="eval_results",
    reuse_traj=True,
    estimator_fn_pfp=estimator_fn_pfp,
    estimator_fn_light_pfp=estimator_fn_light_pfp,
)
```


```
# Display potential energy and volume
display.Image('eval_results/md.png')
```




    
![png](output_45_0.png)
    




```
# Display mean squared displacement
display.Image('eval_results/msd.png')
```




    
![png](output_46_0.png)
    




```
# Display radial distribution function
display.Image('eval_results/rdf.png')
```




    
![png](output_47_0.png)
    





### コマンドラインツール（CLI）

シェルコマンドからすべての評価タスクを実行するためのCLIも提供しています。

```
light-pfp-evaluate  -c assets/how_to_use_lightpfp/evaluate.json
```

詳細については、マニュアルを参照してください。
