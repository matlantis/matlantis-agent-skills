# LightPFP autogen (active learning method) example
本チュートリアルではLightPFPで提供されている light-pfp-autogen パッケージの使用方法について、実際の例を用いて解説します。
light-pfp-autogen パッケージは、 active learning と呼ばれる手法を用いて、目的のMDシミュレーションに適したLightPFPモデルに必要な学習データを自動で収集します。
(active learning の手法の詳細は Matlantis上のドキュメント `LightPFPについて` を参照ください。)

本チュートリアルでは、 `how_to_use_lightpfp` の example で用いられている Ni3Al の系でMDシミュレーションを行うようなケースを想定し、 light-pfp-autogen を用いたLightPFPモデルの構築を以下の3ステップで実践します。


**1. 初期データセットの準備と ActiveLearning オブジェクトの初期化**

active learning を開始するために、まずは初期データセットの準備が必要です。
初期データセットを元に `ActiveLearning` オブジェクトを初期化し、最初のLightPFPモデルの構築を行います。


**2. ActiveLearning オブジェクトを用いたデータ収集**

1で準備した `ActiveLearning` オブジェクトを用い、 active learning 手法に基づいてデータセットの収集を行います。
ここでは実際に行いたいMDシミュレーションと類似のスクリプトを使用することで、目的のMDシミュレーションに必要なデータセットを収集します。


**3. 収集したデータセットを用いての最終学習**

最後に、2で収集したデータセットを用いてLightPFPモデルの最終学習を行います。
active learning の途中でもLightPFPモデルは生成されますが、これらはおおよその精度を計算するために簡易的に学習されたモデルであり、実際にMDシミュレーションで使用するモデルとしての安定性は十分でないことがあります。
よって、実際にMDシミュレーションで用いるLightPFPモデルを構築するためには、収集したデータを用いて最終的な学習を行うことが推奨されます。


このチュートリアルを通して実行することで、 light-pfp-autogen を用いた一連のフローを理解することができます。
以下のセルを実行して必要なライブラリをインストールします。


```python
# Install necessary libraries
%pip install -U light-pfp-client light-pfp-data light-pfp-evaluate light-pfp-autogen matlantis_features
```


## 1: 初期データセットの準備と ActiveLearning オブジェクトの初期化

まずは active learning を開始するための初期データセットを準備します。
ここでは `how_to_use_lightpfp` の example で用いたデータセットの一部を使用します。
以下のコマンドでデータセットを収集できますが、すでに作成したデータセットを `assets/autogen_example/init.h5` に保存しているため、この手順は省略できます。
データセットの詳細については `assets/autogen_example/init_data.json` を参照ください。

```
PFP_DEFAULT_PRIORITY=1 dataset-generation -c assets/autogen_example/init_data.json --log-level info --progress-bar "True" -j 8
```

初期データセットを用いて active learning の設定を行います。


### 注意
本チュートリアルでは処理時間短縮のため、 `sample_config`, `train_config` のパラメータを調整しています。
実際の使用では、これらの設定を変更せずデフォルト値を使用することを推奨します。


```python
from light_pfp_autogen.active_learning import ActiveLearning
from light_pfp_autogen.config import (
    ActiveLearningConfig, 
    MTPConfig,
    SampleConfig,
    TrainConfig,
)

active_learning_config = ActiveLearningConfig(
    task_name="Ni3Al-autogen-example",
    init_dataset=["assets/autogen_example/init.h5"],
    pfp_model_version="v8.0.0",
    pfp_calc_mode="PBE",
    training_time=0.5,
    # NOTE: We don't use pre-trained models in this tutorial to showcase the benefit of the data sampled by active learning.
    # For a real case, we recommend using the default pre-trained model or others like "ALL_ELEMENTS_LARGE_6".
    train_config=TrainConfig(
        mtp_config=MTPConfig(pretrained_model="")
    ),
    # NOTE: This value is set to a quite low value just to show the benefit of active learning in this tutorial.
    # For a real case, we recommend using the default value or a value around 3.0.
    sample_config=SampleConfig(
        dE_min_coef=0.2
    )
)

active_learning = ActiveLearning(active_learning_config)
```



`initialize` メソッドを呼ぶことで `ActiveLearning` オブジェクトの初期化を行います。
初期化処理の中で最初のLightPFPモデルの学習も行うため、以下のセルの実行には30分程度かかることが想定されます。


```python
active_learning.initialize()
```



## 2: ActiveLearning オブジェクトを用いたデータ収集

`ActiveLearning` オブジェクトと、目的のMDシミュレーションに対応するスクリプトを用いてデータ収集を行います。
`how_to_use_lightpfp` の例では NPT を用いて 500K, 1000K でMDシミュレーションを実行した際の結果を評価しているため、同様の設定でMDシミュレーションを実行してデータを収集します。
データセットを収集するMDシミュレーションは`DataCollectionContext` のコンテキスト内で実行します。
なお、 `Atoms` オブジェクトの `calc` 属性は `DataCollectionContext` 内で自動的に設定されるため、MDシミュレーションのスクリプト内で明示的に設定する必要はありません。

また本チュートリアルでは簡易的に、イテレーション数は2回と少なめに設定しています。
`update` メソッドでは収集したデータを用いて次のイテレーションで使用するLightPFPモデルを学習するため、以下のセルの実行には1時間程度かかります。


### Tips
実際にLightPFPを適用するMDシミュレーションは大規模な系であることが想定されます。
active learningでは内部的にPFPも使用するため、メモリや計算時間の観点から、大規模な系のMDシミュレーションをそのまま用いるのではなく、系を小さくしたMDシミュレーションを実行することが推奨されます。


```python
from ase.io import read
from ase.md.npt import NPT
from ase import units

from light_pfp_autogen.context import DataCollectionContext

num_iterations = 2
for _ in range(active_learning.iter, num_iterations):
    print(f"Iteration: {active_learning.iter} starts!")

    atoms = read("assets/autogen_example/AlNi3_mp-2593_conventional_standard.cif") * (3,3,3)

    # Collect training dataset in MD task
    # NOTE: We set the MD steps to 2000, which is quite small for the actual use cases.
    # The intention is just to shorten the execution time, so for the real case, please prolong this value
    # according to your desired simulations.
    md0 = NPT(atoms, timestep=units.fs, temperature_K=500)
    with DataCollectionContext(md=md0, interval=100, max_samples=100):
        md0.run(steps=2000)

    md1 = NPT(atoms, timestep=units.fs, temperature_K=1000)
    with DataCollectionContext(md=md1, interval=100, max_samples=100):
        md1.run(steps=2000)

    # Update the LightPFP model with new dataset
    active_learning.update()
```



## 3: 収集したデータセットを用いてのモデルの学習

ステップ2で収集したデータセットを用いた最終学習を行います。
2回のイテレーションを実行したため、 `autogen_workdir` の中に `dataset_0.h5`, `dataset_1.h5` というデータセットが生成されています。
なお、学習に使用できるデータセットは `ActiveLearning` オブジェクトの `datasets_list` という属性から取得できます。

以下のセルを実行して最終学習を行います。
これでlight-pfp-autogenを用いたLightPFPモデル構築の一連の流れが完了となります。


### 注意
本チュートリアルでは、 active learning で収集したデータセットの効果を示すため、事前学習モデルを用いずにLightPFPモデルを構築しています。
実際の使用では、 `ALL_ELEMENTS_LARGE_6` などの事前学習モデルを用いることでより高精度なLightPFPモデルが構築できます。



```python
from light_pfp_autogen.utils import submit_training_job, check_training_job_status

train_config_dict = {
    "common_config": {
        "total_epoch": 200,
    },
    # NOTE: We don't use pre-trained models in this tutorial to showcase the benefit of the data sampled by active learning.
    # For a real case, we recommend using the default pre-trained model or others like "ALL_ELEMENTS_LARGE_6".
    "mtp_config": {
        "pretrained_model": ""
    },
}

training_config = TrainConfig.from_dict(
    train_config_dict
)

final_model_id = submit_training_job(
    training_config,
    active_learning.datasets_list,
    "Ni3Al-autogen-example-final"
)

status = check_training_job_status(final_model_id)
print(f"Training job {final_model_id} status: {status}")
```



# Appendix
## active learning の効果検証
active learning で生成されたデータセットの有効性を確認します。
PFPとLightPFPの結果を比較し、本チュートリアルで行ったLightPFPモデル構築の一連の流れの効果を検証します。
PFPとLightPFPの結果比較には light-pfp-evaluate パッケージが便利です。

## Tips
active learning の効果は、 Model Library に表示されるMAE値だけでは適切に評価できません。
各イテレーションで学習に使用するデータが異なるため、実際のMDシミュレーションでの挙動を比較することが重要です。


```python
# Please update these values according to your model ids
# You can find them in Model Library
initial_model_id = "initial_model_id"
final_model_id = "final_model_id"

from ase.io import read
from light_pfp_client.estimator_fn import light_pfp_estimator_fn
from light_pfp_evaluate import evaluate_md
from matlantis_features.utils.calculators.pfp_api_calculator import pfp_estimator_fn
from pfp_api_client.pfp.estimator import EstimatorCalcMode


estimator_fn_pfp = pfp_estimator_fn(model_version="v8.0.0", calc_mode=EstimatorCalcMode.PBE)
init_estimator_fn_light_pfp = light_pfp_estimator_fn(model_id=initial_model_id)
final_estimator_fn_light_pfp = light_pfp_estimator_fn(model_id=final_model_id)

atoms = read("assets/autogen_example/AlNi3_mp-2593_conventional_standard.cif")
# Compare initial_model accuracy with PFP
evaluate_md(
    atoms,
    supercell=[5, 5, 5], 
    temperatures=[1000.0], 
    steps=[1000],
    pressures=None,
    save_dir="eval_init_results",
    reuse_traj=True,
    estimator_fn_pfp=estimator_fn_pfp,
    estimator_fn_light_pfp=init_estimator_fn_light_pfp,
)

# Compare final model with PFP
evaluate_md(
    atoms,
    supercell=[5, 5, 5], 
    temperatures=[1000.0], 
    steps=[1000],
    pressures=None,
    save_dir="eval_final_results",
    reuse_traj=True,
    estimator_fn_pfp=estimator_fn_pfp,
    estimator_fn_light_pfp=final_estimator_fn_light_pfp,
)
```


```python
from IPython import display
display.Image('eval_init_results/md.png')
```




    
![png](output_12_0.png)
    




```python
from IPython import display
display.Image('eval_final_results/md.png')
```




    
![png](output_13_0.png)
    





上記の結果から、初期のLightPFPモデルと比較して、active learning で収集したデータを用いて構築したLightPFPモデルの方が高精度であることが分かります。

なお、このチュートリアルはシンプルなMDシミュレーションを用いて評価を行っています。
実際の系に適用する際は、Matlantis内で公開されている以下のような例を参照ください。

* `Amorphous Silicon` (`Matlantis Examples`)
* `light_pfp_SiO2_HF_etching` (`Matlantis Contrib`)
* `light_pfp_SiO2_P2O5_Al2O3_Na2O_glass` (`Matlantis Contrib`)
* `light_pfp_Fe2O3_lubricant` (`Matlantis Contrib`)
* `light_pfp_polymer_ionic_liquid` (`Matlantis Contrib`)



## active learning 中の生成物
light-pfp-autogen のスクリプトを実行すると、各イテレーションごとに `autogen_workdir` の中に以下のファイルが生成されます。
`${i}` はイテレーションのインデックスを表します。

* `dataset_${i}.h5`
* `md_log_${i}`
* `train_log_${i}`

`dataset_${i}.h5` には、各イテレーションで生成されたデータセットが含まれています。

`md_log_${i}` にはデータのサンプリングに使用した計算条件や、どのデータが収集されたのかについての情報が含まれています。
このファイルの冒頭を数行例示すると、以下のようになります。

```
# pfp_model_version v8.0.0
# pfp_calc_model PBE
# model_id xxx 
# dE_min 0.00010776519775390625
# dE_max None
# dF_min None
# dF_max None
# dS_min None
# dS_max None
# Data Collection Section Started
    0     0.000K   -4.9187946   -4.9185892    0.0002055    0.0000016    0.0004036    7.46673g/cm^3 True
  100     0.000K   -4.9187946   -4.9185893    0.0002053    0.0000022    0.0004036    7.46673g/cm^3 True
```

冒頭の数行はデータのサンプリングに使用した計算条件を示しており、 `Data Collection Section Started` 以降には実際に収集されたデータの情報が含まれています。

各カラムはそれぞれ以下の情報に相当します:
1. MDステップ
2. 温度
3. LightPFPのエネルギー
4. PFPのエネルギー
5. エネルギー差
6. 力の差の最大値
7. 応力の差の最大値
8. 密度
9. 訓練データとして選択されたかどうか

これらの情報を元に、 active learning の計算条件をどうチューニングするかの詳細については、 `LightPFPについて` のドキュメント `6.4 FAQ` のセクションを参照ください。


`train_log_${i}` には各イテレーションでのLightPFPモデルの学習についての情報が含まれています。
このログを確認することで、各イテレーションで生成されたLightPFPモデルの model_id および学習が成功したかどうかを確認できます。



## active learning の再開

light-pfp-autogen では前回のイテレーションから active learning を再開する機能を提供しています。
これは以下のような際に便利です。

* active learning が意図せずに停止した場合。
* 指定したイテレーション後にLightPFPモデルが十分に正確ではないと判明した場合。 num_iterations をより大きく設定して、以前の終了地点から追加のイテレーションを続行するためにタスクを再開できます。
* MDスクリプトとサンプリング基準を調整する必要がある場合。active learning を手動で停止し、MDスクリプトと設定を変更してから再開できます。

再開の機能をどのような際に使うべきかの詳細については、 `LightPFPについて` の `6.3 アクティブラーニングの詳細` および `6.4 FAQ` セクションを参照ください。
再開をする際は、active learning の作業ディレクトリの中で以下の手順を実施します。

1. `ActiveLearning` オブジェクトの `initialize` メソッドを呼び出す。これによって `iter` 属性の値が前回の状態まで復元され、 active learning を再開する準備ができます。
2. `iter` 属性を使用してイテレーションを再開し、 `DataCollectionContext` と `update` メソッドを使用してデータを収集する。

例えば、先ほどの active learning から設定を変更し、追加で3回のイテレーションを実行する場合は以下のようなコードを実行します。


```python
active_learning_config = ActiveLearningConfig(
    task_name="Ni3Al-autogen-example",
    init_dataset=["assets/autogen_example/init.h5"],
    pfp_model_version="v8.0.0",
    pfp_calc_mode="PBE",
    training_time=0.5,
    train_config=TrainConfig(
        mtp_config=MTPConfig(pretrained_model="")
    ),
    # Modified the SampleConfig
    sample_config=SampleConfig(
        dE_min_coef=0.5
    )
)

active_learning = ActiveLearning(active_learning_config)

# Please call initialize method to prepare for the resume of active learning
active_learning.initialize()
print(f"Start active learning from {active_learning.iter} iteration")

num_iterations += 3
for i in range(active_learning.iter, num_iterations):
    atoms = read("assets/autogen_example/AlNi3_mp-2593_conventional_standard.cif") * (3,3,3)

    md2 = NPT(atoms, timestep=units.fs, temperature_K=1500)
    with DataCollectionContext(md=md2, interval=100, max_samples=100):
        md2.run(steps=2000)

    # Update the LightPFP model with new dataset
    active_learning.update()
```
