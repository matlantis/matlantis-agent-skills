# Pythonの関数APIを使ったデータセットの生成

LightPFPのチュートリアルでは、JSONコントロールファイルを使用して`dataset-generation` CLIでトレーニングデータセットを生成する方法を紹介しました。しかし、JSONコントロールファイルを使わずに、それぞれのサンプラーを直接呼び出してデータを生成することも可能です。

**重要な注意点**: この例は、各種サンプリング関数の使用方法を示すことを主な目的としています。生成されるデータセットは比較的小規模であり、信頼性の高い LightPFP モデルの学習には不十分な場合があります。実際の使用においては、分子動力学シミュレーションの sampling_steps を増やす、あるいは初期構造の生成数を増やすなど、パラメータを調整し、より網羅的で頑健な学習データセットを作成してください。



## 1. 初期化

* 必要なライブラリを準備する。


```python
!pip install 'light-pfp-data&gt;=1.3.0'
```


```python
from ase.io import read
from ase.visualize import view
from pfp_api_client import ASECalculator, Estimator, EstimatorCalcMode
from light_pfp_data.utils.dataset import H5DatasetWriter
```



* H5DFデータセットを初期化する


```python
dataset = H5DatasetWriter("AlNi3.h5", mode="overwrite")
```



* データセット生成の初期構造を取得する
なお、AlNi3(mp-2593)の構造ファイルは [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) の下で Materials Project から取得しました。



```python
atoms = read("assets/how_to_generate_dataset/AlNi3_mp-2593_conventional_standard.cif")
```


```python
view(atoms*(3,3,3), viewer="ngl")
```



* PFP calculator をロードする。


```python
estimator = Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.PBE)
calc = ASECalculator(estimator)
```



## 2. データセットのサンプリング

訓練データセットを生成するための基本的なサンプリング手法は2つ存在します。1つ目の手法は、既知の初期構造から出発し、分子動力学シミュレーションなどの手法を用いて材料を系統的に変化させて訓練データを作成する方法です。2つ目の手法は、シミュレーションボックス内の元素組成のみに基づいてランダムな構造を生成する方法です。それぞれの手法には明確な利点があります。1つ目の手法は、対象材料によく似た訓練データを生成することでモデルの精度を向上させ、2つ目の手法は、より多様な構造配置を生成することでモデルの頑健性を高めます——特に極端な条件や特殊な条件下での信頼性の高い予測において重要です。



### 2.1. 特定構造からのデータセットサンプリング手法

以下の関数を使用することで、特定の入力構造から異なるタイプの訓練構造を生成できます。これらの関数の使用例は後述します。

* sample_md
* sample_compress
* sample_deformed
* sample_displaced
* sample_vacancy
* sample_surface
* sample_substitution
* sample_rattle_and_relax

サンプリング関数は、`input_structure`、`calculator`、`dataset`、`show_progress_bar`、`num_threads`などの共通パラメータを受け付けます。

* input_structure: トレーニングデータ生成の出発点となる初期の原子構造。
* calculator: トレーニング構造を生成するために使用するPFPカルキュレーター。
* dataset: H5DFデータセット。
* show_progress_bar: プログレスバーの表示を有効/無効にする。
* num_threads: 構造生成時に使用する並列スレッド数。

これらのパラメータは、さまざまなサンプリング関数で共有されています。また、各サンプリング関数には、その特定のサンプリング手法に合わせた固有の引数セットも用意されています。



#### 2.1.1 MDから構造を生成する

* sample_md関数を使って、1つまたは複数のMD軌道からトレーニング構造を収集できます。 この関数の説明については `help(sample_md)` を実行することで確認できます。 

* 以下のコードでは、入力構造から2つのMD計算を開始する。1つは300K、もう1つは500Kで、両方とも1000ステップ実行し、100ステップごとに1つずつトレーニング構造を収集する(2x10 = 20個のトレーニング構造)。

* `sample_md`は最も汎用性の高いサンプリング方法で、ほぼすべての対象材料に適用できます。豊富で多様なトレーニング構造を提供するため、ほとんどのシナリオで基本的な選択肢となります。


```python
from light_pfp_data.sample import sample_md

help(sample_md)
```

    Help on function sample_md in module light_pfp_data.sample.crystal:
    
    sample_md(input_structure: ase.atoms.Atoms, calculator: ase.calculators.calculator.Calculator, dataset: light_pfp_data.utils.dataset.H5DatasetWriter, supercell: Tuple[int, int, int] = (3, 3, 3), sampling_temp: Union[List[float], float] = 500.0, sampling_steps: Union[List[int], int] = 10000, sampling_interval: Union[List[int], int] = 100, sampling_pressure: Union[List[float], float, NoneType] = None, timestep: float = 1.0, structure_type: int = 5, ensemble: str = 'npt', show_progress_bar: bool = False, executor: Optional[concurrent.futures.thread.ThreadPoolExecutor] = None, num_threads: int = 8, max_retries: int = 0, pbar: Optional[tqdm.auto.tqdm] = None) -&gt; List[concurrent.futures._base.Future]
        Generate structures by running MD simulation.
        
        Args:
            input_structure (Atoms): The input structure.
            calculator (Calculator): The ASE calculator for energy calculation.
            dataset (H5DatasetWriter): The dataset.
            supercell (Tuple[int, int, int], optional):
                Expand the input structure into supercell. Defaults to (3, 3, 3).
            sampling_temp (Union[List[float], float], optional):
                MD simulation temperature, in K. Defaults to 500.0.
            sampling_steps (Union[List[int], int], optional): MD steps. Defaults to 10000.
            sampling_interval (Union[List[int], int], optional):
                Collect training structure every this steps. Defaults to 100.
            timestep (float, optional):
                Time step in fs. Defaults to 1.0.
            structure_type (int, optional): The structure type. Defaults to 5.
            ensemble (str, optional):
                The MD ensemble. Supported "nvt" and "npt. Defaults to "npt".
            show_progress_bar (bool, optional):
                Show progress bar. Defaults to False.
            executor (ThreadPoolExecutor, optional):
                Thread pool executor parallel calculation. Defaults to None.
            num_threads (int, optional):
                Max number of threads to use for executor if no executor is passed. Defaults to 8.
            max_retries (int, optional):
                Max retries for PFP calculation. Defaults to 0.
    



```python
sample_md(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    supercell=(3, 3, 3),
    sampling_temp=[300.0, 500.0],
    sampling_steps=[1000, 1000],
    sampling_interval=[100, 100],
    ensemble="npt",
    show_progress_bar=True,
    num_threads=8,
)
```


    MD:   0%|          | 0/20 [00:00&lt;?, ?it/s]





    []





#### 2.1.2 入力構造を均一に圧縮および伸張して構造を取得します

* 関数「sample_compress」は、入力構造の格子を均一に圧縮および伸張し、異なる体積のトレーニング構造を生成できます。さらに、構造最適化と分子動力学を実行して、より多くの構造を取得することもできます。追加してより多くの構造を取得したい場合は、生成された構造から構造最適化やMDを実行することも可能です。 

* 次のコード ブロックでは、次の操作を行います。

    * (1) 入力構造の格子定数を -5%、0%、5% に変更します。(3 つのトレーニング構造)
    
    * (2) 上記の 3 つの構造から固定ボリュームで構造最適化を実行します (3 つのトレーニング構造)
    
    * (3) 上記の (最適化された) 構造のそれぞれから 2 つの MD タスクを開始します。1 つは 300K で、もう 1 つは 500K で実行します。両方の MD タスクは 1000 ステップ実行され、100 MD ステップごとに 1 つのトレーニング構造が収集されます。(3x2x10 = 60 のトレーニング構造)。
    
* `sample_compress`は密度の精度を向上させるのに有用です。異なる体積を持つトレーニング構造を収集し、体積とポテンシャルエネルギーの関係をより明確に学習できるようにします。この方法は、有機液体のシミュレーションなどの材料に特に効果的です。


```python
from light_pfp_data.sample import sample_compress

sample_compress(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    supercell=(1,1,1),
    min_scale=0.95,
    max_scale=1.05,
    interval=0.05,
    opt=True,
    md=True,
    sampling_temp=[300.0, 500.0],
    sampling_steps=[1000, 1000],
    sampling_interval=[100, 100],
    show_progress_bar=True,
    num_threads=8,
)
```


    compress:   0%|          | 0/6 [00:00&lt;?, ?it/s]





    []





#### 2.3 変形構造を生成する
* sample_deformed関数を使って、入力構造の単位格子に歪みテンソルを適用することで、トレーニング構造を生成できます。

* 6つの独立成分(xx、yy、zz、xy、yz、xz)からなる歪みテンソルを適用して、初期構造の格子形状を変形させます。

* 例えば、xy成分に0.01の歪みを加えた場合、strainテンソルは以下のようになります:
$$
\left[
\begin{matrix}
0.0 &amp; 0.01 &amp; 0.0 \\
0.01 &amp; 0.0 &amp; 0.0 \\
0.0 &amp; 0.0 &amp; 0.0 \\
\end{matrix}
\right]
$$

* 以下のコードブロックでは、6つの独立成分に-0.01、-0.005、0.0、0.005、0.01の歪みを適用して30個のトレーニング構造を生成します。

* その後、格子形状を固定して構造最適化を行い、さらに30個のトレーニング構造を生成します。

* `sample_deformed`方法は、ヤング率などの機械的特性を正確に再現するために設計されています。機械的特性が対象材料にとって重要な場合、例えば合金の場合などに、このメソッドを有効にするべきです。


```python
from light_pfp_data.sample import sample_deformed

sample_deformed(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    supercell=(2,2,2),
    min_strain=-0.01,
    max_strain=0.01,
    interval=0.005,
    opt=True,
    show_progress_bar=True,
    num_threads=8,
)
```


    deformed:   0%|          | 0/60 [00:00&lt;?, ?it/s]





    []





#### 2.1.4 1つの原子を移動させて構造を生成します。

* 関数 `sample_displaced` は、x、y、または z 軸に沿って平衡位置から 1 つの原子を移動させた構造を生成できます。

* 以下のコードでは、ランダムに1つの原子を選び、x、y、z方向に0.1 angstromだけ変位させた20個の構造を生成します。

* `sample_displaced`方法は、LightPFPモデルが振動特性をより効果的に学習するのに役立ちます。


```python
from light_pfp_data.sample import sample_displaced

sample_displaced(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    delta=0.1,
    supercell=(3,3,3),
    n_sample=20,
    show_progress_bar=True,
    num_threads=8,
)
```


    displacements:   0%|          | 0/20 [00:00&lt;?, ?it/s]





    []





#### 2.1.5 rattle法を用いたランダム構造の生成
* 関数 `sample_rattle` は、初期構造にランダムな原子変位を加えることによって、多様なトレーニングデータを生成できます。
* 前述の `sample_displaced` 関数とは異なり、この方法では、すべての原子に対して同時に変位を適用し、その変位の方向および大きさはランダムです。
* 以下のコードでは、初期構造から出発して20個のトレーニング構造が生成されます。各原子のx、y、z方向の変位は、標準偏差が0.1の正規分布からサンプルされます。
* 大きな原子変位のため、非常に高いポテンシャルエネルギーや極めて大きな原子間力を持つ非現実的な構造が生成される可能性があります。このようなノイズの多いデータをフィルタリングするために、　`max_forces`　パラメータを利用することができます。
* `sample_rattle`　法は、有機分子の安定性を向上させるのに特に有益です。


```python
from light_pfp_data.sample import sample_rattle

sample_rattle(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    stdev=0.1,
    supercell=(3,3,3),
    n_sample=20,
    max_forces=50.0,
    show_progress_bar=True,
    num_threads=8,
)
```


    rattle:   0%|          | 0/20 [00:00&lt;?, ?it/s]





    []





#### 2.1.6 欠陥構造を生成する
* sample_vacancy関数を使って、いくつかの原子をランダムに削除することで、欠陥構造を生成できます。

* 以下のコードブロックでは、

    * (1) ランダムに1つの原子を削除して3個の欠陥構造を生成する。
    
    * (2) 上記3つの構造について構造最適化を行い、さらに3個のトレーニング構造を得る。
    
    * (3) 最適化された構造それぞれから2つのMD計算を開始する。1つは300K、もう1つは500Kで、両方とも1000ステップ実行し、100ステップごとに1つずつトレーニング構造を収集する(60個のトレーニング構造)。
    
* `sample_vacancy`方法は、LightPFPモデルの空孔形成プロセスの再現能力を向上させます。この方法は、合金などの材料を扱う際に特に有用です。


```python
from light_pfp_data.sample import sample_vacancy

sample_vacancy(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    supercell=(3,3,3),
    n_vacancy=1,
    n_sample=3,
    opt=True,
    md=True,
    sampling_temp=[300.0, 500.0],
    sampling_steps=[1000, 1000],
    sampling_interval=[100, 100],
    show_progress_bar=True,
    num_threads=8,
)
```


    vacancy:   0%|          | 0/6 [00:00&lt;?, ?it/s]





    []





#### 2.1.7 表面構造を生成する

* sample_surface関数を使って、入力構造から表面構造を作成できます。

* 以下のコードブロックでは、

    * (1, 1, 1)、(1, 1, 0)、(1, 0, 0)、(2, 2, 1)、(2, 1, 1)、(2, 1, 0)面に沿って入力構造を切断して6つの表面構造を生成する。表面構造は10 angstromの固体層と10 angstromの真空層から成る。
    
    * 入力構造が高対称であることに注意。同じ物質でも、(1, 0, 0)面と(0, 1, 0)面、(0, 0, 1)面は異なる。
    
    * 上記6つの構造について構造最適化を行い、さらに6個のトレーニング構造を得る。

    * 最適化された構造それぞれから1つのMD計算を開始する。MD計算は1000ステップ実行し、100ステップごとに1つずつトレーニング構造を収集する(60個のトレーニング構造)。
    
* `sample_surface`方法は、結晶表面に関する基本的な情報を提供する簡単なアプローチです。結晶材料を扱い、表面特性がシミュレーションにとって重要な場合に、このメソッドを有効にするべきです。


```python
from light_pfp_data.sample import sample_surface

sample_surface(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    supercell=(1,1,1),
    max_index=2,
    min_slab_size=10.0,
    min_vacuum_size=10.0,
    opt=True,
    md=True,
    sampling_temp=[300.0],
    sampling_steps=[1000],
    sampling_interval=[100],
    show_progress_bar=True,
    num_threads=16,
)
```


    surface:   0%|          | 0/12 [00:00&lt;?, ?it/s]





    []





#### 2.1.8 元素置換による構造生成

* `sample_substitution`メソッドは、初期構造内の原子を指定された確率で1つまたは複数の新しい元素にランダムに置換することで、トレーニング構造を作成します。

* 以下のコードブロックでは：

    * (1) 原子の10%をCrに、5%をMoにランダムに置換することで、3つの新しいトレーニング構造を作成します。（3つのトレーニング構造）
    
    * (2) 上記の各構造について、300Kと500Kで2つのMDタスクを開始します。両MDタスクは1000ステップ実行され、100MDステップごとに1つのトレーニング構造を収集します。（3x2x10 = 60のトレーニング構造）

* `sample_substitution`方法は、特に合金、とりわけハイエントロピー合金に有用です。これは、複数（通常5つ以上）の金属元素が格子位置をランダムに占め、多数の可能な原子配列をもたらす材料を扱うように設計されています。


```python
from light_pfp_data.sample import sample_substitution

sample_substitution(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    n_sample = 3,
    elements = [24, 42],
    possibility = [0.1, 0.05],
    supercell = (4, 4, 4),
    opt=False,
    md=True,
    sampling_temp=[300.0, 500.0],
    sampling_steps=[1000, 1000],
    sampling_interval=[100, 100],
    show_progress_bar=True,
    num_threads=8,
)

```


    substitution:   0%|          | 0/3 [00:00&lt;?, ?it/s]





    []





#### 2.1.9 rattle法と緩和による構造生成

`sample_rattle_and_relax`メソッドは([Gardner 2025](https://arxiv.org/html/2506.10956v1))で導入されています。初期構造に対してランダムな原子変位を適用して新たな配置を生成し、その後緩和させる手法です。このプロセスでは、新たに生成された構造を次のサイクルの開始点として繰り返し使用します。`sample_rattle`メソッドとは異なり、この反復的なアプローチでは原子位置の複数回の修正が可能となり、より広範囲な配置空間を探索できます。

* 以下のコードブロックでは、次の手順で10個の新しい構造を生成します:
    * 初期摂動: 開始構造に対し、標準偏差0.1～0.5の正規分布からサンプリングしたランダムな原子変位を適用します。同時に、格子パラメータも標準偏差0.05の正規分布を用いてランダムに変更します。
    * 緩和処理: 摂動を加えた構造を最適化ステップで緩和させ、生成された構造のコレクションに追加します。
    * 構造選択: `rattle_temp`と`rattle_beta`パラメータに基づいて生成された構造から1つを選択します。`rattle_temp`の値が大きいほど高エネルギー構造が選択されやすくなり、`rattle_beta`の値が大きいほど最近生成された構造が優先されます。
    * 反復処理: 選択された構造を新たな開始点としてラトル法と緩和のプロセスを繰り返し、10個の構造が生成されるまで続けます。



```python
from light_pfp_data.sample import sample_rattle_and_relax

sample_rattle_and_relax(
    input_structure=atoms,
    calculator=calc,
    dataset=dataset,
    rattle_samples=10,
    rattle_sigma_range=(0.1, 0.5),
    rattle_temp=5000.0,
    rattle_beta=0.5,
    rattle_cell_sigma=0.05,
    show_progress_bar=True,
    num_threads=8,
)
```


    rattle_relax:   0%|          | 0/10 [00:00&lt;?, ?it/s]





    []





### 2.2 元素と組成からのデータセットサンプリング

前述の通り、訓練データセットは既定の開始構造を必要とせず、元素組成のみに基づいてランダムな原子構造を生成することで作成できます。これらのランダムに生成された構造は対象材料と大きく異なる場合が多いですが、極端なシミュレーション条件下では重要な利点を提供します。例えば、原子配置が元の訓練データから大きく逸脱する高温分子動力学シミュレーションにおいて、これらの多様なランダム構造はシミュレーションの失敗を防ぎ、計算の安定性を維持する上で不可欠な参照点として機能します。

以下に、ランダムな材料構造を生成するための2つの関数を示します:

* sample_random_structure
* sample_random_composition



#### 2.2.1 特定の組成から構造を生成

`sample_random_structure` メソッドは、シミュレーションボックス内のランダムな位置に指定した数の原子を配置し、物理的に不適切な重なりを避けるために近すぎる原子ペアを調整することで、ランダムな原子構造を生成します。

以下の例は、定義済みの結晶構造に依存せずに Ni₃Al の学習データセットを作成する方法を示しています：

* 初期構造生成: 充填率パラメータによって決定されるボックス寸法内で、ランダムに 30 個の Ni 原子と 10 個の Al 原子を配置します。
* 構造の複製: このランダム配置プロセスを 5 回繰り返し、5 つの学習構造を生成します。
* 分子動力学サンプリング: 各ランダム構造に対して 500 K で 1000 ステップの NVT MD シミュレーションを実行し、100 ステップごとにサンプリングを行うことで、5×10 = 50 個の追加学習構造を取得します。
* Rattle-and-relaxサンプリング: `sample_rattle_and_relax` 関数を各ランダム構造に適用し、反復的な摂動と緩和を通じてさらに 5×10 = 50 個の学習構造を生成します。
このアプローチは、構造ベースのサンプリング手法とは異なり、Ni₃Al の特定の格子構造（すなわち FCC）に制約されない学習データを生成するため、より広範な配置空間を探索できます。


```python
from light_pfp_data.sample import sample_random_structure

sample_random_structure(
    elements_list = ["Al", "Ni"],
    n_atoms=[10, 30],
    calculator=calc,
    dataset=dataset,
    n_sample=5,
    packing_fraction=0.4,
    md=True,
    sampling_temp=[500.0],
    sampling_steps=[500],
    sampling_interval=[100],
    rattle=True,
    rattle_samples=10,
    rattle_sigma_range=(0.05, 0.2),
    rattle_temp=1000.0,
    show_progress_bar=True,
    num_threads=8,
)
```


    random:   0%|          | 0/5 [00:00&lt;?, ?it/s]





    []





#### 2.2.2 元素から構造を生成する

`sample_random_composition` メソッドは、`sample_random_structure` と同様に動作しますが、元素比が事前に決定されておらず、構造生成時にランダムに変動する点が異なります。

以下の例は、組成が可変の Ni-Al 系の学習データセットを生成します。プロセスは前の `sample_random_structure` のアプローチとほぼ同じですが、初期構造生成時に特定の化学量論比を維持せず、ランダムに混合した Ni と Al の原子 40 個をシミュレーションボックス内に配置します。

長時間の分子動力学シミュレーションや極限条件下では、原子の偏析が発生し、局所的な原子組成や構造が初期の結晶配置から大きく逸脱することがあります。固定組成に依存する従来のサンプリング手法ではこのような現象を効果的に捉えられませんが、このランダム組成アプローチはシミュレーションの安定性を維持しつつ、このような困難なシナリオを扱うために必要な構造的多様性を提供します。


```python
from light_pfp_data.sample import sample_random_composition

sample_random_composition(
    elements_list = ["Al", "Ni"],
    n_atoms=40,
    calculator=calc,
    dataset=dataset,
    n_sample=5, 
    packing_fraction=0.4,
    md=True,
    sampling_temp=[500.0],
    sampling_steps=[500],
    sampling_interval=[100],
    rattle=True,
    rattle_samples=10,
    rattle_sigma_range=(0.05, 0.2),
    rattle_temp=1000.0,
    show_progress_bar=True,
    num_threads=8,
)
```


    random:   0%|          | 0/5 [00:00&lt;?, ?it/s]





    []





### 2.3 データセットを閉じる

* データセットの収集が完了したら、H5DFファイルが一つのセッションに占有されないよう、データセットを閉じてください。これにより、他のセッションでもファイルを開くことができます。


```python
dataset.h5.close()
del dataset
```



## 3. 原子構造の確認

* ユーザーはデータセットに保存されたトレーニング構造を個別に確認することができます。
* `H5DatasetWriter.get_atoms`関数からASEのatomsオブジェクトを取得できます。


```python
from h5py import File

h5 = File("AlNi3.h5")
dataset = H5DatasetWriter(h5)
```



* トレーニングポイントのキーを確認します。各トレーニング構造は固有のキーで保存されています。


```python
print(dataset.keys[:20])
```

    ['000000', '000001', '000002', '000003', '000004', '000005', '000006', '000007', '000008', '000009', '000010', '000011', '000012', '000013', '000014', '000015', '000016', '000017', '000018', '000019']




* `key="000250"`のトレーニング構造のASE atomsを取得します。
* `add_calc=True`の場合、エネルギー、力、応力の情報がASE出力atomsに付加されます。


```python
atoms = dataset.get_atoms(key="000250", add_calc=True)
view(atoms, viewer="ngl")
```



* データセットを閉じます。


```python
dataset.h5.close()
del dataset
```



## 4. データセットの統計

* `dataset_dist_analysis`関数を使って、データセット内の位置エネルギー、力、応力の分布を確認できます。


```python
from PIL import Image
from h5py import File
from light_pfp_data.utils.dataset import dataset_dist_analysis
```


```python
dataset_dist_analysis(File("AlNi3.h5"), "AlNi3_dist.jpg")
```

        Total number of structures: 545
        Total number of atoms: 40790
        Elements: Al Cr Ni Mo
        MinEnergy: -4.9219 eV/atom
        MaxEnergy: -0.9408 eV/atom
        MaxForces: 65.979 eV/angstrom
        MinStress: -58.8070 GPa
        MaxStress: 23.0144 GPa



```python
img = Image.open("AlNi3_dist.jpg")
display(img)
```


    
![png](output_50_0.png)
    




## 5. データセットの再計算

* 場合によっては、データサンプリングの全プロセスを繰り返さずに、データセットのラベル(エネルギー、力、応力)を別のPFPモデルバージョンや計算モードで更新したい場合があります。そのような場合は、データセットの`recalculation`メソッドを使用できます。

* 新しいモデルバージョンと計算モードでサンプリング関数(MDシミュレーションなど)を再実行するよりも、サンプリングプロセスをスキップできるため、この方法の方が高速です。

* データセットの再計算時は、既存のラベルをin-placeで置換するため、`mode="append"`でデータセットをロードする必要があります。


```python
from light_pfp_data.utils.dataset import H5DatasetWriter
```


```python
dataset = H5DatasetWriter("AlNi3.h5", mode="append")
```

    CAUTION!!! Data will be appended to 'AlNi3.h5' (using 'append' mode).
    This operation is irreversible! Are you sure you want to continue? [Y/n]  Y



```python
dataset.recalculate(
    model_version = "v8.0.0",
    calc_mode = "PBE_U",
    show_progress_bar = True,
    num_threads = 8
)
```


    recalc:   0%|          | 0/545 [00:00&lt;?, ?it/s]





    []




```python
dataset.h5.close()
del dataset
```
