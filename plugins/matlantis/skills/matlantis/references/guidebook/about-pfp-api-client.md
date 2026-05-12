# pfp-api-clientとは ​

Matlantisでは、第一原理計算並の精度のエネルギー予測を高速に行うことを可能とする機械学習ポテンシャル、 PFP(PreFerred Potential)を提供しています。 `pfp-api-client` とは、Python を用いて PFPを利用するためのライブラリです。

以下では、 `pfp-api-client` の基本的な機能を紹介します。より詳細な仕様に関しては [PFP API Client Reference](/api/docs/ja/pfp-api-client/latest/index.html)、 [エラーの一覧](/api/docs/ja/matlantis-guidebook/errors.html)、 およびTutorialsのHow to use PFPをご参照ください。 これらはLauncherからも開くことができます。

## pfp-api-clientのセットアップ ​

`pfp-api-client` を最新版に更新しておきます。notebook上で以下のコードを実行してください。

default
    
    
    %pip install -U pfp-api-client

ターミナル上でも実行できます。その場合、行頭の%を削除して実行してください。

インストール後、notebookのカーネルを再起動します。notebook上で必要なモジュールをインポートします。

python
    
    
    import pfp_api_client
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

## PFPを使ってみる ​

それでは、さっそく分子のエネルギーをPFPを用いて予測してみましょう。

まず、原子シミュレーション用のライブラリ Atomic Simulation Environment (ASE)から必要なメソッドをインポートします。 ASEとpfp-api-clientを連携させることで、PFPを用いた構造最適化や分子動力学法などを簡単に行うことができます。

ASEが標準で提供しているフラーレンC60の構造を読み込みます。

python
    
    
    from ase.build import molecule
    atoms = molecule("C60")

ASEを用いて計算するには、原子間のポテンシャル等を与えるCalculatorを定義する必要があります。 `pfp-api-client` を用いれば、以下のように簡単にPFPをもとにした `Calculator` を作成し `atoms` にセットすることができます。

python
    
    
    estimator = Estimator(calc_mode = "wb97xd")
    calculator = ASECalculator(estimator)
    atoms.calc = calculator

1行目では、PFPによる予測値を取得するための `Estimator` クラスからインスタンスを作成します。 今回は分子を対象としているため `calc_mode` に `"wb97xd"` を設定しています。

2行目では、 `estimator` をもとにASE用のCalculatorインスタンスを作成します。 3行目では `Calculator` を `Atoms` にセットしています。

それでは、実際にこの系のpotential energyとforceを計算してみましょう。

python
    
    
    print(f"Energy: {atoms.get_potential_energy()}")
    print(f"Force: {atoms.get_forces()}")

Calculatorを経由して計算できる基本物性値と、その計算methodは次の通りです。

  * ポテンシャルエネルギー: get_potential_energy
  * 力: get_forces
  * 応力(ストレス): get_stress
  * 電荷: get_charges



次のmethodもASEに実装されていますが、PFPが対応していないため実行してもエラーが起こります。

  * 磁気モーメント: get_magnetic_moment
  * 双極子モーメント: get_dipole_moment



また、ASEによる構造最適化や分子動力学法の機能を利用することも可能です。以下では構造最適化の例を示します。

python
    
    
    from ase.optimize import BFGS
    
    opt = BFGS(atoms)
    opt.run()
    print(f"Force after structure optimization: {atoms.get_forces()}")

構造が緩和された結果、全体的にforceの絶対値が小さくなっていることが確認できます。

この構造最適化の計算においても、内部的にPFPによる原子間のポテンシャルが呼び出されています。 このように、 `pfp-api-client` を用いることで、 PFPによる予測値を用いながらASEによるさまざまな原子シミュレーションを行うことが可能です。 さらに、 [matlantis-features](./about-matlantis-features.html) を用いればより高度な物性値を手軽に求めることができます。

より詳しい使い方については、Tutorials: How to use PFPをご参照ください。

## PFPのバージョンについて ​

Matlantisではより幅広い物質に対して高精度な予測ができるよう、継続的にPFPのバージョンアップを行っています。 一方で、ユーザーが過去に行った計算を再現できるよう、過去のバージョンのPFPも提供しております。 特に理由がなければ最新版を使うことが推奨されていますが、過去のバージョンのPFPを用いて計算を行うことも可能です。

現在の `estimator` に設定されているPFPのバージョンと、 `pfp-api-client` で利用可能なバージョン一覧を確認してみます。 (PFPのバージョンは `pfp-api-client` のバージョンとは異なることにご注意ください。)

python
    
    
    estimator = Estimator()
    print("Current model version:", estimator.model_version)
    print("Available model versions:", estimator.available_models)

デフォルトでは最新のPFPバージョンが設定されています。 次のように `Estimator` のインスタンス作成時にバージョンを指定したり、 `Estimator.set_model_version()` 関数を用いることで、 任意のバージョンのPFPをセットすることができます。

python
    
    
    # PFP version v0.0.0
    estimator_v0 = Estimator(model_version = "v0.0.0")
    print(f"Current model version:: {estimator_v0.model_version}")
    
    # Change the PFP version to v1.0.0
    estimator_v0.set_model_version("v1.0.0")
    print(f"Current model version:: {estimator_v0.model_version}")

異なるバージョンのPFPを `Estimator` にセットし、potential energyの予測結果を比較してみましょう。

python
    
    
    # PFP version v0.0.0
    estimator.set_model_version("v0.0.0")
    calculator.reset()
    print(f"Model version: {estimator.model_version}, Potential energy: {atoms.get_potential_energy()} eV")
    
    # Latest PFP version
    estimator.set_model_version("latest")
    calculator.reset()
    print(f"Model version: {estimator.model_version}, Potential energy: {atoms.get_potential_energy()} eV")

バージョンごとにPFPの学習に用いたデータが異なるため、予測値も異なります。 この差異は各種物性値の計算にも影響を及ぼす可能性があるので、 Matlantisを用いた結果を比較する際にはPFPバージョンの一貫性にご注意ください。

また、バージョンごとのPFPの対応元素の取得は以下のように行うことができます。

python
    
    
    from pfp_api_client.pfp.estimator import EstimatorElementStatus
    
    # Expected elements of the PFP version set to the current Estimator.
    elements = estimator.supported_elements()
    print(f"Model version: {estimator.model_version}, Supported elements: {elements}")

元素の対応状況には `Expected` と `Experimental` があります。幅広い系での再現が期待される元素を `Expected` 、一部の系で再現が期待される元素を `Experimental` としています。 デフォルトでは、 `supported_elements()` はestimatorにセットされているPFPバージョンの対応元素のうち `Expected` な元素のみを出力するようになっています。

以下のように、PFPのバージョンや対応状況を指定して取得するすることもできます。PFPの初期バージョンと最新バージョンの対応元素を比較してみましょう。

python
    
    
    # Expected + Experimental elements of the specified PFP version.
    model_version = "v0.0.0"
    elements = estimator.supported_elements(model_version = model_version, status=[EstimatorElementStatus.Expected, EstimatorElementStatus.Experimental])
    print(f"Model version: {model_version}, Supported elements: {elements}")
    
    # Expected + Experimental elements of the PFP latest version .
    model_version = "latest"
    elements = estimator.supported_elements(model_version = model_version, status=[EstimatorElementStatus.Expected, EstimatorElementStatus.Experimental])
    print(f"Model version: {model_version}, Supported elements: {elements}")

PFPの仕様についてはHelp内の [About PFP](/api/docs/ja/pfp-description/latest/index.html) を、 各バージョンの対応元素の詳細については [PFP Versions](/api/docs/ja/matlantis-guidebook/pfp_versions.html) の項目をご覧ください。

## PFVM技術について ​

高速かつ大きな系の推論を可能にする PFVM技術を正式版として提供しています。

pfp-api-client を用いた計算はすべて PFVM を用いて実行されます。

  * pfp-api-client v1.11.0 以前のバージョンでは、2023年8月23日以降は `EstimatorMethodType.TORCH` を指定した場合でも PFVM モードが有効になっています。
  * pfp-api-client v1.12.0 以降のバージョンでは、 `EstimatorMethodType.TORCH` は削除されました。 Torch モードを指定した場合はエラーが発生します。
  * pfp-api-client v1.19.0 以降のバージョンでは、 `EstimatorMethodType.PFVM_D3_PFVM` を指定できます。これは PFVM を用いた DFT-D3 補正を行うモードです。



モードを指定する際は、 `Estimator` の `method_type` オプションを使用します。pfp-api-client v1.21.3 以前のバージョンでは、デフォルト値として `EstimatorMethodType.PFVM` が設定されていますが、v1.22.0 以降のバージョンでは、デフォルト値として `EstimatorMethodType.PFVM_D3_PFVM` が設定されています。 以下は `EstimatorMethodType.PFVM_D3_PFVM` を指定する例です。

python
    
    
    from pfp_api_client import Estimator, EstimatorMethodType
    estimator = Estimator(method_type=EstimatorMethodType.PFVM_D3_PFVM)
    # estimator = Estimator(method_type="PFVM_D3_PFVM") # 文字列で指定することもできます。
    # estimator = Estimator(method_type="pfvm_d3_pfvm") # 文字列で指定する場合、大文字小文字は区別されません。

## 入力サイズ上限について ​

推論可能な系のサイズには上限があり、その値をAccount Dialogで確認できます。 PFPのバージョン毎にアーキテクチャが異なり、推論可能な最大近傍数が異なります。 また、D3補正を追加する計算モード(`PBE_PLUS_D3`, `PBE_U_PLUS_D3`, `R2SCAN_PLUS_D3`)を使用する場合には `method_type` によっても推論可能な最大近傍数が異なります。 Dialog画面の表では `v7.0.0` 以降と `v6.0.0` と `v5.0.0` 以前のモデルで推論可能な入力サイズを順に表示しています。 また、表の中で `method_type=EstimatorMethodType.PFVM` を指定した時のD3補正の推論可能な入力サイズはD3 correction (Torch)行に、 `method_type=EstimatorMethodType.PFVM_D3_PFVM` (default) を指定した時のD3補正の推論可能な入力サイズはD3 correction (PFVM)行に記載されています。

PFPおよびD3補正それぞれについて、入力サイズが上限値以下である必要があります。 原子(atom)数の制限が超えた場合には、 `Too many atoms` エラーが、周期境界条件を用いた系において近接セルの仮想原子(ゴースト原子)を含む原子数の制限を超えた場合には `Too many atoms considering periodic boundary conditions` エラーが、近傍(neighbor)数が制限を超えた場合には `Too many neighbors` エラーが発生します。 なお、これらの制約を超えてない入力の推論を保証するものではなく、超えていなくても `Too big for GPU memory` が発生する場合があります。 現在のPFPおよびDFT-D3補正では原子数より近傍数の方がメモリ使用量への影響が大きく、原子数が少ない場合でも近傍数が多いことで `Too big for GPU memory` エラーが発生することがあります。

仮想原子を含む原子数、PFPの入力近傍数、およびDFT-D3補正の入力近傍数は、以下のように `Estimator` の出力 `calc_stats` の中の `n_atoms_plus_ghost_atoms`、 `n_neighbors` (PFP) 、 `d3_n_neighbors` (DFT-D3) の値で確認することができます。

python
    
    
    from ase.build import bulk
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
    
    est = Estimator(calc_mode=EstimatorCalcMode.PBE_PLUS_D3)
    calculator = ASECalculator(est)
    
    si = bulk('Si')
    si.calc = calculator
    si.get_potential_energy()
    
    print(si.calc.results['calc_stats'])
    # {'n_atoms_plus_ghost_atoms': 559, 'n_neighbors': 158, 'd3_n_neighbors': 596}

## 学習データのエネルギー補正値について ​

PFPと第一原理計算ソフトウェアとでは、エネルギーが0となる基準が異なります。 これらの間でエネルギーの値を直接比較するためには、元素毎に定められたエネルギー補正値(shift energy)を用いて値を補正する必要があります。 次の例では、元素ごとのエネルギー補正値をdictとして取得し、水分子での合計補正値を求めています。

python
    
    
    # Dict of atomic number vs. shift energy of a specific model version and calc mode
    shift_energies = estimator.get_shift_energy_table(model_version="v8.0.0", calc_mode="WB97XD")
    
    # If you want to get the shift
    
    print(f"Shift energy of H: {shift_energies[1]} eV.")
    
    # Total shift energy of a molecule
    from ase.build import molecule
    atoms = molecule("H2O")
    total_shift_energy = sum([shift_energies[i] for i in atoms.get_atomic_numbers()])
    print(f"Shift energy of H2O: {total_shift_energy} eV.")

この値はモデルバージョンおよび計算モード (PBE, PBE_U, R2SCAN, WB97XD) でそれぞれ異なります。 PBE, PBE_UおよびR2SCANでの値はVASPと比較する際に、WB97XDでの値はGaussianと比較する際に用いることができます。 shift energyの値はD3補正による影響を受けないため、PBE_PLUS_D3 (PBE_U_PLUS_D3) はPBE (PBE_U) と同じ値を返します。

一部の元素に関してはshift energyの値が0となっており、PFPと第一原理計算ソフトウェアとでエネルギーが0となる基準が同じとなっています。

## トークン制および優先度機能について ​

トークン制では、テナントごとに「トークン」が付与されます。トークンは、 `Estimator` のインスタンスを通じて計算することで消費され、時間経過に応じて回復します。トークンの時間あたりの回復量(以下トークン回復速度)は、ご利用のプランに応じて設定されます。

トークンのご利用状況をご覧になるには、Notebookの右上のPFP Load Statusをクリックしてください。

「Token Consumption Rate」では、ユーザーごとのトークンの消費速度の、トークン回復速度に対する割合が積み上げグラフで表示されています。

  * Token Consumption Rate の合計が100%を超えている場合は、バースト中です。トークン回復速度より多くのトークンを消費されている状態であり、Token Availability は0%になるまで減少していきます。
  * Token Consumption Rate の合計が100%を下回っている場合は、トークン回復速度よりトークン消費速度が小さい状態であり、Token Availability は 100% 近くになるまで増加していきます。



「Token Availability」のグラフでは、テナント全体で現在貯まっているトークンの、トークン容量に対する割合を確認することができます。

  * Token Availability が 0% の場合、トークンを使い切っており、バーストできない状態です。つまり、ちょうどトークン回復速度でしか計算ができない状態です。
  * Token Availability が 100% の場合、トークンが余っている状態です。この場合、トークンの消費速度よりトークン回復速度の方が多かったとしても、トークンはそれ以上回復しません。



トークン容量は、トークン回復速度に比例して設定されます。

### トークンの消費速度を制限する ​

特定のユーザーがテナント全体のトークンを大量に消費することを防ぐために、管理者はユーザーまたはチームごとにトークン消費速度の上限を設定できます。

以下のように、ユーザーの編集画面から「Max Consumption Percentage」を設定します。

Max Consumption Percentage がユーザーに設定されている場合、そのユーザーの Token Consumption Rate が設定された割合を超えないように、計算が適宜保留されます。

同様に、チームの編集画面からも「Max Consumption Percentage」を設定できます。

チームに設定されている場合、そのチームに所属するユーザーの Token Consumption Rate の合計が設定された割合を超えないように、計算が保留されます。

Max Consumption Percentage が設定された複数のチームに所属している場合、最も高い割合が設定されたチームに所属しているものとして扱います。 例えば、チームAに50%、チームBに80%が設定されている場合、チームBのメンバーとしてトークンが消費されます。 どのチームとしてトークンが消費されるかはPFP Load StatusのYour Team Consumption Rateで確認できます。

ユーザーとチームの両方に Max Consumption Percentage が設定されている場合、それぞれの制限が独立して適用されます。以下に例を示します。

  * ユーザーに30%、チームに70%が設定されている場合、チームの Token Consumption Rate に余裕があったとしても、30%を超える計算は実行できません。
  * ユーザーに80%、チームに50%が設定されている場合、先にチームの制限に達するため、50%を超える計算は実行できません。



### 優先度機能 ​

優先度は、1以上100以下の自然数です。これによってToken Availabilityが少ない状況においてどの計算が優先されるかを設定できます。優先度はデフォルトで100です。

  * 優先度として100を設定した場合、Token Availability と関係なく計算を試みます。
  * 優先度として75を設定した場合、Token Availability の残量が 25% 以上になるまで 計算は保留されます。
  * 優先度として25を設定した場合、Token Availability の残量が 75% 以上になるまで 計算は保留されます。
  * 優先度として1を設定した場合、Token Availability の残量が 99% 以上になるまで 計算は保留されます。



計算が保留された場合、Estimator は成功するまで何回でもリトライを行います。タイムアウトは設定されていません。

なおバースト中でも、トークンの消費速度には一定の上限が存在しており、これに到達した場合、計算は保留されます。

優先度の設定方法は複数あります。複数箇所で設定した場合、コンストラクタ > Config > 環境変数の強さで設定が優先されます。

`Estimator` クラスのコンストラクタで優先度の設定が可能です。

python
    
    
    from pfp_api_client.pfp.estimator import Estimator
    estimator = Estimator(priority=80)
    
    assert estimator.priority == 80

以下のように `Config.priority` という変数を設定することでも優先度の設定が可能です。

python
    
    
    from pfp_api_client.pfp.estimator import Estimator, Config
    Config.priority = 80
    estimator = Estimator()
    
    assert estimator.priority == 80

以下のように `PFP_DEFAULT_PRIORITY` という環境変数を設定することでも優先度の設定が可能です。

python
    
    
    from pfp_api_client.pfp.estimator import Estimator, Config
    estimator = Estimator()
    
    print(estimator.priority)

bash
    
    
    $ PFP_DEFAULT_PRIORITY=80 python test_environment_priority.py
    80

何も設定しない場合、優先度はデフォルトで100です。

python
    
    
    from pfp_api_client.pfp.estimator import Estimator
    estimator = Estimator()
    
    assert estimator.priority == 100

## pfp-api-client のロギングについて ​

`pfp-api-client` のログ出力設定は `pfp_api_client.logging` モジュールで変更できます。

`pfp-api-client` のログ出力レベルはデフォルトでは `INFO` に設定されています。 これを変更するには `pfp_api_client.logging.set_verbosity` を使用します。 例えばログ出力レベルを `DEBUG` に変更すると `Estimator` を作成したときの設定が出力されるようになります。

python
    
    
    from pfp_api_client.pfp.estimator import Estimator
    from pfp_api_client.logging import set_verbosity
    Estimator() # => nothing is shown.
    set_verbosity("DEBUG")
    Estimator()
    # => [DEBUG 2025-08-21 12:26:32,507] Created the estimator(model_version=v8.0.0, calc_mode=EstimatorCalcMode.PBE, method_type=EstimatorMethodType.PFVM_D3_PFVM)

### ログのプロパゲーション ​

ユーザーがルートロガーを設定した場合、デフォルトではルートロガーと `pfp-api-client` のロガーの双方でログが出力される場合があります。

python
    
    
    import logging
    logging.basicConfig(
        format="[user-log] %(message)s",
        level="INFO",
    )
    print("Enabled propagation by default:")
    Estimator()
    # => [DEBUG 2025-08-21 12:32:43,196] Created the estimator(model_version=v8.0.0, calc_mode=EstimatorCalcMode.PBE, method_type=EstimatorMethodType.PFVM_D3_PFVM)
    # => [user-log] Created the estimator(model_version=v8.0.0, calc_mode=EstimatorCalcMode.PBE, method_type=EstimatorMethodType.PFVM_D3_PFVM)

ルートロガーから `pfp-api-client` のログを出力させないようにするには `pfp_api_client.logging.disable_propagation` を呼び出します。

python
    
    
    print("Disable propagation:")
    from pfp_api_client.logging import disable_propagation
    disable_propagation()
    estimator.Estimator()
    # => [DEBUG 2025-08-21 12:32:43,299] Created the estimator(model_version=v8.0.0, calc_mode=EstimatorCalcMode.PBE, method_type=EstimatorMethodType.PFVM_D3_PFVM)
