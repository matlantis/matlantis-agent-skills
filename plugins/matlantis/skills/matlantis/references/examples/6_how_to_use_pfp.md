# How to use PFP

ここではPFPの性能を最大限に活用するために、PFPの挙動を制御する設定値の説明や内部動作に関する解説を行います。&lt;br/&gt;
本Tutorial を実行するためには、`pfp_api_client&gt;=1.23.0` のバージョンが必要です。


## Setup

`pfp_api_client`のアップデート。&lt;br/&gt;
一度実行すれば、再度このノートブックを開いた場合や、他のノートブックを利用する場合に再実行する必要はありません。


```python
!pip install 'pfp-api-client&gt;=1.23.0'
```




### 必要な moduleの import

本Tutorialは以下のVersionで動作確認を行っています。


```python
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode, EstimatorElementStatus

print(f"pfp_api_client: {pfp_api_client.__version__}")
```

    pfp_api_client: 2.0.1



```python
from ase.build import molecule

atoms = molecule("H2")
```


## Estimator

`Estimator`は実際に系の構造情報を入力としてNNP (Neural Network Potential)の計算結果であるエネルギーや力を出力とする部分を担当します。&lt;br/&gt;
(実際にはこのNotebook上で計算を行っているわけではなく、内部でGPU Serverに対して構造情報を投げ、計算結果を受け取っています。)


### モデルのバージョンを指定する


`Estimator`をインスタンス化する際の **`model_version`** という引数を用いることで、計算に使用するモデルのバージョン指定ができます。

 - v0.0.0: PFP論文執筆時のモデル
 - v1.0.0: Matlantis製品リリース時のはじめのモデル
 - v1.1.0: v1.0.0に、D3補正を加えた "PBE_U_PLUS_D3"モードをサポートしたモデル
 - v2.0.0: 学習データセットを拡張し、Hubbard 補正のない結晶系モードPBEとPBE_PLUS_D3をサポートしたモデル
 - v3.0.0: 学習データセットを拡張し、対応元素の種類を新たに17元素を追加し、周期表の大部分をカバーする72元素に対応したモデル
 - v4.0.0: 学習データセットを拡張し、Hubbard 補正のない結晶系モードPBEモードをデフォルトとしたモデル
 - v5.0.0: v4.0.0より学習データセットを拡張したモデル
 - v6.0.0: 真空中の分子など疎な環境や高密度な構造での再現性を改善したモデル
 - v7.0.0: 学習データセットを拡張し、対応元素の種類を新たに24元素を追加し、自然界に存在する全ての元素をカバーする96元素に対応したモデル
 - v8.0.0: r2SCANとr2SCAN+D3補正のcalc_modeを追加したモデル


```python
estimator = Estimator(model_version="v8.0.0")
```


`model_version` として指定可能なバージョンは `estimator.available_models` を参照することで確認ができます。


```python
estimator.available_models
```




    ['latest', 'v0.0.0', 'v1.0.0', 'v1.1.0', 'v2.0.0', 'v3.0.0', 'v4.0.0', 'v5.0.0', 'v6.0.0', 'v7.0.0', 'v8.0.0']




デフォルトでは `model_version="latest"` が指定されていて、使用している`pfp-api-client`パッケージリリーズ時の最新モデルが使われる設定となっています。 (`pfp-api-client==1.23.0`では、`"v8.0.0"`が使われます。)&lt;br/&gt;
`available_models` はPFP API Server側で利用可能なモデルを示しているため、`pfp-api-client`のライブラリバージョンが同じでも将来追加される可能性があります。


また、 `estimator.set_model_version` 関数を用いてあとから使用するモデルのVersionを変更することも可能です。&lt;br/&gt;
以下の例では `v8.0.0` と `v7.0.0` で計算を行ってみています。得られる値が少しだけ違うことが確認できます。

※ `calculator.reset()` の動作については後述します。


```python
calculator = ASECalculator(estimator)
atoms.calc = calculator

estimator.set_model_version("v8.0.0")
calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {estimator.model_version}, Potential energy: {Epot} eV")

estimator.set_model_version("v7.0.0")
calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {estimator.model_version}, Potential energy: {Epot} eV")
```

    Model version: v8.0.0, Potential energy: -4.523196788877179 eV
    Model version: v7.0.0, Potential energy: -4.514524186281303 eV



### 計算モードを指定する

PFPでは、PBE汎関数・PAW型擬ポテンシャルを用いたPBE_U Dataset（Hubbard U補正あり）・PBE Dataset（Hubbard U補正なし）、r2SCAN Dataset（Hubbard U補正なし）および、ωB97xd汎関数・6-31G(d)基底を用いた wb97xd Datasetを学習したモデルにより、それぞれのモードでの推論が可能です。&lt;br/&gt;
wb97xd modeは真空中の有機分子のみで使用してください。

計算モードが指定されない場合は、デフォルトの計算モード（v3.0.0 以前では `PBE_U` 、v4.0.0 以降では `PBE` ）が適用されます。詳細に関しては [About PFP](/api/resource/documents/pfp-description/versions/latest/index.html#_22) をご参照ください。

計算可能なモードは `EstimatorCalcMode` のEnumを参照することで確認ができます。


```python
list(EstimatorCalcMode)
```




    [&lt;EstimatorCalcMode.CRYSTAL: 'crystal'&gt;,
     &lt;EstimatorCalcMode.CRYSTAL_PLUS_D3: 'crystal_plus_d3'&gt;,
     &lt;EstimatorCalcMode.CRYSTAL_U0: 'crystal_u0'&gt;,
     &lt;EstimatorCalcMode.CRYSTAL_U0_PLUS_D3: 'crystal_u0_plus_d3'&gt;,
     &lt;EstimatorCalcMode.MOLECULE: 'molecule'&gt;,
     &lt;EstimatorCalcMode.R2SCAN: 'r2scan'&gt;,
     &lt;EstimatorCalcMode.R2SCAN_PLUS_D3: 'r2scan_plus_d3'&gt;,
     &lt;EstimatorCalcMode.PBE: 'pbe'&gt;,
     &lt;EstimatorCalcMode.PBE_PLUS_D3: 'pbe_plus_d3'&gt;,
     &lt;EstimatorCalcMode.PBE_U: 'pbe_u'&gt;,
     &lt;EstimatorCalcMode.PBE_U_PLUS_D3: 'pbe_u_plus_d3'&gt;,
     &lt;EstimatorCalcMode.WB97XD: 'wb97xd'&gt;]




計算モードは `calc_mode` 引数で指定を行う事ができます。


```python
pbe_estimator = Estimator(calc_mode=EstimatorCalcMode.PBE)
pbe_calculator = ASECalculator(pbe_estimator)
```


```python
pbe_u_estimator = Estimator(calc_mode=EstimatorCalcMode.PBE_U)
pbe_u_calculator = ASECalculator(pbe_u_estimator)
```


```python
r2scan_estimator = Estimator(calc_mode=EstimatorCalcMode.R2SCAN)
r2scan_calculator = ASECalculator(r2scan_estimator)
```


```python
wb97xd_estimator = Estimator(calc_mode=EstimatorCalcMode.WB97XD)
wb97xd_calculator = ASECalculator(wb97xd_estimator)
```


```python
wb97xd_estimator.calc_mode
```




    &lt;EstimatorCalcMode.WB97XD: 'wb97xd'&gt;




```python
atoms.calc = pbe_calculator
pbe_calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {pbe_estimator.model_version}, {pbe_estimator.calc_mode.name}  mode, Potential energy: {Epot} eV")

atoms.calc = pbe_u_calculator
pbe_u_calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {pbe_u_estimator.model_version}, {pbe_u_estimator.calc_mode.name}  mode, Potential energy: {Epot} eV")

atoms.calc = r2scan_calculator
r2scan_calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {r2scan_estimator.model_version}, {r2scan_estimator.calc_mode.name}  mode, Potential energy: {Epot} eV")

atoms.calc = wb97xd_calculator
wb97xd_calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {wb97xd_estimator.model_version}, {wb97xd_estimator.calc_mode.name} mode, Potential energy: {Epot} eV")
```

    Model version: v8.0.0, PBE  mode, Potential energy: -4.523197861760831 eV
    Model version: v8.0.0, PBE_U  mode, Potential energy: -4.523196550458586 eV
    Model version: v8.0.0, R2SCAN  mode, Potential energy: -4.527281614584459 eV
    Model version: v8.0.0, WB97XD mode, Potential energy: -4.525407435848635 eV



また、 `estimator.calc_mode` を更新することで後からモデルの計算モードを変更することも可能です。&lt;br/&gt;
以下の例では `wb97xd` mode 、`PBE_U` mode 、`PBE` mode、`R2SCAN` mode で計算を行ってみています。得られる値が少しだけ違うことが確認できます。


```python
estimator = Estimator(model_version="latest")

calculator = ASECalculator(estimator)
atoms.calc = calculator

estimator.calc_mode = EstimatorCalcMode.PBE
calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {estimator.model_version}, {estimator.calc_mode.name}  mode, Potential energy: {Epot} eV")

estimator.calc_mode = EstimatorCalcMode.PBE_U
calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {estimator.model_version}, {estimator.calc_mode.name}  mode, Potential energy: {Epot} eV")

estimator.calc_mode = EstimatorCalcMode.R2SCAN
calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {estimator.model_version}, {estimator.calc_mode.name} mode, Potential energy: {Epot} eV")

estimator.calc_mode = EstimatorCalcMode.WB97XD
calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {estimator.model_version}, {estimator.calc_mode.name} mode, Potential energy: {Epot} eV")
```

    Model version: v8.0.0, PBE  mode, Potential energy: -4.523197444528298 eV
    Model version: v8.0.0, PBE_U  mode, Potential energy: -4.523196848481825 eV
    Model version: v8.0.0, R2SCAN mode, Potential energy: -4.527281614584459 eV
    Model version: v8.0.0, WB97XD mode, Potential energy: -4.525405737116179 eV



なお、wb97xd Datasetは周期境界条件無しで、PBE/PBE_U/r2SCAN Datasetは周期境界条件有りでの教師データが収集されていますが、推論時には周期境界条件有り・無しどちらでも推論することは可能です。


`pfp_api_client&gt;=1.2.2`からは、3方向それぞれの周期境界条件が異なるようなケース(たとえばxy平面に続くSlab `pbc=[True, True, False]`)でも推論ができます。


```python
from ase.build import fcc111

atoms = fcc111("Pt", (3, 3, 3), vacuum=10.0)
atoms.calc = ASECalculator(Estimator())
Epot = atoms.get_potential_energy()
print(f"atoms.pbc = {atoms.pbc}, Epot={Epot:.2f} eV")
```

    atoms.pbc = [ True  True False], Epot=-136.87 eV



v2.0.0以降のモデルはHubbard補正あり(`PBE_U`)となし(`PBE`)の計算モードの両方に対応しています。遷移元素（V、Cr、Mn、Fe、Co、Ni、Cu、Mo、W）を含む場合、それぞれのモードで結果が異なります。以下の例では、`PBE_U` mode と `PBE` mode で計算を行い、得られる結果が異なることが確認しています。


```python
from ase.build import bulk
atoms = bulk("Fe") * (2,2,2)

atoms.calc = pbe_u_calculator
pbe_u_calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {pbe_u_estimator.model_version}, {pbe_u_estimator.calc_mode.name}    mode, Potential energy: {Epot} eV")

atoms.calc = pbe_calculator
pbe_calculator.reset()
Epot = atoms.get_potential_energy()
print(f"Model version: {pbe_estimator.model_version}, {pbe_estimator.calc_mode.name} mode, Potential energy: {Epot} eV")
```

    Model version: v8.0.0, PBE_U    mode, Potential energy: -18.241910249590482 eV
    Model version: v8.0.0, PBE mode, Potential energy: -39.91537051649998 eV



### 対応元素種について

PFPの対応元素種に関して、詳細は "Help" --&gt; "About PFP" をご参照ください。

推論可能な元素種は以下の３つのクラスに分かれています。
 - Expected element: 現時点でPFPが再現する対象として主に注力している元素。
 - Experimental element: 主なターゲットとしていないが、部分的に再現できる可能性があるもの。
 - Unexpected element: 再現できないもの。

Experimental, Unexpectedな元素種を含めて推論を実行した場合もエラーにはならずに推論値が返ってきますが、その精度は保証されないため注意をしてください。

この場合はWARNING Messageが表示されます。実際に例を見てみます。


```python
# Since there are no experimental elements in version 7.0.0 or later,
# we will use version 6.0.0 in this section, in order to display the warning
pfp_estimator = Estimator(model_version="v6.0.0", calc_mode=EstimatorCalcMode.PBE)
pfp_calculator = ASECalculator(pfp_estimator)
```


```python
from ase.build import bulk

# When using an expected element: Warning will not be displayed
atoms = bulk("Pt") * (2, 2, 2)
atoms.calc = pfp_calculator
atoms.get_potential_energy()
```




    -43.54549316209632




```python
# Experimental element: Warning will be displayed
atoms = bulk("Tc") * (2, 2, 2)
atoms.calc = pfp_calculator
atoms.get_potential_energy()
```

    [WARNING 2025-11-04 09:10:45,786] Runtime warning (ID 1001): Experimental element was detected.         You can suppress this message with Estimator.set_message_status().





    -109.91673105729285




```python
from ase.build import bulk

# Unexpected element: Warning will be displayed
atoms = bulk("Eu") * (2, 2, 2)
atoms.calc = pfp_calculator
atoms.get_potential_energy()
```

    [WARNING 2025-11-04 09:10:46,023] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.





    -14.86441255316399




もし複数のクラスの元素種が含まれていた場合は、一番WARNING度の高いメッセージのみが表示されます。


```python
from ase import Atoms

# Unexpected element: Warning will be displayed
atoms = Atoms(["H", "Tc", "Eu"], positions=[[0, 0, 0], [0, 0, 3], [0, 0, 6]])
atoms.calc = pfp_calculator
atoms.get_potential_energy()
```

    [WARNING 2025-11-04 09:10:46,223] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.





    -1.9796038205647146




wb97xd 計算モードを使用した場合は、wb97xd 計算モードでサポート外の元素にはWARNINGが出力されます（注: pfp-api-client 1.13.4以降）。


```python
from ase.build import bulk

# A warning will be displayed since Pt is not included in the wb97xd dataset.
atoms = bulk("Pt") * (2, 2, 2)
pfp_estimator.set_calc_mode(EstimatorCalcMode.WB97XD)
atoms.calc = pfp_calculator
atoms.get_potential_energy()
```

    [WARNING 2025-11-04 09:10:46,483] Runtime warning (ID 1002): Unexpected element was detected. There might be strange behavior.





    -43.673545426413106




これら元素の一覧は、`supported_elements` methodを用いることでの取得も可能です。


```python
# Since there are no experimental elements in version 7.0.0 or later,
# we will use version 6.0.0 in this section, in order to display the warning
pfp_estimator = Estimator(model_version="v6.0.0")

expected_elements = pfp_estimator.supported_elements(status=EstimatorElementStatus.Expected, calc_mode=EstimatorCalcMode.PBE)
experimental_elements = pfp_estimator.supported_elements(status=EstimatorElementStatus.Experimental, calc_mode=EstimatorCalcMode.PBE)

print("Expected    : ", expected_elements)
print("Experimental: ", experimental_elements)
```

    Expected    :  ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Gd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Pb', 'Bi']
    Experimental:  ['Tc', 'Tl']



```python
expected_elements = pfp_estimator.supported_elements(status=EstimatorElementStatus.Expected, calc_mode=EstimatorCalcMode.WB97XD)
experimental_elements = pfp_estimator.supported_elements(status=EstimatorElementStatus.Experimental, calc_mode=EstimatorCalcMode.WB97XD)

print("Expected    : ", expected_elements)
print("Experimental: ", experimental_elements)
```

    Expected    :  ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br']
    Experimental:  []



```python
# Since r2SCAN mode is supported after v8.0.0, we will use version 8.0.0 in this section, in order to display expected elements
pfp_estimator = Estimator(model_version="v8.0.0")

expected_elements = pfp_estimator.supported_elements(status=EstimatorElementStatus.Expected, calc_mode=EstimatorCalcMode.R2SCAN)
experimental_elements = pfp_estimator.supported_elements(status=EstimatorElementStatus.Experimental, calc_mode=EstimatorCalcMode.R2SCAN)

print("Expected    : ", expected_elements)
print("Experimental: ", experimental_elements)
```

    Expected    :  ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']
    Experimental:  []



#### WARNING Message の抑制

上記WARNINGを理解した上でExperimental 元素を使った推論を試したい場合などは、以下の関数でWARNING messageの抑制ができます。


```python
from pfp_api_client.pfp.utils.messages import MessageEnum

# Since there are no experimental elements in version 7.0.0 or later,
# we will use version 6.0.0 in this section, in order to display the warning
estimator = Estimator(model_version="v6.0.0")
estimator.set_message_status(MessageEnum.ExperimentalElementWarning, False)

# Not recommended as it makes it difficult to notice that it contains Unexpected elements carelessly.
# estimator.set_message_status(MessageEnum.UnexpectedElementWarning, False)
```


```python
from ase.build import bulk

# You can see that WARNING is not displayed even when using Unexpected element
atoms = bulk("Tc") * (2, 2, 2)
atoms.calc = ASECalculator(estimator)
atoms.get_potential_energy()
```




    -109.91673291464431




## [補足] ASEのCalculatorに関する説明

本章はMatlantisとは独立した、ASE の一般的な説明となっています。


### AtomsとCalculatorの関係

ASEライブラリでは、原子構造をAtomsクラスで表し、そこにCalculatorをセットすることでそのAtomsに関する各基本物性値（エネルギー・力・ストレス・電荷など）を計算することができます。&lt;br/&gt;
Calculatorをセットする際には `atoms.calc` に対して直接 calculatorをセットします。

**Atoms と Calculatorの関係**

&lt;img src="./assets/tutorial_pfp_usage/atoms-calculator.png" width="350"/&gt;

Calculatorを経由して計算できる基本物性値と、その計算methodは以下のようになっています。

 - ポテンシャルエネルギー: `get_potential_energy`
 - 力: `get_forces`
 - 応力(ストレス): `get_stress`
 - 電荷: `get_charges`
 - マグネティックモーメント: `get_magnetic_moment`
 - ダイポールモーメント: `get_dipole_moment`



```python
calculator = ASECalculator(Estimator())

atoms = bulk("Pt") * (2, 2, 1)
atoms.calc = calculator

E_pot = atoms.get_potential_energy()
charges = atoms.get_charges()
forces = atoms.get_forces()
stress = atoms.get_stress()

print(f"E_pot {E_pot:.2} eV")
print(f"charges {charges} C")
print(f"forces {forces} eV/A")
print(f"stress {stress} eV/A^2")
```

    E_pot -2.2e+01 eV
    charges [-4.87786878e-08  9.12191354e-08  1.55771360e-07 -1.98211922e-07] C
    forces [[ 5.30744388e-07  2.22680500e-07  1.32366059e-08]
     [ 1.08821656e-08  4.76849597e-07 -9.44746923e-07]
     [ 1.50329877e-07 -5.50746567e-07  1.35660978e-06]
     [-6.91956431e-07 -1.48783530e-07 -4.25099462e-07]] eV/A
    stress [-6.11356949e-02 -6.11356781e-02 -6.11356492e-02 -6.10226058e-09
     -3.75040920e-08  5.64582622e-09] eV/A^2



PFP では、すべての方向に周期境界条件が適用されていない場合、応力の計算には対応していません。


```python
# It is possible to be executed depending on the Calculator, but an error occurs on the PFP Calculator.

# atoms.set_pbc([False, True, True]) # unset PBC for x-direction
# stress = atoms.get_stress()
```


PFP では、 magnetic momentや dipole moment の計算には対応していません。


```python
# It is possible to be executed depending on the Calculator, but an error occurs on the PFP Calculator.

# magmom = atoms.get_magnetic_moment()
# dipole = atoms.get_dipole_moment()
```


### Calcultor の計算キャッシュ機能

Calculatorは前に計算を行った際の入力原子構造を `calculator.atoms`に、　またそのときの計算結果を `calculator.results` に保持しています。&lt;br/&gt;
そのため、以前と全く同じ原子構造に対して物性値の `get_XXX` 関数を呼んだ場合、 **再計算をスキップ** するような仕組みになっています。


```python
atoms.calc = calculator
calculator.reset()
```


```python
%time Epot = atoms.get_potential_energy()
%time Epot = atoms.get_potential_energy()
```

    CPU times: user 2.52 ms, sys: 141 μs, total: 2.67 ms
    Wall time: 83.7 ms
    CPU times: user 284 μs, sys: 16 μs, total: 300 μs
    Wall time: 304 μs



上記の例でWall timeを比較すると、１回目の実行ではミリ秒オーダーの時間がかかっており計算が実行されていますが、２回めの実行では計算がスキップされているためミリ秒もかからずにエネルギーを得ることができていることがわかります。


`calculator.reset()` を呼ぶことで明示的に計算結果のキャッシュをクリアして再計算を行うこともできます。以下の例では2回めの計算でもミリ秒の計算時間がかかっていることがわかります。


```python
calculator.reset()

print("----- 1st calc -----")
%time Epot = atoms.get_potential_energy()

print(f"After 1st calc   : {calculator.results}")
calculator.reset()
print(f"After reset      : {calculator.results}")

print("----- 2nd calc -----")
%time Epot = atoms.get_potential_energy()
print(f"After 2nd calc   : {calculator.results}")
```

    ----- 1st calc -----
    CPU times: user 3.19 ms, sys: 0 ns, total: 3.19 ms
    Wall time: 63.1 ms
    After 1st calc   : {'energy': -21.833274701522566, 'forces': array([[-6.04930375e-07,  4.06389801e-07,  1.03574484e-06],
           [ 4.04037239e-07,  2.26034487e-07, -7.55872564e-07],
           [ 1.07867083e-07, -3.65861281e-07, -7.20368305e-07],
           [ 9.30260536e-08, -2.66563007e-07,  4.40496030e-07]]), 'charges': array([ 1.82970552e-07, -1.12249623e-07, -1.45855381e-07,  7.51345439e-08]), 'calc_stats': {'elapsed_usec_preprocess': 2857, 'elapsed_usec': 54670, 'elapsed_usec_infer': 11341, 'n_atoms_plus_ghost_atoms': 648, 'n_neighbors_after_screening': 3541, 'n_neighbors': 400}, 'free_energy': -21.833274701522566, 'stress': array([-6.11354343e-02, -6.11354598e-02, -6.11354023e-02, -2.21590522e-08,
            1.02186774e-09,  9.58753409e-09])}
    After reset      : {}
    ----- 2nd calc -----
    CPU times: user 2.65 ms, sys: 0 ns, total: 2.65 ms
    Wall time: 62.3 ms
    After 2nd calc   : {'energy': -21.833272039686683, 'forces': array([[-4.56774814e-07,  5.43982950e-07, -1.60109391e-07],
           [ 2.13353988e-08, -2.22978808e-07, -2.25485989e-07],
           [-2.14223475e-07,  1.61290689e-07, -5.13372955e-07],
           [ 6.49662890e-07, -4.82294831e-07,  8.98968334e-07]]), 'charges': array([ 2.04857926e-07, -1.45283892e-07, -2.32902153e-07,  1.73328090e-07]), 'calc_stats': {'elapsed_usec_preprocess': 2251, 'elapsed_usec': 51667, 'elapsed_usec_infer': 7774, 'n_atoms_plus_ghost_atoms': 648, 'n_neighbors_after_screening': 2588, 'n_neighbors': 400}, 'free_energy': -21.833272039686683, 'stress': array([-6.11355911e-02, -6.11355642e-02, -6.11355802e-02,  4.55160113e-09,
            1.76677919e-09, -1.91847793e-09])}



なお、原子構造が少しでも変わった場合には `calculator.result`の結果を破棄し、新しく計算がされます。

以下の例では、水素分子の原子間距離を変えて再度 `atoms.get_potential_energy`を呼んでいますが、この場合は入力原子構造が変わったことを `Calculator` が検知して新しく計算を行います。実際にどちらの計算にもミリ秒オーダーのWall timeがかかっていることがわかります。


```python
atoms = Atoms(["H", "H"], positions=[[0, 0, 0], [0, 0, 0.8]])

atoms.calc = calculator
%time E_pot1 = atoms.get_potential_energy()
# --- Here, the calculated atoms are stored in calculator.atoms.
# print(calculator.atoms.positions)

# Changed interatomic distance to 2A
atoms.positions[1, 2] = 2.0
# --- Since the structure of calculator.atoms calculated last time and atoms calculated this time are different, the calculation is executed.
%time E_pot2 = atoms.get_potential_energy()
# print(calculator.atoms.positions)

print(f"E_pot1 {E_pot1:.2f} eV")
print(f"E_pot2 {E_pot2:.2f} eV")
```

    CPU times: user 0 ns, sys: 3.16 ms, total: 3.16 ms
    Wall time: 62.6 ms
    CPU times: user 2.38 ms, sys: 138 μs, total: 2.51 ms
    Wall time: 72.3 ms
    E_pot1 -4.49 eV
    E_pot2 -0.13 eV



## ASECalculator

`pfp-api-client`ライブラリが提供している`ASECalculator`の動作に関する説明を行います。


PFPでは `potential_energy`, `charge` がNNPのForwardで、`forces`, `stress` がNNPのBackwardで計算されるという挙動になっています。

&lt;img src="./assets/tutorial_pfp_usage/nnp-calc.png"/&gt;


`pfp-api-client&gt;=1.2.2` では、常に`energy, charges, forces`は計算されます(※)。このように常に計算されるpropertyはcalculatorの`default_properties` で確認ができます。

以下の例では`get_potential_energy()` を呼んだ際にすでに`forces` が計算されており、`calculator.results`に格納されているため、次の`get_forces()` では計算時間がかかっていないことがWall timeの時間を見ることで確認できます。

※　`pfp-api-client&lt;1.2.2` では`potential_energy` の計算時は`forces`は計算されない仕様でした。


```python
atoms = molecule("H2")
atoms.calc = calculator
calculator.reset()

print(f"default_properties: {calculator.default_properties}")

print("--- calculate potential energy ---")
%time E_pot1 = atoms.get_potential_energy()
print(f"calculator.results = {calculator.results}")

print("--- calculate forces  ------------")
%time forces = atoms.get_forces()
```

    default_properties: ['energy', 'charges', 'forces', 'free_energy']
    --- calculate potential energy ---
    CPU times: user 0 ns, sys: 3.66 ms, total: 3.66 ms
    Wall time: 108 ms
    calculator.results = {'energy': -4.523196550458586, 'forces': array([[-0.        , -0.        ,  0.46052358],
           [-0.        , -0.        , -0.46052358]]), 'charges': array([ 4.26618492e-08, -4.26618492e-08]), 'calc_stats': {'elapsed_usec_preprocess': 3391, 'elapsed_usec': 79620, 'elapsed_usec_infer': 16539, 'n_atoms_plus_ghost_atoms': 2, 'n_neighbors_after_screening': 3788, 'n_neighbors': 1}, 'free_energy': -4.523196550458586}
    --- calculate forces  ------------
    CPU times: user 286 μs, sys: 17 μs, total: 303 μs
    Wall time: 308 μs




### set_default_properties

系のエネルギーだけを計算したいような場合、neural networkのBackward が必要な`forces`の計算を省くことで計算時間を節約することができます。(`forces`は不要で`energy`のみが必要になる用途としては、たとえばモンテカルロ法の各stepでエネルギーのみを評価して構造変化させていく場合などがあります。)

このような場合は、`calculator.set_default_properties` 関数を用いて、常に計算されるpropertyから`forces`を明示的に省くことで対応が可能です。&lt;br/&gt;
`get_potential_energy()` を呼んだ場合の実行時間よりも、 `get_forces` を呼んだ場合の実行時間のほうが長いになることがわかります(おおよそですが2倍程度になることが多いです)。

※ Matlantisでは、複数のGPUを組み合わせたシステムで推論を行っており、実際にどのGPUで推論が行われるかにより推論時間が変わるため、 `get_forces` のほうが早くなることも有りえます。


```python
# We explicitly exclude "forces" to be calculated by default
#calculator.set_default_properties([])  # If you don't require any default property.
calculator.set_default_properties(["energy"])  # Only require "energy" to be calculated as default.
atoms = molecule("H2")
atoms.calc = calculator
calculator.reset()

print(f"default_properties: {calculator.default_properties}")

print("--- calculate potential energy ---")
%time E_pot1 = atoms.get_potential_energy()
print(f"calculator.results = {calculator.results}")

print("--- calculate forces  ------------")
%time forces = atoms.get_forces()
```

    default_properties: ['energy']
    --- calculate potential energy ---
    CPU times: user 3.51 ms, sys: 0 ns, total: 3.51 ms
    Wall time: 110 ms
    calculator.results = {'energy': -4.5231968186795, 'calc_stats': {'elapsed_usec_preprocess': 3322, 'elapsed_usec': 95427, 'elapsed_usec_infer': 10487, 'n_atoms_plus_ghost_atoms': 2, 'n_neighbors_after_screening': 3413, 'n_neighbors': 1}, 'free_energy': -4.5231968186795}
    --- calculate forces  ------------
    CPU times: user 2.8 ms, sys: 0 ns, total: 2.8 ms
    Wall time: 66.4 ms



TIPS: この設定をした際に`energy`, `forces` 双方の値が欲しい場合、 **`get_forces()` を呼んでから `get_potential_energy()` を呼ぶことで**、最初の `get_forces()` で`energy`, `forces` 双方が計算されているため、 **`get_potential_energy()`の計算をスキップすることができ、効率よくPFPを活用することができます。**


```python
atoms = molecule("H2")
atoms.calc = calculator
calculator.reset()

print("--- calculate forces           ---")
%time forces = atoms.get_forces()
print("calculator.results: ", calculator.results)

print("--- calculate potential energy ---")
%time E_pot1 = atoms.get_potential_energy()
```

    --- calculate forces           ---
    CPU times: user 2.26 ms, sys: 3.66 ms, total: 5.92 ms
    Wall time: 102 ms
    calculator.results:  {'energy': -4.523196252435347, 'forces': array([[-0.        , -0.        ,  0.46052292],
           [-0.        , -0.        , -0.46052292]]), 'calc_stats': {'elapsed_usec_preprocess': 3340, 'elapsed_usec': 88436, 'elapsed_usec_infer': 18462, 'n_atoms_plus_ghost_atoms': 2, 'n_neighbors_after_screening': 9808, 'n_neighbors': 1}, 'free_energy': -4.523196252435347}
    --- calculate potential energy ---
    CPU times: user 297 μs, sys: 0 ns, total: 297 μs
    Wall time: 302 μs



TIPS: この設定を行うことで、Defaultで計算される`charges`を出力しないようにすることもできます。ASEのVersionによっては、charge propertyをもつ`atoms`のxyzファイルへの書込・読込で不具合がある場合もあり、`charge`を使用しない場合は明示的に計算しないことで回避することも可能です。
 - https://github.com/matlantis-pfcc/matlantis-contrib/tree/main/matlantis_contrib_examples/ase_read


### [Advanced] Tips: 複数のPropertyを明示的に同時に計算させる方法

エネルギーやchargeなど、複数のproperty を得たいときにそれらの計算を同時に行ってしまいたい場合があります。&lt;br/&gt;
Calculatorの `calculate` method は第２引数に計算したい`properties`を指定しますが、これを用いることで全ての property を同時に計算させることができます。&lt;br/&gt;
つまり、常に計算される`calculator.default_properties` (前述)に`properties`を追加することができます。&lt;br/&gt;
以下の例では１度の`calculate` 実行で、`results`に全ての結果が含まれていることがわかります。


```python
calculator = ASECalculator(Estimator())

atoms = bulk("Pt") * (2, 2, 1)
atoms.calc = calculator

print("--- calculator.calculate ---")
# All properties are calculated at the same time by using the `calculate` method of the calculator.
# You can see that `results` contains all the results.
%time calculator.calculate(atoms, ["energy", "free_energy", "charges", "forces", "stress"])
print(calculator.results)

# Below, you can see that the calculation does not run when each property is acquired, and the execution is completed immediately.
print("--- E_pot ---")
%time E_pot = atoms.get_potential_energy()
print("--- charges ---")
%time charges = atoms.get_charges()
print("--- forces ---")
%time forces = atoms.get_forces()
print("--- stress ---")
%time stress = atoms.get_stress()
```

    --- calculator.calculate ---
    CPU times: user 1.69 ms, sys: 101 μs, total: 1.79 ms
    Wall time: 52.5 ms
    {'energy': -21.83327464366414, 'forces': array([[ 3.61653257e-07, -7.85071957e-07,  2.52180297e-07],
           [ 3.01579026e-07,  8.62054119e-07, -1.72958580e-06],
           [ 1.68411135e-07, -1.12558125e-06,  1.28098664e-06],
           [-8.31643418e-07,  1.04859909e-06,  1.96418869e-07]]), 'charges': array([-8.16645240e-08,  2.22884921e-07,  8.74260948e-08, -2.28646527e-07]), 'calc_stats': {'elapsed_usec_preprocess': 5360, 'elapsed_usec': 45635, 'elapsed_usec_infer': 16265, 'n_atoms_plus_ghost_atoms': 648, 'n_neighbors_after_screening': 62, 'n_neighbors': 400}, 'free_energy': -21.83327464366414, 'stress': array([-6.11363305e-02, -6.11362696e-02, -6.11363612e-02, -8.66598647e-09,
           -1.47810429e-08,  9.76262246e-09])}
    --- E_pot ---
    CPU times: user 265 μs, sys: 16 μs, total: 281 μs
    Wall time: 285 μs
    --- charges ---
    CPU times: user 376 μs, sys: 0 ns, total: 376 μs
    Wall time: 384 μs
    --- forces ---
    CPU times: user 349 μs, sys: 0 ns, total: 349 μs
    Wall time: 357 μs
    --- stress ---
    CPU times: user 354 μs, sys: 0 ns, total: 354 μs
    Wall time: 360 μs



## PFP Advanced Usage

モデリングを行う際、PFPの細かな挙動を理解しておく必要がある場合の解説を行います。


### Cutoff距離について

PFPでは Graph Neural Networkを用いてエネルギーの推論を行いますが、 v6.0.0以降ではその入力グラフの辺の構築のための**cutoff距離は 9.0A** となっています。
ただし、v5.0.0以前ではcutoff距離は6.0Aです。&lt;br/&gt;
これは、原子が 9.0A 以上離れた場所では相互作用を起こさず、エネルギーは変化しないということです。

ここではCa原子とO原子のポテンシャルエネルギーと Bader Charge がCa-O 原子間距離に応じてどのように変わるかを見てみます。


```python
import numpy as np
from ase import Atoms

dists = np.concatenate((np.linspace(1.4, 4.0, 40, endpoint=False), np.linspace(4.0, 10.0, 20)))

atoms_list = []
E_pot_list = []
charges_list = []
with ASECalculator(Estimator(model_version="v8.0.0", calc_mode=EstimatorCalcMode.PBE)) as calculator:
    for dist in dists:
        atoms = Atoms("CaO", positions=[(0, 0, 0), (0, 0, dist)])
        atoms_list.append(atoms)
        atoms.calc = calculator

        E_pot = atoms.get_potential_energy()
        charges = atoms.get_charges()

        E_pot_list.append(E_pot)
        charges_list.append(charges)

E_pot_array = np.array(E_pot_list)
charges_array = np.array(charges_list)
```


```python
# Visualization of the calculated system
from ase.visualize import view

view(atoms_list, viewer="ngl")
```








    HBox(children=(NGLWidget(max_frame=59), VBox(children=(Dropdown(description='Show', options=('All', 'Ca', 'O')…




```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 1, figsize=(10, 8))

axes[0].plot(dists, E_pot_array)
axes[0].vlines([9.0], np.min(E_pot_array), np.max(E_pot_array), color="gray", linestyles="dashed")
axes[0].set_title("Potential energy")
axes[0].set_ylabel("Energy [eV]")

axes[1].plot(dists, charges_array[:, -1], label=atoms.symbols[-1])  # O's charge
axes[1].vlines([9.0], np.min(charges_array), np.max(charges_array), color="gray", linestyles="dashed")
axes[1].legend()
axes[1].set_title("Bader Charges")
axes[1].set_ylabel("Charge [e]")
axes[1].set_xlabel("Distance [A]")

plt.show()
```



![png](output_76_0.png)




以下のように、カットオフ距離を超えたあとのエネルギー、Bader chargeは一定値となっている事がわかります。

入力構造のモデリングをする際には、反応を見たい系がカットオフ距離内でつながっていることをご確認ください。

原子の位置に対して微分可能にするためにエネルギーやBader chargeはカットオフ距離の直前でほとんどゼロになります。(&lt;a href="https://redirect.matlantis.com/documents/detail/documents/pfp-description/versions/2024.01/index.html#_7" target="_blank"&gt;PFPについて&lt;/a&gt;参照).

※ PFP の出力するcharge はBader chargeを学習した値を出力していますが、これはエネルギーを学習する際に補助的に出力ができる値という扱いとなっています。定量的な精度保証をしている物性値ではないので、参考程度にご覧ください。


[For Advanced Users]

現時点でのPFPでは、系全体の電荷が0 (Neutral) な状態のもののみ推論が可能です。&lt;br/&gt;
ただし局所的に電荷が偏った系のエネルギーの推論を行うことは可能です。

アニオン・カチオンなどの分子を含むシミュレーションを行いたい場合、系全体の電荷が釣り合うようにカウンターイオンなどを置くことでシミュレーションを行う事ができますが、その場合はカウンターイオンもカットオフ距離以内に入れて頂く必要があります。&lt;br/&gt;
こういった使い方はまだ深く検証ができていないところにもなりますので、精度に関して不明点などの気になる挙動が見つかった場合には是非フィードバックをください。次期開発に生かしていきたいと思います。
