# Matlantis-features: ガス分子の熱力学特性と生成エンタルピー


## **内容**
* [0. 初期設定](#advanced_0)
* [1. 熱力学的](#advanced_1)
* [2. 生成エンタルピー](#advanced_2)
* [3. 標準生成エンタルピー](#advanced_3)


&lt;a id='advanced_0'&gt;&lt;/a&gt;
## 0. 初期設定

* 必要なライブラリ等の準備を行います。


```python
# If you have already installed matlantis-features, you can skip this.
# Version 0.25.0 or later is required to run this notebook.
!pip install 'matlantis-features&gt;=0.25.0'
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


&lt;a id='advanced_1'&gt;&lt;/a&gt;
## 1. 熱力学的特徴
* ガスの熱力学特性は、振動周波数がわかっている場合、理想気体の仮定に従って計算できます。
* ガスの熱力学特性は`PostVibrationGasThermoFeature`で計算できます。
* 次のプロパティを取得できます。
     * エンタルピー
     * エントロピー
     * ギブズの自由エネルギー


### 1-1 材料の準備
* メタノール（CH3OH）を使用します。
* 初期構造を最適化します。


```python
from ase.build import molecule
from matlantis_features.features.common.opt import FireASEOptFeature

atoms = molecule("CH3OH")
opt = FireASEOptFeature(estimator_fn=estimator_fn)
opt_result = opt(atoms)
```


### 1-2 振動計算
* メタノール分子の振動数は `VibrationFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/vibration.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;)で計算できます。
* 入力には最適化された分子構造を使用します。


```python
from matlantis_features.features.common.vibration import VibrationFeature

vib = VibrationFeature(estimator_fn=estimator_fn)
vib_result = vib(opt_result.atoms)
```


### 1-3 熱力学特性の計算
* 熱力学特性は、 `PostVibrationGasThermoFeature` (&lt;a href="/api/resource/documents/matlantis-guidebook/ja/matlantis-features/vibration_thermo.html" target="_blank"&gt;Guidebookを見る&lt;/a&gt;) を使用して分子の振動周波数から計算されます。
* 温度と圧力を入力として指定する必要があります。
     * 注意：圧力の単位はeV/angstrom^3です。1 eV/angstrom^3 = 160.21766208 GPa.


```python
from ase import units
from matlantis_features.features.common.gas_thermo import PostVibrationGasThermoFeature

thermo = PostVibrationGasThermoFeature()
thermo_result = thermo(
    vib_result,
    temperatures=[200.0, 250.0, 298.15, 300.0, 350.0, 400.0, 450.0, 500.0],
    pressures=[
        1e5 * units.Pascal,
        1e5 * units.Pascal,
        1e5 * units.Pascal,
        1e5 * units.Pascal,
        1e5 * units.Pascal,
        1e5 * units.Pascal,
        1e5 * units.Pascal,
        1e5 * units.Pascal,
    ],
)
```


```python
thermo_csv = thermo_result.to_dataframe()
print(thermo_csv)
```

       temperature (K)  pressure (eV/Angstrom^3)  pressure (Pa)  pressure (GPa)  \
    0           200.00              6.241509e-07       100000.0          0.0001
    1           250.00              6.241509e-07       100000.0          0.0001
    2           298.15              6.241509e-07       100000.0          0.0001
    3           300.00              6.241509e-07       100000.0          0.0001
    4           350.00              6.241509e-07       100000.0          0.0001
    5           400.00              6.241509e-07       100000.0          0.0001
    6           450.00              6.241509e-07       100000.0          0.0001
    7           500.00              6.241509e-07       100000.0          0.0001

       pressure (bar)  entropy (eV/K)  entropy (J/K/mol)  entropy (kJ/K/mol)  \
    0             1.0        0.002300         221.920419            0.221920
    1             1.0        0.002395         231.090419            0.231090
    2             1.0        0.002476         238.875540            0.238876
    3             1.0        0.002479         239.160710            0.239161
    4             1.0        0.002556         246.586804            0.246587
    5             1.0        0.002628         253.598257            0.253598
    6             1.0        0.002698         260.306148            0.260306
    7             1.0        0.002765         266.764130            0.266764

       entropy (kcal/K/mol)  enthalpy (eV)  enthalpy (J/mol)  enthalpy (kJ/mol)  \
    0              0.053040     -21.037905     -2.029849e+06       -2029.849301
    1              0.055232     -21.016583     -2.027792e+06       -2027.792025
    2              0.057093     -20.994501     -2.025661e+06       -2025.661389
    3              0.057161     -20.993617     -2.025576e+06       -2025.576102
    4              0.058936     -20.968624     -2.023165e+06       -2023.164704
    5              0.060611     -20.941389     -2.020537e+06       -2020.536844
    6              0.062215     -20.911853     -2.017687e+06       -2017.687115
    7              0.063758     -20.880070     -2.014621e+06       -2014.620545

       enthalpy (kcal/mol)  gibbs_free_energy (eV)  gibbs_free_energy (J/mol)  \
    0          -485.145626              -21.497914              -2.074233e+06
    1          -484.653926              -21.615354              -2.085565e+06
    2          -484.144691              -21.732652              -2.096882e+06
    3          -484.124307              -21.737235              -2.097324e+06
    4          -483.547969              -21.863117              -2.109470e+06
    5          -482.919896              -21.992733              -2.121976e+06
    6          -482.238794              -22.125901              -2.134825e+06
    7          -481.505866              -22.262478              -2.148003e+06

       gibbs_free_energy (kJ/mol)  gibbs_free_energy (kcal/mol)
    0                -2074.233384                   -495.753677
    1                -2085.564630                   -498.461910
    2                -2096.882131                   -501.166857
    3                -2097.324315                   -501.272542
    4                -2109.470086                   -504.175451
    5                -2121.976147                   -507.164471
    6                -2134.824882                   -510.235392
    7                -2148.002610                   -513.384945


### 1-4 計算結果の可視化


```python
from plotly.subplots import make_subplots
import plotly.graph_objects as go

fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=("Entropy", "Enthalpy &amp; Gibbs free energy")
)
fig.add_trace(
    go.Scatter(x=thermo_result.temperature["K"], y=thermo_result.entropy["J/K/mol"],name= "S (J/K/mol)"), row=1, col=1
)
fig.add_trace(
    go.Scatter(x=thermo_result.temperature["K"], y=thermo_result.enthalpy["kJ/mol"], name="H (kJ/mol)"), row=1, col=2
)
fig.add_trace(
    go.Scatter(x=thermo_result.temperature["K"], y=thermo_result.gibbs_free_energy["kJ/mol"], name="G (kJ/mol)"), row=1, col=2
)

fig.update_xaxes(title="temperature [K]")
fig.update_yaxes(title="entropy [J/K/mol]", row=1, col=1)
fig.update_yaxes(title="thermo properties [kJ/mol]", row=1, col=2)
```


### 1-5 結果と実験値の比較
* メタノールの実験的エントロピーは、標準条件で239.9 J /(mol K) です（https://en.wikipedia.org/wiki/Methanol_(data_page) による）。
* 同じ条件での計算結果は以下のとおりです。


```python
print("Standard molar entropy of methanol:")
print("    Expt: 239.9 J/K/mol)")
print(f"    Calc: {thermo_result.entropy['J/K/mol'][2]:.1f} J/K/mol")

```

    Standard molar entropy of methanol:
        Expt: 239.9 J/K/mol)
        Calc: 238.9 J/K/mol


&lt;a id='advanced_2'&gt;&lt;/a&gt;
## 2. 生成エンタルピー
* 生成エンタルピーは、1 molの特定の材料と、そのすべての構成元素が標準状態にある場合のエンタルピーの差として定義されます。

* matlantis-featuresでは、ガス分子の生成エンタルピーを計算するために`ComplexGasFormationEnthalpyFeature`が提供されています。
* `ComplexGasFormationEnthalpyFeature`は、内部で`FireASEOptFeature`、 `VibrationFeature`、`PostVibrationGasThermoFeature`、 `ForceConstantFeature`、および`PostPhononThermochemistryFeature`を呼び出して、ガス分子とその構成元素のエンタルピーを計算します。
* 現在、この機能はC、H、OとNによって形成されるガス分子にのみ使用できます。
* 構成元素、つまりC（グラファイト）、H（H2ガス）、O（O2ガス）、N（N2ガス）のエンタルピーは内部で計算されます。
    * H、OおよびNのエンタルピーは、内部で `VibrationFeature`および` PostVibrationGasThermoFeature`を使用して計算されます。
    * グラファイトのエンタルピーは、H = U + PV'で近似されます。ここで、Uは`ForceConstantFeature`および` PostPhononThermochemistryFeature`から取得され、V'は0Kおよび0Paでのグラファイトの体積です。熱膨張の影響は計算では無視します。


### 2-1 材料の準備
* メタノール（CH3OH）を使用します。
* `ComplexGasFormationEnthalpyFeature`は入力構造を内部的に最適化するため、構造最適化を実行する必要はありません。


```python
from ase.build import molecule

atoms = molecule("CH3OH")
```


### 2-2 生成エンタルピーの計算
* 生成エンタルピーは `ComplexGasFormationEnthalpyFeature`で簡単に計算できます
* `PostVibrationGasThermoFeature`と同じように、温度の単位はKで、圧力の単位はeV/angstrom^3です。1 eV/angstrom^3 = 160.21766208 GPa.


```python
from ase import units
from matlantis_features.features.common import (
    FireASEOptFeature,
    VibrationFeature,
    ComplexGasFormationEnthalpyFeature,
)

formation_enthalpy = ComplexGasFormationEnthalpyFeature(FireASEOptFeature(estimator_fn=estimator_fn), VibrationFeature(estimator_fn=estimator_fn))
formation_enthalpy_result = formation_enthalpy(atoms, temperatures=[298.15, 500.0], pressures=[1e5 * units.Pascal, 1e5 * units.Pascal])
```


### 2-3 結果と実験値の比較
* メタノールの実験標準生成エンタルピーは -205 ±10.kJ/mol です。（NISTデータベース https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&amp;Units=SI&amp;Mask=1#Thermo-Gas による）。
* 標準状態での計算結果は以下のとおりです。


```python
print("Standard formation enthalpy of methanol:")
print("    Expt: -205. ± 10. kJ/mol")
print(f"    Calc: {formation_enthalpy_result.formation_enthalpy['kJ/mol'][0]:.1f} kJ/mol")
```

    Standard formation enthalpy of methanol:
        Expt: -205. ± 10. kJ/mol
        Calc: -205.5 kJ/mol


&lt;a id='advanced_3'&gt;&lt;/a&gt;
## 3. 標準生成エンタルピー
* 実験的な標準生成エンタルピーとのより良い一致を達成するために、 `GasStandardFormationEnthalpyFeature`がmatlantis-featuresに提供されています。
* `ComplexGasFormationEnthalpyFeature`との違いは、標準状態のC（グラファイト）、H（H2ガス）、O（O2ガス）、およびN（N2ガス）のエンタルピーを実験による標準生成エンタルピーデータセットにフィッティングしていることです。
* データセットはNISTデータベースから取得され、220個のガス分子が含まれています。
* この機能を使用して計算できるのは標準生成エンタルピーのみであり、他の温度と圧力での生成エンタルピーは計算できません。
* **`GasStandardFormationEnthalpyFeature`は現時点ではPFPバージョン2.0.0以上をサポートしています。**使用する前にPFPバージョンを設定してください。


### 3-1 材料の準備
* メタノール（CH3OH）を使用します。
* `ComplexGasFormationEnthalpyFeature`は入力構造を内部的に最適化するため、構造最適化を実行する必要はありません。


```python
from ase.build import molecule

atoms = molecule("CH3OH")
```


### 3-2 標準生成エンタルピーの計算


```python
from matlantis_features.features.common import (
    FireASEOptFeature,
    VibrationFeature,
    GasStandardFormationEnthalpyFeature,
)

standard_formation_enthalpy = GasStandardFormationEnthalpyFeature(FireASEOptFeature(estimator_fn=estimator_fn), VibrationFeature(estimator_fn=estimator_fn))
standard_formation_enthalpy_result = standard_formation_enthalpy(atoms)
```


### 3-3 結果と実験値の比較
* メタノールの実験標準生成エンタルピーは -205 ±10.kJ/mol です。（NISTデータベース https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561&amp;Units=SI&amp;Mask=1#Thermo-Gas による）。
* 標準状態での計算結果は以下のとおりです。


```python
print("Standard formation enthalpy of methanol:")
print("    Expt: -205. ± 10. kJ/mol)")
print(f"    Calc: {standard_formation_enthalpy_result.standard_formation_enthalpy['kJ/mol'][0]:.1f} kJ/mol")
```

    Standard formation enthalpy of methanol:
        Expt: -205. ± 10. kJ/mol)
        Calc: -199.4 kJ/mol

