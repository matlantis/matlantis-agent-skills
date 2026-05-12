# Matlantis CSP について ​

Matlantis CSP はMatlantis上でPFPを用いた結晶構造探索(Crystal Structure Prediction, CSP)を行うための製品です。 Matlantis 上で提供されているPythonライブラリ`mtcsp`を用いて結晶構造探索を行うことができます。

bash
    
    
    pip install -U mtcsp[matlantis]

以下では、Matlantis CSP の提供する結晶構造探索の技術詳細を概説します。 より詳細な技術情報は[テクニカルペーパー](https://arxiv.org/abs/2503.21201)を参照下さい。

## 結晶構造探索 ​

結晶構造探索(Crystal Structure Prediction, CSP)とは、与えられた元素組成のもとでエネルギー的に安定な結晶構造を予測する手法です。 特定の元素系における安定相の全体像を把握したい場合に有効であり、相分離によって現れうる相も含めて、エネルギー的に妥当な候補を体系的に洗い出すことができます。 また、所望の物性をもつ候補物質が、その元素系でどの程度エネルギー的に安定かを定量化する用途にも適しており、既知の結晶構造を母構造として元素置換したとき、相の安定性が損なわれない組成範囲を見つもるためにも用いられます。

## 遺伝的アルゴリズムに基づく結晶構造探索 ​

Matlantis CSP では、遺伝的アルゴリズムに基づく結晶構造探索を採用しています。 遺伝的アルゴリズムの枠組みでは、各世代ごとに構造の集団(population)を用意して、その中から次の世代を生成するために望ましいelite populationを選択します。 選択されたelite populationの構造の中から、交叉(crossover)と突然変異(mutation)を行い、新たな構造を生成します。

### ランダム構造生成 ​

探索の初期集団はランダムに結晶構造を生成して用意します。 ランダム構造はLyakhov et al. 2010 [1]に基づいて、まずセルを複数のサブセルに分割して、そのうちの1つのサブセル内で原子をランダム配置したものを他のサブセルに複製することで生成しています。 サブセル内でのランダム配置には、[PyXtal](https://github.com/MaterSim/PyXtal) [2] を利用しています。 PyXtalによるランダム構造生成では組成と空間群の指定が必要ですが、組成については指定した最大原子数(`max_atoms`)以下までの原子数を一様分布から選び、空間群についてはP1以外の空間群タイプを一様分布から選んでいます。

### 交叉と突然変異 ​

交叉は二つの親構造を組み合わせて子構造を作る操作、突然変異は親構造に変化を与えて子構造を得る操作を指します。 Matlantis CSP では、下表に示すような交叉と突然変異の手法を用いています。

交叉・突然変異の手法| 説明| モジュール  
---|---|---  
`VariableCutAndSplicePairing`| 2つの親構造をランダムな結晶面に沿って切断し、それらを接合して子構造を生成する| `pfn_ase_extras.genetic_operations`  
`RandomCompositionMutationAllowingSingleElement`| 親構造の組成をランダムに変化させて子構造を生成する| `pfn_ase_extras.genetic_operations`  
`PermutationMutationAllowingSingleElement`| 親構造の2つの原子位置をランダムに入れ替える| `pfn_ase_extras.genetic_operations`  
`RemoveRandomAtomMutation`| 親構造の原子をランダムに削除して子構造を生成する| `mtcsp.samplers.genetic_operations`  
`GenerateRandomStructure`| ランダム構造生成を用いて子構造を生成する| `mtcsp.samplers.genetic_operations`  
`ExpandCellAndReplaceAtomsMutation`| 親構造のセルを拡大して一部の原子種を変化させて子構造を生成する| `mtcsp.samplers.genetic_operations`  
`StrainMutation`| 親構造のセルにランダムにひずみを加えて子構造を生成する| `ase_ga.standardmutations`  
`RattleMutation`| 親構造の原子位置にランダムに変位を加えて子構造を生成する| `ase_ga.standardmutations`  
  
これらの交叉・突然変異手法は`ase_ga.offspring_creator.OffspringCreator`の継承クラスとして実装されています。 `ase_ga`パッケージに含まれる手法は[ASE-GA](https://dtu-energy.github.io/ase-ga/)で標準的に提供されている実装です。 [pfn-ase-extrasパッケージ](https://pypi.org/project/pfn-ase-extras/)に含まれる手法はLGPL-2.1-or-laterライセンスで配布されており、ASE-GAのライセンスと互換性があります。 `mtcsp.samplers.genetic_operations`に含まれる手法はMatlantis CSPの内部実装であり、将来のリリースで変更される可能性があります。

## PFPによる形成エネルギー計算 ​

遺伝的アルゴリズムによって生成された構造はPFPを用いてセルと原子位置の両方を構造緩和して、形成エネルギーを計算します。

### 形成エネルギー補正と計算モードの選ばれ方 ​

PFPが複数の計算モード(`EstimatorCalcMode`)をもっていることに対応して、Matlantis CSP では2種類の形成エネルギー補正(`pfp_compatibility`)を用意しています。

#### `pfp_compatibility=pfp2023compatibility` ​

この補正モードでは、[Materials ProjectのGGA/GGA+U mixing 補正およびアニオン補正](https://docs.materialsproject.org/methodology/materials-methodology/calculation-details) [3]に概ね対応する計算モードを選択します。 具体的には、V, Cr, Mn, Fe, Co, Ni, Cu, Mo, W のいずれか遷移金属を含んだ酸化物およびフッ化物に対しては`PBE_U`でエネルギー推論を行い、その他に対しては`PBE`でエネルギー推論を行います。

Cu酸化物はMatlantis CSP での補正がMaterials Projectでの補正と異なる唯一の例外です。 Matlantis CSP ではCu酸化物をHubbard補正ありに相当する`PBE_U`でエネルギー推論しますが、Materials Project ではHubbard補正なしの計算結果を用いています。

#### `pfp_compatibility=r2scan` ​

この補正モードでは、PFPの`EstimatorCalcMode.R2SCAN`を常に用いてエネルギー推論を行います。 R2SCANモードが利用可能なPFP v8.0.0以降のバージョンで利用可能です。

#### `pfp_compatibility=pbe` ​

この補正モードでは、PFPの`EstimatorCalcMode.PBE`を常に用いてエネルギー推論を行います。 アニオン補正やGGA/GGA+U mixing補正はこのモードでは適用されません。

#### `pfp_compatibility=r2scan_plus_d3` ​

この補正モードでは、PFPの`EstimatorCalcMode.R2SCAN_PLUS_D3`を常に用いてエネルギー推論を行います。 単体の基準構造には、R2SCANモードのPFPで最適化した単体構造を用います。

#### `pfp_compatibility=pbe_plus_d3` ​

この補正モードでは、PFPの`EstimatorCalcMode.PBE_PLUS_D3`を常に用いてエネルギー推論を行います。 単体の基準構造には、PBEモードのPFPで最適化した単体構造を用います。

### 単体の基準構造 ​

形成エネルギーの基準となる単体構造は、探索に用いるPFPのバージョン(`pfp_version`)と形成エネルギー補正モード(`pfp_compatibility`)に応じて、事前計算済の同一の構造が用いられます。 これらの単体構造には、Materials Project [4]とAFLOW standard encyclopedia of crystallographic prototypes [5]の構造を候補として、指定したPFPのバージョンと計算モードのもとで構造最適化を行い、エネルギー最小の構造を選択しています。

## サンプルノートブックの概要 ​

Matlantis CSP は、多数のサンプル Jupyter Notebook を提供しており、**Quickstart** 、**Tutorials** 、**How-to ガイド** の 3 つのカテゴリに分類されています。

### Quickstart ​

  * [Quickstart](/api/docs/ja/examples/1_quickstart.html)



Quickstart では、最小構成のエンドツーエンドな結晶構造探索を実行します。 このチュートリアルを完了すると、小規模な CSP 実行結果、保存された候補構造、そして組成空間全体における熱力学的安定性を可視化した相図を得ることができます。

### Tutorials ​

Tutorials では、**mtcsp** の中核機能について、代表的な CSP ワークフローおよび解析タスクを中心に、段階的に解説します。

#### 結晶構造探索（制約あり／なし）と解析 ​

  * [指定した元素系全体の結晶構造探索](/api/docs/ja/examples/2_1_convex_hull_search.html)
  * [指定した元素系全体の結晶構造探索の解析](/api/docs/ja/examples/2_2_convex_hull_analysis.html)
  * [組成指定を含んだ制約付き結晶構造探索](/api/docs/ja/examples/2_3_convex_hull_search_with_constraints.html)



これらのチュートリアルでは、システム全体を対象とした 結晶構造探索と結果解析を行います。 結晶構造探索を完了し、生成された構造を保存するとともに、安定相の相図を作成し、主要な凸包構造を表形式で抽出します。 オプションとして、Materials Project の凸包との比較も紹介します。

さらに、組成空間や酸化数範囲の制限方法、凸包体積の推移を用いた収束判定、結果を再利用するためのエクスポート方法についても学びます。

#### 置換構造探索 ​

  * [指定した母構造からの置換構造探索](/api/docs/ja/examples/3_1_derivative_structure_search.html)
  * [指定した母構造からの置換構造探索の解析](/api/docs/ja/examples/3_2_derivative_structure_analysis.html)
  * [組成とスーパーセルを限定した置換構造探索](/api/docs/ja/examples/3_3_derivative_structure_search_with_composition_and_supercell.html)



これらのチュートリアルでは、母構造から置換構造を生成し、解析する手法を解説します。 置換構造探索を完了し、組成ごとの低エネルギー候補を抽出するとともに、置換組成空間における安定性分布を可視化します。

また、組成比の固定、スーパーセルサイズの選択、大規模な探索空間における列挙法とランダムサンプリングの使い分けなど、探索範囲を制御する手法や、結果のエクスポートについても扱います。

#### 平衡電圧計算 ​

  * [置換構造探索を用いた平衡電位計算](/api/docs/ja/examples/4_1_cell_voltage_search.html)
  * [置換構造探索を用いた平衡電位計算の解析](/api/docs/ja/examples/4_2_cell_voltage_analysis.html)



これらのチュートリアルでは、置換構造探索の結果を用いて平衡電圧プロファイルを計算します。 LiCoO2系を対象とした探索を実行し、組成ごとに低エネルギー構造を選別したうえで、Li 空孔率に対する電圧曲線を生成します。

さらに、置換構造の凸包から安定な終端組成を特定し、詳細解析のために代表構造を抽出する方法も解説します。

### How-to ガイド ​

  * [指定した元素系全体の結晶構造探索(hull search)の再開](/api/docs/ja/examples/5_1_restart_experiment.html)
  * [既存の探索結果を初期構造に追加する](/api/docs/ja/examples/5_2_append_experiment.html)



How-to ガイドでは、**mtcsp** の利用者が必要とする具体的な操作手順を説明します。

  1. **実験の再開** 既存の結晶構造探索結果を読み込み、過去の結果を保持したまま追加の構造計算を実行します。

  2. **既存の探索結果を用いた初期化** 過去の実験で得られた構造を用いて新しい探索を初期化します。




## 参考文献 ​

* * *

  1. A. O. Lyakhov, A. R. Oganov, and M. Valle, [Comput. Phys. Commun. 181, 1623 (2010)](https://www.sciencedirect.com/science/article/abs/pii/S0010465510001840?via%3Dihub). ↩︎

  2. S. Fredericks, K. Parrish, D. Sayre, and Q. Zhu, [Comput. Phys. Commun. 261, 107810 (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0010465520304057?via%3Dihub). ↩︎

  3. A. Wang, R. Kingsbury, M. McDermott, M. Horton, A. Jain, S. P. Ong, S. Dwaraknath, and K. A. Persson, [Sci. Rep. 11, 15496 (2021)](https://www.nature.com/articles/s41598-021-94550-5). ↩︎

  4. A. Jain, S. P. Ong, G. Hautier, W. Chen, W. D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, and K. A. Persson, [APL Mater. 1, 011002 (2013)](https://pubs.aip.org/aip/apm/article/1/1/011002/119685/Commentary-The-Materials-Project-A-materials). ↩︎

  5. M. J. Mehl, D. Hicks, C. Toher, O. Levy, R. M. Hanson, G. Hart, and S. Curtarolo, [Comput. Mater. Sci. 136, S1 (2017)](https://www.sciencedirect.com/science/article/abs/pii/S0927025617300241?via%3Dihub). ↩︎



