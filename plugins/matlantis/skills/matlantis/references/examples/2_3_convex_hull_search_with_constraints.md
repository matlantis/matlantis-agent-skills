# 組成指定を含んだ制約付き結晶構造探索

このチュートリアルでは、組成指定、組成範囲、組成から推定される酸化状態範囲を含む制約付き結晶構造予測(CSP)の例を示します。


## インストール


```python
# Workaround: prevents occasional hang in exception display when using mtcsp
%xmode Plain
%pip install -U mtcsp[matlantis]
```


```python
from fractions import Fraction

from mtcsp.convex_hull_search import SearchConfig
```


## 探索条件の設定

`SearchConfig`を設定することで、探索する組成を制限することができます。
各制約は`SearchConfig`で指定することで、同時に適用することが可能です。
詳細はMTCSPのAPIリファレンスを参照してください。


### 探索する組成範囲の制限

`SearchConfig.composition_range`は、探索する各元素種の組成の下限と上限を指定します。
組成は0から1の間で、全元素の分数組成の和が1になるように指定してください。
端点が許容範囲に含まれるように、また浮動小数点誤差を避けるために、`fractions.Fraction`を使用してください。

デフォルトでは指定されず、組成空間全体を探索します。


```python
search_config = SearchConfig(
    elements=["Ti", "O"],
    experiment_name="dummy-experiment-name",
    # Limit composition Ti_{x}O_{y} to 1/3&lt;= x &lt;= 1/2 and 0 &lt;= y &lt;= 1.
    # Please specify composition by fractions.Fraction to avoid floating point errors.
    composition_range={"Ti": (Fraction(1, 3), Fraction(1, 2)), "O": (Fraction(0), Fraction(1))},
)
```


### 指定した組成間の内挿による探索組成範囲の制限

`SearchConfig.composition_endpoints`は、探索する組成空間の端点を指定します。
例えば、`["SrO", "SrO2", "TiO2"]`と指定すると、SrO-SrO2-TiO2の三角形(2~4価のSrと4価のTiから成る)の組成空間を探索します。

デフォルトでは指定されず、その場合組成空間全体を探索します。


```python
search_config = SearchConfig(
    elements=["Ti", "O"],
    experiment_name="dummy-experiment-name",
    # Limit composition to the triangle SrO-SrO2-TiO2.
    # Please specify compositions as Composition-like strings or pymatgen.Composition.
    composition_endpoints=["SrO", "SrO2", "TiO2"],
)
```


`SearchConfig.composition_endpoints`に以下のように設定することで、探索する組成を直接指定することもできます。
組成を指定すると、交叉などの遺伝的アルゴリズムの操作で目的の組成の構造が生成されない場合があり、探索効率が低下する可能性があることにご注意ください。


```python
search_config = SearchConfig(
    elements=["Sr", "Ti", "O"],
    experiment_name="dummy-experiment-name",
    # Please specify a composition as a Composition-like string or pymatgen.Composition.
    composition_endpoints=["SrTiO3"],
)
```


```python
search_config = SearchConfig(
    elements=["Sr", "Ti", "O"],
    experiment_name="dummy-experiment-name",
    # Limit compositions to those satisfying oxidation state ranges.
    oxidation_state_ranges={"Sr": (2, 2), "Ti": (3, 4), "O": (-2, -2)},
)
```
