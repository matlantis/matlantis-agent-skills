# [Advanced] Nano particle energy

本節は発展編の位置付けとします。次章へ急ぎたい方はスキップし、後から読んでいただいても構いません。

本節では、Nano particle (クラスター構造) の過剰エネルギー(excess energy)計算と、担体上の吸着構造作成を行ってみます。


```python
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pfcc_extras.visualize.view import view_ngl


print(f"pfp_api_client: {pfp_api_client.__version__}")

estimator = Estimator(calc_mode=EstimatorCalcMode.PBE, model_version="v8.0.0")
calculator = ASECalculator(estimator)
```


    


    pfp_api_client: 1.23.1



```python
from ase.cluster import Decahedron, Icosahedron, Octahedron, wulff_construction
from ase import Atoms 
from ase.build import bulk
from ase.constraints import ExpCellFilter, StrainFilter
from ase.optimize import LBFGS
from ase.io import Trajectory
import numpy as np
import pandas as pd
```

## 過剰エネルギー - Excess energy

純金属に1種類以上他の元素を混ぜた材料を**合金 (alloy)** と呼びますが、
金属は2種類以上の元素を組み合わせることで、その性質を大きく変化させるようなものがあります。

例えば、アルミニウムと銅やマグネシウムを混ぜた、[ジュラルミン](https://ja.wikipedia.org/wiki/%E3%82%B8%E3%83%A5%E3%83%A9%E3%83%AB%E3%83%9F%E3%83%B3)はその軽くて硬い性質から航空宇宙機器などに使われています。

2種類の元素が混ざったものを2元系合金、3種類以上の元素が混ざったものは多元系合金と呼ばれます。

多元系合金は、新たな性質を持つような材料の発見も期待されていますが、
元素の種類が増えるほどその組成比の組み合わせが大きくなり、まだまだ解析が進んでいない領域です。&lt;br/&gt;
このような系も[Matlantis](https://matlantis.com)ではその汎用性を活かし、扱うことができます。


こういった合金の原子構造を考える際には、現実ではどのように混ざるのか、そもそもきれいに混ざることができる組み合わせなのかを知る必要があります。

このような合金の安定度合いを測る指標として、過剰エネルギー(excess energy)があります。&lt;br/&gt;
元素Aと元素Bを混ぜた2元系合金 (alloy)の過剰エネルギー(excess energy)は以下のように計算され、その合金のエネルギーが、元素A, Bがそれぞれ混ざっていなかった場合と比べてどのくらい安定かを示します。

$$ E_{\rm{excess}} = \frac{1}{N_{\rm{alloy}}} \left( E_{\rm{alloy}} - \frac{N_{\rm{alloyA}}}{N_{\rm{alloy}}} E_{\rm{A}} - \frac{N_{\rm{alloyB}}}{N_{\rm{alloy}}} E_{\rm{B}} \right) $$

$E_{\rm{alloy}}$が合金のエネルギー、$E_\rm{A}, E_\rm{B}$ は単原子AまたはBで構成された時のエネルギーで、$N_{\rm{alloy}}, N_{\rm{alloyA}}, N_{\rm{alloyB}}$はそれぞれ合金のすべての原子数、元素Aの原子数、元素Bの原子数です。

多元系の場合も同様に定義できます。

ここでは、PtとPdを混ぜた構造に対して、excess energyを評価することでその安定性を評価してみましょう。

参考文献:

 - [Electronic Structure and Phase Stability of PdPt Nanoparticles](https://pubs.acs.org/doi/10.1021/acs.jpclett.5b02753)
 - [Electronic structure and phase stability of Pt3M (M = Co, Ni, and Cu) bimetallic nanoparticles](https://www.sciencedirect.com/science/article/abs/pii/S0927025620303657)
 - [Calculations of Real-System Nanoparticles Using Universal Neural Network Potential PFP](https://arxiv.org/abs/2107.00963)

Pt711 Nano particleに対して、Pdを混ぜた構造に対して、
2元系のNano particleを考えた場合にもその混ざり方は様々な可能性が考えられます。

今回は以下のような構造を作成し、そのexcess energyを評価してみます。

 - PtとPdが均等に混ざる形で存在する構造
 - Ptが内殻に、Pdが外殻に存在する構造
 - Pdが内殻に、Ptが外殻に存在する構造

Core shell構造を作成する関数定義です。読み飛ばしていただいて構いません。


```python
from typing import List, Tuple

from ase import Atoms, neighborlist
from ase.cluster import Cluster
import numpy as np
from ase.data import atomic_numbers


def cluster2atoms(cluster: Cluster) -&gt; Atoms:
    """Convert ASE Cluster to ASE Atoms

    Args:
        cluster (Cluster): input cluster instance

    Returns:
        atoms (Atoms): converted output, atoms instance
    """
    return Atoms(cluster.symbols, cluster.positions, cell=cluster.cell)


def calc_coordination_numbers(atoms: Atoms) -&gt; Tuple[List[int], List[np.ndarray]]:
    """Calculates coordination number

    Args:
        atoms: input atoms

    Returns:
        cns (list): coordination number for each atom
        bonds (list): bond destination indices for each atom
    """
    nl = neighborlist.NeighborList(
        neighborlist.natural_cutoffs(atoms), self_interaction=False, bothways=True
    )
    nl.update(atoms)
    bonds = []
    cns = []
    for i, _ in enumerate(atoms):
        indices, offsets = nl.get_neighbors(i)
        bonds.append(indices)
        cns.append(len(indices))
    return cns, bonds


def substitute(scaffold: Atoms, sites: np.ndarray, elements: List[str]) -&gt; Atoms:
    """Substibute `elements` to `sites` indices of `scaffold` atoms

    Args:
        scaffold (Atoms): Original input atoms
        sites (np.ndarray): site indices to substitute `elements`
        elements (list): elements to be substituted

    Returns:
        substituted (Atoms): substituted ase atoms
    """
    substituted = scaffold.copy()
    for site, element in zip(sites, elements):
        substituted.numbers[site] = atomic_numbers[element]
    return substituted


def make_core_shell(scaffold: Atoms, element: str) -&gt; Atoms:
    """Make core shell structure

    Input `scaffold` element is kept inside core shell,
    and outside surface is replaced by `element`.

    Note that it is not fully tested.
    It was not guaranteed to work on all the cases.

    Args:
        scaffold (Atoms): input cluster structure
        element (str): outside surface element to be replaced

    Returns:
        coreshell (Atoms): core shell structure
    """
    # CN : vertex &lt; edge &lt; surface &lt; core
    cn = np.array(calc_coordination_numbers(scaffold)[0])
    cn_set = np.unique(cn)
    vertexes = None
    surfaces = None
    cores = None
    if len(cn_set) == 1:
        pass
    elif len(cn_set) == 2:
        cores = np.where(cn == cn_set[1])[0]
        surfaces = np.where(cn == cn_set[0])[0]
    elif len(cn_set) == 3 or len(cn_set) == 4:
        cores = np.where(cn == cn_set[-1])[0]
        surfaces = np.where(cn != cn_set[-1])[0]
        vertexes = np.where(cn == cn_set[0])[0]
    else:
        cores = np.where(cn == cn_set[-1])[0]
        surfaces = np.where(cn != cn_set[-1])[0]
        vertexes = np.where(np.isin(cn, cn_set[:-3]))[0]

    if surfaces is None:
        core_shell = scaffold.copy()
    else:
        core_shell = substitute(scaffold, surfaces, [element] * len(surfaces))
    return core_shell

```

まず今回ベースとなるNano particle 骨格構造(scaffold)を作成します。


```python
Pt711 = cluster2atoms(Octahedron("Pt", 11, cutoff=4))
view_ngl(Pt711)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pt'), value='All'), D…



`make_core_shell` 関数を用いて、ベース骨格の内殻のみ "Pd" 元素に置き換えた構造 Pd306Pt405構造を作成します。


```python
Pd306Pt405 = make_core_shell(Pt711, "Pd")
view_ngl(Pd306Pt405)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pd', 'Pt'), value='Al…



同様にして、外殻がPd, 内殻がPtとなる構造Pd405Pt306を作成します。


```python
Pd711 = Octahedron("Pd", 11, cutoff=4)
Pd405Pt306 = make_core_shell(Pd711, "Pt")
view_ngl(Pd405Pt306)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pd', 'Pt'), value='Al…




```python
symbols, counts = np.unique(Pd405Pt306.symbols, return_counts=True)
print(f"Pd405Pt306 structure contains {symbols} with counts {counts}", )

symbols, counts = np.unique(Pd306Pt405.symbols, return_counts=True)
print(f"Pd306Pt405 structure contains {symbols} with counts {counts}", )
```

    Pd405Pt306 structure contains ['Pd' 'Pt'] with counts [405 306]
    Pd306Pt405 structure contains ['Pd' 'Pt'] with counts [306 405]


ここで得られたそれぞれの構造を構造緩和し、エネルギーを求めます。


```python
from ase.optimize import FIRE

def get_opt_energy(atoms, fmax=0.001):    
    opt = FIRE(atoms)
    opt.run(fmax=fmax)
    return atoms.get_total_energy()
```


```python
print("Optimizing Pt711")
Pt711.calc = calculator
E_pt711 = get_opt_energy(Pt711)

print("Optimizing Pd306Pt405")
Pd306Pt405.calc = calculator
E_pd306pt405 = get_opt_energy(Pd306Pt405)

print("Optimizing Pd405Pt306")
Pd405Pt306.calc = calculator
E_pd405pt306 = get_opt_energy(Pd405Pt306)

print("Optimizing Pd711")
Pd711.calc = calculator
E_pd711 = get_opt_energy(Pd711)
```

    Optimizing Pt711
          Step     Time          Energy          fmax
    FIRE:    0 05:06:26    -3609.569290        1.040004
    FIRE:    1 05:06:26    -3610.618085        0.997116
    FIRE:    2 05:06:26    -3612.416157        0.909194
    FIRE:    3 05:06:26    -3613.983257        0.812325
    FIRE:    4 05:06:26    -3615.361402        0.709797
    FIRE:    5 05:06:26    -3616.582632        0.603863
    FIRE:    6 05:06:27    -3617.676292        0.498381
    FIRE:    7 05:06:27    -3618.662917        0.399323
    FIRE:    8 05:06:27    -3619.553256        0.339693
    FIRE:    9 05:06:27    -3620.347276        0.311894
    FIRE:   10 05:06:27    -3621.047967        0.285431
    FIRE:   11 05:06:27    -3621.657066        0.239392
    FIRE:   12 05:06:27    -3622.182567        0.213061
    FIRE:   13 05:06:27    -3622.631155        0.208877
    FIRE:   14 05:06:27    -3623.011035        0.195828
    FIRE:   15 05:06:27    -3623.328874        0.174904
    FIRE:   16 05:06:28    -3623.590156        0.149677
    FIRE:   17 05:06:28    -3623.795950        0.130845
    FIRE:   18 05:06:28    -3623.951072        0.119098
    FIRE:   19 05:06:28    -3624.068673        0.110227
    FIRE:   20 05:06:28    -3624.158420        0.120910
    FIRE:   21 05:06:28    -3624.221578        0.121859
    FIRE:   22 05:06:28    -3624.259647        0.101970
    FIRE:   23 05:06:28    -3624.277939        0.124756
    FIRE:   24 05:06:29    -3624.355435        0.093453
    FIRE:   25 05:06:29    -3624.438236        0.059104
    FIRE:   26 05:06:29    -3624.504409        0.046673
    FIRE:   27 05:06:29    -3624.540790        0.049031
    FIRE:   28 05:06:29    -3624.549812        0.047167
    FIRE:   29 05:06:29    -3624.554610        0.040623
    FIRE:   30 05:06:29    -3624.562142        0.030064
    FIRE:   31 05:06:29    -3624.569143        0.024587
    FIRE:   32 05:06:29    -3624.574742        0.025425
    FIRE:   33 05:06:29    -3624.578666        0.024004
    FIRE:   34 05:06:29    -3624.581829        0.022947
    FIRE:   35 05:06:30    -3624.584106        0.019925
    FIRE:   36 05:06:30    -3624.585299        0.015150
    FIRE:   37 05:06:30    -3624.585321        0.014810
    FIRE:   38 05:06:30    -3624.585917        0.014117
    FIRE:   39 05:06:30    -3624.586723        0.013142
    FIRE:   40 05:06:30    -3624.587187        0.011900
    FIRE:   41 05:06:30    -3624.587886        0.010466
    FIRE:   42 05:06:30    -3624.588935        0.009258
    FIRE:   43 05:06:30    -3624.589417        0.008267
    FIRE:   44 05:06:30    -3624.590290        0.007319
    FIRE:   45 05:06:31    -3624.590656        0.006103
    FIRE:   46 05:06:31    -3624.591165        0.006491
    FIRE:   47 05:06:31    -3624.591811        0.006671
    FIRE:   48 05:06:31    -3624.592166        0.007046
    FIRE:   49 05:06:31    -3624.592468        0.007752
    FIRE:   50 05:06:31    -3624.593492        0.008142
    FIRE:   51 05:06:31    -3624.594433        0.007817
    FIRE:   52 05:06:31    -3624.595124        0.005346
    FIRE:   53 05:06:31    -3624.595783        0.005861
    FIRE:   54 05:06:31    -3624.596135        0.006223
    FIRE:   55 05:06:32    -3624.596072        0.005749
    FIRE:   56 05:06:32    -3624.595873        0.004870
    FIRE:   57 05:06:32    -3624.596144        0.003680
    FIRE:   58 05:06:32    -3624.596099        0.002518
    FIRE:   59 05:06:32    -3624.596321        0.002303
    FIRE:   60 05:06:32    -3624.596206        0.002044
    FIRE:   61 05:06:32    -3624.596292        0.001735
    FIRE:   62 05:06:32    -3624.596192        0.001515
    FIRE:   63 05:06:32    -3624.596448        0.001490
    FIRE:   64 05:06:33    -3624.596260        0.001484
    FIRE:   65 05:06:33    -3624.596439        0.001470
    FIRE:   66 05:06:33    -3624.596298        0.001456
    FIRE:   67 05:06:33    -3624.596388        0.001447
    FIRE:   68 05:06:33    -3624.596480        0.001421
    FIRE:   69 05:06:33    -3624.596424        0.001405
    FIRE:   70 05:06:33    -3624.596384        0.001372
    FIRE:   71 05:06:33    -3624.596381        0.001336
    FIRE:   72 05:06:33    -3624.596433        0.001274
    FIRE:   73 05:06:33    -3624.596414        0.001209
    FIRE:   74 05:06:33    -3624.596202        0.001148
    FIRE:   75 05:06:33    -3624.596449        0.001080
    FIRE:   76 05:06:34    -3624.596358        0.001115
    FIRE:   77 05:06:34    -3624.596519        0.001125
    FIRE:   78 05:06:34    -3624.596378        0.000971
    Optimizing Pd306Pt405
          Step     Time          Energy          fmax
    FIRE:    0 05:06:34    -3125.594099        0.531633
    FIRE:    1 05:06:34    -3125.926152        0.501638
    FIRE:    2 05:06:34    -3126.506870        0.442726
    FIRE:    3 05:06:34    -3127.202026        0.357763
    FIRE:    4 05:06:34    -3127.886976        0.253067
    FIRE:    5 05:06:34    -3128.431623        0.208540
    FIRE:    6 05:06:34    -3128.884657        0.214499
    FIRE:    7 05:06:35    -3129.273722        0.217996
    FIRE:    8 05:06:35    -3129.607956        0.177662
    FIRE:    9 05:06:35    -3129.881719        0.119569
    FIRE:   10 05:06:35    -3130.089639        0.141494
    FIRE:   11 05:06:35    -3130.239066        0.127599
    FIRE:   12 05:06:35    -3130.343764        0.107444
    FIRE:   13 05:06:35    -3130.408191        0.105237
    FIRE:   14 05:06:35    -3130.430606        0.081156
    FIRE:   15 05:06:35    -3130.444893        0.076541
    FIRE:   16 05:06:36    -3130.469668        0.068947
    FIRE:   17 05:06:36    -3130.497958        0.059786
    FIRE:   18 05:06:36    -3130.526015        0.048124
    FIRE:   19 05:06:36    -3130.550095        0.042873
    FIRE:   20 05:06:36    -3130.570362        0.050978
    FIRE:   21 05:06:36    -3130.586748        0.050684
    FIRE:   22 05:06:36    -3130.599216        0.040290
    FIRE:   23 05:06:36    -3130.606180        0.026648
    FIRE:   24 05:06:36    -3130.608459        0.028279
    FIRE:   25 05:06:37    -3130.608876        0.027669
    FIRE:   26 05:06:37    -3130.610401        0.026463
    FIRE:   27 05:06:37    -3130.612210        0.024709
    FIRE:   28 05:06:37    -3130.614231        0.022467
    FIRE:   29 05:06:37    -3130.616692        0.019835
    FIRE:   30 05:06:37    -3130.618739        0.018257
    FIRE:   31 05:06:37    -3130.621115        0.019769
    FIRE:   32 05:06:37    -3130.622977        0.020521
    FIRE:   33 05:06:37    -3130.625344        0.020089
    FIRE:   34 05:06:37    -3130.627259        0.018179
    FIRE:   35 05:06:37    -3130.629322        0.014800
    FIRE:   36 05:06:38    -3130.630452        0.010544
    FIRE:   37 05:06:38    -3130.631953        0.010931
    FIRE:   38 05:06:38    -3130.633008        0.011222
    FIRE:   39 05:06:38    -3130.634638        0.011352
    FIRE:   40 05:06:38    -3130.636388        0.011728
    FIRE:   41 05:06:38    -3130.638124        0.012170
    FIRE:   42 05:06:38    -3130.639831        0.012997
    FIRE:   43 05:06:38    -3130.640126        0.007894
    FIRE:   44 05:06:38    -3130.640162        0.006968
    FIRE:   45 05:06:39    -3130.640492        0.006198
    FIRE:   46 05:06:39    -3130.640475        0.006515
    FIRE:   47 05:06:39    -3130.640574        0.006386
    FIRE:   48 05:06:39    -3130.640730        0.005500
    FIRE:   49 05:06:39    -3130.640497        0.003782
    FIRE:   50 05:06:39    -3130.640700        0.002904
    FIRE:   51 05:06:39    -3130.640720        0.002245
    FIRE:   52 05:06:39    -3130.640654        0.002160
    FIRE:   53 05:06:39    -3130.640812        0.002023
    FIRE:   54 05:06:39    -3130.640926        0.001795
    FIRE:   55 05:06:40    -3130.640969        0.001520
    FIRE:   56 05:06:40    -3130.641039        0.001206
    FIRE:   57 05:06:40    -3130.640825        0.001088
    FIRE:   58 05:06:40    -3130.640900        0.000958
    Optimizing Pd405Pt306
          Step     Time          Energy          fmax
    FIRE:    0 05:06:40    -2936.142464        0.854851
    FIRE:    1 05:06:40    -2937.220114        0.826394
    FIRE:    2 05:06:40    -2938.999905        0.768947
    FIRE:    3 05:06:40    -2940.498574        0.705648
    FIRE:    4 05:06:40    -2941.772741        0.635854
    FIRE:    5 05:06:41    -2942.872919        0.560104
    FIRE:    6 05:06:41    -2943.835201        0.480506
    FIRE:    7 05:06:41    -2944.687657        0.399513
    FIRE:    8 05:06:41    -2945.450229        0.318803
    FIRE:    9 05:06:41    -2946.133509        0.280909
    FIRE:   10 05:06:41    -2946.739569        0.241686
    FIRE:   11 05:06:41    -2947.267283        0.196206
    FIRE:   12 05:06:41    -2947.715539        0.165348
    FIRE:   13 05:06:42    -2948.083806        0.153118
    FIRE:   14 05:06:42    -2948.375531        0.137163
    FIRE:   15 05:06:42    -2948.598662        0.121425
    FIRE:   16 05:06:42    -2948.762756        0.124072
    FIRE:   17 05:06:42    -2948.876736        0.126385
    FIRE:   18 05:06:42    -2948.951986        0.112937
    FIRE:   19 05:06:42    -2948.998015        0.109445
    FIRE:   20 05:06:42    -2949.021608        0.111756
    FIRE:   21 05:06:42    -2949.079012        0.103268
    FIRE:   22 05:06:43    -2949.161898        0.082761
    FIRE:   23 05:06:43    -2949.237292        0.082829
    FIRE:   24 05:06:43    -2949.295640        0.089292
    FIRE:   25 05:06:43    -2949.334601        0.062579
    FIRE:   26 05:06:43    -2949.351760        0.056082
    FIRE:   27 05:06:43    -2949.353941        0.052316
    FIRE:   28 05:06:43    -2949.358608        0.045132
    FIRE:   29 05:06:43    -2949.364378        0.035212
    FIRE:   30 05:06:43    -2949.370473        0.030386
    FIRE:   31 05:06:43    -2949.376524        0.028407
    FIRE:   32 05:06:43    -2949.382239        0.026030
    FIRE:   33 05:06:43    -2949.386982        0.023332
    FIRE:   34 05:06:44    -2949.392044        0.020129
    FIRE:   35 05:06:44    -2949.396432        0.020099
    FIRE:   36 05:06:44    -2949.400091        0.021226
    FIRE:   37 05:06:44    -2949.404231        0.020055
    FIRE:   38 05:06:44    -2949.409252        0.018011
    FIRE:   39 05:06:44    -2949.413977        0.014953
    FIRE:   40 05:06:44    -2949.418989        0.015857
    FIRE:   41 05:06:44    -2949.423180        0.019960
    FIRE:   42 05:06:44    -2949.424781        0.013920
    FIRE:   43 05:06:44    -2949.424999        0.012877
    FIRE:   44 05:06:44    -2949.425377        0.010904
    FIRE:   45 05:06:45    -2949.425272        0.008241
    FIRE:   46 05:06:45    -2949.425685        0.005244
    FIRE:   47 05:06:45    -2949.426045        0.005063
    FIRE:   48 05:06:45    -2949.425682        0.005181
    FIRE:   49 05:06:45    -2949.426244        0.005148
    FIRE:   50 05:06:45    -2949.425904        0.004175
    FIRE:   51 05:06:45    -2949.426189        0.004060
    FIRE:   52 05:06:45    -2949.426053        0.003837
    FIRE:   53 05:06:45    -2949.425944        0.003512
    FIRE:   54 05:06:46    -2949.425973        0.003111
    FIRE:   55 05:06:46    -2949.425993        0.002623
    FIRE:   56 05:06:46    -2949.426294        0.002099
    FIRE:   57 05:06:46    -2949.426012        0.001557
    FIRE:   58 05:06:46    -2949.426477        0.001123
    FIRE:   59 05:06:46    -2949.426254        0.001141
    FIRE:   60 05:06:46    -2949.426264        0.001393
    FIRE:   61 05:06:46    -2949.426410        0.001671
    FIRE:   62 05:06:46    -2949.426095        0.001751
    FIRE:   63 05:06:46    -2949.426225        0.001728
    FIRE:   64 05:06:46    -2949.426248        0.001679
    FIRE:   65 05:06:47    -2949.426392        0.001608
    FIRE:   66 05:06:47    -2949.426136        0.001516
    FIRE:   67 05:06:47    -2949.426391        0.001404
    FIRE:   68 05:06:47    -2949.426468        0.001275
    FIRE:   69 05:06:47    -2949.426149        0.001135
    FIRE:   70 05:06:47    -2949.426182        0.000965
    Optimizing Pd711
          Step     Time          Energy          fmax
    FIRE:    0 05:06:47    -2422.823697        0.360904
    FIRE:    1 05:06:47    -2423.073756        0.333712
    FIRE:    2 05:06:47    -2423.507992        0.304226
    FIRE:    3 05:06:47    -2424.031505        0.272105
    FIRE:    4 05:06:48    -2424.555334        0.229604
    FIRE:    5 05:06:48    -2425.033977        0.189223
    FIRE:    6 05:06:48    -2425.436250        0.208382
    FIRE:    7 05:06:48    -2425.784118        0.192289
    FIRE:    8 05:06:48    -2426.085797        0.151484
    FIRE:    9 05:06:48    -2426.335914        0.115907
    FIRE:   10 05:06:48    -2426.531106        0.128541
    FIRE:   11 05:06:48    -2426.676634        0.120186
    FIRE:   12 05:06:48    -2426.782396        0.104547
    FIRE:   13 05:06:49    -2426.854633        0.113271
    FIRE:   14 05:06:49    -2426.894728        0.094800
    FIRE:   15 05:06:49    -2426.904727        0.084601
    FIRE:   16 05:06:49    -2426.917991        0.075479
    FIRE:   17 05:06:49    -2426.940463        0.063126
    FIRE:   18 05:06:49    -2426.967530        0.059500
    FIRE:   19 05:06:49    -2426.992471        0.053825
    FIRE:   20 05:06:49    -2427.013807        0.047062
    FIRE:   21 05:06:49    -2427.030480        0.049089
    FIRE:   22 05:06:50    -2427.044643        0.046383
    FIRE:   23 05:06:50    -2427.056723        0.045267
    FIRE:   24 05:06:50    -2427.065122        0.029118
    FIRE:   25 05:06:50    -2427.068324        0.029192
    FIRE:   26 05:06:50    -2427.068814        0.027890
    FIRE:   27 05:06:50    -2427.069915        0.025382
    FIRE:   28 05:06:50    -2427.070843        0.021859
    FIRE:   29 05:06:50    -2427.072794        0.017589
    FIRE:   30 05:06:50    -2427.073930        0.015989
    FIRE:   31 05:06:51    -2427.075962        0.014601
    FIRE:   32 05:06:51    -2427.077286        0.013146
    FIRE:   33 05:06:51    -2427.078323        0.012604
    FIRE:   34 05:06:51    -2427.079718        0.012204
    FIRE:   35 05:06:51    -2427.081257        0.009762
    FIRE:   36 05:06:51    -2427.082506        0.010402
    FIRE:   37 05:06:51    -2427.082643        0.010911
    FIRE:   38 05:06:51    -2427.083651        0.008835
    FIRE:   39 05:06:51    -2427.083637        0.009758
    FIRE:   40 05:06:52    -2427.084294        0.009453
    FIRE:   41 05:06:52    -2427.085559        0.008523
    FIRE:   42 05:06:52    -2427.085955        0.010269
    FIRE:   43 05:06:52    -2427.087267        0.008518
    FIRE:   44 05:06:52    -2427.087216        0.005742
    FIRE:   45 05:06:52    -2427.087116        0.004991
    FIRE:   46 05:06:52    -2427.087314        0.003645
    FIRE:   47 05:06:52    -2427.087510        0.003000
    FIRE:   48 05:06:52    -2427.087539        0.002685
    FIRE:   49 05:06:52    -2427.087881        0.003552
    FIRE:   50 05:06:53    -2427.087747        0.003360
    FIRE:   51 05:06:53    -2427.087616        0.002194
    FIRE:   52 05:06:53    -2427.088008        0.001295
    FIRE:   53 05:06:53    -2427.088111        0.002169
    FIRE:   54 05:06:53    -2427.088089        0.002303
    FIRE:   55 05:06:53    -2427.087774        0.001213
    FIRE:   56 05:06:53    -2427.088140        0.001116
    FIRE:   57 05:06:53    -2427.088103        0.001009
    FIRE:   58 05:06:53    -2427.088123        0.000960


それぞれのエネルギーが得られたので、ここからexcess energyを計算することができます。


```python
E_exess_pd306pt405 = (E_pd306pt405 - 306 / 711 * E_pd711 - 405 / 711 * E_pt711) / 711
E_exess_pd405pt306 = (E_pd405pt306 - 405 / 711 * E_pd711 - 306 / 711 * E_pt711) / 711

print(f"Excess energy of Pd306Pt405 = {E_exess_pd306pt405 * 1000:.2f} meV")
print(f"Excess energy of Pd405Pt306 = {E_exess_pd405pt306 * 1000:.2f} meV")
```

    Excess energy of Pd306Pt405 = -30.14 meV
    Excess energy of Pd405Pt306 = -9.78 meV


Excess energyが負となっているPd306Pt405は単独でPd, Pt nanoparticleが存在するよりも混合したほうが安定であるということを示唆しています。&lt;br/&gt;
Pd405Pt306の構造はPd306Pt405と比べると不安定であることになります。

この結果から、**PdとPtの合金微粒子では、Pdが外殻に来るような構造が自然界では安定して存在しそう**だという考察ができます。

最後に、Pt306個、Pd405個がランダムに配置されたような構造のexcess energyを求めてみましょう。


```python
Pd405Pt306_random = Pt711.copy()
# randomly choose 405 atoms to be replaced to Pd
pd_indices = np.random.choice(np.arange(711), 405, replace=False)
Pd405Pt306_random.numbers[pd_indices] = atomic_numbers["Pd"]

print(f"Replaced to Pd with indices {pd_indices}")
view_ngl(Pd405Pt306_random)
```

    Replaced to Pd with indices [557 262 313 188 247  15 220 678 208  53 494 278 384 301 350 168 619 425
      60 245 123 603  45 594 346 638 516  89 404 559   9 334 182 476 394 525
     622 501 221 503 118 106 120 642 698 171 347 423 397 355   0 561  39 543
     196  82 377 645 170 183 222 625 595 502 169 291 211   2 156 530 600 545
     518 692  44 686 424   8 510 623 632 379 212 386 409 531  43  57 635 442
     263 290 460 296 214 505 242 679  75 639  36 145 136 246 283 259 337 508
     288 544  71 491 612 321 429 703 343  55 511 272  64 710 658 327 422 200
     624 489 708 497 426 583 352  63 656 467 325 226 537 215 174 473   4 433
     646 563 108 615 668 388 199 688 702 270 341 382 587 536 441 495 474 651
     150 435 618 370 565 428 192  96 661 165 611 269 655  12 677  34 202 461
     306 224  90 427 414 223 459  94 228   3 693 504 683 593 398 597 447 408
     187 637  97 665 339 621  67 672 237 328 332 164 177 193  59 191 311 111
     670  19 469 641 419 134 371 430 402  74  66 573 261 329 130 614 574 244
      47 472 119 393 493 584 268 417 450 631 648  99 207 216 455 176 304 418
     378  68 517 581 629 556 550  17 186 604 107 345 360 519 650 138  73 310
     636  27 294 206 167  28 488 401 499 369 526 514 218 232   7 500 528  62
     308 303 575 297 571 273 552 131 687 385  81 392 353 161  69  84  54 383
     309  24 399 436 527 569 240 243  10 126 373 653 236 257 660  13 295 101
     451 153 674  65 707  33 694 406 669 319 634 659 478 372 286   1 258 152
     358 453 160 521 381 415  40  52  32 652 689 287 452 463 116 482 485 127
     567 676 178 264 496 522 367 155 564 159 413 317 671 180 354 146 443 462
     103 255  76  87 331 326 533 282 133 534 602 490 626 241 172  78 175 579
     209 590 681 541 132 194  70 157  48]





    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pd', 'Pt'), value='Al…




```python
print("Optimizing Pd405Pt306_random")
Pd405Pt306_random.calc = calculator
E_pd405pt306_random = get_opt_energy(Pd405Pt306_random)
```

    Optimizing Pd405Pt306_random
          Step     Time          Energy          fmax
    FIRE:    0 05:06:54    -2951.463955        0.686472
    FIRE:    1 05:06:54    -2951.883058        0.655774
    FIRE:    2 05:06:54    -2952.612842        0.597436
    FIRE:    3 05:06:54    -2953.475407        0.512045
    FIRE:    4 05:06:54    -2954.220137        0.405601
    FIRE:    5 05:06:54    -2954.781796        0.287330
    FIRE:    6 05:06:54    -2955.206977        0.279969
    FIRE:    7 05:06:54    -2955.536052        0.267057
    FIRE:    8 05:06:54    -2955.794894        0.215046
    FIRE:    9 05:06:54    -2955.996466        0.194621
    FIRE:   10 05:06:54    -2956.149418        0.202107
    FIRE:   11 05:06:55    -2956.267670        0.195909
    FIRE:   12 05:06:55    -2956.365972        0.174524
    FIRE:   13 05:06:55    -2956.449970        0.176056
    FIRE:   14 05:06:55    -2956.515675        0.163480
    FIRE:   15 05:06:55    -2956.559975        0.171774
    FIRE:   16 05:06:55    -2956.589479        0.164634
    FIRE:   17 05:06:55    -2956.639216        0.150769
    FIRE:   18 05:06:55    -2956.693862        0.130952
    FIRE:   19 05:06:55    -2956.742611        0.106435
    FIRE:   20 05:06:55    -2956.779987        0.079620
    FIRE:   21 05:06:55    -2956.806685        0.056871
    FIRE:   22 05:06:56    -2956.823439        0.051401
    FIRE:   23 05:06:56    -2956.833863        0.057812
    FIRE:   24 05:06:56    -2956.836509        0.065774
    FIRE:   25 05:06:56    -2956.838056        0.063241
    FIRE:   26 05:06:56    -2956.840871        0.058339
    FIRE:   27 05:06:56    -2956.844654        0.051365
    FIRE:   28 05:06:56    -2956.849296        0.042768
    FIRE:   29 05:06:56    -2956.853617        0.035873
    FIRE:   30 05:06:56    -2956.857245        0.033567
    FIRE:   31 05:06:56    -2956.860880        0.031002
    FIRE:   32 05:06:57    -2956.864675        0.027350
    FIRE:   33 05:06:57    -2956.867384        0.022381
    FIRE:   34 05:06:57    -2956.870122        0.021284
    FIRE:   35 05:06:57    -2956.872151        0.018652
    FIRE:   36 05:06:57    -2956.874154        0.020530
    FIRE:   37 05:06:57    -2956.875844        0.022471
    FIRE:   38 05:06:57    -2956.877792        0.021527
    FIRE:   39 05:06:57    -2956.880602        0.023119
    FIRE:   40 05:06:57    -2956.883679        0.022468
    FIRE:   41 05:06:57    -2956.887104        0.019191
    FIRE:   42 05:06:57    -2956.890437        0.013375
    FIRE:   43 05:06:58    -2956.892346        0.012415
    FIRE:   44 05:06:58    -2956.892579        0.011058
    FIRE:   45 05:06:58    -2956.892868        0.009238
    FIRE:   46 05:06:58    -2956.893189        0.007834
    FIRE:   47 05:06:58    -2956.893483        0.006969
    FIRE:   48 05:06:58    -2956.893334        0.007206
    FIRE:   49 05:06:58    -2956.893595        0.007273
    FIRE:   50 05:06:58    -2956.893736        0.006222
    FIRE:   51 05:06:58    -2956.894191        0.005410
    FIRE:   52 05:06:58    -2956.893971        0.006673
    FIRE:   53 05:06:59    -2956.894060        0.006728
    FIRE:   54 05:06:59    -2956.894351        0.005678
    FIRE:   55 05:06:59    -2956.894622        0.004578
    FIRE:   56 05:06:59    -2956.894916        0.006321
    FIRE:   57 05:06:59    -2956.895574        0.004792
    FIRE:   58 05:06:59    -2956.895868        0.005141
    FIRE:   59 05:06:59    -2956.896080        0.004266
    FIRE:   60 05:06:59    -2956.896043        0.004952
    FIRE:   61 05:06:59    -2956.896294        0.004693
    FIRE:   62 05:06:59    -2956.896478        0.004576
    FIRE:   63 05:07:00    -2956.896783        0.002752
    FIRE:   64 05:07:00    -2956.896877        0.002213
    FIRE:   65 05:07:00    -2956.896871        0.003885
    FIRE:   66 05:07:00    -2956.896523        0.003477
    FIRE:   67 05:07:00    -2956.896561        0.002706
    FIRE:   68 05:07:00    -2956.896907        0.001703
    FIRE:   69 05:07:00    -2956.896514        0.001313
    FIRE:   70 05:07:00    -2956.896905        0.001457
    FIRE:   71 05:07:00    -2956.896486        0.001937
    FIRE:   72 05:07:01    -2956.896520        0.002071
    FIRE:   73 05:07:01    -2956.896531        0.001899
    FIRE:   74 05:07:01    -2956.896908        0.001363
    FIRE:   75 05:07:01    -2956.896635        0.000957



```python
view_ngl(Pd405Pt306_random)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pd', 'Pt'), value='Al…




```python
E_exess_pd405pt306_random = (E_pd405pt306_random - 405 / 711 * E_pd711 - 306 / 711 * E_pt711) / 711

print(f"Excess energy of Pd306Pt405 random = {E_exess_pd405pt306_random * 1000:.2f} meV")
```

    Excess energy of Pd306Pt405 random = -20.29 meV


この場合も負のexcess energyが得られ、この構造が安定になりうるという結果が得られました。

今回は完全にランダムな配置を行って構造を作成しましたが、実際は物質が混ざる場合には局所的にパターンを取ることでより安定な構造となる可能性もあります。
こういった分析を進める場合はモンテカルロ法を用いて配置パターンを最適化するということがよく行われますがここでは省略します。&lt;br/&gt;
（モンテカルロ法の説明及びその実例については将来後述されるかもしれません。）

このように、様々な構造を作成し、そのexcess energyを計算することで合金の安定構造がどのようになるかを解析する事ができます。

## 担体上への微粒子吸着構造作成

本チュートリアルでは、担体上に微粒子をおいたモデリングの簡単な計算事例を紹介するにとどめ、この構造の現実での妥当性検証までは行いません。

ここでは、担持貴金属担体としてSnO2を、Nano particleとしてはPt を用いて、担体上へPt Nano particleを乗せた構造を作成してみます。

 - [The effect of SnO2(110) supports on the geometrical and electronic properties of platinum nanoparticles](https://link.springer.com/article/10.1007/s42452-019-1478-0)
 - [Calculations of Real-System Nanoparticles Using Universal Neural Network Potential PFP](https://arxiv.org/abs/2107.00963)


```python
# bulk, support material
from ase.io import read
SnO2 = read("../input/SnO2_mp-856_conventional_standard.cif")
view_ngl(SnO2, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'Sn'), value='All…



ここではSlab構造の作成(Bulkからの特定ミラー面の切り出し)にpymatgenの `SlabGenerator` を用いています。&lt;br/&gt;
これは`SlabGenerator` のinstance化をする際の初期化パラメータとして作りたいSlab構造のパラメータを指定し、&lt;br/&gt;
`get_slabs` 関数を呼ぶことで指定された初期化パラメータにマッチするSlab構造が列挙されて得られます。

 - [pymatgen document](https://pymatgen.org/pymatgen.core.surface.html#pymatgen.core.surface.SlabGenerator)

詳しくは、以下のMaterials Project Workshopでも紹介されています。

 - [Working with Surfaces and Interfaces - The Materials Project Workshop](https://workshop.materialsproject.org/lessons/03_heterointerfaces/Main%20Lesson/)


```python
# Generate slab
from ase.build import make_supercell
import numpy as np
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor


slab_gen = SlabGenerator(
    initial_structure=AseAtomsAdaptor.get_structure(SnO2),
    miller_index=[1,1,0],
    min_slab_size=5.0, # ここで層の出方が変わる
    min_vacuum_size=15.0, # nanoparticleを載せるので大きめに
    lll_reduce=False,
    center_slab=True,
    primitive=True,
    max_normal_search=1,
)
slabs = slab_gen.get_slabs(tol=0.3, bonds=None, max_broken_bonds=0, symmetrize=False)
```

得られた `slabs` はpymatgenのinstanceです。これは、`AseAtomsAdaptor` を用いてASEの`Atoms` instanceに変換することができます。

得られた構造をASE `Atoms`に変換して、可視化してみます。


```python
slab_atoms_list = [AseAtomsAdaptor.get_atoms(slab) for slab in slabs]
view_ngl(slab_atoms_list, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(max_frame=1), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'Sn'),…



２つの構造が得られました。

以降では2つ目の構造を使用して、Slabを作成していきます。

ASEでsupercellを作成した後、z軸方向に平行移動して一番下にある原子(Layer)の高さが0になるよう設定しています。


```python
slab = slab_atoms_list[1].copy()

# make supercell: expand to xy-plane
#slab = make_supercell(slab, [[4, 0, 0], [0, 2, 0], [0, 0, 1]])  # this is same with below.
slab = slab * (4, 2, 1)

# shift `slab` to bottom of cell
min_pos_z = np.min(slab.positions, axis=0)[2]
slab.set_positions(slab.positions - [0, 0, min_pos_z])
```


```python
view_ngl(slab, representations=["ball+stick"], replace_structure=True)
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'Sn'), value='All…



担体のSlab構造を作成することができました。

次にこの担体上に乗せるNano particle clusterを作成します。


```python
Pt55 = Octahedron("Pt", 5, cutoff=2)

# cut cluster to make half-sphere
cluster = Pt55.copy()
# Rotate cluster to make triangle surface comes to top.
cluster.rotate([0, 0, 1], [1, 1, 1], center=cluster.get_center_of_mass())
# Cut bottom 2 layers
for _ in range(2):
    target = np.round(cluster.positions, decimals=0)
    del cluster[np.where(target[:,2]==np.min(target[:,2]))[0]]
```


```python
view_ngl(cluster, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'Pt'), value='All'), D…



Nano particle clusterを担体の上に乗せましょう。

ここでは、プログラムを書いて担体真ん中の上にNano particleを乗せていますが、
5章で紹介する`SurfaceEditor`を用いると、InteractiveにNano particleを移動させて初期構造を作成することも可能です。


```python
from ase.data import atomic_numbers, chemical_symbols, covalent_radii


# Put Pt37 on top of SnO2. Both perpendincular lines goes z-axis.
slab_xy_size = np.min(slab.cell.cellpar()[0:2])
cluster_xy_size = np.max(
    (np.max(cluster.positions, axis=0) - np.min(cluster.positions, axis=0))[0:2]
)
# Vacuum size
min_slab_xy_size = cluster_xy_size + 15
for i in range(1, 5):
    print(f"i={i}")
    if slab_xy_size * i &lt; min_slab_xy_size:
        pass
    else:
        slab_sc = make_supercell(slab, [[i, 0, 0], [0, i, 0], [0, 0, 1]])
        break
slab_surface_xy_center = np.append(
    slab_sc.cell.cellpar()[0:2] / 2, np.max(slab_sc.positions, axis=0)[2]
)
cluster_surface_xy_center = np.append(
    np.mean(cluster.positions, axis=0)[0:2], np.min(cluster.positions, axis=0)[2]
)
cluster = Atoms(cluster.get_chemical_symbols(), cluster.positions - cluster_surface_xy_center)
slab_surface_covalent_radii = covalent_radii[
    slab_sc.get_atomic_numbers()[np.argmax(slab_sc.positions, axis=0)[2]]
]
cluster_surface_covalent_radii = covalent_radii[
    cluster.get_atomic_numbers()[np.argmin(cluster.positions, axis=0)[2]]
]
cluster.translate(
    slab_surface_xy_center + [0, 0, slab_surface_covalent_radii + cluster_surface_covalent_radii]
)
supported = slab_sc.copy()
supported += cluster
```

    i=1
    i=2



```python
view_ngl(supported, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'Sn', 'Pt'), valu…




```python
import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from ase.optimize import BFGS, LBFGS, FIRE
from ase.io import Trajectory
print(f"pfp_api_client: {pfp_api_client.__version__}")

estimator = Estimator(calc_mode=EstimatorCalcMode.PBE_U, model_version="v8.0.0")
calculator = ASECalculator(estimator)
supported.calc = calculator
traj = Trajectory("./SnO2_Pt119.traj", "w", supported)
opt = FIRE(supported)
opt.attach(traj.write, interval=1)
opt.run(fmax=0.005)
traj.close()
```

    pfp_api_client: 1.23.1
          Step     Time          Energy          fmax
    FIRE:    0 05:07:02    -1856.170263        2.665578
    FIRE:    1 05:07:02    -1862.262675        2.242122
    FIRE:    2 05:07:02    -1867.449870        2.009662
    FIRE:    3 05:07:02    -1871.758847        1.937803
    FIRE:    4 05:07:02    -1875.235394        1.854351
    FIRE:    5 05:07:02    -1877.972361        1.756824
    FIRE:    6 05:07:02    -1880.108393        1.642316
    FIRE:    7 05:07:02    -1881.825760        1.508070
    FIRE:    8 05:07:02    -1883.341446        1.353713
    FIRE:    9 05:07:02    -1884.849230        1.303777
    FIRE:   10 05:07:02    -1886.443856        1.201917
    FIRE:   11 05:07:03    -1888.094835        0.987404
    FIRE:   12 05:07:03    -1889.686761        0.971872
    FIRE:   13 05:07:03    -1891.098209        0.956527
    FIRE:   14 05:07:03    -1892.274196        0.917137
    FIRE:   15 05:07:03    -1893.250302        0.917339
    FIRE:   16 05:07:03    -1894.129911        1.004214
    FIRE:   17 05:07:03    -1895.013647        1.092757
    FIRE:   18 05:07:03    -1895.940917        1.175548
    FIRE:   19 05:07:03    -1896.872013        1.233957
    FIRE:   20 05:07:03    -1897.722148        1.252479
    FIRE:   21 05:07:03    -1898.419551        1.221698
    FIRE:   22 05:07:03    -1898.950705        1.139072
    FIRE:   23 05:07:04    -1899.368507        1.014610
    FIRE:   24 05:07:04    -1899.740725        0.816489
    FIRE:   25 05:07:04    -1900.097960        0.519997
    FIRE:   26 05:07:04    -1900.418936        0.427894
    FIRE:   27 05:07:04    -1900.666274        0.580771
    FIRE:   28 05:07:04    -1900.841934        0.872977
    FIRE:   29 05:07:04    -1900.999616        0.969969
    FIRE:   30 05:07:04    -1901.182197        0.913152
    FIRE:   31 05:07:04    -1901.366238        0.619061
    FIRE:   32 05:07:04    -1901.500156        0.547311
    FIRE:   33 05:07:04    -1901.597105        0.638849
    FIRE:   34 05:07:04    -1901.703571        0.638583
    FIRE:   35 05:07:05    -1901.812396        0.534441
    FIRE:   36 05:07:05    -1901.872484        0.688772
    FIRE:   37 05:07:05    -1901.995474        0.906043
    FIRE:   38 05:07:05    -1902.037921        0.632184
    FIRE:   39 05:07:05    -1902.096061        0.343188
    FIRE:   40 05:07:05    -1902.118179        0.150365
    FIRE:   41 05:07:05    -1902.121692        0.127778
    FIRE:   42 05:07:05    -1902.127382        0.099333
    FIRE:   43 05:07:05    -1902.133649        0.099688
    FIRE:   44 05:07:05    -1902.139503        0.100180
    FIRE:   45 05:07:05    -1902.145153        0.100802
    FIRE:   46 05:07:05    -1902.151803        0.101519
    FIRE:   47 05:07:06    -1902.159854        0.102313
    FIRE:   48 05:07:06    -1902.170347        0.103250
    FIRE:   49 05:07:06    -1902.182006        0.104338
    FIRE:   50 05:07:06    -1902.193443        0.105569
    FIRE:   51 05:07:06    -1902.205230        0.106905
    FIRE:   52 05:07:06    -1902.218620        0.108270
    FIRE:   53 05:07:06    -1902.232278        0.109647
    FIRE:   54 05:07:06    -1902.245239        0.110953
    FIRE:   55 05:07:06    -1902.259210        0.111549
    FIRE:   56 05:07:06    -1902.274654        0.110601
    FIRE:   57 05:07:06    -1902.292077        0.108556
    FIRE:   58 05:07:06    -1902.313556        0.106564
    FIRE:   59 05:07:07    -1902.339022        0.102928
    FIRE:   60 05:07:07    -1902.369553        0.096015
    FIRE:   61 05:07:07    -1902.403737        0.087521
    FIRE:   62 05:07:07    -1902.441864        0.127309
    FIRE:   63 05:07:07    -1902.480308        0.262127
    FIRE:   64 05:07:07    -1902.489249        0.552484
    FIRE:   65 05:07:07    -1902.546538        0.125462
    FIRE:   66 05:07:07    -1902.548483        0.125572
    FIRE:   67 05:07:07    -1902.551481        0.125704
    FIRE:   68 05:07:07    -1902.554729        0.125682
    FIRE:   69 05:07:07    -1902.558526        0.125287
    FIRE:   70 05:07:07    -1902.563782        0.124341
    FIRE:   71 05:07:07    -1902.570344        0.122752
    FIRE:   72 05:07:07    -1902.577610        0.120500
    FIRE:   73 05:07:07    -1902.585920        0.117343
    FIRE:   74 05:07:08    -1902.596410        0.113322
    FIRE:   75 05:07:08    -1902.609888        0.108598
    FIRE:   76 05:07:08    -1902.626009        0.110980
    FIRE:   77 05:07:08    -1902.646351        0.125224
    FIRE:   78 05:07:08    -1902.672708        0.141231
    FIRE:   79 05:07:08    -1902.706629        0.158127
    FIRE:   80 05:07:08    -1902.752187        0.157374
    FIRE:   81 05:07:08    -1902.813123        0.168110
    FIRE:   82 05:07:08    -1902.893730        0.157536
    FIRE:   83 05:07:08    -1902.993544        0.183617
    FIRE:   84 05:07:08    -1903.089208        0.206657
    FIRE:   85 05:07:09    -1903.183625        0.207540
    FIRE:   86 05:07:09    -1903.276750        0.193371
    FIRE:   87 05:07:09    -1903.364973        0.191780
    FIRE:   88 05:07:09    -1903.429932        0.177726
    FIRE:   89 05:07:09    -1903.459781        0.237578
    FIRE:   90 05:07:09    -1903.471484        0.145258
    FIRE:   91 05:07:09    -1903.486317        0.109691
    FIRE:   92 05:07:09    -1903.496796        0.146517
    FIRE:   93 05:07:09    -1903.507670        0.181363
    FIRE:   94 05:07:09    -1903.515630        0.167505
    FIRE:   95 05:07:09    -1903.526061        0.137261
    FIRE:   96 05:07:09    -1903.537515        0.096263
    FIRE:   97 05:07:09    -1903.549348        0.075398
    FIRE:   98 05:07:10    -1903.562992        0.101925
    FIRE:   99 05:07:10    -1903.577643        0.153720
    FIRE:  100 05:07:10    -1903.587584        0.090648
    FIRE:  101 05:07:10    -1903.585780        0.117460
    FIRE:  102 05:07:10    -1903.588151        0.091638
    FIRE:  103 05:07:10    -1903.591322        0.088607
    FIRE:  104 05:07:10    -1903.593200        0.089134
    FIRE:  105 05:07:10    -1903.593538        0.089907
    FIRE:  106 05:07:10    -1903.593674        0.089956
    FIRE:  107 05:07:10    -1903.594118        0.090052
    FIRE:  108 05:07:10    -1903.594602        0.090202
    FIRE:  109 05:07:10    -1903.595166        0.090401
    FIRE:  110 05:07:10    -1903.595725        0.090647
    FIRE:  111 05:07:11    -1903.596297        0.090951
    FIRE:  112 05:07:11    -1903.596894        0.091304
    FIRE:  113 05:07:11    -1903.597512        0.091748
    FIRE:  114 05:07:11    -1903.598348        0.092302
    FIRE:  115 05:07:11    -1903.599377        0.092970
    FIRE:  116 05:07:11    -1903.600752        0.093768
    FIRE:  117 05:07:11    -1903.602551        0.094725
    FIRE:  118 05:07:11    -1903.604673        0.095912
    FIRE:  119 05:07:11    -1903.607214        0.098675
    FIRE:  120 05:07:11    -1903.610187        0.101455
    FIRE:  121 05:07:11    -1903.613999        0.103958
    FIRE:  122 05:07:11    -1903.618875        0.109311
    FIRE:  123 05:07:11    -1903.624811        0.115539
    FIRE:  124 05:07:11    -1903.632141        0.122548
    FIRE:  125 05:07:12    -1903.641540        0.130080
    FIRE:  126 05:07:12    -1903.653040        0.137014
    FIRE:  127 05:07:12    -1903.667451        0.147126
    FIRE:  128 05:07:12    -1903.685221        0.152102
    FIRE:  129 05:07:12    -1903.705488        0.121521
    FIRE:  130 05:07:12    -1903.726903        0.110721
    FIRE:  131 05:07:12    -1903.749299        0.122232
    FIRE:  132 05:07:12    -1903.772122        0.124668
    FIRE:  133 05:07:12    -1903.790009        0.157177
    FIRE:  134 05:07:12    -1903.798462        0.262995
    FIRE:  135 05:07:12    -1903.803230        0.228290
    FIRE:  136 05:07:12    -1903.807968        0.168494
    FIRE:  137 05:07:13    -1903.814290        0.101565
    FIRE:  138 05:07:13    -1903.820813        0.092466
    FIRE:  139 05:07:13    -1903.825523        0.085685
    FIRE:  140 05:07:13    -1903.831203        0.090986
    FIRE:  141 05:07:13    -1903.835542        0.092721
    FIRE:  142 05:07:13    -1903.840656        0.087556
    FIRE:  143 05:07:13    -1903.846138        0.076264
    FIRE:  144 05:07:13    -1903.852992        0.064538
    FIRE:  145 05:07:13    -1903.860426        0.065194
    FIRE:  146 05:07:13    -1903.870812        0.060389
    FIRE:  147 05:07:13    -1903.883328        0.057322
    FIRE:  148 05:07:13    -1903.897072        0.048027
    FIRE:  149 05:07:14    -1903.909321        0.112484
    FIRE:  150 05:07:14    -1903.881179        0.391966
    FIRE:  151 05:07:14    -1903.928729        0.092869
    FIRE:  152 05:07:14    -1903.930641        0.063874
    FIRE:  153 05:07:14    -1903.932544        0.050678
    FIRE:  154 05:07:14    -1903.932889        0.048923
    FIRE:  155 05:07:14    -1903.932974        0.048788
    FIRE:  156 05:07:14    -1903.933373        0.048527
    FIRE:  157 05:07:14    -1903.933494        0.048128
    FIRE:  158 05:07:14    -1903.933826        0.047604
    FIRE:  159 05:07:14    -1903.934129        0.046952
    FIRE:  160 05:07:14    -1903.934414        0.046169
    FIRE:  161 05:07:15    -1903.934720        0.045284
    FIRE:  162 05:07:15    -1903.935122        0.044551
    FIRE:  163 05:07:15    -1903.935593        0.044376
    FIRE:  164 05:07:15    -1903.936319        0.044188
    FIRE:  165 05:07:15    -1903.937186        0.043992
    FIRE:  166 05:07:15    -1903.938165        0.043787
    FIRE:  167 05:07:15    -1903.939283        0.043597
    FIRE:  168 05:07:15    -1903.940608        0.043466
    FIRE:  169 05:07:15    -1903.942309        0.043438
    FIRE:  170 05:07:15    -1903.944261        0.043602
    FIRE:  171 05:07:15    -1903.946513        0.044094
    FIRE:  172 05:07:15    -1903.949283        0.045032
    FIRE:  173 05:07:15    -1903.952596        0.046501
    FIRE:  174 05:07:16    -1903.956534        0.048471
    FIRE:  175 05:07:16    -1903.961206        0.050850
    FIRE:  176 05:07:16    -1903.966863        0.053290
    FIRE:  177 05:07:16    -1903.973550        0.054873
    FIRE:  178 05:07:16    -1903.981451        0.054489
    FIRE:  179 05:07:16    -1903.990636        0.052538
    FIRE:  180 05:07:16    -1904.000996        0.047763
    FIRE:  181 05:07:16    -1904.011242        0.058041
    FIRE:  182 05:07:16    -1904.011203        0.171578
    FIRE:  183 05:07:16    -1904.025070        0.038943
    FIRE:  184 05:07:17    -1904.025159        0.038934
    FIRE:  185 05:07:17    -1904.025348        0.038920
    FIRE:  186 05:07:17    -1904.025621        0.038909
    FIRE:  187 05:07:17    -1904.025867        0.038903
    FIRE:  188 05:07:17    -1904.026224        0.038902
    FIRE:  189 05:07:17    -1904.026743        0.038919
    FIRE:  190 05:07:17    -1904.027253        0.038945
    FIRE:  191 05:07:17    -1904.027869        0.038990
    FIRE:  192 05:07:17    -1904.028603        0.039035
    FIRE:  193 05:07:17    -1904.029563        0.039062
    FIRE:  194 05:07:17    -1904.030687        0.039018
    FIRE:  195 05:07:18    -1904.032017        0.038839
    FIRE:  196 05:07:18    -1904.033650        0.038529
    FIRE:  197 05:07:18    -1904.035482        0.038197
    FIRE:  198 05:07:18    -1904.037767        0.038016
    FIRE:  199 05:07:18    -1904.040435        0.037982
    FIRE:  200 05:07:18    -1904.043506        0.037759
    FIRE:  201 05:07:18    -1904.046900        0.036600
    FIRE:  202 05:07:18    -1904.050773        0.033697
    FIRE:  203 05:07:18    -1904.054530        0.034787
    FIRE:  204 05:07:18    -1904.055688        0.084179
    FIRE:  205 05:07:18    -1904.059434        0.023902
    FIRE:  206 05:07:18    -1904.056730        0.075488
    FIRE:  207 05:07:19    -1904.057997        0.057194
    FIRE:  208 05:07:19    -1904.059386        0.025089
    FIRE:  209 05:07:19    -1904.059717        0.024690
    FIRE:  210 05:07:19    -1904.059742        0.024685
    FIRE:  211 05:07:19    -1904.059747        0.024661
    FIRE:  212 05:07:19    -1904.059755        0.024641
    FIRE:  213 05:07:19    -1904.059819        0.024604
    FIRE:  214 05:07:19    -1904.059849        0.024551
    FIRE:  215 05:07:19    -1904.059908        0.024492
    FIRE:  216 05:07:19    -1904.059930        0.024423
    FIRE:  217 05:07:19    -1904.059958        0.024352
    FIRE:  218 05:07:20    -1904.060037        0.024244
    FIRE:  219 05:07:20    -1904.060086        0.024127
    FIRE:  220 05:07:20    -1904.060173        0.023984
    FIRE:  221 05:07:20    -1904.060316        0.023822
    FIRE:  222 05:07:20    -1904.060444        0.023628
    FIRE:  223 05:07:20    -1904.060599        0.023374
    FIRE:  224 05:07:20    -1904.060777        0.023082
    FIRE:  225 05:07:20    -1904.061043        0.022717
    FIRE:  226 05:07:20    -1904.061357        0.022261
    FIRE:  227 05:07:20    -1904.061634        0.021675
    FIRE:  228 05:07:20    -1904.062041        0.020952
    FIRE:  229 05:07:21    -1904.062556        0.020086
    FIRE:  230 05:07:21    -1904.063233        0.019073
    FIRE:  231 05:07:21    -1904.063842        0.017974
    FIRE:  232 05:07:21    -1904.064615        0.016796
    FIRE:  233 05:07:21    -1904.065550        0.015527
    FIRE:  234 05:07:21    -1904.066693        0.014106
    FIRE:  235 05:07:21    -1904.067903        0.012550
    FIRE:  236 05:07:21    -1904.069316        0.012666
    FIRE:  237 05:07:21    -1904.070626        0.025486
    FIRE:  238 05:07:21    -1904.070633        0.062790
    FIRE:  239 05:07:21    -1904.072690        0.014748
    FIRE:  240 05:07:21    -1904.071188        0.056420
    FIRE:  241 05:07:22    -1904.071908        0.042272
    FIRE:  242 05:07:22    -1904.072663        0.017555
    FIRE:  243 05:07:22    -1904.072811        0.015096
    FIRE:  244 05:07:22    -1904.072838        0.015095
    FIRE:  245 05:07:22    -1904.072849        0.015106
    FIRE:  246 05:07:22    -1904.072867        0.015116
    FIRE:  247 05:07:22    -1904.072878        0.015130
    FIRE:  248 05:07:22    -1904.072883        0.015143
    FIRE:  249 05:07:22    -1904.072996        0.015167
    FIRE:  250 05:07:22    -1904.072915        0.015180
    FIRE:  251 05:07:22    -1904.072930        0.015206
    FIRE:  252 05:07:22    -1904.072968        0.015248
    FIRE:  253 05:07:23    -1904.072996        0.015290
    FIRE:  254 05:07:23    -1904.073031        0.015339
    FIRE:  255 05:07:23    -1904.073074        0.015394
    FIRE:  256 05:07:23    -1904.073189        0.015479
    FIRE:  257 05:07:23    -1904.073261        0.015582
    FIRE:  258 05:07:23    -1904.073373        0.015719
    FIRE:  259 05:07:23    -1904.073477        0.015886
    FIRE:  260 05:07:23    -1904.073719        0.016108
    FIRE:  261 05:07:23    -1904.073873        0.016406
    FIRE:  262 05:07:23    -1904.074062        0.016785
    FIRE:  263 05:07:23    -1904.074378        0.017284
    FIRE:  264 05:07:24    -1904.074695        0.017950
    FIRE:  265 05:07:24    -1904.075152        0.018819
    FIRE:  266 05:07:24    -1904.075726        0.019866
    FIRE:  267 05:07:24    -1904.076383        0.021015
    FIRE:  268 05:07:24    -1904.077290        0.022206
    FIRE:  269 05:07:24    -1904.078397        0.023734
    FIRE:  270 05:07:24    -1904.079792        0.026171
    FIRE:  271 05:07:24    -1904.081548        0.029605
    FIRE:  272 05:07:24    -1904.082819        0.054769
    FIRE:  273 05:07:24    -1904.084437        0.034138
    FIRE:  274 05:07:24    -1904.083470        0.049763
    FIRE:  275 05:07:24    -1904.084048        0.036849
    FIRE:  276 05:07:25    -1904.084682        0.035197
    FIRE:  277 05:07:25    -1904.084846        0.035429
    FIRE:  278 05:07:25    -1904.084855        0.035454
    FIRE:  279 05:07:25    -1904.085065        0.035489
    FIRE:  280 05:07:25    -1904.084930        0.035546
    FIRE:  281 05:07:25    -1904.084994        0.035621
    FIRE:  282 05:07:25    -1904.085024        0.035726
    FIRE:  283 05:07:25    -1904.085067        0.035841
    FIRE:  284 05:07:25    -1904.085122        0.035985
    FIRE:  285 05:07:25    -1904.085162        0.036161
    FIRE:  286 05:07:25    -1904.085237        0.036380
    FIRE:  287 05:07:26    -1904.085349        0.036648
    FIRE:  288 05:07:26    -1904.085518        0.036976
    FIRE:  289 05:07:26    -1904.085695        0.037391
    FIRE:  290 05:07:26    -1904.085939        0.037911
    FIRE:  291 05:07:26    -1904.086200        0.038550
    FIRE:  292 05:07:26    -1904.086525        0.039360
    FIRE:  293 05:07:26    -1904.086988        0.040375
    FIRE:  294 05:07:26    -1904.087553        0.041666
    FIRE:  295 05:07:26    -1904.088224        0.043343
    FIRE:  296 05:07:26    -1904.089180        0.045567
    FIRE:  297 05:07:26    -1904.090306        0.048571
    FIRE:  298 05:07:26    -1904.091945        0.052665
    FIRE:  299 05:07:27    -1904.094099        0.058527
    FIRE:  300 05:07:27    -1904.097053        0.067188
    FIRE:  301 05:07:27    -1904.101287        0.077279
    FIRE:  302 05:07:27    -1904.107458        0.082007
    FIRE:  303 05:07:27    -1904.116520        0.092922
    FIRE:  304 05:07:27    -1904.129842        0.101447
    FIRE:  305 05:07:27    -1904.149065        0.102895
    FIRE:  306 05:07:27    -1904.172172        0.104890
    FIRE:  307 05:07:27    -1904.173001        0.238219
    FIRE:  308 05:07:27    -1904.201904        0.175324
    FIRE:  309 05:07:27    -1904.204026        0.172856
    FIRE:  310 05:07:27    -1904.206590        0.167927
    FIRE:  311 05:07:27    -1904.208311        0.160482
    FIRE:  312 05:07:28    -1904.210104        0.150994
    FIRE:  313 05:07:28    -1904.213342        0.140133
    FIRE:  314 05:07:28    -1904.217784        0.128385
    FIRE:  315 05:07:28    -1904.222028        0.116125
    FIRE:  316 05:07:28    -1904.226380        0.102499
    FIRE:  317 05:07:28    -1904.232366        0.101920
    FIRE:  318 05:07:28    -1904.240578        0.106826
    FIRE:  319 05:07:28    -1904.249613        0.111331
    FIRE:  320 05:07:28    -1904.261456        0.111254
    FIRE:  321 05:07:28    -1904.275803        0.104447
    FIRE:  322 05:07:28    -1904.291400        0.104734
    FIRE:  323 05:07:28    -1904.308956        0.104343
    FIRE:  324 05:07:28    -1904.327884        0.104614
    FIRE:  325 05:07:28    -1904.345722        0.101730
    FIRE:  326 05:07:29    -1904.359111        0.077625
    FIRE:  327 05:07:29    -1904.362205        0.125875
    FIRE:  328 05:07:29    -1904.370719        0.078520
    FIRE:  329 05:07:29    -1904.368399        0.092028
    FIRE:  330 05:07:29    -1904.370101        0.073619
    FIRE:  331 05:07:29    -1904.372282        0.061998
    FIRE:  332 05:07:29    -1904.373559        0.060091
    FIRE:  333 05:07:29    -1904.373474        0.057306
    FIRE:  334 05:07:29    -1904.373607        0.057166
    FIRE:  335 05:07:29    -1904.373784        0.056892
    FIRE:  336 05:07:29    -1904.373975        0.056481
    FIRE:  337 05:07:29    -1904.374284        0.055927
    FIRE:  338 05:07:29    -1904.374566        0.055222
    FIRE:  339 05:07:29    -1904.374821        0.054366
    FIRE:  340 05:07:30    -1904.375067        0.053347
    FIRE:  341 05:07:30    -1904.375250        0.052051
    FIRE:  342 05:07:30    -1904.375529        0.050431
    FIRE:  343 05:07:30    -1904.375846        0.048468
    FIRE:  344 05:07:30    -1904.376329        0.046129
    FIRE:  345 05:07:30    -1904.376922        0.043362
    FIRE:  346 05:07:30    -1904.377606        0.040096
    FIRE:  347 05:07:30    -1904.378337        0.036249
    FIRE:  348 05:07:30    -1904.379098        0.031770
    FIRE:  349 05:07:30    -1904.379947        0.028296
    FIRE:  350 05:07:30    -1904.381009        0.027023
    FIRE:  351 05:07:30    -1904.382114        0.027127
    FIRE:  352 05:07:30    -1904.383180        0.027873
    FIRE:  353 05:07:30    -1904.384422        0.028466
    FIRE:  354 05:07:31    -1904.385630        0.028484
    FIRE:  355 05:07:31    -1904.386870        0.027816
    FIRE:  356 05:07:31    -1904.388214        0.028658
    FIRE:  357 05:07:31    -1904.389533        0.027402
    FIRE:  358 05:07:31    -1904.390914        0.024117
    FIRE:  359 05:07:31    -1904.392434        0.019369
    FIRE:  360 05:07:31    -1904.393869        0.017925
    FIRE:  361 05:07:31    -1904.395048        0.018804
    FIRE:  362 05:07:31    -1904.395632        0.041300
    FIRE:  363 05:07:31    -1904.396482        0.016707
    FIRE:  364 05:07:31    -1904.396164        0.031320
    FIRE:  365 05:07:31    -1904.396354        0.024830
    FIRE:  366 05:07:31    -1904.396576        0.015147
    FIRE:  367 05:07:32    -1904.396688        0.014801
    FIRE:  368 05:07:32    -1904.396682        0.014296
    FIRE:  369 05:07:32    -1904.396676        0.014270
    FIRE:  370 05:07:32    -1904.396708        0.014222
    FIRE:  371 05:07:32    -1904.396713        0.014159
    FIRE:  372 05:07:32    -1904.396769        0.014071
    FIRE:  373 05:07:32    -1904.396770        0.013951
    FIRE:  374 05:07:32    -1904.396793        0.013825
    FIRE:  375 05:07:32    -1904.396831        0.013806
    FIRE:  376 05:07:32    -1904.396832        0.013775
    FIRE:  377 05:07:32    -1904.396845        0.013743
    FIRE:  378 05:07:33    -1904.396896        0.013678
    FIRE:  379 05:07:33    -1904.396941        0.013605
    FIRE:  380 05:07:33    -1904.397008        0.013508
    FIRE:  381 05:07:33    -1904.397053        0.013382
    FIRE:  382 05:07:33    -1904.397116        0.013227
    FIRE:  383 05:07:33    -1904.397181        0.013046
    FIRE:  384 05:07:33    -1904.397308        0.012846
    FIRE:  385 05:07:33    -1904.397413        0.012608
    FIRE:  386 05:07:33    -1904.397542        0.012314
    FIRE:  387 05:07:33    -1904.397653        0.011936
    FIRE:  388 05:07:33    -1904.397790        0.011430
    FIRE:  389 05:07:33    -1904.397986        0.010770
    FIRE:  390 05:07:33    -1904.398160        0.009937
    FIRE:  391 05:07:34    -1904.398350        0.008944
    FIRE:  392 05:07:34    -1904.398603        0.007678
    FIRE:  393 05:07:34    -1904.398827        0.006196
    FIRE:  394 05:07:34    -1904.399145        0.005507
    FIRE:  395 05:07:34    -1904.399430        0.005894
    FIRE:  396 05:07:34    -1904.399747        0.005839
    FIRE:  397 05:07:34    -1904.400009        0.015864
    FIRE:  398 05:07:34    -1904.399834        0.040805
    FIRE:  399 05:07:34    -1904.400482        0.004770



```python
SnO2_Pt119 = Trajectory("./SnO2_Pt119.traj")
view_ngl(SnO2_Pt119, representations=["ball+stick"])
```




    HBox(children=(NGLWidget(max_frame=399), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'Sn'…




```python
from ase.io import write
from IPython.display import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


write("output/sno2_pt119_opt.png", SnO2_Pt119[-1], rotation="-90x,0y,0z")

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.imshow(mpimg.imread("output/sno2_pt119_opt.png"))
ax.set_axis_off()
ax.set_title("SnO2 - Pt119 optimized structure")
fig.show()
```


    
![png](output_41_0.png)
    


このように、担体上のNano particle 構造を作成することができました。

ここで、論文で示されているようなchargeの可視化も可能です。&lt;br/&gt;
Nano particleの吸着面に近いところでのみ電荷が０からずれていて、核の外側がプラスに、内側がマイナスになっていることが確認できます。


```python
v = view_ngl(SnO2_Pt119[-1], show_charge=True)
v.show_charge_label()
v
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'Sn', 'Pt'), valu…



このように、SnO2担体上にPt Nano particleが吸着している構造を作成することができました。

実際には、本当にこの構造で正しいのか（微粒子の状態で吸着するのか、それとも全体に均質に伸びて吸着するような形のほうが安定なのかなど)の検証が必要となりますが、このTutorialでは省略します。

こういった構造が作成できる事により、ここから更に後述する反応探索やMDなどを組み合わせることで、

 - 様々な面が露出している触媒下での反応
 - 触媒と担体の境界面での反応

なども扱えるようになる可能性があり、より現実に近い状態での触媒反応シミュレーションができる可能性があります。

参考文献

 - [Calculations of Real-System Nanoparticles Using Universal Neural Network Potential PFP](https://arxiv.org/abs/2107.00963)
 - [Structural Stability of Ruthenium Nanoparticles: A Density Functional Theory Study](https://pubs.acs.org/doi/10.1021/acs.jpcc.7b08672)
 - [Electronic Structure and Phase Stability of PdPt Nanoparticles](https://pubs.acs.org/doi/10.1021/acs.jpclett.5b02753)
 - [Electronic structure and phase stability of Pt3M (M = Co, Ni, and Cu) bimetallic nanoparticles](https://www.sciencedirect.com/science/article/abs/pii/S0927025620303657)
 - [The effect of SnO2(110) supports on the geometrical and electronic properties of platinum nanoparticles](https://link.springer.com/article/10.1007/s42452-019-1478-0)
