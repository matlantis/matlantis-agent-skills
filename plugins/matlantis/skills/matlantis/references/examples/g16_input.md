# WB97XDモードに対応するGaussianの入力ファイル
WB97XDモードに対応する、Gaussian16の入力ファイルを作成する例です。PFPによるDFT計算の再現性を確認するためなどに用います。この例を実行することで`opt.single`と`single_point.com`が生成されます。これらをダウンロードして、ご自身の環境のGaussian16でDFT計算を行うことができます。


```python
import ase.build
from ase.calculators.gaussian import Gaussian
```


## スピン多重度
スピン多重度を計算する関数を定義します。singletあるいはdoubletのみを考慮します。


```python
def calc_spin_multiplicity(atoms):
    return sum(atoms.get_atomic_numbers()) % 2 + 1
```


## 計算パラメータ
Gaussianのパラメータを定義します。

### 構造最適化
構造最適化をする場合のパラメータです。`mem`でメモリの最大使用量を、`nprocshared`でコア数を指定します。必要に応じて変更してください。


```python
params_opt = {
    "method": "wb97xd",
    "basis": "6-31g(d)",
    "opt": "maxcycle=256",
    "mem": "4GB",
    "nprocshared": 4,
}
```


### 一点計算
一点計算をする場合のパラメータです。波動関数の対称性が崩れる場合を想定して、非制限DFTを指定し、波動関数のinitial guessの対称性を崩しています。


```python
params_single_point = {
    "method": "uwb97xd",
    "basis": "6-31g(d)",
    "guess": "mix",
    "mem": "4GB",
    "nprocshared": 4,
    "force": None,
}
```


## 入力ファイルの作成
水分子を例として、入力ファイルを作成します。


```python
atoms = ase.build.molecule("H2O")
```


### 構造最適化
以下のcellで入力ファイルを作成できます。


```python
Gaussian(
    label="opt",
    mult = calc_spin_multiplicity(atoms),
    **params_opt,
).write_input(atoms, system_changes=0)
```


`label="opt"`を指定したため、`opt.com`にファイルが出力されます。


```python
!cat opt.com
```

    %mem=4GB
    %nprocshared=4
    #P wb97xd/6-31g(d) ! ASE formatted method and basis
    opt(maxcycle=256)
    
    Gaussian input prepared by ASE
    
    0 1
    O                 0.0000000000        0.0000000000        0.1192620000
    H                 0.0000000000        0.7632390000       -0.4770470000
    H                 0.0000000000       -0.7632390000       -0.4770470000
    
    



### 一点計算
安定構造からずれた構造に対する一点計算のinputを作成します。


```python
atoms = ase.build.molecule("H2O")
atoms.positions *= 1.1
```


```python
Gaussian(
    label="single_point",
    mult = calc_spin_multiplicity(atoms),
    **params_single_point,
).write_input(atoms, system_changes=0)
```


```python
!cat single_point.com
```

    %mem=4GB
    %nprocshared=4
    #P uwb97xd/6-31g(d) ! ASE formatted method and basis
    guess(mix)
    force
    
    Gaussian input prepared by ASE
    
    0 1
    O                 0.0000000000        0.0000000000        0.1311882000
    H                 0.0000000000        0.8395629000       -0.5247517000
    H                 0.0000000000       -0.8395629000       -0.5247517000
    
    



## 計算結果の読み込み
Gaussian16で計算されたファイルを読み込むことができます。例えば`single_point.com`を実行して得られた`single_point.log`をアップロードして以下のように読み込むことで、エネルギーと力を得られます。PFPでのエネルギーの基準と、Gaussian16のエネルギーの基準は異なることに注意してください。


```python
# from ase.io import read
# atoms = read("single_point.log")
# print("Energy: {}".format(atoms.get_potential_energy()))
# print(atoms.get_forces())
```
