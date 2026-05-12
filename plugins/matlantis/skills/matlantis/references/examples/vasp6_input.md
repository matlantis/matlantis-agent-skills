# R2SCANモードに対応するVASPの入力ファイル
`R2SCAN`モードに対応する、VASP6.4.3の入力ファイルを作成する例です。PFPによるDFT計算の再現性を確認するためなどに用います。この例を実行することで`INCAR`, `POSCAR`, `KPOINTS`, `POTCAR`が生成されます。これらのうち、`INCAR`, `POSCAR`, `KPOINTS`をダウンロードし、`POTCAR`をご自身の環境で作成することで、VASPでDFT計算を行うことができます。


```python
import math
import os
from pathlib import Path

import ase.build
import numpy as np
from ase import Atom, Atoms
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.setups import setups_defaults
```


##  準備
### PPの選択
各元素で使用する擬ポテンシャルを定義します。`R2SCAN` 計算モードでは、VASP6.4以降に付属するPBE向けのPAWポテンシャルのバージョン6.4(potpaw_PBE.64)を擬ポテンシャルとして用いています。


```python
setups_defaults["pfp_r2scan"] = {
        # based on ase.calculators.vasp.setups.setups_defaults
        # Alkali and alkali-earth
        "Li": "_sv",
        "Na": "_pv",
        "K": "_sv",
        "Cs": "_sv",
        "Rb": "_sv",
        "Be": "_sv",
        "Mg": "_pv",
        "Ca": "_sv",
        "Sr": "_sv",
        "Ba": "_sv",
        # d-elements, transition metals
        "Sc": "_sv",
        "Y": "_sv",
        "Ti": "_pv",
        "Zr": "_sv",
        "Hf": "_pv",
        "V": "_sv",
        "Nb": "_pv",
        "Ta": "_pv",
        "Cr": "_pv",
        "Mo": "_pv",
        "W": "_sv",
        "Mn": "_pv",
        "Tc": "_pv",
        "Re": "_pv",
        "Fe": "_pv",
        "Co": "",
        "Ni": "_pv",
        "Cu": "_pv",
        "Zn": "",
        "Ru": "_pv",
        "Rh": "_pv",
        "Pd": "",
        "Ag": "",
        "Cd": "",
        "Hg": "",
        "Ir": "",
        "Pt": "",
        "Os": "_pv",
        # Main group
        "Ga": "_d",
        "Ge": "_d",
        "Al": "",
        "As": "",
        "Se": "",
        "Br": "",
        "In": "_d",
        "Sn": "_d",
        "Tl": "_d",
        "Pb": "_d",
        "Bi": "_d",
        # Rare-earth, f-electrons
        "La": "",
        "Ce": "",
    }
```


### ダミーPP
ダミーのPseudo Potentialファイルを用意します。


```python
dummy_pp_path = Path("./dummy_pot")
dummy_pp_pbe_path = dummy_pp_path / "potpaw_PBE"
dummy_pp_pbe_path.mkdir(exist_ok=True, parents=True)

for num in range(Atom("H").number, Atom("Bi").number+1):
    elem = Atom(num).symbol
    pp = setups_defaults["pfp_r2scan"].get(elem, "")
    potcar_path = dummy_pp_pbe_path / f"{elem}{pp}" / "POTCAR"
    potcar_path.parent.mkdir(exist_ok=True)
    with potcar_path.open("w") as f:
        f.write(str(potcar_path.relative_to(dummy_pp_path))+"\n")
```


```python
os.environ["VASP_PP_PATH"] = str(dummy_pp_path)
```


### 関数
k点メッシュを計算するための関数を定義します。


```python
def gen_kpoints_by_spacing(atoms, kspacing=0.22, is_vacuum=(False, False, False)):
    # Generate kpoints based on spacing between k points (kspacing)
    assert len(is_vacuum) == 3

    reciprocal = atoms.cell.reciprocal()
    b = reciprocal.cellpar()[:3]
    num_div = np.ceil(b * 2 * np.pi / kspacing).astype(np.int32)
    # single kpoint in the vacuum direction
    num_div = np.where(is_vacuum, 1, num_div)

    return num_div
```


### パラメータ
VASPでの計算で用いるパラメータを定義します。


```python
vasp_params = {
    "setups": setups_defaults["pfp_r2scan"],
    "xc": "R2SCAN",
    "algo": "all",
    "ediff": 1e-6,
    "enaug": 1360,
    "encut": 680,
    "ismear": 0,
    "ispin": 2,
    "lasph": True,
    "lmaxmix": 6,
    "lmixtau": True,
    "lorbit": 11,
    "lscalu": False,
    "lplane": True,
    "lreal": False,
    "nelm": 200,
    "prec": "Accurate",
    "sigma": 0.05,
    "symprec": 1.0e-8,
    "gamma": True,    
    "ldau_luj": {},
}
```


## 実行例
この例ではFeを対象とします。


```python
atoms = ase.build.bulk('Fe', crystalstructure='bcc', a=2.87, cubic=True)
```


`MAGMOM`の値を定義します。最も安定となる`MAGMOM`を指定してください。Feは強磁性であることが知られているため、今回はこれに対応した値を指定しています。


```python
magmom = [3, 3]
```


k点を計算します。


```python
kpts = gen_kpoints_by_spacing(atoms)
```


真空中の分子や、スラブなど真空層を持つ構造に対していは、真空層がある軸の方向に対して、Γ点のみをとります。これは、a, b, c軸に真空層があるか、否かをis_vacuumを指定することで実現できます。例えば、c軸方向に真空層を持つスラブでは、以下のように指定します。今回は、この条件に当てはまらないため、コメントアウトしています。


```python
# is_vacuum = (False, False, True)
# kpts = gen_kpoints_by_spacing(atoms, is_vacuum=is_vacuum)
```


入力ファイルを出力します。


```python
Vasp(
    directory="./",
    kpts = kpts,
    magmom = magmom,
    **vasp_params
).write_input(atoms)
```


`POSCAR`, `INCAR`, `KPOINTS`, `POTCAR`が出力されます。
VASPでの入力形式に合わせるため、原子のインデックスが元の`atoms`から変わりうることに注意してください。


```python
!cat POSCAR
```

    Fe
     1.0000000000000000
         2.8700000000000001    0.0000000000000000    0.0000000000000000
         0.0000000000000000    2.8700000000000001    0.0000000000000000
         0.0000000000000000    0.0000000000000000    2.8700000000000001
     Fe 
       2
    Cartesian
      0.0000000000000000  0.0000000000000000  0.0000000000000000
      1.4350000000000001  1.4350000000000001  1.4350000000000001



```python
!cat INCAR
```

    ENAUG = 1360.000000
    ENCUT = 680.000000
    SIGMA = 0.050000
    EDIFF = 1.00e-06
    SYMPREC = 1.00e-08
    ALGO = all
    METAGGA = R2SCAN
    PREC = Accurate
    ISMEAR = 0
    ISPIN = 2
    LMAXMIX = 6
    LORBIT = 11
    NELM = 200
    MAGMOM = 2*3.0000
    LASPH = .TRUE.
    LPLANE = .TRUE.
    LSCALU = .FALSE.
    LMIXTAU = .TRUE.
    LREAL = .FALSE.
    LDAU = .TRUE.
    LDAUL = -1
    LDAUU = 0.000
    LDAUJ = 0.000



```python
!cat KPOINTS
```

    KPOINTS created by Atomic Simulation Environment
    0
    Gamma
    10 10 10
    0 0 0



生成される`POTCAR`はダミーであり、実際の`POTCAR`の生成するのに必要なファイルが記述されています。これにもとづき、VASP6.4.3に付属するPPから実際の`POTCAR`を生成してください。


```python
!cat POTCAR
```

    potpaw_PBE/Fe_pv/POTCAR



## 計算結果の読み込み
VASPで計算されたファイルを読み込むことができます。`OUTCAR`をアップロードして以下のように読み込むことで、エネルギーと力を得られます。PFPでのエネルギーの基準と、VASPでのエネルギーの基準は異なることに注意してください。


```python
# import ase.io
# atoms = ase.io.read("OUTCAR")
# print("Energy: {}".format(atoms.get_potential_energy()))
# print(atoms.get_forces())
```
