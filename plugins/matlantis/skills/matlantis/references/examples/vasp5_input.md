# PBE/PBE_Uモードに対応するVASPの入力ファイル
`PBE_U`モード、あるいは、`PBE`モードに対応する、VASP5.4.4の入力ファイルを作成する例です。PFPによるDFT計算の再現性を確認するためなどに用います。この例を実行することで`INCAR`, `POSCAR`, `KPOINTS`, `POTCAR`が生成されます。これらのうち、`INCAR`, `POSCAR`, `KPOINTS`をダウンロードし、`POTCAR`をご自身の環境で作成することで、VASPでDFT計算を行うことができます。


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
各元素で使用する擬ポテンシャルを定義します。`PBE` もしくは `PBE_U` 計算モードでは、VASP5.4以降に付属するPBE向けのPAWポテンシャルのバージョン5.4(potpaw_PBE.54)を擬ポテンシャルとして用いています。


```python
setups_defaults["pfp_pbe"] = {
        # This dictionary is partially based on the one defined in ase.calculators.vasp.setups.setups_defaults (ase==3.22.1, LGPL license).
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
        "Bi": "_d",  # ?
        # Rare-earth, f-electrons
        "La": "",
        "Ce": "",
        "Pr": "_3",
        "Nd": "_3",
        "Pm": "_3",
        "Sm": "_3",
        "Eu": "",
        "Gd": "",
        "Tb": "_3",
        "Dy": "_3",
        "Ho": "_3",
        "Er": "_3",
        "Tm": "_3",
        "Yb": "",
        "Lu": "_3",
        # from VASP recommended
        "Po": "_d",
        "At": "",
        "Rn": "",
        "Fr": "_sv",
        "Ra": "_sv",
        "Ac": "",
        "Th": "",
        "Pa": "",
        "U": "",
        "Np": "",
        "Pu": "",
        "Am": "",
        "Cm": "",
    }
```


### ダミーPP
ダミーのPseudo Potentialファイルを用意します。


```python
dummy_pp_path = Path("./dummy_pot")
dummy_pp_pbe_path = dummy_pp_path / "potpaw_PBE"
dummy_pp_pbe_path.mkdir(exist_ok=True, parents=True)

for num in range(Atom("H").number, Atom("Cm").number+1):
    elem = Atom(num).symbol
    pp = setups_defaults["pfp_pbe"].get(elem, "")
    potcar_path = dummy_pp_pbe_path / f"{elem}{pp}" / "POTCAR"
    potcar_path.parent.mkdir(exist_ok=True)
    with potcar_path.open("w") as f:
        f.write(str(potcar_path.relative_to(dummy_pp_path))+"\n")
```


```python
os.environ["VASP_PP_PATH"] = str(dummy_pp_path)
```


### 関数
k点メッシュと`LMAXMIX`を計算するための関数を定義します。


```python
def loground(v):
    v_ceil = np.ceil(v)
    v_floor = np.max([np.floor(v), np.ones(v.shape)], axis=0)
    is_ceil = np.absolute(np.log(v_ceil) - np.log(v)) &lt; np.absolute(
        np.log(v_floor) - np.log(v)
    )
    v_rounded: np.ndarray = np.where(is_ceil, v_ceil, v_floor)
    return v_rounded

def is_hexagonal(atoms: Atoms, angle_tol=5, length_tol=0.01):
    # based on pymatgen.core.lattice.Lattice.is_hexagonal
    cellpar = atoms.cell.cellpar()
    lengths = cellpar[:3]
    angles = cellpar[3:]

    right_angles = np.absolute(angles - 90) &lt; angle_tol
    hex_angles = np.logical_or(
        np.absolute(angles - 60) &lt; angle_tol, np.absolute(angles - 120) &lt; angle_tol
    )

    if np.sum(right_angles) == 2 and np.sum(hex_angles) == 1:
        right_lengths = lengths[right_angles]
        return bool(np.min(right_lengths) / np.max(right_lengths) &gt;= (1 - length_tol))
    else:
        return False

def gen_kpoints(atoms, kppa=1000):
    # generate kpoints based on kpoints per atom (kppa)
    # based on pymatgen.io.inputs.automatic_density and change to use loground instead of
    # floor to calculate num_div
    cell = atoms.get_cell()
    lengths = cell.lengths()

    if math.fabs((math.floor(kppa ** (1 / 3) + 0.5)) ** 3 - kppa) &lt; 1:
        kppa += kppa * 0.01
    ngrid = kppa / len(atoms)
    mult = (ngrid * np.prod(lengths)) ** (1 / 3)
    div = mult / lengths
    num_div = loground(np.where(div &gt;= 1, div, 1)).astype(np.int32)

    is_hexagonal_ = is_hexagonal(atoms)
    has_odd = any(i % 2 == 1 for i in num_div)
    # vasp strongly suggest to use gamma-centered k-points for hexagonal lattice
    gamma = has_odd or is_hexagonal_

    return num_div, gamma

def gen_lmaxmix(atoms, ldau_luj):
        lmix = []
        ldau_luj = ldau_luj or {}
        for elem in set(atoms.get_chemical_symbols()):
            if elem in ldau_luj and ldau_luj[elem]["L"] &gt;= 0:
                atomno = Atom(elem).number
                if atomno &gt; 56:
                    lmix.append(6)
                elif atomno &gt; 20:
                    lmix.append(4)
                else:
                    lmix.append(2)
            else:
                lmix.append(2)
        return max(lmix)
```


### パラメータ
VASPでの計算で用いるパラメータを定義します。


```python
ldau_luj = {
    "V": {"L": 2, "U": 3.25, "J": 0},
    "Cr": {"L": 2, "U": 3.7, "J": 0},
    "Mn": {"L": 2, "U": 3.9, "J": 0},
    "Fe": {"L": 2, "U": 5.3, "J": 0},
    "Co": {"L": 2, "U": 3.32, "J": 0},
    "Ni": {"L": 2, "U": 6.2, "J": 0},
    "Cu": {"L": 2, "U": 4.0, "J": 0},
    "Mo": {"L": 2, "U": 4.38, "J": 0},
    "W":  {"L": 2, "U": 6.2, "J": 0}
}
```


Hubbard補正なし(`PBE`)に対応するインプットを作成する場合は、下のセルのように`ldau_luj`を空にします。今回は、Hubbard補正あり(`PBE_U`)に対応するインプットを作成するため、コメントアウトしています。


```python
# ldau_luj = {}
```


```python
vasp_params = {
    "setups": setups_defaults["pfp_pbe"],
    "xc": "PBE",
    "algo": "Normal",
    "ediff": 0.0001,
    "encut": 520,
    "ismear": 0,
    "ispin": 2,
    "lorbit": 11,
    "lreal": "Auto",
    "nelm": 200,
    "prec": "Accurate",
    "sigma": 0.05,
    "ldau_luj": ldau_luj
}
```


## 実行例
この例ではFeOを対象とします。


```python
atoms = Atoms(
    "Fe4O4",
    scaled_positions=[
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.5, 0.5],
    ],
    cell=[4.4, 4.4, 4.4],
)
```


`MAGMOM`の値を定義します。最も安定となる`MAGMOM`を指定してください。FeOは強磁性であることが知られているため、今回はこれに対応した値を指定しています。


```python
magmom = [5, 5, 5, 5, 1, 1, 1, 1]
```


k点と`LMAXMIX`を計算します。


```python
kpts, gamma = gen_kpoints(atoms)
lmaxmix = gen_lmaxmix(atoms, ldau_luj)
```


セルのすべての方向に真空層を持つ場合は、Γ点一点計算を指定します。今回の例はこの条件に当てはまらないため、コメントアウトしています。


```python
# kpts = [1,1,1]
# gamma = True
```


入力ファイルを出力します。


```python
Vasp(
    directory="./",
    kpts = kpts,
    gamma = gamma,
    lmaxmix = lmaxmix,
    magmom = magmom,
    **vasp_params
).write_input(atoms)
```


`POSCAR`, `INCAR`, `KPOINTS`, `POTCAR`が出力されます。
VASPでの入力形式に合わせるため、原子のインデックスが元の`atoms`から変わりうることに注意してください。


```python
!cat POSCAR
```

    Fe O 
     1.0000000000000000
         4.4000000000000004    0.0000000000000000    0.0000000000000000
         0.0000000000000000    4.4000000000000004    0.0000000000000000
         0.0000000000000000    0.0000000000000000    4.4000000000000004
     Fe  O  
       4   4
    Cartesian
      0.0000000000000000  0.0000000000000000  0.0000000000000000
      2.2000000000000002  2.2000000000000002  0.0000000000000000
      2.2000000000000002  0.0000000000000000  2.2000000000000002
      0.0000000000000000  2.2000000000000002  2.2000000000000002
      0.0000000000000000  0.0000000000000000  2.2000000000000002
      2.2000000000000002  0.0000000000000000  0.0000000000000000
      0.0000000000000000  2.2000000000000002  0.0000000000000000
      2.2000000000000002  2.2000000000000002  2.2000000000000002



```python
!cat INCAR
```

    ENCUT = 520.000000
    SIGMA = 0.050000
    EDIFF = 1.00e-04
    ALGO = Normal
    GGA = PE
    PREC = Accurate
    ISMEAR = 0
    ISPIN = 2
    LMAXMIX = 4
    LORBIT = 11
    NELM = 200
    MAGMOM = 4*5.0000 4*1.0000
    LREAL = Auto
    LDAU = .TRUE.
    LDAUL = 2 -1
    LDAUU = 5.300 0.000
    LDAUJ = 0.000 0.000



```python
!cat KPOINTS
```

    KPOINTS created by Atomic Simulation Environment
    0
    Gamma
    5 5 5
    0 0 0



生成される`POTCAR`はダミーであり、実際の`POTCAR`の生成するのに必要なファイルが記述されています。これにもとづき、VASP5.4.4に付属するPPから実際の`POTCAR`を生成してください。


```python
!cat POTCAR
```

    potpaw_PBE/Fe_pv/POTCAR
    potpaw_PBE/O/POTCAR



## 計算結果の読み込み
VASPで計算されたファイルを読み込むことができます。`OUTCAR`をアップロードして以下のように読み込むことで、エネルギーと力を得られます。PFPでのエネルギーの基準と、VASPでのエネルギーの基準は異なることに注意してください。


```python
# import ase.io
# atoms = ase.io.read("OUTCAR")
# print("Energy: {}".format(atoms.get_potential_energy()))
# print(atoms.get_forces())
```
