# 形成エネルギーの計算

このノートブックではPFPを用いた形成エネルギーの計算方法を説明します。

形成エネルギーは、化合物がその構成要素から生成されるときのエネルギー差で、一般的には

$$\,E_\mathrm{form} = E_\mathrm{compound} - \sum_i n_i E_i$$

と定義されます。ここで $E_{\mathrm{compound}}$ は化合物の全エネルギー、$E_i$ は元素 $i$ の基準エネルギー、$n_i$ はその元素の化合物中の個数です。形成エネルギーによって、化合物の安定性を定量的に評価できます。

このノートブックは二つの節から構成されています。 

1) 一般的な形成エネルギーの計算方法。PBE calc_modeの例を載せていますが、R2SCAN calc_modeでも同様です。
2) 一部の遷移金属酸化物など、Materials Projectで提案されている混合補正が必要な場合の計算方法。


```python
!pip install ase pymatgen pfp-api-client

```


```python
from ase.build import molecule, bulk
from ase.optimize import LBFGS
from ase.filters import FrechetCellFilter
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility

```


```python
MODEL_VERSION = "v8.0.0"

```


## セクション1: 計算例

ここではNaClの形成エネルギーを計算します。



### 構造最適化

形成エネルギーの計算の前に、構造最適化を行う必要があります。



```python
calc_pbe = ASECalculator(Estimator(model_version=MODEL_VERSION, calc_mode=EstimatorCalcMode.PBE))

def relax_atoms(init_atoms, calc, fmax=0.001, steps=500, cell=False):
    atoms = init_atoms.copy()
    atoms.calc = calc
    opt_target = FrechetCellFilter(atoms) if cell else atoms
    opt = LBFGS(opt_target)
    opt.run(fmax=fmax, steps=steps)
    return atoms.get_potential_energy()

```


```python
E_na = relax_atoms(bulk("Na", crystalstructure="bcc"), calc_pbe, cell=True)
print(f"E_na = {E_na}")

E_nacl = relax_atoms(bulk("NaCl", crystalstructure="rocksalt", a=5.64), calc_pbe, cell=True)
print(f"E_nacl = {E_nacl}")

atoms_chlorine = molecule("Cl2")
atoms_chlorine.set_cell([15, 15, 15])
atoms_chlorine.center()
E_chlorine = relax_atoms(atoms_chlorine, calc_pbe, cell=False)
print(f"E_chlorine = {E_chlorine}")

```

           Step     Time          Energy          fmax
    LBFGS:    0 04:56:47       -1.078337        0.044707
    LBFGS:    1 04:56:47       -1.078419        0.041401
    LBFGS:    2 04:56:47       -1.078906        0.001208
    LBFGS:    3 04:56:48       -1.078907        0.000031
    E_na = -1.078906570233863
           Step     Time          Energy          fmax
    LBFGS:    0 04:56:48       -6.167160        0.118572
    LBFGS:    1 04:56:48       -6.167752        0.110826
    LBFGS:    2 04:56:48       -6.170927        0.025966
    LBFGS:    3 04:56:48       -6.171117        0.001880
    LBFGS:    4 04:56:48       -6.171115        0.000020
    E_nacl = -6.171115245799412
           Step     Time          Energy          fmax
    LBFGS:    0 04:56:48       -2.787997        0.431318
    LBFGS:    1 04:56:48       -2.792039        0.223466
    LBFGS:    2 04:56:48       -2.793467        0.010386
    LBFGS:    3 04:56:48       -2.793470        0.000329
    E_chlorine = -2.7934695616778527



### 形成エネルギーの計算



```python
E_formation = E_nacl - E_na - E_chlorine / 2.0
E_formation_per_atom = E_formation / 2.0

print(f"Formation energy (per atom): {E_formation_per_atom:.6f} eV/atom")

```

    Formation energy (per atom): -1.847737 eV/atom



## セクション2: 遷移金属酸化物の形成エネルギー

ここでは遷移金属酸化物であるFeOの形成エネルギーを計算します。

一部の遷移金属酸化物に関しては、エネルギーを求める際に電子相関に関するHubbard U補正を考慮する必要があります。PFPはU補正の有無で別々のcalc_modeを提供しており、酸化物に関してはU補正あり、その要素となる金属などに対してはU補正なしで計算を行うことができます。

ただし、異なるcalc_modeを用いて算出したエネルギーは原点が異なるため、そのままでは形成エネルギーの計算に使用できません。そのため、Materials Projectが提唱する経験的な補正（アニオン補正やGGA/GGA+Uの混合補正）を適用して、異なるcalc_mode間でエネルギーを比較可能にします。



```python
metal = "Fe"
oxide = "FeO"

```


### 構造最適化

金属FeおよびO2はPBE、酸化物FeOはPBE+Uで最適化します。



```python
calc_pbeu = ASECalculator(Estimator(model_version=MODEL_VERSION, calc_mode=EstimatorCalcMode.PBE_U))

```


```python
atoms_metal = bulk(metal, crystalstructure="bcc")
E_metal = relax_atoms(atoms_metal, calc_pbe, cell=True)
print(f"E_metal = {E_metal}")

atoms_oxide = bulk(oxide, crystalstructure="rocksalt", a=4.33)
E_oxide = relax_atoms(atoms_oxide, calc_pbeu, cell=True)
print(f"E_oxide = {E_oxide}")

atoms_oxygen = molecule("O2")
atoms_oxygen.set_cell([15, 15, 15])
atoms_oxygen.center()
E_oxygen = relax_atoms(atoms_oxygen, calc_pbe, cell=False)
print(f"E_oxygen = {E_oxygen}")

```

           Step     Time          Energy          fmax
    LBFGS:    0 04:56:49       -4.989421        0.325415
    LBFGS:    1 04:56:49       -4.993107        0.203793
    LBFGS:    2 04:56:49       -4.995668        0.029107
    LBFGS:    3 04:56:49       -4.995722        0.002769
    LBFGS:    4 04:56:49       -4.995722        0.000177
    E_metal = -4.995721653298855
           Step     Time          Energy          fmax
    LBFGS:    0 04:56:49       -8.978694        0.203821
    LBFGS:    1 04:56:49       -8.980284        0.159208
    LBFGS:    2 04:56:49       -8.982294        0.049291
    LBFGS:    3 04:56:49       -8.982443        0.009906
    LBFGS:    4 04:56:49       -8.982450        0.000343
    E_oxide = -8.982449645198383
           Step     Time          Energy          fmax
    LBFGS:    0 04:56:49       -5.888375        0.813040
    LBFGS:    1 04:56:50       -5.888474        0.846621
    LBFGS:    2 04:56:50       -5.893245        0.031615
    LBFGS:    3 04:56:50       -5.893250        0.001078
    LBFGS:    4 04:56:50       -5.893250        0.000006
    E_oxygen = -5.893249693761836



### pymatgen ComputedEntryの作成

Materials Projectで提案されている経験的補正はDFTの計算結果に対して提案されている手法です。PFPのエネルギー原点はDFTのエネルギー原点とは異なるため、そのままではこの補正を適用することができません。そのため、事前にエネルギー原点の差であるshift energyを加算してからComputedEntryに渡します。



```python
shift_energies_pbe = Estimator(model_version=MODEL_VERSION, calc_mode=EstimatorCalcMode.PBE).get_shift_energy_table()
shift_energies_pbeu = Estimator(model_version=MODEL_VERSION, calc_mode=EstimatorCalcMode.PBE_U).get_shift_energy_table()

E_metal = E_metal + sum(shift_energies_pbe[i] for i in atoms_metal.get_atomic_numbers())
E_oxygen = E_oxygen + sum(shift_energies_pbe[i] for i in atoms_oxygen.get_atomic_numbers())
E_oxide = E_oxide + sum(shift_energies_pbeu[i] for i in atoms_oxide.get_atomic_numbers())
```


```python
adapter = AseAtomsAdaptor()

def get_clean_structure(ase_atoms):
    atoms_copy = ase_atoms.copy()
    atoms_copy.calc = None
    return adapter.get_structure(atoms_copy)

```


```python
struct_metal = get_clean_structure(atoms_metal)
entry_metal = ComputedEntry(
    composition=struct_metal.composition,
    energy=E_metal,
    parameters={"hubbards": {}, "run_type": "GGA"},
    entry_id=f"{metal}_PBE"
)

struct_oxygen = get_clean_structure(atoms_oxygen)
entry_oxygen = ComputedEntry(
    composition=struct_oxygen.composition,
    energy=E_oxygen,
    parameters={"hubbards": {}, "run_type": "GGA"},
    entry_id="O2_PBE"
)

struct_oxide = get_clean_structure(atoms_oxide)
entry_oxide = ComputedEntry(
    composition=struct_oxide.composition,
    energy=E_oxide,
    parameters={
        "hubbards": {"Fe": 5.3}, 
        "run_type": "GGA+U"
    },
    entry_id=f"{oxide}_PBE+U"
)

```


### Materials Projectの経験的補正の適用



```python
compat = MaterialsProject2020Compatibility(check_potcar=False)

processed_entries = [e for e in compat.process_entries([entry_metal, entry_oxygen, entry_oxide]) if e is not None]

```


### 形成エネルギーの計算

補正後のエネルギーを使って、酸化物の形成エネルギーを求めます。



```python
entry_map = {e.entry_id: e for e in processed_entries}

E_formation = (
    entry_map[f"{oxide}_PBE+U"].energy
    - entry_map[f"{metal}_PBE"].energy
    - entry_map["O2_PBE"].energy / 2.0
)
E_formation_per_atom = E_formation / 2.0

print(f"Formation energy (per atom): {E_formation_per_atom:.6f} eV/atom")

```

    Formation energy (per atom): -1.439854 eV/atom



この値は、Materials Projectに記載されているFeOの形成エネルギー $-1.482$ eV/atom に近い値となっています。



## 参考文献
- [Material Project: Anion and GGA/GGA+U Mixing](https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/anion-and-gga-gga+u-mixing)
- [Material Project: Adding Energy Correction to Custom Entries](https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/adding-energy-corrections-to-custom-entries)
- [Pymatgen document: MaterialsProject2020Compatibility](https://pymatgen.org/pymatgen.entries.html#pymatgen.entries.compatibility.MaterialsProject2020Compatibility)
- [A framework for quantifying uncertainty in DFT energy corrections](https://www.nature.com/articles/s41598-021-94550-5)
- [A flexible and scalable scheme for mixing computed formation energies from different levels of theory](https://www.nature.com/articles/s41524-022-00881-w)

