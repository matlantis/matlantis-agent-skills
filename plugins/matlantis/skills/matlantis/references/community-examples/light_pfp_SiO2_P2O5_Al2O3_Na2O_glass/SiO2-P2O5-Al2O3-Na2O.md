Copyright Preferred Networks inc. as contributors to Matlantis contrib project.

# LightPFP: SiO2-P2O5-Al2O3-Na2O glass

This notebook makes a lightPFP model for SiO2-P2O5-Al2O3-Na2O glass system.
The composition of SiO2 varies between 20.69% and 70.69%
The composition of P2O5 vaires between 0% and 50%

We follow the paper ["Ab initio molecular dynamics simulation of structural and elastic properties of SiO2–P2O5–Al2O3–Na2O glass"](https://doi.org/10.1111/jace.18614)

The total time cost:
* LightPFP model generation (section 1 + 2 + 3 + 4): 25 hours
* PFP reference MD (section 5): about 38 hours depends on PFP Load status
* LightPFP reference MD (section 6): 5 hours

Note: since section 5 (PFP reference MD) is independent from other sections and very time comsuming, I recommend to run it in parallel with other parts.


```python
model_version = "v7.0.0"
calc_mode = "pbe"
```

## Setup


```python
! pip install light-pfp-client==1.0.0 light-pfp-data==1.0.0 light-pfp-evaluate==1.0.0 light-pfp-autogen==0.1.3
```

## 1. Initial structure

* Define the functions to generate the random structure of SiO2-P2O5-Al2O3-Na2O glass networks.
* Si, P, Al, Na, and O atoms are put into a simulation box with random positions.
* Compositions (x 0 ~ 50):
    * SiO2: 79.69-x %
    * P2O5: x %
    * Al2O3: 13.79 %
    * Na2O: 15.52 %
* It outputs small-size (\~250 atoms) and large-size (\~500 atoms) glass structures.

(In glass, the arrangement of atoms should conform to valence bond theory. For example, a Si atom forms a tetrahedral structure with four O atoms. However, I am currently using random initial atomic positions, which clearly violates this principle, making the initial structure very unreasonable. In the future, more suitable methods for generating initial structures could be considered.)


```python
from typing import List
import numpy as np
from pathlib import Path
from pfcc_extras.liquidgenerator.liquid_generator import LiquidGenerator
from ase import Atoms
from ase import units
from ase.io import read
from ase.md.npt import NPT
from ase.md.langevin import Langevin

from IPython.display import clear_output


def get_glass(n_Si, n_P, n_Al, n_Na, density=2.0):
    n_O = n_Si * 2 + n_Al//2 * 3 + n_Na//2 + n_P//2 * 5
    composition = []
    composition += [Atoms(symbols=["Si"])] * n_Si
    composition += [Atoms(symbols=["P"])] * n_P
    composition += [Atoms(symbols=["Al"])] * n_Al
    composition += [Atoms(symbols=["Na"])] * n_Na
    composition += [Atoms(symbols=["O"])] * n_O
    liquid_generator = LiquidGenerator(engine="torch", composition=composition, density=density)
    atoms = liquid_generator.run(epochs=100)
    clear_output()
    return atoms
    

def get_random_glass(size="large"):
    assert size in ["small", "large"]
    if size == "small":
        x = np.random.randint(0, 30)
        n_Al = 16
        n_Na = 18
        n_Si = 41-x
        n_P = 2*x
    elif size == "large":
        x = np.random.randint(0, 60)
        n_Al = 32
        n_Na = 36
        n_Si = 82-2*x
        n_P = 4*x
    density = np.random.uniform(2.0, 2.4)
    atoms = get_glass(n_Si, n_P, n_Al, n_Na, density)
    return atoms

```

## 2. Initial structures
* Download SiO2, Al2O3, Na2O and P2O5 crystal structures from Materials Project for initial dataset generation.
* These structures provides the most energetic stable motif of the glass components.


```python
from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
import os
import numpy as np

def download_mp_structures(api_key, formula, limits=3):
    with MPRester(api_key) as mpr:
        docs = mpr.summary.search(
            formula=formula,
            theoretical=False,
            fields=["material_id", "structure", "energy_above_hull"]
        )
        if len(docs) &gt; limits:
            energy_above_hull = []
            for doc in docs:
                energy_above_hull.append(doc.energy_above_hull)
            ind = np.argsort(energy_above_hull)[:limits]
        else:
            ind = np.arange(len(docs))
        
        atoms_list = []
        mp_id_list = []
        for i in ind:
            doc = docs[i]
            structure = doc.structure
            atoms_list.append(AseAtomsAdaptor.get_atoms(structure))
            mp_id_list.append(doc.material_id)
        return atoms_list, mp_id_list
    

def expand_atoms(atoms, length=10.0):
    par = atoms.get_cell().cellpar()[:3]
    supercell = np.round(length / par, 0).astype(int)
    return atoms * supercell
```


```python
api_key = "" # Please input your api_key, https://next-gen.materialsproject.org/api

structure_dir = Path("structures")
structure_dir.mkdir(parents=True, exist_ok=True)

all_structures = []

if len(list(structure_dir.glob("*.cif"))) &lt; 3:
    atoms_list, mp_id_list = download_mp_structures(api_key, "SiO2")
    all_structures += atoms_list
    for atoms, mp_id in zip(atoms_list, mp_id_list):
        atoms.write(structure_dir / f"SiO2_{mp_id}.cif")
    
    atoms_list, mp_id_list = download_mp_structures(api_key, "Al2O3")
    all_structures += atoms_list
    for atoms, mp_id in zip(atoms_list, mp_id_list):
        atoms.write(structure_dir / f"Al2O3_{mp_id}.cif")
    
    atoms_list, mp_id_list = download_mp_structures(api_key, "Na2O")
    all_structures += atoms_list
    for atoms, mp_id in zip(atoms_list, mp_id_list):
        atoms.write(structure_dir / f"Na2O_{mp_id}.cif")
    
    atoms_list, mp_id_list = download_mp_structures(api_key, "P2O5")
    all_structures += atoms_list
    for atoms, mp_id in zip(atoms_list, mp_id_list):
        atoms.write(structure_dir / f"P2O5_{mp_id}.cif")
else:
    for f in structure_dir.glob("*.cif"):
        all_structures.append(read(f))
```

## 2. Initial dataset
* Make dataset with MD sampling
* Initial dataset includes SiO2, Al2O3, Na2O and P2O5 crystals downloaded from MP and 5 additional random structure generated by `get_random_glass` function.


```python
from light_pfp_data.utils.dataset import H5DatasetWriter

init_dataset_dir = Path("init_dataset")
init_dataset_dir.mkdir(parents=True, exist_ok=True)
initial_dataset = init_dataset_dir / "init.h5"
```


```python
from unittest.mock import patch
from pfp_api_client import Estimator, ASECalculator
from light_pfp_data.sample import sample_md


if initial_dataset.exists():
    print(f"Dataset file {initial_dataset} already exists. Skip generating initial dataset.")
    with patch('builtins.input', return_value='y'):
        dataset = H5DatasetWriter(initial_dataset, mode="append") # automatically input 'y' when running in background
else:
    print(f"Dataset file {initial_dataset} is created. Start generating initial structures.")
    with patch('builtins.input', return_value='y'):
        dataset = H5DatasetWriter(initial_dataset)

    # Initialize estimator and calculator
    estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
    calc = ASECalculator(estimator)
    
    from concurrent.futures import as_completed, ThreadPoolExecutor
    from tqdm.auto import tqdm
    
    futures = []
    pbar = tqdm(desc="Total progress", total=0, leave=True)
    with ThreadPoolExecutor(max_workers=8) as executor:
        for atoms in all_structures:
            futures += sample_md(
                input_structure=expand_atoms(atoms),
                calculator=calc,
                dataset=dataset,
                supercell=(1, 1, 1),
                sampling_temp=[500.0, 1000.0, 1500.0],
                sampling_steps=[2000, 2000, 2000],
                sampling_interval=[100, 100, 100],
                ensemble="nvt",
                executor=executor,
                pbar=pbar
            )
        for _ in range(5):
            futures += sample_md(
                input_structure=get_random_glass("small"),
                calculator=calc,
                dataset=dataset,
                supercell=(1, 1, 1),
                sampling_temp=[500.0, 1000.0, 1500.0],
                sampling_steps=[2000, 2000, 2000],
                sampling_interval=[100, 100, 100],
                ensemble="nvt",
                executor=executor,
                pbar=pbar
            )
```


```python
dataset.h5.close()
```

## 3. Active learning
* Run active learning
* The structure is quite diverse (5 elements, random &amp; unnatural initial structures, variable compositions) which makes the LightPFP models have large Force MAE. Accordingly, the sample ceriterion is changed to "dF_min_coef = 8.0" to avoid over-sampling. (The default value "dF_min=1.0" might be to strict for this situation). 


```python
import pathlib
import logging
from light_pfp_autogen.active_learning import ActiveLearning
from light_pfp_autogen.config import ActiveLearningConfig, TrainConfig, SampleConfig, CommonConfig, MTPConfig


logging.basicConfig(level=logging.INFO)

active_learning_config = ActiveLearningConfig(
    task_name = "glass-test-2",
    work_dir = "./autogen_workdir_2",
    pfp_model_version = model_version,
    pfp_calc_mode = calc_mode,
    init_dataset = [
        "init_dataset/init.h5",
    ],
    train_config = TrainConfig(
        common_config = CommonConfig(max_forces=50.0),
        mtp_config = MTPConfig(pretrained_model="ALL_ELEMENTS_SMALL_NN_6")
    ),
    sample_config = SampleConfig(
        dE_min_coef = 3.0,
        dE_max_coef = 20.0,
        dF_min_coef = 8.0,
        dF_max = 50.0,
        dS_min_coef = 3.0,
        dS_max_coef = 20.0
    )
)

active_learning = ActiveLearning(active_learning_config)
active_learning.initialize()
```

* 0 ~ 4 iteration:
    * Get a random initial structure and just run NVT MD at 300~2000 K for 5 ps.
    * Repeat for 50 time in each iteration.
    * Total MD steps: 250,000
* With the small size initial glass structures ( ~ 250 atoms )


```python
import numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from light_pfp_data.utils.atoms_utils import convert_atoms_to_upper
from light_pfp_autogen.context import DataCollectionContext


class TemperatureControl:
    def __init__(self, md, temperature, cooling_rate):
        # cooling rate: K/steps
        self.md = md
        self.cooling_rate = cooling_rate
        self.init_step = self.md.nsteps
        self.init_temp = temperature
    
    def __call__(self):
        delta = (self.md.nsteps - self.init_step) * self.cooling_rate
        temp = self.init_temp - delta
        self.md.set_temperature(temperature_K=temp)

    
for i in range(active_learning.iter, 5):
    print(f"Current active iteration: {i}")
    for _ in range(50):
        atoms = get_random_glass("small")
        temperature = np.random.uniform(300.0, 2000.0)
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
        clear_output()
        atoms = convert_atoms_to_upper(atoms)

        md = Langevin(atoms, units.fs, temperature_K=temperature, friction=0.1)
        with DataCollectionContext(md=md, interval=100, max_samples=20):
            md.run(steps=5000)

        clear_output()
    active_learning.update()
```

* 5~9 iteration:
    * Get a random initial structure and run NVT MD at 2000 K for 10 ps.
    * Aneal the temperature to 500K with a cooling rate of 0.1K/s
    * Equilibrium the system at 500K for 5 ps.
    * Repeat for 10 time in each iteration.
    * Total MD steps: 300,000
    
* With the large size initial glass structures ( ~ 500 atoms )
* We are getting closer to what we want to simulate, i.e. generate glass structure with melt-quenching method.


```python
for i in range(active_learning.iter, 10):
    print(f"Current active iteration: {i}")
    for _ in range(10):
        atoms = get_random_glass("large")
        MaxwellBoltzmannDistribution(atoms, temperature_K=2000.0)
        clear_output()
        atoms = convert_atoms_to_upper(atoms)

        md = Langevin(atoms, units.fs, temperature_K=2000.0, friction=0.1)
        with DataCollectionContext(md=md, interval=100, max_samples=50):
            md.run(steps=10000)
        
        md = NPT(
            atoms, 
            units.fs, 
            temperature_K=2000.0, 
            externalstress=units.bar,
            ttime=20.0 * units.fs,
            pfactor=2e6 * units.GPa * (units.fs**2)
        )
        with DataCollectionContext(md=md, interval=100, max_samples=50):
            temperature_control = TemperatureControl(
                md, 2000.0, 0.1
            )
            md.attach(temperature_control, interval=10)
            md.run(steps=15000)

        md = NPT(
            atoms, 
            units.fs, 
            temperature_K=500.0, 
            externalstress=units.bar,
            ttime=20.0 * units.fs,
            pfactor=2e6 * units.GPa * (units.fs**2)
        )
        with DataCollectionContext(md=md, interval=100, max_samples=20):
            md.run(steps=5000)
            
        clear_output()
    active_learning.update()
```

* 10~14 iteration:
    * Get a random initial structure and run NVT MD at 2000 K for 10 ps.
    * Aneal the temperature to 500K with a cooling rate of 0.05K/s (slower cooling rate)
    * Equilibrium the system at 500K for 5 ps.
    * Repeat for 5 time in each iteration.
    * Total MD steps: 225,000
    
* With the (2x2x2) supercell of the small size initial glass structures ( ~ 2000 atoms )


```python
for i in range(active_learning.iter, 15):
    print(f"Current active iteration: {i}")
    for _ in range(5):
        atoms = get_random_glass("small") * (2,2,2)  
        MaxwellBoltzmannDistribution(atoms, temperature_K=2000.0)
        clear_output()
        atoms = convert_atoms_to_upper(atoms)

        md = Langevin(atoms, units.fs, temperature_K=2000.0, friction=0.1)
        with DataCollectionContext(md=md, interval=100, max_samples=50):
            md.run(steps=10000)
        
        md = NPT(
            atoms, 
            units.fs, 
            temperature_K=2000.0, 
            externalstress=units.bar,
            ttime=20.0 * units.fs,
            pfactor=2e6 * units.GPa * (units.fs**2)
        )
        with DataCollectionContext(md=md, interval=100, max_samples=50):
            temperature_control = TemperatureControl(
                md, 2000.0, 0.05
            )
            md.attach(temperature_control, interval=10)
            md.run(steps=30000)

        md = NPT(
            atoms, 
            units.fs, 
            temperature_K=500.0, 
            externalstress=units.bar,
            ttime=20.0 * units.fs,
            pfactor=2e6 * units.GPa * (units.fs**2)
        )
        with DataCollectionContext(md=md, interval=100, max_samples=20):
            md.run(steps=5000)
            
        clear_output()
    active_learning.update()
```


```python
active_learning.print_md_statistics()
```

# 4. Post active training
* Train a final LightPFP model with all the dataset collected in active learing.
* Pretrained model: ALL_ELEMENTS_SMALL_NN_6


```python
from light_pfp_autogen.utils import submit_training_job, check_training_job_status, estimate_epoch


epoch = estimate_epoch(active_learning.datasets_list, 2)

train_config_dict = {
    "common_config": {
        "total_epoch": epoch,
        "max_forces": 50.0
    },
    "mtp_config": {
        "pretrained_model": "ALL_ELEMENTS_SMALL_NN_6"
    },
}


train_config = TrainConfig.from_dict(
    train_config_dict
)

model_id = submit_training_job(
    train_config,
    active_learning.datasets_list,
    "glass-test-2-final-small",
)

status = check_training_job_status(model_id)
print(f"Training job {model_id} status: {status}")
```

# 5. Run MD with PFP

**I suggest to run this section in a different notebook because it is independent from the main task**

* To collect evaluate the LightPFP model, we run generate glass structure with 5 different compositions with PFP.
* Melt-quencting method is used.


```python
# from pfp_api_client import Estimator, ASECalculator

# def md_simulation(atoms, name):
#     calc = ASECalculator(Estimator(model_version="v7.0.0", calc_mode="pbe"))
#     atoms = atoms * (3,3,3)
#     atoms.calc = calc
#     MaxwellBoltzmannDistribution(atoms, temperature_K=2000.0)    
#     traj = Trajectory(f"pfp_md/{name}.traj", "w", atoms = atoms)
#     md = Langevin(atoms, units.fs, temperature_K=2000.0, friction=0.1)
#     md.attach(traj, interval=1000)
#     md.run(steps=10000)
    
#     md = NPT(
#         atoms, 
#         units.fs, 
#         temperature_K=2000.0, 
#         externalstress=units.bar,
#         ttime=20.0 * units.fs,
#         pfactor=2e6 * units.GPa * (units.fs**2)
#     )
#     temperature_control = TemperatureControl(
#         md, 2000.0, 0.02
#     )
#     md.attach(temperature_control, interval=10)
#     md.attach(traj, interval=1000)
#     md.run(steps=75000)
    
#     md = NPT(
#         atoms, 
#         units.fs, 
#         temperature_K=500.0, 
#         externalstress=units.bar,
#         ttime=20.0 * units.fs,
#         pfactor=2e6 * units.GPa * (units.fs**2)
#     )
#     md.attach(traj, interval=1000)
#     md.run(steps=20000)
```

**Note: 5 MD tasks are running in parallel and it might cusume a lot of tokens**


```python
# from joblib.parallel import delayed, Parallel


# md_init_dir = Path("md_init_struc")
# md_init_dir.mkdir(exist_ok=True)

# composition_list = [
#     [82, 0, 32, 36], 
#     [78, 8, 32, 36], 
#     [74, 16, 32, 36], 
#     [60, 44, 32, 36], 
#     [36, 92, 32, 36], 
# ]

# md_init_params = []

# for n_Si, n_P, n_Al, n_Na in composition_list:
#     name = f"{n_Si}SiO2-{n_Al//2}Al2O3-{n_Na//2}Na2O-{n_P//2}P2O5"
#     if (md_init_dir / f"{name}.xyz").is_file():
#         atoms = read(md_init_dir / f"{name}.xyz")
#     else:
#         atoms = get_glass(n_Si, n_P, n_Al, n_Na, 2.1)
#         atoms.write(md_init_dir / f"{name}.xyz")
#     md_init_params.append((atoms, name))


# md_dir = Path("pfp_md")
# md_dir.mkdir(exist_ok=True)

# Parallel(n_jobs=5)(delayed(md_simulation)(*params) for params in md_init_params)
```

# 6. Run MD with LightPFP

* Generate the glass with LightPFP. The compositions and MD protocol are the same as above.


```python
import numpy as np
from pathlib import Path
from ase import Atoms
from ase import units
from ase.io import read, Trajectory
from ase.md.langevin import Langevin
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from light_pfp_client.estimator import Estimator
from light_pfp_client.ase_calculator import ASECalculator
from pfcc_extras.liquidgenerator.liquid_generator import LiquidGenerator
from IPython.display import clear_output


def md_simulation(atoms, name):
    MaxwellBoltzmannDistribution(atoms, temperature_K=2000.0)    
    traj = Trajectory(f"light_pfp_md_small/{name}.traj", "w", atoms = atoms)
    md = Langevin(atoms, units.fs, temperature_K=2000.0, friction=0.1)
    md.attach(traj, interval=1000)
    md.run(steps=10000)
    
    md = NPT(
        atoms, 
        units.fs, 
        temperature_K=2000.0, 
        externalstress=units.bar,
        ttime=20.0 * units.fs,
        pfactor=2e6 * units.GPa * (units.fs**2)
    )
    temperature_control = TemperatureControl(
        md, 2000.0, 0.02
    )
    md.attach(temperature_control, interval=10)
    md.attach(traj, interval=1000)
    md.run(steps=75000)
    
    md = NPT(
        atoms, 
        units.fs, 
        temperature_K=500.0, 
        externalstress=units.bar,
        ttime=20.0 * units.fs,
        pfactor=2e6 * units.GPa * (units.fs**2)
    )
    md.attach(traj, interval=1000)
    md.run(steps=20000)
```


```python
md_init_dir = Path("md_init_struc")
md_init_dir.mkdir(exist_ok=True)

composition_list = [
    [82, 0, 32, 36], 
    [78, 8, 32, 36], 
    [74, 16, 32, 36], 
    [60, 44, 32, 36], 
    [36, 92, 32, 36], 
]

md_init_params = []

for n_Si, n_P, n_Al, n_Na in composition_list:
    name = f"{n_Si}SiO2-{n_Al//2}Al2O3-{n_Na//2}Na2O-{n_P//2}P2O5"
    if (md_init_dir / f"{name}.xyz").is_file():
        atoms = read(md_init_dir / f"{name}.xyz")
    else:
        atoms = get_glass(n_Si, n_P, n_Al, n_Na, 2.1)
        atoms.write(md_init_dir / f"{name}.xyz")
    md_init_params.append((atoms, name))
```


```python
md_dir = Path(f"light_pfp_md_small")
md_dir.mkdir(exist_ok=True)

# model_id = "" # The model id of final training
calc = ASECalculator(Estimator(model_id=model_id))

for atoms, name in md_init_params:
    atoms = atoms*(3,3,3)
    atoms.calc = calc
    md_simulation(atoms, name)
```

# 7. Evaluate


```python
from pathlib import Path

eval_dir = Path("evaluate")
eval_dir.mkdir(exist_ok=True)

composition_list = [
    [82, 0, 32, 36], 
    [78, 8, 32, 36], 
    [74, 16, 32, 36], 
    [60, 44, 32, 36], 
    [36, 92, 32, 36], 
]

name_list = [
    f"{n_Si}SiO2-{n_Al//2}Al2O3-{n_Na//2}Na2O-{n_P//2}P2O5" 
    for n_Si, n_P, n_Al, n_Na in composition_list
]

SiO2_P2O5_ratio = [
    (n_Si) / (n_Si + n_P//2)
    for n_Si, n_P, n_Al, n_Na in composition_list
]
```

## Density

* Let's compare the density of glass structures generated by PFP and LightPFP.
* The density is average value of the last 10 ps in the MD simulation.
* My result is as below.

&lt;img src="assets/density.png" width="500"&gt;

* The density agree well with the PFP results:
    * Average error: 0.018g/cm^3
    * Max error: 0.045 g/cm^3


```python
import numpy as np
import matplotlib.pyplot as plt
from ase import units
from ase.io import Trajectory

def get_density(atoms):
    return atoms.get_masses().sum() / units.kg / atoms.get_volume() * 1e27

def get_density_traj(traj, last_n_frames=10):
    return np.mean([get_density(atoms) for atoms in traj[-last_n_frames:]])

```


```python
pfp_density = [get_density_traj(Trajectory(f"pfp_md/{name}.traj")) for name in name_list]
lpfp_density = [get_density_traj(Trajectory(f"light_pfp_md_small/{name}.traj")) for name in name_list]
```


```python
plt.figure()
plt.plot(SiO2_P2O5_ratio, pfp_density, label="PFP", marker="o")
plt.plot(SiO2_P2O5_ratio, lpfp_density, label="LightPFP", marker="o")
plt.xlabel("n_SiO2 / (n_SiO2 + n_P2O5)")
plt.ylabel("density (g/cm3)")
plt.legend()
plt.savefig(eval_dir / "density.png")
```

## RDF

* The radial distribution function of 5 glass structures are calculated and compared.
* My results are as below. The result of LightPFP agree well with PFP.

**Composition 1: 82SiO2-16Al2O3-18Na2O-0P2O5**

&lt;img src="assets/rdf_82SiO2-16Al2O3-18Na2O-0P2O5.png" width="500"&gt;


**Composition 2: 78SiO2-16Al2O3-18Na2O-4P2O5**

&lt;img src="assets/rdf_78SiO2-16Al2O3-18Na2O-4P2O5.png" width="500"&gt;

**Composition 3: 74SiO2-16Al2O3-18Na2O-8P2O5**

&lt;img src="assets/rdf_74SiO2-16Al2O3-18Na2O-8P2O5.png" width="500"&gt;

**Composition 4: 60SiO2-16Al2O3-18Na2O-22P2O5**

&lt;img src="assets/rdf_60SiO2-16Al2O3-18Na2O-22P2O5.png" width="500"&gt;

**Composition 5: 78SiO2-16Al2O3-18Na2O-4P2O5**

&lt;img src="assets/rdf_36SiO2-16Al2O3-18Na2O-46P2O5.png" width="500"&gt;



```python
from light_pfp_evaluate.md import get_rdf, plot_rdf
from ase.io import Trajectory

for name in name_list:
    plot_rdf(
        [500],
        [Trajectory(f"pfp_md/{name}.traj")[-10:]],
        [Trajectory(f"light_pfp_md_small/{name}.traj")[-10:]],
        eval_dir / f"rdf_{name}.png"
    )
```

## Elastic modulus

* The elastic tensor of the glass is calculated with PFP and LightPFP.
* Because it is a little slow, I only calculated one composition, i.e. 74SiO2-16Al2O3-18Na2O-8P2O5.
* My result is

|              |   LightPFP |       PFP |        err |    rel err |
|------------- | ---------- | --------- | ---------- | ---------- |
|C11           |  50.1891   | 51.3725   | 1.18344    | 0.0230365  |
|C12           |  16.9595   | 15.2034   | 1.75608    | 0.115505   |
|C13           |  16.2633   | 15.2597   | 1.00358    | 0.0657666  |
|C14           |   0.689182 |  0.41168  | 0.277501   | 0.674071   |
|C15           |   1.05192  |  0.810074 | 0.24185    | 0.298553   |
|C22           |  49.8913   | 53.1643   | 3.27297    | 0.0615633  |
|C23           |  16.1113   | 15.3702   | 0.741162   | 0.0482208  |
|C24           |   0.906554 |  0.911145 | 0.00459134 | 0.00503909 |
|C25           |   0.532615 |  0.670721 | 0.138106   | 0.205906   |
|C26           |   0.479826 |  0.573146 | 0.09332    | 0.162821   |
|C33           |  48.037    | 52.3248   | 4.28788    | 0.0819473  |
|C35           |   0.670648 |  0.242517 | 0.428131   | 1.76537    |
|C44           |  16.2517   | 18.3673   | 2.11555    | 0.115181   |
|C46           |   0.101012 |  0.22525  | 0.124238   | 0.551555   |
|C55           |  15.9741   | 18.0235   | 2.04932    | 0.113703   |
|C56           |   0.170268 |  0.152172 | 0.0180961  | 0.118918   |
|C66           |  16.359    | 18.5641   | 2.20509    | 0.118783   |
|bulk_modulus  |  27.3836   | 27.5987   | 0.215085   | 0.00779331 |
|shear_modulus |  16.2893   | 18.3781   | 2.08883    | 0.113658   |
|young_modulus |  40.7815   | 45.1193   | 4.33779    | 0.0961405  |
|poisson_ratio |   0.251789 |  0.227527 | 0.0242613  | 0.10663    |



```python
from ase.io import Trajectory
from matlantis_features.utils.calculators.pfp_api_calculator import pfp_estimator_fn
from light_pfp_client.estimator_fn import light_pfp_estimator_fn
from light_pfp_evaluate import evaluate_elastic

# model_id = "" # Please input the model_id of final training
estimator_fn_pfp = pfp_estimator_fn(model_version=model_version, calc_mode=calc_mode)
estimator_fn_light_pfp = light_pfp_estimator_fn(model_id=model_id)

name = "74SiO2-16Al2O3-18Na2O-8P2O5"
atoms = Trajectory(f"pfp_md/{name}.traj")[-1]
evaluate_elastic(atoms, estimator_fn_pfp, estimator_fn_light_pfp, eval_dir / f"elastic_{name}.txt")
```
