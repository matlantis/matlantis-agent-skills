Copyright Preferred Networks inc. and Preferred Computational Chemistry inc. as contributors to Matlantis contrib project.

# LightPFP: dry etching of SiO2 with HF

**Content**

1. Initial structures ' functions
2. Initial dataset
3. Active learning
4. Post training
5. PFP reaction path calculation
6. LightPFP reaction path calculation
7. Real size simulation

Time cost
* LightPFP model generation (section 1 + 2 + 3 + 4): ~ 21 hours
* PFP / LightPFP reference calculation (section  5 + 6) ~ 0.5 hours
* Large-scale simulation (section 7): ~ 45 hours


```python
model_version = "v7.0.0"
calc_mode = "pbe_plus_d3"
```

## Setup


```python
! pip install light-pfp-client==1.0.0 light-pfp-data==1.0.0 light-pfp-evaluate==1.0.0 light-pfp-autogen==0.1.3
```

## 1. Initial structures &amp; functions

* The target materials system contains SiO2 slab, HF etchant and the side products e.g. SiF4.
* Here, we make functions to (a) fill gas molecules into a simulation box, `make_gas` and (b) fill gas molecules on the top of slab.


```python
import random
import numpy as np
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator


def make_gas(mol_list, cell, interval=10.0):
    atoms = Atoms(cell=cell, pbc=True)
    lx = np.linalg.norm(cell[0])
    ly = np.linalg.norm(cell[1])
    lz = np.linalg.norm(cell[2])
    nx = cell[0] / lx
    ny = cell[1] / ly
    nz = cell[2] / lz
    for i in range(int(lx / interval)):
        for j in range(int(ly / interval)):
            for k in range(int(lz / interval)):
                mol = random.choice(mol_list)
                mol_c = mol.copy()
                shift = (i * nx + j * ny + k * nz) * interval
                mol_c.positions += shift.reshape([1,3])
                atoms += mol_c
    return atoms


def make_surface_gas(surface, mol_list, interval=10.0):
    cell = surface.get_cell()[:]
    top_z = surface.positions[:,2].max()
    cell[2,2] -= top_z - 4.0
    gas = make_gas(mol_list, cell, interval)
    gas.positions[:,2] += top_z + 2.0
    surface_gas = surface + gas
    return surface_gas
    
```

* Load SiO2 crystal and surface structures.


```python
from ase.io import read

crystal = read("inputs/SiO2_crystal.cif")
surface = read("inputs/SiO2_surface.xyz")
```

* Make several possible products gas


```python
from ase.build import molecule
from pfcc_extras.structure.ase_rdkit_converter import smiles_to_atoms


product_mol = [
    molecule("HF"),
    molecule("H2"),
    molecule("O2"),
    molecule("F2"),
    molecule("H2O"),
    smiles_to_atoms("F[Si](F)F"),
    smiles_to_atoms("F[Si](F)(F)F"),
    smiles_to_atoms("F[Si](F)(F)(F)F"),
    smiles_to_atoms("F[Si](F)O"),
    smiles_to_atoms("F[Si](F)(F)O"),
    smiles_to_atoms("F[Si](F)(F)(F)O"),
    smiles_to_atoms("F[Si](O)O"),
    smiles_to_atoms("F[Si](F)(O)O"),
    smiles_to_atoms("F[Si](F)(F)(O)O"),
    smiles_to_atoms("F[Si](F)(F)O[Si](F)(F)F"),
]
```


```python
surface_gas = make_surface_gas(surface*(4, 4, 1), product_mol, 10)
```


```python
from ase.visualize import view

view(surface_gas, viewer="ngl")
```

## 2. Initial dataset


```python
from pathlib import Path
from light_pfp_data.utils.dataset import H5DatasetWriter

init_dataset_dir = Path("init_dataset")
init_dataset_dir.mkdir(parents=True, exist_ok=True)
initial_dataset = init_dataset_dir / "init.h5"
```

* Generate the initial dataset from the MD sampling.
* Initial dataset includes: SiO2 crystal, SiO2 slab, gas and SiO2 slab + gas
* time cost: 40 min


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

    # Initialize structures
    input_structures = [
        crystal,
        crystal * (2,2,2),
        crystal * (3,3,3),
        surface,
        surface * (2,2,1),
        surface * (3,3,1),
        make_gas(product_mol, np.eye(3)*20),
        make_gas(product_mol, np.eye(3)*20),
        make_gas(product_mol, np.eye(3)*40),
        make_surface_gas(surface*(2, 2, 1), product_mol),
        make_surface_gas(surface*(3, 3, 1), product_mol),
        make_surface_gas(surface*(4, 4, 1), product_mol),
    ]
    
    from concurrent.futures import as_completed, ThreadPoolExecutor
    from tqdm.auto import tqdm
    
    futures = []
    pbar = tqdm(desc="Total progress", total=0, leave=True)
    
    # Sample structures using various techniques
    with ThreadPoolExecutor(max_workers=16) as executor:
        for atoms in input_structures:
            futures += sample_md(
                input_structure=atoms,
                calculator=calc,
                dataset=dataset,
                supercell=(1, 1, 1),
                sampling_temp=[1000.0, 3000.0, 5000.0],
                sampling_steps=[5000, 5000, 5000],
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
* Run the dry etching script in the active learning framework to collect more training data.
* The dry etching is quite high energy process leads to very disordered structures, thus, we change the sampling criterion to avoid the over sampling. We tried to use a more loose ceriterion by changing the lower bound of energy error range to 5 times of energy MAE, and the lower bound of force error range to 1.5 eV/angstrom.
* To avoid collecting problematic structures (with extremely high energy and forces), we set the `max_forces=50.0` in the train config. It will omit the training structures with maxmium atomic force larger than 50.0 eV/angstrom.


```python
import pathlib
import logging
from light_pfp_autogen.active_learning import ActiveLearning
from light_pfp_autogen.config import ActiveLearningConfig, TrainConfig, SampleConfig, CommonConfig, MTPConfig


logging.basicConfig(level=logging.INFO)

active_learning_config = ActiveLearningConfig(
    task_name = "Si-HF-etching-1",
    work_dir = "./autogen_workdir_1",
    pfp_model_version = model_version,
    pfp_calc_mode = calc_mode,
    init_dataset = ["init_dataset/init.h5"],
    train_config = TrainConfig(
        common_config = CommonConfig(max_forces=50.0),
        mtp_config = MTPConfig(pretrained_model="ALL_ELEMENTS_SMALL_NN_6")
    ),
    sample_config = SampleConfig(
        dE_min_coef = 5.0,
        dE_max_coef = 40.0,
        dF_min = 1.5,
        dF_max = 50.0
    )
)

active_learning = ActiveLearning(active_learning_config)
```


```python
active_learning.initialize()
```

* 0~4 iteration:
  * Simulate the process of 10 HF molecules hitting the surface sequentially, the details are are below:
    * Insert the one HF molecule above the SiO2 surface.
    * Give HF molecule 20~80 eV kinetic energy, and let it collide with the SiO2 surface. The NVE ensemble is used in this step to avoid the unexpected change of initial HF velocity by the thermostate. The timestep is 0.2fs to catch the moment of the impact. 
    * Run NVT MD simulation for 500 fs to cooling down the temperature.
    * Repeat for 10 times to inject 10 HF moleucules.
    * MD steps: 16,000
  * I limited the number of inserted HF molecule to only 10. In the first several iterations, the reliablity of LightPFP model is low, and if the structure becomes unreasonable after first several injection of HF, it is no meaning to collect the following structures.
  * To collect more structures, I repeated the above process for 10 times. The total number of MD steps in one interation is 160,000.


```python
import numpy as np

from ase import units
from ase.io import read
from ase.constraints import FixAtoms
from ase.build import molecule
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from light_pfp_autogen.context import DataCollectionContext


class InsertMol:
    def __init__(self, dyn, incident_energy, init_z):
        self.dyn = dyn
        self.incident_energy = incident_energy
        self.init_z = init_z
        
    def __call__(self):
        atoms = self.dyn.atoms
        if self.dyn.nsteps == 0:
            init_x, init_y = (atoms.cell[0] * np.random.random() + atoms.cell[1] * np.random.random())[:2]
            mol = molecule("HF")
            mol.euler_rotate(np.random.random()*360, np.random.random()*180, np.random.random()*360)
            mol.positions += np.array([[init_x, init_y, self.init_z]])
            vel = (2 * self.incident_energy / mol.get_masses().sum())**0.5
            momenta = np.zeros_like(mol.positions)
            momenta[:, 2] = -mol.get_masses() * vel
            mol.set_momenta(momenta)
            atoms += mol        
        

n_insert = 10
supercell = [2, 2, 1]

for i in range(active_learning.iter, 5):
    for _ in range(10):
        atoms = read("inputs/SiO2_surface.xyz") * supercell
        ind = atoms.positions[:, 2] &lt; 5.0
        fix_atoms = FixAtoms(mask=ind)
        atoms.set_constraint(fix_atoms)
        init_z = 40.0
        MaxwellBoltzmannDistribution(atoms, temperature_K=300.0)
        md0 = NVTBerendsen(atoms, units.fs, temperature_K=300.0, taut=500 * units.fs)
        with DataCollectionContext(md=md0, interval=100, max_samples=5):
            md0.run(1000)
        for i in range(n_insert):
            kinetic_energy = np.random.uniform(20, 80)
            # Hit the surface
            # Time step 0.2 fs
            md1 = VelocityVerlet(atoms, 0.2 * units.fs)
            insert_mol = InsertMol(md1, kinetic_energy, init_z)
            md1.attach(insert_mol, interval=1)
            with DataCollectionContext(md=md1, interval=25, max_samples=20):
                md1.run(1000)
            # Relax
            md2 = NVTBerendsen(atoms, units.fs, temperature_K=300.0, taut=100 * units.fs)
            with DataCollectionContext(md=md2, interval=50, max_samples=5):
                md2.run(500)

    active_learning.update()

```

* 5~9 iterations
  * It is same as above except the number of injected molecule is increased from 10 to 20.
* Total MD steps: 310,000


```python
n_insert = 20
supercell = [2, 2, 1]

for i in range(active_learning.iter, 10):
    for _ in range(10):
        atoms = read("inputs/SiO2_surface.xyz") * supercell
        ind = atoms.positions[:, 2] &lt; 5.0
        fix_atoms = FixAtoms(mask=ind)
        atoms.set_constraint(fix_atoms)
        init_z = 40.0
        MaxwellBoltzmannDistribution(atoms, temperature_K=300.0)
        md0 = NVTBerendsen(atoms, units.fs, temperature_K=300.0, taut=500 * units.fs)
        with DataCollectionContext(md=md0, interval=100, max_samples=5):
            md0.run(1000)
        for i in range(n_insert):
            kinetic_energy = np.random.uniform(20, 80)
            # Hit the surface
            # Time step 0.2 fs
            md1 = VelocityVerlet(atoms, 0.2 * units.fs)
            insert_mol = InsertMol(md1, kinetic_energy, init_z)
            md1.attach(insert_mol, interval=1)
            with DataCollectionContext(md=md1, interval=25, max_samples=20):
                md1.run(1000)
            # Relax
            md2 = NVTBerendsen(atoms, units.fs, temperature_K=300.0, taut=100 * units.fs)
            with DataCollectionContext(md=md2, interval=50, max_samples=5):
                md2.run(500)

    active_learning.update()
```

* 10~14 iterations
  * It is same as above except:
    * The number of injected molecule is increased to 40.
    * The size of SiO2 slab is increased.
* Total MD steps: 305,000


```python
n_insert = 40
supercell = [4, 4, 1]

for i in range(active_learning.iter, 15):
    for _ in range(5):
        atoms = read("inputs/SiO2_surface.xyz") * supercell
        ind = atoms.positions[:, 2] &lt; 5.0
        fix_atoms = FixAtoms(mask=ind)
        atoms.set_constraint(fix_atoms)
        init_z = 40.0
        MaxwellBoltzmannDistribution(atoms, temperature_K=300.0)
        md0 = NVTBerendsen(atoms, units.fs, temperature_K=300.0, taut=500 * units.fs)
        with DataCollectionContext(md=md0, interval=100, max_samples=5):
            md0.run(1000)
        for i in range(n_insert):
            kinetic_energy = np.random.uniform(20, 80)
            # Hit the surface
            # Time step 0.2 fs
            md1 = VelocityVerlet(atoms, 0.2 * units.fs)
            insert_mol = InsertMol(md1, kinetic_energy, init_z)
            md1.attach(insert_mol, interval=1)
            with DataCollectionContext(md=md1, interval=25, max_samples=20):
                md1.run(1000)
            # Relax
            md2 = NVTBerendsen(atoms, units.fs, temperature_K=300.0, taut=100 * units.fs)
            with DataCollectionContext(md=md2, interval=50, max_samples=5):
                md2.run(500)

    active_learning.update()
```


```python
active_learning.print_md_statistics()
```

## 4. Post training 
* Train a final LightPFP model with all the dataset collected in active learing.


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

training_config = TrainConfig.from_dict(
    train_config_dict
)

model_id = submit_training_job(
    training_config,
    active_learning.datasets_list,
    "SiO2-HF-etching-1-final",
)

status = check_training_job_status(model_id)
print(f"Training job {model_id} status: {status}")

# model_id = "jm6mjc5ju4p2qge2"
```

## 5. PFP reaction path calculation

* Use NEB &amp; CINEB method to calculate a chemical reaction of HF and SiO2 cluster. The PFP calculator is used.
* Considering the convergency of NEB is difficult to achieve, we reduce the `fmax` step by step.
* The energy profile of PFP calculation result:

&lt;img src="assets/neb_pfp_19images_fmax0.07.png" width="500"&gt;

* The reaction path is:

&lt;img src="assets/neb_pfp_19images_fmax0.07.gif" width="300"&gt;


```python
from pathlib import Path
from ase import Atoms, units
from ase.io import read, write
from ase.optimize import LBFGS
from ase.data import atomic_numbers, chemical_symbols
from ase.mep import NEB, idpp_interpolate, interpolate, NEBTools

import matplotlib.pyplot as plt

from matlantis_features.utils.calculators import get_calculator, pfp_estimator_fn
from matlantis_features.features.common import FireLBFGSASEOptFeature
from pfp_api_client import Estimator, EstimatorCalcMode, ASECalculator
from pfcc_extras import view_ngl, SurfaceEditor
from pfcc_extras.visualize.povray import traj_to_gif

from IPython.display import Image

estimator_fn = pfp_estimator_fn(model_version='v7.0.0', calc_mode='pbe_plus_d3')
```


```python
inputdir = Path('inputs')
outputdir = Path('./output-pfp')
for _dir in (inputdir, outputdir):
    _dir.mkdir(exist_ok=True)
```


```python
IS = read(inputdir/'IS.cif')
FS = read(inputdir/'FS.cif')
opt = FireLBFGSASEOptFeature(fmax=0.02, estimator_fn=estimator_fn)
IS = opt(IS).atoms.ase_atoms
FS = opt(FS).atoms.ase_atoms
```


```python
num_images = 7
images = [IS.copy() for i in range(num_images-1)] + [FS]
idpp_interpolate(images)
neb = NEB(images, k=5.0)
for atoms in neb.images:
    atoms.calc = get_calculator(estimator_fn)
    
opt = LBFGS(neb)
opt.run(fmax=0.2, steps=1000)
display(view_ngl(neb.images))
nebt = NEBTools(neb.images)
view = nebt.plot_band()
```


```python
neb = NEB(images_new, k=5.0)
for atoms in neb.images:
    atoms.calc = get_calculator(estimator_fn)
    
opt = LBFGS(neb)
opt.run(fmax=0.1, steps=1000)
display(view_ngl(neb.images, ['ball+stick'], replace_structure=True))
nebt = NEBTools(neb.images)
view = nebt.plot_band()
```


```python
fmax = 0.07
k = 5.0
neb = NEB(images_new, k=k)
for atoms in neb.images:
    atoms.calc = get_calculator(estimator_fn)
    
opt = LBFGS(neb)
opt.run(fmax=fmax, steps=1000)
display(view_ngl(neb.images, ['ball+stick'], replace_structure=True))
nebt = NEBTools(neb.images)
view = nebt.plot_band()
```


```python
images = [image.copy() for image in neb.images]
neb_ci = NEB(images, k=k, climb=True)
for atoms in neb_ci.images:
    atoms.calc = get_calculator(estimator_fn)
    
opt = LBFGS(neb_ci)
opt.run(fmax=fmax, steps=2000)
display(view_ngl(neb_ci.images, ['ball+stick'], replace_structure=True))
nebt = NEBTools(neb_ci.images)
view = nebt.plot_band()
```


```python
basename = outputdir/f'neb_{len(neb_ci.images)}images_fmax{fmax}'
trajfile = f'{basename}.traj'
print(trajfile)
write(trajfile, neb_ci.images)
```


```python
pngfile = f'{basename}.png'
print(pngfile)
fig, ax = plt.subplots(figsize=(8,6))
NEBTools(neb_ci.images).plot_band(ax=ax)
fig.savefig(pngfile)
```


```python
giffile = f'{basename}.gif'
images = read(trajfile, index=':')
for a in images:
    a.rotate(140, 'z')
#    a.rotate(-30, 'x')
    a.rotate(100, 'y')

traj_to_gif(
    images,
    gif_filepath=giffile,
    povdir=outputdir / "pov",
    pngdir=outputdir / "png",
    clean=False,
)
Image(giffile)
```

## 6. LightPFP reaction path calculation
* The reaction path is calculated with the same method as above using the LightPFP model.
* The energy profile is

&lt;img src="assets/neb_lpfp_19images_fmax0.07.png" width="500"&gt;

* The reaction path is

&lt;img src="assets/neb_lpfp_19images_fmax0.07.gif" width="300"&gt;

* Comparison of reaction barrier energy between PFP and LightPFP

| | Ef | Er |
|-|-|-|
| PFP | 1.029 | 1.561 |
| LightPFP | 0.844 | 1.560 |

* LightPFP model basically reproduced the same reaction pathway with similar barrier energies, even only MD frames are used for training.


```python
from pathlib import Path
from ase import Atoms, units
from ase.io import read, write
from ase.optimize import LBFGS, FIRE
from ase.data import atomic_numbers, chemical_symbols
from ase.mep import NEB, idpp_interpolate, interpolate, NEBTools

import matplotlib.pyplot as plt

from matlantis_features.utils.calculators import get_calculator
from matlantis_features.features.common import FireLBFGSASEOptFeature
from light_pfp_client.ase_calculator import ASECalculator
from light_pfp_client.estimator import Estimator
from light_pfp_client.estimator_fn import light_pfp_estimator_fn
from pfcc_extras import view_ngl, SurfaceEditor
from pfcc_extras.visualize.povray import traj_to_gif

from IPython.display import Image

estimator_fn = light_pfp_estimator_fn(model_id=model_id)
```


```python
inputdir = Path('inputs')
outputdir = Path('./output-lpfp')
for _dir in (inputdir, outputdir):
    _dir.mkdir(exist_ok=True)
```


```python
IS = read(inputdir/'IS.cif')
FS = read(inputdir/'FS.cif')
opt = FireLBFGSASEOptFeature(fmax=0.02, estimator_fn=estimator_fn)
IS = opt(IS).atoms.ase_atoms
FS = opt(FS).atoms.ase_atoms
```


```python
num_images = 7
images = [IS.copy() for i in range(num_images-1)] + [FS]
idpp_interpolate(images)
neb = NEB(images, k=5.0)
for atoms in neb.images:
    atoms.calc = get_calculator(estimator_fn)

opt = FIRE(neb)
opt.run(fmax=0.2, steps=500)
display(view_ngl(neb.images))
nebt = NEBTools(neb.images)
view = nebt.plot_band()
```


```python
neb = NEB(images_new, k=5.0)
for atoms in neb.images:
    atoms.calc = get_calculator(estimator_fn)
    
opt = LBFGS(neb)
opt.run(fmax=0.1, steps=1000)
display(view_ngl(neb.images, ['ball+stick'], replace_structure=True))
nebt = NEBTools(neb.images)
view = nebt.plot_band()
```


```python
fmax = 0.07
k = 5.0
neb = NEB(images_new, k=k)
for atoms in neb.images:
    atoms.calc = get_calculator(estimator_fn)
    
opt = LBFGS(neb)
opt.run(fmax=fmax, steps=1000)
display(view_ngl(neb.images, ['ball+stick'], replace_structure=True))
nebt = NEBTools(neb.images)
view = nebt.plot_band()
```


```python
images = [image.copy() for image in neb.images]
neb_ci = NEB(images, k=k, climb=True)
for atoms in neb_ci.images:
    atoms.calc = get_calculator(estimator_fn)
    
opt = LBFGS(neb_ci)
opt.run(fmax=fmax, steps=2000)
display(view_ngl(neb_ci.images, ['ball+stick'], replace_structure=True))
nebt = NEBTools(neb_ci.images)
view = nebt.plot_band()
```


```python
basename = outputdir/f'neb_{len(neb_ci.images)}images_fmax{fmax}'
trajfile = f'{basename}.traj'
print(trajfile)
write(trajfile, neb_ci.images)
```


```python
pngfile = f'{basename}.png'
print(pngfile)
fig, ax = plt.subplots(figsize=(8,6))
NEBTools(neb_ci.images).plot_band(ax=ax)
fig.savefig(pngfile)
```


```python
giffile = f'{basename}.gif'
images = read(trajfile, index=':')
for a in images:
    a.rotate(140, 'z')
#    a.rotate(-30, 'x')
    a.rotate(100, 'y')

traj_to_gif(
    images,
    gif_filepath=giffile,
    povdir=outputdir / "pov",
    pngdir=outputdir / "png",
    clean=False,
)
Image(giffile)
```

## 6. Real-size etching simulation

* I use the LightPFP model for simulate the etching process on a large SiO2 slab.
* SiO2 slab size: ~ 10nm x 10nm x 10nm, with 72000 atoms.
* 1000 HF molecule is injected. The kinetic energy is 40 eV.
* HF molecules etch only the central area (2nm x 2nm) of slab.
* The simulation cost about 45 hours.

* Here, the slice of SiO2 slab after etching is presented. A clear etching hole with depth more than 2 nm is observed.

**Depth 0~1 nm**

&lt;img src="assets/top_1.png" width="300"&gt;

**Depth 1~2 nm**

&lt;img src="assets/top_2.png" width="300"&gt;

**Depth 2~3 nm**

&lt;img src="assets/top_3.png" width="300"&gt;




```python
import os
import numpy as np
from time import perf_counter

from ase import units
from ase.io import read
from ase.constraints import FixAtoms
from ase.build import molecule
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from light_pfp_client.estimator import Estimator
from light_pfp_client.ase_calculator import ASECalculator
```


```python
# model_id = ""  # Please use the last LightPFP model generated from active learning
estimator = Estimator(model_id=model_id, book_keeping=True)
calc = ASECalculator(estimator)
```


```python
class InsertMol:
    def __init__(self, dyn, incident_energy, init_z, xrange, yrange):
        self.dyn = dyn
        self.incident_energy = incident_energy
        self.init_z = init_z
        self.xrange = xrange
        self.yrange = yrange
    
    def _init_pos(self):
        init_x = self.xrange[0] + np.random.random() * (self.xrange[1] - self.xrange[0])
        init_y = self.yrange[0] + np.random.random() * (self.yrange[1] - self.yrange[0])
        return init_x, init_y
        
    def __call__(self):
        atoms = self.dyn.atoms
        if self.dyn.nsteps == 0:
            init_x, init_y = self._init_pos()
            mol = molecule("HF")
            mol.euler_rotate(np.random.random()*360, np.random.random()*180, np.random.random()*360)
            mol.positions += np.array([[init_x, init_y, self.init_z]])
            vel = (2 * self.incident_energy / mol.get_masses().sum())**0.5
            momenta = np.zeros_like(mol.positions)
            momenta[:, 2] = -mol.get_masses() * vel
            mol.set_momenta(momenta)
            atoms += mol        
        
```


```python
# From intial state
n_init = -1
atoms = read("inputs/SiO2_surface_20x20x20.xyz")

atoms.calc = calc
```


```python
n_iter = 1000
kinetic_energy = 40.0
init_z = 150.0
xrange = [15, 35]
yrange = [35, 55]
```


```python
snapshot_dir = Path("snapshot")
snapshot_dir.mkdir(exist_ok=True)
```


```python
class PrintDyn:
    def __init__(self, dyn, logfile=None):
        self.dyn = dyn
        self.st = perf_counter()
        self.logfile = logfile
        if self.logfile is not None:
            with open(self.logfile, "w") as fd:
                fd.write("# step E_tot E_pot density T elapsed_time\n")
    def __call__(self):
        dyn = self.dyn
        msg = (
            f"{dyn.get_number_of_steps(): 6d} {atoms.get_total_energy():.3f} {atoms.get_potential_energy():.3f} "
            f"{atoms.get_masses().sum() / units.kg / atoms.get_volume() * 1e27:.5f} "
            f"{atoms.get_temperature():.1f} {perf_counter() - self.st:.2f}"
        )
        print(msg)
        if self.logfile is not None:
            with open(self.logfile, "a") as fd:
                fd.write(msg+"\n")
```


```python
logger = PrintDyn(None, snapshot_dir / "md.log")

# Fix the bottom 10 A
ind = atoms.positions[:, 2] &lt; 10.0
fix_atoms = FixAtoms(mask=ind)
atoms.set_constraint(fix_atoms)

# Initial MD simulation
MaxwellBoltzmannDistribution(atoms, temperature_K=300.0)
md0 = NVTBerendsen(atoms, units.fs, temperature_K=300.0, taut=500 * units.fs)
logger.dyn = md0
md0.attach(logger, interval=100)
md0.run(steps = 5000)

for n in range(n_init+1, n_iter):
    md1 = VelocityVerlet(atoms, 0.5 * units.fs)
    insert_mol = InsertMol(md1, kinetic_energy, init_z, xrange, yrange)
    md1.attach(insert_mol, interval=1)
    logger.dyn = md1
    md1.attach(logger, interval=100)
    md1.run(400)
    md2 = NVTBerendsen(atoms, units.fs, temperature_K=300.0, taut=500 * units.fs)
    logger.dyn = md2
    md2.attach(logger, interval=100)
    md2.run(300)
    atoms.wrap()
    atoms.write(snapshot_dir / f"snapshot_{n}.xyz")
```


```python
init = read("inputs/SiO2_surface_20x20x20.xyz")
final = read("snapshot/snapshot_999.xyz")
z_top = init.positions[:,2].max()
```


```python
from ase.visualize import view

ind = (final.positions[:,2] &lt; z_top+5.0) &amp; (final.positions[:,2] &gt; z_top-10.0)
view(final[ind], viewer="ngl")
```


```python
ind = (final.positions[:,2] &lt; z_top-10.0) &amp; (final.positions[:,2] &gt; z_top-20.0)
view(final[ind], viewer="ngl")
```


```python
ind = (final.positions[:,2] &lt; z_top-20.0) &amp; (final.positions[:,2] &gt; z_top-30.0)
view(final[ind], viewer="ngl")
```
