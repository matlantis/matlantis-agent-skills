Copyright Preferred Networks inc. as contributors to Matlantis contrib project.

# Friction of Fe2O3 surface with lubricant and fatty acid surfactant

In this example, we utilize LightPFP to construct a model that combines metal oxide solids with various lubricants and fatty acid surfactant organic molecules. We employ this model to perform molecular dynamics simulations to observe the behavior of the lubricant and surfactant molecules situated between two solid surfaces during friction. The simulation box consists of a five-layer slab structure: an Fe2O3 slab, a layer of either stearic acid or oleic acid surfactant, a squalane lubricant layer, another layer of either stearic acid or oleic acid surfactant, and another Fe2O3 layer. This configuration results in a complex multilayer structure, making the collection of training data for LightPFP more challenging. This research followes the article [Structure and Friction of Stearic Acid and Oleic Acid Films Adsorbed on Iron Oxide Surfaces in Squalane](https://pubs.acs.org/doi/10.1021/la404024v)

https://www.pure.ed.ac.uk/ws/portalfiles/portal/19397761/langmuir.pdf

The notebook has 5 components:
1. Initial structures
2. Initial dataset
3. Active learning
4. Post training
5. PFP validation MD
6. LightPFP validation MD
7. Large-scale MD with different friction velocities


```python
model_version="v7.0.0"
calc_mode="pbe_plus_d3"
```

## 1. Initial structures

To train a LightPFP model for such a complex system, we need to provide a diverse set of atomic structures that can cover all possible local atomic environments within the aforementioned multilayer structure. Upon analysis, we have identified the following considerations:

a. The multilayer structure is non-uniform, with distinct surfaces and interfaces between each layer.

b. While iron oxide (Fe2O3) has a defined crystal structure and surface configuration, organic molecules do not adhere to a crystalline structure. Their relative positions and angles vary significantly. Therefore, different initial structures should be used to sample a wider range of training data through molecular dynamics simulations.

c. The interfaces between solid and organic molecules, as well as between organic molecules themselves, are crucial.

Based on these assessments, we have implemented four functions:

1. `cut_surface`: This function cuts crystal surfaces from the bulk crystal structure based on the Miller index.

2. `make_bulk_mol`: This function generates the atomic structures of liquid organic molecules, with both the positions and orientations of the molecules being randomized.

3. `make_Fe2O3_mol`: This function creates the interface atomic structures between the solid surface and the liquid organic molecules. The positions and orientations of the organic molecules are randomized.

4. `make_mol_mol_interface`: This function generates the interface structures between two distinct liquid organic molecules. Again, the positions and orientations of the molecules are randomized.


```python
import numpy as np
from ase import units
from ase.build import bulk
from ase.build import diamond111
from ase.cluster import Decahedron
from pfcc_extras.structure.ase_rdkit_converter import atoms_to_smiles, smiles_to_atoms
from pfcc_extras.liquidgenerator.liquid_generator import LiquidGenerator
from density import estimate_density
from ase.io import read
from IPython.display import clear_output
from ase.io import read
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor


def cut_surface(crystal, miller_indices, min_slab_size=10.0, vacuum=20.0):
    # Convert ASE Atoms to Pymatgen Structure
    structure = AseAtomsAdaptor.get_structure(crystal)
    
    # Create slab generator with layer-based thickness
    slabgen = SlabGenerator(
        structure,
        miller_indices,
        min_slab_size=min_slab_size,
        min_vacuum_size=vacuum,
        in_unit_planes=True,
        lll_reduce=True,
    )
    
    # Generate all possible surface terminations
    slabs = slabgen.get_slabs()

    # Return first slab (consider adding termination selection logic if needed)
    slab_atoms = AseAtomsAdaptor.get_atoms(slabs[0])
    del slab_atoms.arrays["bulk_equivalent"]
    del slab_atoms.arrays["bulk_wyckoff"]
    return slab_atoms


def make_bulk_mol(molecule, n_mol=8):
    assert molecule in ["stearic_acid", "oleic_acid", "squalane"]
    # Single molecule structure is loaded from User's file
    if molecule == "stearic_acid":
        smiles = "CCCCCCCCCCCCCCCCCC(=O)O"
    elif molecule == "oleic_acid":
        smiles = "CCCCCCCC\C=C/CCCCCCCC(O)=O"
    else:
        smiles = "CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C"
    mol = smiles_to_atoms(smiles)
    # estimate density of the molecule
    density = estimate_density(smiles)
    # Make a bulk molecule with n_mol molecules
    liquid_generator = LiquidGenerator(engine="torch", composition=[mol] * n_mol, density=density, cubic=True)
    atoms = liquid_generator.run(epochs=100)
    clear_output()
    return atoms


def make_Fe2O3_mol(Fe2O3_surface, molecule, liquid_thickness, slab_interval=2.0, vacuum=50):
    """ Make a Si surface with a liquid layer on top 
    of the surface. The distance between the surface and the liquid layer 
    is `slab_interval`, and there is a vacuum layer with thickness `vacuum`
    on the top of the liquid layer.
    """

    # get the x, y basis vectors of the solid surface
    # so that the liquid layer has the same x, y basis vectors
    cell_surface = Fe2O3_surface.get_cell()[:]
    cell_liquid = cell_surface.copy()
    cell_liquid[2,2] = liquid_thickness
    
    # Load the single molecule structure (provided by user)
    assert molecule in ["stearic_acid", "oleic_acid", "squalane"]
    # Single molecule structure is loaded from User's file
    if molecule == "stearic_acid":
        smiles = "CCCCCCCCCCCCCCCCCC(=O)O"
    elif molecule == "oleic_acid":
        smiles = "CCCCCCCC\C=C/CCCCCCCC(O)=O"
    else:
        smiles = "CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C"
    mol = smiles_to_atoms(smiles)
    density = estimate_density(smiles)

    # Estimate the number of molecules needed to fill the liquid layer
    volume = np.linalg.det(cell_liquid) * 1e-27
    mass = mol.get_masses().sum() / units.kg
    n_mol = int(density * volume / mass)
    
    # Make the liquid layer
    # Use "wall=True" to make the liquid layer with barrier walls on the z-direction boundaries
    liquid_generator = LiquidGenerator(engine="torch", composition=[mol] * n_mol, cell=cell_liquid, wall=True)
    liquid = liquid_generator.run(epochs=100)
    
    # Move the liquid layer to the top of the solid surface
    # Add a margin of slab_interval between the solid surface and the liquid layer
    z_top = Fe2O3_surface.positions[:,2].max()
    liquid.positions[:,2] += z_top + slab_interval
    
    # Set the length of cell in z-direction 
    # to ensure a vacuum layer on the top of the liquid layer
    cell_all = cell_surface.copy()
    cell_all[2,2] = liquid.positions[:,2].max() + vacuum
    Fe2O3_mol = Fe2O3_surface + liquid
    Fe2O3_mol.set_cell(cell_all)
    Fe2O3_mol.set_pbc([True, True, True])
    clear_output()
    return Fe2O3_mol


def make_mol_mol_interface(molecule_1, n1, molecule_2, n2, margin = 1.5):
    if molecule_1 == "stearic_acid":
        smiles = "CCCCCCCCCCCCCCCCCC(=O)O"
    elif molecule_1 == "oleic_acid":
        smiles = "CCCCCCCC\C=C/CCCCCCCC(O)=O"
    else:
        smiles = "CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C"
    mol_1 = smiles_to_atoms(smiles)
    density_1 = estimate_density(smiles)
    
    if molecule_2 == "stearic_acid":
        smiles = "CCCCCCCCCCCCCCCCCC(=O)O"
    elif molecule_2 == "oleic_acid":
        smiles = "CCCCCCCC\C=C/CCCCCCCC(O)=O"
    else:
        smiles = "CC(C)CCCC(C)CCCC(C)CCCCC(C)CCCC(C)CCCC(C)C"
    mol_2 = smiles_to_atoms(smiles)
    density_2 = estimate_density(smiles)
    
    volume_1 = n1 * mol_1.get_masses().sum() / units.kg / density_1 * 10**27
    volume_2 = n2 * mol_2.get_masses().sum() / units.kg / density_2 * 10**27
    volume = volume_1 + volume_2
    x = y = (volume / 2) ** (1/3)
    cell_1 = np.array([[x, 0, 0], [0, y, 0], [0, 0, volume_1 / x / y]])
    liquid_generator = LiquidGenerator(engine="torch", composition=[mol_1]*n1, density=density_1, cell=cell_1, wall=True)
    atoms_1 = liquid_generator.run(epochs=100)
    cell_2 = np.array([[x, 0, 0], [0, y, 0], [0, 0, volume_2 / x / y]])
    liquid_generator = LiquidGenerator(engine="torch", composition=[mol_2]*n2, density=density_2, cell=cell_2, wall=True)
    atoms_2 = liquid_generator.run(epochs=100)
    atoms_2.positions[:,2] += atoms_1.positions[:,2].max() + margin
    z = atoms_2.positions[:,2].max() + margin
    atoms = atoms_1 + atoms_2
    atoms.set_cell([x, y, z])
    atoms.set_pbc([True, True, True])
    clear_output()
    return atoms
```

Based on the functions described above, we generate multiple different initial structures, which are saved in the `structures` folder. These structures will be used for the creation of the initial dataset as outlined below.


```python
from pathlib import Path


structure_dir = Path("structures")
structure_dir.mkdir(parents=True, exist_ok=True)

Fe2O3_crystal_prim = read("assets/Fe2O3.cif")
Fe2O3_crystal_prim.info = {}
del Fe2O3_crystal_prim.arrays["spacegroup_kinds"]

init_structures = []
init_structures.append(Fe2O3_crystal_prim * (3,3,3))
init_structures.append(cut_surface(Fe2O3_crystal_prim, (1, 0, 0), 2, 10.0) * (2, 2, 1))

for molecule in ["stearic_acid", "oleic_acid", "squalane"]:
    for n in [8, 16]:
        if not (structure_dir / f"{molecule}_{n}.xyz").exists():
            atoms = make_bulk_mol(molecule, n)
            atoms.write(structure_dir / f"{molecule}_{n}.xyz")
        else:
            atoms = read(structure_dir / f"{molecule}_{n}.xyz")
        init_structures.append(atoms)

for slab_size in [2, 3]:
    for molecule in ["stearic_acid", "oleic_acid", "squalane"]:
        Fe2O3_100=cut_surface(Fe2O3_crystal_prim, (1, 0, 0), 2, 10.0) * (slab_size, slab_size, 1)
        if not (structure_dir / f"Fe2O3_{slab_size}_{molecule}.xyz").exists():
            atoms = make_Fe2O3_mol(Fe2O3_100, molecule, 15, slab_interval=1.5, vacuum=1.5)
            atoms.write(structure_dir / f"Fe2O3_{slab_size}_{molecule}.xyz", write_info=False)
        else:
            atoms = read(structure_dir / f"Fe2O3_{slab_size}_{molecule}.xyz")
        init_structures.append(atoms)

for molecule in ["stearic_acid", "oleic_acid"]:
    for n in [4, 8]:
        if not (structure_dir / f"squalane_{n}_{molecule}_{n}.xyz").exists():
            atoms = make_mol_mol_interface(molecule, n, "squalane", n)
            atoms.write(structure_dir / f"squalane_{n}_{molecule}_{n}.xyz", write_info=False)
        else:
            atoms = read(structure_dir / f"squalane_{n}_{molecule}_{n}.xyz")
        init_structures.append(atoms)

```

## 2. Initial dataset

For the creation of the initial dataset, we use the dataset generation methods provided by LightPFP, primarily employing Molecular Dynamics (MD) and rattle techniques. 

**1. Molecular Dynamics (MD):**
MD is the fundamental method that provides physically meaningful atomic structures. The MD sampling is conducted using both NVT and NPT ensembles:

- **NVT Ensemble**: Utilizes various high temperatures (500K, 1000K, 1500K) to increase structural diversity.
- **NPT Ensemble**: Utilizes temperatures (300K, 400K, 500K) similar to the target model to ensure accurate density representation at relevant temperatures.

**2. Rattle:**
The rattle method introduces atomic displacement based on a normal distribution to generate more diverse training data, thereby increasing the robustness of the model. 

- Rattle sampling was performed using standard deviations of 0.1 and 0.15 angstroms.

The initial dataset is saved as `init_dataset/init.h5`.


```python
from light_pfp_data.utils.dataset import H5DatasetWriter
from concurrent.futures import as_completed, Future, ThreadPoolExecutor
 
import numpy as np
from h5py import File
from tqdm.auto import tqdm

from pfp_api_client import Estimator, ASECalculator
from light_pfp_data.utils.dataset import H5DatasetWriter
from light_pfp_data.sample.crystal import sample_md, sample_rattle

# Make a folder for initial datasets
init_dataset_dir = Path("init_dataset")
init_dataset_dir.mkdir(parents=True, exist_ok=True)

# Initial dataset file name
initial_dataset = init_dataset_dir / "init.h5"

# Check the initial dataset exist or not
if initial_dataset.exists():
    # If already exist, skip it
    print(f"Dataset file {initial_dataset} already exists. Skip generating initial dataset.")
    dataset = H5DatasetWriter(File(initial_dataset)) 
else:
    # If dataset do not exist, run dataset sampling
    print(f"Dataset file {initial_dataset} is created. Start generating initial structures.")

    # Initialize dataset file
    dataset = H5DatasetWriter(initial_dataset)

    # Initialize PFP estimator and calculator
    estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
    calc = ASECalculator(estimator)
    
    # For Multi-threads running of sampling tasks
    futures = []
    pbar = tqdm(desc="Total progress", total=0, leave=True)
    
    # Sample structures using MD sampling and rattle sampling
    # "initial_structures" contains several LiPF6/EC/DMC mixture atomic structures with various compositions
    with ThreadPoolExecutor(max_workers=16) as executor:
        for atoms in init_structures:
            # Run NVT MD for training structure sampling
            # Temperatures are set as 500K, 1000K and 1500K
            # Run 5000 steps for each temperature
            # Collect 1 training structure every 100 steps
            futures += sample_md(
                input_structure=atoms,
                calculator=calc,
                dataset=dataset,
                supercell=(1, 1, 1),
                sampling_temp=[500.0, 1000.0, 1500.0],
                sampling_steps=[5000, 5000, 5000],
                sampling_interval=[100, 100, 100],
                ensemble="nvt",
                executor=executor,
                pbar=pbar
            )
            # Run NPT MD for training structure sampling
            # Temperatures are set as 300K, 400K and 500K
            # Pressure is set as 1 bar
            # Run 5000 steps for each temperature
            # Collect 1 training structure every 100 steps
            futures += sample_md(
                input_structure=atoms,
                calculator=calc,
                dataset=dataset,
                supercell=(1, 1, 1),
                sampling_temp=[300.0, 400.0, 500.0],
                sampling_pressure=[1.0, 1.0, 1.0],
                sampling_steps=[5000, 5000, 5000],
                sampling_interval=[100, 100, 100],
                ensemble="npt",
                executor=executor,
                pbar=pbar
            )
            # Run random atomic displacement (rattle) to generate training structure
            # Standard deviation is set to 0.1
            # Get 10 random structures
            futures += sample_rattle(
                input_structure=atoms,
                calculator=calc,
                dataset=dataset,
                supercell=(1,1,1),
                stdev=0.1,
                n_sample=10
            )
            # Another random atomic displacement sampling with larger std
            # Standard deviation is set to 0.15
            # Get 10 random structures
            futures += sample_rattle(
                input_structure=atoms,
                calculator=calc,
                dataset=dataset,
                supercell=(1,1,1),
                stdev=0.15,
                n_sample=10
            )

    for f in as_completed(futures):
         _ = f.result()

# Close the dataset after dataset sampling
dataset.h5.close()
```

## 3. Active learning

The initial dataset contains a considerable number of relevant training structures. However, since LightPFP is intended for simulating friction, it's conceivable that molecules will undergo unusual stretching during the simulations. These atomic structures are unlikely to be captured by standard molecular dynamics. Therefore, active learning is essential. Through active learning, we can incorporate molecular dynamics simulations of friction to sample more relevant training structures.

The most crucial aspects are the sampling rules and model update parameters. These need to be defined carefully to ensure that the training data is both relevant and diverse.


```python
import logging
from light_pfp_autogen.active_learning import ActiveLearning
from light_pfp_autogen.config import ActiveLearningConfig, TrainConfig, SampleConfig, CommonConfig, MTPConfig

# To print out the process of active learning
logging.basicConfig(level=logging.INFO)

# Set the parameter of the active learning task
active_learning_config = ActiveLearningConfig(
    task_name = "fe2o3-lubricant-2",
    work_dir = "autogen_workdir_2",
    pfp_model_version = model_version,
    pfp_calc_mode = calc_mode,
    init_dataset = ["init_dataset/init.h5"],
    training_time = 1.0,
    train_config = TrainConfig(
        common_config = CommonConfig(max_forces=50.0),
        mtp_config = MTPConfig(pretrained_model="ORGANIC_SMALL_NN")
    ),
    sample_config = SampleConfig(
        dE_min_coef=3.0,
        dE_max_coef=20.0,
        dF_min_coef=10.0,
        dF_max=50.0,
        dS_min_coef=3.0,
        dS_max_coef=20.0,
    )
)

# Initialize the active learning task
active_learning = ActiveLearning(active_learning_config)

# Start the initial training
active_learning.initialize()
```

To simulate the frictional process, we will define the MD script in the `active_learning_protocol` function:

* Fix the bottom part of the input atomic structure.
* Apply artificial displacements to the top part of the atoms.
* Use NVT MD for the rest of the structure.

We will perform 10 iterations of active learning. Each iteration involves running the above MD script multiple times to ensure that enough training data is collected.

The initial structures for the MD simulations will either be iron oxide-molecule (stearic_acid, oleic_acid, squalane) interfaces or molecule (stearic_acid, oleic_acid)-molecule (squalane) interfaces.
Smaller Initial Structures are used in iterations 1-5, while larger initial structures are used in iterations 6-10.



```python
import numpy as np
from ase import units
from ase.constraints import FixAtoms
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from light_pfp_data.utils.atoms_utils import convert_atoms_to_upper
from light_pfp_autogen.context import DataCollectionContext
from IPython.display import clear_output

# define a MD simulation protocal: Run 5ps NVT + N ps NPT MD with 
# LiPF6/EC/DMC with random composition at 270~450 K
# it is used repeatly in active learning

class ShearDeform:
    def __init__(self, md, ind, velocity=1e-4):
        self.md = md
        self.ind = ind
        self.timestep = self.md.dt / units.fs
        self.velocity = velocity

    def __call__(self):
        v = self.md.atoms.get_cell()[0]
        nv = v / np.linalg.norm(v)
        dr = nv * self.velocity * self.timestep
        self.md.atoms.positions[self.ind] += dr
 

def active_learning_protocol(atoms, steps_npt, steps_shear):
    # set random temperature each MD running to get more diversity
    temperature = np.random.uniform(270, 450)
    velocity = np.random.uniform(1e-4, 1e-3)
    atoms = convert_atoms_to_upper(atoms)
    print(atoms)
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)

    # Initialize NVT MD simulation
    md = Langevin(atoms, units.fs, temperature_K=temperature, friction=0.1)

    # Collect new training structure with DataCollectionContext
    with DataCollectionContext(md=md, interval=100, max_samples=20):
        md.run(steps=5000)

    # Initialize NPT MD simulation
    md = NPT(
        atoms, 
        units.fs, 
        temperature_K=temperature, 
        externalstress=units.bar,
        
        ttime=20.0 * units.fs,
        pfactor=2e6 * units.GPa * (units.fs**2)
    )

    # Collect new training structure with DataCollectionContext
    with DataCollectionContext(md=md, interval=100, max_samples=steps_npt//100//2):
        md.run(steps=steps_npt)

    # Get ready for shear MD
    ## Add 30 A vacuum layer
    cell = atoms.get_cell()[:]
    cell[2,2] += 30
    atoms.set_cell(cell)

    ## Fix some atoms / Shear some atoms
    z = atoms.positions[:,2]
    ind_fix = np.where(z &gt; z.max()-5.0)[0] # get indices for top 5 A
    ind_shear = np.where(z &lt; 5.0)[0] # get indices for bottom 5 A
    ind = np.concatenate([ind_fix, ind_shear])
    atoms.set_constraint(FixAtoms(ind))

    md = Langevin(atoms, units.fs, temperature_K=temperature, friction=0.1)
    shear_deform = ShearDeform(md, ind_shear, velocity)
    md.attach(shear_deform, interval=1)
    with DataCollectionContext(md=md, interval=100, max_samples=steps_shear//100//2):
        md.run(steps=steps_shear)
    
    clear_output()


for i in range(active_learning.iter, 5):
    # First 5 interations, small structure, 15ps NPT MD
    # repeat MD protocol 10 times
    print(f"Current active iteration: {i}")
    for _ in range(6):
        if np.random.random()&lt;0.5:
            # Fe2O3 / molecule interface
            slab_size = np.random.choice([2, 3])
            Fe2O3_100=cut_surface(Fe2O3_crystal_prim, (1, 0, 0), 2, 10.0) * (slab_size, slab_size, 1)
            molecule = np.random.choice(["stearic_acid", "oleic_acid", "squalane"])
            atoms = make_Fe2O3_mol(Fe2O3_100, molecule, 15, slab_interval=1.5, vacuum=1.5)
        else:
            n1 = np.random.randint(1, 5)
            n2 = np.random.randint(1, 5)
            molecule = np.random.choice(["stearic_acid", "oleic_acid"])
            atoms = make_mol_mol_interface(molecule, n1, "squalane", n2)
        active_learning_protocol(atoms, 15000, 20000)
    active_learning.update()


for i in range(active_learning.iter, 10):
    # 5-10 iterations, large structure, 30ps NPT MD
    # repeat MD protocol 5 times
    print(f"Current active iteration: {i}")
    for _ in range(3):
        if np.random.random()&lt;0.5:
            # Fe2O3 / molecule interface
            slab_size = np.random.choice([3, 4])
            Fe2O3_100=cut_surface(Fe2O3_crystal_prim, (1, 0, 0), 2, 10.0) * (slab_size, slab_size, 1)
            molecule = np.random.choice(["stearic_acid", "oleic_acid", "squalane"])
            atoms = make_Fe2O3_mol(Fe2O3_100, molecule, 20, slab_interval=1.5, vacuum=1.5)
        else:
            n1 = np.random.randint(1, 9)
            n2 = np.random.randint(1, 9)
            molecule = np.random.choice(["stearic_acid", "oleic_acid"])
            atoms = make_mol_mol_interface(molecule, n1, "squalane", n2)
        active_learning_protocol(atoms, 15000, 50000)   
    active_learning.update()

```


```python
active_learning.print_md_statistics()
```

## 4. Post training

Due to time constraints during the active learning process, the model updates may not have been fully completed, which could result in a LightPFP model that is not entirely trained. To address this, after the active learning process, we will conduct a final round of model training. This round will utilize all accumulated datasets and will be allocated a longer training duration to ensure that the model is sufficiently trained.


```python
from light_pfp_autogen.utils import submit_training_job, check_training_job_status, estimate_epoch

epoch = estimate_epoch(active_learning.datasets_list, 2)
train_config_dict = {
    "common_config": {
        "total_epoch": epoch,
        "max_forces": 50.0
    },
    "mtp_config": {
        "pretrained_model": "ALL_ELEMENTS_LARGE_7"
    },
}

training_config = TrainConfig.from_dict(
    train_config_dict
)

model_id = submit_training_job(
    training_config,
    active_learning.datasets_list,
    "fe2o3-lubricant-2-final",
)

status = check_training_job_status(model_id)
print(f"Training job {model_id} status: {status}")
```

## 5. PFP validation run

To validate the accuracy of the model, we perform Molecular Dynamics (MD) simulations of the friction process using the PFP (Potential Field Parameters) model. The MD script is similar to the one outlined above. The complete five-layer structure, comprising iron oxide, fatty acid surfactant, and lubricant, will be used for the MD simulation. This initial structure contains 4,604 atoms.

**MD Parameters:**
  - **Friction Velocity:** 30 m/s.
  - **Simulation Duration:** 300 ps (picoseconds).


```python
import time
import numpy as np
from pathlib import Path
from ase import units
from ase.constraints import FixAtoms
from ase.io import read, Trajectory
from ase.md.langevin import Langevin
from pfp_api_client import ASECalculator, Estimator

class PrintDyn:
    def __init__(self, dyn, logfile=None):
        self.dyn = dyn
        self.st = time.perf_counter()
        self.logfile = logfile

    def __call__(self):
        dyn = self.dyn
        msg = (
            f"{dyn.get_number_of_steps(): 6d} {dyn.atoms.get_total_energy():.3f} {dyn.atoms.get_potential_energy():.3f} "
            f"{dyn.atoms.get_masses().sum() / units.kg / dyn.atoms.get_volume() * 1e27:.5f} "
            f"{dyn.atoms.get_temperature():.1f} {time.perf_counter() - self.st:.2f}"
        )

        if self.logfile is not None:
            with open(self.logfile, "a") as fd:
                fd.write(msg+"\n")


class ShearDeform:
    def __init__(self, md, ind, velocity=1e-4):
        self.md = md
        self.ind = ind
        self.timestep = self.md.dt / units.fs
        self.velocity = velocity

    def __call__(self):
        v = self.md.atoms.get_cell()[0]
        nv = v / np.linalg.norm(v)
        dr = nv * self.velocity * self.timestep
        self.md.atoms.positions[self.ind] += dr

```


```python
md_dir = Path("md_pfp")
md_dir.mkdir(parents=True, exist_ok=True)

def md_run(velocity):
    calc = ASECalculator(Estimator(model_version=model_version, calc_mode=calc_mode))
    
    atoms = read("assets/md_init_struc_single.xyz")
    cell = atoms.get_cell()[:]
    cell[2,2] += 50
    atoms.set_cell(cell)
    
    ## Fix some atoms / Shear some atoms
    z = atoms.positions[:,2]
    ind_shear = np.where(z &gt; z.max()-5.0)[0] # get indices for top 5 A
    ind_fix = np.where(z &lt; 5.0)[0] # get indices for bottom 5 A
    ind = np.concatenate([ind_fix, ind_shear])
    atoms.set_constraint(FixAtoms(ind))
    
    atoms.calc = calc

    vel = int(velocity * 1e5)
    traj = Trajectory(md_dir / f"md_{vel}.traj", "w", atoms=atoms)
    
    md = Langevin(atoms, units.fs, temperature_K=300.0, friction=0.1)
    shear_deform = ShearDeform(md, ind_shear, velocity=velocity)
    md.attach(shear_deform, interval=1)
    
    md.attach(traj, interval=1000)
    md.run(steps=300000)

md_run(3e-4)
```

## 6. LightPFP validation run

Run the Light-pfp based MD simulation with the same initial sturcture and MD conditions as above.

**Comparison with PFP Results**

To validate and analyze the effectiveness of the LightPFP model, we compare the results from the final MD frames of the PFP-based simulation with those of the LightPFP model.


1. Molecular Displacement during Friction:

&lt;img src="assets/displacement.png" width="500"&gt;


As expected, the molecules near the moving Fe2O3 layer should show significant displacement.
Both LightPFP and PFP showed similar levels of movement for molecules near the moving layer.

2. Angle between Molecules and Movement Direction

&lt;img src="assets/angle.png" width="500"&gt;

As expected, molecules near the Fe2O3 layer should align their axes with the movement direction, while those further away should be randomly oriented.
The anlge of near Fe2O3 molecules are basically consistent in LightPFP and PFP.



```python
import time
import numpy as np
from pathlib import Path
from ase import units
from ase.constraints import FixAtoms
from ase.io import read, Trajectory
from ase.md.langevin import Langevin
from light_pfp_client import ASECalculator, Estimator


# model_id = "89ogfg0lvv6p8y3j"

class PrintDyn:
    def __init__(self, dyn, logfile=None):
        self.dyn = dyn
        self.st = time.perf_counter()
        self.logfile = logfile

    def __call__(self):
        dyn = self.dyn
        msg = (
            f"{dyn.get_number_of_steps(): 6d} {dyn.atoms.get_total_energy():.3f} {dyn.atoms.get_potential_energy():.3f} "
            f"{dyn.atoms.get_masses().sum() / units.kg / dyn.atoms.get_volume() * 1e27:.5f} "
            f"{dyn.atoms.get_temperature():.1f} {time.perf_counter() - self.st:.2f}"
        )

        if self.logfile is not None:
            with open(self.logfile, "a") as fd:
                fd.write(msg+"\n")


class ShearDeform:
    def __init__(self, md, ind, velocity=1e-4):
        self.md = md
        self.ind = ind
        self.timestep = self.md.dt / units.fs
        self.velocity = velocity

    def __call__(self):
        v = self.md.atoms.get_cell()[0]
        nv = v / np.linalg.norm(v)
        dr = nv * self.velocity * self.timestep
        self.md.atoms.positions[self.ind] += dr

```


```python
md_dir = Path("md_lpfp")
md_dir.mkdir(parents=True, exist_ok=True)
calc = ASECalculator(Estimator(model_id=model_id))

def md_run(velocity):
    atoms = read("assets/md_init_struc_single.xyz")
    cell = atoms.get_cell()[:]
    cell[2,2] += 50
    atoms.set_cell(cell)
    
    ## Fix some atoms / Shear some atoms
    z = atoms.positions[:,2]
    ind_shear = np.where(z &gt; z.max()-5.0)[0] # get indices for top 5 A
    ind_fix = np.where(z &lt; 5.0)[0] # get indices for bottom 5 A
    ind = np.concatenate([ind_fix, ind_shear])
    atoms.set_constraint(FixAtoms(ind))
    
    atoms.calc = calc

    vel = int(velocity * 1e5)
    traj = Trajectory(md_dir / f"md_{vel}.traj", "w", atoms=atoms)
    
    md = Langevin(atoms, units.fs, temperature_K=300.0, friction=0.1)
    shear_deform = ShearDeform(md, ind_shear, velocity=3e-4)
    md.attach(shear_deform, interval=1)
    
    md.attach(traj, interval=1000)
    md.run(steps=300000)

md_run(3e-4)
```

### Comparision of PFP and LightPFP results


```python
import numpy as np
from ase.io import read, Trajectory
import matplotlib.pyplot as plt

atoms_init = read("md_init/md_init_struc_single.xyz")
tags = atoms_init.get_tags()
atoms_lpfp = Trajectory("md_lpfp/md_30.traj")[-1]
atoms_pfp = Trajectory("md_pfp/md_30.traj")[-1]
```


```python
def get_z_pos(atoms_init, atoms_final, tags, ind):
    return (atoms_final[tags==ind].positions[:,2].mean() + atoms_init[tags==ind].positions[:,2].mean()) / 2


def get_displacement(atoms_init, atoms_final, tags, ind):
    axis = atoms_init.get_cell()[0]
    axis_n = axis / np.linalg.norm(axis)
    disp = atoms_final[tags==ind].get_center_of_mass() - atoms_init[tags==ind].get_center_of_mass()
    return np.dot(disp, axis_n)


def get_orientation(atoms_init, atoms_final, tags, ind):
    axis = atoms_init.get_cell()[0]
    axis_n = axis / np.linalg.norm(axis)
    orien = atoms_final[tags==ind].get_moments_of_inertia(vectors=True)[1][0]
    angle = np.arccos(np.dot(axis_n, orien)) * 180 / np.pi
    if angle &gt; 90:
        angle = 180 - angle
    return angle
```


```python
z_lpfp = [get_z_pos(atoms_init, atoms_lpfp, tags, i) for i in range(1, 53)]
z_pfp = [get_z_pos(atoms_init, atoms_pfp, tags, i) for i in range(1, 53)]
displacement_lpfp = [get_displacement(atoms_init, atoms_lpfp, tags, i) for i in range(1, 53)]
displacement_pfp = [get_displacement(atoms_init, atoms_pfp, tags, i) for i in range(1, 53)]
orien_lpfp = [get_orientation(atoms_init, atoms_lpfp, tags, i) for i in range(1, 53)]
orien_pfp = [get_orientation(atoms_init, atoms_pfp, tags, i) for i in range(1, 53)]
```


```python
plt.scatter(z_lpfp, displacement_lpfp, label="LightPFP")
plt.scatter(z_pfp, displacement_pfp, label="PFP")
plt.xlabel("Z position")
plt.ylabel("Displacement along friction direction (A)")
plt.legend()
plt.savefig("displacement.png")
```


```python
plt.scatter(z_lpfp, orien_lpfp, label="LightPFP")
plt.scatter(z_pfp, orien_pfp, label="PFP")
plt.xlabel("Z position")
plt.ylabel("Angle between molecule axis and friction direction")
plt.legend()
plt.savefig("angle.png")
```

## 7. Large-scale MD simulation 

At last, we are able to run long time MD simulation with various friction velocity using a large simulation box.


```python
md_dir = Path("md")
md_dir.mkdir(parents=True, exist_ok=True)
calc = ASECalculator(Estimator(model_id=model_id))

def md_run(velocity):
    atoms = read("assets/md_init_struc_x4.xyz")
    cell = atoms.get_cell()[:]
    cell[2,2] += 50
    atoms.set_cell(cell)
    
    ## Fix some atoms / Shear some atoms
    z = atoms.positions[:,2]
    ind_shear = np.where(z &gt; z.max()-5.0)[0] # get indices for top 5 A
    ind_fix = np.where(z &lt; 5.0)[0] # get indices for bottom 5 A
    ind = np.concatenate([ind_fix, ind_shear])
    atoms.set_constraint(FixAtoms(ind))
    
    atoms.calc = calc

    vel = int(velocity * 1e5)
    traj = Trajectory(md_dir / f"md_{vel}.traj", "w", atoms=atoms)
    
    md = Langevin(atoms, units.fs, temperature_K=300.0, friction=0.1)
    shear_deform = ShearDeform(md, ind_shear, velocity=velocity)
    print_log = PrintDyn(md, logfile=md_dir / f"md_{vel}.log")
    md.attach(shear_deform, interval=1)
    md.attach(print_log, interval=100)
    
    md.attach(traj, interval=100)
    md.run(steps=1000000)
```


```python
for vel in [1e-4, 2e-4, 3e-4, 4e-4, 5e-5]:
    md_run(vel)
```
