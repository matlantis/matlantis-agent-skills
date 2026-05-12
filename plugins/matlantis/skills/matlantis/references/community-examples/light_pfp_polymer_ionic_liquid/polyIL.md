Copyright Preferred Networks inc. as contributors to Matlantis contrib project

# Simulation of Polymer Ionic Liquid with LightPFP

In this instance, we have constructed a model using LightPFP to simulate the properties of polymer ionic liquids (PILs). As an example, we utilize poly(ethyl vinyl imidazolium) coupled with PF6 anions. Polymer ILs exhibit a plethora of unique properties and applications. Specifically, poly(ethyl vinyl imidazolium) forms a polymer network with positive charges, while the PF6- ions are distributed discretely within the interstices of the polymer matrix, carrying negative charges. This research followes the article [Molecular Dynamics Simulations of Polymerized Ionic Liquids: Mechanism of Ion Transport with Different Anions](https://pubs.acs.org/doi/10.1021/acsapm.0c00834)

The notebook has 6 components:

1. Initial structures
2. Initial dataset
3. Active learning
4. Post training
5. PFP validation MD
6. LightPFP validation MD


```python
model_version = "v7.0.0"
calc_mode = "pbe_plus_d3"
```

## 1. Initial structures

The target material system is a uniform. It is the polymer material. To enable LightPFP to learn the features of the target material, it is generally necessary to provide multiple initial structures as starting points for data collection (e.g., MD sampling).

(The Advantages of using multiple initial structures: Starting from a single initial structure and running a prolonged MD simulation can make it challenging to sampling different structural characteristics, such as molecular orientation. By starting from different initial structures, a more diverse training dataset can be more readily acquired.)

The `pfcc_extras` package provides a tool, `LiquidGenerator`, to assemble various molecules into solid or liquid bulk structures. The molecular positions and orientations are randomized, which we find particularly suitable for data collection for LightPFP.

We use monomers and shorter polymer chains of poly(ethyl vinyl imidazolium) as inputs for the `LiquidGenerator` because (1) longer polymer chains are more challenging for the `LiquidGenerator` to handle and tend to result in lower-quality initial structures, and (2) shorter polymer chains still contain sufficient structural information and are approximately equivalent for training LightPFP.

Based on these considerations, we have developed the `make_polyIL_structure` function. This function generates atomic structures of mixtures containing monomers, dimers, trimers, 5-unit chains, and 7-unit chains of poly(ethyl vinyl imidazolium) along with PF6 anions. The number of PF6 is determined based on charge balance.


```python
import numpy as np
from pathlib import Path
from ase.io import read
from pfcc_extras.structure.ase_rdkit_converter import atoms_to_smiles
from pfcc_extras.liquidgenerator.liquid_generator import LiquidGenerator
from density import estimate_density
from IPython.display import clear_output

def make_polyIL_structure(mol_type, n_mol):
    assert mol_type in ["monomer", "dimer", "trimer", "polymer", "polymer_x7"]

    monomer = read(f"assets/monomer.xyz")
    polymer = read(f"assets/{mol_type}.xyz")
    anion = read("assets/PF6.xyz")
    if mol_type == "monomer":
        n_anion = n_mol
    elif mol_type == "dimer":
        n_anion = n_mol * 2
    elif mol_type == "trimer":
        n_anion = n_mol * 3
    elif mol_type == "polymer":
        n_anion = n_mol * 5
    else:
        n_anion = n_mol * 7
    
    polymer.positions -= polymer.get_center_of_mass()
    anion.positions -= anion.get_center_of_mass()
    
    # Estimate density for the polymer based on its SMILES conversion.
    density = estimate_density(atoms_to_smiles(monomer))
    
    # Create a mixture of polymer frames and anions.
    composition = [polymer] * n_mol + [anion] * n_anion
    
    # Generate a bulk random structure using the LiquidGenerator.
    liquid_generator = LiquidGenerator(engine="torch", composition=composition, density=density, cubic=True)
    atoms = liquid_generator.run(epochs=100)
    clear_output()
    return atoms


```

By utilizing the make_polyIL_structure function, we create multiple initial structures and save them in the structures directory. These structures will be used for the subsequent collection of the initial dataset.


```python

# Create a directory to store the generated initial structures.
structure_dir = Path("structures")
structure_dir.mkdir(parents=True, exist_ok=True)

initial_structures = []

for n_mol in [6, 12]:    
    filename = f"monomer_{n_mol}.xyz"
    filepath = structure_dir / filename
    if not filepath.is_file():
        atoms = make_polyIL_structure("monomer", n_mol)
        atoms.write(filepath)
    else:
        atoms = read(filepath)
    initial_structures.append(atoms)

for n_mol in [3, 6]:    
    filename = f"dimer_{n_mol}.xyz"
    filepath = structure_dir / filename
    if not filepath.is_file():
        atoms = make_polyIL_structure("dimer", n_mol)
        atoms.write(filepath)
    else:
        atoms = read(filepath)
    initial_structures.append(atoms)

for n_mol in [2, 3]:    
    filename = f"trimer_{n_mol}.xyz"
    filepath = structure_dir / filename
    if not filepath.is_file():
        atoms = make_polyIL_structure("trimer", n_mol)
        atoms.write(filepath)
    else:
        atoms = read(filepath)
    initial_structures.append(atoms)

for n_mol in [2, 3]:    
    filename = f"polymer_{n_mol}.xyz"
    filepath = structure_dir / filename
    if not filepath.is_file():
        atoms = make_polyIL_structure("polymer", n_mol)
        atoms.write(filepath)
    else:
        atoms = read(filepath)
    initial_structures.append(atoms)
```

## 2. Initial dataset

For the creation of the initial dataset, we use the dataset generation methods provided by LightPFP, primarily employing Molecular Dynamics (MD) and rattle techniques.

**1. Molecular Dynamics (MD)**: MD is the fundamental method that provides physically meaningful atomic structures. The MD sampling is conducted using both NVT and NPT ensembles:

  * NVT Ensemble: Utilizes various high temperatures (500K, 1000K, 1500K) to increase structural diversity.

  * NPT Ensemble: Utilizes temperatures (400K, 500K, 600K) similar to the target model to ensure accurate density representation at relevant temperatures.

**2. Rattle**: The rattle method introduces atomic displacement based on a normal distribution to generate more diverse training data, thereby increasing the robustness of the model.

  * Rattle sampling was performed using standard deviations of 0.1 and 0.15 angstroms.

The initial dataset is saved as init_dataset/init.h5.


```python
import numpy as np
from pathlib import Path
from h5py import File
from tqdm.auto import tqdm
from concurrent.futures import as_completed, ThreadPoolExecutor

from pfp_api_client import Estimator, ASECalculator
from light_pfp_data.utils.dataset import H5DatasetWriter
from light_pfp_data.sample.crystal import sample_md, sample_rattle


# Create folder for the initial dataset
init_dataset_dir = Path("init_dataset")
init_dataset_dir.mkdir(parents=True, exist_ok=True)

# Define the initial dataset file
initial_dataset_file = init_dataset_dir / "init.h5"

if initial_dataset_file.exists():
    print(f"Dataset file {initial_dataset_file} already exists. Skipping dataset generation.")
    dataset = H5DatasetWriter(File(initial_dataset_file, "r+"))
else:
    print(f"Creating dataset file {initial_dataset_file} and starting sampling.")
    dataset = H5DatasetWriter(initial_dataset_file)

    # Initialize the PFP estimator and calculator
    estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
    calc = ASECalculator(estimator)

    # List to store our future tasks
    futures = []
    pbar = tqdm(desc="Total progress", total=0, leave=True)

    # Use ThreadPoolExecutor for multithreading sampling tasks
    with ThreadPoolExecutor(max_workers=16) as executor:
        for atoms in initial_structures:
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
            futures += sample_md(
                input_structure=atoms,
                calculator=calc,
                dataset=dataset,
                supercell=(1, 1, 1),
                sampling_temp=[400.0, 500.0, 600.0],
                sampling_pressure=[1.0, 1.0, 1.0],
                sampling_steps=[5000, 5000, 5000],
                sampling_interval=[100, 100, 100],
                ensemble="npt",
                executor=executor,
                pbar=pbar
            )
            futures += sample_rattle(
                input_structure=atoms,
                calculator=calc,
                dataset=dataset,
                stdev=0.1,
                n_sample=10,
                supercell=(1, 1, 1)
            )
            futures += sample_rattle(
                input_structure=atoms,
                calculator=calc,
                dataset=dataset,
                stdev=0.15,
                n_sample=10,
                supercell=(1, 1, 1)
            )

    for f in as_completed(futures):
        _ = f.result()

dataset.h5.close()
```

## 3. Active learning

While the initial dataset comprises a substantial amount of training data, it is still recommended to employ active learning to further collect training data and enhance the quality of the LightPFP model.

**Advantages of Using Active Learning:**
1. It allows for the verification of the stability of the models trained with the initial dataset during MD simulations.
2. In cases of instability, active learning facilitates the collection of structures with significant errors, thereby improving the stability and performance of the LightPFP model.


```python
import logging
from light_pfp_autogen.active_learning import ActiveLearning
from light_pfp_autogen.config import ActiveLearningConfig, TrainConfig, SampleConfig, CommonConfig, MTPConfig

logging.basicConfig(level=logging.INFO)

# Set hyperparameters for the active learning task
active_learning_config = ActiveLearningConfig(
    task_name="polyIL_diffusion",
    pfp_model_version=model_version,
    pfp_calc_mode=calc_mode,
    init_dataset=["init_dataset/init.h5"],
    work_dir="./autogen_workdir",
    training_time=0.5,
    train_config=TrainConfig(
        common_config=CommonConfig(max_forces=50.0, max_energy=5.0),
        mtp_config=MTPConfig(pretrained_model="ORGANIC_SMALL_NN")
    ),
    sample_config=SampleConfig(
        dE_min_coef=3.0,
        dE_max_coef=20.0,
        dF_min_coef=10.0,
        dF_max=50.0,
        dS_min_coef=3.0,
        dS_max_coef=20.0,
        pfp_fallback_samples=5
    )
)

# Initialize the active learning task
active_learning = ActiveLearning(active_learning_config)

# Start the initial training and active learning process
active_learning.initialize()
```

We have developed a short MD script for active learning. This script is used multiple times in each iteration of active learning. The script includes the following steps:

1. Use the `make_polyIL_structure` function to generate a random initial structure of the polymer ionic liquid.
2. Randomly determine the MD temperature, ranging from 300K to 700K.
3. Run a 5 ps NVT MD simulation, attempting to collect training data every 100 MD steps.
4. Run a 50 ps NPT MD simulation, attempting to collect training data every 100 MD steps.

**Note:**
The MD script starts from a new initial structure each time, further enriching the structural diversity of the dataset.

We will perform 10 iterations of active learning. In each iteration, the above procedure is repeated 5 times.


```python
import numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from IPython.display import clear_output
from light_pfp_autogen.context import DataCollectionContext


# Define the MD simulation protocol for active learning iterations tailored to polyIL diffusion task.
def active_learning_protocol(size, steps):
    temperature = np.random.uniform(300, 700)  # K
    atoms = make_polyIL_structure("polymer_x7", 4)
    
    print(f"Running MD for polyIL with size={len(atoms)}, temperature={temperature:.1f} K")
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
    
    # First stage: Short NVT MD using Langevin thermostat for equilibration.
    md_nvt = Langevin(atoms, units.fs, temperature_K=temperature, friction=0.1)
    with DataCollectionContext(md=md_nvt, interval=100, max_samples=20):
        md_nvt.run(steps=5000)  # short equilibration run
    
    # Second stage: Longer NPT MD to generate diverse training data.
    md_npt = NPT(
        atoms,
        units.fs,
        temperature_K=temperature,
        externalstress=1 * units.bar,
        ttime=20.0 * units.fs,
        pfactor=2e6 * units.GPa * (units.fs**2)
    )
    with DataCollectionContext(md=md_npt, interval=100, max_samples=steps // 100 // 2):
        md_npt.run(steps=steps)
    
    clear_output(wait=True)


for i in range(active_learning.iter, 10):
    print(f"Current active learning iteration: {i} (small structures)")
    for _ in range(5):
        active_learning_protocol(size="small", steps=50000)
    active_learning.update()

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
        "pretrained_model": "ORGANIC_SMALL_NN"
    },
}

training_config = TrainConfig.from_dict(
    train_config_dict
)

model_id = submit_training_job(
    training_config,
    active_learning.datasets_list,
    "polyIL_diffusion_final",
)

status = check_training_job_status(model_id)
print(f"Training job {model_id} status: {status}")
```

## 5. PFP validation run

To validate the performance of LightPFP, we need to compare the differences between LightPFP and PFP in actual MD simulations. For this purpose, we will run molecular dynamics simulations using PFP at temperatures of 400K, 500K, and 600K, and save the MD trajectories for future comparative analysis.


```python
import numpy as np
from ase import units
from ase.io import read, Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.md.npt import NPT


def md_protocol(atoms, temperature, steps, traj):
    print(f"Running MD for polyIL with size={len(atoms)}, temperature={temperature:.1f} K")
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)

    traj = Trajectory(traj, "w", atoms=atoms)
    # First stage: Short NVT MD using Langevin thermostat for equilibration.
    md_nvt = Langevin(atoms, units.fs, temperature_K=temperature, friction=0.1)
    md_nvt.attach(traj.write, interval=100)
    md_nvt.run(steps=5000)  # short equilibration run
    
    # Second stage: Longer NPT MD to generate diverse training data.
    md_npt = NPT(
        atoms,
        units.fs,
        temperature_K=temperature,
        mask=np.eye(3),
        externalstress=1 * units.bar,
        ttime=20.0 * units.fs,
        pfactor=2e6 * units.GPa * (units.fs**2)
    )
    md_npt.attach(traj.write, interval=100)
    md_npt.run(steps=steps)
```


```python
from pathlib import Path

md_dir = Path("pfp_md")
md_dir.mkdir(exist_ok=True)

def md_wrap(t):
    calc = ASECalculator(Estimator(model_version=model_version, calc_mode=calc_mode))
    atoms = read("assets/md_init.xyz")
    atoms.calc = calc
    traj = md_dir / f"md_{t}.traj"
    md_protocol(atoms, t, 50000, traj)
```


```python
from joblib import Parallel, delayed

Parallel(n_jobs=3)(delayed(md_wrap)(t) for t in [400, 500, 600])
```

## 6. LightPFP validation run

To further validate LightPFP's performance and compare it directly against PFP, we will repeat the molecular dynamics simulations using LightPFP under the same initial structure, temperature, and MD conditions as the previous PFP simulations. We will also save the trajectories for comparative analysis.


```python
import numpy as np
from ase import units
from ase.io import read, Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.md.npt import NPT


def md_protocol(atoms, temperature, steps, traj):
    print(f"Running MD for polyIL with size={len(atoms)}, temperature={temperature:.1f} K")
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)

    traj = Trajectory(traj, "w", atoms=atoms)
    # First stage: Short NVT MD using Langevin thermostat for equilibration.
    md_nvt = Langevin(atoms, units.fs, temperature_K=temperature, friction=0.1)
    md_nvt.attach(traj.write, interval=100)
    md_nvt.run(steps=5000)  # short equilibration run
    
    # Second stage: Longer NPT MD to generate diverse training data.
    md_npt = NPT(
        atoms,
        units.fs,
        temperature_K=temperature,
        mask=np.eye(3),
        externalstress=1 * units.bar,
        ttime=20.0 * units.fs,
        pfactor=2e6 * units.GPa * (units.fs**2)
    )
    md_npt.attach(traj.write, interval=100)
    md_npt.run(steps=steps)
```


```python
from light_pfp_client import Estimator, ASECalculator

calc = ASECalculator(Estimator(model_id = model_id))
```


```python
from pathlib import Path

md_dir = Path("light_pfp_md")
md_dir.mkdir(exist_ok=True)

for t in [400, 500, 600]:
    atoms = read("assets/md_init.xyz")
    atoms.calc = calc
    traj = md_dir / f"md_{t}.traj"
    md_protocol(atoms, t, 50000, traj)
```

### 6.1 Analysis of results

To evaluate the performance of LightPFP, we compare the MD trajectories obtained from LightPFP against those from PFP in the following four aspects:

#### A. Density

We measure the density at temperatures of 300K, 400K, and 500K. The densities obtained from LightPFP MD simulations align well with the results from PFP simulations as shown in the following figure.

&lt;img src="./assets/density.png" width="500"/&gt;

**Comparision of the density of polymer ionic liquid at 300K, 400K and 500K.**


```python
import numpy as np
import matplotlib.pyplot as plt
from ase import units
from ase.io import Trajectory

def get_density(atoms):
    return atoms.get_masses().sum() / units.kg / atoms.get_volume() * 1e27

def get_density_traj(traj, last_n_frames=100):
    return np.mean([get_density(atoms) for atoms in traj[-last_n_frames:]])
```


```python
pfp_density_400 = [get_density(atoms) for atoms in Trajectory("pfp_md/md_400.traj")]
pfp_density_500 = [get_density(atoms) for atoms in Trajectory("pfp_md/md_500.traj")]
pfp_density_600 = [get_density(atoms) for atoms in Trajectory("pfp_md/md_600.traj")]
lpfp_density_400 = [get_density(atoms) for atoms in Trajectory("light_pfp_md/md_400.traj")]
lpfp_density_500 = [get_density(atoms) for atoms in Trajectory("light_pfp_md/md_500.traj")]
lpfp_density_600 = [get_density(atoms) for atoms in Trajectory("light_pfp_md/md_600.traj")]
```


```python
plt.plot(np.arange(len(pfp_density_400))*0.1, pfp_density_400, label="PFP 400K", c="r")
plt.plot(np.arange(len(lpfp_density_400))*0.1, lpfp_density_400, label="LightPFP 400K", c="r", linestyle="--")
plt.plot(np.arange(len(pfp_density_500))*0.1, pfp_density_500, label="PFP 500K", c="b")
plt.plot(np.arange(len(lpfp_density_500))*0.1, lpfp_density_500, label="LightPFP 500K", c="b", linestyle="--")
plt.plot(np.arange(len(pfp_density_600))*0.1, pfp_density_600, label="PFP 600K", c="g")
plt.plot(np.arange(len(lpfp_density_600))*0.1, lpfp_density_600, label="LightPFP 600K", c="g", linestyle="--")
plt.xlabel("time (ps")
plt.ylabel("density (g/cm^3)")
plt.legend()
plt.savefig("density.png")
```

#### B. RDF

The Radial Distribution Function is analyzed for the last 10 ps of the MD trajectory. At temperatures of 300K, 400K, and 500K, the RDF results from LightPFP MD closely match those from PFP MD simulations.

&lt;img src="./assets/rdf_400.png" width=500&gt;

**Comparison of radial distribution function at 400K**

&lt;img src="./assets/rdf_500.png" width=500&gt;

**Comparison of radial distribution function at 500K**

&lt;img src="./assets/rdf_600.png" width=500&gt;

**Comparison of radial distribution function at 600K**


```python
from light_pfp_evaluate.md import plot_rdf
from ase.io import Trajectory

for t in [400, 500, 600]:
    plot_rdf(
        [t],
        [Trajectory(f"pfp_md/md_{t}.traj")[-100:]],
        [Trajectory(f"light_pfp_md/md_{t}.traj")[-100:]],
        f"rdf_{t}.png"
    )
```

#### C. MSD

To understand the diffusion properties of anions in the MD simulations, we track the positions of the phosphorus atoms (from PF6 anions) and compute the Mean Squared Displacement. The MSD curves for LightPFP are consistent with those obtained from PFP simulations, as illustrated in the figure.

&lt;img src="./assets/msd.png" width=500&gt;

**Comparision of P atom mean squared displacement at 300K, 400K and 500K**


```python
import numpy as np

def get_msd(traj):
    numbers = traj[0].get_atomic_numbers()
    pos = np.array([atoms[numbers==15].get_positions() for atoms in traj])
    msd = [np.mean(np.sum((pos[i+1:] - pos[:-(i+1)])**2, axis=2)) for i in range(len(pos)-1)]
    return msd

```


```python
pfp_msd_400 = get_msd(Trajectory("pfp_md/md_400.traj")[200:])
pfp_msd_500 = get_msd(Trajectory("pfp_md/md_500.traj")[200:])
pfp_msd_600 = get_msd(Trajectory("pfp_md/md_600.traj")[200:])
lpfp_msd_400 = get_msd(Trajectory("light_pfp_md/md_400.traj")[200:])
lpfp_msd_500 = get_msd(Trajectory("light_pfp_md/md_500.traj")[200:])
lpfp_msd_600 = get_msd(Trajectory("light_pfp_md/md_600.traj")[200:])
```


```python
plt.plot(np.arange(len(pfp_msd_400))*0.1, pfp_msd_400, label="PFP 400K", c="r")
plt.plot(np.arange(len(lpfp_msd_400))*0.1, lpfp_msd_400, label="LightPFP 400K", c="r", linestyle="--")
plt.plot(np.arange(len(pfp_msd_500))*0.1, pfp_msd_500, label="PFP 500K", c="b")
plt.plot(np.arange(len(lpfp_msd_500))*0.1, lpfp_msd_500, label="LightPFP 500K", c="b", linestyle="--")
plt.plot(np.arange(len(pfp_msd_600))*0.1, pfp_msd_600, label="PFP 600K", c="g")
plt.plot(np.arange(len(lpfp_msd_600))*0.1, lpfp_msd_600, label="LightPFP 600K", c="g", linestyle="--")
plt.xlabel("time (ps)")
plt.ylabel("msd")
plt.legend()
plt.savefig("msd.png")
```

#### D. Diffusion active energy

By using the MSD data for P atoms, we further calculate the diffusion coefficients. The molecular diffusion coefficient follows an Arrhenius equation with temperature, allowing us to estimate the diffusion activation energy for P atoms (or PF6 anions). The result from LightPFP is 10.46 kJ/mol, which is very close to the PFP result of 10.82 kJ/mol. 

(Note that due to the limited number of data points, the activation energy here primarily serves to demonstrate the consistency between PFP and LightPFP and is not a reliable absolute result.)


&lt;img src="./assets/arrhenius.png" width=500&gt;

**Arrhenius plot of P diffusion coefficient**


```python
def get_diffusion_coef(msd, time_interval):
    time = np.arange(len(msd)) * time_interval
    D = np.polyfit(time, msd, 1)[0] / 6 *1e-5 # cm^2/s
    return D
```


```python
pfp_d_400 = get_diffusion_coef(pfp_msd_400, 100)
pfp_d_500 = get_diffusion_coef(pfp_msd_500, 100)
pfp_d_600 = get_diffusion_coef(pfp_msd_600, 100)
lpfp_d_400 = get_diffusion_coef(lpfp_msd_400, 100)
lpfp_d_500 = get_diffusion_coef(lpfp_msd_500, 100)
lpfp_d_600 = get_diffusion_coef(lpfp_msd_600, 100)
```


```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

R = 8.314
temperatures = np.array([400, 500, 600])
diff_coeffs_pfp = np.array([pfp_d_400, pfp_d_500, pfp_d_600])
diff_coeffs_lpfp = np.array([lpfp_d_400, lpfp_d_500, lpfp_d_600])
inv_temp = 1 / temperatures
ln_diff_pfp = np.log(diff_coeffs_pfp)
ln_diff_lpfp = np.log(diff_coeffs_lpfp)
activation_energy_pfp = -linregress(inv_temp, ln_diff_pfp)[0] * R / 1000
activation_energy_lpfp = -linregress(inv_temp, ln_diff_lpfp)[0] * R / 1000
plt.plot(inv_temp, ln_diff_pfp, marker="o", label=f'PFP {activation_energy_pfp:4.2f}kJ/mol')
plt.plot(inv_temp, ln_diff_lpfp, marker="o", label=f'LightPFP {activation_energy_lpfp:4.2f}kJ/mol')

# Annotate plot
plt.xlabel('1/T (K⁻¹)')
plt.ylabel('ln(D) (ln(m²/s))')
plt.title('Arrhenius Plot')
plt.legend()
plt.grid(True)
    
plt.tight_layout()
plt.savefig(eval_dir/"arrhenius.png")
```
