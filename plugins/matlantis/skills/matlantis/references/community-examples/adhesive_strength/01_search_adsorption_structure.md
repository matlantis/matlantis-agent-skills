### Copyright Matlantis Corp. as contributors to Matlantis contrib project


##  ポリマーの表面への接着力評価
- シミュレーションは[こちらの論文](https://pubs.acs.org/doi/10.1021/acs.langmuir.5c03183)(K. Kudo and Y. Sumiya, Langmuir (2025))を参考に実装した。
- 本notebookでは表面への分子の吸着構造探索を行う。


```python
import io
import os

from mp_api.client import MPRester
import optuna
from rdkit import Chem

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import LBFGS

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import pfp_api_client
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

from pfcc_extras import show_gui
from pfcc_extras.adsorption.adsorption_structure_search import adstructure_search_for_slab
from pfcc_extras.structure.ase_rdkit_converter import smiles_to_atoms
from pfcc_extras.structure.surface import makesurface


def download_structure(mp_id, api_key):
    """
    Downloads the structure from Materials Project for the given MP-ID.

    Parameters:
        mp_id (str): The Materials Project ID of the structure.
        api_key (str): Your Materials Project API key.

    Returns:
        structure (pymatgen.core.structure.Structure): The pymatgen Structure object.
    """
    # Initialize MPRester with your API key
    
    with MPRester(api_key) as m:
        try:
            # Get the structure by MP-ID
            structure = m.get_structure_by_material_id(mp_id, conventional_unit_cell=True)
            atoms = AseAtomsAdaptor.get_atoms(structure)
            print(f"Successfully downloaded structure for MP-ID: {mp_id}")
            return atoms
        except Exception as e:
            print(f"Error downloading structure for MP-ID {mp_id}: {e}")

def to_conventional(atoms):
    """
    Converts an ASE Atoms object to its conventional cell structure.

    Parameters:
        atoms (ase.atoms.Atoms): The input ASE Atoms object, typically
                                 representing a primitive cell.

    Returns:
        ase.atoms.Atoms: A new ASE Atoms object representing the
                         conventional cell of the input structure.
    """
    structure = AseAtomsAdaptor.get_structure(atoms)
    analyzer = SpacegroupAnalyzer(structure)
    conventional_structure = analyzer.get_conventional_standard_structure()
    return AseAtomsAdaptor.get_atoms(conventional_structure)
```


    



```python
# parameter setting
calc_mode="PBE_PLUS_D3"
model_version="v8.0.0"

api_key = "your_api_key"
mp_id = "mp-81"
```


```python
outpath = f"output/{mp_id}"  
os.makedirs(outpath, exist_ok=True)
```


```python
# create surface
rep_x, rep_y =  3,3
atoms = download_structure(mp_id, api_key)
surface = to_conventional(atoms)
surface = makesurface(surface, miller_indices=(1,1,1), layers=4, rep=[rep_x, rep_y,1], vacuum=23)
surface.set_tags(1)
show_gui(surface, ball_size=1)
```


```python
# set FiaAtoms constraint
thresh = 6.0 # This value was determined by eye from the paper.
constraint = FixAtoms(mask=surface.positions[:, 2] &lt; thresh)
surface.set_constraint(constraint)
write(outpath+"/surface.xyz", surface)
```


```python
# create calculator
estimator = Estimator(calc_mode=calc_mode, model_version=model_version)
calculator = ASECalculator(estimator)

# set calculator and optimize structure
surface.calc = calculator

opt = LBFGS(surface, logfile=None)
opt.run(fmax=0.01)
```


```python
# creare polymer structure (1 unit)
polymer_smiles = "CC(C#N)(C(=O)OCC)CC"
Chem.MolFromSmiles(polymer_smiles)
```


```python
# create polymer atoms and opt the structure
polymer = smiles_to_atoms(polymer_smiles)
polymer.set_tags(2)
polymer.calc = calculator
opt = LBFGS(polymer, logfile=None)
opt.run(fmax=0.01)
```


```python
show_gui(polymer, representations=["ball+stick"])
```


```python
write(outpath+"/polymer.xyz", polymer)
```


```python
# This operation differs from the paper.
# For more accurate results, you can increase the number of trials.
# e.g., n_trials = 500, n_startup_trials = 200
adstructure_search_for_slab(
    calc_mode       = calc_mode,
    model_version   = model_version,
    molec_path      = outpath+"/polymer.xyz",
    slab_path       = outpath+"/surface.xyz",
    z_height_margin = 1.0,    # set margin between slab and molec.
    fix_slab_thickness = thresh,       # set fix slab thickness. 
    TH_max_f        = 5.0,      
    TH_min_f        = 0.0005,   
    tol             = 0.5, 
    rep_x           = rep_x,
    rep_y           = rep_y,
    zgap            = 2.0,    
    fmax            = 0.01,
    optimizer       = LBFGS,
    output_path     = outpath,
    n_trials        = 100,
    njobs           = 20,
    sampler         = optuna.samplers.TPESampler(prior_weight=0.5, n_startup_trials=50),
)
```


```python
# visualize opt history
storage = optuna.storages.JournalStorage(
    optuna.storages.JournalFileStorage(outpath+"/optuna_study.log"),)
study = optuna.load_study(study_name='study',storage=storage)
optuna.visualization.plot_optimization_history(study)
```


```python
# view optimized adsorption structure 
adst = read(io.StringIO(study.best_trial.user_attrs["structure"]), format="json")
show_gui(adst, representations=["ball+stick"])
```


```python
%load_ext watermark
%watermark -n -u -v -iv -w
```

    Last updated: Tue Oct 14 2025
    
    Python implementation: CPython
    Python version       : 3.11.11
    IPython version      : 9.4.0
    
    rdkit         : 2024.9.4
    pfcc_extras   : 0.12.0
    pfp_api_client: 1.24.0
    ase           : 3.25.0
    mp_api        : 0.41.2
    optuna        : 4.5.0
    pymatgen      : 2025.4.10
    
    Watermark: 2.5.0
    

