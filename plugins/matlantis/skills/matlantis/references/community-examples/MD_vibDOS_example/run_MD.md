Copyright ENEOS Corporation as contributors to Matlantis contrib project.


```python
import os
import sys
import glob
import shutil

from pathlib import Path
from io import StringIO
from time import perf_counter

import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ase import Atoms, units
from ase.io import read, write, Trajectory
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
from ase.md.npt import NPT
from ase.md.nptberendsen import NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.logger import MDLogger

from IPython.display import display
```


    



```python
from pfcc_extras.visualize.view import view_ngl
from pfcc_extras.structure.ase_rdkit_converter import smiles_to_atoms
from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
```

    /home/jovyan/.local/lib/python3.7/site-packages/pfp_api_client/__init__.py:36: UserWarning: New version of pfp-api-client is available. Please consider upgrading by `pip install -U pfp-api-client`.
      f"New version of {package_name} is available. Please consider"



```python
MODEL_VERSION='v3.0.0'
CALC_MODE=EstimatorCalcMode.PBE_U_PLUS_D3

def get_calculator(model_version=MODEL_VERSION, calc_mode=CALC_MODE):
    estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
    calculator = ASECalculator(estimator)
    return calculator
    
```


```python
# --- System settings
molecule_name = 'DMSO'
nmols = 50
name_smiles = {'DMSO':'CS(=O)C'}
bulk_name = 'Si'

repeat_mol = [3,3,3]
repeat_solid = [3,3,3]

# --- Initial geometry settings
fmax_opt = 0.03
steps_opt = 200
fmax_isolated = 0.01
fmax_liquid = 0.1
mol_cell = [5.0, 5.0, 5.0]

# --- MD settings ----
# Please change based on configuration
target_temp = 300
target_pressure = 101325 * units.Pascal
#steps_equilib = 1000000
#steps_product = 1000000
steps_equilib = 10000
steps_product = 10000
timestep_solid = 2.0
timestep_mol = 1.0
temperature = 300
log_interval = 10
traj_interval = 1
# --- Parameters
mass_density_factor = units._amu * 1e30 
huge = 999.0
# ------------------
```


```python
inputdir = './input'
outputdir = './output'
os.makedirs(outputdir, exist_ok=True)
```


```python
solid_initial = read(os.path.join(inputdir, 'solid_initial.traj'))
mol_initial = read(os.path.join(inputdir, 'mol_initial.traj'))

mass_molecule = mol_initial.get_masses().sum()
print(f'Number of atoms in the solid unit cell ({solid_initial}): {len(solid_initial)}')
print(f'Number of atoms per molecule ({mol_initial}): {len(mol_initial)}')
```

    Number of atoms in the solid unit cell (Atoms(symbols='Si8', pbc=True, cell=[5.420735127923755, 5.420735031387806, 5.420735088761772], calculator=SinglePointCalculator(...))): 8
    Number of atoms per molecule (Atoms(symbols='CSOCH6', pbc=True, cell=[5.0, 5.0, 5.0], calculator=SinglePointCalculator(...))): 10



```python
v = view_ngl([mol_initial, solid_initial], ['ball+stick'], h=300,replace_structure=True)
display(v)
```


    HBox(children=(NGLWidget(max_frame=1), VBox(children=(Dropdown(description='Show', options=('All', 'H', 'C', '…



```python
mol = mol_initial * repeat_mol
solid = solid_initial * repeat_solid
```


```python
def run_md(atoms,
    time_step    = 2.0,    # fsec
    temperature = 300,    # Kelvin
    num_md_steps = 100000, 
    log_interval = 10,
    traj_interval = 1,
    sigma   = 1.0,     # External pressure in bar
    ttime   = 20.0,    # Time constant in fs
    pfactor = 2e6,     # Barostat parameter in GPa
    prefix = 'md',
    outputdir = './outout',
    print_md=False,
    mask = None,
    print_interval=100,
    fix_shellshape=False,
     ):

    
    # run MD
    print("temperature = ",temperature)
    print(f"sigma = {sigma:.1e} bar")
    print(f"ttime = {ttime:.3f} fs")
    print(f"pfactor = {pfactor:.3e} GPa*fs^2")

    output_filename = os.path.join(outputdir, f"{prefix}_{temperature}K")
    log_filename = output_filename + ".log"
    traj_filename = output_filename + ".traj"

    print("log_filename = ",log_filename)
    print("traj_filename = ",traj_filename)

    atoms = atoms.copy()
    atoms.calc = get_calculator()

    # Set the momenta corresponding to T=300K
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature,force_temp=True)
    Stationary(atoms)

    dyn = NPT(atoms,
          time_step*units.fs,
          temperature_K = temperature,
          externalstress = sigma*units.bar,
          ttime = ttime*units.fs,
          pfactor = pfactor*units.GPa*(units.fs**2),
          trajectory = traj_filename,
          # loginterval argument affects both traj and log interval. Thus, I attach MDLogger later.
#          logfile = log_filename,
#          loginterval=log_interval, 
          logfile = None,
          loginterval=traj_interval,
          mask=mask,
          )
    if fix_shellshape:
        dyn.set_fraction_traceless(0)
        
    # Print statements
    # atttch logger with less/more frequent log_interval than traj_interval
    dyn.attach(MDLogger(dyn, atoms, log_filename, header=True, stress=False, peratom=True, mode="a"), interval=log_interval)
    if print_md:
        def print_dyn():
            imd = dyn.get_number_of_steps()
            etot  = atoms.get_total_energy()
            temp_K = atoms.get_temperature()
            stress = atoms.get_stress(include_ideal_gas=True)/units.GPa
            volume = atoms.get_volume()
            lx, ly, lz, alpha, beta, gamma = atoms.cell.cellpar()

            stress_ave = (stress[0]+stress[1]+stress[2])/3.0 
            elapsed_time = perf_counter() - start_time
            #print(f"  {imd: &gt;3}   {etot:.3f}    {temp_K:.2f}    {stress_ave:.2f}  {stress[0]:.2f}  {stress[1]:.2f}  {stress[2]:.2f}  {stress[3]:.2f}  {stress[4]:.2f}  {stress[5]:.2f}    {elapsed_time:.3f}")
            print(f"  {imd: &gt;3}   {etot:.3f}    {temp_K:.2f} {volume:.2f} {lx:.2f} {ly:.2f} {lz:.2f} {alpha} {beta} {gamma} {elapsed_time:.3f}")
        dyn.attach(print_dyn, interval=print_interval)

    # Run the dynamics
    start_time = perf_counter()
    dyn.run(num_md_steps)

    return atoms, dyn
```


```python
t0 = perf_counter()

timestep = timestep_mol
temperature = target_temp
mask = [1, 1, 1]
print_interval = 10

prefix = 'mol'
prefix_equil = f'{prefix}_equil{int(timestep*steps_equilib/1000)}ps'
prefix_product = f'{prefix}_product{int(timestep*steps_equilib/1000)}ps'

# Equilibration
mol, dyn = run_md(mol, time_step=timestep, temperature=temperature, 
                  num_md_steps=steps_equilib, prefix=prefix_equil, 
                  log_interval=log_interval, traj_interval=traj_interval,
                  outputdir=outputdir, fix_shellshape=True, mask=mask)
t1 = perf_counter()
print(f'Elapsed time (equilib): {t1-t0}')

# Product run
mol, dyn = run_md(mol, time_step=timestep, temperature=temperature, 
                  num_md_steps=steps_product, prefix=prefix_product, 
                  log_interval=log_interval, traj_interval=traj_interval,
                  outputdir=outputdir, fix_shellshape=True, mask=mask)
t2 = perf_counter()
print(f'Elapsed time (product): {t2-t1}')
print(f'Elapsed time (total)  : {t2-t0}')
```

    temperature =  300
    sigma = 1.0e+00 bar
    ttime = 20.000 fs
    pfactor = 2.000e+06 GPa*fs^2
    log_filename =  ./output/mol_equil10ps_300K.log
    traj_filename =  ./output/mol_equil10ps_300K.traj
    Elapsed time (equilib): 1515.3264671945944
    temperature =  300
    sigma = 1.0e+00 bar
    ttime = 20.000 fs
    pfactor = 2.000e+06 GPa*fs^2
    log_filename =  ./output/mol_product10ps_300K.log
    traj_filename =  ./output/mol_product10ps_300K.traj
    Elapsed time (product): 1547.582180276513
    Elapsed time (total)  : 3062.9086474711075



```python
t0 = perf_counter()

timestep = timestep_solid
temperature = target_temp
mask = None
print_interval = 10

prefix = 'solid'
prefix_equil = f'{prefix}_equil{int(timestep*steps_equilib/1000)}ps'
prefix_product = f'{prefix}_product{int(timestep*steps_equilib/1000)}ps'

# Equilibration
solid, dyn = run_md(solid, time_step=timestep, temperature=temperature, 
                    num_md_steps=steps_equilib, prefix=prefix_equil, 
                    log_interval=log_interval, traj_interval=traj_interval,
                    outputdir=outputdir, mask=mask)
t1 = perf_counter()
print(f'Elapsed time (equilib): {t1-t0}')

# Product run
solid, dyn = run_md(solid, time_step=timestep, temperature=temperature, 
                    num_md_steps=steps_product, prefix=prefix_product, 
                    log_interval=log_interval, traj_interval=traj_interval,
                    outputdir=outputdir, mask=mask)
t2 = perf_counter()
print(f'Elapsed time (product): {t2-t1}')
print(f'Elapsed time (total)  : {t2-t0}')
```

    temperature =  300
    sigma = 1.0e+00 bar
    ttime = 20.000 fs
    pfactor = 2.000e+06 GPa*fs^2
    log_filename =  ./output/solid_equil20ps_300K.log
    traj_filename =  ./output/solid_equil20ps_300K.traj
    Elapsed time (equilib): 1401.480807280168
    temperature =  300
    sigma = 1.0e+00 bar
    ttime = 20.000 fs
    pfactor = 2.000e+06 GPa*fs^2
    log_filename =  ./output/solid_product20ps_300K.log
    traj_filename =  ./output/solid_product20ps_300K.traj
    Elapsed time (product): 1573.7805970581248
    Elapsed time (total)  : 2975.261404338293
