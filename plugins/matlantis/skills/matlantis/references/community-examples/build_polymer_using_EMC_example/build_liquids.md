# Build a liquid geometry using EMC.
Copyright ENEOS Corporation as contributors to Matlantis contrib project

Enhanced Monte Carlo\
https://montecarlo.sourceforge.net/emc/Welcome.html \
https://matsci.org/c/emc/50

P.J. in 't Veld and G.C. Rutledge, Macromolecules 2003, 36, 7358

In this example, an amorphous polymer geometry is built using an package Enhanced Monte Calro (EMC), which is distributed under GPL v3 Lincense. 

## Import python packages


```python
import os
import sys
import numpy as np

from ase import Atoms, units
from ase.io import read, write, Trajectory
from ase.data import atomic_masses
from ase.optimize import LBFGS

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
```


```python
from pfcc_extras.visualize.view import view_ngl
from pfcc_extras.structure.ase_rdkit_converter import smiles_to_atoms, atoms_to_smiles
from pfcc_extras.structure.connectivity import get_connected_components
from matlantis_features.utils.calculators import pfp_estimator_fn, get_calculator
```


    


## Import a local package for executing EMC
EMC is a huge package which can be used to setup various kind of molecular simulations. 
For simplicity, this example provides three templates of the input file of EMC: amorphous homopolymer, liquids, and a mixture of homopolymer and liquids. These templates are available through an wrapper package EMC_interface.


```python
sys.path.insert(0, 'EMC_interface/src/emc_interface')
from emc_interface import EMCInterface
```


```python
def set_elements_lammpsdata(atoms):
    mass_number = {int(np.round(m)):i for i, m in enumerate(atomic_masses)}
    masses = atoms.get_masses()
    numbers = [mass_number[int(np.round(m))] for m in masses] 
    atoms.numbers=numbers
    return atoms
```

## System settings


```python
outputdir = './output'
os.makedirs(outputdir, exist_ok=True)
```

The liquid structure is defined by SMILES strings with molecular fractions. `name_smiles_fractions` is a dictionary which stores the names of the components and the tuples of (SMILES, fraction) pairs as the keys and values, respectively.

In this example, the liquid structure is composed of H$_{2}$0, C$_{2}$H$_{5}$OH, and NaCl. The number of atoms in the system is more than 2,000 and the system density is set at 0.85 g/cm^3. The liquid structure is created by Monte Carlo method using the PCFF force field.


```python
name_smiles_fractions = {"water":("O", 85), 
"alcohol":("CCO", 10),
"salt":("[Na+].[Cl-]", 5)}

settings = dict(
            name_smiles_fractions = name_smiles_fractions,
            ntotal = 2000,   # Number of atoms in the cell.
            density = 0.85,  # The mass density of the system [g/cm3]. Please use smaller value compared to the target density.
            field='pcff',    # The force field used in the Monte Carlo method.
            ring_depth = 9,  # Default to 'auto'. The max ring size in the molecules.
            build_dir = './build',
            lammps_prefix = 'liquid', # Prefix for the input files of LAMMPS.
            project='liquid',         # Project name used in the input files of LAMMPS
            seed=12345,               # Random seed for the Monte Carlo simulation.
)
```

### Comfirm the structure of molecules defined by the SMILES strings.


```python
smileses = [s for n,(s,f) in name_smiles_fractions.items()]
mols = [Chem.MolFromSmiles(smiles) for smiles in smileses]
view = Draw.MolsToGridImage(mols)
display(view)
```


    
![png](output_13_0.png)
    


### Confirm the EMC settings and the liquid structure with a small number of atoms.


```python
ntotal, settings['ntotal'] = settings['ntotal'], int(0.05 * settings['ntotal'])

builder = EMCInterface()
builder.verbose  =True
builder.setup('liquid', **settings)
builder.build()

basename =  f'{builder.settings["project"]}'
atoms = read(f'{basename}.pdb')
print(f'Number of atoms: {len(atoms)} ({settings["ntotal"]})')
settings['ntotal'] = ntotal

view_ngl([atoms], ['ball+stick'], replace_structure=True)
```

    {'field': 'pcff', 'ntotal': 100, 'density': 0.85, 'ring_depth': 9, 'build_dir': './build', 'lammps_prefix': 'liquid', 'project': 'liquid', 'seed': 12345, 'emc_execute': 'false', 'groups': 'water           O\nalcohol         CCO\nsalt            [Na+].[Cl-]', 'clusters': 'salt            salt,5\nalcohol         alcohol,10\nwater           water,85'}
    EMC Setup v4.1.5 (March 21, 2023), (c) 2004-2023 Pieter J. in 't Veld
    
    Info: reading script from "./setup.esh"
    Info: phase1 = {salt, alcohol, water}
    Info: project = liquid
    Info: ntotal = 100
    Info: direction = x
    Info: shape = 1
    Info: force field type = "cff"
    Info: force field name = "EMC_interface/EMC/v9.4.4/field/pcff/pcff"
    Info: force field location = "."
    Info: build for LAMMPS script in "./build"
    Info: creating LAMMPS run script "liquid.in"
    Info: adding pressure sampling
    Info: creating EMC build script "build.emc"
    Info: assuming mol fractions
    
    
    
    (* EMC: Enhanced Monte Carlo simulations *)
    
    version 9.4.4, build Aug  2 2023 14:03:34, date Tue Aug 15 12:35:40 2023
    
    valid until Jul 31, 2024
    
    Info: script v1.0 started at Tue Aug 15 12:35:40 2023
    Info: variables = {seed -&gt; -1, ntotal -&gt; 100, fshape -&gt; 1, output -&gt; "liquid",
            field -&gt; "EMC_interface/EMC/v9.4.4/field/pcff/pcff", location1 -&gt;
            "./", nav -&gt; 0.6022141179, temperature -&gt; 300, radius -&gt; 5, nrelax -&gt;
            100, weight_nonbond -&gt; 0.0001, weight_bond -&gt; 0.0001, weight_focus -&gt;
            1, cutoff -&gt; 9.5, charge_cutoff -&gt; 9.5, kappa -&gt; 4, density1 -&gt; 0.85, 
            lprevious -&gt; 0, lphase -&gt; 0, f_salt -&gt; 5, f_alcohol -&gt; 10, f_water -&gt;
            85, chem_water -&gt; "O", chem_alcohol -&gt; "CCO", chem_salt -&gt;
            "[Na+].[Cl-]"}
    Info: output = 
          {detail -&gt; 3, wide -&gt; false, expand -&gt; false, math -&gt; true, reduced -&gt;
            false, info -&gt; true, strict -&gt; true, warning -&gt; true, message -&gt; true,
            debug -&gt; false, exit -&gt; true, ignore -&gt; false}
    Info: field = {id -&gt; pcff, mode -&gt; cff, style -&gt; none, name -&gt;
            {"./EMC_interface/EMC/v9.4.4/field/pcff/pcff.frc", 
              "./EMC_interface/EMC/v9.4.4/field/pcff/pcff_templates.dat"}, 
            compress -&gt; false, error -&gt; true, debug -&gt; false, check -&gt; {atomistic
              -&gt; true, charge -&gt; true}}
    Info: importing './EMC_interface/EMC/v9.4.4/field/pcff/pcff.frc'
    Warning: duplicate field.bond entry 40 [c:c] omitted.
    Warning: duplicate field.bond entry 51 [c:h] omitted.
    Warning: duplicate field.bond entry 62 [c:o_2] omitted.
    Warning: duplicate field.bond entry 176 [c=:c=1] omitted.
    Warning: duplicate field.bond entry 328 [c_1:o_1] omitted.
    Warning: duplicate field.bond entry 364 [cp:cp] omitted.
    Warning: duplicate field.bond entry 365 [cp:cp] omitted.
    Warning: duplicate field.bond entry 367 [cp:h] omitted.
    Warning: duplicate field.bond entry 689 [o:p] omitted.
    Info: groups = {ngroups -&gt; 3, group -&gt; {
              {id -&gt; water, chemistry -&gt; "O", depth -&gt; 6, charges -&gt; forcefield, 
                charge -&gt; 0, terminator -&gt; false}, 
              {id -&gt; alcohol, chemistry -&gt; "CCO", depth -&gt; 6, charges -&gt;
                forcefield, charge -&gt; 0, terminator -&gt; false}, 
              {id -&gt; salt, chemistry -&gt; "[Na+].[Cl-]", depth -&gt; 6, charges -&gt;
                forcefield, charge -&gt; 0, terminator -&gt; false}}, ndeletes -&gt; 0, 
            npolymers -&gt; 0}
    Info: field = {mode -&gt; apply, style -&gt; none, error -&gt; false, debug -&gt; false, 
            check -&gt; {atomistic -&gt; true, charge -&gt; true}}
    Info: applying './EMC_interface/EMC/v9.4.4/field/pcff/pcff.frc'
    Info: variables = {lg_water -&gt; 3, lg_alcohol -&gt; 9, lg_salt -&gt; 2, l_salt -&gt; 2, 
            norm_salt -&gt; 1, l_alcohol -&gt; 9, norm_alcohol -&gt; 1, l_water -&gt; 3, 
            norm_water -&gt; 1, mg_water -&gt; 18.01534, mg_alcohol -&gt; 46.06952, mg_salt
            -&gt; 58.443, m_salt -&gt; 58.443, m_alcohol -&gt; 46.06952, m_water -&gt;
            18.01534, f_salt -&gt; 0.0281690140845, f_alcohol -&gt; 0.253521126761, 
            f_water -&gt; 0.718309859155, norm -&gt; 355, n_salt -&gt; 1, n_alcohol -&gt; 3, 
            n_water -&gt; 24, ntotal -&gt; 0, mtotal -&gt; 0}
    Info: simulation = {units -&gt; {permittivity -&gt; 1, seed -&gt; seed}, types -&gt; {
              coulomb -&gt; {pair -&gt; {active -&gt; true, cutoff -&gt; charge_cutoff}}}}
    Info: clusters = {progress -&gt; none, nbodys -&gt; 0, nclusters -&gt; 3, 
            cluster -&gt; {{id -&gt; salt, system -&gt; main, n -&gt; 1, 
                group -&gt; salt}, {id -&gt; alcohol, system -&gt; main, n -&gt; 3, 
                group -&gt; alcohol}, {id -&gt; water, system -&gt; main, n -&gt; 24, group -&gt;
                water}}, ngrafts -&gt; 0, npolymers -&gt; 0}
          delta = {systems -&gt; 1, bodies -&gt; 0, clusters -&gt; 29, sites -&gt; 101}
    Info: field = {mode -&gt; apply, style -&gt; none, error -&gt; false, debug -&gt; false, 
            check -&gt; {atomistic -&gt; true, charge -&gt; true}}
    Info: applying './EMC_interface/EMC/v9.4.4/field/pcff/pcff.frc'
    Info: types = {cff -&gt; {pair -&gt; {active -&gt; true, mode -&gt; repulsive}}}
    Info: variables = {nphase1 -&gt; 101, mphase1 -&gt; 629.01972, vphase1 -&gt;
            1228.83734872, lbox -&gt; 10.7110357178, lphase1 -&gt; 10.7110357178, lxx -&gt;
            10.7110357178, lyy -&gt; 10.7110357178, lzz -&gt; 10.7110357178, lzy -&gt; 0, 
            lzx -&gt; 0, lyx -&gt; 0, lphase -&gt; 10.7110357178, ntotal -&gt; 101, mtotal -&gt;
            629.01972, vtotal -&gt; 1228.83734872}
    Info: build = {
            system -&gt; {id -&gt; main, temperature -&gt; 300, density -&gt; 0.85, split -&gt;
              false, flag -&gt; {charge -&gt; true, map -&gt; true, pbc -&gt; true, geometry
                -&gt; true}, periodic -&gt; {x -&gt; true, y -&gt; true, z -&gt; true}, geometry
              -&gt; {xx -&gt; 10.7110357178, yy -&gt; 10.7110357178, zz -&gt; 10.7110357178, 
                zy -&gt; 0, zx -&gt; 0, yx -&gt; 0}, deform -&gt; {xx -&gt; 1, yy -&gt; 1, zz -&gt; 1, 
                zy -&gt; 0, zx -&gt; 0, yx -&gt; 0}}, 
            select -&gt; {progress -&gt; list, frequency -&gt; 1, message -&gt; nkt, center -&gt;
              false, origin -&gt; {x -&gt; 0, y -&gt; 0, z -&gt; 0}, order -&gt; random, check -&gt;
              true, nclusters -&gt; 3, cluster -&gt; {salt, alcohol, water}, 
              relax -&gt; {ncycles -&gt; 100, radius -&gt; 5}, name -&gt; "error", 
              grow -&gt; {check -&gt; all, method -&gt; energetic, cutoff -&gt; 0, 
                weight -&gt; {nonbonded -&gt; 0.0001, bonded -&gt; 0.0001, focus -&gt; 1}, 
                theta -&gt; 0, dphi -&gt; 1, nbonded -&gt; 20, ntrials -&gt; 20, niterations
                -&gt; 1000}}}
    Info: building 101 sites.
    
       progress/%      E[0]/nkT
    ---------------------------
                1             0
                3             0
                5             0
                7     0.2799886
                9      0.474497
               11     0.3327791
               13     0.3228805
               15      0.414677
               17     0.1683151
               19     0.5923352
               21     0.4862683
               23     0.5423353
               25     0.7474108
               27     0.6993818
               29     0.7511926
               31      1.018246
               33     0.7573549
               35     0.6389433
               37     0.5542428
               39     0.5539834
               41     0.5619311
               43      0.652573
               45     0.4595564
               47     0.6105327
               49     0.5576211
               50     0.6880128
               51     0.8970456
               52      0.846521
               53     0.7044774
               54     0.5806945
               55      0.622713
               56     0.6248574
               57     0.5818638
               58      0.614266
               59     0.8029058
               60     0.8905107
               61     0.7100114
               62     0.7817141
               63     0.7522645
               64     0.9740345
               65       1.04965
               66     0.9526705
               67     0.8549616
               68      1.079418
               69      1.008996
               70      1.057455
               71     0.9823167
               72     0.9738933
               73     0.9267173
               74     0.9962696
               75     0.9141739
               76     0.9276867
               77     0.9738731
               78       1.11108
               79     0.8390269
               80     0.9575821
               81     0.9998902
               82      1.053398
               83       1.02874
               84        1.0269
               85     0.8892267
               86     0.9049844
               87     0.9346562
               88      1.025412
               89     0.9959041
               90      1.186315
               91       1.10687
               92      1.151216
               93     0.9482834
               94      1.048922
               95      1.286814
               96      1.038365
               97      1.030338
               98     0.9610115
               99     0.8396215
              100     0.9392166
    ---------------------------
          average     0.7796157
    
          niterations = 101, &lt;niterations&gt; = 1
    Info: force = {style -&gt; none, message -&gt; nkt}
    
    (* Energy *)
    
                             E[0]/nkT 
    ---------------------------------
    cff.angle               0.3524173 
    cff.angle_angle      0.0001411964 
    cff.angle_angle_tor  -0.001589655 
    cff.angle_torsion     -0.03160109 
    cff.bond                0.3152275 
    cff.bond_angle       -0.008013846 
    cff.bond_bond        0.0001403451 
    cff.bond_bond_13                0 
    cff.end_bond_torsio   0.002771526 
    cff.improper                    0 
    cff.middle_bond_tor   0.004431474 
    cff.pair                0.4993311 
    cff.torsion            -0.1940392 
    ---------------------------------
    total                   0.9392166 
    
    Info: force = {style -&gt; init, message -&gt; nkt}
    
    (* Energy *)
    
                             E[0]/nkT 
    ---------------------------------
    cff.angle               0.3524173 
    cff.angle_angle      0.0001411964 
    cff.angle_angle_tor  -0.001589655 
    cff.angle_torsion     -0.03160109 
    cff.bond                0.3152275 
    cff.bond_angle       -0.008013846 
    cff.bond_bond        0.0001403451 
    cff.bond_bond_13                0 
    cff.end_bond_torsio   0.002771526 
    cff.improper                    0 
    cff.middle_bond_tor   0.004431474 
    cff.pair                0.4993311 
    cff.torsion            -0.1940392 
    ---------------------------------
    total                   0.9392166 
    
    Info: put = {name -&gt; "liquid", compress -&gt; true, detail -&gt; 3}
    Info: pdb = {name -&gt; "liquid", compress -&gt; false, extend -&gt; false, mode -&gt;
            put, system -&gt; main, length -&gt; auto, forcefield -&gt; cff, atomistic -&gt;
            full, depth -&gt; auto, charges -&gt; false, detect -&gt; false, atom -&gt; index,
            residue -&gt; index, segment -&gt; index, rank -&gt; false, hexadecimal -&gt;
            false, cut -&gt; false, pbc -&gt; true, map -&gt; false, unwrap -&gt; clusters, 
            rigid -&gt; true, fixed -&gt; true, vdw -&gt; false, connectivity -&gt; false, 
            crystal -&gt; false, element -&gt; auto, parameters -&gt; false, flag -&gt; {
              charge -&gt; true, map -&gt; true, pbc -&gt; true, geometry -&gt; true}}
    Info: lammps = {name -&gt; "liquid", compress -&gt; false, system -&gt; main, mode -&gt;
            put, length -&gt; auto, units -&gt; none, forcefield -&gt; cff, shake -&gt; auto, 
            atomistic -&gt; none, detect -&gt; true, cutoff -&gt; false, charges -&gt; true, 
            ewald -&gt; true, ellipsoid -&gt; false, sphere -&gt; false, bonds -&gt; false, 
            types -&gt; false, data -&gt; true, parameters -&gt; true, cross -&gt; false, 
            variables -&gt; true, coefficients -&gt; true, comment -&gt; true, map -&gt; true,
            unwrap -&gt; true, flag -&gt; {charge -&gt; true, map -&gt; true, pbc -&gt; true, 
              geometry -&gt; true}}
    Info: tsc: bonded = {0.2974, 0.07709, 0.06684, 0.09799, 0.0004282, 0.0003834} s
    Info: tsc: simulation = 1.34 s, threads = 0.6436 s
    Info: tsc: thread = 0.5592 s
    Info: script v1.0 finished at Tue Aug 15 12:35:42 2023
    
    Info: Thank you for using EMC v9.4.4
    Info: In any publication of scientific results based in part or
    Info: completely on the use of EMC, please include this reference:
    Info: P.J. in 't Veld and G.C. Rutledge, Macromolecules 2003, 36, 7358
    
    
    
    Number of atoms: 101 (100)



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'C', 'O', 'Na', 'H', '…


## Build a liquid structure


```python
builder = EMCInterface()
builder.setup('liquid', **settings)
```


```python
builder.build()
```


```python
savedir = os.path.join(outputdir, f'build_{builder.settings["project"]}')
builder.savefiles(savedir)
```


```python
basename =  f'{savedir}/{builder.settings["project"]}'
atoms_lammpsdata = read(f'{basename}.data', format='lammps-data')
atoms_lammpsdata = set_elements_lammpsdata(atoms_lammpsdata)
atoms = read(f'{basename}.pdb')

view_ngl([atoms, atoms_lammpsdata], replace_structure=True)
```


    HBox(children=(NGLWidget(max_frame=1), VBox(children=(Dropdown(description='Show', options=('All', 'C', 'O', '…


## Optimize the structure using PFP.
Finally, stable structure is obtaind by optimizing the structure using PFP.


```python
MODEL_VERSION='latest'
CALC_MODE='PBE_PLUS_D3'
estimator_fn = pfp_estimator_fn(model_version=MODEL_VERSION, calc_mode=CALC_MODE)
calculator = get_calculator(estimator_fn)
```


```python
atoms.calc = calculator
opt = LBFGS(atoms)
opt.run(steps=5)
```

           Step     Time          Energy         fmax
    *Force-consistent energies used in optimization.
    LBFGS:    0 12:37:20    -6997.821285*       1.7474
    LBFGS:    1 12:37:21    -7004.949977*       1.3038
    LBFGS:    2 12:37:21    -7019.550723*       2.3027
    LBFGS:    3 12:37:21    -7024.894041*       1.2307
    LBFGS:    4 12:37:21    -7031.139852*       1.7265
    LBFGS:    5 12:37:22    -7040.024333*       3.6010





    False




```python
component_indices = list(get_connected_components(atoms, return_set=True))
components = [opt.atoms[list(indices)] for indices in component_indices]
formula = set([atoms.get_chemical_formula() for atoms in components])
print(formula)
view_ngl([opt.atoms] + components, ['ball+stick'], replace_structure=True)
```

    {'ClNa', 'C2H6NaO', 'C2H10NaO3', 'C2H8NaO2', 'H2NaO', 'H4Na2O2', 'H8NaO4', 'H2Cl2Na2O', 'H4NaO2', 'C6H22Na4O5', 'Cl2Na', 'C2H6O', 'Cl', 'H6ClNaO3', 'H2O', 'C2H8ClNaO2'}



    HBox(children=(NGLWidget(max_frame=540), VBox(children=(Dropdown(description='Show', options=('All', 'C', 'O',…
