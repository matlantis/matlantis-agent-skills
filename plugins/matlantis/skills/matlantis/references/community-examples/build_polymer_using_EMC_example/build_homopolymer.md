# Build an amorphous homopolymer geometry using EMC.
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

The polymer structure is defined by three SMILES strings, `smiles_center`, `smiles_left`, and `smiles_right`. The monomer unit is represented by `smiles_center` and both end of the chain are terminated by the moeclues represented by `smiles_left` and `smiles_right`. Asterisk `*` in the SMILES strings represents the connecting point of the monomer molecule.

In this example, the polymer chain consists of 10 momomers (8 central units terminated at both ends). The number of atoms in the system is more than 2,000 and the system density is set at 0.85 g/cm^3. An amorphous structure is created by Monte Carlo method using the PCFF force field.


```python
smiles_center = '*NCCCCCC(=O)*'
smiles_left = 'NCCCCCC(=O)*'
smiles_right = '*NCCCCCC(=O)O'

settings = dict(
            smiles_center = smiles_center,
            smiles_left = smiles_left,
            smiles_right = smiles_right,
            ntotal = 2000,   # Total number of atoms in the cell.
            density = 0.85,  # The mass density of the system [g/cm3]. Please use smaller value compared to the target density.
            field='pcff',    # The force field used in the Monte Carlo method.
            ring_depth = 9,  # Default to 'auto'. The max ring size in the molecules.
            build_dir = './build',
            lammps_prefix = 'Nylon6',   # Prefix for the input files of LAMMPS.
            project='homopolymer',      # Project name used in the input files of LAMMPS
            seed=12345,       # Random seed for the Monte Carlo simulation.
            repeat_center=8,  # 10 monomers with center + left + right 
)
```

### Comfirm the monomer structures defined by the SMILES strings.


```python
mols = [Chem.MolFromSmiles(smiles) for smiles in [smiles_center, smiles_left, smiles_right]]
view = Draw.MolsToGridImage(mols)
display(view)
```


    
![png](output_13_0.png)
    


### Confirm the EMC settings and the amorpohus structure with a small number of atoms.


```python
ntotal, settings['ntotal'] = settings['ntotal'], int(0.1 * settings['ntotal'])

builder = EMCInterface()
builder.verbose  =True
builder.setup('homopolymer', **settings)
builder.build()

basename =  f'{builder.settings["project"]}'
atoms = read(f'{basename}.pdb')
print(f'Number of atoms: {len(atoms)} (&lt; {settings["ntotal"]})')
settings['ntotal'] = ntotal

view_ngl([atoms], ['ball+stick'], replace_structure=True)
```

    {'center': '*NCCCCCC(=O)*', 'left': 'NCCCCCC(=O)*', 'right': '*NCCCCCC(=O)O', 'field': 'pcff', 'ntotal': 200, 'density': 0.85, 'ring_depth': 9, 'build_dir': './build', 'lammps_prefix': 'Nylon6', 'project': 'homopolymer', 'seed': 12345, 'emc_execute': 'false', 'repeat_center': 8, 'repeat_left': 1, 'repeat_right': 1}
    EMC Setup v4.1.5 (March 21, 2023), (c) 2004-2023 Pieter J. in 't Veld
    
    Info: reading script from "./setup.esh"
    Info: phase1 = {poly}
    Info: project = homopolymer
    Info: ntotal = 200
    Info: direction = x
    Info: shape = 1
    Info: force field type = "cff"
    Info: force field name = "EMC_interface/EMC/v9.4.4/field/pcff/pcff"
    Info: force field location = "."
    Info: build for LAMMPS script in "./build"
    Info: creating LAMMPS run script "homopolymer.in"
    Info: adding pressure sampling
    Info: creating EMC build script "build.emc"
    Info: assuming mol fractions
    
    
    
    (* EMC: Enhanced Monte Carlo simulations *)
    
    version 9.4.4, build Aug  2 2023 14:03:34, date Tue Aug 15 12:38:12 2023
    
    valid until Jul 31, 2024
    
    Info: script v1.0 started at Tue Aug 15 12:38:12 2023
    Info: variables = {seed -&gt; -1, ntotal -&gt; 200, fshape -&gt; 1, output -&gt;
            "homopolymer", field -&gt; "EMC_interface/EMC/v9.4.4/field/pcff/pcff", 
            location1 -&gt; "./", nav -&gt; 0.6022141179, temperature -&gt; 300, radius -&gt;
            5, nrelax -&gt; 100, weight_nonbond -&gt; 0.0001, weight_bond -&gt; 0.0001, 
            weight_focus -&gt; 1, cutoff -&gt; 9.5, charge_cutoff -&gt; 9.5, kappa -&gt; 4, 
            density1 -&gt; 0.85, lprevious -&gt; 0, lphase -&gt; 0, f_poly -&gt; 1, 
            chem_center -&gt; "*NCCCCCC(=O)*", chem_left -&gt; "NCCCCCC(=O)*", 
            chem_right -&gt; "*NCCCCCC(=O)O"}
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
              {id -&gt; center, chemistry -&gt; "*NCCCCCC(=O)*", depth -&gt; 6, charges -&gt;
                forcefield, charge -&gt; 0, terminator -&gt; false, nconnects -&gt; 4, 
                connects -&gt; {
                  {source -&gt; $end1, destination -&gt; {index -&gt; center, site -&gt;
                      $end2}}, 
                  {source -&gt; $end1, destination -&gt; {index -&gt; left, site -&gt;
                      $end1}}, 
                  {source -&gt; $end2, destination -&gt; {index -&gt; center, site -&gt;
                      $end1}}, 
                  {source -&gt; $end2, destination -&gt; {index -&gt; right, site -&gt;
                      $end1}}}}, 
              {id -&gt; left, chemistry -&gt; "NCCCCCC(=O)*", depth -&gt; 6, charges -&gt;
                forcefield, charge -&gt; 0, terminator -&gt; false, nconnects -&gt; 1, 
                connects -&gt; 
                {source -&gt; $end1, destination -&gt; {index -&gt; center, site -&gt;
                    $end1}}}, 
              {id -&gt; right, chemistry -&gt; "*NCCCCCC(=O)O", depth -&gt; 6, charges -&gt;
                forcefield, charge -&gt; 0, terminator -&gt; false, nconnects -&gt; 1, 
                connects -&gt; 
                {source -&gt; $end1, destination -&gt; {index -&gt; center, site -&gt;
                    $end2}}}}, ndeletes -&gt; 0, npolymers -&gt; 0}
    Info: field = {mode -&gt; apply, style -&gt; none, error -&gt; false, debug -&gt; false, 
            check -&gt; {atomistic -&gt; true, charge -&gt; true}}
    Info: applying './EMC_interface/EMC/v9.4.4/field/pcff/pcff.frc'
    Warning: increment pair {na, c_1} not found.
    Info: variables = {lg_center -&gt; 19, lg_left -&gt; 20, lg_right -&gt; 21, norm_poly
            -&gt; 1, l_poly -&gt; 193, mg_center -&gt; 113.16067, mg_left -&gt; 114.16864, 
            mg_right -&gt; 130.16807, m_poly -&gt; 1149.62207, f_poly -&gt; 1, norm -&gt; 193,
            n_poly -&gt; 1, ntotal -&gt; 0, mtotal -&gt; 0}
    Info: simulation = {units -&gt; {permittivity -&gt; 1, seed -&gt; seed}, types -&gt; {
              coulomb -&gt; {pair -&gt; {active -&gt; true, cutoff -&gt; charge_cutoff}}}}
    Info: clusters = {progress -&gt; none, nbodys -&gt; 0, nclusters -&gt; 0, ngrafts -&gt; 0,
            npolymers -&gt; 1, 
            polymer -&gt; 
            {id -&gt; poly, system -&gt; main, type -&gt; alternate, niterations -&gt; 1000, n
              -&gt; 1, groups -&gt; {center, left, right}, weights -&gt; {1, 1, 1}, nrepeat
              -&gt; {8, 1, 1}}}
          delta = {systems -&gt; 1, bodies -&gt; 0, clusters -&gt; 1, sites -&gt; 193}
    Info: field = {mode -&gt; apply, style -&gt; none, error -&gt; false, debug -&gt; false, 
            check -&gt; {atomistic -&gt; true, charge -&gt; true}}
    Info: applying './EMC_interface/EMC/v9.4.4/field/pcff/pcff.frc'
    Warning: increment pair {na, c_1} not found.
    Info: types = {cff -&gt; {pair -&gt; {active -&gt; true, mode -&gt; repulsive}}}
    Info: variables = {nphase1 -&gt; 193, mphase1 -&gt; 1149.62207, vphase1 -&gt;
            2245.87320812, lbox -&gt; 13.0956907686, lphase1 -&gt; 13.0956907686, lxx -&gt;
            13.0956907686, lyy -&gt; 13.0956907686, lzz -&gt; 13.0956907686, lzy -&gt; 0, 
            lzx -&gt; 0, lyx -&gt; 0, lphase -&gt; 13.0956907686, ntotal -&gt; 193, mtotal -&gt;
            1149.62207, vtotal -&gt; 2245.87320812}
    Info: build = {
            system -&gt; {id -&gt; main, temperature -&gt; 300, density -&gt; 0.85, split -&gt;
              false, flag -&gt; {charge -&gt; true, map -&gt; true, pbc -&gt; true, geometry
                -&gt; true}, periodic -&gt; {x -&gt; true, y -&gt; true, z -&gt; true}, geometry
              -&gt; {xx -&gt; 13.0956907686, yy -&gt; 13.0956907686, zz -&gt; 13.0956907686, 
                zy -&gt; 0, zx -&gt; 0, yx -&gt; 0}, deform -&gt; {xx -&gt; 1, yy -&gt; 1, zz -&gt; 1, 
                zy -&gt; 0, zx -&gt; 0, yx -&gt; 0}}, 
            select -&gt; {progress -&gt; list, frequency -&gt; 1, message -&gt; nkt, center -&gt;
              false, origin -&gt; {x -&gt; 0, y -&gt; 0, z -&gt; 0}, order -&gt; random, check -&gt;
              true, nclusters -&gt; 1, cluster -&gt; poly, 
              relax -&gt; {ncycles -&gt; 100, radius -&gt; 5}, name -&gt; "error", 
              grow -&gt; {check -&gt; all, method -&gt; energetic, cutoff -&gt; 0, 
                weight -&gt; {nonbonded -&gt; 0.0001, bonded -&gt; 0.0001, focus -&gt; 1}, 
                theta -&gt; 0, dphi -&gt; 1, nbonded -&gt; 20, ntrials -&gt; 20, niterations
                -&gt; 1000}}}
    Info: building 193 sites.
    
       progress/%      E[0]/nkT
    ---------------------------
                1             0
                2     0.8159381
                3    -0.6204093
                4      -1.41409
                5     -1.169483
                6    -0.7676135
                7     -0.408325
                8    -0.0844594
                9    -0.1862469
               10     0.1226585
               11     0.7067756
               12     0.6886112
               13     0.8714265
               15     0.9424728
               16     0.8928695
               17     0.8106713
               18     0.8021028
               19     0.8404028
               20     0.8546211
               21     0.8474078
               22      1.054581
               23     0.9965999
               24     0.9537725
               25      1.068803
               26     0.9215336
               27     0.8622943
               28      1.152815
               29     0.7743096
               30     0.8832074
               31      1.071096
               32      1.491069
               33      1.412585
               34      1.311339
               35      1.414137
               36      1.318114
               37       1.24802
               38      1.336997
               39      1.213777
               40      1.162637
               41      1.780092
               42      2.358689
               43      2.045759
               44      1.835055
               45      1.918762
               46      1.871039
               47      2.007745
               48      1.884687
               49      1.737559
               50      1.802219
               51      1.715221
               52        1.8334
               53      2.027834
               54      1.814004
               55      2.019107
               56      1.973228
               58      1.945736
               59      1.676051
               60      1.956653
               61      2.039224
               62      1.926788
               63      2.061557
               64      2.201955
               65      2.152751
               66      2.163303
               67      2.118072
               68      2.055814
               69      2.186148
               70      2.278114
               72      2.371516
               73      2.418419
               74      2.292225
               75      2.156477
               76      2.180852
               77      2.125216
               78      2.124903
               79      2.292233
               80       2.31294
               81      2.200027
               82      2.306206
               83      2.283162
               84      2.203002
               85      2.197147
               86      2.276478
               87      2.235909
               88      2.259489
               89      2.251174
               90       2.39275
               91      2.362432
               92      2.380916
               93      2.590187
               94       2.98133
               95       2.97184
               96      3.258321
               97       3.63824
               98       3.71246
               99      4.295616
              100      4.098831
    ---------------------------
          average      1.636597
    
          niterations = 193, &lt;niterations&gt; = 1
    Info: force = {style -&gt; none, message -&gt; nkt}
    
    (* Energy *)
    
                             E[0]/nkT 
    ---------------------------------
    cff.angle                1.903536 
    cff.angle_angle        0.01364377 
    cff.angle_angle_tor    -0.1614076 
    cff.angle_torsion     -0.09981628 
    cff.bond                0.6718893 
    cff.bond_angle       -0.006005401 
    cff.bond_bond        -0.000960316 
    cff.bond_bond_13                0 
    cff.end_bond_torsio  -0.003813308 
    cff.improper           0.03190671 
    cff.middle_bond_tor   -0.02540703 
    cff.pair                 2.590033 
    cff.torsion            -0.8147676 
    ---------------------------------
    total                    4.098831 
    
    Info: force = {style -&gt; init, message -&gt; nkt}
    
    (* Energy *)
    
                             E[0]/nkT 
    ---------------------------------
    cff.angle                1.903536 
    cff.angle_angle        0.01364377 
    cff.angle_angle_tor    -0.1614076 
    cff.angle_torsion     -0.09981628 
    cff.bond                0.6718893 
    cff.bond_angle       -0.006005401 
    cff.bond_bond        -0.000960316 
    cff.bond_bond_13                0 
    cff.end_bond_torsio  -0.003813308 
    cff.improper           0.03190671 
    cff.middle_bond_tor   -0.02540703 
    cff.pair                 2.590033 
    cff.torsion            -0.8147676 
    ---------------------------------
    total                    4.098831 
    
    Info: variables = {nl_poly -&gt; 1}
    Info: put = {name -&gt; "homopolymer", compress -&gt; true, detail -&gt; 3}
    Info: pdb = {name -&gt; "homopolymer", compress -&gt; false, extend -&gt; false, mode
            -&gt; put, system -&gt; main, length -&gt; auto, forcefield -&gt; cff, atomistic
            -&gt; full, depth -&gt; auto, charges -&gt; false, detect -&gt; false, atom -&gt;
            index, residue -&gt; index, segment -&gt; index, rank -&gt; false, hexadecimal
            -&gt; false, cut -&gt; false, pbc -&gt; true, map -&gt; false, unwrap -&gt; clusters,
            rigid -&gt; true, fixed -&gt; true, vdw -&gt; false, connectivity -&gt; false, 
            crystal -&gt; false, element -&gt; auto, parameters -&gt; false, flag -&gt; {
              charge -&gt; true, map -&gt; true, pbc -&gt; true, geometry -&gt; true}}
    Info: lammps = {name -&gt; "homopolymer", compress -&gt; false, system -&gt; main, mode
            -&gt; put, length -&gt; auto, units -&gt; none, forcefield -&gt; cff, shake -&gt;
            auto, atomistic -&gt; none, detect -&gt; true, cutoff -&gt; false, charges -&gt;
            true, ewald -&gt; true, ellipsoid -&gt; false, sphere -&gt; false, bonds -&gt;
            false, types -&gt; false, data -&gt; true, parameters -&gt; true, cross -&gt;
            false, variables -&gt; true, coefficients -&gt; true, comment -&gt; true, map
            -&gt; true, unwrap -&gt; true, flag -&gt; {charge -&gt; true, map -&gt; true, pbc -&gt;
              true, geometry -&gt; true}}
    Info: tsc: bonded = {3.216, 1.164, 0.5445, 1.352, 0.01103, 0.008832} s
    Info: tsc: simulation = 9.8 s, threads = 5.704 s
    Info: tsc: thread = 5.308 s
    Info: script v1.0 finished at Tue Aug 15 12:38:23 2023
    
    Info: Thank you for using EMC v9.4.4
    Info: In any publication of scientific results based in part or
    Info: completely on the use of EMC, please include this reference:
    Info: P.J. in 't Veld and G.C. Rutledge, Macromolecules 2003, 36, 7358
    
    
    
    Number of atoms: 193 (&lt; 200)



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'N', 'O', 'C', 'H'), v…


## Build an amorphous polymer structure


```python
builder = EMCInterface()
builder.setup('homopolymer', **settings)
```


```python
builder.build()
```


```python
savedir = os.path.join(outputdir, f'build_{builder.settings["project"]}')
print(savedir)
builder.savefiles(savedir)
```

    ./output/build_homopolymer



```python
basename =  f'{savedir}/{builder.settings["project"]}'
atoms_lammpsdata = read(f'{basename}.data', format='lammps-data')
atoms_lammpsdata = set_elements_lammpsdata(atoms_lammpsdata)
atoms = read(f'{basename}.pdb')

view_ngl([atoms, atoms_lammpsdata], replace_structure=True)
```


    HBox(children=(NGLWidget(max_frame=1), VBox(children=(Dropdown(description='Show', options=('All', 'N', 'O', '…


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
    LBFGS:    0 12:40:56    -8497.392387*       7.0775
    LBFGS:    1 12:40:57    -8529.433837*       5.9746
    LBFGS:    2 12:40:57    -8546.926021*       2.1862
    LBFGS:    3 12:40:57    -8560.358824*       4.4393
    LBFGS:    4 12:40:58    -8570.197520*       4.6909
    LBFGS:    5 12:40:58    -8576.513506*       3.0999





    False




```python
component_indices = list(get_connected_components(atoms, return_set=True))
components = [opt.atoms[list(indices)] for indices in component_indices]
formula = set([atoms.get_chemical_formula() for atoms in components])
print(formula)
view_ngl([opt.atoms] + components, ['ball+stick'], replace_structure=True)
```

    {'C60H112N10O11'}



    HBox(children=(NGLWidget(max_frame=10), VBox(children=(Dropdown(description='Show', options=('All', 'N', 'O', …
