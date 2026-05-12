# Build a polymer/liquid mixture geometry using EMC.
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

The liquid structure is defined by SMILES strings with molecular fractions. `name_smiles_fractions` is a dictionary which stores the names of the components and the tuples of (SMILES, fraction) pairs as the keys and values, respectively.

In this example, the polymer chain consists of 10 momomers (8 central units terminated at both ends), and the liquid structure is composed of H$_{2}$0 and C$_{2}$H$_{5}$OH. The number of atoms in the system is more than 2,000 and the system density is set at 0.85 g/cm^3. The structure is created by Monte Carlo method using the PCFF force field.


```python
name_smiles_fractions = {"water":("O", 1), 
"alcohol":("CCO", 1),
# "salt":("[Na+].[Cl-]", 5)
}

smiles_center = '*NCCCCCC(=O)*'
smiles_left = 'NCCCCCC(=O)*'
smiles_right = '*NCCCCCC(=O)O'

settings = dict(
            name_smiles_fractions = name_smiles_fractions,
            smiles_center = smiles_center,
            smiles_left = smiles_left,
            smiles_right = smiles_right,
            polymer_fraction=1,
            ntotal = 2000,   # Total number of atoms in the cell.
            density = 0.85,  # The mass density of the system [g/cm3]. Please use smaller value compared to the target density.
            field='pcff',    # The force field used in the Monte Carlo method.
            ring_depth = 9,  # Default to 'auto'. The max ring size in the molecules.
            build_dir = './build',
            lammps_prefix = 'mixture',  # Prefix for the input file of LAMMPS.
            project='mixture',          # Project name used in the input files of LAMMPS.
            seed=12345,       # Random seed for the Monte Carlo simulation.
            repeat_center=8,  # 10 monomers with center + left + right 
)
```

### Comfirm the molecular structures defined by the SMILES strings.


```python
smileses = [s for n,(s,f) in name_smiles_fractions.items()]
smileses.extend([smiles_center, smiles_left, smiles_right])

mols = [Chem.MolFromSmiles(smiles) for smiles in smileses]
view = Draw.MolsToGridImage(mols)

display(view)
```


    
![png](output_13_0.png)
    


### Confirm the EMC settings and the structure of the mixture with a small number of atoms.


```python
ntotal, settings['ntotal'] = settings['ntotal'], int(0.2 * settings['ntotal'])

builder = EMCInterface()
builder.verbose  =True
builder.setup('mixture', **settings)
builder.build()

basename =  f'{builder.settings["project"]}'
atoms = read(f'{basename}.pdb')
print(f'Number of atoms: {len(atoms)} (&lt; {settings["ntotal"]})')
settings['ntotal'] = ntotal

view_ngl([atoms], ['ball+stick'], replace_structure=True)
```

    {'center': '*NCCCCCC(=O)*', 'left': 'NCCCCCC(=O)*', 'right': '*NCCCCCC(=O)O', 'field': 'pcff', 'ntotal': 400, 'density': 0.85, 'ring_depth': 9, 'build_dir': './build', 'lammps_prefix': 'mixture', 'project': 'mixture', 'seed': 12345, 'emc_execute': 'false', 'polymer_fraction': 1, 'repeat_center': 8, 'repeat_left': 1, 'repeat_right': 1, 'groups': 'water           O\nalcohol         CCO', 'clusters': 'water           water,1\nalcohol         alcohol,1\npoly            alternate,1'}
    EMC Setup v4.1.5 (March 21, 2023), (c) 2004-2023 Pieter J. in 't Veld
    
    Info: reading script from "./setup.esh"
    Info: phase1 = {water, alcohol, poly}
    Info: project = mixture
    Info: ntotal = 400
    Info: direction = x
    Info: shape = 1
    Info: force field type = "cff"
    Info: force field name = "EMC_interface/EMC/v9.4.4/field/pcff/pcff"
    Info: force field location = "."
    Info: build for LAMMPS script in "./build"
    Info: creating LAMMPS run script "mixture.in"
    Info: adding pressure sampling
    Info: creating EMC build script "build.emc"
    Info: assuming mol fractions
    
    
    
    (* EMC: Enhanced Monte Carlo simulations *)
    
    version 9.4.4, build Aug  2 2023 14:03:34, date Tue Aug 15 12:58:03 2023
    
    valid until Jul 31, 2024
    
    Info: script v1.0 started at Tue Aug 15 12:58:03 2023
    Info: variables = {seed -&gt; -1, ntotal -&gt; 400, fshape -&gt; 1, output -&gt;
            "mixture", field -&gt; "EMC_interface/EMC/v9.4.4/field/pcff/pcff", 
            location1 -&gt; "./", nav -&gt; 0.6022141179, temperature -&gt; 300, radius -&gt;
            5, nrelax -&gt; 100, weight_nonbond -&gt; 0.0001, weight_bond -&gt; 0.0001, 
            weight_focus -&gt; 1, cutoff -&gt; 9.5, charge_cutoff -&gt; 9.5, kappa -&gt; 4, 
            density1 -&gt; 0.85, lprevious -&gt; 0, lphase -&gt; 0, f_water -&gt; 1, f_alcohol
            -&gt; 1, f_poly -&gt; 1, chem_center -&gt; "*NCCCCCC(=O)*", chem_left -&gt;
            "NCCCCCC(=O)*", chem_right -&gt; "*NCCCCCC(=O)O", chem_water -&gt; "O", 
            chem_alcohol -&gt; "CCO"}
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
    Info: groups = {ngroups -&gt; 5, group -&gt; {
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
                    $end2}}}, 
              {id -&gt; water, chemistry -&gt; "O", depth -&gt; 6, charges -&gt; forcefield, 
                charge -&gt; 0, terminator -&gt; false}, 
              {id -&gt; alcohol, chemistry -&gt; "CCO", depth -&gt; 6, charges -&gt;
                forcefield, charge -&gt; 0, terminator -&gt; false}}, ndeletes -&gt; 0, 
            npolymers -&gt; 0}
    Info: field = {mode -&gt; apply, style -&gt; none, error -&gt; false, debug -&gt; false, 
            check -&gt; {atomistic -&gt; true, charge -&gt; true}}
    Info: applying './EMC_interface/EMC/v9.4.4/field/pcff/pcff.frc'
    Warning: increment pair {na, c_1} not found.
    Info: variables = {lg_center -&gt; 19, lg_left -&gt; 20, lg_right -&gt; 21, lg_water -&gt;
            3, lg_alcohol -&gt; 9, l_water -&gt; 3, norm_water -&gt; 1, l_alcohol -&gt; 9, 
            norm_alcohol -&gt; 1, norm_poly -&gt; 1, l_poly -&gt; 193, mg_center -&gt;
            113.16067, mg_left -&gt; 114.16864, mg_right -&gt; 130.16807, mg_water -&gt;
            18.01534, mg_alcohol -&gt; 46.06952, m_water -&gt; 18.01534, m_alcohol -&gt;
            46.06952, m_poly -&gt; 1149.62207, f_water -&gt; 0.0146341463415, f_alcohol
            -&gt; 0.0439024390244, f_poly -&gt; 0.941463414634, norm -&gt; 205, n_water -&gt;
            2, n_alcohol -&gt; 2, n_poly -&gt; 2, ntotal -&gt; 0, mtotal -&gt; 0}
    Info: simulation = {units -&gt; {permittivity -&gt; 1, seed -&gt; seed}, types -&gt; {
              coulomb -&gt; {pair -&gt; {active -&gt; true, cutoff -&gt; charge_cutoff}}}}
    Info: clusters = {progress -&gt; none, nbodys -&gt; 0, nclusters -&gt; 2, 
            cluster -&gt; {{id -&gt; water, system -&gt; main, n -&gt; 2, 
                group -&gt; water}, {id -&gt; alcohol, system -&gt; main, n -&gt; 2, group -&gt;
                alcohol}}, ngrafts -&gt; 0, npolymers -&gt; 1, 
            polymer -&gt; 
            {id -&gt; poly, system -&gt; main, type -&gt; alternate, niterations -&gt; 1000, n
              -&gt; 2, groups -&gt; {center, left, right}, weights -&gt; {1, 1, 1}, nrepeat
              -&gt; {8, 1, 1}}}
          delta = {systems -&gt; 1, bodies -&gt; 0, clusters -&gt; 6, sites -&gt; 410}
    Info: field = {mode -&gt; apply, style -&gt; none, error -&gt; false, debug -&gt; false, 
            check -&gt; {atomistic -&gt; true, charge -&gt; true}}
    Info: applying './EMC_interface/EMC/v9.4.4/field/pcff/pcff.frc'
    Warning: increment pair {na, c_1} not found.
    Info: types = {cff -&gt; {pair -&gt; {active -&gt; true, mode -&gt; repulsive}}}
    Info: variables = {nphase1 -&gt; 410, mphase1 -&gt; 2427.41386, vphase1 -&gt;
            4742.13560739, lbox -&gt; 16.8005947534, lphase1 -&gt; 16.8005947534, lxx -&gt;
            16.8005947534, lyy -&gt; 16.8005947534, lzz -&gt; 16.8005947534, lzy -&gt; 0, 
            lzx -&gt; 0, lyx -&gt; 0, lphase -&gt; 16.8005947534, ntotal -&gt; 410, mtotal -&gt;
            2427.41386, vtotal -&gt; 4742.13560739}
    Info: build = {
            system -&gt; {id -&gt; main, temperature -&gt; 300, density -&gt; 0.85, split -&gt;
              false, flag -&gt; {charge -&gt; true, map -&gt; true, pbc -&gt; true, geometry
                -&gt; true}, periodic -&gt; {x -&gt; true, y -&gt; true, z -&gt; true}, geometry
              -&gt; {xx -&gt; 16.8005947534, yy -&gt; 16.8005947534, zz -&gt; 16.8005947534, 
                zy -&gt; 0, zx -&gt; 0, yx -&gt; 0}, deform -&gt; {xx -&gt; 1, yy -&gt; 1, zz -&gt; 1, 
                zy -&gt; 0, zx -&gt; 0, yx -&gt; 0}}, 
            select -&gt; {progress -&gt; list, frequency -&gt; 1, message -&gt; nkt, center -&gt;
              false, origin -&gt; {x -&gt; 0, y -&gt; 0, z -&gt; 0}, order -&gt; random, check -&gt;
              true, nclusters -&gt; 3, cluster -&gt; {water, alcohol, poly}, 
              relax -&gt; {ncycles -&gt; 100, radius -&gt; 5}, name -&gt; "error", 
              grow -&gt; {check -&gt; all, method -&gt; energetic, cutoff -&gt; 0, 
                weight -&gt; {nonbonded -&gt; 0.0001, bonded -&gt; 0.0001, focus -&gt; 1}, 
                theta -&gt; 0, dphi -&gt; 1, nbonded -&gt; 20, ntrials -&gt; 20, niterations
                -&gt; 1000}}}
    Info: building 410 sites.
    
       progress/%      E[0]/nkT
    ---------------------------
                0             0
                1     0.1359911
                2     0.3156861
                3     0.6649138
                4     0.2466837
                5     0.5200354
                6   -0.08975174
                7     0.1184721
                8     -0.255089
                9    -0.1870438
               10   -0.03276051
               11  -0.004297709
               12    -0.1805734
               13     0.1208399
               14    0.06030153
               15     0.4810551
               16      0.558814
               17      1.184378
               18      1.244127
               19     0.9412187
               20      1.140481
               21      1.015785
               22      1.161644
               23      1.153356
               24      1.313417
               25      1.877579
               26      1.757096
               27      1.849571
               28      1.773361
               29      1.653729
               30      1.808911
               31      1.825831
               32      1.864993
               33       1.78609
               34      1.829452
               35      1.878162
               36      1.929598
               37      1.938106
               38      1.857655
               39      1.909555
               40        1.8176
               41      1.953736
               42      1.973238
               43      1.908306
               44      1.886531
               45      2.031299
               46      1.871511
               47      1.890593
               48      1.858203
               49      1.819135
               50      1.841864
               51      1.726088
               52      1.694557
               53       1.77785
               54      1.789856
               55      1.802545
               56      1.922006
               57      1.843414
               58       1.87664
               59      1.911094
               60      1.861114
               61        1.9178
               62      1.886437
               63      1.905704
               64      1.979776
               65      2.004301
               66      1.955082
               67      1.880253
               68      1.907425
               69      1.997224
               70      2.012778
               71      2.094785
               72      2.155713
               73      2.177148
               74      2.155081
               75      2.181065
               76      2.143043
               77      2.226853
               78      2.157804
               79      2.054323
               80      2.130594
               81      2.339544
               82      2.521075
               83      2.493252
               84      2.707073
               85      2.718689
               86      2.557206
               87      2.547981
               88      2.584635
               89       2.73994
               90      2.659722
               91      3.087891
               92      3.136939
               93      3.349444
               94      3.571546
               95      3.471617
               96      3.680436
               97      3.779361
               98      3.757244
               99      3.800831
              100      3.619079
    ---------------------------
          average      1.784864
    
          niterations = 410, &lt;niterations&gt; = 1
    Info: force = {style -&gt; none, message -&gt; nkt}
    
    (* Energy *)
    
                             E[0]/nkT 
    ---------------------------------
    cff.angle                1.829341 
    cff.angle_angle        0.01121326 
    cff.angle_angle_tor     -0.179041 
    cff.angle_torsion      -0.1221533 
    cff.bond                 0.635214 
    cff.bond_angle        -0.02739145 
    cff.bond_bond        0.0009878045 
    cff.bond_bond_13                0 
    cff.end_bond_torsio    0.00104441 
    cff.improper           0.05152967 
    cff.middle_bond_tor  -0.009794768 
    cff.pair                 2.164808 
    cff.torsion            -0.7366786 
    ---------------------------------
    total                    3.619079 
    
    Info: force = {style -&gt; init, message -&gt; nkt}
    
    (* Energy *)
    
                             E[0]/nkT 
    ---------------------------------
    cff.angle                1.829341 
    cff.angle_angle        0.01121326 
    cff.angle_angle_tor     -0.179041 
    cff.angle_torsion      -0.1221533 
    cff.bond                 0.635214 
    cff.bond_angle        -0.02739145 
    cff.bond_bond        0.0009878045 
    cff.bond_bond_13                0 
    cff.end_bond_torsio    0.00104441 
    cff.improper           0.05152967 
    cff.middle_bond_tor  -0.009794768 
    cff.pair                 2.164808 
    cff.torsion            -0.7366786 
    ---------------------------------
    total                    3.619079 
    
    Info: variables = {nl_poly -&gt; 2}
    Info: put = {name -&gt; "mixture", compress -&gt; true, detail -&gt; 3}
    Info: pdb = {name -&gt; "mixture", compress -&gt; false, extend -&gt; false, mode -&gt;
            put, system -&gt; main, length -&gt; auto, forcefield -&gt; cff, atomistic -&gt;
            full, depth -&gt; auto, charges -&gt; false, detect -&gt; false, atom -&gt; index,
            residue -&gt; index, segment -&gt; index, rank -&gt; false, hexadecimal -&gt;
            false, cut -&gt; false, pbc -&gt; true, map -&gt; false, unwrap -&gt; clusters, 
            rigid -&gt; true, fixed -&gt; true, vdw -&gt; false, connectivity -&gt; false, 
            crystal -&gt; false, element -&gt; auto, parameters -&gt; false, flag -&gt; {
              charge -&gt; true, map -&gt; true, pbc -&gt; true, geometry -&gt; true}}
    Info: lammps = {name -&gt; "mixture", compress -&gt; false, system -&gt; main, mode -&gt;
            put, length -&gt; auto, units -&gt; none, forcefield -&gt; cff, shake -&gt; auto, 
            atomistic -&gt; none, detect -&gt; true, cutoff -&gt; false, charges -&gt; true, 
            ewald -&gt; true, ellipsoid -&gt; false, sphere -&gt; false, bonds -&gt; false, 
            types -&gt; false, data -&gt; true, parameters -&gt; true, cross -&gt; false, 
            variables -&gt; true, coefficients -&gt; true, comment -&gt; true, map -&gt; true,
            unwrap -&gt; true, flag -&gt; {charge -&gt; true, map -&gt; true, pbc -&gt; true, 
              geometry -&gt; true}}
    Info: tsc: bonded = {7.077, 2.521, 1.202, 2.993, 0.02548, 0.01749} s
    Info: tsc: simulation = 21.45 s, threads = 12.75 s
    Info: tsc: thread = 11.97 s
    Info: script v1.0 finished at Tue Aug 15 12:58:27 2023
    
    Info: Thank you for using EMC v9.4.4
    Info: In any publication of scientific results based in part or
    Info: completely on the use of EMC, please include this reference:
    Info: P.J. in 't Veld and G.C. Rutledge, Macromolecules 2003, 36, 7358
    
    
    
    Number of atoms: 410 (&lt; 400)



    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H', 'N', 'C'), v…


## Build the mixture of liquid molecules and amorphous polymers.


```python
builder = EMCInterface()
builder.setup('mixture', **settings)
```


```python
builder.build()
```


```python
savedir = os.path.join(outputdir, f'build_{builder.settings["project"]}')
print(savedir)
builder.savefiles(savedir)
```

    ./output/build_mixture



```python
basename =  f'{savedir}/{builder.settings["project"]}'
atoms_lammpsdata = read(f'{basename}.data', format='lammps-data')
atoms_lammpsdata = set_elements_lammpsdata(atoms_lammpsdata)
atoms = read(f'{basename}.pdb')

view_ngl([atoms, atoms_lammpsdata], replace_structure=True)
```


    HBox(children=(NGLWidget(max_frame=1), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H', '…


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
    LBFGS:    0 13:00:37    -8957.630457*       7.0484
    LBFGS:    1 13:00:37    -8991.848565*       6.2971
    LBFGS:    2 13:00:37    -9009.746541*       1.7679
    LBFGS:    3 13:00:38    -9023.342350*       4.4884
    LBFGS:    4 13:00:38    -9033.152521*       3.6167
    LBFGS:    5 13:00:39    -9038.953378*       2.9377





    False




```python
component_indices = list(get_connected_components(atoms, return_set=True))
components = [opt.atoms[list(indices)] for indices in component_indices]
formula = set([atoms.get_chemical_formula() for atoms in components])
print(formula)
view_ngl([opt.atoms] + components, ['ball+stick'], replace_structure=True)
```

    {'H2O', 'C2H6O', 'C60H112N10O11'}



    HBox(children=(NGLWidget(max_frame=30), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H', …
