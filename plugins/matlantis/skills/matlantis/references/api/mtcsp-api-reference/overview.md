  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * Overview



# Overview#

MTCSP is a Python package for performing Crystal Structure Prediction on Matlantis. This package consists of three main categories of APIs:

  1. High-level search APIs: wraps the overall CSP workflow implemented by the subsequent low-level APIs.

  2. Low-level search APIs: implements the core functionalities such as structure generation, sampling, and evaluation.

  3. Utility APIs: provides helper functions for post analysis, visualization, and data export.




The MTCSP’s modules are categorized with this three layers in the following diagram.


## High-level search APIs#

MTCSP currently provides the following two search types as high-level search APIs:

  * [mtcsp.convex_hull_search](mtcsp.convex_hull_search.html#api-convex-hull-search): corresponds to CSP in entire compositional space.

  * [mtcsp.derivative_structure](mtcsp.derivative_structure.html#api-derivative-structure): corresponds to CSP with a given prototype structure.




## Low-level search APIs#

MTCSP is modularly designed so that each functionality is implemented in a separate module:

  * A database during CSP search is managed by [`mtcsp.experiment.Experiment()`](mtcsp.experiment.html#mtcsp.experiment.Experiment "mtcsp.experiment.Experiment").

  * Each trial (structure sampling and its evaluation) is serialized as [`mtcsp.crystal_structure.CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure").

  * Configurations of crystal structures (such as lattice, positions, and atomic species) are separately stored against the database file and managed by [`mtcsp.atoms.AtomsStore`](mtcsp.atoms.html#mtcsp.atoms.AtomsStore "mtcsp.atoms.AtomsStore").




The following component diagram illustrates the relationship among these low-level search APIs.
    
    
            ---
    title: MTCSP Low-level Search APIs Overview
    ---
    flowchart TB
       experiment["mtcsp.experiment
       .Experiment
       Manages DB in CSP search"]
    
       frozen["mtcsp.experiment
       .FrozenExperiment
       Query DB after CSP search"]
    
       atoms["mtcsp.atoms
       .AtomsStore
       Store Atoms in JSON format"]
    
       sampler["CSPSampler
       Generate structures"]
    
       db[("Database
       *.journal or *.sqlite3")]
    
       fs[("atoms_store_dir
       Directory to store Atoms")]
    
       experiment --> atoms
       experiment --> sampler
       experiment --> db
       frozen --> atoms
       frozen --> db
       atoms --> fs
        

Besides, MTCSP is built on top of [Optuna](https://optuna.org/) to manage the CSP workflow. If you are familiar with Optuna, the following mapping table between Optuna and MTCSP concepts would be helpful to understand the MTCSP low-level search APIs.

Optuna Concept | MTCSP Concept  
---|---  
[`optuna.study.Study`](https://optuna.readthedocs.io/en/stable/reference/generated/optuna.study.Study.html#optuna.study.Study "\(in Optuna v4.7.0\)") w/ sampler | [`mtcsp.experiment.Experiment()`](mtcsp.experiment.html#mtcsp.experiment.Experiment "mtcsp.experiment.Experiment") `mtcsp.experiment.DerivativeExperiment()`  
[`optuna.study.Study`](https://optuna.readthedocs.io/en/stable/reference/generated/optuna.study.Study.html#optuna.study.Study "\(in Optuna v4.7.0\)") w/o sampler | [`mtcsp.experiment.FrozenExperiment()`](mtcsp.experiment.html#mtcsp.experiment.FrozenExperiment "mtcsp.experiment.FrozenExperiment")  
[`optuna.trial.Trial`](https://optuna.readthedocs.io/en/stable/reference/generated/optuna.trial.Trial.html#optuna.trial.Trial "\(in Optuna v4.7.0\)") | [`mtcsp.crystal_structure.CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")  
[`optuna.artifacts.FileSystemArtifactStore`](https://optuna.readthedocs.io/en/stable/reference/artifacts.html#optuna.artifacts.FileSystemArtifactStore "\(in Optuna v4.7.0\)") | [`mtcsp.atoms.FileSystemAtomsStore`](mtcsp.atoms.html#mtcsp.atoms.FileSystemAtomsStore "mtcsp.atoms.FileSystemAtomsStore")  
  
[ previous MTCSP Reference ](mtcsp.html "previous page") [ next mtcsp.analysis package ](mtcsp.analysis.html "next page")

On this page 

  * High-level search APIs
  * Low-level search APIs


