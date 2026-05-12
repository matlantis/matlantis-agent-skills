# matlantis_features.ase_ext.optimize.FIRELBFGS#

_class _matlantis_features.ase_ext.optimize.FIRELBFGS(_atoms : Union[Atoms, Filter, NEB]_, _logfile : Optional[str] = '-'_, _trajectory : Optional[str] = None_, _master : Optional[bool] = None_, _append_trajectory : bool = False_, _force_consistent : Optional[bool] = None_, _switch : float = 0.05_, _switch_decrease_rate : float = 0.9_, _switch_patience_fire : int = 10_, _switch_patience_lbfgs : int = 10_, _maxstep_fire : float = 0.2_, _maxstep_lbfgs : float = 0.2_)[[source]](../_modules/matlantis_features/ase_ext/optimize/fire_lbfgs.html#FIRELBFGS)#
    

Bases: `Optimizer`

A optimization method that combine FIRE and LBFGS algorithm.

In this method, a threshold value is define. The FIRE (LBFGS) optimization is used when the maximum force is larger (smaller) than the switch hold. If the force and potential energy increase is observed after switching to LBFGS, the optimization methid is be switch back to FIRE and the threshold value will be reduced by certain ratio.

Methods

`__init__`(atoms[, logfile, trajectory, ...]) | Initialize an instance.  
---|---  
`attach`(function[, interval]) | Attach callback function.  
`call_observers`() |   
`close`() |   
`closelater`(fd) |   
`converged`(gradient) | Did the optimization converge?  
`dump`(data) |   
`get_number_of_steps`() |   
`initialize`() |   
`insert_observer`(function[, position, interval]) | Insert an observer.  
`irun`([fmax, steps]) | Run optimizer as generator.  
`load`() |   
`log`(gradient) | a dummy function as placeholder for a real logger, e.g.  
`openfile`(file, comm[, mode]) |   
`read`() |   
`run`([fmax, steps]) | Run optimizer.  
`step`() | Take a single step.  
`todict`() |   
  
__init__(_atoms : Union[Atoms, Filter, NEB]_, _logfile : Optional[str] = '-'_, _trajectory : Optional[str] = None_, _master : Optional[bool] = None_, _append_trajectory : bool = False_, _force_consistent : Optional[bool] = None_, _switch : float = 0.05_, _switch_decrease_rate : float = 0.9_, _switch_patience_fire : int = 10_, _switch_patience_lbfgs : int = 10_, _maxstep_fire : float = 0.2_, _maxstep_lbfgs : float = 0.2_) → None[[source]](../_modules/matlantis_features/ase_ext/optimize/fire_lbfgs.html#FIRELBFGS.__init__)#
    

Initialize an instance.

Parameters
    

  * **atoms** (_AtomsLike_) – The input structure.

  * **logfile** (_str_ _or_ _None_ _,__optional_) – Filename for log file. If ‘-’ is used, stdout is used. Defaults to ‘-‘.

  * **trajectory** (_str_ _or_ _None_ _,__optional_) – Filename for trajectory file. If None, no trajectory file will be saved. Defaults to None.

  * **master** (_bool_ _or_ _None_ _,__optional_) – If True, the current rank will save files. If None, only rank 0 will save files. Defaults to None.

  * **append_trajectory** (_bool_ _,__optional_) – If True and trajectory file exists, trajectory will be appended to the existing file. Defaults to False. Defaults to False.

  * **force_consistent** (_bool_ _or_ _None_ _,__optional_) – Use force-consistent energy calls. Defaults to False.

  * **switch** (_float_ _,__optional_) – The initial value of algorithm switch force threshold value. Defaults to 0.05.

  * **switch_decrease_rate** (_float_ _,__optional_) – The ratio of threshold value decreasing. Defaults to 0.9.

  * **switch_patience_fire** (_int_ _,__optional_) – Number of optimization steps without switching to LBFGS after maximum force lower than the threshold value. Defaults to 10.

  * **switch_patience_lbfgs** (_int_ _,__optional_) – Number of optimization steps without switching to FIRE after the energy increase is observed. Defaults to 10.

  * **maxstep_fire** (_float_ _,__optional_) – The maximum distance an atom can move per iteration when using FIRE. Defaults to 0.2.

  * **maxstep_lbfgs** (_float_ _,__optional_) – The maximum distance an atom can move per iteration when using LBFGS. Defaults to 0.2.




attach(_function_ , _interval =1_, _* args_, _** kwargs_)#
    

Attach callback function.

If _interval > 0_, at every _interval_ steps, call _function_ with arguments _args_ and keyword arguments _kwargs_.

If _interval <= 0_, after step _interval_ , call _function_ with arguments _args_ and keyword arguments _kwargs_. This is currently zero indexed.

call_observers()#
    

close()#
    

closelater(_fd_)#
    

converged(_gradient_)#
    

Did the optimization converge?

dump(_data_)#
    

get_number_of_steps()#
    

initialize()#
    

insert_observer(_function_ , _position =0_, _interval =1_, _* args_, _** kwargs_)#
    

Insert an observer.

This can be used for pre-processing before logging and dumping.

Examples
    
    
    >>> from ase.build import bulk
    >>> from ase.calculators.emt import EMT
    >>> from ase.optimize import BFGS
    ...
    ...
    >>> def update_info(atoms, opt):
    ...     atoms.info["nsteps"] = opt.nsteps
    ...
    ...
    >>> atoms = bulk("Cu", cubic=True) * 2
    >>> atoms.rattle()
    >>> atoms.calc = EMT()
    >>> with BFGS(atoms, logfile=None, trajectory="opt.traj") as opt:
    ...     opt.insert_observer(update_info, atoms=atoms, opt=opt)
    ...     opt.run(fmax=0.05, steps=10)
    True
    

irun(_fmax =0.05_, _steps =100000000_)#
    

Run optimizer as generator.

Parameters
    

  * **fmax** (_float_) – Convergence criterion of the forces on atoms.

  * **steps** (_int_ _,__default=DEFAULT_MAX_STEPS_) – Number of optimizer steps to be run.



Yields
    

**converged** (_bool_) – True if the forces on atoms are converged.

load()#
    

log(_gradient_)#
    

a dummy function as placeholder for a real logger, e.g. in Optimizer

openfile(_file_ , _comm_ , _mode ='w'_)#
    

read()#
    

run(_fmax =0.05_, _steps =100000000_)#
    

Run optimizer.

Parameters
    

  * **fmax** (_float_) – Convergence criterion of the forces on atoms.

  * **steps** (_int_ _,__default=DEFAULT_MAX_STEPS_) – Number of optimizer steps to be run.



Returns
    

**converged** – True if the forces on atoms are converged.

Return type
    

bool

step() → None[[source]](../_modules/matlantis_features/ase_ext/optimize/fire_lbfgs.html#FIRELBFGS.step)#
    

Take a single step.

todict()#
    

[ __ previous Structure Optimization API ](../matlantis_features.ase_ext.optimize.html "previous page") [ next matlantis_features.ase_ext.optimize __](matlantis_features.ase_ext.optimize.html "next page")
