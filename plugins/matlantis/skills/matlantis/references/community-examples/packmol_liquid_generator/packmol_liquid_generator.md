### Fast Liquid Genaration by Packmol
[LiquidGenerator](https://github.com/matlantis-pfcc/matlantis-contrib/tree/main/matlantis_contrib_examples/liquid_generator) was provided to create the initial liquid structure.
PackmolLiquidGenerator uses Packmol as its engine, making it possible to generate liquid structures more than 10 times faster, even for large models.

### Install packmol
Clone the source code from git and build packmol.


```python
!git clone https://github.com/m3g/packmol.git
!cd packmol &amp;&amp; ./configure &amp;&amp; make
```

    fatal: destination path 'packmol' already exists and is not an empty directory.
    Setting compiler to /usr/bin/gfortran
     ------------------------------------------------------ 
     Compiling packmol with /usr/bin/gfortran 
     Flags: -O3 --fast-math -march=native -funroll-loops 
     ------------------------------------------------------ 
     ------------------------------------------------------ 
     Packmol succesfully built.
     ------------------------------------------------------ 


### `PackmolLiquidGenerator` class
Although packmol can create and run its own scripts and xyz files, the `PackmolLiquidGenerator` class can be used to create ASE Atoms from multiple SMILES.


```python
import os
import subprocess as sub
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, TypedDict, Tuple
import warnings

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.units import _Nav
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt


class Composition(TypedDict):
    name: str
    smiles: str
    number: int


class Cell(TypedDict):
    lx: float
    ly: float
    lz: float

@dataclass
class PackmolLiquidGenerator:
    outputname: str = "packmol.xyz"
    outputdir: Path = "output"
    composition: List[Composition] = field(default_factory=list)
    density: float = 0.5
    """ target density in g/cm^3 """
    pbc: Tuple[bool, bool, bool] = (True, True, True)
    packmol_bin: Path = Path("packmol/packmol")
    inputname: str = "packmol.inp"
    filetype: str = "xyz"
    """ same as filetype in packmol. Please see https://m3g.github.io/packmol/userguide.shtml """
    tolerance: float = 2.0
    """ same as tolerance in packmol. Please see https://m3g.github.io/packmol/userguide.shtml """
    margin: float = 1.0
    """ margin for box and atom """
    cell: Optional[Cell] = None
    verbose: bool = False

    def __post_init__(self) -&gt; None:
        if not os.path.isdir(self.outputdir):
            os.makedirs(self.outputdir)

        if not self.cell:
            total_mass = 0
            for props in self.composition:
                mol = Chem.MolFromSmiles(props["smiles"])
                mol = AllChem.AddHs(mol)
                total_mass += props["number"] * ExactMolWt(mol)
            total_mass /= _Nav
            l = round((total_mass / self.density) ** (1 / 3) * 10**8, 2)
            self.cell = {"lx": l, "ly": l, "lz": l}

    def _smiles_to_atoms(self, smiles: str) -&gt; Atoms:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        atoms = Atoms(
            positions=mol.GetConformer().GetPositions(),
            numbers=np.array([a.GetAtomicNum() for a in mol.GetAtoms()]),
        )
        return atoms

    def run(self) -&gt; Atoms:
        assert self.cell is not None
        assert self.packmol_bin.exists()

        with open(os.path.join(self.outputdir, self.inputname), "w") as f:
            print(
                f"tolerance {self.tolerance}",
                f"filetype {self.filetype}",
                f"output {os.path.join(self.outputdir, self.outputname)}",
                sep="\n",
                end="\n\n",
                file=f,
            )
            for props in self.composition:
                fname = os.path.join(self.outputdir, props["name"] + "." + self.filetype)
                atoms = self._smiles_to_atoms(props["smiles"])
                write(fname, atoms)
                print(
                    f"structure {fname}",
                    f"  number {props['number']}",
                    f"  inside box {self.margin} {self.margin} {self.margin} {self.cell['lx']-self.margin} {self.cell['ly']-self.margin} {self.cell['lz']-self.margin}",
                    "end structure",
                    sep="\n",
                    end="\n\n",
                    file=f,
                )

        cmd = [self.packmol_bin]
        stdout = None if self.verbose else sub.DEVNULL

        with open(os.path.join(self.outputdir, self.inputname), "rt") as f:
            sub.run(cmd, text=True, stdout=stdout, stdin=f)

        packed_atoms = read(os.path.join(self.outputdir, self.outputname))
        packed_atoms.set_cell([self.cell["lx"], self.cell["ly"], self.cell["lz"]])
        packed_atoms.set_pbc(self.pbc)

        write(os.path.join(self.outputdir, self.outputname), packed_atoms)

        return packed_atoms
```

### Example: Mixture of 1000 water moleculs and 1000 ethanol molecules


```python
liq_generator = PackmolLiquidGenerator(
    density=1.0,
    composition=[
        {"name": "water", "smiles": "O", "number": 1000},
        {"name": "ehtanol", "smiles": "CCO", "number": 1000},
    ],
    outputname="mixture.xyz",
    outputdir="output"
)
atoms = liq_generator.run()
```


```python
atoms
```




    Atoms(symbols='C2000H8000O2000', pbc=True, cell=[47.38, 47.38, 47.38])




```python
from ase.visualize import view
view(atoms, viewer="ngl")
```


    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'O', 'H', 'C'), value=…

