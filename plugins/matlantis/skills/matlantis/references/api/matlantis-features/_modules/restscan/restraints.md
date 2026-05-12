# Source code for restscan.restraints
    
    
    import math
    from typing import Any, Dict, Mapping, Optional, Tuple
    
    import numpy as np
    from ase import Atoms
    from ase.constraints import slice2enlist
    from ase.geometry import (
        get_angles,
        get_angles_derivatives,
        get_dihedrals,
        get_dihedrals_derivatives,
        get_distances,
        get_distances_derivatives,
    )
    
    from restscan.utils import NDArrayFloat, max_norm
    
    
    def sin_sigmoid(x: float) -> float:
        if x <= -0.5:
            return 0.0
        elif 0.5 <= x:
            return 1.0
        else:
            return 0.5 * (math.sin((x * math.pi)) + 1)
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint)class ScanRestraint:
        """Super class for scan restraints.
    
        Parameters:
    
            direction: 1 for scan towords large value. -1 for scan for towards small value.
    
            destination: The scan destination.
    
            exceed: In this class, scan force is turned off when it reaches the destination.
                    This is a value indicating how far the scan can exceed the destination.
                    This value is defined to smooth out the force field that Optimizer feels.
    
            fmax: Rsstraint scan force.
    
        """
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint.__init__)    def __init__(
            self, direction: int, destination: float, exceed: float, fmax: float
        ) -> None:
            super().__init__()
            if direction not in (-1, 1):
                raise ValueError(f"direction should be -1 or 1 ({direction}).")
            if exceed <= 0:
                raise ValueError(f"exceed should be larger than 0 ({exceed}).")
            if fmax <= 0:
                raise ValueError(f"fmax should be larger than 0 ({fmax}).")
            self.direction = direction
            self.destination = destination
            self.exceed = exceed
            self.fmax = fmax
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint.get_removed_dof)    def get_removed_dof(self, atoms: Atoms) -> int:
            return 0
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint.get_coefficient)    def get_coefficient(self, atoms: Atoms) -> float:
            colvar = self.get_colvar(atoms)
            return sin_sigmoid(
                -1
                * self.direction
                * (colvar - (self.destination + self.direction * self.exceed * 0.5))
                / self.exceed
            )
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint.adjust_positions)    def adjust_positions(self, atoms: Atoms, positions: NDArrayFloat) -> None:
            pass
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint.adjust_momenta)    def adjust_momenta(self, atoms: Atoms, momenta: NDArrayFloat) -> None:
            pass
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint.adjust_forces)    def adjust_forces(self, atoms: Atoms, forces: NDArrayFloat) -> None:
            original_forces: NDArrayFloat = atoms.get_forces(apply_constraint=False)
            colvar_forces: NDArrayFloat = (
                self.get_colvar_forces(atoms) * self.direction * -1
            )
            normalized_colvar_forces = colvar_forces / float(np.linalg.norm(colvar_forces))
            c_cancel = -1 * np.vdot(normalized_colvar_forces, original_forces)
            c_forces = self.fmax / max_norm(normalized_colvar_forces)
    
            external_forces = normalized_colvar_forces * max(c_forces + c_cancel, 0)
            alpha = self.get_coefficient(atoms)
            forces[:, :] += external_forces * alpha
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint.get_colvar)    def get_colvar(self, atoms: Atoms) -> float:
            raise NotImplementedError()
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint.get_colvar_forces)    def get_colvar_forces(self, atoms: Atoms) -> NDArrayFloat:
            raise NotImplementedError()
    
    
    
    
    def calculate_distance(atoms: Atoms, indices: Tuple[int, int]) -> float:
        positions = atoms.get_positions()
        i1, i2 = indices
        p1 = positions[i1]
        p2 = positions[i2]
        _, rij = get_distances(p1=p1, p2=p2, cell=atoms.cell, pbc=atoms.pbc)
        return float(rij)
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDistance.html#matlantis_features.features.reaction.rest_scan.ScanDistance)class ScanDistance(ScanRestraint):
        """Scan for distance.
    
        Parameters:
    
            direction: 1 for scan towords large distance. -1 for scan for towards small distance.
    
            destination: The scan finish when the distance reach this value.
    
            exceed: In this class, scan force is turned off when it reaches the destination.
                    This is a value indicating how far the scan can exceed the destination.
                    This value is defined to smooth out the force field that Optimizer feels.
    
            fmax: Rsstraint scan force.
        """
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDistance.html#matlantis_features.features.reaction.rest_scan.ScanDistance.__init__)    def __init__(
            self,
            indices: Tuple[int, int],
            direction: int,
            destination: float,
            exceed: float = 0.1,
            fmax: float = 0.05,
        ) -> None:
            super().__init__(
                direction=direction, destination=destination, exceed=exceed, fmax=fmax
            )
            self.indices = indices
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDistance.html#matlantis_features.features.reaction.rest_scan.ScanDistance.get_colvar)    def get_colvar(self, atoms: Atoms) -> float:
            return calculate_distance(atoms, self.indices)
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDistance.html#matlantis_features.features.reaction.rest_scan.ScanDistance.get_colvar_forces)    def get_colvar_forces(self, atoms: Atoms) -> NDArrayFloat:
            positions = atoms.get_positions()
            i1, i2 = self.indices
            p1 = positions[i1]
            p2 = positions[i2]
            vec_, _ = get_distances(p1=p1, p2=p2, cell=atoms.cell, pbc=atoms.pbc)
            vec = vec_[0]
            der = get_distances_derivatives(v0=vec, cell=atoms.cell, pbc=atoms.pbc)[0]
            forces: NDArrayFloat = np.zeros_like(positions)
            forces[i1] = -der[0]
            forces[i2] = -der[1]
            return forces
    
    
    
        def __repr__(self) -> str:
            return "ScanDistance(%d, %d)" % self.indices
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDistance.html#matlantis_features.features.reaction.rest_scan.ScanDistance.todict)    def todict(self) -> Dict[str, Any]:
            dct = {
                "name": "ScanDistance",
                "kwargs": {
                    "indices": self.indices,
                    "direction": self.direction,
                    "destination": self.destination,
                    "exceed": self.exceed,
                    "fmax": self.fmax,
                },
            }
            return dct
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDistance.html#matlantis_features.features.reaction.rest_scan.ScanDistance.index_shuffle)    def index_shuffle(self, atoms: Atoms, ind: NDArrayFloat) -> None:
            newa = [-1, -1]  # Signal error
            for new, old in slice2enlist(ind, len(atoms)):
                for i, a in enumerate(self.indices):
                    if old == a:
                        newa[i] = new
            if any(v == -1 for v in newa):
                raise IndexError("Constraint not part of slice")
            self.indices = (newa[0], newa[1])
    
    
    
    
    def calculate_angle(atoms: Atoms, indices: Tuple[int, int, int]) -> float:
        positions = atoms.get_positions()
        i1, i2, i3 = indices
        p1 = positions[i1]
        p2 = positions[i2]
        p3 = positions[i3]
        v0 = [p1 - p2]
        v1 = [p3 - p2]
        a = get_angles(v0=v0, v1=v1, cell=atoms.cell, pbc=atoms.pbc)
        return float(a)
    
    
    class ScanAngle(ScanRestraint):
        """Scan for distance.
    
        Parameters:
    
            direction: 1 for scan towords large distance. -1 for scan for towards small distance.
    
            destination: The scan finish when the distance reach this value.
    
            exceed: In this class, scan force is turned off when it reaches the destination.
                    This is a value indicating how far the scan can exceed the destination.
                    This value is defined to smooth out the force field that Optimizer feels.
    
            fmax: Rsstraint scan force.
        """
    
        def __init__(
            self,
            indices: Tuple[int, int, int],
            direction: int,
            destination: float,
            exceed: float = 0.1,
            fmax: float = 0.05,
        ) -> None:
            super().__init__(
                direction=direction, destination=destination, exceed=exceed, fmax=fmax
            )
            self.indices = indices
    
        def get_colvar(self, atoms: Atoms) -> float:
            return calculate_angle(atoms, self.indices)
    
        def get_colvar_forces(self, atoms: Atoms) -> NDArrayFloat:
            positions = atoms.get_positions()
            i1, i2, i3 = self.indices
            p1 = positions[i1]
            p2 = positions[i2]
            p3 = positions[i3]
            v0 = [p1 - p2]
            v1 = [p3 - p2]
            der = get_angles_derivatives(v0=v0, v1=v1, cell=atoms.cell, pbc=atoms.pbc)[0]
            forces: NDArrayFloat = np.zeros_like(positions)
            forces[i1] = -der[0]
            forces[i2] = -der[1]
            forces[i3] = -der[2]
            return forces
    
        def __repr__(self) -> str:
            return "ScanAngle(%d, %d, %d)" % self.indices
    
        def todict(self) -> Dict[str, Any]:
            dct = {
                "name": "ScanAngle",
                "kwargs": {
                    "indices": self.indices,
                    "direction": self.direction,
                    "destination": self.destination,
                    "exceed": self.exceed,
                    "fmax": self.fmax,
                },
            }
            return dct
    
        def index_shuffle(self, atoms: Atoms, ind: NDArrayFloat) -> None:
            newa = [-1, -1, -1]  # Signal error
            for new, old in slice2enlist(ind, len(atoms)):
                for i, a in enumerate(self.indices):
                    if old == a:
                        newa[i] = new
            if any(v == -1 for v in newa):
                raise IndexError("Constraint not part of slice")
            self.indices = (newa[0], newa[1], newa[2])
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.calculate_dihedral.html#matlantis_features.features.reaction.rest_scan.calculate_dihedral)def calculate_dihedral(atoms: Atoms, indices: Tuple[int, int, int, int]) -> float:
        i0, i1, i2, i3 = indices
        positions = atoms.get_positions()
        p0 = positions[i0]
        p1 = positions[i1]
        p2 = positions[i2]
        p3 = positions[i3]
        v0 = [p1 - p0]
        v1 = [p2 - p1]
        v2 = [p3 - p2]
        dijkl = get_dihedrals(v0=v0, v1=v1, v2=v2, cell=atoms.cell, pbc=atoms.pbc)
        return float(dijkl)
    
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDihedral.html#matlantis_features.features.reaction.rest_scan.ScanDihedral)class ScanDihedral(ScanRestraint):
        """Scan for dihedral.
    
        Parameters:
    
            direction: 1 for scan towords large angle. -1 for scan for towards small angle.
    
            destination: The scan will stop when the angle is rotated by this value.
    
            exceed: In this class, scan force is turned off when it reaches the destination.
                    This is a value indicating how far the scan can exceed the destination.
                    This value is defined to smooth out the force field that Optimizer feels.
    
            fmax: Rsstraint scan force.
    
            current: The value required to restart scan. Normally, you can ignore this value.
        """
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDihedral.html#matlantis_features.features.reaction.rest_scan.ScanDihedral.__init__)    def __init__(
            self,
            indices: Tuple[int, int, int, int],
            direction: int,
            destination: float,
            exceed: float = 1.0,
            fmax: float = 0.05,
            current: Optional[float] = None,
        ) -> None:
            self.current: Optional[float] = current
            super().__init__(
                direction=direction, destination=destination, exceed=exceed, fmax=fmax
            )
            self.indices = indices
    
    
    
        def _min_change_angle(self, x: float) -> float:
            if self.current is None:
                return x
            else:
                dx = x - self.current
                while dx < -180.0:
                    dx = dx + 360.0
                while dx > 180.0:
                    dx = dx - 360.0
                return self.current + dx
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDihedral.html#matlantis_features.features.reaction.rest_scan.ScanDihedral.adjust_forces)    def adjust_forces(self, atoms: Atoms, forces: NDArrayFloat) -> None:
            self.observe(atoms)
            super().adjust_forces(atoms, forces)
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDihedral.html#matlantis_features.features.reaction.rest_scan.ScanDihedral.get_colvar)    def get_colvar(self, atoms: Atoms) -> float:
            return self._min_change_angle(calculate_dihedral(atoms, self.indices))
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDihedral.html#matlantis_features.features.reaction.rest_scan.ScanDihedral.get_colvar_forces)    def get_colvar_forces(self, atoms: Atoms) -> NDArrayFloat:
            i0, i1, i2, i3 = self.indices
            positions = atoms.get_positions()
            p0 = positions[i0]
            p1 = positions[i1]
            p2 = positions[i2]
            p3 = positions[i3]
            v0 = [p1 - p0]
            v1 = [p2 - p1]
            v2 = [p3 - p2]
            der = get_dihedrals_derivatives(
                v0=v0, v1=v1, v2=v2, cell=atoms.cell, pbc=atoms.pbc
            )[0]
            forces: NDArrayFloat = np.zeros_like(positions)
            forces[i0] = -der[0]
            forces[i1] = -der[1]
            forces[i2] = -der[2]
            forces[i3] = -der[3]
            return forces
    
    
    
        def __repr__(self) -> str:
            return "ScanDihedral(%d, %d, %d, %d)" % self.indices
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDihedral.html#matlantis_features.features.reaction.rest_scan.ScanDihedral.todict)    def todict(self) -> Dict[str, Any]:
            dct = {
                "name": "ScanDihedral",
                "kwargs": {
                    "indices": self.indices,
                    "direction": self.direction,
                    "destination": self.destination,
                    "exceed": self.exceed,
                    "fmax": self.fmax,
                    "current": self.current,
                },
            }
            return dct
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDihedral.html#matlantis_features.features.reaction.rest_scan.ScanDihedral.index_shuffle)    def index_shuffle(self, atoms: Atoms, ind: NDArrayFloat) -> None:
            newa = [-1, -1, -1, -1]  # Signal error
            for new, old in slice2enlist(ind, len(atoms)):
                for i, a in enumerate(self.indices):
                    if old == a:
                        newa[i] = new
            if any(v == -1 for v in newa):
                raise IndexError("Constraint not part of slice")
            self.indices = (newa[0], newa[1], newa[2], newa[3])
    
    
    
    
    
    [[docs]](../../generated/matlantis_features.features.reaction.rest_scan.ScanDihedral.html#matlantis_features.features.reaction.rest_scan.ScanDihedral.observe)    def observe(self, atoms: Atoms) -> None:
            if self.current is None:
                self.current = self.get_colvar(atoms)
                # self.destination = self.current + self.direction * self.destination
            else:
                self.current = self.get_colvar(atoms)
    
    
    
    
    scan_classes = {"ScanDistance": ScanDistance, "ScanDihedral": ScanDihedral}
    
    
    def get_scan(atoms: Atoms, dct: Mapping[str, Any]) -> ScanRestraint:
        """Get and initialize your scan class."""
        name = dct["name"]
        kwargs = dct["kwargs"]
        return scan_classes[name](atoms=atoms, **kwargs)  # type: ignore
    
