from __future__ import annotations
from ase import io, Atoms
import numpy as np
import os
from pathlib import Path
from typing import Union

from mydftlib.config.config import Config as DefaultPath

from ase.build import add_vacuum

def gamma_struct(
    structure: Union[str, Atoms],
    struct_format: str = 'vasp',
    interlayer_shift: float = 0.0,
    increments: int = 10,
    duplicate_z: int = 10,
    shift_plane: float = 0.5,
    write_file: bool = True,
    dft_workflow_root: Path = DefaultPath.DFT_WORKFLOW_ROOT,
    write_dir: Union[str, Path] = "GAMMA_SURF",
) -> list[tuple[float, float, float]]:
    """
    Generalized Gamma surface POSCAR generator

    Parameters
    ----------
    structure : str or ASE Atoms
        CIF file path or ASE Atoms object.
    interlayer_shift : float
        Additional z-coordinate shift for the sliding layer in Angstrom
    increments : float
        Number of increments along x and y displacement (default: 10).
    duplicate_z : int
        Number of repetitions along z.
    shift_plane : float
        Fractional height (0–1) of the plane to shift atoms above.
    write_file : bool
        If True, writes displaced structures to file.
    write_dir : str or Path
        Directory to save files if write_file is True.

    Returns
    -------
    list of (x, y, energy)
        Interpolated GSFE surface data points.
    """

    # Load structure
    atoms = io.read(structure,format = struct_format) if isinstance(structure, (str,Path)) else structure.copy()
    atoms = atoms.repeat((1, 1, duplicate_z))
    original_positions = atoms.get_positions()
    lattice_vectors = atoms.cell

    limits = atoms.cell.lengths()
    z_mid = limits[2] * shift_plane + interlayer_shift
    areaxy = np.linalg.norm(np.cross(lattice_vectors[0] , lattice_vectors[1]))    

    # Displacement ranges
    x_range = np.linspace(0,1,increments)
    y_range = np.linspace(0,1,increments)

    add_vacuum(atoms,vacuum=10)
    for y in y_range:
        for x in x_range:
            atoms.set_positions(original_positions)
            burger = lattice_vectors[0] * x + lattice_vectors[1] * y
            for atom in atoms:
                if atom.position[2] > z_mid:
                    atom.position += burger
            atoms.wrap()
            
            
            if write_file:
                out_file = dft_workflow_root / write_dir / f"poscar_x_{round(x, 2)}_y_{round(y, 2)}.vasp"
                out_file.parent.mkdir(parents=True, exist_ok=True)  # Make sure the directory exists
                io.write(out_file, atoms, format="cif")



def gsfe_struct(
    structure: Union[str, Atoms],
    struct_format: str = 'vasp',
    interlayer_shift: float = 0.0,
    increments: int = 10,
    xmax : float = 1.0,
    ymax : float = 0.0,
    duplicate_z: int = 10,
    shift_plane: float = 0.5,
    write_file: bool = True,
    dft_workflow_root: Path = DefaultPath.DFT_WORKFLOW_ROOT,
    write_dir: Union[str, Path] = 'GSFE',
) -> list[tuple[float, float, float]]:
    """
    Generalized GSFE POSCAR generator

    Parameters
    ----------
    structure : str or ASE Atoms
        CIF file path or ASE Atoms object.
    interlayer_shift : float
        Additional z-coordinate shift for the sliding layer.
    increments : float
        Number of increments along x and y displacement (default: 10).
    xmax : float
        The final displacement length along the x direction (default: 1.0 Angstrom)
    ymax : float
        The final displacement length along the y direction (default: 0.0 Angstrom)
    duplicate_z : int
        Number of repetitions along z.
    shift_plane : float
        Fractional height (0–1) of the plane to shift atoms above.
    write_file : bool
        If True, writes displaced structures to file.
    write_dir : str or Path
        Directory to save files if write_file is True.

    Returns
    -------
    list of (x, y, energy)
        Interpolated GSFE surface data points.
    """

    # Load structure
    atoms = io.read(structure,format = struct_format) if isinstance(structure, (str,Path)) else structure.copy()
    atoms = atoms.repeat((1, 1, duplicate_z))
    original_positions = atoms.get_positions()
    lattice_vectors = atoms.cell

    limits = atoms.cell.lengths()
    z_mid = limits[2] * shift_plane + interlayer_shift
    areaxy = np.linalg.norm(np.cross(lattice_vectors[0] , lattice_vectors[1]))    

    # Displacement ranges
    x_range = np.linspace(0,xmax/np.linalg.norm(lattice_vectors[0]),increments)
    y_range = np.linspace(0,ymax/np.linalg.norm(lattice_vectors[1]),increments)

    add_vacuum(atoms,vacuum=10)
    for x, y in zip(x_range,y_range):
        atoms.set_positions(original_positions)
        burger = lattice_vectors[0] * x + lattice_vectors[1] * y
        for atom in atoms:
            if atom.position[2] > z_mid:
                atom.position += burger
        atoms.wrap()

        if write_file:
            out_file = dft_workflow_root / write_dir / f"poscar_x_{round(x, 2)}_y_{round(y, 2)}.vasp"
            out_file.parent.mkdir(parents=True, exist_ok=True)  # Make sure the directory exists
            io.write(out_file, atoms,format="cif")

