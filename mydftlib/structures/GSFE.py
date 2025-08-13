from __future__ import annotations
from ase import io, Atoms
import numpy as np
import os
from pathlib import Path
from typing import Union
from mydftlib.core.POSCAR import POSCAR
from ase.io import read,write
from mydftlib.config.config import Config as DefaultPath
import itertools
from ase.build import add_vacuum

import numpy as np
import itertools
import subprocess
import sys

from ase import io, visualize, Atoms
from ase.data import atomic_numbers, covalent_radii
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier
import atomman as am
import atomman.unitconvert as uc

def slicer(structure_file,tol = 0):
    # Load your atomic structure
    pipeline = import_file(structure_file)  # or .xyz, .vasp, etc.

    # Set up the Voronoi analysis modifier.
    voro = VoronoiAnalysisModifier(
        compute_indices = True,
        use_radii = True,
        edge_threshold = 0.1
    )
    pipeline.modifiers.append(voro)

    # Let OVITO compute the results.
    data = pipeline.compute()

    # Access computed Voronoi indices.
    # This is an (N) x (M) array, where M is the maximum face order.
    voro_indices = data.particles['Voronoi Index']
    voro_volume = data.particles['Atomic Volume']
    positions = np.array(data.particles.positions)

    df = pd.DataFrame({'position': positions[:,2],'volume': voro_volume})
    df_sorted = df.sort_values(by='position')

    temp_vol = 0
    prev_layer = "x"
    tol=0
    volume_at_each_layer = []
    volumes_of_layers = []
    for layer, vol in zip(df_sorted.position,df_sorted.volume):
        if layer == prev_layer:
            temp_vol += vol
        if layer != prev_layer:
            if prev_layer != "x":
                volume_at_each_layer.append([(layer - prev_layer)/2 + prev_layer, round(vol - temp_vol,tol)])
            volumes_of_layers.append([prev_layer,round(temp_vol,tol)])
            temp_vol = vol
            prev_layer = layer

    volume_at_each_layer = np.array(volume_at_each_layer)
    unique_layer = np.array([])
    unique_pos = np.array([])
    layer_diff = np.array([])
    
    for pos,item in zip(volume_at_each_layer[:,0],volume_at_each_layer[:,1]):
        if item in unique_layer:
            continue
        unique_layer = np.append(unique_layer, item)
        unique_pos = np.append(unique_pos, pos)
        if len(unique_layer) > 1:
            layer_diff = np.append(layer_diff, (unique_layer[-1]-unique_layer[-2]))
    return unique_pos


def normalize_miller_indices(indices):
    """Normalizes a set of Miller indices to ensure uniqueness."""
    indices = np.array(indices)
    norms = np.linalg.norm(indices, axis=1)
    normalized = indices / norms[:, None]

    return normalized

def unique_slip_directions(miller_indices,miller_type = "hexagonal"):
    """Eliminates duplicate slip directions based on normalized direction."""

    if miller_type == "hexagonal":
        normalized_indices = normalize_miller_indices(miller_indices)
        unique_dirs = []
        seen = list()

        for idx, norm_vec in zip(miller_indices, normalized_indices):
            key = tuple(np.round(norm_vec, decimals=6))
            key_elements = [round(x, 6) for x in key]
            key_together = key_elements + [-x for x in key_elements]
            all_combinations = list(itertools.permutations(key_together, 4))

            flattened_combinations = [tuple(x) for x in all_combinations]

            if key not in seen:
                seen= seen + flattened_combinations
                unique_dirs.append(idx)

        return unique_dirs

    normalized_indices = normalize_miller_indices(miller_indices)
    unique_dirs = []
    seen = list()

    for idx, norm_vec in zip(miller_indices, normalized_indices):
        key = tuple(np.round(norm_vec, decimals=6))
        key_elements = [round(x, 6) for x in key]
        key_together = key_elements + [-x for x in key_elements]
        all_combinations = list(itertools.permutations(key_together, 3))

        flattened_combinations = [tuple(x) for x in all_combinations]

        if key not in seen:
            seen= seen + flattened_combinations
            unique_dirs.append(idx)

    return unique_dirs

def miller_indices_creator(integ, miller_type="hexagonal"):
    if miller_type != "hexagonal":
        raise NotImplementedError("Only hexagonal system is implemented.")

    # Define the range of values for h, k, l
    vals = list(range(-integ, integ + 1))
    
    miller_indices = []

    for h, k, l in itertools.product(vals, repeat=3):
        i = -(h + k)
        # Ensure i is in range
        if i in vals:
            if [h, k, i, l] != [0, 0, 0, 0]:
                miller_indices.append([h, k, i, l])

    # Remove duplicates if any
    unique_indices = [list(x) for x in set(tuple(m) for m in miller_indices)]
    
    return sorted(unique_indices)

def orthonormal_basis_gen_hexagonal(hkil, a=1.0, c=1.0):
    """
    Generate orthonormal basis for a hexagonal (hkil) direction.
    """
    # Convert to 3D Cartesian vector
    normal = miller_bravais_to_cartesian(hkil, a, c)
    normal = normal / np.linalg.norm(normal)

    # Choose arbitrary vector not parallel to normal
    if np.allclose(normal, [1, 0, 0]):
        arbitrary = np.array([0, 1, 0])
    else:
        arbitrary = np.array([1, 0, 0])
    
    # First perpendicular vector
    orth1 = np.cross(normal, arbitrary)
    orth1 /= np.linalg.norm(orth1)

    # Second perpendicular vector
    orth2 = np.cross(normal, orth1)
    orth2 /= np.linalg.norm(orth2)
    
    vector_cartesian = np.array([orth1, orth2, normal])
    orth1_vec4 = am.tools.miller.vector3to4(orth1)
    orth2_vec4 = am.tools.miller.vector3to4(orth2)
    normal_vec4 = am.tools.miller.plane3to4(normal)
    
    miller_bravais = [orth1_vec4,orth2_vec4,normal_vec4]
    normalized_miller_bravais = None
    # normalized_orth1 = normalize_miller_bravais(orth1_vec4, tol=1e-5)
    # normalized_orth2 = normalize_miller_bravais(orth2_vec4, tol=1e-5)
    # normalized_normal = normalize_miller_bravais(normal_vec4, tol=1e-5)
    
    # normalized_miller_bravais = [normalized_orth1,normalized_orth2,normalized_normal]
    return normalized_miller_bravais, miller_bravais

import numpy as np

import numpy as np

def format_miller_bravais_vector(vec, tol=1e-5):
    """
    Formats a 4-index vector (hkil) into a string like [1 -2 1 0],
    ensuring the constraint h + k + i = 0.
    """

    # Try to find integer multiple to scale to whole numbers
    for s in range(1, 1000):
        scaled = vec * s
        if np.allclose(scaled, np.round(scaled), atol=tol):
            scaled = np.round(scaled).astype(int)
            break
    else:
        raise ValueError("Cannot find integer scaling")

    # Enforce i = - (h + k)
    h, k, _, l = scaled
    i = - (h + k)
    scaled[2] = i  # enforce

    # Reduce by GCD
    nonzero = scaled[scaled != 0]
    gcd = np.gcd.reduce(nonzero) if nonzero.size > 0 else 1
    scaled //= gcd

    # Format with minus sign
    return "[" + " ".join(str(x) for x in scaled) + "]"


def normalize_miller_bravais(vectors, tol=1e-5):
    """
    Normalize Miller-Bravais indices to smallest integer representation.

    Parameters:
    - vectors: (N, 4) array-like of float Miller-Bravais directions
    - tol: tolerance for rounding to integers

    Returns:
    - normalized_vectors: (N, 4) array of integer Miller-Bravais indices
    """
    vectors = np.array(vectors)

    # Flatten all values to find a common multiple
    flat = vectors.flatten()
    nonzero = flat[np.abs(flat) > tol]

    # Get common multiple using inverse of GCD-like method
    lcm_factor = np.lcm.reduce(np.round(1 / nonzero).astype(int))
    
    # Try a safer way by scaling until everything is close to int
    scale = 1
    for s in range(1, 1000):
        test = np.round(vectors * s)
        if np.allclose(vectors * s, test, atol=tol):
            scale = s
            break

    scaled = np.round(vectors * scale).astype(int)

    # Optionally reduce each vector by its GCD
    def reduce_gcd(vec):
        gcd = np.gcd.reduce(vec[vec != 0]) if np.any(vec != 0) else 1
        return vec // gcd

    return np.array([reduce_gcd(v) for v in scaled])


def miller_bravais_to_cartesian(hkil, a=1.0, c=1.0):
    h, k, i, l = hkil
    assert h + k + i == 0, "Invalid Miller-Bravais indices: h + k + i must be 0"

    x = a * (3/2) * h
    y = a * (np.sqrt(3)/2) * h + a * np.sqrt(3) * k
    z = c * l
    return np.array([x, y, z])

def cartesian_to_hkil(vector, a=1.0, c=1.0, tol=1e-6):
    x, y, z = vector

    # Hexagonal lattice vectors for directions
    a1 = np.array([a, 0, 0])
    a2 = np.array([a/2, a*np.sqrt(3)/2, 0])
    c_vec = np.array([0, 0, c])

    # Build inverse matrix
    B = np.column_stack([a1, a2, c_vec])  # 3x3 matrix
    coeffs = np.linalg.solve(B, vector)  # solves B @ [u, v, l] = vector

    u, v, l = coeffs
    u_r = int(np.round(u))
    v_r = int(np.round(v))
    w_r = -u_r - v_r
    l_r = int(np.round(l))

    # Check error
    recon = miller_bravais_to_cartesian([u_r, v_r, w_r, l_r], a, c)
    residual = np.linalg.norm(vector - recon)

    if residual > tol:
        raise ValueError(f"Conversion error too large: residual = {residual:.3e}")

    return [u_r, v_r, w_r, l_r]

def format_miller_indices(array):
    array = np.round(array).astype(int)

    formatted = []
    for vec in array:
        parts = []
        for n in vec:
            if n < 0:
                parts.append(f"-{abs(n)}")
            else:
                parts.append(str(n))
        formatted.append("".join(parts))
    
    return (" ".join(formatted))



def gamma_struct(
    structure: Union[str, Atoms],
    struct_format: str = 'vasp',
    interlayer_shift: float = 0.0,
    increments: int = 10,
    vacuum : float = 10.0,
    duplicate_z: int = 10,
    shift_plane: float = 0.5,
    write_file: bool = True,
    dft_workflow_root: Path = DefaultPath.DFT_WORKFLOW_ROOT,
    write_dir: Union[str, Path] = Path("GAMMA_SURF"),
) -> list[tuple[float, float, float]]:
    """
    Generalized Gamma surface POSCAR generator

    Parameters
    ----------
    structure : str or ASE Atoms
        CIF file path or ASE Atoms object.
    interlayer_shift : float
        z-coordinate of the sliding layer in Angstrom
    increments : float
        Number of increments along x and y displacement (default: 10).
    duplicate_z : int
        Number of repetitions along z.
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
    z_mid = interlayer_shift
    areaxy = np.linalg.norm(np.cross(lattice_vectors[0] , lattice_vectors[1]))    

    # Displacement ranges
    x_range = np.linspace(0,1,increments)
    y_range = np.linspace(0,1,increments)

    add_vacuum(atoms,vacuum=vacuum)
    for y in y_range:
        for x in x_range:
            atoms.set_positions(original_positions)
            burger = lattice_vectors[0] * x + lattice_vectors[1] * y
            for atom in atoms:
                if atom.position[2] > z_mid:
                    atom.position += burger
            atoms.wrap()
            
            
            if write_file:
                out_file =  write_dir / f"poscar_x_{round(x, 2)}_y_{round(y, 2)}.vasp"
                out_file.parent.mkdir(parents=True, exist_ok=True)  # Make sure the directory exists
                io.write(out_file, atoms, format=struct_format)



def gsfe_struct(
    structure: Union[str, Atoms],
    struct_format: str = 'vasp',
    interlayer_shift: float = 0.0,
    increments: int = 10,
    vacuum :float = 10.0,
    xmax : float = 1.0,
    ymax : float = 0.0,
    duplicate_z: int = 10,
    freeze_layers : int = 6,
    write_file: bool = True,
    dft_workflow_root: Path = DefaultPath.DFT_WORKFLOW_ROOT,
    write_dir: Union[str, Path] = Path('GSFE'),
) -> list[tuple[float, float, float]]:
    """
    Generalized GSFE POSCAR generator

    Parameters
    ----------
    structure : str or ASE Atoms
        CIF file path or ASE Atoms object.
    interlayer_shift : float
        The z-coordinate for the sliding layer.
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
    write_dir = Path(write_dir)
    # Load structure
    atoms = io.read(structure,format = struct_format) if isinstance(structure, (str,Path)) else structure.copy()
    atoms = atoms.repeat((1, 1, duplicate_z))
    original_positions = atoms.get_positions()
    lattice_vectors = atoms.cell
    
    atoms = apply_selectivedynamics(
        structure= atoms, 
        freeze_top_n_layers=freeze_layers, 
        freeze_bottom_n_layers=freeze_layers, 
        relax_flag='F F T') 

    z_mid = interlayer_shift
    areaxy = np.linalg.norm(np.cross(lattice_vectors[0] , lattice_vectors[1]))    

    # Displacement ranges
    x_range = np.linspace(0,xmax/np.linalg.norm(lattice_vectors[0]),increments)
    y_range = np.linspace(0,ymax/np.linalg.norm(lattice_vectors[1]),increments)

    add_vacuum(atoms,vacuum=vacuum)
    for x, y in zip(x_range,y_range):
        atoms.set_positions(original_positions)
        burger = lattice_vectors[0] * x + lattice_vectors[1] * y
        for atom in atoms:
            if atom.position[2] > z_mid:
                atom.position += burger
        atoms.wrap()

        if write_file:
            print(f"writing poscar_x_{round(x, 2)}_y_{round(y, 2)}.vasp in {write_dir}")
            out_file = write_dir / f"poscar_x_{round(x, 2)}_y_{round(y, 2)}.vasp"
            out_file.parent.mkdir(parents=True, exist_ok=True)  # Make sure the directory exists
            io.write(out_file, atoms,format=struct_format)
    return out_file

def create_needle(
    structure: Union[str, Atoms],
    struct_format: str = 'vasp',
    duplicate_z: int = 10,
    vacuum :float = 10.0,
    write_file: bool = True,
    write_dir: Union[str, Path] = Path('INI'),
):
    atoms = io.read(structure,format = struct_format) if isinstance(structure, (str,Path)) else structure.copy()
    atoms = atoms.repeat((1, 1, duplicate_z))
    add_vacuum(atoms,vacuum=vacuum)
    
    if write_file:
        out_file = write_dir / f"POSCAR"
        out_file.parent.mkdir(parents=True, exist_ok=True)  # Make sure the directory exists
        io.write(out_file, atoms,format=struct_format)
    return atoms


def apply_selectivedynamics(structure, freeze_top_n_layers, freeze_bottom_n_layers, relax_flag):
    
    if isinstance(structure, (str,Path)):
        poscar = POSCAR(structure) 
    else:
        write('POSCAR', structure, format='vasp', vasp5=True, direct=False)
        poscar = POSCAR('POSCAR')
    
    poscar.pos_data, poscar.fixation_data  = poscar.setup_layers(freeze_top_n_layers, freeze_bottom_n_layers, relax_flag)
    poscar.write_out('POSCAR')
    atoms = read("POSCAR", format="vasp")
    return atoms
    
