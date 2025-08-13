from __future__ import annotations

import os
import subprocess
from pathlib import Path
import yaml
from mydftlib.config.config import Config as DefaultPath
from ase.io import read,write
# from mydftlib.config import hexagonal_orth_slip_systems
from typing import Union 

import yaml
from pathlib import Path

# other upcoming functions: import from MaterialsProject

def create_structure_set(
    slip_systems_yaml : Union[str,list] = 'hexagonal_orth_slip_systems.yaml',
    lattice_type: str = None,
    lattice_constants : str = None,
    element_list : str = None,
    suffix_name : str = None,
    duplicates_xyz : str = '1 1 1',
    orthogonal : bool = False, 
):
    if isinstance(slip_systems_yaml,str):
        with open(slip_systems_yaml, "r") as f:
            slip_systems = yaml.safe_load(f)
    elif isinstance(slip_systems_yaml,list):
        slip_systems = slip_systems_yaml
    else:
        raise ValueError("slip_systems_yaml must be a file path or a loaded YAML list ")
    
    # Access or loop through slip systems
    for system in slip_systems:
        print(
        f"{ system['system']} → Plane: {''.join(str(x) for x in system['plane']) }, "
        f" Direction1: {''.join(str(x) for x in system['direction1']) }"
        f" Direction2: {''.join(str(x) for x in system['direction2']) }"
    )

    # Define the input command for Atomsk
    input_command = f"{lattice_type} {lattice_constants} {element_list} "  #
    duplicate_x = duplicates_xyz.split()[0]
    duplicate_y = duplicates_xyz.split()[1]
    duplicate_z = duplicates_xyz.split()[2]

    for system in slip_systems:
        direction1 = ''.join(str(x) for x in system['direction1'])
        direction2 = ''.join(str(x) for x in system['direction2'])
        direction3 = ''.join(str(x) for x in system['plane'])
        
        # Optionally specify output filename
        output_file = f"{suffix_name}_x{direction1}_y{direction2}_z{direction3}_x{duplicate_x}y{duplicate_y}z{duplicate_z}.lmp"
        vasp_file =f"{suffix_name}_x{direction1}_y{direction2}_z{direction3}_x{duplicate_x}y{duplicate_y}z{duplicate_z}.vasp"
        
        folder = 'atomsk_cells'
        output_root = Path(folder)
        output_root.mkdir(parents=True, exist_ok=True)
        
        # Call the function
        result_path = run_atomsk(
            input_command=input_command,
            orientation=f"[{direction1}] [{direction2}] [{direction3}]",
            output_filename=output_file,
            duplication = duplicates_xyz,
            subfolder = 'lmps',
            root = output_root,
            orthogonal = orthogonal,
        )
        
        # write to vasp folder
        vasp_subfolder = output_root / "vasp"
        vasp_subfolder.mkdir(parents=True,exist_ok=True)
        if result_path is None or not os.path.isfile(result_path):
            print(f"Output file not found: {result_path}")
            continue
        else:
            atoms = read(result_path,format='lammps-data')
            write(vasp_subfolder / vasp_file,images=atoms)
            print("written as vasp file:", vasp_file)
    return

# mydftlib/utils/loader.py

import yaml
from pathlib import Path

def load_hexagonal_slip_systems():
    # Find root path of the mydftlib package
    root_path = Path(__file__).resolve().parents[1]  # points to mydftlib/
    yaml_path = root_path / "config" / "hexagonal_orth_slip_systems.yaml"

    with open(yaml_path, "r") as f:
        return yaml.safe_load(f)


def run_atomsk(
    input_command: str,
    atomsk_path: str = "/home/wz300646/Installations/atomsk/src/atomsk",
    root: Union[str,Path] = Path('.'),
    orientation: str = "[100] [010] [001]",
    output_filename: str = "reoriented.cif",
    subfolder: str = "rotated_cells",
    duplication: str = "1 1 1",
    orthogonal : bool = False,
) -> Path:
    """
    Run Atomsk to generate a structure with given parameters and save it to a subfolder under DFT_Workflow.

    Parameters
    ----------
    input_command : str
        Atomsk input command string (e.g., "fcc 3.6 Cu").
    orientation : str
        Orientation string for Atomsk (e.g., "[100] [010] [001]") (default: "[100] [010] [001]").
    output_filename : str (default: "reoriented.cif")
        The output file name (e.g., "Cu_oriented.lmp").
    atomsk_path : str
        Full path to the Atomsk binary (provided externally).
    dft_workflow_root : Path
        The root path to your DFT_Workflow directory.
    subfolder : str
        Name of the subdirectory inside DFT_Workflow to store output (default: "oriented_cells").
    duplication : str
        Duplication factor string (default: "1 1 10").

    Returns
    -------
    Path
        Full path to the generated structure file.
    """
    # Resolve output directory inside DFT_Workflow
    root = Path(root)
    output_dir = root / subfolder
    output_dir.mkdir(parents=True, exist_ok=True)

    if orthogonal is True:
        orthogonal = '-orthogonal-cell'
    else:
        orthogonal = ''
    # Build command
    create_cmd = f"{atomsk_path} --create {input_command} orient {orientation} {orthogonal} -duplicate {duplication} {output_filename}"
    move_cmd = f"mv {output_filename} {output_dir}"

    try:
        subprocess.run(create_cmd, shell=True, check=True, input='n\n', text=True)
        subprocess.run(move_cmd, shell=True, check=True)
        print(f"Atomsk ran successfully. File saved to: {output_dir / output_filename}")
    except subprocess.CalledProcessError as e:
        print(f"Atomsk failed: {e}")
        return None

    return output_dir / output_filename
