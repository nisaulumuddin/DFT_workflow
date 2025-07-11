from __future__ import annotations

import os
import subprocess
from pathlib import Path

from config.config import Config as DefaultPath

# other upcoming functions: import from MaterialsProject

def run_atomsk(
    input_command: str,
    atomsk_path: str,
    dft_workflow_root: Path = DefaultPath.DFT_WORKFLOW_ROOT,
    orientation: str = "[100] [010] [001]",
    output_filename: str = "reoriented.cif",
    subfolder: str = "oriented_cells",
    duplication: str = "1 1 1"
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
    output_dir = dft_workflow_root / subfolder
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build command
    create_cmd = f"{atomsk_path} --create {input_command} orient {orientation} -duplicate {duplication} {output_filename}"
    move_cmd = f"mv {output_filename} {output_dir}"

    try:
        subprocess.run(create_cmd, shell=True, check=True)
        subprocess.run(move_cmd, shell=True, check=True)
        print(f"Atomsk ran successfully. File saved to: {output_dir / output_filename}")
    except subprocess.CalledProcessError as e:
        print(f"Atomsk failed: {e}")
        return None

    return output_dir / output_filename
