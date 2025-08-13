from __future__ import annotations
from pathlib import Path
import shutil
import os
from typing import Union

def wipe_to_raw(*filenames, raw_dir="raw_inputs", src_dir="."):
    """
    Copies specified VASP input files into a raw input directory.

    Parameters:
    - *filenames (str): Any number of filenames (e.g., 'INCAR', 'POSCAR')
    - raw_dir (str or Path): Target directory to store raw copies.
    - src_dir (str or Path): Source directory where the files are located.
    """
    
    src_dir = Path(src_dir)
    raw_dir = src_dir / raw_dir
    raw_dir.mkdir(parents=True, exist_ok=True)

    for name in filenames:
        src_file = src_dir / name
        dst_file = raw_dir / name
        if src_file.exists():
            shutil.move(src_file, dst_file)
            print(f"[✓] Copied {src_file} → {dst_file}")
        else:
            print(f"[!] Warning: {src_file} not found.")
            



def distribute_vasp_inputs(
    incar: str,
    kpoints: str,
    potcar: str,
    poscar_glob: str,
    output_root: Path | str,
):
    """
    Distribute VASP input files to directories corresponding to POSCAR variants.

    Parameters:
    - incar: Path to INCAR file
    - kpoints: Path to KPOINTS file
    - potcar: Path to POTCAR file
    - poscar_glob: Glob pattern to match POSCAR variants (e.g. 'poscar_x_*')
    - output_root: Path to root directory where subfolders will be created
    """
        
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    # Find all matching POSCAR variants
    for poscar_path in output_root.glob(poscar_glob):
        folder_name = poscar_path.name.replace('poscar_', '')  # e.g. poscar_x_001 → x_001
        target_dir = output_root  / folder_name
        target_dir.mkdir(parents=True, exist_ok=True)

        # Copy shared files
        shutil.copy(output_root / incar, target_dir / 'INCAR')
        shutil.copy(output_root / kpoints, target_dir / 'KPOINTS')
        shutil.copy(output_root /  potcar, target_dir / 'POTCAR')

        # Rename and copy POSCAR variant
        shutil.copy(poscar_path, target_dir / 'POSCAR')
        print(f"copied INCAR, KPOINTS, POTCAR, POSCAR to  {target_dir}")
