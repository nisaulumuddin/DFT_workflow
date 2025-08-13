#!/usr/bin/env python3
from __future__ import annotations

# === Standard Library ===
import os
from pathlib import Path
from typing import Union

# === Third-Party Libraries ===
import numpy as np
from ase import io, Atoms
from ase.io import write

# === Local Modules ===
from mydftlib.config.config import Config as DefaultPath
from mydftlib.core.POSCAR import POSCAR
from mydftlib.core.INCAR import INCAR
from mydftlib.core.KPOINTS import KPOINTS
from mydftlib.core.POTCAR import create_POTCAR_from_poscar
from mydftlib.structures.GSFE import create_needle, apply_selectivedynamics
from mydftlib.utilities.setupcalc import  wipe_to_raw
import yaml
from pathlib import Path

def main(
    config_file : str = None,
    ):
    
    cfg = yaml.safe_load(open(config_file))
    ini_dir = Path(cfg["folder"])
    ini_dir.mkdir(exist_ok=True)

    atoms = create_needle(
        structure=cfg['structure'],
        struct_format='vasp',
        duplicate_z=cfg['duplicate_z'],
        vacuum=cfg['vacuum'],
        write_dir=ini_dir,
        write_file=False,
    )

    atoms = apply_selectivedynamics(
        structure=atoms, 
        freeze_top_n_layers=cfg['freeze_top_n_layers'], 
        freeze_bottom_n_layers=cfg['freeze_bottom_n_layers'], 
        relax_flag=cfg['relax_flag']
    )
    poscar_path = ini_dir / 'POSCAR'
    write(poscar_path, atoms, format='vasp', vasp5=True, direct=False)

    incar = INCAR.generate_default_incar()
    
    for tag, value in cfg["incar"].items():
        incar.set_tag(tag, str(value))
    incar.write_out("INCAR", folder=ini_dir)

    create_POTCAR_from_poscar(poscar_path=ini_dir / 'POSCAR', target_dir=ini_dir)

    kpoints = KPOINTS.generate_default_kpoints()
    kpoints.set_centering(cfg['kcentering'])
    kpoints.set_grid(cfg['kgrid'])
    kpoints.write_out('KPOINTS', folder=ini_dir)

    wipe_to_raw('INCAR','KPOINTS','POSCAR')
    print(f"Folder {ini_dir} is now ready for VASP")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: ./GSFE1_iniopt.py config.yaml")
        sys.exit(1)
    main(sys.argv[1])
