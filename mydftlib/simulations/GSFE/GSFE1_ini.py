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

def main():
    ini_dir = Path("INI")


    # Step 1: Generate POSCAR
    atoms = create_needle(
        structure='CONTCAR_Ta6Fe7_fo_unit',
        struct_format='vasp',
        duplicate_z=4,
        write_dir=ini_dir,
        write_file=False,
    )

    atoms = apply_selectivedynamics(
        structure= atoms, 
        freeze_top_n_layers=6, 
        freeze_bottom_n_layers=6, 
        relax_flag='T T T') 
    
    # Now write POSCAR with selective dynamics flags included
    write('INI/POSCAR', atoms, format='vasp', vasp5=True, direct=False)
    
    # Step 3: Create INCAR
    incar = INCAR.generate_default_incar()
    incar.set_tag('ISPIN', '1')
    # incar.set_tag('MAGMOM', '18*0  21*3  18*0  21*3  18*0  21*3  18*0  21*3')
    incar.set_tag('NPAR', '16')
    incar.write_out('INCAR', folder=ini_dir)

    # Step 4: Create POTCAR
    create_POTCAR_from_poscar(poscar_path=ini_dir / 'POSCAR', target_dir=ini_dir)

    # Step 5: Create KPOINTS
    kpoints = KPOINTS.generate_default_kpoints()
    kpoints.set_centering('Gamma')
    kpoints.set_grid('5 5 1')
    kpoints.write_out('KPOINTS', folder=ini_dir)

    print(f"Folder {ini_dir} is now ready for VASP")
    
    # Step 6: Tidy up files
    wipe_to_raw('INCAR','KPOINTS','POSCAR')
    
    
if __name__ == "__main__":
    main()
