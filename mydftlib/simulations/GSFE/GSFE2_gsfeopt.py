#!/usr/bin/env python3
import yaml
from pathlib import Path

from mydftlib.config.config import Config as DefaultPath
from mydftlib.structures import GSFE
from mydftlib.core.POTCAR import create_POTCAR_from_poscar
from mydftlib.core.INCAR import INCAR
from mydftlib.core.KPOINTS import KPOINTS
from mydftlib.utilities.setupcalc import distribute_vasp_inputs, wipe_to_raw

def main(config_file):
    cfg = yaml.safe_load(open(config_file))

    poscar = cfg["poscar"]
    folder = Path(cfg["folder"])
    folder.mkdir(exist_ok=True)

    # --- Step 1: Create GSFE POSCAR files ---
    gsfe_results = GSFE.gsfe_struct(
        structure=poscar,
        struct_format=cfg["gsfe"]["struct_format"],
        interlayer_shift=cfg["gsfe"]["interlayer_shift"],
        increments=cfg["gsfe"]["increments"],
        vacuum=cfg["gsfe"]["vacuum"],
        xmax=cfg["gsfe"]["xmax"],
        ymax=cfg["gsfe"]["ymax"],
        duplicate_z=cfg["gsfe"]["duplicate_z"],
        freeze_layers=cfg["gsfe"]["freeze_layers"],
        write_dir=folder,
    )

    # --- Step 2: Create POTCAR ---
    create_POTCAR_from_poscar('POSCAR', target_dir=folder)

    # --- Step 3: Create INCAR (fully generic) ---
    incar_obj = INCAR.generate_default_incar()
    for tag, value in cfg.get("incar", {}).items():
        incar_obj.set_tag(tag, str(value))
    incar_obj.write_out("INCAR", folder=folder)
    print("INCAR created with tags:", cfg.get("incar", {}))

    # --- Step 4: Create KPOINTS (generic) ---
    kpoints_obj = KPOINTS.generate_default_kpoints()
    if "centering" in cfg.get("kpoints", {}):
        kpoints_obj.set_centering(cfg["kpoints"]["centering"])
    if "grid" in cfg.get("kpoints", {}):
        kpoints_obj.set_grid(cfg["kpoints"]["grid"])
    kpoints_obj.write_out("KPOINTS", folder=folder)

    # --- Step 5: Distribute files for each GSFE structure ---
    distribute_vasp_inputs(
        incar='INCAR',
        kpoints='KPOINTS',
        potcar='POTCAR',
        poscar_glob=cfg["distribution"]["poscar_glob"],
        output_root=folder
    )

    created_poscar_files = [p.name for p in Path(folder).glob(cfg["distribution"]["poscar_glob"])]
    print("Created POSCAR variants:", created_poscar_files)

    # --- Step 6: Wipe to raw directory ---
    wipe_to_raw(
        'INCAR', 'KPOINTS', 'POTCAR', *created_poscar_files,
        raw_dir=cfg["distribution"]["raw_dir"], src_dir=folder
    )

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: ./GSFE2_gsfe.py config.yaml")
        sys.exit(1)
    main(sys.argv[1])
