#!/usr/bin/env python3

from mydftlib.config.config import Config as DefaultPath
from mydftlib.structures import GSFE
from mydftlib.core.POTCAR import create_POTCAR_from_poscar
from mydftlib.core.INCAR import INCAR
from mydftlib.core.KPOINTS import KPOINTS
from mydftlib.utilities.setupcalc import distribute_vasp_inputs, wipe_to_raw
import glob
from pathlib import Path

poscar = 'CONTCAR_Fe7Ta6_x1y1z5_nm'
folder = Path('GSFE_dir11-20')

# create POSCAR 
gsfe_results = GSFE.gsfe_struct(
    structure= poscar,
    struct_format = 'vasp',
    interlayer_shift = 70.4,
    increments= 10,
    vacuum =0.0,
    xmax = 4.8576272320999996,
    ymax = 0.0,
    duplicate_z = 1,
    freeze_layers= 6,
    write_dir= folder,
)

# create POTCAR 
create_POTCAR_from_poscar( 'POSCAR', target_dir=folder)

# create INCAR 
incar_obj = INCAR.generate_default_incar() 
# user needs to set these up if applicable. else, remove. 
incar_obj.set_tag('NPAR','16')
incar_obj.set_tag('AMIN','0.01')
incar_obj.write_out('INCAR', folder=folder)
print(incar_obj)
# create KPOINTS 
kpoints_obj = KPOINTS.generate_default_kpoints()
kpoints_obj.set_centering('Gamma')
kpoints_obj.set_grid('5 5 1')
kpoints_obj.write_out('KPOINTS',folder=folder)


# creating paths 
distribute_vasp_inputs(
    incar= 'INCAR',
    kpoints= 'KPOINTS',
    potcar= 'POTCAR',
    poscar_glob = "poscar_x_*",
    output_root= folder
)
created_poscar_files = [p.name for p in Path(folder).glob("poscar_x_*")]
print(created_poscar_files)

wipe_to_raw('INCAR','KPOINTS','POTCAR' ,*created_poscar_files, raw_dir= "raw_inputs", src_dir=folder)
