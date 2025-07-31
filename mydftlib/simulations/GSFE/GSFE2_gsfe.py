#!/usr/bin/env python3

from mydftlib.config.config import Config as DefaultPath
from mydftlib.structures import GSFE
from mydftlib.core.POTCAR import create_POTCAR
from mydftlib.core.INCAR import INCAR
from mydftlib.core.KPOINTS import KPOINTS
from mydftlib.utilities.setupcalc import distribute_vasp_inputs, wipe_to_raw
import glob

poscar = 'POSCAR'
folder = 'GSFE'

# create POSCAR 
gsfe_results = GSFE.gsfe_struct(
    structure=poscar,
    struct_format='vasp',  
    interlayer_shift=60.6888,
    increments=10,
    xmax=4.8487695842000003/2,
    duplicate_z=1,
)

# create POTCAR 
create_POTCAR(gsfe_results)

# create INCAR 
incar_obj = INCAR.generate_default_incar() 
# user needs to set these up if applicable. else, remove. 
incar_obj.set_tag('ISPIN','2')
incar_obj.set_tag('MAGMOM','18*0  21*3  18*0  21*3  18*0  21*3  18*0  21*3')
incar_obj.set_tag('NPAR','16')
incar_obj.set_tag('AMIN','0.01')
incar_obj.write_out('INCAR')

# create KPOINTS 
kpoints_obj = KPOINTS.generate_default_kpoints()
kpoints_obj.set_centering('Gamma')
kpoints_obj.set_grid('5 5 1')
kpoints_obj.write_out('KPOINTS')

# creating paths 
distribute_vasp_inputs(
    incar= 'INCAR',
    kpoints= 'KPOINTS',
    potcar= 'POTCAR',
    poscar_glob = "poscar_x_*",
    output_root= folder
)

created_poscar_files = glob.glob("poscar_x_*")

wipe_to_raw('INCAR','KPOINTS','POSCAR','POTCAR',*created_poscar_files)
