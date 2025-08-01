from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy

from pathlib import Path
from mydftlib.config.config import Config as DefaultPath
from mydftlib.core.POSCAR import POSCAR

class POTCAR:
    all_elements = [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
        'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
        'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re',
        'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
        'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm'
    ]

    std_potcar = [
        'H', 'He', 'Li_sv', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na_pv', 'Mg', 'Al', 'Si',
        'P', 'S', 'Cl', 'Ar', 'K_sv', 'Ca_sv', 'Sc_sv', 'Ti_sv', 'V_sv', 'Cr_pv', 'Mn_pv', 'Fe', 'Co', 'Ni',
        'Cu', 'Zn', 'Ga_d', 'Ge_d', 'As', 'Se', 'Br', 'Kr', 'Rb_sv', 'Sr_sv', 'Y_sv', 'Zr_sv', 'Nb_sv',
        'Mo_sv', 'Tc_pv', 'Ru_pv', 'Rh_pv', 'Pd', 'Ag', 'Cd', 'In_d', 'Sn_d', 'Sb', 'Te', 'I', 'Xe',
        'Cs_sv', 'Ba_sv', 'La', 'Ce', 'Pr_3', 'Nd_3', 'Pm_3', 'Sm_3', 'Eu_2', 'Gd_3', 'Tb_3', 'Dy_3',
        'Ho_3', 'Er_3', 'Tm_3', 'Yb_2', 'Lu_3', 'Hf_pv', 'Ta_pv', 'W_sv', 'Re', 'Os', 'Ir', 'Pt', 'Au',
        'Hg', 'Tl_d', 'Pb_d', 'Bi_d', 'Po_d', 'At', 'Rn', 'Fr_sv', 'Ra_sv', 'Ac', 'Th', 'Pa', 'U', 'Np',
        'Pu', 'Am', 'Cm'
    ]

    potcar_map = dict(zip(all_elements, std_potcar))

    def __init__(self, symbols, potcar_root):
        self.symbols = symbols
        self.potcar_root = Path(potcar_root)
        self.mapped = [self.potcar_map[symbol] for symbol in symbols]
        self._binary_data = []

    @classmethod
    def from_poscar(cls, poscar_obj, potcar_root):
        try:
            symbols = poscar_obj.atom_header[0]
        except AttributeError:
            raise ValueError("POSCAR object must have `atom_header[0]` containing element symbols.")
        return cls(symbols, potcar_root)

    def read_potcars(self):
        self._binary_data = []
        for name in self.mapped:
            potcar_path = self.potcar_root / name / "POTCAR"
            if potcar_path.exists():
                with potcar_path.open("rb") as f:
                    self._binary_data.append(f.read())
            else:
                print(f"Warning: {potcar_path} not found.")
                self._binary_data.append(b'')

    def write_out(self, output_path="POTCAR"):
        if not self._binary_data:
            self.read_potcars()
        with open(output_path, "wb") as f:
            for data in self._binary_data:
                f.write(data)

def create_POTCAR_from_poscar(poscar_path, target_dir, potcar_root=DefaultPath.POTCAR_ROOT):
    """
    Creates a POTCAR file in `target_dir` using the element list from `poscar_path`.
    
    Parameters:
    - poscar_path (Path or str): Path to POSCAR file.
    - target_dir (Path or str): Directory to write POTCAR into.
    - potcar_root (Path or str): Root directory where pseudopotential folders are stored.
    """
    # Load POSCAR
    poscar = POSCAR(poscar_path)
    

    # Create POTCAR object from POSCAR
    potcar = POTCAR.from_poscar(poscar, potcar_root=potcar_root)

    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    # Write to target directory
    output_path = Path(target_dir) / 'POTCAR'

    
    potcar.write_out(output_path)
    return output_path