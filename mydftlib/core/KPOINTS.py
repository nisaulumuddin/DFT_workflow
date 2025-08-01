from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy
from pathlib import Path

class KPOINTS:
    def __init__(self, filepath=None, *, is_str_not_file=False):
        if filepath is None:
            self.kpoints_data = [
                ["Automatic mesh"],
                ["0"],
                ["Gamma"],  # default centering
                ["1", "1", "1"],  # default grid
                ["0", "0", "0"]  # default shift
            ]
        else:
            if is_str_not_file:
                self.kpoints_data = [line.split() for line in filepath.splitlines()]
            else:
                with open(filepath, 'r') as f:
                    self.kpoints_data = [line.split() for line in f.read().splitlines()]

    @classmethod
    def generate_default_kpoints(cls):
        return cls()

    def set_centering(self, centering):
        centering = centering.strip().capitalize()
        if centering not in ["Gamma", "Monkhorst-pack"]:
            raise ValueError("Centering must be 'Gamma' or 'Monkhorst-Pack'")
        self.kpoints_data[2] = [centering]

    def set_grid(self, grid_string):
        grid_parts = grid_string.strip().split()
        if len(grid_parts) != 3:
            raise ValueError("Grid must be a string like '5 5 1'")
        self.kpoints_data[3] = grid_parts

    def set_shift(self, shift_string="0 0 0"):
        shift_parts = shift_string.strip().split()
        if len(shift_parts) != 3:
            raise ValueError("Shift must be a string like '0 0 0'")
        self.kpoints_data[4] = shift_parts


    def write_out(self, filename="KPOINTS", folder="."):
        """
        Write the KPOINTS data to a file.

        Parameters:
        - filename (str): Name of the file to write (default: "KPOINTS")
        - folder (str or Path): Output directory (default: current directory)
        """
        
        if folder is None:
            folder = Path('.')  # current directory
        else:
            folder = Path(folder)
            folder.mkdir(parents=True, exist_ok=True)

        filepath = folder / filename

        with open(filepath, 'w') as f:
            for line in self.kpoints_data:
                f.write(" ".join(line) + "\n")


