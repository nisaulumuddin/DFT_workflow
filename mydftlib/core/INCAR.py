from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy
from pathlib import Path

class INCAR:
        def __init__(self,filepath=None):

                if filepath is None:
                        self.tags = INCAR.default_tags()
                        return  # prevent further execution
                
                with open(filepath,'r') as f:
                        self.file_data = [i.split() for i in f.read().splitlines()]
                

                for fields in self.file_data:
                        if len(fields) != 0:
                                try:
                                        if len(" ".join(fields).split("!")[0]) != 0:
                                                line = " ".join(fields).split("!")[0]
                                                tagname = line.split("=")[0].split()[0]
                                                tagdat = line.split("=")[1].split()
                                                self.tags[tagname] = tagdat

                                except:
                                        continue

                

        def delete_tag(self, tag):
                if tag in self.tags:
                        del self.tags[tag]
                        print(f"Tag '{tag}' has been deleted.")
                else:
                        print(f"Tag '{tag}' does not exist.")
                        
        def set_tag(self,tag,value):
                self.tags[tag] = [value]         

        def get_sorted_tags_string(self):
                taglist = '\n'.join([str(tag) + ' = ' + ' '.join([str(i) for i in self.tags[tag]]) for tag in self.tags])
                return taglist

        def write_out(self,filename, folder=None):
        
                if folder is None:
                        folder = Path('.')  # current directory
                else:
                        folder = Path(folder)
                        folder.mkdir(parents=True, exist_ok=True)
                
                filepath = folder / filename
                if os.path.isfile(filepath):
                        os.remove(filepath)

                with open(filepath,'a') as f:
                        f.write(self.get_sorted_tags_string())
        
        @classmethod
        def generate_default_incar(cls):
                obj = cls()
                obj.tags = cls.default_tags()
                return obj

        @staticmethod
        def default_tags():
                return {
                        'ISTART': ['0'],
                        'ICHARG': ['2'],
                        'ENCUT': ['550'],
                        'PREC': ['Accurate'],
                        'ISIF': ['2'],
                        'EDIFF': ['1E-6'],
                        'EDIFFG': ['-0.01'],
                        'ISMEAR': ['1'],
                        'SIGMA': ['0.1'],
                        'IBRION': ['2'],
                        'NSW': ['100'],
                        'POTIM': ['0.1'],
                        'LREAL': ['Auto'],
                        'ALGO': ['Fast'],
                        'LWAVE': ['F'],
                        'LCHARG': ['F'],
                        'LVTOT': ['F'],
                        'LORBIT': ['10'],
                        'NELM': ['200'],
                        'NELMIN': ['6'],
                        'LASPH': ['T'],
                        'ISPIN': ['1'],
                        'NPAR': ['16'],
                        'AMIN': ['0.01'],}