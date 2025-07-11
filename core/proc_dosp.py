from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy
import sys

# forgot what this is for
class proc_dosp:
        def __init__(self,filepath,*,is_str_not_file=False):

                if is_str_not_file:
                        file_data = [i.split() for i in filepath.splitlines()]
                else:
                        with open(filepath,'r') as f:
                                file_data = [i.split() for i in f.read().splitlines()]
                self.file_data = []
                for f in range(1,len(file_data)):
                    self.file_data.append([float(i) for i in file_data[f]])
                self.header = file_data[0]
                self.colnum = len(self.file_data[1])
                self.nedos = len(self.file_data)
                
                
        def energy(self):
            # Extracts the energy column from the dosp file 
            self.energy = np.zeros(self.nedos)
            count = 0
            for line in self.file_data[1:self.nedos]:
                self.energy[count] = line[0]
                count+=1
            return self.energy


        def select_pdos(self,index):
            # Select the column PDOS we are interested in
            pdos = np.zeros(self.nedos)
            count = 0
            for line in self.file_data[1:self.nedos]:
                pdos[count] = line[index]
                count+=1
            
            return pdos
             
        

        def nearest_value(self,user_energy, Energy):
            # Finds the nearest energy value and its index within the array, within the dosp file, given a user-specified energy 
            diff = [abs(user_energy - float(i)) for i in Energy]

            search_val_upp = user_energy + min(diff)
            search_val_low = user_energy - min(diff)

            if search_val_upp in Energy:
                nearE = deepcopy(search_val_upp)
                nearest_row = np.where(Energy == search_val_upp)

            elif search_val_low in Energy:
                nearE = deepcopy(search_val_low)
                nearest_row = np.where(Energy == search_val_low)

            return nearE, nearest_row
        
        def select_section(self,lower_row, upper_row,pdos):
            sect_energy = self.energy[lower_row:upper_row]
            sect_pdos = pdos[lower_row:upper_row]
            return sect_energy, sect_pdos

        def integrate_trapezoid(self,sect_energy,sect_pdos):
            # Calculates the area under the DOS curve using trapezoid method

            if (len(sect_energy) != len(sect_pdos)):
                print("Input error. The energy and PDOS arrays needs to be the same size")
                sys.exit()

            area = 0.0

            for i in range(len(sect_energy) - 1):
                area += 0.5 * (sect_pdos[i] + sect_pdos[i + 1]) * (sect_energy[i + 1] - sect_energy[i])

            return area 


                
    
