from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy

#This is the output file of vdos_split.pl. DOSCAR file split into projected DOS. 
class dosp:
        def __init__(self,filepath,*,is_str_not_file=False):

                if is_str_not_file:
                        self.file_data = [i.split() for i in filepath.splitlines()]
                else:
                        with open(filepath,'r') as f:
                                self.file_data = [i.split() for i in f.read().splitlines()]

                self.colnum = len(self.file_data[1])
                self.nedos = len(self.file_data)
                self.case = 0
                self.energy = [] 
                self.noncollinear = 'n' #str(input('is this calculation non-collinear? y/n:   '))
 
                if self.colnum <=5:
                        self.case = 1 #ispin=1,lorbit=10
                        self.s = []
                        self.p = []
                        self.d = []
                        self.f = []
                        for line in self.file_data[0:self.nedos]:
                                self.energy.append(line[0])
                                if self.colnum >= 5:
                                        self.f.append(line[4])
                                if self.colnum >= 4:
                                        self.d.append(line[3])
                                if self.colnum >= 3:
                                        self.p.append(line[2])
                                if self.colnum >= 2:
                                        self.s.append(line[1])
                                else:
                                        print('dosp file cannot be read')
                


                elif self.colnum <= 9:
                        self.case = 2 #ispin=2,lorbit=10
                        self.s_u = []
                        self.s_d = []
                        self.p_u = []
                        self.p_d = []
                        self.d_u = []
                        self.d_d = []
                        self.f_u = []
                        self.f_d = []
                        for line in self.file_data[0:self.nedos]:
                                self.energy.append(line[0])
                                if self.colnum >= 9:
                                        self.f_u.append(line[7]) 
                                        self.f_d.append(line[8])
                                if self.colnum >= 7:
                                        self.d_u.append(line[5])
                                        self.d_d.append(line[6]) 
                                if self.colnum >= 5:
                                        self.p_u.append(line[3])
                                        self.p_d.append(line[4])
                                if self.colnum >= 3:
                                        self.s_u.append(line[1])
                                        self.s_d.append(line[2])

                elif self.colnum <= 17 and self.noncollinear == 'n':
                        self.case = 3 #ispin=1,lorbit=11

                elif self.colnum <=17 and self.noncollinear == 'y':
                        self.case = 5 #ispin=1, lorbit = 10, non-collinear 
                        self.s = []
                        self.smx = []
                        self.smy = []
                        self.smz = []
                        self.p = []
                        self.pmx = []
                        self.pmy = []
                        self.pmz = []
                        self.d = []
                        self.dmx = []
                        self.dmy = []
                        self.dmz = []
                        self.f = []
                        self.fmx = []
                        self.fmy = []
                        self.fmz = []
                        for line in self.file_data[0:self.nedos]:
                                self.energy.append(line[0])
                                if self.colnum >= 5: 
                                        self.s.append(line[1])
                                        self.smx.append(line[2])
                                        self.smy.append(line[3])
                                        self.smz.append(line[4])
                                if self.colnum >= 9: 
                                        self.p.append(line[5])
                                        self.pmx.append(line[6])
                                        self.pmy.append(line[7])
                                        self.pmz.append(line[8])
                                if self.colnum >= 13:
                                        self.d.append(line[9])
                                        self.dmx.append(line[10])
                                        self.dmy.append(line[11])
                                        self.dmz.append(line[12])
                                if self.colnum >= 17:
                                        self.f.append(line[13])
                                        self.fmx.append(line[14])
                                        self.fmy.append(line[15])
                                        self.fmz.append(line[16])


                        
                elif self.colnum <= 33:
                        self.case = 4 #ispin=2,lorbit=11, collinear
                        self.s = []
                        self.s_u = []
                        self.s_d = []
                        self.py_u  = []
                        self.py_d = []
                        self.pz_u = []
                        self.pz_d = []
                        self.px_u = []
                        self.px_d = []
                        self.dxy_u = []
                        self.dxy_d = []
                        self.dyz_u = []
                        self.dyz_d = []
                        self.dz2_r2_u = []
                        self.dz2_r2_d = []
                        self.dxz_u = []
                        self.dxz_d = []
                        self.dx2_y2_u = []
                        self.dx2_y2_d = []
                        self.energy = []
                        self.p = []
                        self.d = []
                        self.f = [] # 14 columns
                        # for the f orbital, not sure the sequence in which it is decomposed

                        for line in self.file_data[0:self.nedos]:
                                self.energy.append(line[0])

                                if self.colnum >= 33:
                                        self.f.append(line[19])

                                if self.colnum>=19:
                                        self.dxy_u.append(line[9])
                                        self.dxy_d.append(line[10])
                                        self.dyz_u.append(line[11])
                                        self.dyz_d.append(line[12])
                                        self.dz2_r2_u.append(line[13])
                                        self.dz2_r2_d.append(line[14])
                                        self.dxz_u.append(line[15])
                                        self.dxz_d.append(line[16])
                                        self.dx2_y2_u.append(line[17])
                                        self.dx2_y2_d.append(line[18])

                                if self.colnum >= 9:
                                        self.py_u.append(line[3])
                                        self.py_d.append(line[4])
                                        self.pz_u.append(line[5])
                                        self.pz_d.append(line[6])
                                        self.px_u.append(line[7])
                                        self.px_d.append(line[8])

                                if self.colnum >=3:
                                        self.s_u.append(line[1])
                                        self.s_d.append(line[2])


                else:
                        print("case not covered")