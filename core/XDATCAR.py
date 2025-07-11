from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy

class XDATCAR:
        def __init__(self,filepath,*,is_str_not_file=False):

                if is_str_not_file:
                        self.file_data = [i.split() for i in filepath.splitlines()]
                else:
                        with open(filepath,'r') as f:
                                self.file_data = [i.split() for i in f.read().splitlines()]

                self.name = ' '.join(self.file_data[0])

                # Get the lattice parameter data (different frames may have updated lattice parameter)
                self.lattice_scalar = float(self.file_data[1][0])
                #self.factored_lattice_parameters = np.array(self.file_data[2:5]).astype(np.float64)

                self.factored_lattice_parameters = []
                for num,line in enumerate(self.file_data):
                        if 'configuration' in ' '.join(line):
                                self.factored_lattice_parameters.append(np.array(self.file_data[num-5:num-5+3]).astype(np.float64))
                #print(self.factored_lattice_parameters)
                #self.lattice_parameters = self.lattice_scalar * self.factored_lattice_parameters
                self.lattice_parameters = [self.lattice_scalar * i for i in self.factored_lattice_parameters]
                #print(self.lattice_parameters)
                #self.volume = np.linalg.det(self.lattice_parameters)
                self.volume = [np.linalg.det(i) for i in self.lattice_parameters]
                print(self.volume)

                # Get atom header data
                self.atom_header = [self.file_data[5],[int(i) for i in self.file_data[6]]]
                self.atom_list = []
                for speciesnum,species in enumerate(self.atom_header[0]):
                        for occur in range(self.atom_header[1][speciesnum]):
                                self.atom_list.append(species)

                self.atom_list = np.array(self.atom_list)

                # Get atomic positions by finding lines with 'configuration' in them (must be Direct, think this is normal anyhow)
                self.pos_data = []
                for num,line in enumerate(self.file_data):
                        if 'configuration' in ' '.join(line):
                                self.pos_data.append(self.file_data[num+1:num+len(self.atom_list)+1])

                self.pos_data = np.array(self.pos_data).astype(float)
                print(np.shape(self.pos_data))
        def export_poscar(self,filepath,frame):
                
                with open(filepath,'w') as f:
                        f.write(self.name + '\n')
                        f.write(str(self.lattice_scalar) + '\n')
                        for row in self.factored_lattice_parameters[frame]:
                                f.write('   '.join([f'{i:15.10f}' for i in row]) + '\n')

                        for row in self.atom_header:
                                f.write('  '.join([f'{i:>4}' for i in row]) + '\n')

                        f.write('Direct\n')
                        pos_frame = self.pos_data[frame]
                        
                        for rownum,row in enumerate(pos_frame):
                                f.write('  '.join([f'{i:20.13f}' for i in row]))
                                f.write('\n')
 

        # Loads the positions of a POSCAR into the XDATCAR (needs same lattice and number of atoms)
        def load_poscar_frame(self,poscar_object):
                if not len(self.pos_data) == 0:
                        self.pos_data = np.append(self.pos_data,np.array([poscar_object.get_dir_positions()]),axis=0)
                else:
                        self.pos_data = np.array([poscar_object.get_dir_positions()])

        # Scan through all frames of an atom looking for movements across boundaries and change the atom position if they exist
        # A direct coordinate movement threshold of 0.75 is set as a default
        def handle_periodicity(self,*,threshold=0.75):
                # For each atom in all frames (except the first one) and in all three directions
                for atomnum in range(len(self.atom_list)):
                        for frame_num,frame in enumerate(self.pos_data):
                                if not frame_num == 0:
                                        for direction_num,direction in enumerate(frame[atomnum]):
                                                # Calculate the movement of the atom in a direction between steps
                                                displacement = direction - self.pos_data[frame_num-1][atomnum][direction_num]

                                                # If the atom moved across the negative boundary
                                                if displacement > threshold:
                                                        self.pos_data[frame_num][atomnum][direction_num] -= 1

                                                # If the atom moved across the positive boundary
                                                elif displacement < -threshold:
                                                        self.pos_data[frame_num][atomnum][direction_num] += 1

        # Perform cubic interpolation of XDATCAR frames with a specified number of frames between.
        # Hermite smoothing means that the atoms smoothly stop at the original XDATCAR frames before moving on.
        # Hermite fixed frames is a list of boolean values for each frame in order on whether its derivatives are set to zero.
        # Note: This operation requires boundary-periodic interactions to be accounted for by creating direct coordinates less than 0 or greater than 1.
        def smooth_trajectories(self,frames_between,*,hermite_smooth=False,hermite_fixed_frames=[],periodicity_threshold=0.75):
                self.handle_periodicity(threshold=periodicity_threshold)

                # Convert all of the coordinates to Cartesian
                self.pos_data = np.array([np.matmul((self.lattice_scalar*self.factored_lattice_parameters),frame.T).T for frame in self.pos_data])

                # Rearrange the positions so that its atom -> frame -> position
                self.pos_data = np.array([self.pos_data[:,atom_num] for atom_num,atom in enumerate(self.atom_list)])

                def get_3D_spline(tdata,xyzdata):

                        xspline = si.CubicSpline(x=tdata,y=xyzdata[:,0])
                        yspline = si.CubicSpline(x=tdata,y=xyzdata[:,1])
                        zspline = si.CubicSpline(x=tdata,y=xyzdata[:,2])

                        if hermite_smooth:
                                # We want the velocities of all ions to be zero at the POSCARs
                                dudt = np.array([xspline(tdata,1),yspline(tdata,1),zspline(tdata,1)]).T

                                dudt[np.array(hermite_fixed_frames),:] = 0

                                xspline = si.CubicHermiteSpline(x=tdata,y=xyzdata[:,0],dydx=dudt[:,0])
                                yspline = si.CubicHermiteSpline(x=tdata,y=xyzdata[:,1],dydx=dudt[:,1])
                                zspline = si.CubicHermiteSpline(x=tdata,y=xyzdata[:,2],dydx=dudt[:,2])
                                return lambda t: np.array([xspline(t),yspline(t),zspline(t)])

                        else:
                                return lambda t: np.array([xspline(t),yspline(t),zspline(t)])


                frame_function = lambda t: np.array([get_3D_spline(tdata=np.linspace(0,1,endpoint=True,num=len(self.pos_data[0])),xyzdata=atom_pos_data)(t) for atom_pos_data in self.pos_data])
                number_of_frames = len(self.pos_data[0] - 1) * (frames_between + 1) + 1
                sample_time = np.linspace(start=0,stop=1,endpoint=True,num=number_of_frames)

                self.pos_data = np.array([frame_function(t) for t in sample_time])
                self.pos_data = np.array([np.matmul(np.linalg.inv(self.lattice_scalar*self.factored_lattice_parameters),frame.T).T for frame in self.pos_data])

        # Write out the XDATCAR file
        def write_out(self,filepath):
                if os.path.isfile(filepath):
                        os.remove(filepath)

                with open(filepath,'a') as f:
                        f.write(self.name + '\n')
                        f.write(str(self.lattice_scalar) + '\n')
                        for row in self.factored_lattice_parameters:
                                f.write('  ' + ' '.join([f'{i:11.6f}' for i in row]) + '\n')

                        for row in self.atom_header:
                                f.write(' ' + '  '.join([f'{i:>4}' for i in row]) + '\n')

                        for frame_num,frame in enumerate(self.pos_data):
                                f.write(f'Direct configuration={str(frame_num+1):>6}\n')
                                for rownum,row in enumerate(frame):
                                        f.write('  ' + ' '.join([f'{i:>11.8f}' for i in row]))
                                        f.write('\n')

