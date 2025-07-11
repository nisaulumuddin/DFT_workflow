from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy

class POSCAR:
        def __init__(self,filepath,*,is_str_not_file=False):
                
                if is_str_not_file:
                        self.file_data = [i.split() for i in filepath.splitlines()]
                else:
                        with open(filepath,'r') as f:
                                self.file_data = [i.split() for i in f.read().splitlines()]
                
                # Check to see if the file has the 0s left over from CONTCARs - this section removes the break characters and zeros at the end of the file until real data is seen again
                for revline in reversed(self.file_data):
                        if len(revline) == 0 or revline[0] == '0.00000000E+00' or revline[0] == 'NaN':
                                self.file_data.pop(-1)
                        else:
                                break

                self.name = ' '.join(self.file_data[0])

                # Get the lattice parameter data
                self.lattice_scalar = float(self.file_data[1][0])
                self.factored_lattice_parameters = np.array(self.file_data[2:5]).astype(np.float64)
                self.lattice_parameters = self.lattice_scalar * self.factored_lattice_parameters
                self.volume = np.linalg.det(self.lattice_parameters)
                


                self.a = np.linalg.norm(self.lattice_parameters[0])
                self.b = np.linalg.norm(self.lattice_parameters[1])
                self.c = np.linalg.norm(self.lattice_parameters[2])


                self.alpha = np.arccos(np.dot(self.lattice_parameters[1],self.lattice_parameters[2])/(self.b*self.c))
                self.beta  = np.arccos(np.dot(self.lattice_parameters[0],self.lattice_parameters[2])/(self.a*self.c))
                self.gamma = np.arccos(np.dot(self.lattice_parameters[0],self.lattice_parameters[1])/(self.a*self.b))


                self.alpha_deg = self.alpha*180/np.pi
                self.beta_deg = self.beta*180/np.pi
                self.gamma_deg = self.gamma*180/np.pi

                # Get atom header data
                self.atom_header = [self.file_data[5],[int(i) for i in self.file_data[6]]]
                self.atom_list = []
                for speciesnum,species in enumerate(self.atom_header[0]):
                        for occur in range(self.atom_header[1][speciesnum]):
                                self.atom_list.append(species)

                self.atom_list = np.array(self.atom_list)
                self.atom_count = len(self.atom_list)

                # Check whether the POSCAR has selective dynamics on
                startstr = ' '.join(self.file_data[7])[0]
                if startstr == 's' or startstr == 'S': # If the first letter is s or S its selective dynamics
                        self.selective_dynamics = True
                else:
                        self.selective_dynamics = False


                # If there is sel. dyn. then the atom position type (cart or direct) will be on the
                if self.selective_dynamics:
                        startstr = ' '.join(self.file_data[8])[0]
                else:
                        startstr = ' '.join(self.file_data[7])[0]
                if startstr == 'd' or startstr == 'D':
                        self.coordinate_type = 'dir'
                else:
                        self.coordinate_type = 'cart'
                
                # Get atomic positions and fixation data if it exists
                if self.selective_dynamics:
                        self.pos_data = np.array(self.file_data[9:])[:,:3].astype(np.float64)
                        self.fixation_data = np.array(self.file_data[9:])[:,3:]
                else:
                        self.pos_data = np.array(self.file_data[8:]).astype(np.float64)
                        self.fixation_data = None
        
                
                if self.selective_dynamics:
                        count_t = 0
                        count_f = 0
                        for f in self.fixation_data:
                                fixation = f[1]
                                if fixation == 'T' or fixation == 't': # if the atom is relaxed, i.e. T for True
                                        count_t += 1
                                else:
                                        count_f += 1
                        self.count_true  = count_t
                        self.count_false = count_f
                else:
                        elem_count = [int(x) for x in self.file_data[6]]
                        self.count_true = np.sum(elem_count)
                        self.count_false = 0

        # Reorder the atom coordinates in ascending/descending way per element on the POSCAR
        def reorder_coord(self, order):
                orig_posn = deepcopy(self.pos_data)
                atom_list = self.atom_list
                elm_count = deepcopy(self.atom_header[1])
                fix_data = deepcopy(self.fixation_data)

                #elm_count.insert(0,0) 
                elm = 0
                new_posn = np.zeros(np.shape(orig_posn))
                new_fix = [] 
                endline = 0
                
                while elm < len(elm_count):
                        line_a = elm_count[elm]
                        
                        elm_posn = np.array(orig_posn[endline:line_a+endline])
                        
                        if order.startswith('d') or order.startswith('D') : #descending
                                elm_posn = elm_posn[np.lexsort(( elm_posn[:,0],  elm_posn[:,1], elm_posn[:,2]))][::-1]
                                new_posn[endline:line_a+endline] = elm_posn

                        elif order.startswith('a') or order.startswith('A') : #ascending 
                                elm_posn = elm_posn[np.lexsort(( elm_posn[:,0]  , elm_posn[:,1], elm_posn[:,2]))]
                                new_posn[endline:line_a+endline] = elm_posn
                        else:
                                new_posn[endline:line_a+endline] =  np.array(elm_posn)
                        endline += line_a
                        elm += 1
                #new_posn = np.concatenate(np.array(new_posn).ravel())
                reordered_indices = []
                for ind,line in enumerate(orig_posn) :
                        indices = [i for i in range(0,len(new_posn)) if np.array_equal(new_posn[i], line)]
                        reordered_indices.append(indices)                

                try:
                        new_fix = fix_data[reordered_indices].reshape(np.shape(fix_data))
                except:
                        new_fix = []

                self.pos_data = new_posn
                self.fixation_data = new_fix
                return self.pos_data , self.fixation_data

        # Fix and relax user-specified atomic layers within a structure (number = number of layers, direction: upper/lower/both,  setup : relax/fix/custom) 
        def setup_layers(self,  number, direction, setup):
                 
                orig_posn = deepcopy(self.pos_data)
                atom_list = self.atom_list
                elm_count = deepcopy(self.atom_header[1])
                
                reordered_posn, reordered_fix  = deepcopy(self.reorder_coord('d'))

                if len(reordered_fix) == 0 : # if selective dynamics was turned off, then by default we freeze all the other atoms
                        atoms_quant = np.sum(self.atom_header[1])
                        fix = ['F','F','F']
                        reordered_fix = []
                        for line in range(0,atoms_quant):
                                reordered_fix.append(fix)
                        reordered_fix = np.array(reordered_fix)
                        reordered_fix = reordered_fix.reshape(np.shape(reordered_posn))
                        
                ref_order = orig_posn[orig_posn[:,2].argsort()][::-1]
                unique_z = np.unique(ref_order[:,2])

                if direction.startswith('u') or direction.startswith('U'): # if upper layers need to be fixed:
                        unique_z = unique_z[unique_z.argsort()]
                        ref_z = unique_z[0:number]
                elif direction.startswith('l') or direction.startswith('L'):  # if lower layers need to be fixed
                        unique_z = unique_z[unique_z.argsort()][::-1]  
                        ref_z = unique_z[0:number]
                elif direction.startswith('b') or direction.startswith('B'):  # if both lower and upper layers need to be fixed
                        unique_z = unique_z[unique_z.argsort()]
                        ref_z.append(unique_z[0:number])
                        unique_z = unique_z[unique_z.argsort()][::-1]
                        ref_z.append( unique_z[0:number])
                else:
                        print("incomplete arguments")

                if setup.startswith('f') or setup.startswith('F'):
                        fix = ['F','F','F']
                elif setup.startswith('r') or setup.startswith('R'):
                        fix = ['T','T','T']
                elif setup.startswith('c') or setup.startswith('C'):
                        inp = input('please specify F/T in the x,y,z direction: ') 
                        fix = np.array(inp)
                else: 
                        print("incomplete arguments")   


                for ind, posn in enumerate(reordered_posn):
                        if posn[2] in ref_z:
                                reordered_fix[ind] = fix


                

                self.fixation_data = reordered_fix
                self.pos_data = reordered_posn 
                self.selective_dynamics = True

                return self.pos_data, self.fixation_data                 

        # Return the cartesian coordinates for the atoms
        def get_cart_positions(self):
                if not self.coordinate_type == 'cart':
                        return np.matmul((self.lattice_scalar*self.factored_lattice_parameters),self.pos_data.T).T
                else:
                        return self.pos_data
        
        # Return the fractional coordinates for the atoms
        def get_dir_positions(self):
                if not self.coordinate_type == 'dir':
                        return np.matmul(np.linalg.inv(self.lattice_scalar*self.factored_lattice_parameters),self.pos_data.T).T
                else:
                        return self.pos_data


        # Refresh the atom header based on the curret atomlist
        def remake_atom_header(self):
                temp_atom_header = []
                current_species = ''
                current_species_count = 0
                for atom in self.atom_list:
                        if current_species != atom:
                                if current_species != '':
                                        temp_atom_header.append([current_species,current_species_count])
                                current_species = atom
                                current_species_count = 1
                        else:
                                current_species_count += 1

                temp_atom_header.append([current_species,current_species_count])
                self.atom_header = np.array(temp_atom_header).T.tolist()

        # Return the index of the atom with the coordinates closes to coord (a list of x,y,z)
        def find_closest_atom(self,coord):
                sse = np.asarray(self.pos_data).astype(float)
                for axis in range(3):
                        sse[:,axis] = np.square(sse[:,axis] - coord[axis])

                return (sse[:,0] + sse[:,1] + sse[:,2]).argmin()

        # Fix all of the atoms above a specific c cutoff coordinate
        def affix_above(self,cutoff):
                # If selective dynamics was not enabled, copy over the pos_data array to the fixation_data (has the same shape as we need but will be overwitten)
                if not self.selective_dynamics:
                        self.selective_dynamics = True
                        self.fixation_data = deepcopy(self.pos_data.astype(str))
                
                for rownum, row in enumerate(self.fixation_data):
                        for colnum, col in enumerate(row):
                                if self.pos_data[rownum][2] > cutoff:
                                        self.fixation_data[rownum][colnum] = 'T'
                                else:
                                        self.fixation_data[rownum][colnum] = 'F'

        # Creates supercells of the current lattice - this operation can only be performed once
        def replicate(self,xrep,yrep,zrep):
                
                # Because the math is incredibly hard to do with cartesian coordinates, we first convert over to fractional
                self.pos_data = self.get_dir_positions()

                # Make sure that we have a copy of the coordinates that we can work off of
                original_positions = deepcopy(self.pos_data)
                original_atom_list = deepcopy(self.atom_list)
                original_fixation_data = deepcopy(self.fixation_data)

                # Now we need to loop through all of the replications that we need to make
                # Keep in mind that a rep=1 means that no replication in that direction happens (think of them as multipliers almost)
                for xnum in range(xrep):
                        for ynum in range(yrep):
                                for znum in range(zrep):
                                        
                                        # Make sure that we're not going to replicate the positions that alreadu exist
                                        if not(xnum + ynum + znum == 0):
                                                replicated_positions = deepcopy(original_positions)
                                                
                                                # Add the replicated positions
                                                replicated_positions[:,0] += xnum
                                                replicated_positions[:,1] += ynum
                                                replicated_positions[:,2] += znum
                                                
                                                # Append the atoms
                                                self.pos_data = np.append(self.pos_data,replicated_positions,axis=0)
                                                self.atom_list = np.append(self.atom_list,original_atom_list)

                                                # If we have selective dynamics, make sure to add the flags alongside the position data
                                                if self.selective_dynamics:
                                                        self.fixation_data = np.append(self.fixation_data,original_fixation_data,axis=0)
                
                # Increase the lengths of the lattice vectors accordingly
                self.factored_lattice_parameters[0] = self.factored_lattice_parameters[0] * xrep
                self.factored_lattice_parameters[1] = self.factored_lattice_parameters[1] * yrep
                self.factored_lattice_parameters[2] = self.factored_lattice_parameters[2] * zrep

                # Make sure that the un-factored lattice parameters are updated as well
                self.lattice_parameters = self.lattice_scalar * self.factored_lattice_parameters

                # Recalculate the direct coordinates to be between 0 and 1 (we already adjusted the lattice vectors)
                self.pos_data[:,0] = self.pos_data[:,0] / xrep
                self.pos_data[:,1] = self.pos_data[:,1] / yrep
                self.pos_data[:,2] = self.pos_data[:,2] / zrep

                # If we originally had cartesian coordinates, revert them back to them
                if self.coordinate_type == 'cart':
                        self.pos_data = self.get_cart_positions()

                # Remake the atom header using all of the new atoms
                self.remake_atom_header()

        # Create the ICORE POSCAR for the atom in index atomnum
        def create_XPS_POSCAR(self,atomnum):
                if self.selective_dynamics:
                        atom_selective = self.fixation_data[atomnum]
                        self.fixation_data = np.delete(self.fixation_data,atomnum,0)
                        self.fixation_data = np.append(self.fixation_data,[atom_selective],axis=0)
                        
                        
                atomtype = self.atom_list[atomnum]
                atompos = self.pos_data[atomnum]
                
                self.pos_data = np.delete(self.pos_data,atomnum,0)
                self.atom_list = np.delete(self.atom_list,atomnum,0).tolist()
                
                self.remake_atom_header()
                
                self.atom_list = np.append(self.atom_list,atomtype)
                self.pos_data = np.append(self.pos_data,[atompos],axis=0)
                
                self.atom_header = np.append(np.array(self.atom_header),[[atomtype],[1]],axis=1).tolist()

        # Return a list of indicies of the atoms of 'atomtype' in the atom_list
        def get_atom_type_match_indicies(self,atomtype):
                return np.argwhere(self.atom_list == atomtype).T[0].tolist()

        def get_lattice_parameters(self,filepath):
                a = np.linalg.norm(self.lattice_parameters[0])
                b = np.linalg.norm(self.lattice_parameters[1])
                c = np.linalg.norm(self.lattice_parameters[2])


                alpha = np.arccos(np.dot(self.lattice_parameters[1],self.lattice_parameters[2])/(b*c))
                beta  = np.arccos(np.dot(self.lattice_parameters[0],self.lattice_parameters[2])/(a*c))
                gamma = np.arccos(np.dot(self.lattice_parameters[0],self.lattice_parameters[1])/(a*b))


                alpha_deg = alpha*180/np.pi
                beta_deg = beta*180/np.pi
                gamma_deg = gamma*180/np.pi

                vol = self.volume
                area_prod = np.cross(self.lattice_parameters[0] , self.lattice_parameters[1])
                surf_area = np.linalg.norm(area_prod)

                with open(filepath,'w') as f:
                        f.write(f' Lattice parameters \n a = {a} , b = {b} , c = {c} \n  (in radians)  alpha = {alpha} , beta = {beta} , gamma = {gamma} \n (in degrees)   alpha = {alpha_deg},  beta = {beta_deg} ,   gamma = {gamma_deg} \n volume = {vol}  \n  surface_area (axb) = {surf_area} '  )
                        
                with open(filepath,'r') as f:
                        print(f.read())

        # Write out information with the atom_order style atomic numbering for lammps
        def write_out_lammps(self,filepath,atom_order):
                if os.path.isfile(filepath):
                        os.remove(filepath)
                
                a = np.linalg.norm(self.lattice_parameters[0])
                b = np.linalg.norm(self.lattice_parameters[1])
                c = np.linalg.norm(self.lattice_parameters[2])

                # Lattice angles in radians
                alpha = np.arccos(np.dot(self.lattice_parameters[1],self.lattice_parameters[2])/(b*c))
                beta  = np.arccos(np.dot(self.lattice_parameters[0],self.lattice_parameters[2])/(a*c))
                gamma = np.arccos(np.dot(self.lattice_parameters[0],self.lattice_parameters[1])/(a*b))

                cart_pos = self.get_cart_positions()

                lx = a
                xy = b*np.cos(gamma)
                xz = c*np.cos(beta)
                ly = np.sqrt(b**2-xy**2)
                yz = (b*c*np.cos(alpha)-xy*xz)/ly
                lz = np.sqrt(c**2-xz**2-yz**2)

                with open(filepath,'a') as f:
                        # Write out the lattice information
                        f.write(f'\n{len(self.atom_list)} atoms\n{len(atom_order)} atom types\n   0.00000   {lx:>10.5f} xlo xhi\n   0.00000   {ly:>10.5f} ylo yhi\n   0.00000   {lz:>10.5f} zlo zhi\n{xy:>10.5f}   {xz:>10.5f}   {yz:>10.5f} xy xz yz\n\nMasses\n\n')
                        
                        # Get the dictionary of mass data
                        mass_dict = getMassDict()
                        
                        for atomnum,atom in enumerate(atom_order):
                                # Calcualte the mass of the atom in amu
                                atom_mass = mass_dict[atom]*6.0221409E23*1000

                                # Write the atom mass
                                f.write(f'{atomnum+1:>6} {atom_mass:>10.5f}\n')
                        
                        f.write(f'\nAtoms\n\n')

                        for atomnum,atom in enumerate(self.atom_list):
                                atom_car_pos = cart_pos[atomnum]
                                f.write(f'  {atomnum+1:>4}  {atom_order.index(atom)+1:>4} {atom_car_pos[0]:>15.10f} {atom_car_pos[1]:>15.10f} {atom_car_pos[2]:>15.10f}\n')

        def write_out_potcar(self,filepath,potcar_root):
                
                # Make sure that potcar_root ends with / (that way it's a directory)
                if potcar_root[-1] != '/':
                        potcar_root += '/'

                potcar = ''
                for atom in self.atom_header[0]:
                        with open(potcar_root + atom + '/POTCAR','r') as f:
                                potcar += f.read()

                with open(filepath,'a') as f:
                        f.write(potcar)

        # Write out the POSCAR file
        def write_out(self,filepath):
                if os.path.isfile(filepath):
                        os.remove(filepath)
                
                with open(filepath,'a') as f:
                        f.write(self.name + '\n')
                        f.write(str(self.lattice_scalar) + '\n')
                        for row in self.factored_lattice_parameters:
                                f.write('   '.join([f'{i:15.10f}' for i in row]) + '\n')
                        
                        for row in self.atom_header:
                                f.write('  '.join([f'{i:>4}' for i in row]) + '\n')
                        
                        if self.selective_dynamics:
                                f.write('Selective Dynamics\n')
                        
                        if self.coordinate_type == 'dir':
                                f.write('Direct\n')
                        else:
                                f.write('Cartesian\n')
                        if np.shape(self.pos_data) == np.shape([1,2,3]):
                                f.write('   '.join(map(str,self.pos_data)))
                        else:
                                for rownum,row in enumerate(self.pos_data):
                                        f.write('  '.join([f'{i:20.13f}' for i in row]))
                                        if self.selective_dynamics:
                                                f.write('   ' + ' '.join(self.fixation_data[rownum]) + '\n')
                                        else:
                                                f.write('\n')