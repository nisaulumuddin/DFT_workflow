#!/usr/bin/env python3

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy

class OSZICAR:
        def __init__(self,filepath):
                self.mag = []

                with open(filepath,'r') as f:
                        self.file_data = f.read().splitlines()

                        for num,line in enumerate(self.file_data):
                                if 'mag=' in line:
                                        self.mag.append(float(line.split()[9]))

class OUTCAR:
        def __init__(self,filepath):
                # Set some constants that will be checked against in the file read
                self.is_converged = False
                self.volume = None
                self.ionic_steps = 0
                self.selective_dynamics = False
                self.elapsed_time = None
                self.e_sigzero = []
                self.forces = []
                self.mag_forces = []
                self.is_neb = False
                self.tangent_vectors = []
                self.scf_cycle = 0
                self.is_static = False
                self.magmom = []
                self.indmagmom = []
                self.finmagmom = 0
                self.nions = 0
                self.atoms = []

                # By only reading the OUTCARs once, the operation time is cut significantly (althought the string searching takes a while)
                with open(filepath,'r') as f:
                        self.file_data = f.read().splitlines()
                        for num,line in enumerate(self.file_data):
                                # This indicates a converged file
                                if 'reached required accuracy - stopping structural energy minimisation' in line:
                                        self.is_converged = True
                                
                                elif 'NIONS =' in line:
                                        self.nions = line.split()[11]

                                elif 'NELM' in line:
                                        nelm = line.split()[2]
                                        self.nelm = nelm.replace(';','')
                                        self.nelm = int(self.nelm)

                                elif 'number of electron ' in line:
                                        try:
                                                totmagmom = line.split()[5]
                                                self.magmom.append(float(totmagmom))
                                                self.finmagmom = totmagmom
                                        except:
                                                self.magmom = 0
                                                self.finmagmom = 0 

 
                                elif 'IBRION =' in line:
                                        self.ibrion = int(line.split()[2])
                                        
                                # Grab the volume of the cell
                                elif 'volume of cell' in line:
                                        self.volume = float(line.split()[4])
                                
                                # This only appears at the end of an ionic step, also shows the energy after the step
                                elif 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM' in line:
                                        self.ionic_steps += 1
                                        self.e_sigzero.append(float(self.file_data[num+4].split('=')[-1].split()[0]))
                                
                                # This counts the total number of SCF cycles during the simulation
                                elif 'Iteration' in line:
                                        self.scf_cycle += 1
                                
                                elif 'magnetization (x)' in line:
                                        self.indmagmom = []
                                        start_line = num+3
                                        for atom in range(1,int(self.nions)+1):
                                                atom_mag =  self.file_data[start_line+atom]
                                                self.indmagmom.append(atom_mag.split()[4])

                                # List of atoms in the same order as the POSCAR
                                elif 'TITEL' in line:
                                        self.atoms.append(line.split()[3])

                                       
                                # The time elapsed
                                elif 'Elapsed time (sec)' in line:
                                        self.elapsed_time = float(line.split()[3])
                                
                                # The forces listed at the end of each ionic step
                                elif 'TOTAL-FORCE' in line:
                                        forces = []
                                        for force_line in self.file_data[num+2:]:
                                                # There is a nice line drawn at the end of the forces, stop looping once you see it
                                                if '------------' in force_line:
                                                        break
                                                else:
                                                        forces.append([float(i) for i in force_line.split()[3:]])
                                        
                                        # We want the forces of all of the atoms
                                        self.forces.append(forces)

                                        # Get the magnitude of all of the force vectors
                                        self.mag_forces.append([np.linalg.norm(i) for i in self.forces[-1]])
                                
                                # NPAR value VASP 5
                                elif 'CORES_PER_BAND=' in line:
                                        self.npar = line.split()[7]
                                
                                # NPAR value VASP 6
                                elif 'one band on NCORE=' in line:
                                        self.npar = line.split()[5]

                                # The total number of cores
                                elif 'running on' in line: 
                                        self.cores = line.split()[2]
                                        if self.cores == "running":
                                                self.cores = line.split()[4]

                                
                                
                                # The maximum memory used (kb)
                                elif 'Maximum memory used (kb):' in line:
                                        self.totmem = line.split()[4]

                                # Check if the OUTCAR is for an NEB
                                elif 'CHAIN: Running the NEB' in line:
                                        self.is_neb = True
                                
                                # Grab NEB tangent band direction unit vectors
                                elif 'NEB: Tangent' in line:
                                        tangent_vectors = []
                                        for tangent_vector_line in self.file_data[num+2:]:
                                                # There is a blank line at the end of the tangent section, stop looking once you see it
                                                if len(tangent_vector_line.split()) == 0:
                                                        break
                                                else:
                                                        # Turn the vector elements into a float
                                                        vector = np.array([float(i) for i in tangent_vector_line.split()])

                                                        # Turn the vector into an actual unit vector by dividing the original magnitude
                                                        if not np.linalg.norm(vector) == 0:
                                                                vector = vector / np.linalg.norm(vector)
                                                        else:
                                                                vector = np.array([0,0,0])
                                                        
                                                        # The way that numpy handles appending confuses me - I like doing it all in python lists
                                                        tangent_vectors.append(vector.tolist())

                                        # We want the vectors of all of the atoms
                                        self.tangent_vectors.append(tangent_vectors)
                               
                                # Prints how many atoms type there are in the POSCAR  

                                elif 'ions per type' in line:
                                        ionspertyp = []
                                        iontyps = len(line.split()[4:])
                                        self.iontyps = float(iontyps)

                                        for typ in line.split()[4:]:
                                                ionspertyp.append(typ)
                                        self.ionspertyp = np.array(ionspertyp).astype(int)
                                
                                # Prints the the line number of where the total charge of the first atom is displayed
                                elif 'LOOP+' in line:
                                        self.occupationline = int(num) + 9
                                        
                                        div = self.file_data[int(num) + 7 ].split()
                                        # also prints the list of l orbitals that is modeled as the valence 
                                        for l in div[3:-1]:
                                                self.l_orb = l
                                
                                elif 'LDAUL' in line:
                                        self.ldaul = np.array(line.split()[7:]).astype(int)


                # If forces were found, turn them into arrays
                if len(self.forces) > 0:
                        self.forces = np.array(self.forces)
                        self.mag_forces = np.array(self.mag_forces)

                # If NEB tangent vectors were found, turn them into arrays
                if len(self.tangent_vectors) > 0:
                        self.tangent_vectors = np.array(self.tangent_vectors)
        
                # This returns if the static calculation has converged to electronic minimum
                if (self.ibrion == int(-1) and self.nelm >= self.scf_cycle and self.elapsed_time):
                        self.is_static = True

                # This returns the list of atoms in the same order as the POSCAR
                i=0
                atomlists = []
                for ind in self.ionspertyp:
                        concatlist = ' '.join(ind*[str(self.atoms[i])])
                        atomlists.append(concatlist.split(' '))
                        i+=1
                self.atomlist = [item for row in atomlists for item in row]
                        

                        
        def ion_occupancies(self, ldau, alpha):
                # Storing the converged ion occupancy as an array
                        # alpha: the perturbation value inserted as LDAUU when LDATYPE = 3. set as 0 if none.
                        # ldau: the 'ldaul' tag in INCAR lets us know which atom occupancies should be taken into account if they are not set to -1. use angular quantum number: (s = 0 , p = 1, d = 2 , f = 3) 
                totalatoms = np.sum(np.array(self.ionspertyp).astype(int))
               
                ion_occ = [] 
                for line in self.file_data[self.occupationline: self.occupationline+totalatoms ]:
                        ion_occ.append(line.split())
               
                ion_occ = np.array(ion_occ).astype('float')
                ang_occ = [float(alpha)]
                # Stores the ion occupancy for a specific angular momentum (angmom) 
                
                ldaumatrix = []
                for i in range(len(self.ionspertyp)):
                        ldaumatrix += self.ionspertyp[i] * [ldau[i]]
                 

                count = 0
                for ion in ion_occ:
                        if int(ldaumatrix[count]) >= 0:
                                col = int(ldaumatrix[count]) + 1 
                                ang_occ.append(ion[col])
                        count += 1
                
                ang_occ_f = np.array(ang_occ).astype(float)
                ion_occ_f = np.array(ion_occ).astype(float)

                return ion_occ_f, ang_occ_f

        def select_occupancies(self, ldaul, alpha, sel_atoms):
                # Storing select ion occupancies as an array based on sel_atoms array 
                        # sel_atoms: the selected atoms 
                        # alpha : the peturbation value inserted as LDAUU when LDATYPE = 3. set as 0 if none.
                        # ldaul : user-specified angular quantum numbers (s = 0 , p = 1 , d = 2, f = 3) 

                totalatoms = np.sum(np.array(self.ionspertyp).astype(int))
                ion_occ = []

                for line in self.file_data[self.occupationline: self.occupationline+totalatoms ]:
                        ion_occ.append(line.split())

                ion_occ = np.array(ion_occ).astype('float')              
                sel_occ = [np.float(alpha)]

                count = 0
                sel_atoms = np.asarray(sel_atoms).flatten()
                ldaul = np.asarray(ldaul).flatten()


                for at in sel_atoms:
                        sel_occ.append(ion_occ[int(at)-1][ldaul[count]+1])
                        count+=1 
                        
                sel_occ_f = np.array(sel_occ).astype(np.float)
                return sel_occ_f




        def LinResponseMatrix(self, *ang_occ_f):
                size = np.shape(ang_occ_f)

                return RespMatrix

               

class INCAR:
        def __init__(self,filepath):

                with open(filepath,'r') as f:
                        self.file_data = [i.split() for i in f.read().splitlines()]
                
                self.tags = {}
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

                                #if len(" ".join(fields).split("!")[0]) != 0:
                                #        line = " ".join(fields).split("!")[0]
                                #        tagname = line.split("=")[0].split()[0]
                                #        tagdat = line.split("=")[1].split()
                                #        self.tags[tagname] = tagdat

                

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

        def write_out(self,filepath):
                if os.path.isfile(filepath):
                        os.remove(filepath)

                with open(filepath,'a') as f:
                        f.write(self.get_sorted_tags_string())
        
                


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


                
    
