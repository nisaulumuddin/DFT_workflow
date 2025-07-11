from __future__ import annotations

import numpy as np
#import scipy.interpolate as si
import os
from copy import deepcopy

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
