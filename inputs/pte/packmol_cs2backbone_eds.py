from peptidesim import PeptideSim
import textwrap, sys, re, os
import sys
import dill as pickle
from shutil import copyfile, move

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
seq = sys.argv[1]
name = sys.argv[2]
debug = False
pickle_name = name + '.pickle'
MPI_NP = 8
peptide_copies=4

#try to reload                                                                                
if(os.path.exists(pickle_name)):
    print 'loading restart'
    with open(pickle_name, 'r') as f:
        ps = pickle.load(f)
else:
    ps = PeptideSim(name, [seq], [peptide_copies], job_name='2mer_{}'.format(name))
    ps.mdrun_driver='gmx_mpi'
ps.forcefield='amber99sb'
#ps.water='tip4pew'                                                                           
ps.water='tip4p'
ps.peptide_density = 0.009#mg/ml                                                              
ps.ion_concentration=0.001#10mM                                                               
ps.initialize()

ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**2,'rcoulomb':1}, mpi_np=MPI_NP)
file0=ps.pdb_file
def total_aa(file1):
    with open(file1, "r") as f:
        lines=f.readlines()
        i=0
        for line in lines:
            if(line.startswith('TER')==True):
                break
            i=i+1
        last_line=lines[i-2].strip()
        last_line=last_line.split()
        print last_line
        return int(last_line[1])+1,int(last_line[5])
total_no_atoms,number_chains=total_aa(file0)
atoms_in_chain=total_no_atoms/peptide_copies
print atoms_in_chain
with open(file0,'r') as file_read:
    file_data=file_read.read()
file_data=file_data.replace('HIS','HIE')
with open(file0,'w') as file_write:
    file_write.write(file_data)

def pdbfile_generator_w_chain_id(number_of_chains,atoms_in_chain,first_atom_index,output_pdbfile,input_pdbfile):
    '''a funtion that takes an old pdbfile that has all hydrogens but without unique chain IDs 
        and without terminii of chains indicated and generates a new pdbfile with uniqe chain id's
        and terminii indicated. 
       
       var:
       
           number_of_chains-an int that holds copies of the same chain
           atoms_in_chain-an int that holds number of atoms in each chain
           first_atom_index-an index of the line containing first atom 
                            in the old pdbfile (line indexing starts at 0)
           output_pdbfile-a string containing the name of the output pdbfile
           input_pdbfile-a string containing the input pdbfile
        
        returns:
        
            newp pdb file
        '''
    from string import ascii_uppercase
    with open(input_pdbfile,'r')as f:                       #opens the old pdbfile
        lines=f.readlines()                                 #reads all the lines 
        beginning=lines[:first_atom_index]                  #saves the first useless lines that don't contain conf info
        lines=lines[first_atom_index:]                      #gets rid of those lines from readlines
        with open(output_pdbfile,'w') as f:
            for index in np.arange(len(beginning)):         #iterates through useless lines
                f.write("{}\n".format(beginning[index].strip())) #writes the useless lines into the new file   
            for i in np.arange(number_of_chains):           #iterates through the copies of chains  
                for j in np.arange(atoms_in_chain):         #itarates through the atoms in the chain
                    a=lines[i*atoms_in_chain+j]             #reads the old pdbfile info pertaining to the atoms of interest  
                    a=list(a)                               #converts the string into a list of characters 
                    a[21]=ascii_uppercase[i]                #replaces the repetittive chain ID with unique chain ID chosen from 
                    a="".join(a)                            #puts the list of characters back to a string that contains line info
                    f.write('{}'.format(a))                 #writes the line into the new pdbfile
                    if((i*atoms_in_chain+j+1)==(atoms_in_chain*i+atoms_in_chain)): #checks whether the atoms is at the end of the chain
                
                        f.write('TER\n')                    #puts ter at the end of each chain
            f.write('ENDMDL\n')                             #puts finish touches
            f.close()                                       #done
    return output_pdbfile 
file3='template.pdb'
new_pdb_file=pdbfile_generator_w_chain_id(peptide_copies,atoms_in_chain,2,'new_pdb.pdb',file0)

directory='data'
if not os.path.exists(directory):
    os.mkdir(directory)
template_pdb=os.path.join(directory,file3)

camshift_src='/home/damirkul/simulation_software/simulations/AAA/data/camshift.db'
gromacs_a03_mdb_src='/home/damirkul/simulation_software/simulations/AAA/data/a03_gromacs.mdb'
h_src='/scratch/damirkul/pep_simulations/sequences_amyloid/amyloid_21_30/data/Hshifts.dat'
ha_src='/scratch/damirkul/pep_simulations/sequences_amyloid/amyloid_21_30/data/HAshifts.dat'
#c_src='/scratch/damirkul/pep_simulations/sequences_amyloid/amyloid/data_old/Cshifts.dat'     
#ca_src='/scratch/damirkul/pep_simulations/sequences_amyloid/amyloid/data_old/CAshifts.dat'   
#cb_src='/scratch/damirkul/pep_simulations/sequences_amyloid/amyloid/data_old/CBshifts.dat'   
#n_src='/scratch/damirkul/pep_simulations/sequences_amyloid/amyloid/data_old/Nshifts.dat'     
gromacs_a03_mdb_dst=os.path.join(directory,'a03_gromacs.mdb')
camshift_dst=os.path.join(directory,'camshift.db')
template_pdb=os.path.join(directory,file3)
h_dst=os.path.join(directory,'Hshifts.dat')
ha_dst=os.path.join(directory,'HAshifts.dat')
#c_dst=os.path.join(directory,'Cshifts.dat')                                                  
#ca_dst=os.path.join(directory,'CAshifts.dat')                                                
#cb_dst=os.path.join(directory,'CBshifts.dat')                                                
#n_dst=os.path.join(directory,'Nshifts.dat')                                                  
copyfile(h_src,h_dst)
copyfile(ha_src,ha_dst)
#copyfile(c_src,c_dst)                                                                        
#copyfile(ca_src,ca_dst)                                                                      
#copyfile(cb_src,cb_dst)                                                                      
#copyfile(n_src,n_dst)                                                                     
move(new_pdb_file,template_pdb)
copyfile(camshift_src,camshift_dst)
copyfile(gromacs_a03_mdb_src, gromacs_a03_mdb_dst)

center=''


def data_folder(number_amino_acids,name,copies_chains):
    ''' a function that takes total number of amino acids in the pdbfile,name of the chemical shifts file and how many copies of a peptide are present and puts chemical shifta file into the data folder        
    Var:
       number_amino_acids-total number of amino acids in the pdbfile 
       name-name of the chemical shifts file                                                   
       copies_chains-the number of peptide copies
    '''
    new_file=os.path.join(directory,name)                                  #creates a chemical shift file in data folder
    if name is not('HAshifts.dat' or 'Hshifts.dat'):                       #checks whether the exp chem shifts exist, for now only H shifts are available
         with open(new_file,'w') as f:                                     #creates a new file to copy over the chemical shufts
            beg_end_chains=[]                                              #amino acid residue indeces for beginning and end of each chain of peptides                                                                  
            number_chains=number_amino_acids/copies_chains                 #number of amino acids in one chain                                                                                                                    
            for i in range(copies_chains):                                 #iterates throug a number of copies
                beg_end_chains.append(i*number_chains+1)                   #adds the begining amino acid index of the chain i, first index starts at 1
                beg_end_chains.append(i*number_chains+number_chains)       #adds the amino acid ending index of the chain
                for j in range(number_chains):                             #iterates through the amino acid indeces in the chain 
                    j=j+1                                                  #amino acid index also starts at 1
                    if (j in beg_end_chains):                              #checks whether the amino acid index is at the beginning or ending of the chain                        
                        f.write('#{} 1\n'.format(number_chains*i+j))       #puts a hash sign for beginning and amino acid indeces in the chain and 1 for lack of exp chem shift
                    else:
                        f.write('{} 1\n'.format(number_chains*i+j))        #puts 1 for lack of exp chem shifts
    else:
        with open(new_file,'w') as f:
            beg_end_chains=[]                                              #amino acid residue indeces for beginning and end of each chain of peptides                                                                                                         
            number_chains=number_amino_acids/copies_chains                 #number of amino acids in one chain                                                                                                                    
            chemical_shifts=np.genfromtxt('/scratch/damirkul/pep_simulations/sequences_amyloid/amyloid_21_30/data/{}'.format(name))
            for i in range(copies_chains):
                beg_end_chains.append(i*number_chains+1)                   #adds the begining amino acid index of the chain i, first index starts at 1
                beg_end_chains.append(i*number_chains+number_chains)       #adds the amino acid ending index of the chain
                for j in range(number_chains):                             #iterates through the amino acid indeces in the chain
                    j=j+1                                                  #amino acid index also starts at 1
                    if (j in beg_end_chains):                              #checks whether the amino acid index is at the beginning or ending of the chain
                        f.write('#{} {}\n'.format(number_chains*i+j,chemical_shifts[j-1][1]))       #puts a hash sign for beginning and amino acid indeces in the chain
                    else:
                        f.write('{} {}\n'.format(number_chains*i+j,chemical_shifts[j-1][1]))       #copies the chem shifts from file with chemical shifts for the correspoding amino acid

filenames=['Cshifts.dat','CAshifts.dat','HAshifts.dat','Hshifts.dat','CBshifts.dat', 'Nshifts.dat']
for i in filenames:
    data_folder(number_chains,i,peptide_copies)
ps.add_file(directory)

plumed_input=textwrap.dedent(
    '''                                                                                                            
    peptide: GROUP ATOMS=1-{}                                                                                      
    WHOLEMOLECULES ENTITY0=peptide                                                                                 
    cs: CS2BACKBONE ATOMS=peptide DATA=data NRES={} NOPBC                                                          
    eds: EDS CENTER=4 ARG=cs.hn_2 PERIOD=10 TEMP=287                                                               
    phi_27ASN_28LYS: TORSION ATOMS=86,88,90,108#phi 27ASN to 28LYS bend                                            
    psi_28LYS_29GLY: TORSION ATOMS=88,90,108,110#psi 28 LYS and 29GLY                                              
    phi_ALA1_GLU2: TORSION ATOMS=11,13,14,26#phi ALA and GLU                                                       
    psi_ALA1_GLU2: TORSION ATOMS=1,5,11,13#psi ALA and GLU                                                         
    com_glu_22: COM ATOMS=23,24,25#coo                                                                             
    com_lys_28: COM ATOMS=104,105,106,107#nh3                                                                      
    com_asp_23: COM ATOMS=35,36,37#coo                                                                             
    com_val_24: COM ATOMS=44,45,46,47,48,49,50,51,52,53#ch2cch3                                                    
    com_gly_25: COM ATOMS=58,59,60#ch2                                                                             
    distance_val24_K28: DISTANCE ATOMS=com_asp_23,com_val_24#val24 and asp23                                       
    distance_K28_asp23: DISTANCE ATOMS=com_lys_28,com_asp_23#K28 and asp23                                         
    distance_K28_E22: DISTANCE ATOMS=com_lys_28,com_glu_22#K28 and glu22                                           
    distance_SER26HN_asp23: DISTANCE ATOMS=71,com_asp_23#SER26HN and asp23                                         
    distance_GLY25_asp23: DISTANCE ATOMS=com_gly_25,com_asp_23#GLY25 and asp23                                     
    PRINT ARG=phi_27ASN_28LYS,psi_28LYS_29GLY,phi_ALA1_GLU2,psi_ALA1_GLU2,distance_val24_K28,distance_K28_asp23,distance_K28_E22,distance_SER26HN_asp23,distance_GLY25_asp23 FILE={} STRIDE=10                                       
    PRINT ARG=(cs\.exphn_.*),(cs\.expha_.*),(cs\.hn_.*),(cs\.ha_.*) FILE=CS_H_shifts STRIDE=10                     
    PRINT ARG=eds.bias,eds.force2 FILE=eds_output.dat STRIDE=10                                                    
    #PRINT ARG=anti_beta_1,alpha_1 FILE=SECONDARY_STRUCTURE STRIDE=10                                              
    #PRINT ARG=ab FILE=alpha_beta_similarity STRIDE=10                                                             
    ENDPLUMED'''.format(total_no_atoms,number_chains,'OUTPUT_COLVAR_EDS'))
with open('plumed_distance_torsion_eds_copies.dat', 'w') as f:
    f.write(plumed_input)
ps.add_file('plumed_distance_torsion_eds_copies.dat')
ps.run(mdpfile='peptidesim_anneal.mdp',tag='annealing',mdp_kwargs={'nsteps':int(1 * 5*10**2)},mpi_np=MPI_NP, pickle_name=pickle_name )
ps.run(mdpfile='peptidesim_npt.mdp', tag='equil_npt', mdp_kwargs={'nsteps': int(1 * 5*10**2),'ref_p':0.061}, mpi_np=MPI_NP, pickle_name=pickle_name)

final_time=int(1*5*10**2)

ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_tune',  mdp_kwargs={'nsteps':final_time,'ref_t':278 }, mpi_np=MPI_NP,run_kwargs={'plumed':'plumed_distance_torsion_eds_copies.dat'}, pickle_name=pickle_name)
#ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_tune',  mdp_kwargs={'nsteps':final_time,'ref_t':278 }, mpi_np=MPI_NP, pickle_name=pickle_name)#                                                                                       
#short simulation                                                                                                  
#ps.analyze()                                                                                                      
print ps.box_size_angstrom



