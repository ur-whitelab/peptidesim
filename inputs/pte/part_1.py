from peptidesim import PeptideSim
import textwrap, sys, re, os
import sys
import gromacs
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
MPI_NP = 4
peptide_copies=int(sys.argv[3])

#try to reload                                                                                
if(os.path.exists(pickle_name)):
    print 'loading restart'
    with open(pickle_name, 'r') as f:
        ps = pickle.load(f)
else:
    ps = PeptideSim(name, [seq], [peptide_copies], job_name='{}'.format(name))
    ps.mdrun_driver='gmx_mpi'
    ps.forcefield='amber99sb'
    ps.water='tip4p'
    ps.peptide_density = 0.008#mg/ml                                                              
    ps.ion_concentration=0.001#10mM                                                 

ps.initialize()
file00=ps.pdb_file
cwd=os.getcwd()
gromacs.pdb2gmx(f=file00,o="file_pdb_output.pdb",water=ps.water,ff=ps.forcefield)
file00="{}/{}".format(cwd,"file_pdb_output.pdb")
print file00,"file_pdb_output.pdb", cwd
def total_aa(file1,output_file):
    print ps.sims[-1].location
    output=open(output_file,'w')
    with open(file1, "r") as f:
        lines=f.readlines()
        i=4
     
        for line in lines[:4]:
            output.write(line)
            
     
        for line in lines[4:]:
            old_line=line
            line=line.split()
            if(len(line)>1 and line[3]=='SOL'):
                
                output.write('TER\n')
                output.write('ENDMDL\n')      
                output.close()
                break
            i=i+1
            output.write(old_line)
        last_line=[]
        print lines[i-1],lines[i-2], lines[i-3]
        if (lines[i-1].startswith('END')==True):
                
            last_line=lines[i-3].strip()
            print "ends with END"
        else:
            last_line=lines[i-1].strip()
            print "doesnt end with END"

        last_line=last_line.split()

        output.close()
        if (last_line[4].isdigit()==True):
            return int(last_line[1]),int(last_line[4]),output_file
        else:
            return int(last_line[1]),int(last_line[5]),output_file
print total_aa(file00,'template.pdb')

def number_all_atoms_chains():
    if (os.path.isdir("{}/".format(os.getcwd)+"data")==True):
        template_file="{}/".format(os.getcwd)+"data/template.pdb"
        
        with open(template_file, 'r') as f:
            lines=f.readlines()
            lines=lines[-3].strip()
            lines=lines.split()
            return lines[1],last_lines[5]
total_no_atoms,number_chains,file0=total_aa(file00,'template.pdb')
atoms_in_chain=total_no_atoms/peptide_copies


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
        
           newp pdb file'''
    from string import ascii_uppercase
    print input_pdbfile,number_of_chains,atoms_in_chain
    with open(input_pdbfile, 'r') as f:
        lines=f.readlines()
        beginning=lines[:first_atom_index]                  #saves the first useless lines that don't contain conf info
        lines=lines[first_atom_index:]                      #gets rid of those lines from readlines
        first_line=lines[0].strip()
        first_line=first_line.split()
        with open(output_pdbfile,'w') as f:
            for index in np.arange(len(beginning)):         #iterates through useless lines
                f.write("{}\n".format(beginning[index].strip())) #writes the useless lines into the new file   
            for i in np.arange(number_of_chains):           #iterates through the copies of chains  
                for j in np.arange(atoms_in_chain):         #itarates through the atoms in the chain
                    print i*(atoms_in_chain+1)+j,i*(atoms_in_chain)+j, len(lines), i,j
                    a=lines[i*(atoms_in_chain)+j]             #reads the old pdbfile info pertaining to the atoms of interest  
                    a=list(a)                               #converts the string into a list of characters 
                    if (len(a)>=21):
                        a[21]=ascii_uppercase[i]                #replaces the repetittive chain ID with unique chain ID chosen from 
                        
                    a="".join(a)
                            #puts the list of characters back to a string that contains line info
                    f.write('{}'.format(a))                 #writes the line into the new pdbfile

                    #if(int(i*(atoms_in_chain-1)+i+j)==int((atoms_in_chain-1)*i+atoms_in_chain-1+i) and peptide_copies!=1): #checks whether the atoms is at the end of the chain    
                    if(j==int(atoms_in_chain-1) and peptide_copies!=1):   #print int(i*atoms_in_chain+j),int((atoms_in_chain*i+atoms_in_chain-1))
                        f.write('TER\n')                    #puts ter at the end of each chain
                    

            if (peptide_copies!=1):
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
    if (name!='HAshifts.dat' and name!='Hshifts.dat'):                       #checks whether the exp chem shifts exist, for now only H shifts are available

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
            chemical_shifts=np.genfromtxt('/scratch/damirkul/pep_simulations/sequences_amyloid/amyloid_21_30/data/{}'.format(name),comments='!#')
#            print chemical_shifts
            for i in range(copies_chains):
                beg_end_chains.append(i*number_chains+1)                   #adds the begining amino acid index of the chain i, first index starts at 1
                beg_end_chains.append(i*number_chains+number_chains)       #adds the amino acid ending index of the chain
                #print beg_end_chains
                for j in range(number_chains):                             #iterates through the amino acid indeces in the chain
                    j=j+1                                                  #amino acid index also starts at 1
                    if (j in beg_end_chains):                              #checks whether the amino acid index is at the beginning or ending of the chain
                        #print i,j,number_chains*i+j,chemical_shifts[j-1][1],'a' 
                        f.write('#{} {}\n'.format(number_chains*i+j,chemical_shifts[j-1][1]))       #puts a hash sign for beginning and amino acid indeces in the chain
                    else:
                        #print i,j,number_chains*i+j,chemical_shifts[j-1][1]
                        f.write('{} {}\n'.format(number_chains*i+j,chemical_shifts[j-1][1]))       #copies the chem shifts from file with chemical shifts for the correspoding amino acid
                        
filenames=['Cshifts.dat','CAshifts.dat','HAshifts.dat','Hshifts.dat','CBshifts.dat', 'Nshifts.dat']
for i in filenames:
    data_folder(number_chains,i,peptide_copies)
ps.add_file(directory)
with open(ps.pickle_name, 'w') as f:
    pickle.dump(ps, file=f)

print ps.pickle_name, 'picklename1'
with open(ps.pickle_name, 'w') as f:
    pickle.dump(ps, file=f)

ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 8*10**2,'rcoulomb':1}, mpi_np=MPI_NP)

print ps.pickle_name, 'picklename2'
with open(ps.pickle_name, 'w') as f:
    pickle.dump(ps, file=f)

ps.run(mdpfile='peptidesim_anneal.mdp',tag='annealing',mdp_kwargs={'nsteps':int(200* 5*10**2)},mpi_np=MPI_NP)#change the time step to 2 ns

print ps.pickle_name, 'picklename3'
with open(ps.pickle_name, 'w') as f:
    pickle.dump(ps, file=f)

ps.run(mdpfile='peptidesim_npt.mdp', tag='equil_npt', mdp_kwargs={'nsteps': int(200 * 5*10**2),'ref_t':278}, mpi_np=MPI_NP)

print ps.pickle_name, 'picklename3'
with open(ps.pickle_name, 'w') as f:
    pickle.dump(ps, file=f)
            

