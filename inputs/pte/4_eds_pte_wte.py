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
MPI_NP = 16
peptide_copies=1

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
ps.peptide_density = 0.008#mg/ml                                                              
ps.ion_concentration=0.001#10mM                                                 
                                                                                       
ps.initialize()
replicas=16
file00=ps.pdb_file
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
            
        last_line=lines[i-1].strip()
        last_line=last_line.split()
        output.close()
        return int(last_line[1]),int(last_line[4]),output_file

total_no_atoms,number_chains,file0=total_aa(file00,'template.pdb')
print total_no_atoms, "total_no_atoms",file0
atoms_in_chain=total_no_atoms/peptide_copies
print atoms_in_chain, "atoms_in_chain"
eds_period=250
remd_exchange_period=300
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
    print input_pdbfile,number_of_chains,atoms_in_chain
    #with open(input_pdbfile,'r',buffering=-1)as f:                       #opens the old pdbfile
    with open(input_pdbfile, 'r') as f:
        lines=f.readlines()
        #print lines                         #reads all the lines 
        beginning=lines[:first_atom_index]                  #saves the first useless lines that don't contain conf info
        #print beginning
        lines=lines[first_atom_index:]                      #gets rid of those lines from readlines
        with open(output_pdbfile,'w') as f:
            for index in np.arange(len(beginning)):         #iterates through useless lines
                f.write("{}\n".format(beginning[index].strip())) #writes the useless lines into the new file   
            for i in np.arange(number_of_chains):           #iterates through the copies of chains  
                for j in np.arange(atoms_in_chain):         #itarates through the atoms in the chain
                    a=lines[i*atoms_in_chain+j]             #reads the old pdbfile info pertaining to the atoms of interest  
                    a=list(a)                               #converts the string into a list of characters 
         #           print a
                    a[21]=ascii_uppercase[i]                #replaces the repetittive chain ID with unique chain ID chosen from 
#                    print i,j,a
                    a="".join(a)
                    #print j,i*atoms_in_chain+j,atoms_in_chain*i+atoms_in_chain,i
                            #puts the list of characters back to a string that contains line info
                    f.write('{}'.format(a))                 #writes the line into the new pdbfile
                    

                    if(int(i*atoms_in_chain+j)==int((atoms_in_chain*i+atoms_in_chain-1))): #checks whether the atoms is at the end of the chain    
                        print int(i*atoms_in_chain+j),int((atoms_in_chain*i+atoms_in_chain-1))
                        f.write('TER\n')                    #puts ter at the end of each chain
            f.write('ENDMDL\n')                             #puts finish touches
            f.close()                                       #done
    return output_pdbfile 
#file0=ps.pdb_file
file3='template.pdb'
new_pdb_file=pdbfile_generator_w_chain_id(peptide_copies,atoms_in_chain,4,'new_pdb.pdb',file0)

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
        #print name
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


'''
ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**2,'rcoulomb':1}, mpi_np=MPI_NP,pickle_name=pickle_name)
with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)
ps.run(mdpfile='peptidesim_anneal.mdp',tag='annealing',mdp_kwargs={'nsteps':int(2000* 5*10**2)},mpi_np=MPI_NP, pickle_name=pickle_name )#change the time step to 2 ns
with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)
ps.run(mdpfile='peptidesim_npt.mdp', tag='equil_npt', mdp_kwargs={'nsteps': int(2000 * 5*10**2),'ref_t':278}, mpi_np=MPI_NP, pickle_name=pickle_name)
with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)
'''


ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**2,'rcoulomb':1}, mpi_np=MPI_NP,pickle_name=pickle_name)
with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)
ps.run(mdpfile='peptidesim_anneal.mdp',tag='annealing',mdp_kwargs={'nsteps':int(200* 5*10**2)},mpi_np=MPI_NP, pickle_name=pickle_name )#change the time step to 2 ns
with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)
ps.run(mdpfile='peptidesim_npt.mdp', tag='equil_npt', mdp_kwargs={'nsteps': int(200 * 5*10**2),'ref_t':278}, mpi_np=MPI_NP, pickle_name=pickle_name)
with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)

def get_replex_e(ps, replica_number):
    with open(ps.sims[-1].location + '/' + ps.sims[-1].metadata['md-log']) as f:
        p1 = re.compile('Repl  average probabilities:')       
        p2 = re.compile('Repl\s*' + ''.join(['([0-9\.]+\s*)' for _ in range(replica_number - 1)]) + '$')
        ready = False
        for line in f:
            if not ready and p1.findall(line):
                ready = True
            elif ready:
                match = p2.match(line)
                if match:
                    return [float(s) for s in match.groups()]





#pte_result = ps.pte_replica(mpi_np=MPI_NP, max_tries=15,min_iters=10,
                              #mdp_kwargs={'nsteps':int(150 * 5*10**2) }, replicas=replicas,
                              #hot=400,hill_height=2, eff_threshold=0.3,cold=278,exchange_period=50)
pte_result = ps.pte_replica(mpi_np=MPI_NP, max_tries=4,min_iters=2,
                              mdp_kwargs={'nsteps':int(150 * 5*10**2) }, replicas=replicas,
                              hot=400,hill_height=2, eff_threshold=0.3,cold=278,exchange_period=200)
with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)
pte_plumed_script=pte_result['plumed']
replica_temps=pte_result['temperatures']
temps = ','.join(str(e) for e in replica_temps)
kwargs = [ {'ref_t': ti} for ti in replica_temps]
#now reload PTE_WTE hills file and run eds with cs2backbone to generate the eds parameters with multiple replicas
plumed_input0=textwrap.dedent(
    '''

    peptide: GROUP ATOMS=1-{}                                                                                      
    WHOLEMOLECULES ENTITY0=peptide                                                                                 
    cs: CS2BACKBONE ATOMS=peptide DATA=data NRES={}                                                           

    #bias the simalation with EDS and chem chifts until one gets good eds convergence
     
    PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=200                     
    exp_mean: MATHEVAL ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_7,cs.exphn_8,cs.expha_1,cs.expha_2,cs.expha_3,cs.expha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    calc_mean: MATHEVAL ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    num1: MATHEVAL ARG=cs.exphn_1,exp_mean,cs.hn_1,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm1: MATHEVAL ARG=cs.exphn_1,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num2: MATHEVAL ARG=cs.exphn_2,exp_mean,cs.hn_2,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm2: MATHEVAL ARG=cs.exphn_2,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num3: MATHEVAL ARG=cs.exphn_3,exp_mean,cs.hn_3,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm3: MATHEVAL ARG=cs.exphn_3,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num4: MATHEVAL ARG=cs.exphn_4,exp_mean,cs.hn_4,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm4: MATHEVAL ARG=cs.exphn_4,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num5: MATHEVAL ARG=cs.exphn_7,exp_mean,cs.hn_7,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm5: MATHEVAL ARG=cs.exphn_7,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num6: MATHEVAL ARG=cs.exphn_8,exp_mean,cs.hn_8,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm6: MATHEVAL ARG=cs.exphn_8,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num7: MATHEVAL ARG=cs.expha_1,exp_mean,cs.ha_1,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm7: MATHEVAL ARG=cs.expha_1,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num8: MATHEVAL ARG=cs.expha_2,exp_mean,cs.ha_2,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm8: MATHEVAL ARG=cs.expha_2,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num9: MATHEVAL ARG=cs.expha_3,exp_mean,cs.ha_3,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm9: MATHEVAL ARG=cs.expha_3,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num10: MATHEVAL ARG=cs.expha_5,exp_mean,cs.ha_5,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm10: MATHEVAL ARG=cs.expha_5,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    sum_num: COMBINE ARG=num1,num2,num3,num4,num5,num6,num7,num8,num9,num10 POWERS=1,1,1,1,1,1,1,1,1,1 PERIODIC=NO
    sum_dm: COMBINE ARG=dm1,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10 POWERS=1,1,1,1,1,1,1,1,1,1 PERIODIC=NO
    slope_cumul: MATHEVAL ARG=sum_num,sum_dm FUNC=x/y PERIODIC=NO
    int_cumul: MATHEVAL ARG=slope_cumul,exp_mean,calc_mean FUNC=z-x*y PERIODIC=NO
    slope: AVERAGE ARG=slope_cumul
    int: AVERAGE ARG=int_cumul 
    PRINT ARG=slope_cumul,int_cumul,slope,int FILE=slope_intercept.dat 
    cal_cs1: MATHEVAL ARG=slope,int,cs.exphn_1 FUNC=x*z+y PERIODIC=NO
    cal_cs2: MATHEVAL ARG=slope,int,cs.exphn_2 FUNC=x*z+y PERIODIC=NO
    cal_cs3: MATHEVAL ARG=slope,int,cs.exphn_3 FUNC=x*z+y PERIODIC=NO
    cal_cs4: MATHEVAL ARG=slope,int,cs.exphn_4 FUNC=x*z+y PERIODIC=NO
    cal_cs5: MATHEVAL ARG=slope,int,cs.exphn_7 FUNC=x*z+y PERIODIC=NO
    cal_cs6: MATHEVAL ARG=slope,int,cs.exphn_8 FUNC=x*z+y PERIODIC=NO
    cal_cs7: MATHEVAL ARG=slope,int,cs.expha_1 FUNC=x*z+y PERIODIC=NO
    cal_cs8: MATHEVAL ARG=slope,int,cs.expha_2 FUNC=x*z+y PERIODIC=NO
    cal_cs9: MATHEVAL ARG=slope,int,cs.expha_3 FUNC=x*z+y PERIODIC=NO
    cal_cs10: MATHEVAL ARG=slope,int,cs.expha_5 FUNC=x*z+y PERIODIC=NO
    eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 CENTER_ARG=cal_cs1,cal_cs2,cal_cs3,cal_cs4,cal_cs5,cal_cs6,cal_cs7,cal_cs8,cal_cs9,cal_cs10 PERIOD={} TEMP=@replicas:{{{}}} COVAR OUT_RESTART=restart_pt_wte.dat MULTI_PROP=0.1 RANGE=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1   
    PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat
    PRINT ARG=cal_cs1,cal_cs2,cal_cs3,cal_cs4,cal_cs5,cal_cs6,cal_cs7,cal_cs8,cal_cs9,cal_cs10 STRIDE=200 FILE=calibrated_shifts.dat
    phi_27ASN_28LYS: TORSION ATOMS=86,88,90,108#phi 27ASN to 28LYS bend                                            
    psi_28LYS_29GLY: TORSION ATOMS=88,90,108,110#psi 28 LYS and 29GLY                                              

    distance_SER26HN_asp23: DISTANCE ATOMS=com_ser_26,com_asp_23#SER26HN and asp23

    distance_ASN27_glu22: DISTANCE ATOMS=com_asn_27,com_glu_22#asn27 and glu22                                     
    
    coor_number_K28_E22: COORDINATION GROUPA=23-25 GROUPB=104-107 R_0=0.8
    coor_number_K28_D23: COORDINATION GROUPA=35-37 GROUPB=104-107 R_0=1.4
    coor_number_N27_E22: COORDINATION GROUPA=22-25 GROUPB=83-85 R_0=1.4
    coor_number_N27_D23: COORDINATION GROUPA=35-37 GROUPB=83-85 R_0=1.5        
    PRINT ARG=phi_27ASN_28LYS,psi_28LYS_29GLY,distance_ASN27_asp23,distance_ASN27_glu22,coor_number_K28_E22,coor_number_K28_D23,pmf_coor_number_N27_D23,coor_number_N27_E22,distance_val24_A30 FILE={} STRIDE=250     

    ENDPLUMED
    '''.format(total_no_atoms,number_chains,eds_period,temps))
plumed_input=pte_plumed_script+plumed_input0
with open('plumed_eds_conver_pt_wte_metad.dat', 'w') as f:
    f.write(plumed_input)
ps.add_file('plumed_eds_conver_pt_wte_metad.dat')
time_ns=0.10#10
replex_eff = 0
for kw in kwargs:
    kw['nsteps']=  int(time_ns * 5 * 10 ** 5)
max_iterations=6
min_iterations=4
for i in range(max_iterations):
    with open(pickle_name, 'w') as f:
        pickle.dump(ps, file=f)
    ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_conver_eds_{}'.format(i),  mdp_kwargs=kwargs,run_kwargs={'plumed':'plumed_eds_conver_pt_wte_metad.dat', 'replex': remd_exchange_period}, pickle_name=pickle_name,mpi_np=MPI_NP)
    
    replex_eff = min(get_replex_e(ps, replicas))
    if (replex_eff >= 0.28 and i>=min_iterations):
        print 'Reached replica exchange efficiency of {}. Continuing to production'.format(replex_eff)
        break
    else:
        print 'Replica exchange efficiency of {}. Continuing simulation'.format(replex_eff)


eds_conver_ptwte_folder=ps.sims[-1].location
print eds_conver_ptwte_folder,'eds_conver_ptwte_folder'


colvar_file='{}/restart_pt_wte.0.dat'.format(os.path.abspath(eds_conver_ptwte_folder))

plumed_input=textwrap.dedent(
'''
    peptide: GROUP ATOMS=1-{}                                                                                      
    WHOLEMOLECULES ENTITY0=peptide                                                                                 
    cs: CS2BACKBONE ATOMS=peptide DATA=data NRES={}                                                           

    #bias the simalation with EDS and chem chifts until one gets good eds convergence
     
    PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=CS_shifts STRIDE=200                     
    exp_mean: MATHEVAL ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_7,cs.exphn_8,cs.expha_1,cs.expha_2,cs.expha_3,cs.expha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    calc_mean: MATHEVAL ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    num1: MATHEVAL ARG=cs.exphn_1,exp_mean,cs.hn_1,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm1: MATHEVAL ARG=cs.exphn_1,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num2: MATHEVAL ARG=cs.exphn_2,exp_mean,cs.hn_2,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm2: MATHEVAL ARG=cs.exphn_2,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num3: MATHEVAL ARG=cs.exphn_3,exp_mean,cs.hn_3,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm3: MATHEVAL ARG=cs.exphn_3,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num4: MATHEVAL ARG=cs.exphn_4,exp_mean,cs.hn_4,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm4: MATHEVAL ARG=cs.exphn_4,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num5: MATHEVAL ARG=cs.exphn_7,exp_mean,cs.hn_7,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm5: MATHEVAL ARG=cs.exphn_7,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num6: MATHEVAL ARG=cs.exphn_8,exp_mean,cs.hn_8,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm6: MATHEVAL ARG=cs.exphn_8,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num7: MATHEVAL ARG=cs.expha_1,exp_mean,cs.ha_1,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm7: MATHEVAL ARG=cs.expha_1,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num8: MATHEVAL ARG=cs.expha_2,exp_mean,cs.ha_2,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm8: MATHEVAL ARG=cs.expha_2,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num9: MATHEVAL ARG=cs.expha_3,exp_mean,cs.ha_3,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm9: MATHEVAL ARG=cs.expha_3,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num10: MATHEVAL ARG=cs.expha_5,exp_mean,cs.ha_5,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm10: MATHEVAL ARG=cs.expha_5,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    sum_num: COMBINE ARG=num1,num2,num3,num4,num5,num6,num7,num8,num9,num10 POWERS=1,1,1,1,1,1,1,1,1,1 PERIODIC=NO
    sum_dm: COMBINE ARG=dm1,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10 POWERS=1,1,1,1,1,1,1,1,1,1 PERIODIC=NO
    slope_cumul: MATHEVAL ARG=sum_num,sum_dm FUNC=x/y PERIODIC=NO
    int_cumul: MATHEVAL ARG=slope_cumul,exp_mean,calc_mean FUNC=z-x*y PERIODIC=NO
    slope: AVERAGE ARG=slope_cumul CLEAR=200#moving average of last 200
    int: AVERAGE ARG=int_cumul CLEAR=200#moving average of last 200
    cal_cs1: MATHEVAL ARG=slope,int,cs.exphn_1 FUNC=x*z+y PERIODIC=NO
    cal_cs2: MATHEVAL ARG=slope,int,cs.exphn_2 FUNC=x*z+y PERIODIC=NO
    cal_cs3: MATHEVAL ARG=slope,int,cs.exphn_3 FUNC=x*z+y PERIODIC=NO
    cal_cs4: MATHEVAL ARG=slope,int,cs.exphn_4 FUNC=x*z+y PERIODIC=NO
    cal_cs5: MATHEVAL ARG=slope,int,cs.exphn_7 FUNC=x*z+y PERIODIC=NO
    cal_cs6: MATHEVAL ARG=slope,int,cs.exphn_8 FUNC=x*z+y PERIODIC=NO
    cal_cs7: MATHEVAL ARG=slope,int,cs.expha_1 FUNC=x*z+y PERIODIC=NO
    cal_cs8: MATHEVAL ARG=slope,int,cs.expha_2 FUNC=x*z+y PERIODIC=NO
    cal_cs9: MATHEVAL ARG=slope,int,cs.expha_3 FUNC=x*z+y PERIODIC=NO
    cal_cs10: MATHEVAL ARG=slope,int,cs.expha_5 FUNC=x*z+y PERIODIC=NO
    eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 CENTER_ARG=cal_cs1,cal_cs2,cal_cs3,cal_cs4,cal_cs5,cal_cs6,cal_cs7,cal_cs8,cal_cs9,cal_cs10 TEMP=278 IN_RESTART={} MULTI_PROP=0.1 RANGE=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 FREEZE MEAN PERIOD=50000
    phi_27ASN_28LYS: TORSION ATOMS=86,88,90,108#phi 27ASN to 28LYS bend                                            
    psi_28LYS_29GLY: TORSION ATOMS=88,90,108,110#psi 28 LYS and 29GLY                                              
    phi_ALA1_GLU2: TORSION ATOMS=11,13,14,26#phi ALA and GLU                                                       
    psi_ALA1_GLU2: TORSION ATOMS=1,5,11,13#psi ALA and GLU                                                         
    com_glu_22: COM ATOMS=23,24,25#coo                                                                             
    com_lys_28: COM ATOMS=104,105,106,107#nh3                                       group_lys_28: ATOMS=104,105,106,107#nh3                                                                      
    com_asp_23: COM ATOMS=35,36,37#coo                                                                             
    com_val_24: COM ATOMS=44,45,46,47,48,49,50,51,52,53#ch2cch3                                                    
    com_asn_27: COM ATOMS=81,82,83,84,85#conh2
    com_ser_26: COM ATOMS=71,70  #OH       
    com_ala_30: COM ATOMS=121,122,123,124                      
    distance_val24_A30: DISTANCE ATOMS=com_ala_30,com_val_24#val24 and ala30
    distance_SER26HN_glu22: DISTANCE ATOMS=com_ser_26,com_glu_22#SER26HN and glu22

    distance_SER26HN_asp23: DISTANCE ATOMS=com_ser_26,com_asp_23#SER26HN and asp23

    distance_ASN27_glu22: DISTANCE ATOMS=com_asn_27,com_glu_22#asn27 and glu22                                     
    distance_ASN27_asp23: DISTANCE ATOMS=com_asn_27,com_asp_23#asn27 and asp23
    coor_number_K28_E22: COORDINATION GROUPA=23-25 GROUPB=104-107 R_0=0.8
    coor_number_K28_D23: COORDINATION GROUPA=35-37 GROUPB=104-107 R_0=1.4
    coor_number_N27_E22: COORDINATION GROUPA=22-25 GROUPB=83-85 R_0=1.4
    coor_number_N27_D23: COORDINATION GROUPA=35-37 GROUPB=83-85 R_0=1.5        
    PRINT ARG=phi_27ASN_28LYS,psi_28LYS_29GLY,phi_ALA1_GLU2,psi_ALA1_GLU2,distance_SER26HN_glu22,distance_SER26HN_asp23,distance_ASN27_asp23,distance_ASN27_glu22,coor_number_K28_E22,coor_number_K28_D23,pmf_coor_number_N27_D23,coor_number_N27_E22,distance_val24_A30 FILE={} STRIDE=250                                       
    PRINT ARG=eds.bias,eds.force2 FILE=eds_output.dat STRIDE=100
    PRINT ARG=cal_cs1,cal_cs2,cal_cs3,cal_cs4,cal_cs5,cal_cs6,cal_cs7,cal_cs8,cal_cs9,cal_cs10 STRIDE=200 FILE=calibrated_shifts.dat
    ENDPLUMED'''.format(total_no_atoms,number_chains,colvar_file,'OUTPUT_COLVAR_EDS'))


with open('plumed_eds_colvars.dat', 'w') as f:
    f.write(plumed_input)
ps.add_file('plumed_eds_colvars.dat')
final_time_eds=int(0.1*0.05*10**5)#int(40000*0.05*10**5)
ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod_eds_colvar',  mdp_kwargs={'nsteps': final_time_eds, 'ref_t': 278},run_kwargs={'plumed':'plumed_eds_colvars.dat'},mpi_np=MPI_NP, pickle_name=pickle_name)
#finally:
with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)
  
print ps.box_size_angstrom, replica_temps

print 'done'

