import numpy as np
import matplotlib.pyplot as plt
from peptidesim import PeptideSim, utilities
import textwrap
import sys
import re
import os
import sys
import gromacs
import dill as pickle
from shutil import copyfile, move

import matplotlib
matplotlib.use('Agg')

name = sys.argv[2]
debug = False
MPI_NP = 1
data_dir = sys.argv[4]
if sys.argv[3][0]=='[':
    seq=sys.argv[1].strip('[').strip(']').split(',')
    copies = sys.argv[3].strip('[').strip(']').split(',')
    peptide_copies =[int(i) for i in copies]
else:
    peptide_copies=[int(sys.argv[3])]
    seq=[sys.argv[1]]

ps = PeptideSim(name, seq, peptide_copies, job_name='{}'.format(name))
ps.mdrun_driver = 'gmx_mpi'
ps.forcefield = 'amber99sb'
ps.water = 'tip4p'
ps.peptide_density = 0.008  # mg/ml
ps.ion_concentration = 0.001  # 10mM
ps.initialize()

file00 = ps.pdb_file
cwd = os.getcwd()
# gromacs.pdb2gmx(f=file00,o="file_pdb_output.pdb",water=ps.water,ff=ps.forcefield)
# file00="{}/{}".format(cwd,"file_pdb_output.pdb")

def atoms_in_each_chain(input_file, output_file='template.pdb'):
    output = open(output_file, 'w')
    with open(input_file, "r") as f:
        lines = f.readlines()
        i = 4 
        for line in lines[:4]:
            output.write(line)

        for line in lines[4:]:
            old_line = line
            line = line.split()
            i = i + 1
            if(len(line) > 1 and line[3] == 'SOL'):
                output.write('TER final\n')
                output.write('ENDMDL final\n')
                output.close()
                break
            # i=i+1
            output.write(old_line)
        last_line = []
        print(lines[i - 1], lines[i - 2], lines[i - 3])
        if (lines[i - 1].startswith('END')):

            last_line = lines[i - 3].strip()
            print("ends with END")
        else:
            last_line = lines[i - 1].strip()
            print("doesnt end with END")

        last_line = last_line.split()

        output.close()
        if (last_line[4].isdigit()):
            return int(last_line[1]), output_file
        else:
            return int(last_line[1]), output_file

total_no_atoms, file0 = total_aa(file00, 'template.pdb')
#find the total number of chains
number_chains = int(sum(peptide_copies)
atoms_in_chain = int(total_no_atoms / number_chains)

# replace 'HIS' as 'HIE' because plumed likes it that way
with open(file0, 'r') as file_read:
    file_data = file_read.read()
file_data = file_data.replace('HIS', 'HIE')
with open(file0, 'w') as file_write:
    file_write.write(file_data)

file3 = 'template.pdb'
print('file0 is:{}'.format(file0))
new_pdb_file = utilities.pdb_for_plumed(input_pdbfile=file0,
                                        peptide_copies=peptide_copies,
                                        atoms_in_chain=atoms_in_chain,
                                        first_atom_index=5,
                                        output_pdbfile=file3)

directory = 'data'
if not os.path.exists(directory):
    os.mkdir(directory)
template_pdb = os.path.join(directory, file3)

camshift_src = os.path.join(data_folder + 'camshift.db')
gromacs_a03_mdb_src = os.path.join(data_folder + 'a03_gromacs.mdb')
h_src = os.path.join(data_folder + 'Hshifts.dat')
ha_src = os.path.join(data_folder + 'HAshifts.dat')
# c_src=os.path.join(data_folder + 'Cshifts.dat')
# ca_src=os.path.join(data_folder + 'CAshifts.dat')
# cb_src=os.path.join(data_folder + 'CBshifts.dat')
# n_src=os.path.join(data_folder + 'Nshifts.dat')
gromacs_a03_mdb_dst = os.path.join(directory, 'a03_gromacs.mdb')
camshift_dst = os.path.join(directory, 'camshift.db')
template_pdb = os.path.join(directory, file3)
h_dst = os.path.join(directory, 'Hshifts.dat')
ha_dst = os.path.join(directory, 'HAshifts.dat')
# c_dst=os.path.join(directory,'Cshifts.dat')
# ca_dst=os.path.join(directory,'CAshifts.dat')
# cb_dst=os.path.join(directory,'CBshifts.dat')
# n_dst=os.path.join(directory,'Nshifts.dat')
copyfile(h_src, h_dst)
copyfile(ha_src, ha_dst)
# copyfile(c_src,c_dst)
# copyfile(ca_src,ca_dst)
# copyfile(cb_src,cb_dst)
# copyfile(n_src,n_dst)
move(new_pdb_file, template_pdb)
copyfile(camshift_src, camshift_dst)
copyfile(gromacs_a03_mdb_src, gromacs_a03_mdb_dst)

center = ''

def data_folder(number_amino_acids, name, copies_chains):
    ''' a function that takes total number of amino acids in the pdbfile,
    name of the chemical shifts file and how many copies of a peptide are
    present and puts chemical shifta file into the data folder
    Var:
       number_amino_acids-total number of amino acids in the pdbfile
       name-name of the chemical shifts file
       copies_chains-the number of peptide copies
       '''
    new_file = os.path.join(
        directory,
        name)  # creates a chemical shift file in data folder
    # checks whether the exp chem shifts exist, for now only H shifts are
    # available
    if (name != 'HAshifts.dat' and name != 'Hshifts.dat'):

        # creates a new file to copy over the chemical shifts
        with open(new_file, 'w') as f:
            # amino acid residue indeces for beginning and end of
            # each chain of peptides
            beg_end_chains = []
            # number of amino acids in one chain
            number_chains = int(number_amino_acids / copies_chains)
            for i in range(
                    int(copies_chains)):  # iterates throug a number of copies
                # adds the begining amino acid index of the chain i, first
                # index starts at 1
                beg_end_chains.append(i * number_chains + 1)
                # adds the amino acid ending index of the chain
                beg_end_chains.append(i * number_chains + number_chains)
                # iterates through the amino acid indeces in the chain
                for j in range(int(number_chains)):
                    j = j + 1  # amino acid index also starts at 1
                    # checks whether the amino acid index is at the beginning
                    # or ending of the chain
                    if (j in beg_end_chains):
                        # puts a hash sign for beginning and amino acid indeces
                        # in the chain and 1 for lack of exp chem shift
                        f.write('#{} 1\n'.format(number_chains * i + j))
                    else:
                        # puts 1 for lack of exp chem shifts
                        f.write('{} 1\n'.format(number_chains * i + j))
    else:
        with open(new_file, 'w') as f:
            # amino acid residue indeces for beginning and end of each chain of
            # peptides
            beg_end_chains = []
            # number of amino acids in one chain
            number_chains = int(number_amino_acids / copies_chains)
            chemical_shifts = np.genfromtxt(
                os.path.join(data_folder, '{}'.format(name)),
                comments='!#')
#            print chemical_shifts
            for i in range(copies_chains):
                # adds the begining amino acid index of the chain i, first
                # index starts at 1
                beg_end_chains.append(i * number_chains + 1)
                # adds the amino acid ending index of the chain
                beg_end_chains.append(i * number_chains + number_chains)
                # print beg_end_chains
                # iterates through the amino acid indeces in the chain
                for j in range(number_chains):
                    j = j + 1  # amino acid index also starts at 1
                    # checks whether the amino acid index is at
                    # the beginning or ending of the chain
                    if (j in beg_end_chains):
                        # print
                        # i,j,number_chains*i+j,chemical_shifts[j-1][1],a
                        # puts a hash sign for beginning and amino acid indeces
                        # in the chain
                        f.write('#{} {}\n'.format(
                            number_chains * i + j,
                            chemical_shifts[j - 1][1]))
                    else:
                        # print i,j,number_chains*i+j,chemical_shifts[j-1][1]
                        # copies the chem shifts from file with chemical shifts
                        # for the correspoding amino acid
                        f.write('{} {}\n'.format(
                            number_chains * i + j,
                            chemical_shifts[j - 1][1]))


filenames = [
    'Cshifts.dat',
    'CAshifts.dat',
    'HAshifts.dat',
    'Hshifts.dat',
    'CBshifts.dat',
    'Nshifts.dat']
for i in filenames:
    data_folder(int(number_chains), i, int(peptide_copies))
ps.add_file(directory)
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

ps.run(
    mdpfile='peptidesim_emin.mdp',
    tag='init_emin',
    mdp_kwargs={
        'nsteps': 8 *
        10**2,
        'rcoulomb': 1},
    mpi_np=MPI_NP)

ps.run(
    mdpfile='peptidesim_anneal.mdp',
    tag='annealing',
    mdp_kwargs={
        'nsteps': int(
            5 * 5 * 10**2)},
    run_kwargs={
        'cpt': 5},
    mpi_np=MPI_NP)  # change the time step to 2 ns

ps.run(
    mdpfile='peptidesim_npt.mdp',
    tag='equil_npt',
    mdp_kwargs={
        'nsteps': int(
            5 * 5 * 10**2),
        'ref_t': 278},
    run_kwargs={
        'cpt': 5},
    mpi_np=MPI_NP)
