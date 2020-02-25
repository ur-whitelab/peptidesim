import numpy as np
import matplotlib.pyplot as plt
from peptidesim import PeptideSim
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
seq = sys.argv[1]
name = sys.argv[2]
debug = False
pickle_name = name + '.pickle'
MPI_NP = 4
peptide_copies = int(sys.argv[3])
data_folder = sys.argv[4]

# try to reload
if(os.path.exists(pickle_name)):
    print('loading restart')
    with open(pickle_name, 'rb') as f:
        ps = pickle.load(f)
else:
    ps = PeptideSim(name, [seq], [peptide_copies], job_name='{}'.format(name))
    ps.mdrun_driver = 'gmx_mpi'
    ps.forcefield = 'amber99sb'
    ps.water = 'tip4p'
    ps.peptide_density = 0.008  # mg/ml
    ps.ion_concentration = 0.001  # 10mM
    ps.initialize()
    with open(ps.pickle_name, 'wb') as f:
        pickle.dump(ps, file=f)

file00 = ps.pdb_file
cwd = os.getcwd()
# gromacs.pdb2gmx(f=file00,o="file_pdb_output.pdb",water=ps.water,ff=ps.forcefield)
# file00="{}/{}".format(cwd,"file_pdb_output.pdb")
# print(file00,"file_pdb_output.pdb", cwd)


def total_aa(file1, output_file):
    print(ps.sims[-1].location)
    output = open(output_file, 'w')
    with open(file1, "r") as f:
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
            return int(last_line[1]), int(last_line[4]), output_file
        else:
            return int(last_line[1]), int(last_line[5]), output_file
# print total_aa(file00,'template.pdb')


def number_all_atoms_chains():
    if (os.path.isdir("{}/".format(os.getcwd) + "data")):
        template_file = "{}/".format(os.getcwd) + "data/template.pdb"
        # template_file="{}/".format(os.getcwd)+"data/template.pdb"
        with open(template_file, 'r') as f:
            lines = f.readlines()
            lines = lines[-3].strip()
            lines = lines.split()
            return int(lines[1]), int(last_lines[5])


total_no_atoms, number_chains, file0 = total_aa(file00, 'template.pdb')
number_chains = int(number_chains)
total_no_atoms = int(total_no_atoms)
peptide_copies = int(peptide_copies)
atoms_in_chain = int(total_no_atoms / peptide_copies)


with open(file0, 'r') as file_read:
    file_data = file_read.read()
file_data = file_data.replace('HIS', 'HIE')
with open(file0, 'w') as file_write:
    file_write.write(file_data)


def pdbfile_generator_w_chain_id(
        number_of_chains,
        atoms_in_chain,
        first_atom_index,
        output_pdbfile,
        input_pdbfile):
    '''a funtion that takes an old pdbfile that has all hydrogens
    but without unique chain IDs and without terminii of chains
    indicated and generates a new pdbfile with uniqe chain id's
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
    print(input_pdbfile, number_of_chains, atoms_in_chain)
    with open(input_pdbfile, 'r') as f:
        lines = f.readlines()
        # saves the first useless lines that don't contain conf info
        beginning = lines[:first_atom_index]
        # gets rid of those lines from readlines
        lines = lines[first_atom_index:]
        first_line = lines[0].strip()
        first_line = first_line.split()
        with open(output_pdbfile, 'w') as f:
            for index in np.arange(
                    len(beginning)):  # iterates through useless lines
                # writes the useless lines into the new file
                f.write("{}\n".format(beginning[index].strip()))
            # iterates through the copies of chains
            for i in np.arange(number_of_chains):
                # itarates through the atoms in the chain
                for j in np.arange(atoms_in_chain):
                    print(i * (atoms_in_chain + 1) + j, i *
                          (atoms_in_chain) + j, len(lines), i, j)
                    # reads the old pdbfile info pertaining to the atoms of
                    # interest
                    a = lines[i * (atoms_in_chain) + i + j]
                    # converts the string into a list of characters
                    a = list(a)
                    if (len(a) >= 21):
                        # replaces the repetittive chain ID with unique chain
                        # ID chosen from
                        a[21] = ascii_uppercase[i]
                    # puts the list of characters back to a
                    # string that contains line info
                    a = "".join(a)
                    # writes the line into the new pdbfile
                    f.write('{}'.format(a))

                    # if(int(i*(atoms_in_chain-1)+i+j)==int((atoms_in_chain-1)*i+atoms_in_chain-1+i)
                    # and peptide_copies!=1): #checks whether the atoms is at
                    # the end of the chain
                    # print
                    # int(i*atoms_in_chain+j),int((atoms_in_chain*i+atoms_in_chain-1))
                    if(j == int(atoms_in_chain - 1)):
                        #    print "wrote TER"
                        # puts ter at the end of each chain
                        f.write('TER second\n')
                if (peptide_copies == i + 1):
                    f.write('ENDMDL second\n')  # puts finish touches
                    f.close()  # done
    return output_pdbfile


file3 = 'template.pdb'
new_pdb_file = pdbfile_generator_w_chain_id(
    peptide_copies, atoms_in_chain, 5, 'new_pdb.pdb', file0)

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

number_chains = int(number_chains)
peptide_copies = int(peptide_copies)


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

print(ps.pickle_name, 'picklename1')
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

# print ps.pickle_name, 'picklename2'
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)
ps.run(
    mdpfile='peptidesim_anneal.mdp',
    tag='annealing',
    mdp_kwargs={
        'nsteps': int(
            5 * 5 * 10**2)},
    run_kwargs={
        'cpt': 5},
    mpi_np=MPI_NP)  # change the time step to 2 ns
# ps.run(mdpfile='peptidesim_anneal.mdp',tag='annealing',mdp_kwargs={'nsteps':int(2000*
# 5*10**2)},run_kwargs={'cpt':5},mpi_np=MPI_NP)#change the time step to 2
# ns


print(ps.pickle_name, 'picklename3')
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

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

print(ps.pickle_name, 'picklename3')
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)
