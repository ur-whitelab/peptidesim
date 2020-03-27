import signal
import functools
import os
import pandas as pd
import numpy as np

class TimeoutError(Exception):
    pass


def timeout(seconds, error_message='Function call timed out'):
    def decorated(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return functools.wraps(func)(wrapper)

    return decorated


# http://stackoverflow.com/questions/3812849/how-to-check-whether-a-directory-is-a-sub-directory-of-another-directory
def in_directory(file, directory, allow_symlink=False):
    # make both absolute
    directory = os.path.abspath(directory)
    file = os.path.abspath(file)

    # check whether file is a symbolic link, if yes, return false if they are not allowed
    if not allow_symlink and os.path.islink(file):
        return False

    # return true, if the common prefix of both is equal to directory
    # e.g. /a/b/c/d.rst and directory is /a/b, the common prefix is /a/b
    return os.path.commonprefix([file, directory]) == directory


def validity_check(df, shifts):
    Inval_index = []
    Shif_index = []
    count = 0
    for i in range(len(shifts)):
        Index = shifts[i]
        Shif_index.append(Index+1)
        if df.values[Index][1] == 1:
            Inval_index.append(Index+1)
            count += 1
    return (count, Inval_index, Shif_index)


def cs_validity(plumed_dat, exec_dir):
    '''plumed_dat = The path to the plumed file (include the plumed data filename and extension)
       exec_dir = The path to the cs_shifts dat files)'''
    with open(plumed_dat, 'r') as file_plumed:
        shifts_plumed = parser(file_plumed)
    keys = ['CAshifts', 'CBshifts', 'Cshifts',
            'HAshifts', 'Hshifts', 'Nshifts']
    shifts_plumed_dict = dict.fromkeys(keys)
    dat_files = [f for f in os.listdir(exec_dir) if f.endswith('.dat')]
    for i in range(len(shifts_plumed)):
        shifts_dat = exec_dir + dat_files[i]
        df_shifts_dat = pd.read_csv(shifts_dat, sep=' ', header=None)
        shifts_plumed_dict[keys[i]] = shifts_plumed[i]
        if shifts_plumed[i] != []:
            count, Invalid_index, Shifts_index = validity_check(
                df_shifts_dat, shifts_plumed_dict[keys[i]])
            if count != 0:
                raise ValueError('Simulation uses data {} for {}. Invalid EDS chemical shift found for data {}!'.format(
                    Shifts_index, keys[i], Invalid_index))
    return True


def parser(plumed_file):
    # Finding plumed data
    string = 'eds: EDS ARG=cs'
    for line in plumed_file:
        if string in line:
            txt = line
    # Skipping 4 letters to account for 'ARG=cs'
    txt = txt.split(' ')
    txt = txt[2][4:]
    txt_list = txt.split(',')
    hn_shifts = []
    nh_shifts = []
    ca_shifts = []
    cb_shifts = []
    c_shifts = []
    ha_shifts = []
    for i in range(len(txt_list)):
        cs_type = txt_list[i][0:5]
        if cs_type == 'cs.hn':
            hn_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == 'cs.nh':
            nh_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == 'cs.ca':
            ca_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == 'cs.cb':
            cb_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == 'cs.c_':
            c_shifts.append(int(txt_list[i][5:])-1)
        elif cs_type == 'cs.ha':
            ha_shifts.append(int(txt_list[i][6:])-1)
    return (ca_shifts, cb_shifts, c_shifts, ha_shifts, hn_shifts, nh_shifts)


def pdb_for_plumed(input_file, peptide_copies,
                   atoms_in_chain, first_atom_index,
                   output_file):
    ''' Funtion that takes an old pdbfile that has all hydrogens
        but without unique chain IDs and without terminii of chains
        indicated and generates a new pdbfile with unique chain IDs
        and terminii indicated.

        Parameters
        ----------
        input_file: input pdb filename
        number_of_chains: list of number of copies of each sequence
        atoms_in_chain: list of number of atoms in each sequence
        first_atom_index: the line number containing first atom
                        in the old pdbfile (line indexing starts at 0)
        output_file: output pdb filename

        Returns
        -------
        output_pdbfile

        Example function call
        -------------
        pdb_for_plumed(input_file='template.pdb',
                        peptide_copies=[6,2],
                        atoms_in_chain=[35,43],
                        first_atom_index=5,
                        output_file='new.pdb')
    '''
    from string import ascii_uppercase

    # read the pdb file
    with open(input_pdbfile, 'r') as f:
        lines = f.readlines()
        # grab lines that don't matter
        beginning = lines[:first_atom_index]
        # grab lines that need to be changed
        lines = lines[first_atom_index:]

        with open(output_pdbfile, 'w') as f:
            # iterate through first few lines and write them as is
            for index in np.arange(len(beginning)):
                f.write("{}\n".format(beginning[index].strip()))

            skip_lines = 0
            atoms_scanned = 0
            # unique chain ID
            letter = 0
            # iterate through the number of different sequences
            for seq in np.arange(len(peptide_copies)):
                if seq != 0:
                    letter += 1
                # iterate through the copies of that sequence
                for copy in np.arange(peptide_copies[seq]):
                    if copy != 0:
                        letter += 1
                    # iterate through the atoms in the chain
                    for atom in np.arange(atoms_in_chain[seq][copy]):
                        current_line = lines[skip_lines + copy *
                                             (atoms_in_chain[seq][copy]) + atom]
                        # converts the string into a list of characters
                        split_line = list(current_line)
                        # unique ID on at position 21 of current_line
                        if (len(split_line) >= 21):
                            split_line[21] = ascii_uppercase[letter]
                            current_line = "".join(split_line)
                            f.write('{}'.format(current_line))
                        if(atom == int(atoms_in_chain[seq][copy]-1)):
                            f.write('TER second\n')
                            skip_lines += 1
                            # puts ter at the end of each chain
                    atoms_scanned += atoms_in_chain[seq][copy]
                    print(atoms_scanned, skip_lines)
                skip_lines += atoms_scanned
            f.write('ENDMDL second\n')
            f.close()

    return output_pdbfile

def get_atoms_in_chains(input_file):
    ''' This function returns the number of atoms in each peptide
    chain. Note that the solvent molecules are excluded from counting.

    Parameters
    ----------
    input_pdbfile: input pdb filename to use to count atoms

    Returns
    -------
    atoms_in_chains

    Example function call
    -------------
    get_atoms_in_chains(input_pdbfile='template.pdb')
    '''
    with open(input_file, "r") as f:
        lines = f.readlines()
    # use rest of the lines to get atoms in chain
    atoms_in_chains = []
    old_line = lines[4]
    for line in lines[4:]:
        if line == 'TER\n':
            split = old_line.split()
            atoms_in_chains.append(int(split[1])-sum(atoms_in_chains))
        # do not count solvent atoms
        # solvent molecules are always at the end
        if 'SOL' in line:
            break
        old_line = line
    return atoms_in_chains

def remove_solvent_from_pdb(input_file,
                            output_file='template_no_solvent.pdb'):
    ''' Remove solvent molecules from the pdb file

    Parameters
    ----------
    input_file: input pdb filename
    output_file: output pdb filename

    Returns
    -------
    output_file

    Example function call
    -------------
    remove_solvent_from_pdb(input_file='template.pdb',
                            output_file='no_solvent.pdb')
    '''
    with open(input_file, "r") as f:
        lines = f.readlines()
    with open(output_file, 'w') as output:
        # write back first few lines as is 
        for line in lines[:4]:
            output.write(line)
        # Then check if there is solvent
        for line in lines[4:]:
            if('SOL' in line):
                output.write('TER final\n')
                output.write('ENDMDL final\n')
                output.close()
                break
            output.write(line)
    return output_file