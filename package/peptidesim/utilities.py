import signal
import functools
import os
import pandas as pd
import pkg_resources

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
        elif cs_type == 'cs.c-':
            c_shifts.append(int(txt_list[i][5:])-1)
        elif cs_type == 'cs.ha':
            ha_shifts.append(int(txt_list[i][6:])-1)
    return (ca_shifts, cb_shifts, c_shifts, ha_shifts, hn_shifts, nh_shifts)



def load_eds_restart(filename):
    with open(filename, 'r') as f:
        header = f.readline().split()[2:]
    # initial read of fields
    data = pd.read_table(filename, sep=r'\s+', comment='#', names=header)
    return data


def plot_couplings(eds_filename, output_plot='couplings.png'):
    import matplotlib.pyplot as plt
    data = load_eds_restart(eds_filename)
    # get cv names
    cv_names = set()
    for n in data.columns:
        sn = n.split('_')
        if len(sn) > 1:
            cv_names.add(sn[0])
    fig, ax = plt.subplots(nrows=len(cv_names) // 2, ncols=2, figsize=(12, 8))
    index = 0
    cv_names = list(cv_names)
    cv_names.sort()
    for i in range(len(cv_names) // 2):
        for j in range(2):
            cv = cv_names[index]
            ax[i, j].plot(data.time, data[f'{cv}_coupling'])
            ax[i, j].set_title(cv)
            index += 1
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)


def pdb_for_plumed(input_file, output_file):
    ''' Funtion that takes an old .pdb file that has all hydrogens
        but without unique chain IDs and without terminii of chains
        indicated and generates a new .pdb file with unique chain IDs
        and terminii indicated.

        Parameters
        ----------
        input_file: input pdb filename
        output_file: output pdb filename

        Returns
        -------
        number of atom records written
    '''
    from string import ascii_uppercase
    last_cid = ''
    cindex = 0
    atom_number = 0
    with open(input_file, 'r') as f, open(output_file, 'w') as o:
        for line in f.readlines():
            if line.startswith('ATOM'):
                atom_number += 1
                # add chain id and atom type to 77
                cid = line[21]
                if last_cid != cid:
                    if last_cid:
                        o.write('TER\n')
                        cindex += 1
                    last_cid = cid
                o.write(line[:21] + ascii_uppercase[cindex] + line[22:77] + line[13] + line[78:])
            else:
                o.write(line)
    return atom_number


def prepare_cs_data(ps, shift_dict=None, pte_reweight=False):
    '''Prepare a directory for adding chemical shifts.

    Parameters
    ----------
    ps: peptidesim object
    shift_dict: dictionary containing shifts.
      Keys should be '[peptide id, from 0]-[resid, from 1]-[res code, one character]-[atom name]' and
      value is shift. For example: 4-G-HA: 4.2. Can be None
    pte_reweight: should pte-rewighting be used

    Returns
    ---------
    dictionary containing following keys:
        data_dir: path to directory
        shift_dict: the given shift dictionary
        plumed: the header necessary to use cs2backbone in scripts
        cs2_names: the names for the given shifts (if averaging)
        cs2_values: the given experimental values (if averaging)
    '''
    data_dir = os.path.join(ps.dir_name, 'camshift_data')
    os.makedirs(data_dir, exist_ok=True)
    template_pdb = os.path.join(data_dir, 'template.pdb')
    atom_number = pdb_for_plumed(input_file=ps.pdb_file, output_file=template_pdb)
    # find which nuclei to consider
    if shift_dict is None:
        shift_dict = {}
    shift_files = set()
    for k in shift_dict.keys():
        n = k.split('-')[-1]
        shift_files.add(n)
    # check to make sure all shifts have been seen
    seen_shifts = set()
    cs2_name_conversions = {'h': 'hn', 'n': 'nh', 'c': 'co'}
    cs2_names = dict()
    for rn in list(shift_files):
        rindex = 1
        cindex = 0
        with open(os.path.join(data_dir, f'{rn}shifts.dat'), 'w') as f:
            for i, s in enumerate(ps.sequences):
                for j in range(ps.counts[i]):
                    for k in range(len(s)):
                        shift = 0.0
                        # check for match
                        key = f'{i}-{k+1}-{s[k]}-{rn}'
                        if key in shift_dict:
                            shift = shift_dict[key]
                            seen_shifts.add(key)
                            if key not in cs2_names:
                                cs2_names[key] = []
                            # y tho
                            rn_csname = rn.lower()
                            if rn_csname in cs2_name_conversions:
                                rn_csname = cs2_name_conversions[rn.lower()]
                            cs2_names[key].append(f'cs.{rn_csname}-{cindex}-{rindex}')
                        if k == 0 or k == len(s) - 1:
                            f.write(f'#{rindex} {shift}\n')
                        else:
                            f.write(f'{rindex} {shift}\n')
                        rindex += 1
                    cindex += 1
    # check for missed shifts
    given_shifts = set(shift_dict.keys())
    if given_shifts != seen_shifts:
        raise ValueError('Not all given shifts were assigned. Could not assign' + str(given_shifts - seen_shifts))

    # add extra files
    with open(os.path.join(data_dir, 'camshift.db'), 'wb') as f:
        f.write(pkg_resources.resource_string(__name__, 'templates/' + 'camshift.db'))
    with open(os.path.join(data_dir, 'a03_gromacs.mdb'), 'wb') as f:
        f.write(pkg_resources.resource_string(__name__, 'templates/' + 'a03_gromacs.mdb'))

    plumed_script = f'peptide: GROUP ATOMS=1-{atom_number}\nWHOLEMOLECULES ENTITY0=peptide\n'
    plumed_script += f'cs: CS2BACKBONE ATOMS=peptide DATADIR={os.path.abspath(data_dir)}\n'

    cs2_values = []
    cs2_avg_names = []
    print_arg_list = []

    # now add averaging across chains
    for k in shift_dict.keys():
        n = k.split('-')[-1]
        i = k.split('-')[0]
        # this is for biasing
        arg_list = ','.join(cs2_names[k])
        plumed_script += f'avg-{k}: COMBINE ARG={arg_list} PERIODIC=NO NORMALIZE\n'
        # for computing running means
        plumed_script += f'all-avg-{k}: AVERAGE ARG=avg-{k}'
        if pte_reweight:
            plumed_script += ' LOGWEIGHTS=pte-lw\n'
        else:
            plumed_script += '\n'
        cs2_values.append(shift_dict[k])
        cs2_avg_names.append(f'avg-{k}')
        exp_name = cs2_names[k][0].replace('cs.', 'cs.exp')
        print_arg_list.append(f'avg-{k},all-avg-{k},{exp_name}')

    plumed_script += f'PRINT FILE=cs_shifts.dat ARG={",".join(print_arg_list)} STRIDE=500\n'

    return {'data_dir': data_dir, 'shift_dict': shift_dict, 'plumed': plumed_script, 'cs2_names': cs2_avg_names, 'cs2_values': cs2_values}
