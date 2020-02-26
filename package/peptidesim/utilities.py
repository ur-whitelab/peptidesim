import signal
import functools
import os


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
    """plumed_dat = The path to the plumed file (include the plumed data filename and extension)
       exec_dir = The path to the cs_shifts data files)"""
    file_plumed = open(plumed_dat, 'r')
    shifts_plumed = parser(file_plumed)
    file_plumed.close()
    keys = ['CAshifts', 'CBshifts', 'Cshifts',
            'HAshifts', 'Hshifts', 'Nshifts']
    shifts_plumed_dict = dict.fromkeys(keys)
    dat_files = [f for f in os.listdir(exec_dir) if f.endswith('.dat')]
    for i in range(len(shifts_plumed)):
        shifts_dat = exec_dir + dat_files[i]
        df_shifts_dat = pd.read_csv(shifts_dat, sep=" ", header=None)
        shifts_plumed_dict[keys[i]] = shifts_plumed[i]
        if shifts_plumed[i] != []:
            count, Invalid_index, Shifts_index = validity_check(
                df_shifts_dat, shifts_plumed_dict[keys[i]])
            if count != 0:
                raise ValueError('Simulation uses data {} for {}. Invalid EDS chemical shift found for data {}!'.format(
                    Shifts_index, keys[i], Invalid_index))
    print('All chemical shifts are valid'.)


def parser(plumed_file):
    # Finding plumed data
    string = 'eds: EDS ARG=cs'
    for line in plumed_file:
        if string in line:
            txt = line
    # Skipping 4 letters to account for "ARG=cs"
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
        if cs_type == "cs.hn":
            hn_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == "cs.nh":
            nh_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == "cs.ca":
            ca_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == "cs.cb":
            cb_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == "cs.c_":
            c_shifts.append(int(txt_list[i][5:])-1)
        elif cs_type == "cs.ha":
            ha_shifts.append(int(txt_list[i][6:])-1)
    return (ca_shifts, cb_shifts, c_shifts, ha_shifts, hn_shifts, nh_shifts)
