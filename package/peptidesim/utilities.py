import matplotlib.pyplot as plt
import numpy as np
import signal
import functools
import os
import pandas as pd
import pkg_resources
import subprocess
from string import ascii_uppercase


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


def plot_couplings(eds_filename, output_plot='couplings.png', weights_file=None):
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
            if index == len(cv_names):
                break
            cv = cv_names[index]
            ax[i, j].plot(data.time, data[f'{cv}_coupling'])
            ax[i, j].set_title(cv)
            index += 1
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)


def plot_pte_overlap(ps, production_sim_name, output_plot='wte_overlap.png'):
    '''
    This will show how much overlap there is between the PTE bias and the production simulation
    '''
    # sum bias from pte step
    hills_file_loc = os.path.join(ps.dir_name, 'pte_hills', 'HILLS_PTE.0')
    bias_file_loc = os.path.join(ps.dir_name, 'pte_hills', 'BIAS.0')
    if not os.path.exists(hills_file_loc):
        raise FileNotFoundError('Could not find PT-WTE hills file')
    result = subprocess.call(f'plumed sum_hills --bin 1000 '
                             f'--hills {hills_file_loc} --outfile {bias_file_loc}',
                             shell=True)
    if result != 0:
        raise ValueError('Failed to run sum_hills')

    # load bias
    bias = np.loadtxt(bias_file_loc)

    # now load details from simulation
    s = ps.get_simulation(production_sim_name)
    weights_file = os.path.join(s.location, s.metadata['multi-dirs'][0], 'weights.0.dat')
    if not os.path.exists(weights_file):
        raise FileNotFoundError('Could not find weights file')

    data = np.loadtxt(weights_file)

    # now plot
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].set_title('Bias')
    ax[0].plot(bias[:, 0], bias[:, 1])
    ax[1].set_title('PE')
    hist, bin_edges = np.histogram(data[:, 2], bins=250)
    ax[1].plot((bin_edges[1:] + bin_edges[:-1]) / 2, hist)
    plt.tight_layout()
    plt.savefig(output_plot)

    plt.figure()
    plt.title('weights')
    plt.plot(data[:, 0], data[:, 1])
    plt.savefig('weights.png')




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

def pretty_cslabel(cs_name):
    chain, resid, name, atom = cs_name.split('-')
    atom_labels = {'CB': r'C$\beta$', 'C': r'C', 
                   'CA': r'C$\alpha$', 'H': r'H',
                   'HA': r'H$\alpha$'}
    if atom in atom_labels:
        atom = atom_labels[atom]
    
    return 'Chain ' +  ascii_uppercase[int(chain)] + ' ' + name + resid + ' ' + atom


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


def load_pte_weights(sim):
    weights_filename = sim.get_file('weights.dat')
    return np.loadtxt(weights_filename)

def plot_cs(sim_list, use_weights=True, output_plot='cs.png'):
    '''Will plot chemical shift and its average over time across multiple simulations. Can also reweight
    '''
    from tqdm import tqdm

    if type(sim_list) != list:
        sim_list = [sim_list]
    sim_data = []
    weights = []
    for sim in sim_list:
        cs_filename = sim.get_file('cs_shifts.dat')
        with open(cs_filename, 'r') as f:
            header = f.readline().split()[2:]
        # get target values using kind of hack
        targets = dict()
        for i, h in enumerate(header):
            if 'exp' in h:
                n = header[i - 1].split('all-avg-')[-1]
                targets[n] = h
        data = pd.read_table(cs_filename, sep=r'\s+', comment='#', names=header)
        # drop zeroth row
        data = data.drop(data.index[0])
        sim_data.append(data)
        # load weights
        if use_weights:
            wdata = load_pte_weights(sim)
            # drop zeroth row
            weights.append(wdata[1:, 1])
        else:
            weights.append(np.ones(data.count()))

    # get shift names
    shift_names = set()
    for n in sim_data[-1].columns:
        if 'all-avg' in n:
            sn = n.split('all-avg-')[-1]
            shift_names.add(sn)

    fig, ax = plt.subplots(nrows=len(shift_names) // 2, ncols=2, figsize=(6, 4), sharex=True)
    index = 0
    shift_names = list(shift_names)
    shift_names.sort()
    for i in range(len(shift_names) // 2):
        for j in range(2):
            if index == len(shift_names):
                break
            s = shift_names[index]
            last_time = 0
            for k, data in enumerate(sim_data):
                if k > 0:
                    ax[i, j].axvline(last_time, linestyle='--', color='gray')
                ax[i, j].plot(data.time / 1000 + last_time, data[f'avg-{s}'], color='C0', alpha=0.5, label='Instantaneous' if k == 0 else None)
                # compute weighted running average
                widx = 0  # data.time.count() // 2
                time = data.time[widx + 1:] / 1000
                wmean = []
                shifts = data[f'avg-{s}'].to_numpy()
                wmax = np.max(weights[k])
                w = np.exp(weights[k] - wmax)
                for wi in tqdm(range(widx + 1, data.time.count())):
                    wmean.append(np.sum(shifts[widx:wi] * w[widx:wi]) / np.sum(w[widx:wi]))
                ax[i, j].plot(time + last_time, wmean, color='C0', alpha=1.0, label='Run. Avg.' if k == 0 else None)
                if use_weights:
                    ax[i, j].plot(data.time[1:] / 1000 + last_time,
                                  [np.mean(shifts[:i]) for i in range(1, data.time.count())],
                                  color='C2', alpha=1.0, label='Unw. Run. Avg.' if k == 0 else None, linestyle='--')
                ax[i, j].plot(data.time / 1000 + last_time, data[targets[s]], color='C1', alpha=1.0, label='Experiment' if k == 0 else None)
                last_time = data.time.max() / 1000
            ax[i, j].set_title(pretty_cslabel(s))
            ax[i, j].set_xlabel('Time [ns]')
            ax[i, j].set_ylabel('Shift [ppm]')
            if i == 0 and j == 1:
                ax[i,j].legend(loc='upper left', bbox_to_anchor=(1.05, 1))
            index += 1
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)

def plot_pte_ramachandran(ps, sim, temperature, pte_stride=250, traj_stride=1000, units='kcal/mol', output_name='ramachandran_data'):
    
    import matplotlib as mpl

    mpl.rc('text', usetex=True)
    

    # line-up weights with trajectory steps
    lcm = np.lcm(pte_stride, traj_stride)
    pte_every = lcm // pte_stride
    traj_every = lcm // traj_stride
    print(f'Based on strides, will only compute Ramachandran at steps {lcm}')

    # set energy
    plumed_script = f'UNITS ENERGY={units}\n'

    # make pdb for plumed
    data_dir = os.path.join(ps.dir_name, output_name)
    os.makedirs(data_dir, exist_ok=True)
    template_pdb = os.path.join(data_dir, 'template.pdb')
    atom_number = pdb_for_plumed(input_file=ps.pdb_file, output_file=template_pdb)
    plumed_script += f'MOLINFO STRUCTURE={template_pdb}\n'

    # get PTE info
    weights_filename = sim.get_file('weights.dat')
    # I'm pretty sure I did this correctly, but I cannot get plumed to believe me
    plumed_script += f'pte-lw: READ FILE={weights_filename} VALUES=pte-lw STRIDE={traj_every} EVERY={pte_every} IGNORE_TIME\n'
    
    # make plumed script
    names =[]
    for i, s in enumerate(ps.sequences):
        index = 0
        for j in range(ps.counts[i]):
            for k in range(len(s)):
                if k == 0 or k == len(s) - 1:
                    continue
                n = f'{i}-{k+1}-{s[k]}'
                names.append(n)
                output_grid = os.path.join(data_dir, f'fes-{n}.dat')
                plumed_script += f'phi-{n}: TORSION ATOMS=@phi-{ascii_uppercase[index]}{k+1}\n'\
                                  f'psi-{n}: TORSION ATOMS=@psi-{ascii_uppercase[index]}{k+1}\n'\
                                  f'HISTOGRAM ...\n'\
                                  f'ARG=phi-{n},psi-{n}\n'\
                                  f'GRID_MIN=-3.14,-3.14\n'\
                                  f'GRID_MAX=3.14,3.14\n'\
                                  f'GRID_BIN=200,200\n'\
                                  f'BANDWIDTH=0.05,0.05\n'\
                                  f'LABEL=hrr-{n}\n'\
                                  f'STRIDE={traj_every}\n'\
                                  f'... HISTOGRAM\n'\
                                  f'fes-{n}: '\
                                  f'CONVERT_TO_FES GRID=hrr-{n} TEMP={temperature}\n'\
                                  f'DUMPGRID GRID=fes-{n} FILE={output_grid}\n'
            index += 1

    plumed_file = os.path.join(data_dir, 'compute_ramachandran.dat')
    with open(plumed_file, 'w') as f:
        f.write(plumed_script)
    if 'multi-dirs' in sim.metadata:
        trj = os.path.join(sim.location, 
                            sim.metadata['multi-dirs'][0], 
                            sim.metadata['traj'])
    else:
        path = os.path.join(sim.location, sim.metadata['traj'])

    if True:
       result = subprocess.call(f'plumed driver --mf_trr {trj} '
                                f'--plumed {plumed_file}',
                                shell=True)
       if result != 0:
           raise ValueError('Failed to run plumed_driver')

    # now we plot
    for n in names:
        data  = np.loadtxt(os.path.join(data_dir, f'fes-{n}.dat'))
        imdata = data[:,2].reshape(200,200)
        plt.figure(figsize=(4,3))
        #plt.title(n)
        # pick reasonable scaling
        result = np.quantile(imdata[imdata != np.inf], [0.01, 0.99])
        # set background color
        plt.gca().set_facecolor(mpl.cm.get_cmap('viridis')(0))
        # plot
        plt.imshow(imdata - result[1], cmap='viridis_r', origin='lower', vmin=result[0] - result[1], vmax=0)
        plt.xticks((50, 100, 150), (r'$-\frac{\pi}{2}$', '0', r'$\frac{\pi}{2}$'))
        plt.yticks((50, 100, 150), (r'$-\frac{\pi}{2}$', '0', r'$\frac{\pi}{2}$'))
        plt.xlabel('$\phi$ [Rad]')
        plt.ylabel('$\psi$ [Rad]')
        cbar = plt.colorbar()
        cbar.ax.set_ylim(0, result[0] - result[1])
        cbar.ax.set_ylabel(f'$\Delta$ A[{units}]')
        plt.tight_layout()
        plt.savefig(n + '.png', dpi=300)
            
