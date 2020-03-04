import numpy as np
import matplotlib.pyplot as plt
from peptidesim import PeptideSim
import textwrap
import sys
import re
import os
import sys
import dill as pickle
from shutil import copyfile, move

import matplotlib
matplotlib.use('Agg')

name = sys.argv[2]
debug = False
pickle_name = name + '.pickle'
MPI_NP = 4
peptide_copies = int(sys.argv[3])
replicas = 4
# try to reload
total_no_atoms = 0  # init
number_chains = 0  # init
ps = 3  # initialize
pte_tune_time = 0.005  # in ns will run at least 4 time, total time will be 40 ns
pt_wte_eds_time_ns = 0.02  # 10 #in ns will run at least 4 times,resulting in 40 ns sim
final_eds_time_ns = 0.05  # 40 in ns will run once resulting in 40ns simulation
replica_eff_threshold = 0.01  # minimum thereshold for both pt-wte and eds+pt-wte
cpt_time = 1  # 60                   #checkpoint file output interval in mins

remd_exchange = 10  # 250
total_no_atoms, number_chains = 0, 0

if(os.path.exists(pickle_name)):
    print('loading restart')
    with open(pickle_name, 'rb') as f:
        ps = pickle.load(f)
        # ps.pickle_name=pickle_name
        print(os.getcwd())
        # ps.rel_dir_name='.'

        def number_all_atoms_chains():
            if (os.path.isdir("{}/data/".format(os.getcwd()))):
                template_file = "{}/data/template.pdb".format(os.getcwd())
                with open(template_file, 'r') as f:
                    lines = f.readlines()
                    lines = lines[-3].strip()
                    lines = lines.split()
                    print(lines)
                    return int(lines[1]), int(lines[5])
            else:
                print("there is no data directory")
        total_no_atoms, number_chains = number_all_atoms_chains()

old_gro_file = ps.gro_file
mdp_kwargs = {'nsteps': int(400 * 5 * 10**2)}
replica_kwargs_0 = [mdp_kwargs.copy() for _ in range(4)]
for k in range(4):
    replica_kwargs_0[k]['ref_t'] = 278

#ps.run(mdp_kwargs=replica_kwargs, run_kwargs={'cpt':cpt_time, 'cpnum':'yes',  'replex':200, 'multi':4,'plumed':'plumed_restraint.dat'},mdpfile='peptidesim_nvt.mdp', tag='restraint',mpi_np=4)
ps.gro_file = old_gro_file

eds_period = 10
remd_exchange_period = 10


def get_replex_e(ps, replica_number):
    with open(ps.sims[-1].location + '/' + ps.sims[-1].metadata['md-log']) as f:
        p1 = re.compile('Repl  average probabilities:')
        p2 = re.compile(
            r'Repl\s*' + ''.join([r'([0-9\.]+\s*)' for _ in range(replica_number - 1)]) + '$')
        ready = False
        # answer=-1
        for line in f:
            if not ready and p1.findall(line):
                ready = True
            elif ready:
                match = p2.match(line)
                if match:
                    return [float(s) for s in match.groups()]


for i in ps.sims:
    print(i.name)
# ps.remove_simulation('pte_tune-0-9852615d')
pte_result = ps.pte_replica(
    mpi_np=MPI_NP,
    max_tries=2,
    min_iters=1,
    mdp_kwargs={
        'nsteps': int(
            pte_tune_time * 5 * 10**5)},
    replicas=replicas,
    hills_file_location=os.getcwd(),
    run_kwargs={
        'cpt': cpt_time},
    hot=279,
    hill_height=1,
    sigma=500,
    bias_factor=20,
    eff_threshold=replica_eff_threshold,
    cold=278,
    exchange_period=2)
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

pte_plumed_script = pte_result['plumed']
replica_temps = pte_result['temperatures']
temps = ','.join(str(e) for e in replica_temps)
kwargs = [{'ref_t': ti} for ti in replica_temps]
print(total_no_atoms, number_chains, eds_period, temps)
# now reload PTE_WTE hills file and run eds with cs2backbone to generate
# the eds parameters with multiple replicas

iteration_number = 0
list_of_indeces = [0]
index_of_eds_file = 0
for any_file in ps.sims:

    if any_file.name.startswith('nvt_conver_eds_2'):

        for file_name in ps.file_list:
            if 'number' in file_name:
                f2 = open(file_name, 'r')
                iteration_number = int(f2.readlines()[0])
            list_of_indeces.append(iteration_number)
    index_of_eds_file = max(list_of_indeces)
    with open('number_{}'.format(index_of_eds_file + 1), 'w') as f:
        f.write('{}'.format(index_of_eds_file + 1))
        f.close()
    ps.add_file('number_{}'.format(index_of_eds_file + 1))
# if the eds is being called for the first time then the eds file index zero,
# hence we do not need call IN_Restart argument in EDS line
# because we not reading any eds bias.

plumed_input0 = textwrap.dedent(
    '''

    peptide: GROUP ATOMS=1-{}
    WHOLEMOLECULES ENTITY0=peptide
    cs: CS2BACKBONE ATOMS=peptide DATADIR=data

    #bias the simalation with EDS and chem chifts until one gets good eds convergence

    PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=200
    #exp_mean: MATHEVAL ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_7,cs.exphn_8,cs.expha_1,cs.expha_2,cs.expha_3,cs.expha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    #calc_mean: MATHEVAL ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO


    eds: EDS ARG=cs.hn-0-5,cs.hn-0-2,cs.hn-0-3,cs.hn-0-4,cs.hn-0-7,cs.hn-0-8,cs.ha-0-4,cs.ha-0-2,cs.ha-0-3,cs.ha-0-6 CENTER_ARG=cs.exphn-0-5,cs.exphn-0-2,cs.exphn-0-3,cs.exphn-0-4,cs.exphn-0-7,cs.exphn-0-8,cs.expha-0-4,cs.expha-0-2,cs.expha-0-3,cs.expha-0-6 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.9 PERIOD={} TEMP=@replicas:{} OUT_RESTART={}/eds_bias_restart_correct0.dat


    ENDPLUMED
    '''.format(total_no_atoms, eds_period, '{' + temps + '}', os.getcwd()))

# if the eds is being called for not the first time then the eds file index is not  zero,
# hence we do need call IN_Restart argument in EDS line because we  reading any eds bias.
# THE index iterates through how many times, it was restarted.
# We also need create a list of eds files that need to be read during restart
eds_file_string = ''
for replica in range(replicas):
    eds_file_string += '{}/eds_bias_restart_correct{}{}.dat,'.format(
        os.getcwd(), index_of_eds_file, replica)
# remove the comma after the last file name
eds_file_string = eds_file_string[:-1]
plumed_input1 = textwrap.dedent(
    '''

    peptide: GROUP ATOMS=1-{}
    WHOLEMOLECULES ENTITY0=peptide
    cs: CS2BACKBONE ATOMS=peptide DATADIR=data

    #bias the simalation with EDS and chem chifts until one gets good eds convergence

    PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=200
    #exp_mean: MATHEVAL ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_7,cs.exphn_8,cs.expha_1,cs.expha_2,cs.expha_3,cs.expha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    #calc_mean: MATHEVAL ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO


    eds: EDS ARG=cs.hn-0-5,cs.hn-0-2,cs.hn-0-3,cs.hn-0-4,cs.hn-0-7,cs.hn-0-8,cs.ha-0-4,cs.ha-0-2,cs.ha-0-3,cs.ha-0-6 CENTER_ARG=cs.exphn-0-5,cs.exphn-0-2,cs.exphn-0-3,cs.exphn-0-4,cs.exphn-0-7,cs.exphn-0-8,cs.expha-0-4,cs.expha-0-2,cs.expha-0-3,cs.expha-0-6 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.9 PERIOD={} TEMP=@replicas:{} OUT_RESTART={}/eds_bias_restart_correct{}.dat IN_RESTART=@replicas:{}


    ENDPLUMED
    '''.format(total_no_atoms, eds_period, '{' + temps + '}', os.getcwd(), index_of_eds_file + 1, '{' + eds_file_string + '}'))
plumed_file = 'plumed_eds_conver_pt_wte' + \
    '{}'.format(index_of_eds_file) + '.dat'
if index_of_eds_file == 0:
    plumed_input = pte_plumed_script + plumed_input0
    with open(plumed_file, 'w') as f:
        f.write(plumed_input)
    ps.add_file(plumed_file)
else:
    plumed_input = pte_plumed_script + plumed_input1
    with open(plumed_file, 'w') as f:
        f.write(plumed_input)
    ps.add_file(plumed_file)
time_ns = .40
replex_eff = 0
for kw in kwargs:
    kw['nsteps'] = int(pt_wte_eds_time_ns * 5 * 10 ** 5)


ps.run(
    mdpfile='peptidesim_nvt.mdp',
    tag='nvt_conver_eds_{}'.format(2),
    mdp_kwargs=kwargs,
    run_kwargs={
        'plumed': plumed_file,
        'replex': remd_exchange,
        'cpt': cpt_time},
    mpi_np=MPI_NP)


with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

eds_conver_ptwte_folder = ps.sims[-1].location
print(eds_conver_ptwte_folder)
colvar_file = '{}/restart_pt_wte.0.dat'.format(
    os.path.abspath(eds_conver_ptwte_folder))

plumed_input = textwrap.dedent(
    '''
    peptide: GROUP ATOMS=1-{}
    WHOLEMOLECULES ENTITY0=peptide
    cs: CS2BACKBONE ATOMS=peptide DATADIR=data

    #bias the simalation with EDS and chem chifts until one gets good eds convergence


eds: EDS ARG=cs.hn-0-5,cs.hn-0-2,cs.hn-0-3,cs.hn-0-4,cs.hn-0-7,cs.hn-0-8,cs.ha-0-4,cs.ha-0-2,cs.ha-0-3,cs.ha-0-6 CENTER_ARG=cs.exphn-0-5,cs.exphn-0-2,cs.exphn-0-3,cs.exphn-0-4,cs.exphn-0-7,cs.exphn-0-8,cs.expha-0-4,cs.expha-0-2,cs.expha-0-3,cs.expha-0-6 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 PERIOD=2 TEMP=@replicas:{} OUT_RESTART={}/{}

    ENDPLUMED'''.format(total_no_atoms, '{' + temps + '}', os.getcwd(), 'OUTPUT_COLVAR_EDS'))


with open('plumed_eds_colvars_new.dat', 'w') as f:
    f.write(plumed_input)
ps.add_file('plumed_eds_colvars_new.dat')
# ps.req_files.append('plumed_eds_colvars_new.dat')
