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
import subprocess
import matplotlib
matplotlib.use('Agg')
seq = 'AAA'
name = 'metad1'
debug = False
pickle_name = name + '.pickle'
MPI_NP = 4
peptide_copies = 1


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

ps.run(
    mdpfile='peptidesim_emin.mdp',
    tag='init_emin',
    mdp_kwargs={
        'nsteps': 8 *
        10**2,
        'rcoulomb': 1},
    mpi_np=MPI_NP)

# print ps.pickle_name, 'picklename2'

ps.run(
    mdpfile='peptidesim_anneal.mdp',
    tag='annealing',
    mdp_kwargs={
        'nsteps': int(
            0.1 * 5 * 10**2)},
    run_kwargs={
        'cpt': 5},
    mpi_np=MPI_NP)


ps.run(
    mdpfile='peptidesim_npt.mdp',
    tag='equil_npt',
    mdp_kwargs={
        'nsteps': int(
            0.1 * 5 * 10**2),
        'ref_t': 278},
    run_kwargs={
        'cpt': 5},
    mpi_np=MPI_NP)
TEMP = 278
STRIDE = 100  # the number frames at which the trajecotry is written and biases are added and then are reweighted

plumed_input = textwrap.dedent(
    '''
    gyration: GYRATION ATOMS=1-3
    distance: DISTANCE ATOMS=5,8
    METAD ...
    LABEL=metad TEMP={}
    ARG=gyration,distance SIGMA=0.05,0.05 HEIGHT=0.3 PACE={} FILE=HILLS2
    ... METAD

    PRINT ARG=metad.bias,distance,gyration FILE=COLVAR_OUTPUT_metad STRIDE={}
    ENDPLUMED''').format(TEMP, STRIDE, STRIDE)

# IMPORTANT the file COLVAR_OUTPUT_Metad has the biad being added under
# metad.bias column for each CV that was  used a CV to be biased. the
# column can be used to measure the bias for other CV that need to be
# measured after the simulation is over.
with open('plumed_metad.dat', 'w') as f:
    f.write(plumed_input)
ps.add_file('plumed_metad.dat')

ps.run(
    mdpfile='peptidesim_nvt.mdp',
    tag='meat2',
    mdp_kwargs={
        'nsteps': int(
            100 * 5 * 10**2),
        'nstxout': int(STRIDE),
        'ref_t': TEMP
    },
    run_kwargs={
        'cpt': 5, 'plumed': 'plumed_metad.dat'},
    mpi_np=MPI_NP)
