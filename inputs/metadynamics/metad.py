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
# TO find the weights for other CV's. Assuming you did not print
# metad.bias previously during simulation, you can compute a new CV and
# measure the weights for those CVs using the following plumed files and
# command line executable called 'plumed driver'
plumed_input_reweight = textwrap.dedent(
    '''
    # after the simulation is done we rememeber to compute these CVs because there was metad applied
    # we need to reweights the biases for those CVs as well.

    gyration: GYRATION ATOMS=1-3
    distance: DISTANCE ATOMS=5,8

    gyration1: GYRATION ATOMS=5-8 # new CV
    distance1: DISTANCE ATOMS=1,3 # new CV

    # restarting metad but we not addinh anymore biases shown by big PACE and zero for hill height
    # Very important: if you are using the plumed driver, the RESTART=YES in metad line
    # also importnt to specify the temperature at the trajectory was generated

    METAD ...
    LABEL=metad TEMP={} RESTART=YES
    ARG=gyration,distance SIGMA=0.05,0.05 HEIGHT=0.0 PACE=10000000 FILE={}/{}/HILLS2
    ... METAD

    # again metad.bias column contains weights, make sure the Stride is equal to 1 if you
    # do not want to lose the hard-earned trajectory points
    # stride=1 means that it will compute the CV and its weights at every frame
    # at the trr file was outputted, aka nstfout in the mdp file that was used
    # to make the trajectory file
    PRINT ARG=metad.bias STRIDE=1 FILE=WEIGHTS
    PRINT ARG=metad.bias,distance,gyration,distance1,gyration1 FILE=COLVAR_OUTPUT_metad_reweight STRIDE=1

    ENDPLUMED''').format(TEMP, os.getcwd(), ps.sims[-1].location)

with open('plumed_metad_reweight.dat', 'w') as f:
    f.write(plumed_input_reweight)
ps.add_file('plumed_metad_reweight.dat')
# command line excutable used in reweighting new CV's after the simulation
# is completed
driver = 'plumed driver'

# the following command is run the command line
# timestep is the time step in picesecond used during generation of the trajectory file. In peptidesim the default is 2 fs.
# the trajectory stride is equal to the nstfxout in mdp file as mentioned
# before.


pwd = os.getcwd()
os.chdir(ps.sims[-1].location)
# we are using method mentioned in https://www.plumed.org/doc-master/user-doc/html/lugano-3.html
# Exercise 4: reweighting

result = subprocess.call(
    '{} --plumed {}/{} --mf_trr traj.trr --timestep 0.002 --trajectory-stride {}'.format(
        driver, pwd, 'plumed_metad_reweight.dat', STRIDE), shell=True)
print(result)
