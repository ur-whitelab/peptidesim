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
MPI_NP = 8
peptide_copies=3

#try to reload                                                                                
if(os.path.exists(pickle_name)):
    print('loading restart')
    with open(pickle_name, 'rb') as f:
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
ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**2,'rcoulomb':1}, mpi_np=MPI_NP)
ps.run(mdpfile='peptidesim_anneal.mdp',tag='annealing',mdp_kwargs={'nsteps':int(1 * 5*10**2)},mpi_np=MPI_NP, pickle_name=pickle_name )#change the time step to 2 ns
ps.run(mdpfile='peptidesim_npt.mdp', tag='equil_npt', mdp_kwargs={'nsteps': int(1 * 5*10**2)}, mpi_np=MPI_NP, pickle_name=pickle_name)
filename=ps.pte_replica(mpi_np=MPI_NP,final_time=int(1* 5*10**2),replicas=4)
print((ps.sim_dicts), filename)
