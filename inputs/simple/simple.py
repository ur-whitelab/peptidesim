from peptidesim import PeptideSim
import sys

seq = sys.argv[1]
name = sys.argv[2]
MPI_NP = 4

ps = PeptideSim(name, [seq], [1], job_name='2mer_{}'.format(name))
ps.peptide_density = 0.005
ps.water = 'spce'
ps.initialize()

ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**5}, mpi_np=MPI_NP)
ps.run(mdpfile='peptidesim_anneal.mdp', tag='anneal_nvt', mpi_np=MPI_NP)
ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod', mdp_kwargs={'nsteps': int(3 * 5*10**5), 'constraints': 'h-bonds'}, mpi_np=MPI_NP)
