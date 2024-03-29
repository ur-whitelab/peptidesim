import numpy as np
from peptidesim import PeptideSim, utilities
import sys
import os
import sys
import gromacs
from shutil import copyfile, move
import fire


def run(name, seq, peptide_copies, peptide_density, shift_dict, eds_reweight=True, debug=False, temperature=4 + 273.15):
    ps = PeptideSim(name, seq, peptide_copies, job_name='{}'.format(name))
    ps.mdrun_driver = 'gmx'
    ps.run_kwargs = {'nt': 4}
    # To use an arbitrary/new/your force field,
    # add the forcefield files to peptidesim list of files,
    # and add it to the peptide_structures folder.
    # This will ensure that the ff files are copied over to every simulation directory.
    # ps.add_file('path/to/forcefield.ff')
    # ps._put_in_dir('peptide_structures')
    ps.forcefield = 'charmm27'
    ps.water = 'tip3p'
    # For your own water model, use water='select' and the input argument for pdb2gmx_args
    ps.peptide_density = peptide_density  # mg/ml
    ps.ion_concentration = 0.008  # 10mM
    ps.initialize()

    def ns_ts(ns, scale=None):
        '''nano seconds to timesteps. Use scale to turn down times for debugging'''
        if scale is None and debug:
            scale = 0.001
        else:
            scale = 1.0
        return str(max(1000, int(ns * 10**3 * 10**3 / 2 * scale)))

    ps.mdrun_driver = 'gmx_mpi'
    ps.run_kwargs = {}
    ps.mpi_np = 1


    ps.run(mdpfile='peptidesim_anneal.mdp', tag='annealing', mdp_kwargs={'nsteps': ns_ts(0.5)})
    ps.run(mdpfile='peptidesim_npt.mdp', tag='equil_npt', mdp_kwargs={'nsteps': ns_ts(2), 'ref_t': temperature})

    ps.mpi_np = None

    pteinfo = ps.pte_replica(cold=temperature, eff_threshold=0.2 if not debug else 0.0,
                             hill_height=1., hot=375, sigma=250, bias_factor=20,
                             max_tries=10 if not debug else 5, mdp_kwargs={'nsteps': ns_ts(0.04)})

    kwargs = [{'ref_t': ti} for ti in pteinfo['temperatures']]

    # set-up camshift
    csinfo = utilities.prepare_cs_data(ps, shift_dict, True)


    # now run a short simulation to confirm everything works
    plumed_file = 'cs2.dat'
    with open(plumed_file, 'w') as f:
        f.write(pteinfo['plumed'])
        f.write(csinfo['plumed'])
    ps.add_file(plumed_file)

    for kw in kwargs:
        kw['nsteps'] = ns_ts(1)

    ps.run(mdpfile='peptidesim_nvt.mdp',
           tag='cs-unbiased',
           run_kwargs={'plumed': plumed_file, 'replex': pteinfo['replex']},
           mdp_kwargs=kwargs)

    # run with camshift shifts w EDS
    for kw in kwargs:
        kw['nsteps'] = ns_ts(40)
    plumed_file = 'eds.dat'
    plumed_file_restart = 'eds-restart.dat'
    with open(plumed_file, 'w') as f, open(plumed_file_restart, 'w') as rf:
        f.write(pteinfo['plumed'])
        f.write(csinfo['plumed'])
        rf.write(pteinfo['plumed'])
        rf.write(csinfo['plumed'])
        eds_str = 'eds: EDS ARG=' + ','.join(csinfo['cs2_names'])
        eds_str += ' CENTER=' + ','.join([str(v) for v in csinfo['cs2_values']])
        eds_str += ' RANGE=' + ','.join(['0.1' for _ in csinfo['cs2_values']])
        eds_str += ' MULTI_PROP=0.2'
        if eds_reweight:
            eds_str += ' LOGWEIGHTS=pte-lw'
        eds_str += ' PERIOD=500 LM TEMP=@replicas:{{{}}}'.format(','.join([str(t) for t in pteinfo['temperatures']]))
        eds_str += ' OUT_RESTART=checkpoint.eds'
        f.write(eds_str + '\n')
        eds_str += ' IN_RESTART=@replicas:{{{}}} RESTART=YES\n'.format(','.join(['checkpoint.{}.eds'.format(i) for i in range(len(pteinfo['temperatures']))]))
        rf.write(eds_str)
    ps.add_file(plumed_file)
    ps.add_file(plumed_file_restart)

    restarted = ps.sims[-1].short_name == 'eds-equil'

    # remove extra headers from EDS file if present
    if restarted:
        for i, md in enumerate(ps.sims[-1].metadata['multi-dirs']):
            with open(os.path.join(ps.sims[-1].location, md, 'checkpoint.{}.eds'.format(i)), 'r') as f:
                restart_lines = f.readlines()
            with open(os.path.join(ps.sims[-1].location, md, 'checkpoint.{}.eds'.format(i)), 'w') as f:
                header_complete = False
                for line in restart_lines:
                    if line[0] != '#':
                        header_complete = True
                        f.write(line)
                    if not header_complete:
                        f.write(line)

    ps.run(mdpfile='peptidesim_nvt.mdp',
           tag='eds-equil',
           run_kwargs={'plumed': plumed_file_restart if restarted else plumed_file, 'replex': pteinfo['replex']},
           mdp_kwargs=kwargs)

    final_dir = os.path.join(ps.sims[-1].location, ps.sims[-1].metadata['multi-dirs'][0])

    plot_couplings(os.path.join(final_dir, 'checkpoint.0.eds'))

if __name__ == '__main__':
    fire.Fire(run)
