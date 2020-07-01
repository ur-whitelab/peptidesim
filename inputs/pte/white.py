import numpy as np
from peptidesim import PeptideSim, utilities
import sys
import os
import sys
import gromacs
from shutil import copyfile, move
import fire


def run(name, seq, peptide_copies, peptide_density, shift_dict, temperature=4 + 273.15, debug=False):
    ps = PeptideSim(name, seq, peptide_copies, job_name='{}'.format(name))
    ps.mdrun_driver = 'gmx'
    ps.run_kwargs = {'nt': 1}
    ps.forcefield = 'amber99sb'
    ps.water = 'tip4p'
    ps.peptide_density = peptide_density  # mg/ml
    ps.ion_concentration = 0.001  # 10mM
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
    ps.mpi_np = None


    ps.run(mdpfile='peptidesim_anneal.mdp', tag='annealing', mdp_kwargs={'nsteps': ns_ts(0.5)})
    ps.run(mdpfile='peptidesim_npt.mdp', tag='equil_npt', mdp_kwargs={'nsteps': ns_ts(2), 'ref_t': temperature})
    pteinfo = ps.pte_replica(cold=temperature, eff_threshold=0.2 if not debug else 0.0, hill_height=1.2, hot=375,
                             max_tries=100 if not debug else 5, mdp_kwargs={'nsteps': ns_ts(0.01)})

    kwargs = [{'ref_t': ti} for ti in pteinfo['temperatures']]

    # set-up camshift
    csinfo = utilities.prepare_cs_data(ps, shift_dict, True)


    # now run a short simulation to confirm everything works
    plumed_file = 'cs2.dat'
    with open(plumed_file, 'w') as f:
        f.write(pteinfo['plumed'])
        f.write(csinfo['plumed'])
    ps.add_file(plumed_file)

    ps.run(mdpfile='peptidesim_nvt.mdp',
           tag='cs-unbiased',
           run_kwargs={'plumed': plumed_file, 'replex': pteinfo['replex']},
           mdp_kwargs=kwargs)

    # run with camshift shifts w EDS
    for kw in kwargs:
        kw['nsteps'] = ns_ts(100)
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


if __name__ == '__main__':
    fire.Fire(run)
