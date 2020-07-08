import numpy as np
from peptidesim import PeptideSim, utilities
import sys
import os
import sys
import gromacs
from shutil import copyfile, move
import fire


def run(dir_name, seqs):
    ps = PeptideSim(dir_name, seqs)
    print(ps.sims_dict.keys())

    eds_sim = ps.get_simulation('eds-equil')
    final_dir = os.path.join(eds_sim.location, eds_sim.metadata['multi-dirs'][0])

    utilities.plot_couplings(os.path.join(final_dir, 'checkpoint.0.eds'))
    utilities.plot_wte_overlap(ps, 'eds-equil')
    cs_sim = ps.get_simulation('cs-unbiased')
    utilities.plot_cs([cs_sim, eds_sim])


if __name__ == '__main__':
    fire.Fire(run)
