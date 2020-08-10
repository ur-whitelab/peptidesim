'''Simulate a peptide with a defined sequence and conditions.

Example
-------
Here's an example for creating a ``PeptideSim`` object with the peptide AEAE using the default configuration, saving in the current directory. ::

    p = PeptideSim( dir_name = '.', seqs = ['AEAE'], counts = [1] )

Here's an example showing **one** AEAE peptide and **two** LGLG peptides, saving in the current directory. ::

    #counts in order of the list of peptides
    p = PeptideSim( dir_name = '.', seqs = ['AEAE', 'LGLG'], counts = [1,2])
'''


from functools import reduce
import numpy as np
import logging
import os
import shutil
import datetime
import subprocess
import re
import textwrap
import sys
import pkg_resources
import contextlib
import uuid
import json
import ast
import requests
import signal
import PeptideBuilder
import Bio.PDB
import glob
import dill
from math import *
from .utilities import *
from .version import __version__
from traitlets.config import Configurable, Application, PyFileConfigLoader
from traitlets import Int, Float, Unicode, Bool, List, Instance, Dict
import gromacs
import gromacs.cbook
gromacs.environment.flags['capture_output'] = True


class SimulationInfo(object):
    '''A class that stores information about a simulation. Only for keeping history of simulation runs
    '''

    def __init__(self, name, short_name):
        self.run_fxn = None
        self.run_kwargs = None
        self.name = ''
        self.location = ''
        self.short_name = ''
        self.restart_count = 0
        self.complete = False
        self.metadata = dict()
        self.name = name
        self.short_name = short_name

    def __str__(self):
        metadata_str = json.dumps(self.metadata, indent=12, sort_keys=True, default=str)
        return textwrap.dedent(f'''
        name          : {self.name}
        location      : {self.location}
        restart_count : {self.restart_count}
        run_fxn       : {self.run_fxn}
        run_kwargs    : {self.run_kwargs}
        metadata      : {metadata_str}''')[1:]

    def run(self, run_fxn=None, run_kwargs=None):
        if(self.restart_count > 0):
            if (run_fxn is not None and run_fxn != self.run_fxn) or (run_kwargs is not None and run_kwargs != self.run_kwargs):
                raise ValueError(
                    'Name collision in simulation. You tried to repeat a non-identical simulation')
        else:
            self.run_fxn = run_fxn
            self.run_kwargs = run_kwargs

        # have to add now in case we die
        self.restart_count += 1
        result = self.run_fxn(**self.run_kwargs)

        self.complete = True
        return result

    def get_file(self, filename, replica_index=0):
        if 'multi-dirs' in self.metadata:
            path = os.path.join(self.location,
                                self.metadata['multi-dirs'][replica_index],
                                '{0:}.{2:}.{1:}'.format(*filename.split('.'), replica_index))
        else:
            path = os.path.join(self.location, filename)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f'Could not find file {path} in simulation {self.name}')
        return path



class PeptideSim(Configurable):
    '''Class to use for conducting simulations.

    '''

    sim_name = Unicode('peptidesim',
                       help='The name for the type of simulation job (e.g., NVE-equil-NVT-prod)',
                       )

    config_file = Unicode('peptidesim_config.py',
                          help='The config file to load',
                          )

    req_files = List(
        help='List of files required for simulation. For example, restraints or plumed input'
    )

    pressure = Float(0,
                     help='Barostat pressure. Ignored if not doing NPT'
                     )

    peptide_density = Float(20,
                            help='The density of the peptides in milligrams / milliliter',
                            )

    ion_concentration = Float(0.002,
                              help='The concentration of sodium chloride to add in moles / liter'
                              )

    log_file = Unicode('simulation.log',
                       help='The location of the log file. If relative path, it will be in simulation directory.',
                       )
    packmol_exe = Unicode('packmol',
                          help='The command to run the packmol program.'
                          )
    demux_exe = Unicode('demux',
                        help='The command to demux the replica temperatures.'
                        )

    forcefield = Unicode('charmm27',
                         help='The gromacs syntax forcefield',
                         )
    water = Unicode('tip3p',
                    help='The water model to use',
                    )

    pdb2gmx_args = Dict(dict(ignh=True),
                        help='Any additional special arguments to give to pdb2gmx, aside from force-field and water which are separately specified.',
                        )

    mdp_directory = Unicode('.',
                            help='The directory to find gromacs MDP files'
                            )
    mdp_base = Unicode('peptidesim_base.mdp',
                       help='The MDP file containing basic forcefield parameters'
                       )

    mdp_emin = Unicode('peptidesim_emin.mdp',
                       help='The energy miniziation MDP file. Built from mdp_base. Used specifically for adding ions'
                       )

    mpiexec = Unicode('mpiexec',
                      help='The MPI executable'
                      )

    run_kwargs = Dict(dict(),
                      help='mdrun kwargs'
                      ).tag(config=False)

    mpi_np = Int(1,
                 allow_none=True,
                 help='Number of mpi processes. If not set, it is not specified'
                 )

    mdrun_driver = Unicode(None,
                           allow_none=True,
                           help='An override command for mdrun. Replaces gromacswrapper.cfg prefix (e.g., gmx)'
                           )

    @property
    def peptidesim_version(self):
        ''' the version of peptidesim used when object was instantiated
        '''
        return self._peptidesim_version
    @property
    def box_size_angstrom(self):
        ''' a property that returns a box size in Angstroms
        '''
        return self._box_size

    @box_size_angstrom.setter
    def box_size_angstrom(self, v):
        assert(len(v) == 3)
        self._box_size[:] = v[:]

    @property
    def box_size_nm(self):
        '''  a property that returns a box size in nm
        '''
        return [x / 10 for x in self._box_size]

    @box_size_nm.setter
    def box_size_nm(self, v):
        assert(len(v) == 3)
        self._box_size = [x * 10 for x in v]

    @property
    def file_list(self):
        ''' a property that returns the list of files needed for a simulations
        '''
        result = []
        result.extend(self._file_list)
        # for now, we'll keep these as paths relative to basedirectory. So we won't have local copies everywhere
        # result.extend([self.pdb_file, self.gro_file, self.top_file, self.tpr_file])
        return result

    @property
    def pdb_file(self):
        ''' a property that returns the latest generated pdb file
        '''
        if len(self._pdb) > 0:
            return os.path.normpath(os.path.join(self._rel_dir_name, self._pdb[-1]))
        elif self.gro_file is not None:
            output = '{}.pdb'.format(os.path.basename(
                self.gro_file).split('.gro')[0])
            gromacs.editconf(f=self.gro_file, o=output)
            return os.path.normpath(os.path.join(self._rel_dir_name, output))
        return None

    @pdb_file.setter
    def pdb_file(self, f):
        self._pdb.append(self._convert_path(f))

    @property
    def gro_file_list(self):
        ''' The latest gro file(s) as a list
        '''
        if len(self._gro) == 0:
            return None
        if type(self._gro[-1]) is list:
            result = self._gro[-1]
        else:
            result = [self._gro[-1]]
        return [os.path.normpath(os.path.join(self._rel_dir_name, ri)) for ri in result]

    @property
    def gro_file(self):
        '''The latest gro file
        '''
        if len(self._gro) == 0:
            return None
        if type(self._gro[-1]) is list:
            result = self._gro[-1][0]
        else:
            result = self._gro[-1]
        return os.path.normpath(os.path.join(self._rel_dir_name, result))

    @gro_file.setter
    def gro_file(self, f):
        if type(f) is list:
            f = [self._convert_path(fi) for fi in f]
        else:
            f = self._convert_path(f)
        self._gro.append(f)

    @property
    def top_file(self):
        '''a property that returns the latest topology file
        '''
        if len(self._top) == 0:
            return None
        return os.path.normpath(os.path.join(self._rel_dir_name, self._top[-1]))

    @top_file.setter
    def top_file(self, f):
       self._top.append(self._convert_path(f))

    @property
    def tpr_file(self):
        ''' a property that returns the latest tpr file
        '''
        if len(self._tpr) == 0:
            return None
        return os.path.normpath(os.path.join(self._rel_dir_name, self._tpr[-1]))

    @tpr_file.setter
    def tpr_file(self, f):
        self._tpr.append(self._convert_path(f))

    @property
    def traj_file(self):
        if len(self._trr) == 0:
            return None
        return os.path.normpath(os.path.join(self._rel_dir_name, self._trr[-1]))

    @traj_file.setter
    def traj_file(self, f):
        ''' a property that returns the latest trr file
        '''

        self._trr.append(self._convert_path(f))

    @property
    def ndx(self):
        '''a property that returns the latest ndx file
        '''
        if self._ndx_file is None:
            return None
        n = gromacs.fileformats.NDX()
        n.read(os.path.normpath(os.path.join(self._rel_dir_name, self._ndx_file)))
        return n

    @ndx.setter
    def ndx(self, n):
        self._ndx_file = self._convert_path(n)
        self._file_list.append(self._ndx_file)

    @property
    def sims(self):
        '''a property that returns a list of simulation objects'''
        return self._sim_list

    @property
    def sims_dict(self):
        '''a property that returns a dictionary of simulation names'''
        return self._sims

    def __new__(cls, dir_name=None, seqs=None, counts=None, config_file=None, job_name=None, *args, **kwargs):
        if dir_name is not None:
            if job_name is None:
                _job_name = os.path.split(dir_name)[-1]
            else:
                _job_name = job_name
            pickle_file = os.path.join(
                dir_name, '{}-restart.pickle'.format(_job_name))
            if os.path.exists(pickle_file):
                with open(pickle_file, 'r+b') as f:
                    inst = dill.load(f)
                if not isinstance(inst, cls):
                    raise TypeError(
                        'Unpickled object is not of type {}'.format(cls))
                inst.log.info(
                    'Successfully loaded from restart pickle {}'.format(pickle_file))
                inst.log.info('Current state: \n{}'.format(
                    textwrap.indent(str(inst), '------   ')))
                inst._pickle_loaded = True
                return inst
        inst = super(PeptideSim, cls).__new__(
            cls, dir_name, seqs, counts, config_file, job_name, *args, **kwargs)
        inst._pickle_loaded = False
        return inst

    def __init__(self, dir_name, seqs, counts=None, config_file=None, job_name=None):
        '''This is an initiator that takes the arguments from the command
        line and creates the class simulation.

        Parameters
        ----------
        dir_name : str
            name of the directory where your simulation should be saved. Can have multiple levels. The last directory on path will be taken as the name of the simulation for logging purposes if not overridden by job_name
        seqs : List[str]
            A list of amino acid sequences.
        counts : List[int]
            A list of the number of occurrences of each amino acid, in order.
        '''

        # detect if loaded from pickle file
        if self._pickle_loaded:
            return

        # We have to declare class attributes at the instance level
        # for pickling purposes. Silly, I know.

        # Keep a chain of all files created. Hide behind properties
        self._top = []
        self._gro = []
        self._pdb = []
        self._tpr = []
        self._rmsd = []
        self._trr = []
        self._file_list = []
        self._ndx_file = None
        self._peptidesim_version = __version__
        self._initialized = False

        # set-up saving
        self._save_count = 0
        self._save_directory = os.path.join(os.path.abspath(dir_name), 'restarts')
        os.makedirs(self._save_directory, exist_ok=True)

        # keep track of simulations in case of restart needs
        self._sims = dict()
        self._sim_list = []

        # other variables
        self._box_size = [0, 0, 0]  # box size in angstroms

        # Set-up directory to begin with
        self.job_name = job_name
        if job_name is None:
            self.job_name = os.path.split(dir_name)[-1]

        self.dir_name = dir_name
        self._rel_dir_name = '.'

        if not os.path.exists(self.dir_name):
            os.mkdir(self.dir_name)

        for f in self.req_files:
            if os.path.isdir(f):
                shutil.copytree(f, self._convert_path(f))
            else:
                shutil.copy2(f, self._convert_path(f))
            self._file_list.append(os.path.basename(f))

        self._start_logging()

        # Note: use load_pyconfig_files to merge them. Useful in future
        # load in configuration file
        config = config_file
        if config_file is not None:
            self.log.info('Loading config file {}'.format(config_file))
            config = PyFileConfigLoader(config_file).load_config()
        else:
            self.log.info('Using default configuration'.format(config_file))

        self.log.info('Loaded {}:'.format(str(config)))
        super(PeptideSim, self).__init__(config=config, parent=None)

        # store passed parameters
        self.sequences = seqs
        # store the copies we'd like to have of each sequence
        if counts is None:
            counts = [1 for s in seqs]
        self.counts = counts
        self.log.info(
            'Have {} many of these sequences {}:'.format(counts, seqs))

    def _start_logging(self):
        # check if logger is relative
        # split path and see if folder is empty
        if self.log_file == os.path.basename(self.log_file):
            self.log_file = os.path.join(self.dir_name, self.log_file)

        # don't know how we got here, so we'll just add our logger
        file_handler = logging.FileHandler(self.log_file)
        formatter = logging.Formatter(
            '%(asctime)s [%(filename)s, %(lineno)d, %(funcName)s]: %(message)s (%(levelname)s)')
        file_handler.setFormatter(formatter)

        self.log = logging.getLogger('peptidesim:{}'.format(self.job_name))
        self.log.addHandler(file_handler)

        self.log_handler = file_handler
        self.log.setLevel(logging.DEBUG)
        self.log.info('Started logging for PeptideSim...')

    def __str__(self):
        data = {}
        for k, v in self.traits().items():
            data[k] = v.get(self)
        for k, v in self.__dict__.items():
            if k not in data and type(v) in [str, int, float, list, dict, tuple, str] and k[0] != '_':
                data[k] = v
        # list of properties I think are important
        data['box_size_nm'] = self.box_size_nm
        data['file_list'] = self.file_list
        data['current_tpr'] = self.tpr_file
        data['current_gro'] = self.gro_file
        data['current_traj'] = self.traj_file
        data['current_top'] = self.top_file
        data['current_pdb'] = self.pdb_file
        data['initialized'] = self._initialized
        data['peptidesim_version'] = self.peptidesim_version
        keys = list(data.keys())
        keys.sort()
        result = ['==========PeptideSim Object==========']
        for k in keys:
            if k == '_sim_list':
                continue
            v = data[k]
            result.append('{:14s} : {}'.format(k, v))
        result.append('------------Simulations--------------')
        for i, s in enumerate(self._sim_list):
            result.append('Simulation {}'.format(i))
            result.append(textwrap.indent(str(s), '  '))
        result.append('=====================================')
        return '\n'.join(result)

    def initialize(self, force=False, acetylate=False, amidate=False):
        '''Build PDB files, pack them, convert to gmx, add water and ions

        This method accomplishes the following steps:
          1. Use Bio Python to convert sequences into PDB files
          2. Combine the PDB files using packmol
          3. Convert them into gmx files using the pdb2gmx command and configuration parameters
          4. Add water using the editconf/genbox
          5. Add ions
          6. Compile our energy-minimization tpr file for purposes of adding ions
        '''
        # generate pdbs from sequences and store their extents
        if not force and self._initialized:
            self.log.warn('Ignoring request to re-initialize. Pass force=True to initialize to force')
            return
        self._structure_extents = []
        self.peptide_mass = []
        self.peptide_pdb_files = []
        for i, s in enumerate(self.sequences):
            structure, minmax, mass = self._pdb_file_generator(
                s, 'seq_' + str(i), amidate=amidate, acetylate=acetylate)
            self.peptide_pdb_files.append(structure)
            self._structure_extents.append(minmax)
            self.peptide_mass.append(mass)

        # pack the peptides together into an initial structure
        self._packmol()

        # now get gromcas files
        self._pdb2gmx(amidate=amidate, acetylate=acetylate)

        # center the peptides
        self._center()

        # add solvent
        self._solvate()

        # add ions
        self._add_ions()

        # energy minimize it
        self.run(mdpfile='peptidesim_emin.mdp', tag='initialize-emin', mdp_kwargs={'nsteps': 500})

        self._initialized = True
        self.save('initialized')
        self.log.info('Completed Initialization')

    def demux(self, path_to_replica_dir):
        with self._put_in_dir('demux'):

            class Demux(gromacs.core.Command):
                command_name = self.demux_exe
            cmd = Demux()
            result = cmd('{}/md0.log'.format(path_to_replica_dir))
            if result[0] != 0:
                self.log.error(
                    'Demux failed with retcode {}. Out: {} Err: {} '.format(*result))
            else:
                self.log.info(
                    'Demux succeeded with retcode {}'.format(*result))

                assert os.path.exists(
                    'replica_temp.xvg'), 'Demux succeeded with retcode {} but has no output. Out: {} Err: {}'.format(*result)
        current_dir = os.getcwd()
        return '{}/replica_temp.xvg'.format(current_dir)

    def pte_replica(self, tag='pte_tune', mpi_np=None, replicas=8, max_tries=30, min_iters=4,
                    mdp_kwargs=None, run_kwargs=None, force=False,
                    cold=300.0, hot=400.0, eff_threshold=0.3,
                    hill_height=1.2, sigma=140.0, bias_factor=10,
                    exchange_period=25, dump_signal=None):
        '''
        Runs NVT tuning with plumed PTE and replica exchange to obtain high replica efficiency.

        Parameters
        ----------
        tag : str
            name of run
        mpi_np: int
            number of mpi processes. Must be multiple of replicas
        replicas : int
            Number of replicas
        min_iters: int
            minimum number of simulations required to run
        max_tries : int
            number of simulations to attempt for achieving eff_threshold
        mdp_kwargs : dict
            additional gromacs mdp arguments
        hot : float
            hottest replica temperatrue
        force : boolean
            Accept final replica exchange efficiency, no matter the value
        cold : float
            coldest replica temperatrue
        eff_threshold : float
            The target replica exchange efficiency
        hill_height : float
             hill height for pt-wte
        sigma: float
            the width of the gaussian for pte-wte
        bias_factor: int or float
            biasfactor to add to pte_wte
        exchange_period : int
            the number of steps between exchange attempts

        Returns
        ----------
        A dictionary with:
        plumed : str
            A plumed input section that will utilize the PTE bias
        efficiency : float
            The final efficiency
        temperature : list
            The list of temperatures used for the biases
        '''

        if mdp_kwargs is None:
            mdp_kwargs = {}
        if run_kwargs is None:
            run_kwargs = {}
        def get_replex_e(self, replica_number):
            pte_sims = []
            index = 0
            for simulation in self.sims:
                if simulation.name.startswith(tag):
                    pte_sims.append(index)
                index = index+1
            with open(self.sims[pte_sims[-1]].location + '/' + self.sims[pte_sims[-1]].metadata['md-log']) as f:
                p1 = re.compile('Repl  average probabilities:')
                p2 = re.compile(
                    r'Repl\s*' + ''.join([r'([0-9\.]+\s*)' for _ in range(replica_number - 1)]) + '$')
                ready = False
                for line in f:
                    if not ready and p1.findall(line):
                        ready = True
                    elif ready:
                        match = p2.match(line)
                        if match:
                            return[float(s) for s in match.groups()]
            raise RuntimeError(
                'Unable to parse replica exchange efficiency. Probably simulation was incomplete')

        if min_iters >= max_tries:
            raise ValueError(
                'minimum number of simulations should be greater than the maximum number of iterations')
        hills_file_location = os.path.join(os.path.abspath(self.dir_name), 'pte_hills')
        if not os.path.exists(hills_file_location):
            os.mkdir(hills_file_location)

        # replica temperatures
        for i in range(max_tries):

            #replex eff initiated
            replex_eff = 0
            replica_temps = [
                cold * (hot / cold) ** (float(m) / (replicas - 1)) for m in range(replicas)]
            #plumed input for WT-PTE
            temps = ','.join(str(e) for e in replica_temps)
            plumed_input_name = self._convert_path(os.path.join(hills_file_location, 'plumed_wte.dat'))
            plumed_input = textwrap.dedent(
                '''
                RESTART
                ene: ENERGY
                METAD ...
                LABEL=METADPT
                ARG=ene
                SIGMA={}
                HEIGHT={}
                PACE=250
                TEMP=@replicas:{{{}}}
                FILE={}/HILLS_PTE
                BIASFACTOR={}
                ... METAD
                PRINT ARG=ene STRIDE=250 FILE={}/COLVAR_PTWTE
                '''.format(sigma, hill_height, temps, hills_file_location, bias_factor, hills_file_location))
            # putting the above text into a file
            with open(plumed_input_name, 'w') as f:
                f.write(plumed_input)

            self.add_file(plumed_input_name)

            # arguments for WT-PTE
            replica_kwargs = [mdp_kwargs.copy() for _ in range(replicas)]
            for k, ti in enumerate(replica_temps):
                replica_kwargs[k]['ref_t'] = ti

            plumed_output_script = None
            run_kwargs = {}
            run_kwargs['plumed'] = 'plumed_wte.dat'
            run_kwargs['replex'] = exchange_period

            self.run(tag=tag + '-{}'.format(i), mdpfile='peptidesim_nvt.mdp',
                     mdp_kwargs=replica_kwargs, mpi_np=mpi_np,
                     run_kwargs=run_kwargs,
                     dump_signal=dump_signal)

            replex_eff = min(get_replex_e(self, replicas))
            if (replex_eff >= eff_threshold and i >= (min_iters-1)) or (i == (max_tries - 1) and force):
                self.log.info('Completed the simulation. Reached replica exchange efficiency of {}.'
                              ' The replica temperatures were {}. The name of the plumed input scripts is'
                              '\'plumed_wte.dat\'. Continuing to production'.format(replex_eff, replica_temps))
                plumed_output_script = textwrap.dedent(
                    '''
                 ene: ENERGY
                 METAD ...
                 RESTART=YES
                 LABEL=METADPT
                 ARG=ene
                 SIGMA={}
                 HEIGHT={}
                 PACE=99999999
                 TEMP=@replicas:{{{}}}
                 FILE={}/HILLS_PTE
                 BIASFACTOR={}
                 ... METAD
                 pte-lw: REWEIGHT_BIAS TEMP=@replicas:{{{}}} ARG=METADPT.bias
                 PRINT FILE=weights.dat ARG=pte-lw,ene STRIDE=250
                        '''.format(sigma, hill_height, temps, hills_file_location, bias_factor, temps))
                break
            else:
                self.log.info('Did not complete the simulation. Replica exchange efficiency of {}.'
                              ' The replica temperatures were {}. The name of the plumed input scripts'
                              ' is \'plumed_wte.dat\''.format(replex_eff, replica_temps))
                continue

        if plumed_output_script is None:
            raise RuntimeError('Did not reach high enough efficiency')
        else:
            #if(hills_file_location != os.getcwd() and hills_file_location is not None):
            #    for f in glob.glob('{}/HILLS_PTE*'.format(hills_file_location)):
            #        self.add_file(f)
            self.save('pte-complete')
            self.remove_file(plumed_input_name)
            return {'plumed': plumed_output_script, 'efficiency': replex_eff, 'temperatures': replica_temps, 'replex': exchange_period}

    def get_simulation(self, sim_tag):
        '''
        give most recent simulation with the same tag
        '''
        for s in self.sims:
            if s.name.startswith(sim_tag):
                return s
        raise KeyError()


    def remove_simulation(self, sim_name, debug=None):
        ''' method that takes a name of the simulation to be removed and deletes the anything related to that simulation

        Parameters:
        sim_name - name of the simulation a.k.a 'tag'
        '''
        sim_index = 0
        full_name = None
        for sim in self.sims:
            if sim.name.startswith(sim_name):
                full_name = sim.name
                if (debug is not None):
                    full_name = debug
                sim_location = sim.location
                del self.sims[sim_index]
                del self._sims[full_name]
                break
            sim_index = sim_index+1
        if (debug is not None):
            full_name = debug
        all_tpr_file_index = 0
        replica_tpr = 0
        tpr_file_index = []
        for tpr_file in self._tpr:
            if type(tpr_file) == list:
                file_path = tpr_file[0].split('/')
                all_tpr_file_index += 1
                if file_path[1] == full_name:
                    tpr_file_index.append(all_tpr_file_index)
                    continue
                elif(file_path[1] != full_name):
                    continue
            else:
                file_path = tpr_file.split('/')
                all_tpr_file_index += 1
                if file_path[1] == full_name:
                    replica_tpr += 1
                    tpr_file_index.append(all_tpr_file_index)
                    continue

        all_gro_file_index = 0
        replica_gro = 0
        gro_file_index = []
        for gro_file in self._gro:
            if type(gro_file) == list:
                file_path = gro_file[0].split('/')
                all_gro_file_index += 1
                if file_path[1] == full_name:
                    gro_file_index.append(all_gro_file_index)
                    replica_gro += 1
                    continue

                elif(file_path[1] != full_name):
                    continue
            else:
                file_path = gro_file.split('/')
                all_gro_file_index += 1
                if (file_path[1] == full_name):
                    gro_file_index.append(all_gro_file_index)
                    continue

        if len(self._gro) < len(gro_file_index):
            raise ValueError('Gro file index is larger than the array size')
        if len(self._tpr) < len(tpr_file_index):
            raise ValueError('Tpr file index is larger than tpr file array')
        if len(self._sims) < sim_index:
            raise ValueError('Simulation index is larger than sims array')
        if type(full_name) != str and full_name is None:
            raise ValueError(
                'name of the simulation to be removed should be a string or the simulation does not exist anymore')
        if os.path.exists(sim_location):
            dir_name = sim_location.split('/')[-1]
            dir_path = sim_location.rpartition('/')[0] + '/'
            bck_count = len(glob.glob(dir_path + dir_name + '-backup-*'))
            shutil.move(sim_location, dir_path + dir_name + '-backup-'+str(bck_count))
        for i in range(len(gro_file_index)):
            del self._gro[gro_file_index[i] - i - 1]

        for i in range(len(tpr_file_index)):
            del self._tpr[tpr_file_index[i] - i - 1]

    def run(self, mdpfile, tag='', repeat=False, mpi_np=None, mdp_kwargs=None, run_kwargs=None, metadata=None, dump_signal=None):
        '''Run a simulation with the given mdpfile

        The name of the simulation will be the name of the mpdfile
        plus the tag name plus a hash of the current gro/top file. If
        the tag/mdp file are the same then the simulation will be restarted

        Parameters
        ----------
        mdpfile : str
            name of the mdp file
        tag : str
            A tag to add onto the simulation.
            Probably good idea if you're adding additional arguments.
        repeat : bool
            This will prevent the simulation for becoming a new name/hash/directory. Useful for continuing a simulation.
        mdp_kwargs : dict or list
            Additional arguments that will be added to the mdp file. Can be list of dicts, which indicates replica exchange
        run_kwargs : dict
            Additional arguments that will be converted to mdrun flags
        metadata : dict
            Metadata which will be saved with the SimulationInfo Object
        dump_signal : signal
            This is the signal that will trigger saving the simulation, the pickle file, and gracefully shutdown

        '''
        if mdp_kwargs is None:
            mdp_kwargs = {}
        if run_kwargs is None:
            run_kwargs = {}
        if metadata is None:
            metadata = {}
        if mpi_np is None:
            mpi_np = self.mpi_np
        if tag == '':
            tag = os.path.basename(mdpfile).split('.')[0]
        # combine kwargs with my config kwargs - using given once as more important
        for k, v in self.run_kwargs.items():
            # avoid override
            if k not in run_kwargs:
                run_kwargs[k] = v

        # make name contain info about mdp and tag
        # remove nsteps, because extending simulations is not something
        # that should lead to new one
        # mdp_kwargs_copy = {k: v for k,v in mdp_kwargs.items() if k != 'nsteps'}
        # file_hash = uuid.uuid5(uuid.NAMESPACE_DNS, mdpfile + str(mdp_kwargs_copy))
        # the hash is huge. Take the first few chars

        simname = tag  # + '-' + str(file_hash)[:8]

        # sick of the stupid hashing. It just breaks stuff
        with self._simulation_context(simname, tag, dump_signal, repeat=repeat) as ec:
            self.log.info('Running simulation with name {}'.format(ec.name))
            ec.metadata.update(metadata)
            self.save('{}-pre'.format(ec.name))
            self._run(mpi_np, mdpfile, ec, mdp_kwargs, run_kwargs)
            self.save('{}-post'.format(ec.name))

    def __del__(self):
        self._stop_logging()

    def _stop_logging(self):
        # gracefully stop logging
        self.log_handler.close()
        self.log.removeHandler(self.log_handler)

    def __getstate__(self):
        odict = self.__dict__.copy()
        del odict['log']
        del odict['log_handler']
        return odict

    def __setstate__(self, dict):
        self.__dict__.update(dict)
        self._start_logging()
        self.log.info('Restarted from pickled state')
        if self.peptidesim_version != __version__:
            self.log.warning('Loaded pickle file is from peptidesim version {} but you are running version {}'.format(
                self.peptidesim_version, __version__))

    def _convert_path(self, p):
        '''Converts path to be local to our working directory'''
        if(os.path.exists(p)):
            return os.path.relpath(os.path.abspath(p), os.path.abspath(self._rel_dir_name))
        else:
            # join(where you are relative to root, filename)
            return os.path.normpath(os.path.join(os.path.relpath(os.getcwd(), os.path.abspath(self._rel_dir_name)), p))

    @contextlib.contextmanager
    def _simulation_context(self, simname, short_name, dump_signal, repeat):
        '''This context will handle restart and keeping a history of simulations performed.
        TODO: Maybe have this context submit a simulation job?
        '''
        # construct signal handler
        def handler(signum, frame):
            self.log.warning(
                'Simulation being interrupted by signal {}'.format(signum))
            rel_dir_name = self._rel_dir_name
            self._rel_dir_name = '.'  # so that we restart from where we started
            self.save('{}-interrupt'.format(simname))
            os.chdir(rel_dir_name)  # put us cleanly into the correct place
            raise KeyboardInterrupt()  # Make sure we do actually end
        if dump_signal is None:
            dump_signal = [
                signal.SIGHUP, signal.SIGINT, signal.SIGQUIT,
                signal.SIGTRAP, signal.SIGABRT, signal.SIGKILL,
                signal.SIGALRM, signal.SIGPIPE, signal.SIGSTOP,
                signal.SIGTSTP, signal.SIGXCPU, signal.SIGVTALRM,
                signal.SIGXFSZ, signal.SIGPROF, signal.SIGTERM]
        else:
            dump_signal = [dump_signal]
        cached = {}
        for ds in dump_signal:
            # cache existing
            oh = signal.getsignal(ds)
            # set new one
            try:
                signal.signal(ds, handler)
                # do not cache, since we cannot set it anyway
                self.log.debug('Will intercept signal ' + str(ds) + ' for saving from simulation')
                cached[ds] = oh
            except (OSError, RuntimeError) as m:
                # skip signals we cannot intercept
                pass

        if simname in self._sims:
            si = self._sims[simname]
        else:
            si = SimulationInfo(simname, short_name)
            if repeat:
                si.location = self._convert_path(os.path.join(
                    self.dir_name, self._sim_list[-1].name))
            else:
                si.location = self._convert_path(
                    os.path.join(self.dir_name, si.name))
            self._sims[simname] = si
            self._sim_list.append(si)
        try:
            yield si
        finally:
            #reset signal handlers
            for k,v in cached.items():
                self.log.debug('Resetting signal handler' + str(k))
                signal.signal(k, v)

    @contextlib.contextmanager
    def _put_in_dir(self, dirname):
        ''' This is a context that will wrap a block of code so that it is executed inside of a particular directory.

        Parameters
        ----------
        dirname : str
            The name of the directory which the function should be executed within. This is relative
            to the dir_name of the PeptideSim object.
        '''
        # check if the path contains our dir_name already
        if(in_directory(dirname, self.dir_name)):
            d = dirname
        else:
            d = self._convert_path(os.path.join(self.dir_name, dirname))
        if not os.path.exists(d):
            os.mkdir(d)
        # bring files
        for f in self.file_list:
            # make sure it exists, is not none and needs to be copied
            if(f is not None and os.path.exists(f) and not os.path.exists(os.path.join(d, os.path.basename(f)))):
                if os.path.isdir(f):
                    shutil.copytree(f, os.path.join(d, os.path.basename(f)))
                else:
                    shutil.copy2(f, os.path.join(d, os.path.basename(f)))
            try:
                # TODO: Why this separate logic?
                if(f.startswith('plumed') and os.path.exists(f)):
                    shutil.copy2(f, os.path.join(d, os.path.basename(f)))
            except shutil.SameFileError:
                pass
        curdir = os.getcwd()
        # keep path to original directory
        self._rel_dir_name = os.path.relpath(curdir, d)
        os.chdir(d)

        try:
            yield

        finally:
            os.chdir(curdir)

            # update path to original directory
            self._rel_dir_name = '.'

            # bring back files
            # TODO: Why? Is this ever necessary
            for f in self.file_list:
                if(f is not None and os.path.exists(os.path.join(d, f)) and not os.path.exists(f)):
                    if os.path.isdir(f):
                        shutil.copytree(os.path.join(d, f), f)
                    else:
                        shutil.copy2(os.path.join(d, f), f)

    def add_file(self, f):
        if not os.path.exists(f):
            raise FileNotFoundError('Could not find requested add_file file {}'.format(f))
        if f != self._convert_path(f):
            try:
                shutil.copyfile(f, self._convert_path(f))
            except shutil.SameFileError:
                pass
        self._file_list.append(self._convert_path(f))

    def remove_file(self, f):
        if not os.path.exists(f):
            raise FileNotFoundError('Could not find requested remove_file file {}'.format(f))
        if f != self._convert_path(f):
            try:
                os.remove(self._convert_path(f))
            except shutil.SameFileError:
                pass
        self._file_list.remove(self._convert_path(f))


    def get_mdpfile(self, f):

        mdpfile = None

        # check if mdp file exist in current directory, on the file list or in dir_name, or in parent of dir_name
        for d in ['.', self.dir_name, os.path.join('..', self.dir_name)]:
            if(os.path.exists(os.path.join(d, f))):
                mdpfile = os.path.join(d, f)

        # now check if it's in our mdp file path
        if(os.path.exists(os.path.join(self.mdp_directory, f))):
            mdpfile = os.path.join(self.mdp_directory, f)

        # now check if it's in our package resource
        if(pkg_resources.resource_exists(__name__, 'templates/' + f)):
            # if so, copy it to here
            self.log.info(
                'Could not locate MDP file {} locally, so using it from package resource.'.format(f))
            with open(f, 'wb') as newf:
                newf.write(pkg_resources.resource_string(
                    __name__, 'templates/' + f))
            mdpfile = f

        if mdpfile is None:
            raise IOError('Could not find MDP file called {}'.format(f))

        return mdpfile

    def _pdb_file_generator(self, sequence, name, amidate=False, acetylate=False):
        '''Generates PDB file from  sequence using BioPython + PeptideBuilder

           Parameters
           ----------
           sequence : str
               The peptide sequence

           Returns
           -------
           tuple
               The pdb file, a list containing the min and maximum extents of the structure, and the molecular weight of the resulting peptide

        '''

        from Bio.PDB import PDBIO
        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        # add glycine if we're acetylating
        if acetylate:
            sequence = 'G' + sequence

        # sets the pdbfilename after the sequence
        pdbfile = name+'.pdb'
        structure = PeptideBuilder.initialize_res(sequence[0])
        for s in sequence[1:]:
            structure = PeptideBuilder.add_residue(structure, s)

        # extract minmax
        smax = [-10**10, -10**10, -10**10]
        smin = [10**10, 10**10, 10**10]
        for a in structure.get_atoms():
            for i in range(3):
                smax[i] = a.coord[i] if a.coord[i] > smax[i] else smax[i]
                smin[i] = a.coord[i] if a.coord[i] < smin[i] else smin[i]

        with self._put_in_dir('peptide_structures'):
            out = PDBIO()

            out.set_structure(structure)
            # adds aminoacids one at a time and generates a pdbfile
            # process amidation/acetyltion
            class Modifier(Bio.PDB.Select):

                def accept_atom(self, atom):
                    rid = atom.get_parent().get_id()[1]
                    if amidate:
                        # remove terminal oxygen
                        if rid == len(sequence) and atom.get_name() == 'O':
                            return True
                    if acetylate:
                        if rid == 1 and atom.get_name() == 'N':
                            return False
                    return True

            out.save(pdbfile, select=Modifier())

            # seems to be the easiest way.....lol?
            if acetylate:
                with open(pdbfile) as f:
                    lines = f.readlines()
                with open(pdbfile, 'w') as f:
                    lines[0] = lines[0].replace('CA ', 'CH3', 1)
                    for i in range(3):
                        lines[i] = lines[i].replace('GLY', 'ACE')
                    f.writelines(lines)


            # get molecular weight
            p = ProteinAnalysis(sequence)

            return (pdbfile, [smin, smax], p.molecular_weight())

    def _packmol(self, output_file='dry_packed.pdb'):
        '''This function takes multiple pdbfiles and combines them into one pdbfile
        '''

        import subprocess

        # compute volume based on density
        mass = sum([c * m for c, m in zip(self.counts, self.peptide_mass)])
        vol = mass / self.peptide_density
        # vol is in L / mol, need to convert to angstrom per molecule
        vol *= 1661

        # sum volumes and get longest dimension
        long_dim = 0
        for e in self._structure_extents:
            diff = [smax - smin for smax, smin in zip(e[0], e[1])]
            long_dim = max(max(diff), long_dim)

        proposed_box_dim = vol**(1/3.)
        proposed_box_dim = max(long_dim, proposed_box_dim)
        self.box_size_angstrom = [proposed_box_dim, sqrt(
            vol / proposed_box_dim), sqrt(vol / proposed_box_dim)]

        if min(self.box_size_angstrom) / 2 < 10.0:
            raise ValueError('The boxsize is too small ({}) to conduct a simulation.'
                             ' Try decreasing the peptide density. Half the shortest box'
                             ' length must be greater than cutoff.'.format(self.box_size_angstrom))

        self.log.info('Final box size to achieve density {} mg/mL: {} Angstrom'.format(
                        self.peptide_density, self.box_size_angstrom))



        # build input text
        input_string = textwrap.dedent(
            '''
                tolerance 2.0
                filetype pdb
                output {}
                '''.format(output_file))
        for f, c in zip(self.peptide_pdb_files, self.counts):
            input_string += textwrap.dedent(
                '''
                    structure {}
                      number {}
                      #make chains be different between molecules
                      changechains
                      #residue numbers are all unique between chains
                      resnumbers 2
                      add_amber_ter
                      inside box 0 0 0 {} {} {}
                    end structure

                    '''.format(f, c, *self.box_size_angstrom))

        with self._put_in_dir('peptide_structures'):

            self.pdb_file = output_file
            input_file = 'packmol.inp'

            with open(input_file, 'w') as f:
                f.write(input_string)

            # pack up packmol into a gromacs command

            result = subprocess.call('{} < {}'.format(
                self.packmol_exe, input_file), shell=True)
            if result != 0:
                self.log.error('Packmol failed with retcode {}. Out: {} Err: {} Input: {input}'.format(
                    result, input=input_string))
            else:
                self.log.info(
                    'Packmol succeeded with retcode {}'.format(result))

            assert os.path.exists(output_file), 'Packmol succeeded with retcode {} but has no output. Input: {input}'.format(
                *result, input=input_string)

    def _pdb2gmx(self, amidate=False, acetylate=False):

        with self._put_in_dir('prep'):

            # only support one kind
            if amidate or acetylate:
                if amidate and acetylate and self.forcefield == 'charmm27':
                    self.pdb2gmx_args['input'] = tuple(sum(self.counts) * ['3', '2'])
                    self.pdb2gmx_args['ter'] = True
                else:
                    raise NotImplementedError('Have not made implemented this combination of FF and ters')

            output = 'dry_mixed.gro'
            topology = 'dry_topology.top'
            self.log.info('Attempting to convert {} to {} with pdb2gmx. Args {}'.format(
                self.pdb_file, output, self.pdb2gmx_args))
            gromacs.pdb2gmx(f=self.pdb_file, o=output, p=topology,
                            water=self.water, ff=self.forcefield, **self.pdb2gmx_args)
            self.gro_file = output
            self.top_file = topology
            # also make pdb file so we have explicit hydrogens
            # cannot use editconf, otherwise chains do not appear
            pdb = 'dry_mixed.pdb'
            gromacs.pdb2gmx(f=self.pdb_file, o=pdb, p=topology,
                            water=self.water, ff=self.forcefield, **self.pdb2gmx_args)
            self.pdb_file = pdb

    def calc_rmsd(self):
        with self._put_in_dir('analysis'):
            self.log.info('Plotting the RMSD of your simulation')
            noe_data = 'noe.dat'
            rmsd_output = 'distrmsd.xvg'
            gromacs.g_rms(f=self.traj_file, s=self.tpr_file,
                          input=['Backbone', 4], o=rmsd_output)
            return rmsd_output

    def calc_sham(self):
        with self._put_in_dir('analysis'):
            self.log.info(
                'Generating free energy landscape of your simulation by boltzman inversion')
            gromacs.g_sham(f='distrmsd.xvg', notime=True, bin='bindex.ndx', lp='probability.xpm',
                           ls='gibbs.xpm', histo='histogram.xvg', lsh='enthalpy.xpm', lss='entropy.xpm')
            return 'energy.xvg', 'bindex.ndx', 'probability.xpm', 'entropy.xpm', 'enthalpy.xpm', 'gibbs.xpm', 'histogram.xvg'

    def _center(self):
        with self._put_in_dir('prep'):
            output = 'dry_centered.gro'
            self.log.info('Centering molecule')
            gromacs.editconf(f=self.gro_file, o=output, c=True,
                             box=self.box_size_nm, bt='cubic')
            self.gro_file = output

    def _solvate(self):

        with self._put_in_dir('prep'):
            output = 'wet_mixed.gro'
            water = self.water + '.gro'

            if self.water == 'spce' or self.water == 'spc' or self.water == 'tip3p':
                # swtich to spc
                water = 'spc216.gro'

            gromacs.solvate(cp=self.gro_file, cs=water, o=output,
                            p=self.top_file, box=self.box_size_nm)
            self.gro_file = output

    def _add_ions(self):

        with self._put_in_dir('prep'):

            # need a TPR file to add ions
            self.log.info('Building first TPR file for adding ions')
            ion_tpr = 'ion.tpr'
            ion_mdp = 'ion.mdp'
            ion_gro = 'prepared.gro'

            mdp_base = gromacs.fileformats.mdp.MDP(
                self.get_mdpfile(self.mdp_base))
            mdp_sim = gromacs.fileformats.mdp.MDP(
                self.get_mdpfile(self.mdp_emin))
            mdp_sim.update(mdp_base)
            mdp_sim.write(ion_mdp)

            gromacs.grompp(f=ion_mdp, c=self.gro_file,
                           p=self.top_file, o=ion_tpr, maxwarn=1)
            # adding max warn because we want to avoid the 'ssytem has net charge' warningbefore output from this grompp command
            self.tpr_file = ion_tpr

            self.log.info('Preparing NDX file')
            ndx_file = 'index.ndx'

            # build ndx input to get an index group for each peptide
            # The one-liner below just explodes the list of peptides/repeats into a
            # list of all peptides

            # get current list of ndx groups
            _, out, _ = gromacs.make_ndx(
                f=self.gro_file, o=ndx_file, input=('', 'q'))
            groups = gromacs.cbook.parse_ndxlist(out)

            input_str = []
            ri = 1  # residue index counter
            name_i = len(groups)  # which group we're naming

            for i, pi in enumerate(
                    reduce(
                        lambda x, y: x + y,
                        [
                            [si for _ in range(ni)]
                            for si, ni in zip(self.sequences, self.counts)
                        ]
                    )):

                input_str.append('r {}-{}'.format(ri, ri + len(pi) - 1))
                ri += len(pi)
                input_str.append('name {} peptide_{}'.format(name_i, i))
                name_i += 1

                input_str.append('{} & a CA'.format(name_i - 1))
                input_str.append('name {} peptide_CA_{}'.format(name_i, i))
                name_i += 1

            # cause it to ouput the final list
            input_str.append('')
            input_str.append('q')

            # update indices
            _, t, _ = gromacs.make_ndx(
                f=self.gro_file, o=ndx_file, input=tuple(input_str))
            self.log.debug('make_ndx output')
            for ti in t.split('\n'):
                self.log.debug('ndx_output: {}'.format(ti))

            solvent_index = -1
            for g in groups:
                if g['name'] == 'SOL':
                    solvent_index = g['nr']
                    break
            assert solvent_index >= 0, 'Problem with making index file and finding solvent'
            self.log.info(
                'Identified {} as the solvent group'.format(solvent_index))

            self.log.info('Adding Ions...')
            gromacs.genion(s=ion_tpr, conc=self.ion_concentration, neutral=True,
                           o=ion_gro, p=self.top_file, input=('', solvent_index))

            # now we need to remove all the include stuff so we can actually pass the file around if needed
            self.log.info('Resovling include statements via GromacsWrapper...')
            self.top_file = gromacs.cbook.create_portable_topology(
                self.top_file, ion_gro)

            self.gro_file = ion_gro
            self.tpr_file = ion_tpr
            self.ndx = ndx_file

    def _run(self, mpi_np, mdpfile, sinfo, mdp_kwargs, run_kwargs):

        with self._put_in_dir(sinfo.location):

            exit_code = 0

            # make this out of restart/no restart logic so we can check for success
            gro = sinfo.short_name + '.gro'
            if type(mdp_kwargs) is list:
                gro = [os.path.join(
                    'multi-{:04d}'.format(i), sinfo.short_name + '.gro') for i in range(len(mdp_kwargs))]

            # check if it's a restart
            if(sinfo.restart_count > 0):
                self.log.info(
                    'Found existing information about this simulation. Using restart')
                # yup, no prep needed
                if(sinfo.complete):
                    self.log.info('Simulation was completed already. Skipping')
                else:
                    # add restart string if this is our first
                    if(sinfo.restart_count == 1):
                        sinfo.run_kwargs['args'] += ' -cpi state.cpt'
                    # allow basic updating of kwargs too, in case they are different
                    for k, v in run_kwargs.items():
                        if ('-' + k) in sinfo.run_kwargs['args'].split():
                            sargs = sinfo.run_kwargs['args'].split()
                            i = sargs.index('-' + k)
                            self.log.info('Noticed different run_kwarg in restart, '
                                          'updating {} with {} (previously, {}'.format(
                                              k, v, sargs[i + 1]))
                            sargs[i + 1] = str(v)
                            sinfo.run_kwargs['args'] = ' '.join(sargs)
                    exit_code = sinfo.run()
            else:
                # need to prepare for simulation
                # preparing tpr file
                self.log.info(
                    'Compiling TPR file for simulation {}'.format(sinfo.name))
                final_mdp = sinfo.short_name + '.mdp'
                mdp_base = gromacs.fileformats.mdp.MDP(
                    self.get_mdpfile(self.mdp_base))
                mdp_sim = gromacs.fileformats.mdp.MDP(
                    self.get_mdpfile(mdpfile))
                # make sure sim has higher priority than base
                mdp_base.update(mdp_sim)
                mdp_sim = mdp_base
                # check if we're doing multiple mdp files
                if type(mdp_kwargs) is list:
                    if not isinstance(mdp_kwargs[0], dict):
                        raise RuntimeError(
                            'To make multiple tpr files, must pass in list of dicts')
                    final_mdp = []
                    mdp_data = []
                    # make a bunch of tpr files and put them into a subdirectories
                    # make a note about how we're making the tprs
                    if(len(self.gro_file_list) == len(mdp_kwargs)):
                        self.log.info('Using list of Gro files from previous  \
                                        run since the number matches the number pf replicas here')
                    multidirs = []
                    for i, mk in enumerate(mdp_kwargs):
                        mdp_temp = gromacs.fileformats.mdp.MDP()
                        mdp_temp.update(mdp_sim)
                        mdp_temp.update(mk)
                        final_mdp.append(sinfo.short_name + str(i) + '.mdp')
                        mdp = gromacs.fileformats.mdp.MDP()
                        mdp.update(mdp_temp)
                        mdp.write(final_mdp[i])
                        mdp_data.append(dict(mdp))
                        tpr = os.path.join(
                            '{}-{:04d}.tpr'.format(sinfo.short_name, i))
                        # get gro file, which may come from a list of same length
                        c = self.gro_file
                        if(len(self.gro_file_list) == len(mdp_kwargs)):
                            c = self.gro_file_list[i]
                        gromacs.grompp(
                            f=final_mdp[i], c=c, p=self.top_file, o=tpr)
                        # move to sub-directory
                        subdir = 'multi-{:04d}'.format(i)
                        multidirs.append(subdir)
                        if not os.path.exists(subdir):
                            os.mkdir(subdir)
                        # feels like overkill, but I guess we copy all files
                        # TPR
                        shutil.copy2(tpr, os.path.join(
                            subdir, sinfo.short_name + '.tpr'))
                        # copy everything(!) on files list
                        for f in self.file_list:
                            # put_in_dir should have already copied it to run directory
                            f = os.path.basename(f)
                            # make sure it exists, is not none and needs to be copied
                            if(f is not None and os.path.exists(f) and not os.path.exists(os.path.join(subdir, f))):
                                if os.path.isdir(f):
                                    shutil.copytree(f, os.path.join(
                                        subdir, os.path.basename(f)))
                                else:
                                    shutil.copy2(f, os.path.join(
                                        subdir, os.path.basename(f)))

                            if f.startswith(
                                    run_kwargs['plumed']) and os.path.exists(
                                        os.path.join(
                                            subdir,
                                            os.path.basename(f))) and os.path.exists(f):
                                shutil.copy2(
                                    f, os.path.join(
                                        subdir, os.path.basename(f)))
                            elif(f.startswith(run_kwargs['plumed']) and not os.path.exists(os.path.join(subdir, os.path.basename(f))) and os.path.exists(f)):
                                shutil.copy2(
                                    f, os.path.join(
                                        subdir, os.path.basename(f)))
                    # keep a reference to current topology. Use 0th since it will exist
                    self.tpr_file = os.path.join(
                        '{}-{:04d}.tpr'.format(sinfo.short_name, 0))
                    # add the multi option
                    run_kwargs.update(dict(multidir=' '.join(multidirs)))
                    sinfo.metadata['md-log'] = os.path.join(
                        multidirs[0], 'md.log')
                    sinfo.metadata['mdp-data'] = mdp_data
                    sinfo.metadata['multi-dirs'] = multidirs
                    tpr = sinfo.short_name + '.tpr'
                else:
                    tpr = sinfo.short_name + '.tpr'
                    mdp = gromacs.fileformats.mdp.MDP()
                    mdp.update(mdp_sim)
                    # make passed highest priority
                    mdp.update(mdp_kwargs)
                    mdp.write(final_mdp)
                    mdp_data = dict(mdp)
                    gromacs.grompp(f=final_mdp, c=self.gro_file,
                                   p=self.top_file, o=tpr)
                    self.tpr_file = tpr
                    sinfo.metadata['md-log'] = 'md.log'
                    sinfo.metadata['mdp-data'] = mdp_data

                # update metadata
                sinfo.metadata['mdp-name'] = mdpfile
                # store reference to trajectory
                if 'o' in run_kwargs:
                    sinfo.metadata['traj'] = run_kwargs['o']
                else:
                    trajectory_file = 'traj.trr'
                    sinfo.metadata['traj'] = trajectory_file
                    run_kwargs['o'] = trajectory_file
                    self.traj_file = trajectory_file
                run_kwargs.update(dict(s=tpr, c=sinfo.short_name + '.gro'))
                sinfo.metadata['run-kwargs'] = run_kwargs
                # add mpiexec to command
                np_str = '' if mpi_np is None else '-np {}'.format(mpi_np)
                # store original driver and prepend mpiexec to it
                temp = gromacs.mdrun.driver
                # also add custom mdrun executable if necessary
                if(self.mdrun_driver is not None):
                    gromacs.mdrun.driver = ' '.join(
                        [self.mpiexec, np_str, self.mdrun_driver])
                else:
                    gromacs.mdrun.driver = ' '.join(
                        [self.mpiexec, np_str, temp])
                self.log.info('Starting simulation {} in directory {}...'.format(
                    sinfo.name, os.getcwd()))
                cmd = gromacs.mdrun._commandline(**run_kwargs)
                gromacs.mdrun.driver = temp  # put back the original command
                self.log.info('Simulation command:' + ' '.join(map(str, cmd)))
                # make it run in shell
                exit_code = sinfo.run(subprocess.call, {
                          'args': ' '.join(map(str, cmd)), 'shell': True})

            # check if the output file was created and exit code
            if exit_code != 0 or ((not os.path.exists(gro[0])) if type(gro) is list else (not os.path.exists(gro))):
                # open the md log and check for error message
                if os.path.exists(sinfo.metadata['md-log']):
                    with open(sinfo.metadata['md-log']) as f:
                        s = f.read()
                        m = re.search(gromacs.mdrun.gmxfatal_pattern,
                                    s, re.VERBOSE | re.DOTALL)
                        if(m is None):
                            self.log.error(
                                'Gromacs simulation (args: {}) failed for unknown reason. Unable to locate output gro file ({})'.format(run_kwargs, gro))
                            self.log.error('Found ({})'.format(
                                [f for f in os.listdir('.') if os.path.isfile(f)]))
                        else:
                            self.log.error('SIMULATION FAILED:')
                            for line in m.group('message'):
                                self.log.error('SIMULATION FAILED: ' + line)
                self.log.error('Removing simulation from restarts')
                self.remove_simulation(sinfo.short_name)
                raise RuntimeError('Failed to complete simulation')
            else:
                self.log.info('{} simulation succeeded'.format(sinfo.name))
            # finished, store any info needed
            self.gro_file = gro

    def save(self, name):

        # we now do not overwrite previously saved pickle files
        #name = '{:05d}-{}.pickle'.format(self._save_count, name)
        name = '{}.pickle'.format(name)
        pickle_path = os.path.join(self._save_directory, name)
        if os.path.exists(pickle_path):
            return
        with open(pickle_path, 'w+b') as f:
            dill.dump(self, file=f)
        shutil.copy2(pickle_path, os.path.join(self._save_directory,
                                               '..', '{}-restart.pickle'.format(self.job_name)))
        with open(os.path.join(self._save_directory, '..', '{}-status.txt'.format(self.job_name)), 'w') as f:
            f.write(str(self) + '\n')
        self._save_count += 1
