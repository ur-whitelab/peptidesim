'''Simulate a peptide with a defined sequence and conditions.

Example
-------

Here's an example for creating a ``PeptideSim`` object with the peptide AEAE using the default configuration, saving in the current directory. ::

    p = PeptideSim( dir_name = ".", seqs = ['AEAE'], counts = [1] )

Here's an example showing **one** AEAE peptide and **two** LGLG peptides, saving in the current directory. ::

    p = PeptideSim( dir_name = ".", seqs = ['AEAE', 'LGLG'], counts = [1,2]) #counts in order of the list of peptides
'''
import numpy as np 

import logging, os, shutil, datetime, subprocess, re, textwrap, sys, pkg_resources, contextlib, uuid, json, ast, requests, signal

import PeptideBuilder 
import Bio.PDB
from math import *
from .utilities import *
import dill

from traitlets.config import Configurable, Application, PyFileConfigLoader
from traitlets import Int, Float, Unicode, Bool, List, Instance, Dict

import gromacs
gromacs.environment.flags['capture_output'] = True


class SimulationInfo(object):
    '''A class that stores information about a simulation. Only for keeping history of simulation runs
    '''

    def __init__(self,name, short_name):
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

    def run(self,run_fxn=None, run_kwargs=None):
        if(self.restart_count > 0):            
            if (run_fxn is not None and run_fxn != self.run_fxn) or (run_kwargs is not None and run_kwargs != self.run_kwargs):
                raise ValueError('Name collision in simulation. You tried to repeat a non-identical simulation')
        else:
            self.run_fxn = run_fxn
            self.run_kwargs = run_kwargs

        #have to add now in case we die            
        self.restart_count += 1 
        result = self.run_fxn(**self.run_kwargs)

        self.complete = True
        return result
        

    


class PeptideSim(Configurable):
    '''PeptideSim class. Class to use for conducting simulations.

    '''

    sim_name           = Unicode(u'peptidesim',
                                 help='The name for the type of simulation job (e.g., NVE-equil-NVT-prod)',
                                ).tag(config=True)

    config_file       = Unicode(u'peptidesim_config.py',
                                help="The config file to load",
                               ).tag(config=True)

    req_files         = List(
                             help='List of files required for simulation. For example, restraints or plumed input'
                            ).tag(config=True)

    pressure          = Float(0,
                              help='Barostat pressure. Ignored if not doing NPT'
                             ).tag(config=True)
    
    peptide_density   = Float(0.02,
                              help='The density of the peptides in milligrams / milliliter',
                             ).tag(config=True)

    ion_concentration = Float(0.002,
                              help='The concentration of sodium chloride to add in moles / liter'
                              ).tag(config=True)

    log_file          = Unicode(u'simulation.log',
                                 help='The location of the log file. If relative path, it will be in simulation directory.',
                                ).tag(config=True)
    packmol_exe       = Unicode(u'packmol',
                                help='The command to run the packmol program.'
                                ).tag(config=True)

    forcefield        = Unicode(u'amber99sb-ildn',
                                help='The gromacs syntax forcefield',
                                ).tag(config=True)
    water             = Unicode(u'tip3p',
                               help='The water model to use',
                               ).tag(config=True)

    pdb2gmx_args      = Dict(dict(ignh=None),
                             help='Any additional special arguments to give to pdb2gmx, aside from force-field and water which are separately specified.',
                             ).tag(config=True)

    mdp_directory     = Unicode(u'.',
                                help='The directory to find gromacs MDP files'
                            ).tag(config=True)
    mdp_base          = Unicode(u'peptidesim_base.mdp',
                                help='The MDP file containing basic forcefield parameters'
                       ).tag(config=True)
    
    mdp_emin          = Unicode(u'peptidesim_emin.mdp',
                                help='The energy miniziation MDP file. Built from mdp_base. Used specifically for adding ions'
                       ).tag(config=True)

    host              = Unicode(u'http://52.71.14.39',
                                help='The host address for the Redis database'
                        ).tag(config=True)
    
    post_address      = Unicode(u'/insert/simulation',
                                help='The extension to post simulation data'
                        ).tag(config=True)


                                  
    @property
    def box_size_angstrom(self):
        return self._box_size

    @box_size_angstrom.setter
    def box_size_angstrom(self, v):
        assert(len(v) == 3)
        self._box_size[:] = v[:]        

    @property
    def box_size_nm(self):        
        return [x / 10 for x in self._box_size]

    @box_size_nm.setter
    def box_size_nm(self, v):
        assert(len(v) == 3)
        self._box_size = [x * 10 for x in v]
    

    @property
    def file_list(self):
        result = []
        result.extend(self._file_list)
        #for now, we'll keep these as paths relative to basedirectory. So we won't have local copies everywhere
        #result.extend([self.pdb_file, self.gro_file, self.top_file, self.tpr_file])
        return result

    @property
    def pdb_file(self):
        if(len(self._pdb) == 0):
            return None
        return os.path.normpath(os.path.join(self.rel_dir_name, self._pdb[-1]))


    @pdb_file.setter
    def pdb_file(self, f):
        self._pdb.append(self._convert_path(f))

    @property
    def gro_file(self):
        if(len(self._gro) == 0):
            return None
        return os.path.normpath(os.path.join(self.rel_dir_name, self._gro[-1]))


    @gro_file.setter
    def gro_file(self, f):
        self._gro.append(self._convert_path(f))

    @property
    def top_file(self):
        if(len(self._top) == 0):
            return None
        return os.path.normpath(os.path.join(self.rel_dir_name, self._top[-1]))


    @top_file.setter
    def top_file(self, f):
        self._top.append(self._convert_path(f))

    @property
    def tpr_file(self):
        if(len(self._tpr) == 0):
            return None
        return os.path.normpath(os.path.join(self.rel_dir_name, self._tpr[-1]))

    @tpr_file.setter
    def tpr_file(self, f):
        self._tpr.append(self._convert_path(f))        
 
    @property
    def ndx(self):
        if(self._ndx_file is None):
            return None
        n = gromacs.fileformats.NDX()
        n.read(os.path.normpath(os.path.join(self.rel_dir_name, self._ndx_file)))
        return n

    @property 
    def sims(self):
        return self._sim_list

    @ndx.setter
    def ndx(self, n):
        self._ndx_file = self._convert_path(n)
        self._file_list.append(self._ndx_file)

    def __init__(self,dir_name,seqs,counts=None,config_file=None, job_name=None):        
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

        #We have to declare class attributes at the instance level
        #for pickling purposes. Silly, I know.

        #Keep a chain of all files created. Hide behind properties
        self._top = []    
        self._gro = []
        self._pdb = []
        self._tpr = []
        self._file_list = []
        self._ndx_file = None
        
        #keep track of simulations in case of restart needs        
        self._sims = dict()
        self._sim_list = []
        
        #other variables
        self._box_size = [0,0,0] #box size in angstroms



        #Set-up directory to begin with
        self.job_name = job_name
        if job_name is None:
            self.job_name = os.path.split(dir_name)[-1]

        self.dir_name = dir_name
        self.rel_dir_name = '.'
        
        if not os.path.exists(self.dir_name):
            os.mkdir(self.dir_name)

        for f in self.req_files:
            shutil.copyfile(f, self._convert_path(f))
            self._file_list.append(os.path.basename(f))

        self._start_logging()

        
        #Note: use load_pyconfig_files to merge them. Useful in future
        #load in configuration file
        config = config_file
        if config_file is not None:
            self.log.info('Loading config file {}'.format(config_file))
            config = PyFileConfigLoader(config_file).load_config()
        else:
            self.log.info('Using default configuration'.format(config_file))
            
        self.log.debug('Loaded {}:'.format(str(config)))            
        super(PeptideSim, self).__init__(config=config, parent=None)

        #store passed parameters
        self.sequences = seqs
        #store the copies we'd like to have of each sequence
        if counts is None:
            counts = [1 for s in seqs]
        self.counts=counts
        self.log.info('Have {} many of these sequences {}:'.format(counts, seqs))

    def _start_logging(self):
        #check if logger is relative
        #split path and see if folder is empty
        if(self.log_file == os.path.basename(self.log_file)):
            self.log_file = os.path.join(self.dir_name, self.log_file)

        #don't know how we got here, so we'll just add our logger
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s [%(filename)s, %(lineno)d, %(funcName)s]: %(message)s (%(levelname)s)")
        file_handler.setFormatter(formatter)
        
        self.log = logging.getLogger('peptidesim:{}'.format(self.job_name))
        self.log.addHandler(file_handler)
        
        self.log_handler = file_handler
        self.log.setLevel(logging.DEBUG)
        self.log.info('Started logging for PeptideSim...')

    def store_data(self):
        '''Writes the instance's Traits to a json file
        '''
        if not os.path.exists('data'):
            os.mkdir('data')
        with open('data/simdata.json', 'w' ) as f:
            data = {}
            for k, v in self.traits().iteritems():
                if type(v.default_value) not in [unicode,int, float]:
                    data[k] = ast.literal_eval(v.default_value_repr())
                else:
                    data[k] = v.default_value
            for k,v in self.__dict__.iteritems():
                if k not in data and type(v) in [unicode, int, float, list, dict,tuple, str] and k[0] != '_':
                    data[k] = v
            f.write(json.dumps(data))

            #put information to database

            #properties that we want to store
            properties = ['peptide_density', 'ion_concentration']
            url = (self.host + self.post_address).encode('utf-8')

            i = 0
            for seq in data['sequences']:
                for prop in properties:
                    payload = {'sim_name': self.sim_name, 'property':prop, 'property_value':data[prop]}
                    r = requests.put(url, payload)
                payload = {'sim_name': self.sim_name, 'property':'counts', 'property_value':data['counts'][i]}
                r = requests.put(url, payload)
                i+=1    


    def initialize(self):
        '''Build PDB files, pack them, convert to gmx, add water and ions

        This method accomplishes the following steps:
          1. Use Bio Python to convert sequences into PDB files
          2. Combine the PDB files using packmol
          3. Convert them into gmx files using the pdb2gmx command and configuration parameters
          4. Add water using the editconf/genbox
          5. Compile our energy-minimization tpr file for purposes of adding ions
          6. Add ions

        '''
        #generate pdbs from sequences and store their extents
        self.structure_extents = []
        self.peptide_mass = []
        self.peptide_pdb_files = []
        for i, s in enumerate(self.sequences):
            structure, minmax, mass = self._pdb_file_generator(s,'seq_' + str(i))
            self.peptide_pdb_files.append(structure)
            self.structure_extents.append(minmax)
            self.peptide_mass.append(mass)


        #pack the peptides together into an initial structure
        self._packmol()

        #now get gromcas files
        self._pdb2gmx()

        #Add solvent
        self._solvate()

        #add ions
        self._add_ions()


        self.log.info('Completed Initialization')

    def run(self, mdpfile, tag='', mpi_np=1, mdp_kwargs=dict(), run_kwargs=dict(), metadata=dict(), pickle_name=None, dump_signal=signal.SIGTERM):
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
        mdp_kwargs : dict or list
            Additional arguments that will be added to the mdp file. Can be list of dcits, which indicates replica exchange
        run_kwargs : dict
            Additional arguments that will be convreted to mdrun flags
        '''
        if pickle_name is None:
            pickle_name = self.job_name + '.pickle'
        with self._simulation_context(os.path.basename(mdpfile).split('.')[0] + '-' + tag, pickle_name, dump_signal) as ec:
            self.log.info('Running simulation with name {}'.format(ec.name))
            ec.metadata.update(metadata)
            self._run(mpi_np, mdpfile, ec, mdp_kwargs, run_kwargs)


    def energy_minimize(self, tag, steps=1000):
        '''Energy minimize the system. Can be called anytime after initialize.        
        '''

            

    def __del__(self):
        self._stop_logging()

    def _stop_logging(self):
        #gracefully stop logging
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
        self.log.info('Restarting from pickled state')

        
    def _convert_path(self, p):
        '''Converts path to be local to our working directory'''
        if(os.path.exists(p)):
            return os.path.relpath(os.path.abspath(p), os.path.abspath(self.rel_dir_name))
        else:
            #join(where you are relative to root, filename)
            return os.path.normpath(os.path.join(os.path.relpath(os.getcwd(), os.path.abspath(self.rel_dir_name)), p))


    @contextlib.contextmanager
    def _simulation_context(self, name, pickle_name, dump_signal = signal.SIGTERM):
        '''This context will handle restart and keeping a history of simulations performed.

        TODO: Maybe have this context submit a simulation job?
        '''

        #construct signal handler
        def handler(signum, frame):
            with open(os.path.join(self.rel_dir_name,pickle_name), 'w') as f:
                dill.dump(self, file=f)
            os.chdir(self.rel_dir_name) #put us cleanly into the correct place
            raise KeyboardInterrupt #Make sure we do actually end
        #cache existing
        oh = signal.getsignal(dump_signal)
        #set new one
        signal.signal(dump_signal, handler)

        #construct name and add to simulation infos
        file_hash = uuid.uuid5(uuid.NAMESPACE_DNS, self.top_file + self.gro_file + self.pdb_file)
        simname = name + '-' + str(file_hash)

        if simname in self._sims:
            si = self._sims[simname]
        else:
            si = SimulationInfo(simname, name)            

        self._sims[simname] = si
        self._sim_list.append(si)

        try:
            yield  si    
        finally:
            #reset signal handler
            signal.signal(dump_signal, oh)

        
            

    
    @contextlib.contextmanager
    def _put_in_dir(self, dirname):
        ''' This is a context that will wrap a block of code so that it is executed inside of a particular directory.
	
        Parameters
        ----------
        dirname : str
            The name of the directory which the function should be exectued within. This is relative
            to the dir_name of the PeptideSim object.
        '''
        d = self._convert_path(os.path.join(self.dir_name, dirname))
        if not os.path.exists(d):
            os.mkdir(d)
        #bring files
        for f in self.file_list:
            if(f is not None and os.path.exists(f)):
                try:
                    shutil.copyfile(f, os.path.join(d, os.path.basename(f)))
                except shutil.Error:
                    pass #same file. Stupid that this is necessary
                
        #go there
        curdir = os.getcwd()
        #keep path to original directory
        self.rel_dir_name = os.path.relpath(curdir, d)
        os.chdir(d)


        try:
            yield

        finally:
            os.chdir(curdir)

            #update path to original directory
            self.rel_dir_name  = '.'
            
            #bring back files
            for f in self.file_list:
                if(f is not None and os.path.exists(os.path.join(d, f))):
                    try:
                        shutil.copyfile(os.path.join(d, f),f)
                    except shutil.Error:
                        pass #same file


    def add_file(self, f):
        if f != self._convert_path(f):        
            shutil.copyfile(f, self._convert_path(f))
        self._file_list.append(os.path.basename(f))

                    
    def get_mdpfile(self, f):

        mdpfile = None
        
        #check if mdp file exist in current directory, on the file list or in dir_name, or in parent of dir_name
        for d in ['.', self.dir_name, os.path.join('..', self.dir_name)]:
            if(os.path.exists(os.path.join(d, f))):
                mdpfile =  os.path.join(d, f)

        #now check if it's in our mdp file path
        if(os.path.exists(os.path.join(self.mdp_directory, f))):
            mdpfile =  os.path.join(self.mdp_directory, f)
        

        #now check if it's in our package resource
        if(pkg_resources.resource_exists(__name__, 'templates/' + f)):
            #if so, copy it to here
            self.log.info('Could not located MDP file {} locally, so using it from package resource.'.format(f))
            with open(f, 'wb') as newf:
                newf.write(pkg_resources.resource_string(__name__, 'templates/' + f))
            mdpfile = f

        if mdpfile is None:
            raise IOError('Could not find MDP file called {}'.format(f))

        return mdpfile
                
    def analysis(self):
            """This function analyzes the output of the simulation. 

            It reads md.log file specified in configuration. It extracts the total
            energy and temperature at each timestep and creates respective
            histograms. The function creates a folder called Images and saves the
            historams as TotalEnergyHist.png and TemperatureHist.png.
                            
            Returns
            -------
            list
                The relative path to the *.png images.

            """
            
            pass
        

    def equilibration(self,dirname):
        '''This is the method of equiliberation which contains five steps.

        First step: Adding water into the system and the output goes to file named peptide_addwater_gro.
        
        Second step: Adding ions into the system and the output goes to file named peptide_addwater_addions_gro. 

        Third step: Getting energy minimization  
        
        Fourth step: Annealing for simulated annealing 

        Fifth step: Equil run

        Parameters
        -------------
        dirname : str 
            the name of directory that has output files of the equiliberation simulation

        Returns                                                                                                          
        -------
        '''
        pass

    def _pdb_file_generator(self, sequence, name):
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
        
        pdbfile=name+'.pdb'# sets the pdbfilename after the sequence
        structure = PeptideBuilder.initialize_res(sequence[0])
        for s in sequence[1:]:
            structure = PeptideBuilder.add_residue(structure, s)
        
        #extract minmax
        smax = [-10**10, -10**10, -10**10]
        smin = [10**10, 10**10, 10**10]
        for a in structure.get_atoms():
            for i in range(3):
                smax[i] =  a.coord[i] if a.coord[i] >  smax[i] else smax[i]
                smin[i] =  a.coord[i] if a.coord[i] <  smin[i] else smin[i]

        with self._put_in_dir('peptide_structures'):
            out = PDBIO()
            out.set_structure(structure)
            out.save( pdbfile ) #adds aminoacids one at a time and generates a pdbfile

            #get molecular weight
            p = ProteinAnalysis(sequence)        

            return (pdbfile, [smin, smax], p.molecular_weight())
        
    def _packmol(self, output_file='dry_packed.pdb'):
        '''This function takes multiple pdbfiles and combines them into one pdbfile
        '''

        import subprocess

        #compute volume based on density
        mass = sum([c * m for c,m in zip(self.counts, self.peptide_mass)])
        vol = mass / self.peptide_density
        
        #sum volumes and get longest dimension
        long_dim = 0
        for e in self.structure_extents:
            diff = [smax - smin for smax,smin in zip(e[0], e[1])]
            long_dim = max(max(diff), long_dim)

        proposed_box_dim = vol**(1/3.)
        proposed_box_dim = max(long_dim, proposed_box_dim)
        self.box_size_angstrom = [proposed_box_dim, sqrt(vol / proposed_box_dim), sqrt(vol / proposed_box_dim)]
        
        
        #build input text
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
                      inside box 0 0 0 {} {} {}
                    end structure
                    '''.format(f, c, *self.box_size_angstrom))


        with self._put_in_dir('peptide_structures'):

            self.pdb_file = output_file
        

            #pack up packmol into a gromacs command
            class Packmol(gromacs.core.Command):
                command_name = self.packmol_exe            
            cmd = Packmol()
            result = cmd(input=input_string)
        
            if result[0] != 0:
                self.log.error('Packmol failed with retcode {}. Out: {} Err: {} Input: {input}'.format(*result, input=input_string))
            else:
                self.log.info('Packmol succeeded with retcode {}'.format(*result))

            assert os.path.exists(output_file), 'Packmol claimed to succeed but no output file found'

    def _pdb2gmx(self):

        with self._put_in_dir('prep'):

            output = 'dry_mixed.gro'
            topology = 'dry_topology.top'
            self.log.info('Attempting to convert {} to {} with pdb2gmx'.format(self.pdb_file, output))
            gromacs.pdb2gmx(f=self.pdb_file, o=output, p=topology, water=self.water, ff=self.forcefield, **self.pdb2gmx_args)
            self.gro_file = output
            self.top_file = topology        

    def _solvate(self):
        
        with self._put_in_dir('prep'):
            output = 'wet_mixed.gro'
            water = self.water + '.gro'

            if self.water == 'tip3p':
                #swtich to spc
                water = 'spc216.gro'
            gromacs.solvate(cp=self.gro_file, cs=water, o=output, p=self.top_file, box=self.box_size_nm)
            self.gro_file = output


    def _add_ions(self):

        with self._put_in_dir('prep'):
        
            #need a TPR file to add ions
            self.log.info('Building first TPR file for adding ions')
            ion_tpr = 'ion.tpr'
            ion_mdp = 'ion.mdp'
            ion_gro = 'prepared.gro'
            
            
            #NOTE: edit_mdp rewrites, it does not load
            mdp_base = gromacs.cbook.read_mdp(self.get_mdpfile(self.mdp_base))
            gromacs.cbook.edit_mdp(self.get_mdpfile(self.mdp_emin),new_mdp=ion_mdp, **mdp_base)
            
            gromacs.grompp(f=ion_mdp, c=self.gro_file, p=self.top_file, o=ion_tpr)
            self.tpr_file = ion_tpr
            
            
            self.log.info('Preparing NDX file')
            ndx_file = 'index.ndx'
            
            #build ndx input to get an index group for each peptide
            #The one-liner below just explodes the list of peptides/repeats into a
            #list of all peptides
            
            #get current list of ndx groups
            _,out,_  = gromacs.make_ndx(f=self.gro_file, o=ndx_file, input=('', 'q'))
            groups = gromacs.cbook.parse_ndxlist(out)
            
            input_str = []
            ri = 1 #residue index counter
            name_i =  len(groups) #which group we're naming

            for i, pi in enumerate(                                                    \
                reduce(                                                                \
                    lambda x,y: x + y,                                                 \
                    [                                                                  \
                        [si for _ in xrange(ni)]                                       \
                        for si, ni in zip(self.sequences,self.counts)                  \
                    ]                                                                  \
                    )):                                                                
                
                
                input_str.append('r {}-{}'.format(ri, ri + len(pi) - 1))
                ri += len(pi)
                input_str.append('name {} peptide_{}'.format(name_i, i))
                name_i += 1
                
                input_str.append('{} & a CA'.format(name_i - 1))
                input_str.append('name {} peptide_CA_{}'.format(name_i, i))
                name_i += 1

            #cause it to ouput the final list
            input_str.append('')                            
            input_str.append('q')

            #update indices
            _,t,_  = gromacs.make_ndx(f=self.gro_file, o=ndx_file, input=tuple(input_str))
            self.log.debug('make_ndx output')
            for ti in t.split('\n'):
                self.log.debug('ndx_output: {}'.format(ti))
            
            solvent_index = -1
            for g in groups:
                if g['name'] == 'SOL':
                    solvent_index = g['nr']
                    break
            assert solvent_index >= 0, 'Problem with making index file and finding solvent'
            self.log.info('Identified {} as the solvent group'.format(solvent_index))
        

            self.log.info('Adding Ions...')
            gromacs.genion(s=ion_tpr, conc=self.ion_concentration, neutral=True, o=ion_gro, p=self.top_file, input=('', solvent_index))
            self.log.info('...OK')
            
            #now we need to remove all the include stuff so we can actually pass the file around if needed
            self.log.info('Resovling include statements via GromacsWrapper...')
            output = self.top_file = gromacs.cbook.create_portable_topology(self.top_file, ion_gro)
            self.log.info('...OK')
            
            self.gro_file = ion_gro
            self.tpr_file = ion_tpr
            self.top_file = output
            self.ndx = ndx_file

    def _run(self, mpi_np, mdpfile, sinfo, mdp_kwargs, run_kwargs):

        if(mpi_np == 1):
            #assuming debug. Fixing the openmp number
            run_kwargs.update({'nt': '1'})
        
        with self._put_in_dir(sinfo.name):

            #make the simulation info an absolute path
            sinfo.location = os.path.abspath(os.getcwd())

            #make this out of restart/no restart logic so we can check for success
            gro = sinfo.short_name + '.gro'
            if isinstance(mdp_kwargs,list):                    
                gro = sinfo.short_name + '0.gro'

            #check if it's a restart
            if(sinfo.restart_count > 0):
                self.log.info('Found existing information about this simulation. Using restart')
                #yup, no prep needed
                if(sinfo.complete):
                    self.log.info('Simulation was completed already. Skipping')
                else:
                    #add restart string if this is our first
                    if(sinfo.restart_count == 1):
                        sinfo.run_kwargs['args'] += ' -cpi state.cpt'
                    sinfo.run()
            else:
                #need to prepare for simulation                
                #Preparing emin tpr file        
                self.log.info('Compiling TPR file for simulation {}'.format(sinfo.name))
                final_mdp = sinfo.short_name + '.mdp'
                
                mdp_base = gromacs.cbook.read_mdp(self.get_mdpfile(self.mdp_base))

                #check if we're doing multiple mdp files
                if isinstance(mdp_kwargs,list):                    
                    assert isinstance(mdp_kwargs[0], dict), 'To make multiple tpr files, must pass in list of dicts'

                    final_mdp = []
                    mdp_data = []
                    top_dir = 'TOPOL'
                    
                    if not os.path.exists(top_dir):
                        os.mkdir(top_dir)

                    #make a bunch of tpr files and put them into a subdirectory (TOPOL)
                    for i, mk in enumerate(mdp_kwargs):
                        mdp_temp = mdp_base.copy()
                        mdp_temp.update(mk)
                        final_mdp.append(sinfo.short_name + str(i) + '.mdp')
                        gromacs.cbook.edit_mdp(self.get_mdpfile(mdpfile), new_mdp=final_mdp[i], **mdp_temp)
                        mdp_data.append(gromacs.cbook.read_mdp(final_mdp[i]))
                        tpr = os.path.join(top_dir, sinfo.short_name + str(i) + '.tpr')
                        gromacs.grompp(f=final_mdp[i], c=self.gro_file, p=self.top_file, o=tpr)
                    tpr = os.path.join(top_dir, sinfo.short_name)

                    #keep a reference to current topology. Use 0th since it will exist
                    self.tpr_file = tpr + '0.tpr'
                    
                    #add the multi option
                    run_kwargs.update(dict(multi=len(mdp_kwargs)))

                    sinfo.metadata['md-log'] = 'md0.log'
                    sinfo.metadata['mdp-data'] = mdp_data
                        

                else:
                    tpr = sinfo.short_name + '.tpr'
                    mdp_base.update(mdp_kwargs)
                    gromacs.cbook.edit_mdp(self.get_mdpfile(mdpfile), new_mdp=final_mdp, **mdp_base)
                    mdp_data = gromacs.cbook.read_mdp(final_mdp)
                    gromacs.grompp(f=final_mdp, c=self.gro_file, p=self.top_file, o=tpr)
                    self.tpr_file = tpr                    
                    sinfo.metadata['md-log'] = 'md.log'
                    sinfo.metadata['mdp-data'] = mdp_data

                #update metadata
                sinfo.metadata['mdp-name'] = mdpfile

                
                run_kwargs.update(dict(s=tpr, c=sinfo.short_name + '.gro', dds=0.5))

                sinfo.metadata['run-kwargs'] = run_kwargs

                #add mpiexec to command                
                #store original driver and prepend mpiexec to it
                temp = gromacs.mdrun.driver
                gromacs.mdrun.driver = ' '.join(['mpiexec.hydra', '-np {}'.format(mpi_np), temp])

                self.log.info('Starting simulation...'.format(sinfo.name))
                cmd = gromacs.mdrun._commandline(**run_kwargs)
                gromacs.mdrun.driver = temp #put back the original command
                self.log.info(cmd)
                self.log.info(' '.join(map(str, cmd)))
                #make it run in shell
                sinfo.run(subprocess.call, {'args': ' '.join(map(str,cmd)), 'shell':True})
                #sinfo.run(gromacs.mdrun, run_kwargs)
                #check if the output file was created                
            if(not os.path.exists(gro)):
                #open the md log and check for error message
                with open(sinfo.metadata['md-log']) as f:
                    s = f.read()
                    m = re.search(gromacs.mdrun.gmxfatal_pattern, s, re.VERBOSE | re.DOTALL)
                    if(m is None):                        
                        self.log.error('Gromacs simulation failed for unknown reason. Unable to locate output gro file ({})'.format(gro))
                        self.log.error('Found ({})'.format([f for f in os.listdir('.') if os.path.isfile(f)]))
                    else:
                        self.log.error('SIMULATION FAILED:') 
                        for line in m.group('message'):
                            self.log.error('SIMULATION FAILED: ' + line)
            else:
                self.log.info('...done'.format(sinfo.name))            
                #finished, store any info needed
                self.store_data()
                self.gro_file = gro
            
            


