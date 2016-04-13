'''Simulate a peptide with a defined sequence and conditions.

    Configuration File:
    -------------------

    config_name: The name of the config file to look at for the sim parameters to use
    
    See the template directory for predefined config files. Run ::

        $ peptidesim --config <config_name> 

    to generate a config file in the current directory based on the config templates provided, or use
    ``default`` to generate the default configuration.

    Module Documentation
    --------------------
    '''
import numpy as np 
import logging, os, shutil, datetime, subprocess, re, textwrap, sys   
import gromacs.tools as tools
import PeptideBuilder 
import Bio.PDB
from .version import __version__
from math import *

from traitlets.config import Configurable, Application, PyFileConfigLoader
from traitlets import Int, Float, Unicode, Bool, List, Instance

class PeptideSim(Configurable):
    '''PeptideSim    

    Example
    -------
    
    Here's an example for creating a ``PeptideSim`` object with the peptide AEAE using the default configuration, saving in the current directory. ::

        p = PeptideSim( dir_name = ".", seqs = ['AEAE'], counts = [1] )

    Here's an example showing **one** AEAE peptide and **two** LGLG peptides, saving in the current directory. ::

        p = PeptideSim( dir_name = ".", seqs = ['AEAE', 'LGLG'], counts = [1,2]) #counts in order of the list of peptides

    '''

    sim_name           = Unicode(u'peptidesim',
                                 help='The name for the type of simulation job (e.g., NVE-equil-NVT-prod)',
                                ).tag(config=True)

    config_file       = Unicode('peptidesim_config.py',
                                help="The config file to load",
                               ).tag(config=True)

    req_files         = List(
                             help='List of files required for simulation. For example, restraints or plumed input'
                            ).tag(config=True)

    pressure          = Float(0,
                              help='Barostat pressure. Ignored if not doing NPT'
                             ).tag(config=True)
    
    peptide_density   = Float(0.2,
                              help='The density of the peptides in milligrams / milliliter',
                             ).tag(config=True)

    log_file          = Unicode(u'simulation.log',
                                 help='The location of the log file. \
                                 If relative path, it will be in \
                                 simulation directory.',
                                ).tag(config=True)


    #Keep a chain of all files created
    _top = []    
    _gro = []
    _pdb = []
    _tpr = []
    _file_list = []

    @property
    def file_list(self):
        result = self._file_list[:]
        result.extend([self.pdb_file, self.gro_file, self.top_file, self.tpr_file])
        return result

    @property
    def pdb_file(self):
        if(len(self._pdb) == 0):
            return None


    @pdb_file.setter
    def pdb_file(self, f):
        self._pdb.append(os.path.abspath(f))

    @property
    def gro_file(self):
        if(len(self._gro) == 0):
            return None


    @gro_file.setter
    def gro_file(self, f):
        self._gro.append(os.path.abspath(f))

    @property
    def top_file(self):
        if(len(self._top) == 0):
            return None


    @top_file.setter
    def top_file(self, f):
        self._top.append(os.path.abspath(f))

    @property
    def tpr_file(self):
        if(len(self._tpr) == 0):
            return None
        return self._tpr[-1]

    @tpr_file.setter
    def tpr_file(self, f):
        self._tpr.append(os.path.abspath(f))
        
    
                     
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

        #Set-up directory to begin with
        self.job_name = job_name
        if job_name is None:
            self.job_name = os.path.split(dir_name)[-1]

        self.dir_name = dir_name
        
        if not os.path.exists(self.dir_name):
            os.mkdir(self.dir_name)
        self.dir_name = os.path.abspath(self.dir_name)
        for f in self.req_files:
            shutil.copyfile(f, self._convert_path(f))
            self._file_list.append(os.path.basename(f))


        #check if logger is relative
        #split path and see if folder is empty
        if(len(os.path.split(dir_name)[0]) == 0):
            self.log_file = self._convert_path(self.log_file)

        #don't know how we got here, so we'll just add our logger
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s [%(filename)s, %(lineno)d, %(funcName)s]: %(message)s (%(levelname)%s)")
        file_handler.setFormatter(formatter)
        
        self.log = logging.getLogger()
        self.log.addHandler(file_handler)
        
        self.log_handler = file_handler
        self.log.info('Started logging for PeptideSim...{}')

        
        #Note: use load_pyconfig_files to merge them. Useful in future
        #load in configuration file
        config = config_file
        if config_file is not None:
            config = PyFileConfigLoader(config_file).load_config()
        super(PeptideSim, self).__init__(config=config)

        #generate pdbs from sequences and store their extents
        self.structure_extents = []
        self.peptide_mass = []
        self.peptide_pdb_files = []
        for i, sequence in enumerate(seqs):
            print(i, sequence)
            structure, minmax, mass = self._pdb_file_generator(sequence,'seq_' + str(i))
            self.peptide_pdb_files.append(structure)
            self.structure_extents.append(minmax)
            self.peptide_mass.append(mass)

        #store the copies we'd like to have of each sequence
        if counts is None:
            counts = [1 for s in seqs]
        self.counts=counts

    def __del__(self):
        #gracefully stop logging
        self.log.removeHandler(self.log)


    def _convert_path(self, p, dir='.'):
        '''Converts path and optional subdirectory to be local to our working directory'''
        return os.path.join(self.dir_name, dir, os.path.basename(p))
    
    def _put_in_dir(dirname):
        ''' This is a decorator that will wrap a function so that it is executed inside of a particular directory.
	
        Parameters
        ----------
        dirname : str
            The name of the directory which the function should be exectued within. This is relative
            to the dir_name of the PeptideSim object.
        Returns
        -------
        fxn
            A new wrapped function (usually used as annotation)
        '''
        
        def wrap(fxn):
            def mod_f(self, *args, **kwargs):
                d = self._convert_path(dirname)
                if not os.path.exists(d):
                    os.mkdir(d)
                #bring files
                for f in self.file_list:
                    if(f is not None and os.path.exists(f)):
                        shutil.copyfile(f, os.path.join(d, os.path.basename(f)))           
                    #go there
                    curdir = os.getcwd()
                    os.chdir(d)
                    try:
                        return fxn(self, *args, **kwargs)
                    finally:
                        #make sure we leave
                        os.chdir(curdir)
                        #bring back files
                        for f in self.file_list:
                            if(f is not None and os.path.exists(os.path.join(d, f))):
                                shutil.copyfile(os.path.join(d, f),f)
            return mod_f
        return wrap

        
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

        Returns                                                                                                                         ------                                                                                                                             
        '''
        pass

    @_put_in_dir('peptide_structures')
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
                
        out = PDBIO()
        out.set_structure(structure)
        out.save( pdbfile ) #adds aminoacids one at a time and generates a pdbfile

        #get molecular weight
        p = ProteinAnalysis(sequence)        
        
        return (pdbfile, [smin, smax], p.molecular_weight())
        
    @_put_in_dir('packing')
    def _packmol(self, output_file='dry_packed.pdb'):
        '''This function takes multiple pdbfiles and combines them into one pdbfile
        '''


        #compute volume based on density
        mass = sum([c * m for c,m in zip(self.counts, self.peptide_mass)])
        vol = mass / self.peptide_density
        
        #sum volumes and get longest dimension
        long = 0
        for e in self.structure_extents:
            diff = [max - min for max,min in e]
            long = max(max(diff), long)

        proposed_box_dim = vol**(1/3.)
        proposed_box_dim = max(long, proposed_box_dim)
        box_size = [proposed_box_dim, sqrt(vol / proposed_box_dim), sqrt(vol / proposed_box_dim)]
        
        
        #build input file
        input_file = 'input.inp'
        with open(input_file, 'w') as f:
            f.write(textwrap.dedent(
                '''
                tolerance 2.0
                filetype pdb 
                output {}
                '''.format(self.dry_packed_file)))
            for fname in self.peptide_pdb_files:#iterate through peptide pdb structures
                f.write(textwrap.dedent(
                    '''
                    structure {} 
                      number {}
                      inside box 0 0 0. {}. {}. {}. 
                    end structure
                    '''.format(fname, self.counts, *box_size))) 
        self.pdb_file = output_file

        #pack up packmol into a gromacs command 
        

class PeptideSimConfigurator(Application):
    #information for running peptidesim from command line
    name = 'peptidesim'
    version = __version__
    classes = [PeptideSim]
    description = '''The peptidesim application.

    Currently, this application only generates a configuration
    '''
    

    #command line flags
    aliases = {'f': 'PeptideSim.config_file',
               'config': 'PeptideSim.generate_config'}
    
    config_file = Unicode('peptidesim_config.py',
        help="The config file to load",
    ).tag(config=True)
                   
    generate_config = Bool(True,
        help="Generate default config file",
    ).tag(config=True)


    def write_config_file(self):
        '''Write our default config to a .py config file'''
        if os.path.exists(self.config_file):
            answer = ''
            def ask():
                prompt = "Overwrite {} with default config? [y/N]".format(self.config_file)
                try:
                    return raw_input(prompt).lower() or 'n'
                except KeyboardInterrupt:
                    print('') # empty line
                    return 'n'
            answer = ask()
            while not answer.startswith(('y', 'n')):
                print("Please answer 'yes' or 'no'")
                answer = ask()
            if answer.startswith('n'):
                return

        config_text = self.generate_config_file()
        if isinstance(config_text, bytes):
            config_text = config_text.decode('utf8')
        print("Writing default config to: %s" % self.config_file)
        with open(self.config_file, mode='w') as f:
            f.write(config_text)
            






def main():
    p = PeptideSimConfigurator.instance()
    p.write_config_file()
    
if __name__ == "__main__":
    main()
