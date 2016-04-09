'''Simulate a peptide with a defined sequence and conditions.

    Configuration File:
    -------------------

    config_name: The name of the config file to look at for the sim parameters to use
    
    See the template directory for predefined config files. Run ::

        $ peptidesim config <config_name> 

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

from traitlets.config.configurable import Configurable
from traitlets import Int, Float, Unicode, Bool, List

PDB2GMX='gmx pdb2gmx'
GMXSOLVATE='gmx solvate'

class PeptideSim(Configurable):
    '''PeptideSim    

    Example
    -------
    
    Here's an example for creating a ``PeptideSim`` object with the peptide AEAE using the default configuration, saving in the current directory. ::

        p = PeptideSim( dir_name = ".", seqs = ['AEAE'], counts = [1] )

    Here's an example showing **one** AEAE peptide and **two** LGLG peptides, saving in the current directory. ::

        p = PeptideSim( dir_name = ".", seqs = ['AEAE', 'LGLG'], counts = [1,2]) #counts in order of the list of peptides
'''

    name = Unicode(u'peptidesim',
                   help='The name for the type of simulation job (e.g., NVE-equil-NVT-prod)'
                   ).tag(config=True)
                   
    topol = Unicode(u'topology.top',
                    help='Gromacs topology file'
                    ).tag(config=True)
    gro = Unicode(u'protein.gro',
                  help='The Gromcas structure file'
                  ).tag(config=True)
    pressure = Float(0,
                     help='Barostat pressure. Ignored if not doing NPT'
                     ).tag(config=True)
    file_list = List(help = 'Files to carry over for each simulation'
                     ).tag(config=True)
    pdbfiles = List(help = 'Initial PDB files for peptides'
                     ).tag(config=True)

    dry_packed_file = Unicode(u'dry_mix.pdb',
                              help = 'The name of the combined peptides without water'
                              ).tag(config=True)

    peptide_density = Float(0.2,
                          help='The density of the peptides in milligrams / milliliter'
                          ).tag(config=True)
    
                     
    def __init__(self,job_name,seqs,counts=None):        
        '''This is an initiator that takes the arguments from the command
line and creates the class simulation.

        Parameters
        ----------
        job_name : str
            name of the directory where your simulation should be saved, and 
        seqs : List[str]
            A list of amino acid sequences.
        counts : List[int]
            A list of the number of occurrences of each amino acid, in order.
        '''
        self.job_name = job_name
        #generate pdbs from sequences and store their extents
        self.structure_extents = []
        self.peptide_mass = []
        for sequence in seqs:
<<<<<<< HEAD
            self.pdbfile=self._pdb_file_generator(sequence)#sequence will be converted to a pdb file
            self.pdbfiles.append(self.pdbfile)#copies the pdbfiles
            self.files_tocopy.append(self.pdbfile)#copies the files needed to start the simulation
        self._files_to_take()#this actually copies and carries the files around
        self._setup_directory(*self._files_to_take())# sets up a directory with given name and files
        self.copies1=1#number of copies of pdb file 
        self.x_dim_box,self.y_m_box,self.z_dim_box=40,40,40#box legth in Angstroms
        self.copies2=1#number of copies of pdb file 
    

    def change_density(self, new_density):
        '''This will change the density of the simulation.

        Density is given by the equation:
            rho = m / v
        where m is mass and v is volume. This function
        will...


        Parameters
        -----------
        new_density : float
            The density at which the simulation should be do'''
    def equiliberation(self,dirname):
        '''This is the method of equiliberation which contains five steps.

        First step: Adding water into the system and the output goes to file named peptide_addwater_gro.
        
        Second step: Adding ions into the system and the output goes to file named peptide_addwater_addions_gro. 

        Third step: Getting energy minimization  
        
        Fourth step: Annealing for simulated annealing 

        Fifth step: Equil run

        Parameters
        -------------
        dirname: the name of directory that has output files of the equiliberation simulation
=======
            structure, minmax, mass = self._pdb_file_generator(sequence)
            self.pdbfiles.append(structure)
            self.structure_extents.append(minmax)
            self.peptide_mass.append(mass)
            
        #we'll carry around the pdb files
        self.file_list.append(self.pdbfiles)        
        self._setup_directory(*self.file_list)# sets up a directory with given name and files

        #store the copies we'd like to have of each sequence
        if counts is None:
            counts = [1 for s in seqs]
        self.counts=counts

        #pack initial structure
        
       
>>>>>>> 6add3844966b37560aefa36cdff1544d083239bc
        
        .....                                                                                                                              
        Returns                                                                                                                            
        ------                                                                                                                             
        res : res is the path to md.log file

        '''
        pass

        '''

    def _putInDir(dirname):
        ''' This is a decorator that will wrap a function so that it is executed inside of a particular directory.

        Args:
            dirname: The name of the directory which the function should be exectued within.
        returns:
            A new wrapped function (usually used as annotation)
        '''
        def wrap(fxn):
            def mod_f(self, *args):
                if not os.path.exists(dirname):
                    os.mkdir(dirname)
                    #bring files
                for f in self.file_list:
                    if(f is not None and os.path.exists(f)):
                        shutil.copyfile(f, os.path.join(dirname, os.path.basename(f)))            
                #go there
                os.chdir(dirname)
                try:
                    fxn(self, *args)
                finally:
                    #make sure we leave
                    os.chdir( self.dir)
                    #bring back files
                    for f in self.file_list:
                        if(f is not None and os.path.exists(os.path.join(dirname, f))):
                            shutil.copyfile(os.path.join(dirname, f),f)
            return mod_f
        return wrap


    def _exec_log(self, string, arg_dic=None, input=None):
        '''
        Smartly executes the given string and logs the results
        Args: 
             string: a string of command line commands
             arg_dic: dictionary of arguments needed to execute the command
             input: input string if a command needs and input from a user
        '''
        if(arg_dic is not None):
            for k,v in arg_dic.iteritems():
                string += ' -{} {}'.format(k,v)
        logging.info('{}'.format(string))
        try:
            process = subprocess.Popen(string, shell=True,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       stdin=subprocess.PIPE)
                           
            if(input is None):
                out, err = process.communicate()
            else:
                out, err = process.communicate(input)
            retcode = process.returncode
            logging.debug(out)
            if(len(err) > 0):
                logging.warning(err)
                logging.debug(out)
            if retcode < 0:
                logging.error('{} failed with {}. stderr follow'.format(string, retcode))
                logging.error(err)
                raise OSError('{} failed with {}. stderr follow'.format(string, retcode))
        except OSError as e:
            raise OSError('Execution of {} failed:'.format(string), e)

    def _setup_directory(self, *to_copy):
        '''builds directory, starts log
           *to_copy: takes a list of files to be copied, can take as many files needed 
           returns: nothing
        '''

        if not os.path.exists(self.name):
            os.mkdir(self.name)
        self.dir = os.path.abspath(self.name)
        for f in to_copy:
            shutil.copyfile(f, os.path.join(self.dir, os.path.basename(f)))
        os.chdir(self.dir)
        #set-up logging
        logging.basicConfig(filename='simulation.log',level=logging.DEBUG)
        logging.info('Beginning simulation {}'.format(datetime.datetime.now().strftime('%b-%d-%I%M%p-%G')))

    def _files_to_take(self):
        '''This function returns an array of files to be copied
           args:none
           returns: an array of files to be copied
        '''
        return self.file_list#takes the files that are needed for the simulation
    def _pdb_file_generator(self, sequence):        
        '''This function generates a pdbfile using peptide builder from a string of Amino Acids.
           Args:
               sequence-a string of amino acids
           returns: a pdbfile   
        '''

        from Bio.PDB import PDBIO
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        
        pdbfile=sequence+'.pdb'# sets the pdbfilename after the sequence
        structure = PeptideBuilder.initialize_res(sequence[0])
        for s in sequence[1:]:
            structure = PeptideBuilder.add_residue(structure, s)
        
        #extract minmax
        smax = [-10**10, -10**10, -10**10]
        smin = [10**10, 10**10, 10**10]
        for a in structure.get_atoms():
            for i in range(3):
                smax[i] =  a.coords[i] if a.coords[i] >  smax[i] else smax[i]
                smin[i] =  a.coords[i] if a.coords[i] <  smin[i] else smin[i]
                
        out = PDBIO()
        out.set_structure(structure)
        out.save( pdbfile ) #adds aminoacids one at a time and generates a pdbfile

        #get molecular weight
        p = ProteinAnalysis(sequence)        
        
        return pdbfile, [smin, smax], p.molecular_weight()
        
        
    @_putInDir('init_setup')
    def packmol(self):
        '''This function takes multiple pdbfiles and combines them into one pdbfile
        '''

        #compute the box size we'll need
        #sum volumes and get longest dimension
        vol = 0
        long = 0
        for e in self.structure_extents:
            diff = [max - min for max,min in e]
            long = max(max(diff), long)
            vol += reduce(lambda x,y: x * y, diff)


        vol *= self.volume_margin
        proposed_box_dim = vol**(1/3.)
        proposed_box_dim = max(long, proposed_box_dim)
        
        
        #build input file
        input_file = 'input.inp'
        with open(input_file, 'w') as f:
            f.write(textwrap.dedent(
                '''
                tolerance 1.0
                filetype pdb 
                output {}
                '''.format(self.dry_packed_file)))
            for i in self.pdbfiles:#iterate through peptide pdb structures
                f.write(textwrap.dedent(
                    '''
                    structure {} 
                      number {}
                      inside box 0 0 0. {}. {}. {}. 
                    end structure
                    '''.format(i, self.counts, self.x_dim_box, self.y_dim_box, self.z_dim_box))) 
        PACKMOL1='packmol < {}'.format(input_file) #packmol command that runs the input file and generates a pdbfile with multiple amino acid sequences.
        self._exec_log(PACKMOL1)#executes the packmol command
        self.file_list.append(self.pdbfile)
        self.file_list.append(input_file)
        
    @_putInDir('initial_housekeeping')        
    def initial_setup(self):
        '''Prepares the components before the simulation is run.This is where the type of water and type of force fields are chosen. Generates a box with the amino acid sequence. It generates gro file with amino acids in a box and itp files. 
        Args: none
        Returns:none
        '''
        self.file_list
        self.gro='initial.gro' #gro file that contains the amino acid.
        self.water_model=['spc', 'spce', 'tip3p', 'tip4p', 'tip5p']#types of water models available
        self.water_model_number=1#index of the water model to be used
        self.box_protein='box_protein.gro'#gro file that contains the amino acid in a box.
        
        force_fields={'AMBER03':1, 'AMBER94':2, 'AMBER96':3, 'AMBER99':4, 'AMBER99SB':5,'AMBER99SB_ILDN':6, 'AMBERGS':7, 'CHARMM27':8, 'GROMOS96_43a1':9, 'GROMOS96_43a2':10, 'GROMOS96_5a3':11, 'GROMOS96_53a5':12, 'GROMOS96_53a6':13, 'GROMOS96_54a7':14, 'OPLS_AAL':15, 'all_atom_solvent':16, 'all_atom_vac':17, 'gromacs':18, 'gromacs_nmr_H':19}# types of force fields to choose for intial setup of the simulation
       
        arg_force={'f':self.pdbfile,'o':self.gro, 'p':self.topol, 'water':self.water_model[self.water_model_number]}
        self._exec_log(PDB2GMX, arg_dic=arg_force, input='{}'.format(force_fields['GROMOS96_43a1'])+'\n')#executes the conversion of pdbfile to a gro file with a given water model and force field
        
        box_protein=tools.Editconf(f=self.gro, c=[], o=self.box_protein, box=[self.x_dim_box/10, self.y_dim_box/10, self.z_dim_box/10])#runs an editconf command from gromacs_wrapper
        box_protein.run()#runs the editconf command
      
        self.file_list.append(self.topol)
        self.file_list.append(self.box_protein)
        self._itp()#finds and takes itp files
    @_putInDir('energy_min_solvate')    
    def solvate_enem(self):
        ''' This function will solvate the box. Runs energy minimization and adds ions into the system. Generates a lot of files but takes in topology and gro file after the finish of the simulation.
        Args: none
        returns: none
        '''
        self.file_list
        self.solvated_gro="solvated.gro" #gro file that will contain solvated box of amino acid sequences
        topol_em_tpr='topol_em.tpr'# name of the tpr file generated by initial energy minimization
       
        self.conf_ion_gro="ions.gro"#gro file after ions are added
        topol_em_tpr2="topol_ions.tpr"#name of the tpr file generated after  energy minimization after adding ions
        energy_mdrun_trr='energy.trr'
        if(self.water_model_number==0 or self.water_model_number==1 or self.water_model_number==2):
            self.water_model_gro='spc216.gro'
        else:
            self.water_model_gro=self.water_model[self.water_model_number]+'.gro'#defines the gro file for corresponding water models defined in initial_setup
            
        args_solvated={'cp':self.box_protein, 'cs':self.water_model_gro, 'o':self.solvated_gro, 'p':self.topol, 'box':'{} {} {}'.format(self.x_dim_box/10, self.y_dim_box/10, self.z_dim_box/10)}#water will be added to the box
        self._exec_log(GMXSOLVATE, arg_dic=args_solvated) #solvates the box
        energy_min=tools.Grompp( f=self._energy(), c=self.solvated_gro, p=self.topol, o=topol_em_tpr, pp='processed1.top')#arguments for running energy minimization
        energy_min.run()#runs the energy minimization from gromacs wrapper commands
        genion=tools.Genion(s=topol_em_tpr, p=self.topol, o=self.conf_ion_gro, conc=0.001)#arguments needed to add ions to the system
        genion.run()#adds ions
        energy_min2=tools.Grompp( f=self._energy(), c=self.conf_ion_gro, p=self.topol, o=topol_em_tpr2)#arguments needed to run grommp commands
        energy_min2.run()#second energy minimization after adding ions is running
        energy_mdrun=tools.Mdrun(s=topol_em_tpr2, c=self.conf_ion_gro, o=energy_mdrun_trr)#arguments needed to run mdrun command
        energy_mdrun.run()#runs the mdrun command
        self.file_list.append(self.conf_ion_gro)#takes the output of the last command
        self.file_list.append(self.topol)

    @_putInDir('equilibration')
    def equilibrate(self):
        ''' Equilibrates the system. Copies the topolgy and gro file after the end of the simulation.
        Args: none
        returns: none
        '''
        self.file_list#copies files needed to run this equilibration simulation
        self.equilib_tpr='equilib.tpr'#tpr file generated after running the equilibration
        self.equil_top='equil.top'#name of the topology file after equilibration step
        self.equilib_mdrun_trr='equilib.trr'#trr file generated after running equilibration
        self.equil_gro='equil.gro' #gro file generated after running equilibration
        equilibrate_grommp=tools.Grompp(f=self._equilib_file(), c=self.conf_ion_gro, p=self.topol,pp=self.equil_top, o=self.equilib_tpr)
        equilibrate_grommp.run()#reads the equilibration file and runs the grommp command
        equilibrate_mdrun=tools.Mdrun(s=self.equilib_tpr, c=self.equil_gro, o=self.equilib_mdrun_trr )
        equilibrate_mdrun.run()#run the equilibration mdrun
        self.file_list.append(self.equil_gro)
        self.file_list.append(self.equil_top)
    @_putInDir('nvt')
    def nvt_grompp_mdrun(self):
        '''Runs the simulation in nvt ensemble
        Args:
        Returns: 
        '''
        self.file_list
        self.nvt_tpr='nvt.tpr'#tpr file generated by grompp command
        self.nvt_top='nvt.top'#top file
        self.nvt_trr='nvt.trr'#trr file
        self.nvt_gro='nvt.gro'#gro file
        nvt_grompp=tools.Grompp(f=self._nvt_mdp(), c=self.equil_gro, p=self.equil_top, pp=self.nvt_top, o=self.nvt_tpr)
        nvt_grompp.run()#runs grompp
        nvt_mdrun=tools.Mdrun(s=self.nvt_tpr, c=self.nvt_gro, o=self.nvt_trr)
        nvt_mdrun.run()#runs mdrun
        self.file_list.append(self.nvt_gro)
        self.file_list.append(self.nvt_trr)
    def _equilib_file(self):
        '''mdp file needed for equilibration
        Args: none
        returns: 
          input_file-an mdp input file for running equilbration 
        '''
        self.file_list
        input_file='equilib.mdp'
        
        with open(input_file, 'w') as f:
            f.write(textwrap.dedent('''
            ; RUN CONTROL PARAMETERS = 
            integrator               = sd
            nsteps                   = {time}  
            dt                       = 0.0005
            
            
            
            ; Output frequency for coords (x), velocities (v) and forces (f) = 
            nstxout                  = 1000
            nstvout                  = 0
            nstfout                  = 0
            
            ; Output frequency for energies to log file and energy file = 
            nstlog                   = 1000
            nstenergy                = 0
            
                       
            ;Neighbor searching neighbor
            cutoff-scheme		=Verlet

            ; OPTIONS FOR ELECTROSTATICS AND VDW = 
            ; Method for doing electrostatics = 
            coulombtype     = PME
            nstcalclr	=1 ;Controls the period between calculations of long-range forces when using the group cut-off scheme.
            rlist           = 0.01
            optimize-fft    = yes
            fourierspacing  = 0.12
            pme-order       = 4
            ewald-rtol      = 1e-5
            dispcorr                 = ener
            ;Temperature
            tcoupl                   = v-rescale ;Temperature coupling using velocity rescaling with a stochastic term (JCP 126, 014101). This thermostat is similar to Berendsen coupling, with the same scaling using tau-t, but the stochastic term ensures that a proper canonical ensemble is generated. The random seed is set with ld-seed. This thermostat works correctly even for tau-t=0. For NVT simulations the conserved energy quantity is written to the energy and log file.
            tau_t                    = 2 ;[ps] time constant for coupling (one for each group in tc-grps), -1 means no temperature coupling
            ref_t                    = 300 ;[K] reference temperature for coupling (one for each group in tc-grps)
            tc-grps                  = System ;groups to couple separately to temperature bath
            
            ;Pressure            
            pcoupl                   = berendsen ;Exponential relaxation pressure coupling with time constant tau-p [ps]. The box is scaled every timestep. It has been argued that this does not yield a correct thermodynamic ensemble, but it is the most efficient way to scale a box at the beginning of a run.
            tau_p                    = 250 ;[ps] time constant for coupling
            compressibility          = 4.5e-5
            ref_p                    = {pressure}
            ;Constraints
            constraints              = h-bonds ;Convert the bonds with H-atoms to constraints.
            lincs_iter               = 2 ;Number of iterations to correct for rotational lengthening in LINCS. For normal runs a single step is sufficient, but for NVE runs where you want to conserve energy accurately or for accurate energy minimization you might want to increase it to 2.
            lincs_order              = 6 ;Highest order in the expansion of the constraint coupling matrix. When constraints form triangles, an additional expansion of the same order is applied on top of the normal expansion only for the couplings within such triangles. For ``normal'' MD simulations an order of 4 usually suffices, 6 is needed for large time-steps with virtual sites or BD. For accurate energy minimization an order of 8 or more might be required. With domain decomposition, the cell size is limited by the distance spanned by lincs-order+1 constraints. When one wants to scale further than this limit, one can decrease lincs-order and increase lincs-iter, since the accuracy does not deteriorate when (1+lincs-iter)*lincs-order remains constant.
            '''.format(pressure=self.pressure * 1.01325, 
                       time=self.equil_time * 10**6 / 0.5)))
            return input_file
        
    def _itp(self):
        '''Finds files that have an itp extension and appends them to files to be taken. 
        Args: none
        returns: none
        '''
        self.itp_files = []
        with open(self.topol) as f:
            dir = os.path.dirname(self.topol)
            for line in f:
                m = re.match(r'#include "(\S+)\.itp"', line)
                if(m):
                    file = '{}.itp'.format(m.group(1))
                    if(os.path.exists(os.path.join(dir, file))):
                        self.itp_files.append(file)
                        #print self.itp_files
                        self.file_list.append(file)

                        
    def _energy(self):
        '''
        Energy minimize the system. 
        Args: none
        return: mdp_file(an mdp file needed to run an energy minimization)
        '''
        self.energy_mdp_file= 'emin.mdp'
        with open(self.energy_mdp_file, 'w') as f:
            f.write(textwrap.dedent('''
            ; RUN CONTROL PARAMETERS = 
            integrator               = steep
            ;A steepest descent algorithm for energy minimization. The maximum step size is emstep [nm], the tolerance is emtol [kJ mol-1 nm-1].
            emtol                    = 0.01 ; in [kJ mol-1 nm-1]the minimization            ;is converged when the maximum force is smaller than this value  
            
            nsteps                   = 100; max steps
            
            
            ; Output frequency for coords (x), velocities (v) and forces (f) = 
            nstxout                  = 5
            nstvout                  = 0
            nstfout                  = 0
            
            ; Output frequency for energies to log file and energy file = 
            nstlog                   = 1000
            nstenergy                = 0
            
           
            
            ;Neighbor searching neighbor
            cutoff-scheme		=Verlet;Generate a pair list with buffering. The buffer size is automatically set based on verlet-buffer-tolerance, unless this is set to -1, in which case rlist will be used. This option has an explicit, exact cut-off at rvdw=rcoulomb. Currently only cut-off, reaction-field, PME electrostatics and plain LJ are supported. Some mdrun functionality is not yet supported with the Verlet scheme, but grompp checks for this. Native GPU acceleration is only supported with Verlet. With GPU-accelerated PME or with separate PME ranks, mdrun will automatically tune the CPU/GPU load balance by scaling rcoulomb and the grid spacing. This can be turned off with -notunepme. Verlet is faster than group when there is no water, or if group would use a pair-list buffer to conserve energy.

            ; OPTIONS FOR ELECTROSTATICS AND VDW = 
            ; Method for doing electrostatics = 
            coulombtype     = PME ;Fast smooth Particle-Mesh Ewald (SPME) electrostatics. Direct space is similar to the Ewald sum, while the reciprocal part is performed with FFTs. Grid dimensions are controlled with fourierspacing and the interpolation order with pme-order. With a grid spacing of 0.1 nm and cubic interpolation the electrostatic forces have an accuracy of 2-3*10-4. Since the error from the vdw-cutoff is larger than this you might try 0.15 nm. When running in parallel the interpolation parallelizes better than the FFT, so try decreasing grid dimensions while increasing interpolation.

            rlist           = 0.01 ; [nm]
            ;Cut-off distance for the short-range neighbor list. With cutoff-scheme=Verlet, this is by default set by the verlet-buffer-tolerance option and the value of rlist is ignored.
            optimize-fft    = yes
            fourierspacing  = 0.12 ;For ordinary Ewald, the ratio of the box dimensions and the spacing determines a lower bound for the number of wave vectors to use in each (signed) direction. For PME and P3M, that ratio determines a lower bound for the number of Fourier-space grid points that will be used along that axis. In all cases, the number for each direction can be overridden by entering a non-zero value for fourier_n[xyz]. For optimizing the relative load of the particle-particle interactions and the mesh part of PME, it is useful to know that the accuracy of the electrostatics remains nearly constant when the Coulomb cut-off and the PME grid spacing are scaled by the same factor.
            pme-order       = 4 ;Interpolation order for PME. 4 equals cubic interpolation. You might try 6/8/10 when running in parallel and simultaneously decrease grid dimension.
            ewald-rtol      = 1e-5 ;The relative strength of the Ewald-shifted direct potential at rcoulomb is given by ewald-rtol. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum. 
            dispcorr        = ener ;apply long range dispersion corrections for Energy only
            ;constraints    = h-bonds

            '''.format()))
        
        return  self.energy_mdp_file
    
    def _nvt_mdp(self):
        '''generates an mdp file for running an nvt ensemble simulation
        Args: none
        returns: input_ file(an mdp file needed to run the nvt ensemble simulation)
        '''
        input_file = 'prod.mdp'
        with open(input_file, 'w') as f:
            f.write(textwrap.dedent('''
            ; RUN CONTROL PARAMETERS = 
            integrator               = sd
            nsteps                   = {time}
            dt                       = 0.002
            
            
            
            ; Output frequency for coords (x), velocities (v) and forces (f) = 
            nstxout                  = 1000
            nstvout                  = 0
            nstfout                  = 0
            
            ; Output frequency for energies to log file and energy file = 
            nstlog                   = 1000
            nstenergy                = 0
            
            
            ; Cut-off scheme:
	    cutoff-scheme		=Verlet
           
            ; OPTIONS FOR ELECTROSTATICS AND VDW = 
            ; Method for doing electrostatics = 
            coulombtype     = PME
            rlist           = 0.01            
            optimize-fft    = yes
            fourierspacing  = 0.12
            pme-order       = 4
            ewald-rtol      = 1e-5
            dispcorr                 = ener
            ;Temperature
            tcoupl                   = v-rescale
            tau_t                    = 2
            ref_t                    = 300
            tc-grps                  = System
            
            ;Constraints
            constraints              = h-bonds
            lincs_iter               = 2
            lincs_order              = 6
            '''.format(pressure=self.pressure * 1.01325, 
                       time=self.prod_time * 10**6 / 2.)))
            return input_file

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
