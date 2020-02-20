# PeptideSim
​
## Goal

We want to be able to post-process explore/analyze the data. We need to know the following then:

   1. The peptides/amounts, density, etc (traits)
   2. The mdp parameters
   3. the trajectory locations
   4. metadata system of simulations

## Installation

### Bluehive Install

1. Load required modules
```bash
module load anaconda3 packmol libmatheval gromacs-plumed/2018.3/b4 openblas openmpi
```

2. Create a virtual environment
```bash
conda create -n yourenvname python=2.7
```
Make sure you specify _python 2.7_ for installing your conda environment as PeptideSim is not compatible with the newer versions of python. Once you create your environment, you should activate it for the next steps:
```bash
source activate yourenvname
```
**Note**: Avoid using `conda activate yourenvname` as it might cause some problems in installing different python modules using pip. Check on pip and python and make sure they are sourced to your environment directory:
```bash
which pip #or which python
# The output should look like this:
# /.conda/envs/yourenvname/bin/pip
```
**Note**: In case the environment is not properly activated, i.e., `which pip` outputs a path not as mentioned above, try `conda init`. Close and reopen a new terminal and do `conda activate <yourenvname>`.

3. You will need to install the GromacsWrapper module. You can
download the [source code](https://github.com/Becksteinlab/GromacsWrapper/releases) for version 0.8.0 (latest) and install using the following command:
```bash
pip install GromacsWrapper-release-0.8.0.tar.gz
```

4. Setup a config file for Gromacswrapper to be able to find the gromacs installation. Use interactive python to do this:
```python
import gromacs
gromacs.config.setup()
```
This will create a `.gromacswrapper.cfg` file in your home directory.
Using any text editor, make the changes to that file and add the GMXRC location. Your file should look like this:
```text
 [ DEFAULT ]
 qscriptdir = %( configdir )s/qscripts
 templatesdir = %( configdir )s/templates
 configdir = ~/.gromacswrapper
​
 [ Gromacs ]
 release = 2018.3
 gmxrc = /software/gromacs-plumed/2018.3/b4/bin/GMXRC
 extra =
 tools = gmx gmx_mpi
 groups = tools
​
 [ Logging ]
 logfilename = gromacs .log
 loglevel_console = INFO
 loglevel_file = DEBUG
```

5. For installing PeptideSim, you need to clone the package from github.
```bash
git clone https://github.com/ur-whitelab/peptidesim.git
```
Change directory to `package` and install the requirements and the module:
```bash
cd package
pip install -r requirements.txt
pip install .
```

### Developer Environment

First, prepare the docker image in the test-docker folder by running
the build script. Then, run the test script in the root directory. It
is only necessary to rebuild the docker script when newer gromacs,
gromacswrapper packages are available.

## Typical Workflow

### Regular Simulation/Multi_chain
Refer to the scripts under peptidesim/inputs/simple/ as an example.
1. Create a Python script simple.py and import required module.
```python
from pepsidesim import PeptideSim
```
2. Specify input condition and initialize
```python
seq = 'EKEKEKEKEKEK' #input the sequence in one-letter-code
name = 'EK6' #naming the sequence
pep_copies = 1 #specifiy the number of copies, increase the number if running multiple chains to study self assembly and etc.
MPI_NP = 4
ps = PeptideSim(name, [seq], [pep_copies], job_name='2mer_{}'.format(name)) #input to PeptideSim
ps.peptide_density = 0.008 #g/mol
ps.water = 'spce' #specify water model
ps.mpi_np = MPI_NP
ps.initialize()
```
3. Specify energy minimization, annealing and NVT tuning input parameters
```python
ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**5})
ps.run(mdpfile='peptidesim_anneal.mdp', tag='anneal_nvt')
ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod', mdp_kwargs={'nsteps': int(3 * 5*10**5), 'constraints': 'h-bonds'})
```
4. Create a bash script and run it
```bash
#!/bin/sh
#SBATCH --partition=standard --time=120:00:00 --output=EK6.txt
#SBATCH --mail-type ALL
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH --mem=24gb

export OMP_NUM_THREADS=2
python simple.py
```
5. Adding pickle

Simluation may takes a long time depending on the number of steps inputed and resource requrested and it may be the case that the job won't be able to finish even with the upper limit time of Bluehive. You can use pickle module to help resuming the job. Firs step is to  add if statement before initializing:
```python
if(os.path.exists(pickle_name)):
    print('loading restart')
    with open(pickle_name, 'rb') as f:
        ps = pickle.load(f)
else:	
    #step 2.
```
It is recommended to dump information into pickle after each step for example:
```python
ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**5})
print(ps.pickle_name, 'picklename1')
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

ps.run(mdpfile='peptidesim_anneal.mdp', tag='anneal_nvt')
print(ps.pickle_name, 'picklename2')
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod', mdp_kwargs={'nsteps': int(3 * 5*10**5), 'constraints': 'h-bonds'})
print(ps.pickle_name, 'picklename3')
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

```
Note that for Python3, pickle opened with binary mode must be specified. Otherwise error will occur.

### PT-WTE/EDS
Refer to the scripts under peptidesim/inputs/pte/ as example
1. Preparation

Create a Python script for the preparation step. It is almost identical to the regular simulation. (Check the part1.py)

If you wish to add EDS, you will need additional code for the preparation. Check part1.py as well as https://www.plumed.org/doc-v2.5/user-doc/html/_m_e_t_a_d.html.

2. PT-WTE

Create a second script and import all required modules. 
Specify input. Note that the name of the job must be consistant with that in part 1.

```python
name = nameinpart1
pickle_name = name + '.pickle'
MPI_NP = 16
peptide_cpoies = 1 #number of peptide per replica
replicas = 16 #number of replicas
remd_exhcange_period = 250
#reload
ps=3#initialize
if(os.path.exists(pickle_name)):
    print 'loading restart'
    with open(pickle_name, 'r') as f:
        ps = pickle.load(f)
        #ps.pickle_name=pickle_name
        print os.getcwd()
        ps.rel_dir_name='.'
```

Write a function for getting replica exchange efficiency, this will be used in determining whether to end the simulation
```python
def get_replex_e(ps, replica_number):
    with open(ps.sims[-1].location + '/' + ps.sims[-1].metadata['md-log']) as f:
        p1 = re.compile('Repl  average probabilities:')
        p2 = re.compile('Repl\s*' + ''.join(['([0-9\.]+\s*)' for _ in range(replica_number - 1)]) + '$')
        ready = False
        answer=-1
        for line in f.readlines():
            if not ready and p1.findall(line):
                ready = True
            elif ready:
                match = p2.match(line)
                
                if match:
                    answer= [float(s) for s in match.groups()]
                    break
            
        return answer
```

Specify PTE inputs and generate plumed file
```python
pte_result = ps.pte_replica(mpi_np=MPI_NP, max_tries=3,min_iters=1, mdp_kwargs={'nsteps':int(400* 5*10**2), }, replicas=replicas,hills_file_location=os.getcwd(),hot=400,hill_height=0.6,sigma=250,bias_factor=16, eff_threshold=0.20,cold=278,exchange_period=200)
pte_plumed_script=pte_result['plumed']
replica_temps=pte_result['temperatures']
with open('plumed_pte.dat', 'w') as f:
    f.write(pte_plumed_script)
ps.add_file('plumed_pte.dat')
temps = ','.join(str(e) for e in replica_temps)
kwargs = [ {'ref_t': ti} for ti in replica_temps]
time_ns=5.0 
replex_eff = 0.2
for kw in kwargs:
    kw['nsteps']=  int(time_ns * 5 * 10 ** 5)
max_iterations=20
min_iterations=10
```

If you are running EDS or require additional output from plumed, you can add additional plumed script by
```python
plumed_input0=textwrap.dedent(
    '''
    plumed script
    ''')
plumed_input=pte_plumed_script+plumed_input0
with open('plumed_eds_conver_pt_wte_metad.dat', 'w') as f:
    f.write(plumed_input)
ps.add_file('plumed_eds_conver_pt_wte_metad.dat')
```

Code for running simulation.
```python
for i in range(max_iterations):
    ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_conver_{}'.format(i),  
           mdp_kwargs=kwargs, run_kwargs={'plumed':'plumed.dat','replex':remd_exhcange_period}, mpi_np=MPI_NP)
    with open(ps.pickle_name, 'wb') as f:
        pickle.dump(ps, file=f)

    rep_eff_1 = get_replex_e(ps, replicas)
    if rep_eff_1 ==-1:
        ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_conver_{}'.format(i),  mdp_kwargs=kwargs,mpi_np=MPI_NP, run_kwargs={'plumed':'plumed.dat','replex':remd_exhcange_period},mpi_np=MPI_NP)
        with open(ps.pickle_name, 'wb') as f:
            pickle.dump(ps, file=f)
        
    elif (min(rep_eff_1) >= replex_eff and i>=min_iterations):
        print 'Reached replica exchange efficiency of {}. Continuing to production'.format(rep_eff_1)
        break

        with open(ps.pickle_name, 'wb') as f:
            pickle.dump(ps, file=f)
    else:
        print 'Replica exchange efficiency of {}. Continuing simulation'.format(rep_eff_1)


final_time_eds=int(0.040*5*10**5)
ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod',  mdp_kwargs={'nsteps': final_time_eds, 'ref_t': 278},mpi_np=MPI_NP)
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)
#finally:
print('done')
```

Write a bash script and run.

**Note**:MPI_NP must be consistant with Ntask in bash script. Restarts should also have the same number of MPI processes otherwise it will result in PTE tuning log file missing error.


## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate when making changes.
​

&copy; Andrew White at University of Rochester
