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
module load anaconda3 packmol libmatheval gromacs-plumed/2019.4/b2 openblas git
```

2. Create a virtual environment
```bash
conda create -n yourenvname python=3.7
```
Peptidesim is no longer available for python2 versions. Once you create your environment, you should activate it for the next steps:
```bash
source activate yourenvname
```
Check on pip and python and make sure they are sourced to your environment directory:
```bash
which pip #or which python
# The output should look like this:
# /.conda/envs/yourenvname/bin/pip
```
**Note**: In case the environment is not properly activated, i.e., `which pip` outputs a path not as mentioned above, try `conda init`. Close and reopen a new terminal and do `conda activate <yourenvname>`.

3. You will need to install the GromacsWrapper module separately due to a pending bug fix.
```bash
pip install --no-cache-dir git+https://github.com/whitead/GromacsWrapper
```

4. For installing PeptideSim, you need to clone the package from github.
```bash
git clone https://github.com/ur-whitelab/peptidesim.git
```
Change directory to `package` and install the requirements and the module:
```bash
cd package
pip install -e .
```

5. Setup a config file for Gromacswrapper to be able to find the gromacs installation. You can use the one given below and save it to `~/.gromacswrapper.cfg`


`.gromacswrapper.cfg`
```text
 [ DEFAULT ]
 qscriptdir = %( configdir )s/qscripts
 templatesdir = %( configdir )s/templates
 configdir = ~/.gromacswrapper
​
 [ Gromacs ]
 release = 2019.4
 gmxrc = /software/gromacs-plumed/2019.4/b2/bin/GMXRC
 extra =
 tools = gmx gmx_mpi
 groups = tools
​
 [ Logging ]
 logfilename = gromacs .log
 loglevel_console = INFO
 loglevel_file = DEBUG
```

 Alternatively, on **a non-head node** use interactive python to do this:
```python
import gromacs
gromacs.config.setup()
```


# Developer Test Environment

## Creating Docker Image

Load the plumed gromacs docker image from dockerhub:

```sh
docker pull whitelab/plumed-gromacs
```

Now we build the peptidesim testing image

```sh
cd test-docker
docker build -t peptidesim/test .
```

These two steps gather the plumed and gromacs version. Generally,
you do not need to re-run them.

## Running Unit Tests

From the repo root directory:

```sh
docker run -it --rm -v [path_to_peptidesim_root]:/home/whitelab/peptidesim peptidesim/test
```

This will run all tests and clean-up.

## Running Unit Tests Interactively

If you want to leave all test files around and have python access to troubleshoot,
including an editable install so code you change is reflected, use:

```sh
docker run --rm -it -v [path_to_peptidesim_root]:/home/whitelab/peptidesim peptidesim/test bash ../interact.sh
python -m pytest -x ../peptidesim/package/tests/
```

Prepare the docker image in the test-docker folder by running
the build script. You may need to have `sudo` depending on your docker configuration.

```bash
cd test-docker && [sudo] ./build.sh
```
It is only necessary to rebuild the docker script when newer gromacs,
gromacswrapper packages are available.


#### Running Tests with Docker

From the repo root directory:

```bash
[sudo] ./test.sh
```

 You may need to have `sudo` depending on your docker configuration.

 #### Interacting with Environment

 To have access to the docker environment with the built package
 so you can debug/develop, execute this commend from the repo
 root directory:

 ```bash
[sudo] ./interact.sh
 ```

You may need to have `sudo` depending on your docker configuration.
Type `exit` to leave the docker environment. See instructions that are printed after
running the command for how to interact/use the environment.

## Example Workflows

### NVT simulations with multiple peptides
A complete example can be found in `peptidesim/inputs/simple`.

#### 1. Imports
Create a Python script simple.py and import the PeptideSim class
```python
from pepsidesim import PeptideSim
```
#### 2. Input Conditions
Specify peptides, conditions, and initialize. Note
that all properties of the PeptideSim have defaults,
so you do not necessarily need to specify a concentration
or water model.

```python
seq = 'EKEKEKEKEKEK' # input the sequence in one-letter-code
name = 'EK6' # naming the sequence
pep_copies = 1 # specify the number of copies, increase the number if running multiple chains to study self assembly and etc.
ps = PeptideSim(name, [seq], [pep_copies], job_name='2mer_{}'.format(name)) # input to PeptideSim
ps.peptide_density = 0.008 # g/mol
ps.water = 'spce' # specify water model
ps.mpi_np = 4 # Number of MPI processes to use
ps.initialize()
```

#### 3. Simulation Steps
Here we do energy minimization, annealing and NVT equilibration. Note that
we can pass specific gromacs `mdp` arguments as python objects. The `tag` is
the short name of the simulation we are doing. The `mdpfile` names used here
are templates found in the peptidesim package. You can provide your own files
if desired.
```python
ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**5})
ps.run(mdpfile='peptidesim_anneal.mdp', tag='anneal_nvt')
ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod', mdp_kwargs={'nsteps': int(3 * 5*10**5), 'constraints': 'h-bonds'})
```

#### 4. Run the script

If you are using a slurm job script, you
could use this example:
```bash
#!/bin/sh
#SBATCH --partition=standard --time=120:00:00 --output=EK6.txt
#SBATCH --mail-type ALL
#SBATCH -N 2

python simple.py
```

#### 5. Saving at intermediate steps

Simulations may takes a long time depending on the number of steps chosen and resource requested. It is a good idea
to periodically save your python `PeptideSim` object so that you can restart. First, you'll want to check
in your script if a `pickle` file is found and restart rather than creating a new one.
```python
pickle_name = ...
if(os.path.exists(pickle_name)):
    print('loading restart')
    with open(pickle_name, 'rb') as f:
        ps = pickle.load(f)
else:
    # ...code from step 2
```

It is recommended to save the pickle file after each simulation:
```python
ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 10**5})
print('Pickled filename': ps.pickle_name)
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

ps.run(mdpfile='peptidesim_anneal.mdp', tag='anneal_nvt')
print('Pickled filename': ps.pickle_name)
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod', mdp_kwargs={'nsteps': int(3 * 5*10**5), 'constraints': 'h-bonds'})
print('Pickled filename': ps.pickle_name)
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)

```
Note that for python3 pickle files should always be opened in binary mode.

### Enhanced Sampling (PT-WTE)/ Experiment Directed simulation (EDS)

A complete example can be found in `peptidesim/inputs/pte`.

#### 1. Preparation
Create a Python script for the preparation step. It is similiar to NVT simulation script. 

If you are using experimental bias, additional code will be needed for the preparation. Check part1.py, [plumed METAD documentation](https://www.plumed.org/doc-v2.5/user-doc/html/_m_e_t_a_d.html) and [plumed EDS documentation](https://www.plumed.org/doc-v2.5/user-doc/html/_e_d_s.html).

#### 2. PT-WTE

Create a second script and import all required modules.
Specify peptides, replica number and exchange period. Note that the name of the job must be consistant with that in part 1.

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
    print('loading restart')
    with open(pickle_name, 'rb') as f:
        ps = pickle.load(f)
        #ps.pickle_name=pickle_name
        print(os.getcwd())
        ps.rel_dir_name='.'
```

Specify conditions for replica exchange. The `ps.pte_replica` function automatically generates plumed scripts based on the arguments inputted and run NVT tuning with PTE. Note the function has default input, you don't have to specify all of the inputs.

```python
pte_result = ps.pte_replica(tag='pte_tune', mpi_np=None, replicas=8, max_tries=30, min_iters=4, mdp_kwargs=dict(), run_kwargs=dict(), hills_file_location=None,cold=300.0, hot=400.0, eff_threshold=0.3,hill_height=1.2, sigma=140.0, bias_factor=10,exchange_period=25,dump_signal=signal.SIGTERM)
```


Proceed to NVT_production step after replica exchange.
```python
plumed_input0=textwrap.dedent(
    '''
    plumed script
    ''')

with open('plumed_prod', 'wb') as f:
    f.write(plumed_input0)
ps.add_file('plumed_prod.dat')
final_time_eds=int(0.040*5*10**5)
ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod',  mdp_kwargs={'nsteps': final_time_eds, 'ref_t': 278},run_kwargs={'plumed':'plumed_prod.dat'},mpi_np=MPI_NP)
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)
```

Write a bash script and run.

**Note**: During replica exchange simulations, MPI_NP must be consistant with N (number of nodes), ntask_per_node and c (cpus_per_task) in bash script. Restarts should also have the same number of MPI processes otherwise it will result in PTE tuning log file missing error.


## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate when making changes.
​

&copy; Andrew White at University of Rochester
