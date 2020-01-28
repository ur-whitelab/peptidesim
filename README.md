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

### Regular Simulation

1. Create a Python script simple.py and import PeptideSim
```python
from pepsidesim import PeptideSim
```
2. Specify input condition and initialize
```python
seq = 'EKEKEKEKEKEK' #input the sequence in one-letter-code
name = 'EK6' #naming the sequence
pep_copies = 1 #specifiy the number of copies
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
4. Run the script
```python
python simple.py
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate when making changes.
​

&copy; Andrew White at University of Rochester
