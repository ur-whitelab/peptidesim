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
module load anaconda3 packmol libmatheval gromacs-plumed/2018.3/b4 openblas git
```

2. Create a virtual environment
```bash
conda create -n yourenvname python=3.7
```
Peptidesim is no longer available for python2 versions. Once you create your environment, you should activate it for the next steps:
```bash
conda activate yourenvname
```
Check on pip and python and make sure they are sourced to your environment directory:
```bash
which pip #or which python
# The output should look like this:
# /.conda/envs/yourenvname/bin/pip
```

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

 Alternatively, on **a non-head node** use interactive python to do this:
```python
import gromacs
gromacs.config.setup()
```


# Developer Test Environment

## Creating Docker Image

Load the plumed gromacs docker image from dockerhub:

```sh
[sudo] docker pull whitelab/plumed-gromacs
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
docker run -it --rm -v "`pwd`:/home/whitelab/peptidesim" peptidesim/test
```

This will run all tests and clean-up.

## Running Unit Tests Interactively

If you want to leave all test files around and have python access to troubleshoot,
including an editable install so code you change is reflected, use:

```sh
docker run --rm -it -v `pwd`:/home/whitelab/peptidesim peptidesim/test bash ../interact.sh
python -m pytest -x ../peptidesim/package/tests/
```

After running these two lines, you will be in the docker container and the tests will be run. Use `exit`
to leave the container.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate when making changes.
​

&copy; Andrew White at University of Rochester