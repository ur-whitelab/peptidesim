This is a script to run a GROMACS 5.0.6 simulation a given a sequence of amino acids. Run.py is the file that actually calls the commands defined in a file called Simulation.py which has all of the necessary functions and steps needed to perform a gromacs simulation. Run.py script takes in a list of sequences of Amino acids and converts it to a pdb file using Peptide Builder. It then combines the pdb files into one using Packmol. It then calls gromacs 5 commands using GromacsWrapper. The commands will generate 4 folders: init_setup(where Packmol, peptidebuilder, generating box, solvate are called), energy_min_solvate(where energy minimization and solvation, and addition of ions occur), equilibration(equilibrates the system) and nvt(runs the simulation in nvt ensemble). The command that runs this script is the following:
'python Run.py folder_name sequence1 sequence2 sequence3 sequence_etc '
To run the script one needs:
gromacs-5.0.x, packmol, PeptideBuilder, gromacs_wrapper be installed and be on the path of the computer.
Installing and downloading instructions are in given below:
Peptide Builder: 'https://github.com/mtien/PeptideBuilder'
Packmol:'http://www.ime.unicamp.br/~martinez/packmol/'
GromacsWrapper:'https://pypi.python.org/pypi/GromacsWrapper'
Gromacs:'http://www.gromacs.org/Documentation/Installation_Instructions_5.0'
