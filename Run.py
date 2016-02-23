from Simulation import Simulation
import sys
t=Simulation(*sys.argv[1:])
t.packmol()
t.initial_setup()
t.solvate_enem()
t.equilibrate()
t.nvt_grompp_mdrun()
