'''
Configuration File:
-------------------

config_name: The name of the config file to look at for the sim parameters to use

See the template directory for predefined config files. Run ::

$ peptidesim --config <configuration file name>

to generate a config file in the current directory based on the config templates provided, or use
``default`` to generate the default configuration.
'''

#make it so we can write in Python3


from .version import __version__

from peptidesim.config import *

def main():
    generate_config()


if __name__ == '__main__':
    main()
else:
    from peptidesim.peptidesim import SimulationInfo, PeptideSim
