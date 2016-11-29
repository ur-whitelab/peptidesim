from traitlets.config import Configurable, Application, PyFileConfigLoader
from traitlets import Int, Float, Unicode, Bool, List, Instance, Dict
from .peptidesim import PeptideSim
from .version import __version__

import os


class PeptideSimConfigurator(Application):
    #information for running peptidesim from command line
    name = 'peptidesim'
    version = __version__
    classes = List([PeptideSim])
    description = '''The peptidesim application.

    Currently, this application only generates a configuration
    '''
    

    #command line flags
    aliases = Dict({'config:':'PeptideSimConfigurator.config_file'})
    
    config_file = Unicode('peptidesim_config.py',
        help="The config file to load",
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

        config_text = PeptideSim.class_config_section()
        if isinstance(config_text, bytes):
            config_text = config_text.decode('utf8')
        print("Writing default config to: %s" % self.config_file)
        with open(self.config_file, mode='w') as f:
            f.write(config_text)
            



def generate_config():
    p = PeptideSimConfigurator.instance()
    p.initialize()
    p.write_config_file()
    
