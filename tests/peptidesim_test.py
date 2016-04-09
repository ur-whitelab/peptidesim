from unittest import TestCase
from peptidesim import *
import shutil, os, textwrap

class TestPeptideSim(TestCase):
    def test_sanity(self):
        self.assertTrue(True)


class TestConfig(TestCase):
    def test_config_setname(self):
        name = "test_config_setname.py"
        try:
            with open(name, "w") as f:
                f.write(textwrap.dedent('''
                    c = get_config()
                    c.PeptideSim.sim_name = 'NVT'
                '''))
            p = PeptideSim('testconfig', ['ALA'], [1], name)
            self.assertTrue(p.sim_name == 'NVT')
        finally:
            os.remove(name)
            try:
                shutil.rmtree('testconfig')
            except OSError as e:
                pass
            
            
