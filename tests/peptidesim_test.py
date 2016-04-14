from unittest import TestCase


import os.path, time

from peptidesim import *
import shutil, os, textwrap

class TestPeptideSimInit(TestCase):
    def setUp(self):
        self.p = PeptideSim('psim_test', ['AA', 'RE'], [3, 1], job_name='testing')
        
    def test_init(self):
        self.assertTrue(os.path.exists('psim_test'))

    def test_filechain(self):
        self.p.gro_file = 'a'
        self.p.gro_file = 'b'
        self.assertEqual(len(self.p._gro), 2)
        
    def test_logging_started(self):
        log_file = self.p.log_file
        self.assertTrue(os.path.exists(log_file))
        self.assertTrue(os.stat(log_file).st_size > 0)

    def test_packmol_success(self):
        self.p.initialize()    
        output_file = self.p.pdb_file
        self.assertIsNotNone(output_file)
        self.assertTrue(os.path.exists(output_file))

    def test_pdb2gmx_success(self):
        self.p.initialize()    
        self.assertTrue(os.path.exists(self.p.top_file))
        self.assertTrue(os.path.exists(self.p.gro_file))
        self.assertTrue(os.stat(self.p.gro_file).st_size > 0)

        
    def tearDown(self):
        shutil.rmtree('psim_test')

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

