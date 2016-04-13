from unittest import TestCase


import os.path, time

from peptidesim import *
import shutil, os, textwrap

class EasyPeptideSimTests(TestCase):
    def test_tests_working(self):
        self.assertTrue(True)

    def test_init(self):
        p = PeptideSim('example', ['A'], [2])
        self.assertTrue(os.path.exists('example'))
        shutil.rmtree('example')



class TestPeptideSim(TestCase):
    def setUp(self):
        self.p = PeptideSim('psim_test', ['AA', 'RE'], [3, 1], job_name='testing')
    def test_logging_started(self):
        log_file = self.p.log_file
        print(log_file)
        self.assertTrue(os.path.exists(log_file))

    def test_packmol_success(self):
        output_file = self.p.pdb_file
        self.assertIsNotNone(output_file)
        self.assertTrue(os.path.exists(output_file))

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

