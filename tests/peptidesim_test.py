from unittest import TestCase


import os.path, time

from peptidesim import *
import shutil, os, textwrap

import json, requests
import signal
import functools

class TimeoutError(Exception): pass

def timeout(seconds, error_message = 'Function call timed out'):
    def decorated(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return functools.wraps(func)(wrapper)

    return decorated

class TestPeptideSimSimple(TestCase):
    def setUp(self):
        self.p = PeptideSim('psim_test', ['AA', 'RE'], [3, 1], job_name='testing')
        
    def test_init(self):
        self.assertTrue(os.path.exists('psim_test'))

    def test_filechain(self):
        start = len(self.p._gro)
        self.p.gro_file = 'a'
        self.p.gro_file = 'b'
        self.assertEqual(len(self.p._gro), 2 + start)  

    def test_logging_started(self):
        log_file = self.p.log_file
        self.assertTrue(os.path.exists(log_file))
        self.assertTrue(os.stat(log_file).st_size > 0)
    
    def test_req_files(self):
        with open('test.txt', 'w') as f:
            f.write('fdsa\n')
        self.p.add_file('test.txt')
        #make sure it's on required files list. Use string find since other path info will be present
        self.assertTrue(self.p._file_list[-1].find('test.txt') != -1, 'Adding required file did not put it on file_list. Found {}'.format(self.p._file_list[-1]))
        #make sure it's present now in directory
        self.p.initialize()
        self.assertTrue(os.path.exists('psim_test/prep/test.txt'))
        os.remove('test.txt')

    def tearDown(self):
       shutil.rmtree('psim_test')
    


class TestPeptideSimInitialize(TestCase):

    def setUp(self):
        self.p = PeptideSim('pinit_test', ['AA', 'REE'], [3, 1], job_name='testing')
        self.p.initialize()
    
    def test_packmol_success(self):
        output_file = self.p.pdb_file
        self.assertIsNotNone(output_file)
        self.assertTrue(os.path.exists(output_file))

    def test_pdb2gmx_success(self):
        self.assertTrue(os.path.exists(self.p.top_file))
        self.assertTrue(os.path.exists(self.p.gro_file))
        self.assertTrue(os.stat(self.p.gro_file).st_size > 0)

    def test_include_resolution(self):
        if '#include' in open(self.p.top_file).read():
            self.fail('There are unresolved #include directives in the topology file')

    def test_pickle(self):
        import dill as pickle
        from cStringIO import StringIO
        phash = self.p.top_file + self.p.gro_file +  self.p.pdb_file
        string = pickle.dumps(self.p)
        #need to delete old object so we don't get duplicate logging
        del self.p
        new_p = pickle.load(StringIO(string))        
        self.assertEqual(phash, new_p.top_file + new_p.gro_file +  new_p.pdb_file)

    def test_data_stored(self):
        #tests if the data was correctly written to the json file
        self.p.store_data()
        self.assertTrue(os.path.exists('data/simdata.json'))
        with open('data/simdata.json', 'r') as f:
            data = json.load(f)
            self.assertEqual(data['sim_name'], self.p.sim_name)

    def test_to_database(self):
        #verifies that data can be properly sent to redis database
        with open('data/simdata.json', 'r') as f:
            data = json.load(f)
            prop = 'peptide_density'
            url = 'http://52.71.14.39/insert/simulation'
            payload = {'sim_name':data['sim_name'], 'property':prop, 'property_value': data[prop]}
            r = requests.put(url, payload)
            self.assertEqual(r.status_code, 200)
            
            



    def test_ndx(self):
        #just see if the list item exists
        #later can write test to check for atom numbers
        self.assertIsNotNone(self.p.ndx['peptide_3'])


    def test_neutral(self):
        pass
        

    def tearDown(self):
        shutil.rmtree('pinit_test')
        

class TestPeptideEmin(TestCase):
    def setUp(self):
        self.p = PeptideSim('pemin_test', ['RE', 'DC'], [1, 1], job_name='testing-emin')
        self.p.initialize()

    def test_short_emin(self):
        start_gro = self.p.gro_file
        self.p.run(mdpfile='peptidesim_emin.mdp', tag='test', mdp_kwargs={'steps':10})
        self.assertTrue(start_gro != self.p.gro_file)

    def test_restart_emin(self):
        
        start_gro = self.p.gro_file

        #call and interrupt the function
        @timeout(3,'')
        def wrap(this=self):
            self.p.run(mdpfile='peptidesim_emin.mdp', tag='timeout', mdp_kwargs={'steps':10**10})

        try:
            wrap()
        except TimeoutError:
            pass

        #attempt restart
        try:
            wrap()
        except TimeoutError:
            pass

        for k,v in self.p._sims.iteritems():
            if(k.startswith('emin-timeout')):
                self.assertTrue(v.restart_count == 2)

    def tearDown(self):
        shutil.rmtree('pemin_test')



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
                pass
                #shutil.rmtree('testconfig')
            except OSError as e:
                pass

