from unittest import TestCase, skip
import os.path
import time
from peptidesim import *
import shutil
import os
import textwrap
import dill as pickle
from io import BytesIO
import json
import requests
import signal


SIM_KWARGS = {'nt': 1}

class TestPeptideSimSimple(TestCase):
    def setUp(self):
        self.p = PeptideSim('psim_test', ['AA', 'RE'], [
                            3, 1], job_name='testing')
        self.p.run_kwargs = SIM_KWARGS

    def test_init(self):
        self.assertTrue(os.path.exists('psim_test'))
        self.assertEqual(self.p.run_kwargs, SIM_KWARGS)

    def test_filechain(self):
        start = len(self.p._gro)
        self.p.gro_file = 'a'
        self.p.gro_file = 'b'
        self.assertEqual(len(self.p._gro), 2 + start)

    def test_gro_file_list(self):
        self.p.gro_file = ['a', 'b']
        self.assertEqual(self.p.gro_file, 'a')
        self.assertEqual(self.p.gro_file_list[1], 'b')

    def test_logging_started(self):
        log_file = self.p.log_file
        self.assertTrue(os.path.exists(log_file))
        self.assertTrue(os.stat(log_file).st_size > 0)

    def test_req_files(self):
        with open('test.txt', 'w') as f:
            f.write('fdsa\n')
        self.p.add_file('test.txt')
        # make sure it's on required files list. Use string find since other path info will be present
        self.assertTrue(self.p._file_list[-1].find('test.txt') != -1,
                        'Adding required file did not put it on file_list. Found {}'.format(self.p._file_list[-1]))
        # make sure it's present now in directory
        self.p.initialize()
        self.assertTrue(os.path.exists('psim_test/prep/test.txt'))
        os.remove('test.txt')

    def tearDown(self):
       shutil.rmtree('psim_test')


class TestPeptideSimInitialize(TestCase):

    def setUp(self):
        self.p = PeptideSim('pinit_test', ['AA', 'REE'], [
                            3, 1], job_name='testing')
        self.p.run_kwargs = SIM_KWARGS
        self.p.initialize()

    def test_packmol_success(self):
        output_file = self.p.pdb_file
        self.assertIsNotNone(output_file)
        self.assertTrue(os.path.exists(output_file))

    def test_pdb2gmx_success(self):
        self.assertTrue(os.path.exists(self.p.top_file))
        self.assertTrue(os.path.exists(self.p.gro_file))
        self.assertTrue(os.stat(self.p.gro_file).st_size > 0)

    def test_pdbfile(self):
        self.assertTrue(os.path.exists(self.p.pdb_file))

    def test_include_resolution(self):
        if '#include' in open(self.p.top_file).read():
            self.fail(
                'There are unresolved #include directives in the topology file')

    def test_pickle(self):

        phash = self.p.top_file + self.p.gro_file + self.p.pdb_file
        string = pickle.dumps(self.p)
        # need to delete old object so we don't get duplicate logging
        del self.p
        new_p = pickle.load(BytesIO(string))
        self.assertEqual(phash, new_p.top_file +
                         new_p.gro_file + new_p.pdb_file)

    def test_ndx(self):
        # just see if the list item exists
        # later can write test to check for atom numbers
        self.assertIsNotNone(self.p.ndx['peptide_3'])

    def test_restore(self):
        '''
        Assert that we don't repeat the initialization simulations
        and that no errors occur
        '''
        self.p.run(mdpfile='peptidesim_emin.mdp', mdp_kwargs={'nsteps': 50})
        self.p.run(mdpfile='peptidesim_nvt.mdp', mdp_kwargs={'nsteps': 25})
        n = len(self.p.sims)
        string = pickle.dumps(self.p)
        del self.p
        new_p = pickle.load(BytesIO(string))
        # make sure on pickling we still have our history intact
        self.assertEqual(len(new_p.sims), n)
        # now actually ensure that we don't repeat the initialization
        new_p.initialize()
        self.assertEqual(len(new_p.sims), n)

    def tearDown(self):
        shutil.rmtree('pinit_test')


class TestFileTransfer(TestCase):

    def setUp(self):
        self.p = PeptideSim('file_trans', ['AA', 'REE'], [3, 1])
        self.p.run_kwargs = SIM_KWARGS
        # make some files to move around
        if(not os.path.exists('data')):
           os.makedirs('data')

        with open('data/file.txt', 'w') as f:
            f.write('test\n')

    def test_move_directory(self):
        self.p.add_file('data')
        self.p.initialize()
        self.assertTrue(os.path.exists(os.path.join(
            self.p.sims[-1].location, 'data/file.txt')))

    def tearDown(self):
        shutil.rmtree('file_trans')
        shutil.rmtree('data')


class TestPeptideRestart(TestCase):

    def test_can_save(self):
        p = PeptideSim('can-save-test', ['AA'], [2])
        p.marker = True
        p.save('foo')
        g = PeptideSim('can-save-test', ['AA'], [2])
        self.assertTrue(
            g.marker, msg='Failed to reload all attributes in restarting from pickle')

    def test_can_skip(self):
        p = PeptideSim('can-skip-test', ['AA'], [2])
        p.initialize()
        p.run(mdpfile='peptidesim_emin.mdp',
              tag='restart-test', mdp_kwargs={'nsteps': 25})
        g = PeptideSim('can-skip-test', ['AA'], [2])
        self.assertEqual(len(p.sims), len(g.sims))
        g.initialize()
        g.run(mdpfile='peptidesim_emin.mdp',
              tag='restart-test', mdp_kwargs={'nsteps': 25})
        self.assertEqual(len(p.sims), len(g.sims))

    def tearDown(self):
        shutil.rmtree('can-skip-test', ignore_errors=True)
        shutil.rmtree('can-save-test', ignore_errors=True)

class TestPeptideStability(TestCase):
    def test_dipeptides(self):
        for a in ['GG', 'VV', 'EE', 'SS', 'YY', 'SS', 'PP']:
            p = PeptideSim('dipeptide', [a], [1], job_name='dipeptide')
            p.run_kwargs = SIM_KWARGS
            p.peptide_density = 0.005
            p.initialize()
            p.run(mdpfile='peptidesim_emin.mdp', mdp_kwargs={'nsteps': 50})
            p.run(mdpfile='peptidesim_nvt.mdp', mdp_kwargs={'nsteps': 25})
            shutil.rmtree('dipeptide')


class TestPTE(TestCase):
    def test_pte(self):
        # run a pte to get plumed output
        p = PeptideSim('pte_test', ['AA'], [1], job_name='test-pte')
        print(p.run_kwargs)
        p.mpi_np = 1
        print(p.run_kwargs)
        p.mdrun_driver = 'gmx_mpi'
        p.peptide_density = 0.005
        print(p.run_kwargs)
        p.initialize()
        print(p.run_kwargs)
        p.run(mdpfile='peptidesim_emin.mdp', mdp_kwargs={'nsteps': 100})
        p.run(mdpfile='peptidesim_nvt.mdp', mdp_kwargs={'nsteps': 100})
        pte_result = ''
        with self.assertRaises(RuntimeError) as cm:
            p.pte_replica(mpi_np=2, max_tries=2, min_iters=1, mdp_kwargs={
                          'nsteps': 250}, replicas=2, hot=315, eff_threshold=0.99)
        # make the plumed input file
        pte_result = p.pte_replica(mpi_np=2, max_tries=5, mdp_kwargs={
                                   'nsteps': 400}, replicas=2, hot=315, eff_threshold=0.00)
        with open('plumed.dat', 'w') as f:
            f.write(pte_result['plumed'])
        p.add_file('plumed.dat')

        # now try running it with PTE
        p.run(tag='pte_check', mdpfile='peptidesim_nvt.mdp', mpi_np=2,
              mdp_kwargs=[{'nsteps': 100, 'ref_t': ti}
                          for ti in pte_result['temperatures']],
              run_kwargs={'plumed': 'plumed.dat', 'replex': 25})
        shutil.rmtree('pte_test')

        shutil.rmtree('pte_test')

    def test_pte_restart(self):
        # run a pte to get plumed output
        p = PeptideSim('pte_test_restart', ['AA'], [1], job_name='test-pte-restart')
        p.mpi_np = 1
        p.mdrun_driver = 'gmx_mpi'
        p.peptide_density = 0.005
        p.initialize()
        p.run(mdpfile='peptidesim_emin.mdp', mdp_kwargs={'nsteps': 100})
        p.run(mdpfile='peptidesim_nvt.mdp', mdp_kwargs={'nsteps': 100})
        #test pickle on signal
        signal.alarm(1)
        try:
            p.pte_replica(mpi_np=2, tag='pte_tune_test', max_tries=5, mdp_kwargs={
                          'nsteps': 250}, replicas=2, hot=315, eff_threshold=0.01, min_iters=1, dump_signal=signal.SIGALRM)
        except KeyboardInterrupt:
            pass

        del p

        new_p = PeptideSim('pte_test_restart', ['AA'], [1], job_name='test-pte-restart')

        # make sure there is one simulation in history with pte

        self.assertTrue(new_p.sims[-1].short_name.startswith('pte_tune_test'))

        # try to restart it
        new_p.pte_replica(mpi_np=2, tag='pte_tune_test', max_tries=5, mdp_kwargs={
                          'nsteps': 250}, replicas=2, hot=315, min_iters=1, eff_threshold=0.01, dump_signal=signal.SIGALRM)
        shutil.rmtree('pte_test')

        shutil.rmtree('pte_test_restart')

class TestRemoveSimulation(TestCase):
    def test_remove(self):
        # run a pte to get plumed output
        p = PeptideSim('test_remove', ['AA'], [1], job_name='remove')
        p.run_kwargs = SIM_KWARGS
        p.peptide_density = 0.005
        p.initialize()
        p.run(mdpfile='peptidesim_emin.mdp', mdp_kwargs={'nsteps': 100})
        p.run(tag='nvt_check', mdpfile='peptidesim_nvt.mdp',
              mdp_kwargs={'nsteps': 100})
        old_gro_files_number = len(p._gro)
        old_tpr_files_number = len(p._tpr)
        old_sim_files_number = len(p._sims)

        # now try running it with PTE
        p.remove_simulation('nvt_check')
        new_gro_len = len(p._gro)
        new_tpr_len = len(p._tpr)
        new_sim_len = len(p._sims)
        self.assertGreaterEqual(old_gro_files_number, new_gro_len)
        self.assertGreaterEqual(old_tpr_files_number, new_tpr_len)
        self.assertGreaterEqual(old_sim_files_number, new_sim_len)

    def test_remove_restart(self):
        # run a pte to get plumed output
        p = PeptideSim('test_remove_restart', ['AA'], [1], job_name='test-remove')
        p.run_kwargs = SIM_KWARGS
        p.peptide_density = 0.005
        p.initialize()
        p.run(tag='eminiiii', mdpfile='peptidesim_emin.mdp',
              mdp_kwargs={'nsteps': 100})
        p.run(mdpfile='peptidesim_nvt.mdp', mdp_kwargs={'nsteps': 100})
        # test pickle on signal

        # make sure there is one simulation

        # try to restart it
        p.remove_simulation('eminiiii')
        with self.assertRaises(ValueError) as cm:
            p.remove_simulation('wrong_sim_name')
        with self.assertRaises(TypeError) as cm:
            p.remove_simulation(None)

class TestPeptideEmin(TestCase):
    def setUp(self):
        self.p = PeptideSim('pemin_test', ['VV'], [1], job_name='testing-emin')
        self.p.run_kwargs = SIM_KWARGS
        self.p.peptide_density = 0.005
        self.p.initialize()
        # do short emin
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='set-up-emin-constraints', mdp_kwargs={'nsteps': 25})

    def test_short_emin(self):
        start_gro = self.p.gro_file
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='short-test', mdp_kwargs={'nsteps': 10})
        self.assertTrue(start_gro != self.p.gro_file)

    def test_mpinp_emin(self):
        start_gro = self.p.gro_file
        self.p.mpi_np = 1
        self.p.run_kwargs = {}
        self.p.mdrun_driver = 'gmx_mpi'
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='short-test', mdp_kwargs={'nsteps': 10})
        self.assertTrue(start_gro != self.p.gro_file)

    def test_notag(self):
        start_gro = self.p.gro_file
        self.p.run(mdpfile='peptidesim_emin.mdp', mdp_kwargs={'nsteps': 1})
        self.assertTrue(start_gro != self.p.gro_file)

    def test_emin_mdp_combine(self):
        self.p.run(mdpfile='peptidesim_nvt.mdp',
                   tag='test-mdp',  mdp_kwargs={'nsteps': 2})
        self.assertIn('constraints', self.p.sims[-1].metadata['mdp-data'])

    def test_emin_mdp_kwargs(self):
        self.p.run(mdpfile='peptidesim_nvt.mdp',
                   tag='test-mdp',  mdp_kwargs={'nsteps': 7})
        self.assertEqual(
            str(self.p.sims[-1].metadata['mdp-data']['nsteps']), '7')

    def test_emin_metadata(self):
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='test', mdp_kwargs={'nsteps': 10})
        self.assertTrue('md-log' in self.p.sims[-1].metadata)

    def test_emin_metadata_multiple(self):
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='test2', mdp_kwargs={'nsteps': 10})
        self.assertTrue('md-log' in self.p.sims[-1].metadata)
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='test1', mdp_kwargs={'nsteps': 10})
        self.assertTrue('md-log' in self.p.sims[0].metadata)
        self.assertTrue('md-log' in self.p.sims[1].metadata)

    def test_emin_metadata_frompickle(self):
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='test1', mdp_kwargs={'nsteps': 10})
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='timeout', mdp_kwargs={'nsteps': 10})
        string = pickle.dumps(self.p)
        # need to delete old object so we don't get duplicate logging
        del self.p
        new_p = pickle.load(BytesIO(string))
        string = pickle.dumps(new_p)
        del new_p
        new_p = pickle.load(BytesIO(string))
        # I do it twice because I'm paranoid

        self.assertTrue('md-log' in new_p.sims[0].metadata)
        self.assertTrue('md-log' in new_p.sims[1].metadata)

    def test_pickle_emin(self):
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='test', mdp_kwargs={'nsteps': 10})

        string = pickle.dumps(self.p)
        # need to delete old object so we don't get duplicate logging
        del self.p
        new_p = pickle.load(BytesIO(string))
        string = pickle.dumps(new_p)
        del new_p
        new_p = pickle.load(BytesIO(string))

        self.assertTrue(len(new_p.sims) > 0)

    def test_restart_emin(self):

        # call and interrupt the function
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='timeout', mdp_kwargs={'nsteps': 10})
        for k, v in self.p._sims.items():
            if(k.startswith('emin-timeout')):
                self.assertTrue(v.restart_count == 2)

        # check metadata is intact
        self.assertTrue('md-log' in self.p.sims[-1].metadata)

    def test_continue_emin(self):
        # call and interrupt the function
        self.p.run(mdpfile='peptidesim_emin.mdp',
                   tag='repeat', mdp_kwargs={'nsteps': 10})
        self.p.run(mdpfile='peptidesim_emin.mdp', tag='repeat-1',
                   mdp_kwargs={'nsteps': 10}, repeat=True)
        # check locations are the same
        self.assertEqual(self.p.sims[-1].location, self.p.sims[-2].location)

    def test_signal_restart_emin(self):
        '''
        Test that if the simulation is killed exeternally, it can still be pickled and recovered
        '''
        # test pickle on signal
        signal.alarm(1)
        try:
            self.p.run(mdpfile='peptidesim_emin.mdp', tag='timeout-signal',
                       mdp_kwargs={'nsteps': 2500}, dump_signal=signal.SIGALRM)
        except KeyboardInterrupt:
            pass

        del self.p
        # rely on automatic restarting from pickle
        new_p = PeptideSim('pemin_test', ['VV'], [1], job_name='testing-emin')

        for k, v in new_p._sims.items():
            if(k.startswith('emin-timeout')):
                self.assertTrue(v.restart_count == 2)

        # check metadata is intact
        self.assertTrue('md-log' in new_p.sims[-1].metadata)

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
                shutil.rmtree('testconfig')
            except OSError as e:
                pass


class TestRestartPlumed(TestCase):
    def test_plumed_restart(self):
        # run a pte to get plumed output
        p = PeptideSim('plumed_test', ['AA'], [1], job_name='test-plumed')
        p.run_kwargs = SIM_KWARGS
        p.peptide_density = 0.005
        p.initialize()
        p.run(mdpfile='peptidesim_emin.mdp', mdp_kwargs={'nsteps': 50})
        p.run(mdpfile='peptidesim_nvt.mdp', mdp_kwargs={'nsteps': 50})
        # the name of the plumed file should always start with a plumed string
        # and end with .dat
        plumed_test_name = 'plumed_test.dat'

        # plumed input file contents
        text_old = textwrap.dedent('''
        DISTANCE ATOMS=1,2 LABEL=d1
        PRINT ARG=d1 FILE=colvar STRIDE=10
        ''')

        with open(plumed_test_name, 'w') as f0:
            f0.write(text_old)

        p.add_file(plumed_test_name)
        # the new content of the plumed file
        signal.alarm(1)
        try:
            p.run(
                mdpfile='peptidesim_nvt.mdp', tag='nvt', mdp_kwargs={
                    'nsteps': 250}, run_kwargs={
                        'plumed': plumed_test_name}, dump_signal=signal.SIGALRM)
        except KeyboardInterrupt:
            pass

        del p

        with open('test-plumed.pickle', 'r+b') as f:
            new_p = pickle.load(f)

        # now we modify the content of the same plumed file and add it to the
        # list of required files
        text_new = textwrap.dedent('''
        DISTANCE ATOMS=1,3 LABEL=d1
        PRINT ARG=d1 FILE=colvar STRIDE=10
        ''')
        # test pickle on signal
        f1 = open(plumed_test_name, 'w')
        f1.write(text_new)
        f1.close()
        new_p.add_file(plumed_test_name)
        signal.alarm(1)
        try:
            new_p.run(
                mdpfile='peptidesim_nvt.mdp', tag='nvt', mdp_kwargs={
                    'nsteps': 250}, run_kwargs={
                    'plumed': plumed_test_name}, dump_signal=signal.SIGALRM)
        except KeyboardInterrupt:
            pass

        reloaded_plumed = open(
            './' + new_p.sims[-1].location + '/' + plumed_test_name, 'r')
        f2 = open('./' + plumed_test_name, 'r')
        self.assertEqual(reloaded_plumed.readlines(), f2.readlines())
        f2.close()
        reloaded_plumed.close()
        os.remove('test-plumed.pickle')
        os.remove(plumed_test_name)
        shutil.rmtree('plumed_test')


if __name__ == '__main__':
    import unittest
    unittest.main()
