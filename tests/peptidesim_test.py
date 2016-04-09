from unittest import TestCase


import os.path, time

from peptidesim import *
import shutil, os, textwrap


class TestPeptideSim(TestCase):
    def test_sanity(self):
        self.assertTrue(True)

        
    def test_analysis_returns_images(self):

        prea_time = time.ctime()
        ret = p.analysis()
        
        self.assertTrue(len(ret) > 0)
        
        #check if the images exist
        oldest = time.ctime()
        for path in ret:
            self.assertTrue(os.path.exists(path))
            #get oldest
            oldest = min(oldest, time.ctime(os.path.getctime(path)))
        #check that the oldest image was created after analysis call
        self.assertTrue(prea_time < oldest)


    def equil_top_sanity(self):
        self.assertTrue(topol.exists())
        try:
            self.assertTrue(topol.read()>=0)
        finally:
            topol.close()
    def equil_addwater(self):
        self.assertTrue(peptide_addwater_gro.exists())
        try:
            self.assertTrue(peptide_addwater_gro.read()>=0)
        finally:
            topol.close()

    def equil_add_ions(self):
        self.assertTrue(peptide_addwater_addions_gro.exists())
        try:
            self.assertTrue(peptide_addwateraddions_gro.read()>=0)
        finally:
            peptide_addwateraddions_gro.close()



    def equil_added_ions(self):
        try:
            self.assertTrue(peptide_addwater_gro.read()!=peptide_addwater_addions_gro.read())
        finally:
            peptide_addwater_gro.close()
            peptide_addwateraddions_gro.close()

    def equil_mdp(self):
        self.assertTrue(mdp_equil.exists())
        try:
            self.assertTrue(mdp_equil.read()>=0)
        finally:
            mdp_equil.close()
            
    def equil_trr(self):
        self.assertTrue(trr.exists())
        try:
            self.assertTrue(trr.read()>=0)
        finally:
            trr.close()
    def eauil_tpr(self):
        self.assertTrue(tpr.exists())
        try:
            self.assertTrue(tpr.read()>=0)
        finally:
            tpr.close()
    def equil_success(self):
        self.assertTrue(equil_success_value)
   

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

