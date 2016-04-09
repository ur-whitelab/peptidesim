from unittest import TestCase
import peptidesim
import os.path, time

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
