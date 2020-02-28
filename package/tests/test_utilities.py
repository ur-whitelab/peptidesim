import pandas as pd
import os
import sys
import unittest
from peptidesim import utilities


class Test(unittest.TestCase):
    '''These tests can identify if EDS chemical shifts are valid based on plumed data. The test is an example of how the chemical shift assertion works.
        '''

    def test_cs_shifts_example(self):
        txt_plumed = ('PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=250\n' +
                      'eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_10,cs.nh_14,cs.nh_15,cs.nh_16,cs.nh_22,cs.cb_1,cs.cb_11,cs.cb_16,cs.ha_1,cs.ha_10,cs.ca_1,cs.ca_2,cs.ca_5,cs.c_3,cs.c_5,cs.c_7 CENTER_ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_6,cs.exphn_7,cs.exphn_8,cs.exphn_9,cs.exphn_10,cs.exphn_13,cs.exphn_14,cs.exphn_15,cs.exphn_16,cs.exphn_18,cs.exphn_19,cs.exphn_20,cs.exphn_21,cs.exphn_22,cs.expnh_1,cs.expnh_2,cs.expnh_3,cs.expnh_4,cs.expnh_6,cs.expnh_7,cs.expnh_8,cs.expnh_9,cs.expnh_10,cs.expnh_13,cs.expnh_14,cs.expnh_15,cs.expnh_16,cs.expnh_18,cs.expnh_19,cs.expnh_20,cs.expnh_21,cs.expnh_22 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 PERIOD=250 TEMP=@replicas:{278.0,284.825679594,291.818948762,298.9839223,306.324816031,313.84594929,321.551747461,329.446744586,337.535586031,345.823031217,354.313956423,363.013357653,371.926353579,381.05818855,390.414235678,400.0} OUT_RESTART=eds_bias_restart_correct3.dat\n' +
                      'PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat')
        plumed_dat = 'Plumed.dat'
        with open(plumed_dat, 'w') as f1:
            f1.write(txt_plumed)
        keys = ['CAshifts', 'CBshifts', 'Cshifts',
                'HAshifts', 'Hshifts', 'Nshifts']
        df_ca = pd.DataFrame({'0': ['#1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '#12', '#13', '14', '15', '16', '17', '18',
                                    '19', '20', '21', '22', '23', '#24'], '1': [1.5, 11.5, 1, 1, 1.5, 1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})
        df_cb = pd.DataFrame({'0': ['#1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '#12', '#13', '14', '15', '16', '17', '18',
                                    '19', '20', '21', '22', '23', '#24'], '1': [1.5, 1, 1.5, 1, 1.5, 1, 1.5, 1, 1, 1, 1.2, 1, 1, 1, 1, 1.2, 1, 1, 1, 1, 1, 1, 1, 1]})
        df_c = pd.DataFrame({'0': ['#1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '#12', '#13', '14', '15', '16', '17', '18',
                                   '19', '20', '21', '22', '23', '#24'], '1': [1, 1, 1.5, 1.5, 1.2, 1, 1.6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})
        df_ha = pd.DataFrame({'0': ['#1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '#12', '#13', '14', '15', '16', '17', '18',
                                    '19', '20', '21', '22', '23', '#24'], '1': [1.5, 1.2, 1, 1, 1, 1, 1, 1, 1, 1.3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})
        df_h = pd.DataFrame({'0': ['#1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '#12', '#13', '14', '15', '16', '17', '18',
                                   '19', '20', '21', '22', '23', '#24'], '1': [1.5, 1.5, 8.51, 8.39, 1, 8.5, 8.46, 8.52, 8.62, 8.02, 8.44, 8.14, 1, 1, 8.51, 8.39, 8.25, 1, 8.46, 8.52, 8.62, 8.02, 8.44, 8.14]})
        df_n = pd.DataFrame({'0': ['#1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '#12', '#13', '14', '15', '16', '17', '18',
                                   '19', '20', '21', '22', '23', '#24'], '1': [1, 1, 120.7, 122.2, 123.8, 121.4, 122, 116.9, 111, 120.4, 122.7, 121.5, 1, 10, 120.7, 122.2, 123.8, 121.4, 122.2, 123.8, 1111, 120.4, 122.7, 121.5]})
        df = [df_ca, df_cb, df_c, df_ha, df_h, df_n]
        for i in range(len(keys)):
            shifts_dat = keys[i] + '.dat'
            df[i].to_csv(shifts_dat, index=False,
                         index_label=False, sep=' ', header=False)

        with open(plumed_dat, 'r') as file_plumed:
            shifts_plumed = utilities.parser(file_plumed)
        os.remove(plumed_dat)
        exec_dir = os.getcwd()
        dat_files = [f for f in os.listdir(exec_dir) if f.endswith('.dat')]
        shifts_plumed_dict = dict.fromkeys(keys)
        for i in range(len(shifts_plumed)):
            shifts_dat = exec_dir + '//' + dat_files[i]
            df_shifts_dat = pd.read_csv(shifts_dat, sep=' ', header=None)
            os.remove(shifts_dat)
            shifts_plumed_dict[keys[i]] = shifts_plumed[i]
            if shifts_plumed[i] != []:
                count, Invalid_index, Shifts_index = utilities.validity_check(
                    df_shifts_dat, shifts_plumed_dict[keys[i]])
                #print (df_shifts_dat.values[4][1],shifts_plumed_dict[keys[0]],count)
                self.assertTrue(count == 0, msg=(
                    'Simulation uses data {} for {}. Invalid EDS chemical shift found for data {}!'.format(Shifts_index, keys[i], Invalid_index)))


if __name__ == '__main__':
    unittest.main()
