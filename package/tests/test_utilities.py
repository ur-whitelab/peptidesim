import pandas as pd
import os
import shutil
import sys
import unittest
import peptidesim
from peptidesim import utilities, PeptideSim
import textwrap


class CSTests(unittest.TestCase):
    '''These tests can identify if EDS chemical shifts are valid based on plumed data. The test is an example of how the chemical shift assertion works.
        '''

    def test_cs_data_prepare(self):
        p = PeptideSim('cs-test', ['AA', 'GGR'], [2, 1])
        p.run_kwargs = {'nt': 1}
        p.initialize()
        # Just make sure it runs
        utilities.prepare_cs_data(p, {'0-1-A-H': 4.4})
        # make sure it fails with bad h shift
        with self.assertRaises(ValueError):
            utilities.prepare_cs_data(p, {'0-1-A-H': 4.4, '4-3-A-CA': 44})
        # check with pte reweight
        utilities.prepare_cs_data(p, {'0-1-A-H': 4.4, '1-3-R-CA': 44}, True)
        # check without reweight
        utilities.prepare_cs_data(p, {'0-1-A-H': 4.4, '1-3-R-CA': 44}, False)
        shutil.rmtree('cs-test')


    def test_cs_shifts_example(self):
        txt_plumed = ('PRINT ARG=(cs\.cb-.*),(cs\.ca-.*),(cs\.hn-.*),(cs\.ha-.*),(cs\.nh-.*),(cs\.co-.*) FILE=WT_PTE_CS_shifts STRIDE=250\n' +
                      'eds: EDS ARG=cs.hn-1,cs.hn-2,cs.hn-3,cs.hn-10,cs.nh-14,cs.nh-15,cs.nh-16,cs.nh-22,cs.cb-1,cs.cb-11,cs.cb-16,cs.ha-1,cs.ha-10,cs.ca-1,cs.ca-2,cs.ca-5,cs.c-3,cs.c-5,cs.c-7 CENTER_ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_6,cs.exphn_7,cs.exphn_8,cs.exphn_9,cs.exphn_10,cs.exphn_13,cs.exphn_14,cs.exphn_15,cs.exphn_16,cs.exphn_18,cs.exphn_19,cs.exphn_20,cs.exphn_21,cs.exphn_22,cs.expnh_1,cs.expnh_2,cs.expnh_3,cs.expnh_4,cs.expnh_6,cs.expnh_7,cs.expnh_8,cs.expnh_9,cs.expnh_10,cs.expnh_13,cs.expnh_14,cs.expnh_15,cs.expnh_16,cs.expnh_18,cs.expnh_19,cs.expnh_20,cs.expnh_21,cs.expnh_22 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 PERIOD=250 TEMP=@replicas:{278.0,284.825679594,291.818948762,298.9839223,306.324816031,313.84594929,321.551747461,329.446744586,337.535586031,345.823031217,354.313956423,363.013357653,371.926353579,381.05818855,390.414235678,400.0} OUT_RESTART=eds_bias_restart_correct3.dat\n' +
                      'PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat')
        plumed_dat = 'Plumed.dat'
        with open(plumed_dat, 'w') as f1:
            f1.write(txt_plumed)
        keys = ['CAshifts', 'CBshifts', 'Cshifts',
                'HAshifts', 'Hshifts', 'Nshifts']
        df_ca = pd.DataFrame({'0': ['#1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '#12', '#13', '14', '15', '16', '17', '18',
                                    '19', '20', '21', '22', '23', '#24'], '1': [1.5, 11.5, 1, 1, 1.5, 1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})
        df_cb = pd.DataFrame({'0': ['#1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '#12', '#13', '14', '15', '16', '17', '18',
                                    '19', '20', '21', '22', '23', '#24'], '1': [1.5, 1, 1.5, 1, 1.5, 1, 1.5, 1, 1, 1, 1.5, 1, 1, 1, 1, 1.2, 1, 1, 1, 1, 1, 1, 1, 1]})
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
                if count != 0:
                    with self.assertRaises(AssertionError):
                        print ('Simulation uses data {} for {}. Invalid EDS chemical shift found for data {}!'.format(Shifts_index, keys[i], Invalid_index))
                        self.assertTrue(count == 0)


class EDSTests(unittest.TestCase):
    def test_load_eds(self):
        with open('eds_out', 'w') as f:
            f.write(textwrap.dedent('''
            #! FIELDS time cn0.mean_center cn0.mean_set cn0.mean_target cn0.mean_coupling cn0.mean_maxrange cn0.mean_maxgrad cn0.mean_accum cn0.mean_mean cn0.mean_pseudovirial cn0.mean_std cn1.mean_center cn1.mean_set cn1.mean_target cn1.mean_coupling cn1.mean_maxrange cn1.mean_maxgrad cn1.mean_accum cn1.mean_mean cn1.mean_pseudovirial cn1.mean_std cn2.mean_center cn2.mean_set cn2.mean_target cn2.mean_coupling cn2.mean_maxrange cn2.mean_maxgrad cn2.mean_accum cn2.mean_mean cn2.mean_pseudovirial cn2.mean_std cn3.mean_center cn3.mean_set cn3.mean_target cn3.mean_coupling cn3.mean_maxrange cn3.mean_maxgrad cn3.mean_accum cn3.mean_mean cn3.mean_pseudovirial cn3.mean_std
            #! SET adaptive  1
            #! SET update_period  25
            #! SET seed  0
            #! SET kbt  2.494339e+00
            2.625000e-01 2.881080e+00 0.000000e+00 0.000000e+00 0.000000e+00 6.235847e+02 6.235847e+02 0.000000e+00 3.346075e+00 0.000000e+00 1.724248e-06 8.228687e+00 0.000000e+00 0.000000e+00 0.000000e+00 6.235847e+02 6.235847e+02 0.000000e+00 9.035338e+00 0.000000e+00 1.089532e-05 2.365638e+01 0.000000e+00 0.000000e+00 0.000000e+00 6.235847e+02 6.235847e+02 0.000000e+00 2.452908e+01 0.000000e+00 6.367754e-05 6.856322e+01 0.000000e+00 0.000000e+00 0.000000e+00 6.235847e+02 6.235847e+02 0.000000e+00 6.707487e+01 0.000000e+00 3.344642e-04
            2.750000e-01 2.881080e+00 -6.235847e+02 0.000000e+00 -5.986413e+02 6.235847e+02 6.235847e+02 7.539987e+03 0.000000e+00 0.000000e+00 0.000000e+00 8.228687e+00 -6.235847e+02 0.000000e+00 -5.986413e+02 6.235847e+02 6.235847e+02 1.503917e+02 0.000000e+00 0.000000e+00 0.000000e+00 2.365638e+01 -6.235847e+02 0.000000e+00 -5.986413e+02 6.235847e+02 6.235847e+02 3.223598e+00 0.000000e+00 0.000000e+00 0.000000e+00 6.856322e+01 -6.235847e+02 0.000000e+00 -5.986413e+02 6.235847e+02 6.235847e+02 7.675593e-02 0.000000e+00 0.000000e+00 0.000000e+00
            2.880000e-01 2.881080e+00 -6.235847e+02 0.000000e+00 -6.235847e+02 6.235847e+02 6.235847e+02 7.539987e+03 3.522762e+00 -3.075063e-01 3.184181e-05 8.228687e+00 -6.235847e+02 0.000000e+00 -6.235847e+02 6.235847e+02 6.235847e+02 1.503917e+02 9.574345e+00 1.312518e-01 2.969051e-04 2.365638e+01 -6.235847e+02 0.000000e+00 -6.235847e+02 6.235847e+02 6.235847e+02 3.223598e+00 2.613460e+01 -5.923110e-01 2.822439e-03 6.856322e+01 -6.235847e+02 0.000000e+00 -6.235847e+02 6.235847e+02 6.235847e+02 7.675593e-02 7.176664e+01 -4.892661e+00 2.698455e-02
            3.005000e-01 2.881080e+00 -5.959157e+00 0.000000e+00 -3.066418e+01 6.235847e+02 6.235847e+02 3.963978e+05 0.000000e+00 -3.075063e-01 0.000000e+00 8.228687e+00 -1.205515e-01 0.000000e+00 -2.505912e+01 6.235847e+02 6.235847e+02 3.890082e+05 0.000000e+00 1.312518e-01 0.000000e+00 2.365638e+01 -2.584716e-03 0.000000e+00 -2.494587e+01 6.235847e+02 6.235847e+02 3.888610e+05 0.000000e+00 -5.923110e-01 0.000000e+00 6.856322e+01 -6.154410e-05 0.000000e+00 -2.494345e+01 6.235847e+02 6.235847e+02 3.888579e+05 0.000000e+00 -4.892661e+00 0.000000e+00
            '''))
        peptidesim.plot_couplings('eds_out', 'plot.png')
        os.remove('eds_out')
        os.remove('plot.png')


if __name__ == '__main__':
    unittest.main()
