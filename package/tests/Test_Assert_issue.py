import numpy as np
import pandas as pd
import os,sys
import unittest


def validity_check(df, shifts):
            Invalid_index = []
            Shifts_index = []
            for i in range(len(shifts)):
                count = 0
                Index = shifts[i]
                Shifts_index.append(Index+1)
                if df.values[Index][1] ==1:
                    Invalid_index.append(Index+1)
                    count +=1
            return (count,Invalid_index,Shifts_index)

def parser(plumed_file):
    ## Finding plumed data
    string = 'eds: EDS ARG=' 
    for line in plumed_file: 
        if string in line: 
            txt = line
    ##skipping 4 letters to account for "ARG="
    txt = txt.split(' ')
    txt = txt[2][4:]
    txt_list = txt.split(',')
    hn_shifts = []
    nh_shifts = []
    ca_shifts = []
    cb_shifts = []
    for i in range(len(txt_list)):
        cs_type = txt_list[i][0:5]
        if cs_type == "cs.hn":
            hn_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == "cs.nh":
            nh_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == "cs.ca":
            ca_shifts.append(int(txt_list[i][6:])-1)
        elif cs_type == "cs.cb":
            cb_shifts.append(int(txt_list[i][6:])-1)
    return (hn_shifts,nh_shifts,ca_shifts,cb_shifts)

class Test(unittest.TestCase):
    """These tests can identify if EDS chemical shifts are valid based on plumed data. test_parser verifies that the parser fuction works.
       Other tests are designed to evaluate how validity_check fuction works with different chemical shifts.
        """
    def test_parser(self):
        txt_plumed = ("PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=250\n"+                     
        "eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_10,cs.hn_13,cs.hn_19,cs.hn_20,cs.hn_21,cs.nh_14,cs.nh_15,cs.nh_16,cs.nh_22,cs.cb_11,cs.cb_16 CENTER_ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_6,cs.exphn_7,cs.exphn_8,cs.exphn_9,cs.exphn_10,cs.exphn_13,cs.exphn_14,cs.exphn_15,cs.exphn_16,cs.exphn_18,cs.exphn_19,cs.exphn_20,cs.exphn_21,cs.exphn_22,cs.expnh_1,cs.expnh_2,cs.expnh_3,cs.expnh_4,cs.expnh_6,cs.expnh_7,cs.expnh_8,cs.expnh_9,cs.expnh_10,cs.expnh_13,cs.expnh_14,cs.expnh_15,cs.expnh_16,cs.expnh_18,cs.expnh_19,cs.expnh_20,cs.expnh_21,cs.expnh_22 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 PERIOD=250 TEMP=@replicas:{278.0,284.825679594,291.818948762,298.9839223,306.324816031,313.84594929,321.551747461,329.446744586,337.535586031,345.823031217,354.313956423,363.013357653,371.926353579,381.05818855,390.414235678,400.0} OUT_RESTART=eds_bias_restart_correct3.dat\n"+
        "PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat")
        plumed_dat = "Plumed.dat"
        with open(plumed_dat,"w") as f1:
            f1.write(txt_plumed)
        f1.close()
        file_plumed = open(plumed_dat, 'r')
        try: 
            parsed_inputs = parser(file_plumed)
        finally:
            file_plumed.close()
        True_parsed_inputs = ([0, 1, 2, 9, 12, 18, 19, 20], [13, 14, 15, 21],[], [10, 15])
        
        os.remove(plumed_dat)
        
        self.assertEqual(parsed_inputs,True_parsed_inputs)

    def test_EDS_cb(self):
        df_cb = pd.DataFrame({'0':['#1','2','3','4','5','6','7','8','9','10','11','#12','#13','14','15','16','17','18','19','20','21','22','23','#24'],'1':[1.5,2,1.5,2,1,1,1.5,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]})
        shifts_dat = "Shifts.dat"
        df_cb.to_csv(shifts_dat,index=False,index_label=False)
        df = pd.read_csv(shifts_dat)
        txt_plumed = ("PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=250\n"+                     
        "eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_10,cs.hn_13,cs.hn_19,cs.hn_20,cs.hn_21,cs.nh_14,cs.nh_15,cs.nh_16,cs.nh_22,,cs.cb_1,cs.cb_11,cs.cb_16 CENTER_ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_6,cs.exphn_7,cs.exphn_8,cs.exphn_9,cs.exphn_10,cs.exphn_13,cs.exphn_14,cs.exphn_15,cs.exphn_16,cs.exphn_18,cs.exphn_19,cs.exphn_20,cs.exphn_21,cs.exphn_22,cs.expnh_1,cs.expnh_2,cs.expnh_3,cs.expnh_4,cs.expnh_6,cs.expnh_7,cs.expnh_8,cs.expnh_9,cs.expnh_10,cs.expnh_13,cs.expnh_14,cs.expnh_15,cs.expnh_16,cs.expnh_18,cs.expnh_19,cs.expnh_20,cs.expnh_21,cs.expnh_22 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 PERIOD=250 TEMP=@replicas:{278.0,284.825679594,291.818948762,298.9839223,306.324816031,313.84594929,321.551747461,329.446744586,337.535586031,345.823031217,354.313956423,363.013357653,371.926353579,381.05818855,390.414235678,400.0} OUT_RESTART=eds_bias_restart_correct3.dat\n"+
        "PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat")
        plumed_dat = "Plumed.dat"
        with open(plumed_dat,"w") as f1:
            f1.write(txt_plumed)
        file_plumed = open(plumed_dat, 'r')
        hn_shifts,nh_shifts,ca_shifts,cb_shifts = parser(file_plumed)
        file_plumed.close()
        if cb_shifts ==[]:
                raise ValueError('Desired chemical shift is unavailable in the provided plumed data!')
        os.remove(shifts_dat)
        count, Invalid_index, Shifts_index = validity_check(df,cb_shifts)
        self.assertTrue(count==0,msg=('Simulation uses chemical shifts {}. Invalid EDS Chemical shift found for data {}!'.format(Shifts_index,Invalid_index)))
        
    def test_EDS_ca(self):
        df_ca = pd.DataFrame({'0':['#1','2','3','4','5','6','7','8','9','10','11','#12','#13','14','15','16','17','18','19','20','21','22','23','#24'],'1':[1.5,2,1.5,2,1,1,1.5,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]})
        shifts_dat = "Shifts.dat"
        df_ca.to_csv(shifts_dat,index=False,index_label=False)
        df = pd.read_csv(shifts_dat)
        txt_plumed = ("PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=250\n"+                     
        "eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_10,cs.hn_13,cs.hn_19,cs.hn_20,cs.hn_21,cs.nh_14,cs.nh_15,cs.nh_16,cs.nh_22,,cs.cb_1,cs.cb_11,cs.cb_16 CENTER_ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_6,cs.exphn_7,cs.exphn_8,cs.exphn_9,cs.exphn_10,cs.exphn_13,cs.exphn_14,cs.exphn_15,cs.exphn_16,cs.exphn_18,cs.exphn_19,cs.exphn_20,cs.exphn_21,cs.exphn_22,cs.expnh_1,cs.expnh_2,cs.expnh_3,cs.expnh_4,cs.expnh_6,cs.expnh_7,cs.expnh_8,cs.expnh_9,cs.expnh_10,cs.expnh_13,cs.expnh_14,cs.expnh_15,cs.expnh_16,cs.expnh_18,cs.expnh_19,cs.expnh_20,cs.expnh_21,cs.expnh_22 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 PERIOD=250 TEMP=@replicas:{278.0,284.825679594,291.818948762,298.9839223,306.324816031,313.84594929,321.551747461,329.446744586,337.535586031,345.823031217,354.313956423,363.013357653,371.926353579,381.05818855,390.414235678,400.0} OUT_RESTART=eds_bias_restart_correct3.dat\n"+
        "PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat")
        plumed_dat = "Plumed.dat"
        with open(plumed_dat,"w") as f1:
            f1.write(txt_plumed)
        file_plumed = open(plumed_dat, 'r')
        hn_shifts,nh_shifts,ca_shifts,cb_shifts = parser(file_plumed)
        file_plumed.close()
        if ca_shifts ==[]:
                raise ValueError('Desired chemical shift is unavailable in the provided plumed data!')
        os.remove(shifts_dat)
        count = validity_check(df,ca_shifts) 
        count, Invalid_index, Shifts_index = validity_check(df,ca_shifts)
        self.assertTrue(count==0,msg=('Simulation uses chemical shifts {}. Invalid EDS Chemical shift found for data {}!'.format(Shifts_index,Invalid_index)))
        

    def test_EDS_nh(self):
        df_nh = pd.DataFrame({'0':['#1','2','3','4','5','6','7','8','9','10','11','#12','#13','14','15','16','17','18','19','20','21','22','23','#24'],'1':[1.5,2,1.5,2,1,1,1.5,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]})
        shifts_dat = "Shifts.dat"
        df_nh.to_csv(shifts_dat,index=False,index_label=False)
        df = pd.read_csv(shifts_dat)
        txt_plumed = ("PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=250\n"+                     
        "eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_10,cs.hn_13,cs.hn_19,cs.hn_20,cs.hn_21,cs.nh_14,cs.nh_15,cs.nh_16,cs.nh_22,,cs.cb_1,cs.cb_11,cs.cb_16 CENTER_ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_6,cs.exphn_7,cs.exphn_8,cs.exphn_9,cs.exphn_10,cs.exphn_13,cs.exphn_14,cs.exphn_15,cs.exphn_16,cs.exphn_18,cs.exphn_19,cs.exphn_20,cs.exphn_21,cs.exphn_22,cs.expnh_1,cs.expnh_2,cs.expnh_3,cs.expnh_4,cs.expnh_6,cs.expnh_7,cs.expnh_8,cs.expnh_9,cs.expnh_10,cs.expnh_13,cs.expnh_14,cs.expnh_15,cs.expnh_16,cs.expnh_18,cs.expnh_19,cs.expnh_20,cs.expnh_21,cs.expnh_22 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 PERIOD=250 TEMP=@replicas:{278.0,284.825679594,291.818948762,298.9839223,306.324816031,313.84594929,321.551747461,329.446744586,337.535586031,345.823031217,354.313956423,363.013357653,371.926353579,381.05818855,390.414235678,400.0} OUT_RESTART=eds_bias_restart_correct3.dat\n"+
        "PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat")
        plumed_dat = "Plumed.dat"
        with open(plumed_dat,"w") as f1:
            f1.write(txt_plumed)
        file_plumed = open(plumed_dat, 'r')
        hn_shifts,nh_shifts,ca_shifts,cb_shifts = parser(file_plumed)
        file_plumed.close()
        if nh_shifts ==[]:
                raise ValueError('Desired chemical shift is unavailable in the provided plumed data!')
        os.remove(shifts_dat)
        count, Invalid_index, Shifts_index = validity_check(df,nh_shifts)
        self.assertTrue(count==0,msg=('Simulation uses chemical shifts {}. Invalid EDS Chemical shift found for data {}!'.format(Shifts_index,Invalid_index)))
        

    def test_EDS_hn(self):
        df_hn = pd.DataFrame({'0':['#1','2','3','4','5','6','7','8','9','10','11','#12','#13','14','15','16','17','18','19','20','21','22','23','#24'],'1':[1.5,2,1.5,2,1,1,1.5,1,1,1,1.5,1,1,1.2,1,1,1,1,1,1.2,1.2,1,1,1]})
        shifts_dat = "Shifts.dat"
        df_hn.to_csv(shifts_dat,index=False,index_label=False)
        df = pd.read_csv(shifts_dat)
        txt_plumed = ("PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=250\n"+                     
        "eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_10,cs.hn_13,cs.hn_19,cs.hn_20,cs.hn_21,cs.nh_14,cs.nh_15,cs.nh_16,cs.nh_22,,cs.cb_1,cs.cb_11,cs.cb_16 CENTER_ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_6,cs.exphn_7,cs.exphn_8,cs.exphn_9,cs.exphn_10,cs.exphn_13,cs.exphn_14,cs.exphn_15,cs.exphn_16,cs.exphn_18,cs.exphn_19,cs.exphn_20,cs.exphn_21,cs.exphn_22,cs.expnh_1,cs.expnh_2,cs.expnh_3,cs.expnh_4,cs.expnh_6,cs.expnh_7,cs.expnh_8,cs.expnh_9,cs.expnh_10,cs.expnh_13,cs.expnh_14,cs.expnh_15,cs.expnh_16,cs.expnh_18,cs.expnh_19,cs.expnh_20,cs.expnh_21,cs.expnh_22 MULTI_PROP=0.4 RANGE=0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 PERIOD=250 TEMP=@replicas:{278.0,284.825679594,291.818948762,298.9839223,306.324816031,313.84594929,321.551747461,329.446744586,337.535586031,345.823031217,354.313956423,363.013357653,371.926353579,381.05818855,390.414235678,400.0} OUT_RESTART=eds_bias_restart_correct3.dat\n"+
        "PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat")
        plumed_dat = "Plumed.dat"
        with open(plumed_dat,"w") as f1:
            f1.write(txt_plumed)
            
        file_plumed = open(plumed_dat, 'r')
        hn_shifts,nh_shifts,ca_shifts,cb_shifts = parser(file_plumed)
        file_plumed.close()
        if hn_shifts ==[]:
                raise ValueError('Desired chemical shift is unavailable in the provided plumed data!')
        os.remove(shifts_dat)
        count, Invalid_index, Shifts_index = validity_check(df,hn_shifts)
        self.assertTrue(count==0,msg=('Simulation uses chemical shifts {}. Invalid EDS Chemical shift found for data {}!'.format(Shifts_index,Invalid_index)))
        

if __name__ == "__main__":
    unittest.main()
    
    
