from peptidesim import PeptideSim
import textwrap, sys, re, os
import sys
import dill as pickle
from shutil import copyfile, move

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

name = sys.argv[1]
debug = False
pickle_name = name + '.pickle'
MPI_NP = 16
peptide_copies=int(sys.argv[2])
replicas=16
#try to reload                                                            
total_no_atoms=0#init
number_chains=0#init
ps=3#initialize
if(os.path.exists(pickle_name)):
    print('loading restart')
    with open(pickle_name, 'rb') as f:
        ps = pickle.load(f)
        #ps.pickle_name=pickle_name
        print(os.getcwd())
        ps.rel_dir_name='.'
def number_all_atoms_chains():
    if (os.path.isdir("{}/data/".format(os.getcwd()))==True):
        template_file="{}/data/template.pdb".format(os.getcwd())
        with open(template_file, 'r') as f:
            lines=f.readlines()
            lines=lines[-3].strip()
            lines=lines.split()
            print(lines)
            return lines[1],lines[5]
    else:
        print("there is no data directory")
total_no_atoms,number_chains=number_all_atoms_chains()



eds_period=250
remd_exchange_period=250
def get_replex_e(ps, replica_number):
    with open(ps.sims[-1].location + '/' + ps.sims[-1].metadata['md-log']) as f:
        p1 = re.compile('Repl  average probabilities:')       
        p2 = re.compile('Repl\s*' + ''.join(['([0-9\.]+\s*)' for _ in range(replica_number - 1)]) + '$')
        ready = False
        answer=-1
        for line in f:
            if not ready and p1.findall(line):
                ready = True
            elif ready:
                match = p2.match(line)
                if match:
                    answer= [float(s) for s in match.groups()]
                    break
        return answer
    
                



pte_result = ps.pte_replica(mpi_np=MPI_NP, max_tries=3,min_iters=1,
                              mdp_kwargs={'nsteps':int(200* 5*10**2) }, replicas=replicas,hills_file_location=os.getcwd(),
                              hot=400,hill_height=0.6,sigma=250,bias_factor=16, eff_threshold=0.14,cold=278,exchange_period=200)

pte_plumed_script=pte_result['plumed']
replica_temps=pte_result['temperatures']
temps = ','.join(str(e) for e in replica_temps)
kwargs = [ {'ref_t': ti} for ti in replica_temps]
print(total_no_atoms, number_chains, eds_period, temps)
#now reload PTE_WTE hills file and run eds with cs2backbone to generate the eds parameters with multiple replicas
plumed_input0=textwrap.dedent(
    '''

    peptide: GROUP ATOMS=1-{}                                                                                      
    WHOLEMOLECULES ENTITY0=peptide                                                                                 
    cs: CS2BACKBONE ATOMS=peptide DATA=data NRES={}                                                           

    #bias the simalation with EDS and chem chifts until one gets good eds convergence
     
    PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=WT_PTE_CS_shifts STRIDE=200                     
    exp_mean: MATHEVAL ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_7,cs.exphn_8,cs.expha_1,cs.expha_2,cs.expha_3,cs.expha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    calc_mean: MATHEVAL ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    num1: MATHEVAL ARG=cs.exphn_1,exp_mean,cs.hn_1,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm1: MATHEVAL ARG=cs.exphn_1,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num2: MATHEVAL ARG=cs.exphn_2,exp_mean,cs.hn_2,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm2: MATHEVAL ARG=cs.exphn_2,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num3: MATHEVAL ARG=cs.exphn_3,exp_mean,cs.hn_3,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm3: MATHEVAL ARG=cs.exphn_3,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num4: MATHEVAL ARG=cs.exphn_4,exp_mean,cs.hn_4,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm4: MATHEVAL ARG=cs.exphn_4,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num5: MATHEVAL ARG=cs.exphn_7,exp_mean,cs.hn_7,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm5: MATHEVAL ARG=cs.exphn_7,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num6: MATHEVAL ARG=cs.exphn_8,exp_mean,cs.hn_8,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm6: MATHEVAL ARG=cs.exphn_8,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num7: MATHEVAL ARG=cs.expha_1,exp_mean,cs.ha_1,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm7: MATHEVAL ARG=cs.expha_1,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num8: MATHEVAL ARG=cs.expha_2,exp_mean,cs.ha_2,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm8: MATHEVAL ARG=cs.expha_2,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num9: MATHEVAL ARG=cs.expha_3,exp_mean,cs.ha_3,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm9: MATHEVAL ARG=cs.expha_3,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num10: MATHEVAL ARG=cs.expha_5,exp_mean,cs.ha_5,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm10: MATHEVAL ARG=cs.expha_5,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    sum_num: COMBINE ARG=num1,num2,num3,num4,num5,num6,num7,num8,num9,num10 POWERS=1,1,1,1,1,1,1,1,1,1 PERIODIC=NO
    sum_dm: COMBINE ARG=dm1,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10 POWERS=1,1,1,1,1,1,1,1,1,1 PERIODIC=NO
    slope_cumul: MATHEVAL ARG=sum_num,sum_dm FUNC=x/y PERIODIC=NO
    int_cumul: MATHEVAL ARG=slope_cumul,exp_mean,calc_mean FUNC=z-x*y PERIODIC=NO
    slope: AVERAGE ARG=slope_cumul
    int: AVERAGE ARG=int_cumul 
    PRINT ARG=slope_cumul,int_cumul,slope,int FILE=slope_intercept.dat 
    cal_cs1: MATHEVAL ARG=slope,int,cs.exphn_1 FUNC=x*z+y PERIODIC=NO
    cal_cs2: MATHEVAL ARG=slope,int,cs.exphn_2 FUNC=x*z+y PERIODIC=NO
    cal_cs3: MATHEVAL ARG=slope,int,cs.exphn_3 FUNC=x*z+y PERIODIC=NO
    cal_cs4: MATHEVAL ARG=slope,int,cs.exphn_4 FUNC=x*z+y PERIODIC=NO
    cal_cs5: MATHEVAL ARG=slope,int,cs.exphn_7 FUNC=x*z+y PERIODIC=NO
    cal_cs6: MATHEVAL ARG=slope,int,cs.exphn_8 FUNC=x*z+y PERIODIC=NO
    cal_cs7: MATHEVAL ARG=slope,int,cs.expha_1 FUNC=x*z+y PERIODIC=NO
    cal_cs8: MATHEVAL ARG=slope,int,cs.expha_2 FUNC=x*z+y PERIODIC=NO
    cal_cs9: MATHEVAL ARG=slope,int,cs.expha_3 FUNC=x*z+y PERIODIC=NO
    cal_cs10: MATHEVAL ARG=slope,int,cs.expha_5 FUNC=x*z+y PERIODIC=NO
    eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 CENTER_ARG=cal_cs1,cal_cs2,cal_cs3,cal_cs4,cal_cs5,cal_cs6,cal_cs7,cal_cs8,cal_cs9,cal_cs10 PERIOD={} TEMP=@replicas:{{{}}} COVAR OUT_RESTART={}/restart_pt_wte.dat MULTI_PROP=0.1 RANGE=1,1,1,1,1,1,1,1,1,1   
    PRINT ARG=eds.bias,eds.force2 STRIDE=200 FILE=eds.dat
    PRINT ARG=cal_cs1,cal_cs2,cal_cs3,cal_cs4,cal_cs5,cal_cs6,cal_cs7,cal_cs8,cal_cs9,cal_cs10 STRIDE=200 FILE=calibrated_shifts.dat
    phi_27ASN_28LYS: TORSION ATOMS=86,88,90,108#phi 27ASN to 28LYS bend                                            
    psi_28LYS_29GLY: TORSION ATOMS=88,90,108,110#psi 28 LYS and 29GLY              com_glu_22: COM ATOMS=23,24,25#coo                        
    com_lys_28: COM ATOMS=104,105,106,107#nh3                                  
    com_glu_22: COM ATOMS=23,24,25#coo                                                                          

    com_asp_23: COM ATOMS=35,36,37#coo                                                                             
    com_val_24: COM ATOMS=44,45,46,47,48,49,50,51,52,53#ch2cch3                                                    
    com_asn_27: COM ATOMS=81,82,83,84,85#conh2
    com_ser_26: COM ATOMS=71,70  #OH       
    com_ala_30: COM ATOMS=121,122,123,124                      

    distance_SER26HN_asp23: DISTANCE ATOMS=com_ser_26,com_asp_23#SER26HN and asp23

    distance_ASN27_glu22: DISTANCE ATOMS=com_asn_27,com_glu_22#asn27 and glu22                                     
    
    coor_number_K28_E22: COORDINATION GROUPA=23-25 GROUPB=104-107 R_0=0.8
    coor_number_K28_D23: COORDINATION GROUPA=35-37 GROUPB=104-107 R_0=1.4
    coor_number_N27_E22: COORDINATION GROUPA=22-25 GROUPB=83-85 R_0=1.4
    coor_number_N27_D23: COORDINATION GROUPA=35-37 GROUPB=83-85 R_0=1.5        
    PRINT ARG=phi_27ASN_28LYS,psi_28LYS_29GLY,distance_ASN27_glu22,coor_number_K28_E22,coor_number_K28_D23,coor_number_N27_D23,coor_number_N27_E22 FILE={} STRIDE=250     

    ENDPLUMED
    '''.format(total_no_atoms,number_chains,eds_period,temps, os.getcwd(),'COLVAR_OUTPUT_EDS_PTE'))
plumed_input=pte_plumed_script+plumed_input0
with open('plumed_eds_conver_pt_wte_metad.dat', 'w') as f:
    f.write(plumed_input)
ps.add_file('plumed_eds_conver_pt_wte_metad.dat')
time_ns=0.050
replex_eff = 0
for kw in kwargs:
    kw['nsteps']=  int(time_ns * 5 * 10 ** 5)
max_iterations=2
min_iterations=1
for i in range(max_iterations):
    ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_conver_eds_{}'.format(i),  mdp_kwargs=kwargs,run_kwargs={'plumed':'plumed_eds_conver_pt_wte_metad.dat', 'replex': remd_exchange_period},mpi_np=MPI_NP)
    with open(ps.pickle_name, 'wb') as f:
        pickle.dump(ps, file=f)
    
    replex_eff = min(get_replex_e(ps, replicas))
    if (replex_eff >= 0.024 and i>=min_iterations):
        print('Reached replica exchange efficiency of {}. Continuing to production'.format(replex_eff))
        break
    elif(replex_eff==-1):
        ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_conver_eds_{}'.format(i),  mdp_kwargs=kwargs,run_kwargs={'plumed':'plumed_eds_conver_pt_wte_metad.dat', 'replex': remd_exchange_period},mpi_np=MPI_NP)
        with open(ps.pickle_name, 'wb') as f:
            pickle.dump(ps, file=f)
    else:
        print('Replica exchange efficiency of {}. Continuing simulation'.format(replex_eff))

print(ps.sim_dicts)
nvt_conver_names=[]
k=0
indeces=[]
for i in ps.sim_dicts:
    
    if i.startswith("nvt_conver_eds"):
        nvt_conver_names.append(i)
        indeces.append(k)
        print(i, k)
    k=k+1    
x=0
adresses=[]
for i in  ps.sims:
    if ("nvt_conver_eds"  in "{}".format(i.location)):
        adresses.append(x)
        print(i.location,x)
    x=x+1

print(indeces[-1], ps.sims[indeces[len(indeces)-1]].location)
eds_conver_ptwte_folder=ps.sims[adresses[-1]].location
print(eds_conver_ptwte_folder)
colvar_file='{}/restart_pt_wte.0.dat'.format(os.path.abspath(eds_conver_ptwte_folder))

plumed_input=textwrap.dedent(
'''
    peptide: GROUP ATOMS=1-{}                                                                                      
    WHOLEMOLECULES ENTITY0=peptide                                                                                 
    cs: CS2BACKBONE ATOMS=peptide DATA=data NRES={}                                                           

    #bias the simalation with EDS and chem chifts until one gets good eds convergence
     
    PRINT ARG=(cs\.cb_.*),(cs\.ca_.*),(cs\.hn_.*),(cs\.ha_.*),(cs\.nh_.*),(cs\.co_.*) FILE=CS_shifts STRIDE=200                     
    exp_mean: MATHEVAL ARG=cs.exphn_1,cs.exphn_2,cs.exphn_3,cs.exphn_4,cs.exphn_7,cs.exphn_8,cs.expha_1,cs.expha_2,cs.expha_3,cs.expha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    calc_mean: MATHEVAL ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 VAR=x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 FUNC=(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10 PERIODIC=NO
    num1: MATHEVAL ARG=cs.exphn_1,exp_mean,cs.hn_1,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm1: MATHEVAL ARG=cs.exphn_1,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num2: MATHEVAL ARG=cs.exphn_2,exp_mean,cs.hn_2,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm2: MATHEVAL ARG=cs.exphn_2,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num3: MATHEVAL ARG=cs.exphn_3,exp_mean,cs.hn_3,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm3: MATHEVAL ARG=cs.exphn_3,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num4: MATHEVAL ARG=cs.exphn_4,exp_mean,cs.hn_4,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm4: MATHEVAL ARG=cs.exphn_4,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num5: MATHEVAL ARG=cs.exphn_7,exp_mean,cs.hn_7,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm5: MATHEVAL ARG=cs.exphn_7,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num6: MATHEVAL ARG=cs.exphn_8,exp_mean,cs.hn_8,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm6: MATHEVAL ARG=cs.exphn_8,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num7: MATHEVAL ARG=cs.expha_1,exp_mean,cs.ha_1,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm7: MATHEVAL ARG=cs.expha_1,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num8: MATHEVAL ARG=cs.expha_2,exp_mean,cs.ha_2,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm8: MATHEVAL ARG=cs.expha_2,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num9: MATHEVAL ARG=cs.expha_3,exp_mean,cs.ha_3,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm9: MATHEVAL ARG=cs.expha_3,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    num10: MATHEVAL ARG=cs.expha_5,exp_mean,cs.ha_5,calc_mean VAR=x,y,z,d FUNC=(x-y)*(z-d) PERIODIC=NO
    dm10: MATHEVAL ARG=cs.expha_5,exp_mean FUNC=(x-y)*(x-y) PERIODIC=NO
    sum_num: COMBINE ARG=num1,num2,num3,num4,num5,num6,num7,num8,num9,num10 POWERS=1,1,1,1,1,1,1,1,1,1 PERIODIC=NO
    sum_dm: COMBINE ARG=dm1,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10 POWERS=1,1,1,1,1,1,1,1,1,1 PERIODIC=NO
    slope_cumul: MATHEVAL ARG=sum_num,sum_dm FUNC=x/y PERIODIC=NO
    int_cumul: MATHEVAL ARG=slope_cumul,exp_mean,calc_mean FUNC=z-x*y PERIODIC=NO
    slope: AVERAGE ARG=slope_cumul #CLEAR=200 moving average of last 200
    int: AVERAGE ARG=int_cumul #CLEAR=200 moving average of last 200
    PRINT ARG=slope_cumul,int_cumul,slope,int FILE=slope_intercept.dat 

    cal_cs1: MATHEVAL ARG=slope,int,cs.exphn_1 FUNC=x*z+y PERIODIC=NO
    cal_cs2: MATHEVAL ARG=slope,int,cs.exphn_2 FUNC=x*z+y PERIODIC=NO
    cal_cs3: MATHEVAL ARG=slope,int,cs.exphn_3 FUNC=x*z+y PERIODIC=NO
    cal_cs4: MATHEVAL ARG=slope,int,cs.exphn_4 FUNC=x*z+y PERIODIC=NO
    cal_cs5: MATHEVAL ARG=slope,int,cs.exphn_7 FUNC=x*z+y PERIODIC=NO
    cal_cs6: MATHEVAL ARG=slope,int,cs.exphn_8 FUNC=x*z+y PERIODIC=NO
    cal_cs7: MATHEVAL ARG=slope,int,cs.expha_1 FUNC=x*z+y PERIODIC=NO
    cal_cs8: MATHEVAL ARG=slope,int,cs.expha_2 FUNC=x*z+y PERIODIC=NO
    cal_cs9: MATHEVAL ARG=slope,int,cs.expha_3 FUNC=x*z+y PERIODIC=NO
    cal_cs10: MATHEVAL ARG=slope,int,cs.expha_5 FUNC=x*z+y PERIODIC=NO
    eds: EDS ARG=cs.hn_1,cs.hn_2,cs.hn_3,cs.hn_4,cs.hn_7,cs.hn_8,cs.ha_1,cs.ha_2,cs.ha_3,cs.ha_5 CENTER_ARG=cal_cs1,cal_cs2,cal_cs3,cal_cs4,cal_cs5,cal_cs6,cal_cs7,cal_cs8,cal_cs9,cal_cs10 TEMP=278 IN_RESTART={}/restart_pt_wte.0.dat MULTI_PROP=0.1 RANGE=1,1,1,1,1,1,1,1,1,1 FREEZE MEAN PERIOD=5000000
    phi_27ASN_28LYS: TORSION ATOMS=86,88,90,108#phi 27ASN to 28LYS bend                                            
    psi_28LYS_29GLY: TORSION ATOMS=88,90,108,110#psi 28 LYS and 29GLY                                              
    phi_ALA1_GLU2: TORSION ATOMS=11,13,14,26#phi ALA and GLU                                                       
    psi_ALA1_GLU2: TORSION ATOMS=1,5,11,13#psi ALA and GLU
    com_glu_22: COM ATOMS=23,24,25#coo
                                                                             
   
    com_lys_28: COM ATOMS=104,105,106,107#nh3                                      
    
    com_asp_23: COM ATOMS=35,36,37#coo                                                                             
    com_val_24: COM ATOMS=44,45,46,47,48,49,50,51,52,53#ch2cch3                                                    
    com_asn_27: COM ATOMS=81,82,83,84,85#conh2
    com_ser_26: COM ATOMS=71,70  #OH       
    com_ala_30: COM ATOMS=121,122,123,124                      
    distance_val24_A30: DISTANCE ATOMS=com_ala_30,com_val_24#val24 and ala30
    distance_SER26HN_glu22: DISTANCE ATOMS=com_ser_26,com_glu_22#SER26HN and glu22

    distance_SER26HN_asp23: DISTANCE ATOMS=com_ser_26,com_asp_23#SER26HN and asp23

    distance_ASN27_glu22: DISTANCE ATOMS=com_asn_27,com_glu_22#asn27 and glu22                                     
    distance_ASN27_asp23: DISTANCE ATOMS=com_asn_27,com_asp_23#asn27 and asp23
    coor_number_K28_E22: COORDINATION GROUPA=23-25 GROUPB=104-107 R_0=0.8
    coor_number_K28_D23: COORDINATION GROUPA=35-37 GROUPB=104-107 R_0=1.4
    coor_number_N27_E22: COORDINATION GROUPA=22-25 GROUPB=83-85 R_0=1.4
    coor_number_N27_D23: COORDINATION GROUPA=35-37 GROUPB=83-85 R_0=1.5        
    PRINT ARG=phi_27ASN_28LYS,psi_28LYS_29GLY,phi_ALA1_GLU2,psi_ALA1_GLU2,distance_SER26HN_glu22,distance_SER26HN_asp23,distance_ASN27_asp23,distance_ASN27_glu22,coor_number_K28_E22,coor_number_K28_D23,coor_number_N27_D23,coor_number_N27_E22,distance_val24_A30 FILE={} STRIDE=250                                       
    PRINT ARG=eds.bias,eds.force2 FILE=eds_output.dat STRIDE=100
    PRINT ARG=cal_cs1,cal_cs2,cal_cs3,cal_cs4,cal_cs5,cal_cs6,cal_cs7,cal_cs8,cal_cs9,cal_cs10 STRIDE=200 FILE=calibrated_shifts.dat
    ENDPLUMED'''.format(total_no_atoms,number_chains,os.getcwd(),'OUTPUT_COLVAR_EDS'))


with open('plumed_eds_colvars_new.dat', 'w') as f:
    f.write(plumed_input)
ps.add_file('plumed_eds_colvars_new.dat')
final_time_eds=int(0.040*5*10**5)
ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod_eds_colvar',  mdp_kwargs={'nsteps': final_time_eds, 'ref_t': 278},run_kwargs={'plumed':'plumed_eds_colvars_new.dat'},mpi_np=MPI_NP)
with open(ps.pickle_name, 'wb') as f:
    pickle.dump(ps, file=f)
#finally:
  
print(ps.box_size_angstrom, replica_temps)

print('done')

