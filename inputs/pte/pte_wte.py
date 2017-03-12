from peptidesim import PeptideSim
import textwrap, sys, re, os,json
import sys
import dill as pickle
seq1 = sys.argv[1]
seq2=sys.argv[2]
name = sys.argv[3]
data_file=sys.argv[4]
#configure=sys.argv[5]
debug = False
pickle_name = name + '.pickle'

com_data=json.load(open(data_file))
#try to reload 
if(os.path.exists(pickle_name)):
    print 'loading restart'
    with open(pickle_name, 'r') as f:
        ps = pickle.load(f)
else:
    ps = PeptideSim(name, [seq1,seq2], [1,1], job_name='2mer_{}'.format(name))#config_file=configure)
ps.peptide_density = 0.02
ps.mdrun_driver='gmx_mpi'
ps.mpi_np = com_data['mpi_np']
ps.forcefield='oplsaa'
ps.initialize()



def make_ladder(hot, N, cold=300.):
    return [cold * (hot / cold) ** (float(i) / N) for i in range(N)]

def get_replex_e(ps, replica_number):
    with open(ps.sims[-1].location + '/' + ps.sims[-1].metadata['md-log']) as f:
        p1 = re.compile('Repl  average probabilities:')        
        p2 = re.compile('Repl\s*' + ''.join(['([0-9\.]+\s*)' for _ in range(replica_number - 1)]) + '$')
        ready = False
        for line in f:
            if not ready and p1.findall(line):
                ready = True
            elif ready:
                match = p2.match(line)
                if match:
                    return [float(s) for s in match.groups()]
                
        

#get ndx indices
i0 = ps.ndx['peptide_CA_0']
i1 = ps.ndx['peptide_CA_1']
i0 = [str(x) for x in i0]
i1 = [str(x) for x in i1]

center_of_mass=com_data['{}{}'.format(seq1,seq2)]
atoms1=center_of_mass[0]
atoms2=center_of_mass[1]
hill_height=com_data['hill_height']

#make plumed wte files
plumed_input = textwrap.dedent(
    '''
    RESTART
    ene: ENERGY
    
    
    METAD ...
    LABEL=METADPT
    ARG=ene
    SIGMA=100.0
    HEIGHT={}
    PACE=250
    TEMP=300
    FILE=HILLS_PTWTE
    BIASFACTOR=10
    ... METAD
    
    
    PRINT ARG=ene STRIDE=250 FILE=COLVAR_PTWTE
    '''.format(hill_height)
with open('plumed_wte.dat', 'w') as f:
    f.write(plumed_input)            
ps.add_file('plumed_wte.dat')

    

#make plumed wte + metad file
plumed_input = textwrap.dedent(
    '''
    RESTART
    WHOLEMOLECULES STRIDE=1 ENTITY0={c1} ENTITY1={c2}
    
    c1: COM ATOMS={c1}
    c2: COM ATOMS={c2}
    
    d: DISTANCE ATOMS=c1,c2
    ene: ENERGY

    #load hills but don't add to bias
    METAD ...
    LABEL=METADPT
    ARG=ene
    SIGMA=150.0
    HEIGHT=1.0
    PACE=99999999
    TEMP=300
    FILE=HILLS_PTWTE
    BIASFACTOR=10
    ... METAD

    
    
    METAD ...
    LABEL=METAD
    ARG=d
    SIGMA=0.1
    HEIGHT=1.0
    PACE=500
    TEMP=300
    GRID_MIN=0
    GRID_MAX=10
    GRID_BIN=250
    GRID_WFILE=BIAS
    GRID_WSTRIDE=1000
    #BIASFACTOR=10
    ... METAD
    
    
    PRINT ARG=d,ene STRIDE=10 FILE=COLVAR
    '''.format(c1=atoms1,c2=atoms2))
with open('plumed.dat', 'w') as f:
    f.write(plumed_input)            
ps.add_file('plumed.dat')



ps.run(mdpfile='peptidesim_emin.mdp', tag='init_emin', mdp_kwargs={'nsteps': 0.00001*10**5})
ps.run(mdpfile='peptidesim_anneal.mdp', tag='anneal_nvt')

#equilibrate parallel tempering metadynamics -> add hills until replica exchange efficiency is high enough
    time_ns = com_data['final_time']
if debug:
    time_ns = 0.00001
replicas = com_data['replica_number']
hot = com_data['hot_temperature']
replex_eff = 0
max_iters = com_data['max_iterations']
if debug:
    max_iters = 2
kwargs = [{'nsteps': int(time_ns * 5 * 10 ** 5), 'ref_t': ti} for ti in make_ladder(hot, replicas)]

#take our hill files with us
for i in xrange(replicas):
    ps.add_file('HILLS_PTWTE.{}'.format(i))

for i in range(max_iters):
    with open(pickle_name, 'w') as f:
        pickle.dump(ps, file=f)
    ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_pte_tune_{}'.format(i),  mdp_kwargs=kwargs, run_kwargs={'plumed':'plumed_wte.dat', 'replex': 25}, pickle_name=pickle_name)
    replex_eff = min(get_replex_e(ps, replicas))
    if replex_eff >= 0.3:
        print 'Reached replica exchange efficiency of {}. Continuing to production'.format(replex_eff)
        break
    else:
        print 'Replica exchange efficiency of {}. Continuing simulation'.format(replex_eff)

with open(pickle_name, 'w') as f:
    pickle.dump(ps, file=f)

#now do production metadynamics

if debug:
    time = 0.00005
for kw in kwargs:
    kw['nsteps']=  int(time_ns * 5 * 10 ** 5)
try:
    ps.run(mdpfile='peptidesim_nvt.mdp', tag='nvt_prod'.format(i),  mdp_kwargs=kwargs, run_kwargs={'plumed':'plumed.dat', 'replex': 100}, pickle_name=pickle_name)
finally:
    with open(pickle_name, 'w') as f:
        pickle.dump(ps, file=f)
    
