"""
scp -rp aidanmac@ebi-003.ebi.ac.uk:/gpfs/nobackup/saezgrp/aidan/chem_dimer_paper/[name]/ save/
"""

import pip3bayes as bayes
import pickle
import numpy as np
import sys
import pylab as pl
sys.path.append('examples/user_models/')
import random
from bayessb import convergence

iters = [0,1,2]
path = 'save/'
h2o2_compare = False
test_convergence = False

# get model
model = bayes.model_1
pip3_max = 1
pi3k_max = 1
stim_frame = 20
inhib_frame = 164 

# all data
d = bayes.Data('test', bayes.get_data(), pip3_max, pi3k_max, stim_frame, inhib_frame)
d.scale()

# skip data - poor data
keep_data = np.array((
        [1,1],
        [1,2],
        [1,3],
        [1,4],
        [1,5],
        [2,1],
        [2,3],
        [2,4],
        [2,5],
        [2,6],
        [3,1],
        [3,7],
        [3,9],
        [3,11]
#         [2,2],
#         [3,2],
#         [3,3],
#         [3,4],
#         [3,5],
#         [3,6],
#         [3,8],
#         [3,10],                                       
    ))

# TOADD y = np.where((skip_data==x).all(axis=1)); y.size = 0/1

if h2o2_compare:
    shift_params = bayes.get_scaling_factor(d, inhib_frame, stim_frame)
    h = bayes.Data('test', bayes.get_data('H2O2_data_all_wide_corrected.csv'), pip3_max, pi3k_max, stim_frame, inhib_frame)
    h.scale(is_h2o2=True, time_add=27.5*60, shift_params=shift_params)
    exper_1 = bayes.get_scaling_factor(d, inhib_frame, stim_frame, return_max=False)[1]
    exper_2 = h.pip3_pre_shift
    exper_2_shift = h.pip3_post_shift
    allowed_diff = 0.05

# setup model options and model
m_opts_1 = bayes.ModelOptions('model_1_scale_free_fixed')
m_1 = bayes.Model(m_opts_1, model)

if h2o2_compare:
    m_opts_2 = bayes.ModelOptions('pip3bayes.model_2_scale_free')
    m_opts_2.time_start = 0
    m_opts_2.time_end = 2850
    m_opts_2.mod_sim = {
        'time' : [1800,2100], \
        'species' : [ \
                        ['h2o2(b=None)'], \
                        ['sh2_p110(b=None)', 'sh2(b=None)'] \
                    ], \
        'value' : [[10],[0,0]] \
    }
    m_opts_2.mod_sim['time'].append(m_opts_2.time_end+m_opts_2.gap)
    m_2 = bayes.Model(m_opts_2, model)

converge_res = []
likelihoods = []
for key in d.information.keys():
    values = list(d.information[key])
    
    for j in range(len(d.information[key])):
        root = '%s_%s_%s' % (key, values[j], m_opts_1.name)
        chain_set = []
        
        for k in iters:
            file_name = '%s%s/%s/%s_iter_%s.pck' % (path, m_opts_1.name, root, root, str(k))
            try:
                mcmc = pickle.load(open(file_name,'rb'))
                likelihoods.append(mcmc.likelihoods[-1])
                if test_convergence:
                    chain_set.append(10 ** mcmc.positions[50000:,:])
            except IOError:
                print "Can't get positions at %s!" % file_name
                continue
            data_exper_1 = d.single_cell(key, values[j])
            m_1.simulate(plot=False, position=mcmc.cur_params(), data=data_exper_1, exper=key, run=values[j], save=True)
            if h2o2_compare:
                pre_ss = exper_1.loc[(values[j], key)]
                post_ss = exper_2_shift[(exper_2_shift>pre_ss-allowed_diff) & (exper_2_shift<pre_ss+allowed_diff)]
                run_h2o2, exper_h2o2, _ = random.sample(post_ss.index, 1)[0]
                data_exper_2 = h.single_cell(exper_h2o2, run_h2o2, plot=False)
                m_2.simulate(plot=False, position=mcmc.cur_params(), data=data_exper_2, exper=exper_h2o2, run=run_h2o2, save=False)
        
        W = convergence.within_chain_variances(chain_set)
        var = convergence.parameter_variance_estimates(chain_set)
        converge_res.append(np.sqrt(var/W))


# get distributions/simulations and save to file

# simulation part 1
f = open('R/figure_data/figure_14.csv', 'a')
f.write('Count,Experiment,Run,Data,Sim,Best\n')
count = 1
for i in range(len(keep_data)):
    min_likelihood = []
    root = '%s_%s_%s' % (str(keep_data[i][0]), str(keep_data[i][1]), m_opts_1.name)
    
    for j in iters:
        file_name = '%s%s/%s/%s_iter_%s.pck' % (path, m_opts_1.name, root, root, str(j))
        mcmc = pickle.load(open(file_name,'rb'))
        min_likelihood.append(mcmc.likelihoods[-1])
    
    l = min_likelihood.index(min(min_likelihood))
    file_name = '%s%s/%s/%s_iter_%s.pck' % (path, m_opts_1.name, root, root, str(l))
    mcmc = pickle.load(open(file_name,'rb'))
    data = d.single_cell(str(keep_data[i][0]), str(keep_data[i][1]))
    sim_1 = m_1.simulate(plot=False, position=mcmc.cur_params(), data=data, save=False)['PH']
    m1 = np.hstack(( \
        np.repeat(count,len(data)).reshape((181,1)), \
        np.repeat(keep_data[i][0],len(data)).reshape((181,1)), \
        np.repeat(keep_data[i][1],len(data)).reshape((181,1)), \
        data['PH'].reshape((181,1)), \
        sim_1.reshape((181,1)), \
        np.repeat(min_likelihood[l],len(data)).reshape((181,1)), \
    ))
    np.savetxt(f, m1, delimiter=',')
    count = count+1
f.close()

# simulation part 2
f = open('R/figure_data/figure_15.csv', 'a')
f.write('Count,Experiment,Run,Data,Sim,Best\n')
count = 1
for i in range(len(keep_data)):
    min_likelihood = []
    root = '%s_%s_%s' % (str(keep_data[i][0]), str(keep_data[i][1]), m_opts_1.name)
    
    for j in iters:
        file_name = '%s%s/%s/%s_iter_%s.pck' % (path, m_opts_1.name, root, root, str(j))
        mcmc = pickle.load(open(file_name,'rb'))
        min_likelihood.append(mcmc.likelihoods[-1])
    
    l = min_likelihood.index(min(min_likelihood))
    file_name = '%s%s/%s/%s_iter_%s.pck' % (path, m_opts_1.name, root, root, str(l))
    mcmc = pickle.load(open(file_name,'rb'))
    
    pre_ss = exper_1.loc[(str(keep_data[i][1]), str(keep_data[i][0]))]
    post_ss = exper_2_shift[(exper_2_shift>pre_ss-allowed_diff) & (exper_2_shift<pre_ss+allowed_diff)]
    run_h2o2, exper_h2o2, _ = random.sample(post_ss.index, 1)[0]
    data_2 = h.single_cell(exper_h2o2, run_h2o2, plot=False)
    sim_2 = m_2.simulate(plot=False, position=mcmc.cur_params(), data=data_2, exper=exper_h2o2, run=run_h2o2, save=False)['PH']
    nan_array = np.empty((110,))
    nan_array[:] = np.NAN
    
    m1 = np.hstack(( \
        np.repeat(count,len(sim_2)).reshape((191,1)), \
        np.repeat(keep_data[i][0],len(sim_2)).reshape((191,1)), \
        np.repeat(keep_data[i][1],len(sim_2)).reshape((191,1)), \
        np.hstack((nan_array,data_2['PH'])).reshape((191,1)), \
        sim_2.reshape((191,1)), \
        np.repeat(min_likelihood[l],len(sim_2)).reshape((191,1)), \
    ))
    np.savetxt(f, m1, delimiter=',')
    count = count+1
f.close()

# distributions
f = open('R/figure_data/figure_20.csv', 'a')
f.write('Count,Experiment,Run,Likelihood,Value,Best\n')
count = 1
n = 10000
param = 0
for i in range(len(keep_data)):
    min_likelihood = []
    root = '%s_%s_%s' % (str(keep_data[i][0]), str(keep_data[i][1]), m_opts_1.name)
    
    for j in iters:
        file_name = '%s%s/%s/%s_iter_%s.pck' % (path, m_opts_1.name, root, root, str(j))
        mcmc = pickle.load(open(file_name,'rb'))
        min_likelihood.append(mcmc.likelihoods[-1])
    
    l = min_likelihood.index(min(min_likelihood))
    file_name = '%s%s/%s/%s_iter_%s.pck' % (path, m_opts_1.name, root, root, str(l))
    mcmc = pickle.load(open(file_name,'rb'))

    m1 = np.hstack(( \
        np.repeat(count,n).reshape((n,1)), \
        np.repeat(keep_data[i][0],n).reshape((n,1)), \
        np.repeat(keep_data[i][1],n).reshape((n,1)), \
        np.repeat(min_likelihood[l],n).reshape((n,1)), \
        10 ** mcmc.positions[-n:,param].reshape((n,1)), \
        np.repeat(mcmc.cur_params()[param],n).reshape((n,1)) \
    ))

    np.savetxt(f, m1, delimiter=',')
    count = count+1
f.close()
