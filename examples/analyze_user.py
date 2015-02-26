"""
scp -rp aidanmac@ebi-003.ebi.ac.uk:/gpfs/nobackup/saezgrp/aidan/chem_dimer_paper/[name]/ save/
"""

import pip3bayes as bayes
import pickle
import numpy as np
import sys
import csv
import os.path
sys.path.append('examples/user_models/')

iters = range(0,3)
path = 'save/'

# get model
from user_1 import *

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
    
# setup model options and model
m_opts = bayes.ModelOptions("user_1_scale_free_pip3")

likelihoods = []
for key in d.information.keys():
    values = list(d.information[key])
    for j in range(len(d.information[key])):
        root = '%s_%s_%s' % (key, values[j], m_opts.name)
        for k in iters:
            file_likelihood = '%s%s/%s/%s_iter_%s_likelihood.pck' % (path, m_opts.name, root, root, str(k))
            file_name = '%s%s/%s/%s_iter_%s_positions.pck' % (path, m_opts.name, root, root, str(k))
            try:
                positions = pickle.load(open(file_name,'rb'))
                l = pickle.load(open(file_likelihood,'rb'))
                likelihoods.append(l)
            except IOError:
                print "Can't get positions at %s!" % file_name
                continue
            data = d.single_cell(key, values[j])
            m = bayes.Model(m_opts, type="user", user_params=user_params, user_odes=user_odes, init_cond=init_cond, obs_idx=obs_idx, data=data, observables=observables)
            position = np.power(10, positions[-1,:])
            if k==0:
                m.simulate(plot=True, position=position, data=data, exper=key, run=values[j], save=False)

out = open("out.csv", 'w')
wr = csv.writer(out, quoting=csv.QUOTE_ALL)
wr.writerow(likelihoods)
out.close()

# sanity check: make sure of order of experiments/runs
for key in d.information.keys():
    values = list(d.information[key])
    for j in range(len(d.information[key])):
        for k in range(0,3):
            print values[j]

# get fit
f = open('R/figure_data/figure_4.csv', 'a')
f.write('Count,Experiment,Run,Data,Sim,Best\n')
count = 1
for i in range(len(keep_data)):
    min_likelihood = []
    root = '%s_%s_%s' % (str(keep_data[i][0]), str(keep_data[i][1]), m_opts.name)
    
    for j in iters:
        file_name = '%s%s/%s/%s_iter_%s.pck' % (path, m_opts.name, root, root, str(j))
        if os.path.isfile(file_name) is False:
            print 'This file (%s) does not exist ...' % file_name
        else:
            mcmc = pickle.load(open(file_name,'rb'))
            min_likelihood.append(mcmc.likelihoods[-1])
    
    l = min_likelihood.index(min(min_likelihood))
    file_name = '%s%s/%s/%s_iter_%s_positions.pck' % (path, m_opts.name, root, root, str(l))
    positions = pickle.load(open(file_name,'rb'))
    data = d.single_cell(str(keep_data[i][0]), str(keep_data[i][1]))
    m = bayes.Model(m_opts, type="user", user_params=user_params, user_odes=user_odes, init_cond=init_cond, obs_idx=obs_idx, data=data, observables=observables)
    position = np.power(10, positions[-1,:])
    sim = m.simulate(plot=False, position=position, data=data, save=False)['PH']
    m1 = np.hstack(( \
        np.repeat(count,len(data)).reshape((181,1)), \
        np.repeat(keep_data[i][0],len(data)).reshape((181,1)), \
        np.repeat(keep_data[i][1],len(data)).reshape((181,1)), \
        data['PH'].reshape((181,1)), \
        sim.reshape((181,1)), \
        np.repeat(min_likelihood[l],len(data)).reshape((181,1)), \
    ))
    np.savetxt(f, m1, delimiter=',')
    count = count+1
f.close()
            