"""
Check the model predictions against the H2O2 data

Cluster instructions:

    1. source pip3bayesVE/bin/activate
    2. unset PYTHONPATH
    3. Copy pip3bayes to cluster
    4. python setup.py install
    5. Run cluster_setup.py
"""

import pip3bayes as bayes
import sys
import scipy.stats as stats
import pylab as pl
import numpy as np
import random

# pick experiment and replicate
exper = 1
run = 1
iter = 0
path = './save/model_2_scale_free_fixed/'
# exper = sys.argv[1]
# run = sys.argv[2]
# iter = sys.argv[3]
# path = sys.argv[4]

model = bayes.model_2
pip3_max = 1
pi3k_max = 1
stim_frame = 20 # frame where dimerizer is added TODO check this
inhib_frame = 164 # frame where competitive dimerizer is added TODO check this

# experiment 1
d = bayes.Data('test', bayes.get_data('SingleCell_corrected.csv'), pip3_max, pi3k_max, stim_frame, inhib_frame)
d.scale()
data_exper_1 = d.single_cell(exper, run, plot=False)

# experiment 2
shift_params = bayes.get_scaling_factor(d, inhib_frame, stim_frame)
h = bayes.Data('test', bayes.get_data('H2O2_data_all_wide_corrected.csv'), pip3_max, pi3k_max, stim_frame, inhib_frame)
h.scale(is_h2o2=True, time_add=27.5*60, shift_params=shift_params)

# compare ph steady-state distributions for both experiments
exper_1 = bayes.get_scaling_factor(d, inhib_frame, stim_frame, return_max=False)[1]
exper_2 = h.pip3_pre_shift
exper_2_shift = h.pip3_post_shift
bins = np.linspace(0, 2, 100)
# pl.hist(exper_1, bins, alpha=0.5, label='exper_1')
# pl.hist(exper_2, bins, alpha=0.5, label='exper_2')
# pl.hist(exper_2_shift, bins, alpha=0.5, label='exper_2_shift')
# pl.legend(loc='upper right')
# pl.show()

# choose the H2O2 experiment that best matches experiment 1
allowed_diff = 0.05
pre_ss = exper_1.loc[(str(run), str(exper))]
post_ss = exper_2_shift[(exper_2_shift>pre_ss-allowed_diff) & (exper_2_shift<pre_ss+allowed_diff)]
run_h2o2, exper_h2o2, _ = random.sample(post_ss.index, 1)[0]
data_exper_2 = h.single_cell(exper_h2o2, run_h2o2, plot=False) # TODO scale iSH2 per experiment?

# model
m_opts_1 = bayes.ModelOptions('model_2_scale_free_fixed')
m_1 = bayes.Model(m_opts_1, model)
# m_1.simulate(plot=True)

m_opts_2 = bayes.ModelOptions('model_2_scale_free_fixed')
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
# m_2.simulate(plot=True)

# optimize
o_opts = bayes.OptimizeOptions(exper, run, m_1, data_exper_1)
o_opts.weights = np.ones((len(data_exper_1)+len(data_exper_2), len(o_opts.fitted_obs)))
o_opts.weights[o_opts.weight_range,:] = 5
o_opts.params = m_1.model.parameters[0:5] # what params to estimate
o_opts.prior_mean = [i.value for i in o_opts.params]
o_opts.prior_var = 10
o_opts.nsteps = 20000
o_opts.exp_var = stats.norm.fit(data_exper_1['PH'][-10:])[1]
o_opts.iter = iter
o_opts.path = path
o_opts.fitted_obs = ['PH']
mcmc = bayes.optimize(m_1, data_exper_1, o_opts, write=False, m_2=m_2, data_2=data_exper_2)

# plot result
m_1.simulate(plot=True, position=mcmc.cur_params(), data=data_exper_1)
m_2.simulate(plot=True, position=mcmc.cur_params(), data=data_exper_2)
# bayes.scatter(mcmc, [1,2])
