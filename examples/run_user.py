"""
This script passes a user-defined set of ODEs to BayesSB

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
sys.path.append('examples/user_models/')

# pick experiment and replicate
exper = 1
run = 1
iter = 0
path = './save/user_1_scale_free/'
# exper = sys.argv[1]
# run = sys.argv[2]
# iter = sys.argv[3]
# path = sys.argv[4]

# get model
from user_3 import *

# all data
d = bayes.Data('test', bayes.get_data(), pip3_max, pi3k_max, stim_frame, inhib_frame)
d.scale()

# get single cell data
data = d.single_cell(exper, run)
# data = data[0:inhib_frame-stim_frame]

# setup model options and model
m_opts = bayes.ModelOptions("user_3_scale_free_pip3")
# m_opts.time_end = 2145
m = bayes.Model(m_opts, type="user", user_params=user_params, user_odes=user_odes, init_cond=init_cond, obs_idx=obs_idx, data=data, observables=observables)
m.simulate(plot=True, data=None, exper=exper, run=run, save=False)

# optimize
o_opts = bayes.OptimizeOptions(exper, run, m, data)
o_opts.weights[o_opts.weight_range,:] = 5
o_opts.use_pysb = False
o_opts.params = m.model.parameters
o_opts.prior_mean = user_params.values()
o_opts.prior_var = 20
o_opts.nsteps = 1000
o_opts.exp_var = stats.norm.fit(data['PH'][-10:])[1]
o_opts.fitted_obs = ['PH']
o_opts.iter = iter
o_opts.path = path
mcmc = bayes.optimize(m, data, o_opts, write=False)

# plot
m.simulate(plot=True, position=mcmc.cur_params(), data=data)
bayes.scatter(mcmc, [0,1])
