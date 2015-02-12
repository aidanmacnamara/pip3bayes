"""
"""

import pip3bayes as bayes
import sys
import scipy.stats as stats

# pick experiment and replicate
exper = 1
run = 1
iter = 0
path = './save/model_1_scale_free_pip3_fit/'
# exper = sys.argv[1]
# run = sys.argv[2]
# iter = sys.argv[3]
# path = sys.argv[4]

model = bayes.model_1
pip3_max = 1
pi3k_max = 1
stim_frame = 20 # frame where dimerizer is added TODO check this
inhib_frame = 164 # frame where competitive dimerizer is added TODO check this

# all data
d = bayes.Data('test', bayes.get_data('SingleCell_corrected.csv'), pip3_max, pi3k_max, stim_frame, inhib_frame)
d.scale()
# plot(d.raw_long)

# get single cell data
data = d.single_cell(exper, run)

# model
m_opts = bayes.ModelOptions('model_3_scale_free_pip3_fit')
m = bayes.Model(m_opts, model)
# m.simulate(plot=True)

# optimize
o_opts = bayes.OptimizeOptions(exper, run, m, data)
o_opts.weights[o_opts.weight_range,:] = 5
o_opts.params = m.model.parameters[0:3] # what params to estimate
o_opts.prior_mean = [i.value for i in o_opts.params]
o_opts.prior_var = 10
o_opts.nsteps = 20000
o_opts.exp_var = stats.norm.fit(data['PH'][-10:])[1]
o_opts.iter = iter
o_opts.path = path
# o_opts.fitted_obs = ['PH', 'iSH2']
o_opts.fitted_obs = ['PH']
mcmc = bayes.optimize(m, data, o_opts, write=False, init_random=True)

# plot
# m.simulate(plot=True, position=mcmc.cur_params(), data=data)
# bayes.scatter(mcmc, [1,2])
