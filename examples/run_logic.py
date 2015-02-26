from pip3bayes import *

# TODO
# 1. import sif and midas
# 2. make a ModelSIF() class that returns a model instance with a simulate function
# 3. should be able to pass above (+ midas) to optimize()

# pick experiment and replicate
exper = 1
run = 1
model = model_1
name = 'model_1'

# all data
d = Data('test', get_data())
d.scale()
# plot(d.raw_long)

# get single cell data
data = d.single_cell(exper, run)

# model
m_opts = ModelOptions(name)
m = Model(model, m_opts)
m.simulate(plot=True)

# optimize
o_opts = OptimizeOptions(exper, run, m, data)
mcmc = optimize(m, data, o_opts)

# plot
m.simulate(plot=True, position=mcmc.cur_params(), data=data)
