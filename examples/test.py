"""
To incorporate the H2O2 data, a scaling factor is needed, which is derived from the original data
(pre-inhibition)
"""

import pip3bayes as bayes
import pylab as pl
import numpy as np

pip3_max = 0.5
pi3k_max = 0.08
stim_frame = 20 # frame where dimerizer is added TODO check this
inhib_frame = 164 # frame where competitive dimerizer is added TODO check this

# all data
d = bayes.Data('test', bayes.get_data('SingleCell_corrected.csv'), pip3_max, pi3k_max, stim_frame, inhib_frame)
d.scale()
data = d.single_cell(1, 1, plot=True)

# confirm that the relationship is linear
times = d.raw['time']
x = d.data_scaled
xx = x.loc[(x['readout']=='PH') & (x['time']<times[inhib_frame-stim_frame]), 'value_raw']
yy = x.loc[(x['readout']=='PH') & (x['time']<times[inhib_frame-stim_frame]), 'value_post_inhib_scaled']
# pl.plot(range(len(xx)), xx, '--')
# pl.plot(range(len(yy)), yy)
# pl.show()
pl.plot(xx,yy)
pl.show()

# get a scaling factor
m = (yy[2]-yy[1]) / (xx[2]-xx[1])

# get distributions, take the n points pre fk506 addition
n = 10
t_range = range(inhib_frame-n, inhib_frame)
xx = x.loc[(x['readout']=='PH') & (x['time']>=times[inhib_frame-n]) & (x['time']<times[inhib_frame]), 'value_post_inhib_scaled']
dist = stats.norm.fit(xx)
