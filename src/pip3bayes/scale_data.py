import numpy as np
import pandas as pd
import re
import ggplot as gg
import pylab as pl

__all__ = ["Data", "plot", "get_scaling_factor"]

"""
Specific scaling code for the chemical dimerization data
"""

def get_scaling_factor(d, inhib_frame, stim_frame, n=10, return_max=True):
    """
    To incorporate the H2O2 data, a scaling factor is needed, which is derived from the original data
(pre-inhibition)
    Additionally this function returns the maximum steady-state value for the initial conditions of the H2O2 data
    """
    times = d.raw['time']
    x = d.data_scaled
    xx = x.loc[(x['readout']=='PH') & (x['time']<times[inhib_frame-stim_frame]), 'value_raw']
    yy = x.loc[(x['readout']=='PH') & (x['time']<times[inhib_frame-stim_frame]), 'value_post_inhib_scaled']
    # pl.plot(xx,yy)
    # get a scaling factor
    m = (yy[2]-yy[1]) / (xx[2]-xx[1])

    def mean_by_group(group):
        ss_data = group.loc[(x['time']>=times[inhib_frame-stim_frame-n]) \
        & (x['time']<times[inhib_frame-stim_frame]), 'value_post_inhib_scaled']
        return(np.mean(ss_data))
        
    pip3 = x.loc[(x['readout']=='PH')].groupby(['run','experiment']).apply(mean_by_group)
    pi3k = x.loc[(x['readout']=='iSH2')].groupby(['run','experiment']).apply(mean_by_group)
    
    if return_max:
        return (m, np.max(pip3), np.max(pi3k))
    else:
        return (m, pip3, pi3k)

def plot(data_long, y='value'):
    """
    Plots all the data in a grid format
        'data_long': A Pandas data frame in long-format
        'y': What column to use for the y-axis
    """
    p = gg.ggplot(gg.aes('time', y, color='readout'), data=data_long) + \
        gg.geom_line() + gg.facet_grid("experiment", "run", scales="free")
    return p

def get_raw_data(file_name):
    """
    Get the raw chemical dimerization data and edit (used internally only in 'scale_data' module)
    """

    data = pd.read_csv(file_name, index_col=0)
    # correct the header information
    # headers = data.columns.values.tolist()
    # replace points with underscore
    # headers = [re.sub("[\.]", "_", elem) for elem in headers]
    # correct experiment 1
    # headers = [re.sub('(fk506_(?![0-9]))', 'fk506_1_', elem) for elem in headers]
    # rename
    # data.columns = headers
    return data

def long_format(data):
    """
    Transform the data to long-format (used internally only in 'scale_data' module)
    """
    
    cols = ["time", "value", "readout", "run", "experiment"]
    data_long = pd.DataFrame(index=range(data.shape[0]*(data.shape[1]-1)), columns=cols)
    tExper = data['time'].tolist()
    pattern = '^([\w_]+)_([0-9]+)_([0-9]+)$'
    exper_sets = {} # dictionary to store experiment/replicate names
    
    for i in range(1,data.shape[1]):
        run = re.sub(pattern, '\\3', data.columns.values[i])
        exper = re.sub(pattern, '\\2', data.columns.values[i])
        readout = re.sub(pattern, '\\1', data.columns.values[i])
        exper_sets.setdefault(exper,[]).append(run)
            
        for j in range(data.shape[0]):
            data_long.iloc[data.shape[0]*(i-1)+j,:] = [tExper[j], data.iloc[j,i], readout, run, exper]
            
    for i, j in exper_sets.iteritems():
           exper_sets[i] = set(j)
    return (data_long, exper_sets)

class Data:

    """
    This will read the data, store and plot transformations
        'name': Give the object a name
        'file_name': The file name of the single-cell data
        'pip3_max': What os the maximum PIP3 amounts across all experiments (arbitrary units)?
        'pi3k_max': As above, for PI3K
    """

    def __init__(self, name, file_name, pip3_max, pi3k_max, stim_frame, inhib_frame):
        self.name = 'all data'
        self.file_name = file_name
        self.raw = get_raw_data(self.file_name)
        self.raw_long, self.information = \
        long_format(self.raw)
        self.pip3_max = pip3_max
        self.pi3k_max = pi3k_max
        self.stim_frame = stim_frame
        self.inhib_frame = inhib_frame
        
    def single_cell(self, exper, run, plot=False):
        """
        Return single cell data by experiment and replicate
        """
        
        # make sure arguments are strings
        exper = str(exper)
        run = str(run)
        
        s = self.data_scaled
        s_select = s.loc[(s['run']==run) & (s['experiment']==exper)]
        s_select = s_select.pivot(index='time',columns='readout',values='value_post_inhib_scaled')
        # below is necessary as a structured array is needed to pick columns by name
        # this is important for matching data to simulation correctly and flexibly in likelihood function
        names = s_select.columns.values.tolist()
        formats = ['float64'] * len(names)
        s_return = np.core.records.fromarrays(s_select.values.transpose(), names=names, formats=formats)
        
        if plot:
            pl.ion()
            pl.figure()
            for i in range(len(names)):
                pl.plot(s_select.index, s_return[names[i]], '--', label=names[i])
            pl.legend()
            pl.xlabel("Time (s)")
            pl.ylabel("Molecules/cell")
            pl.title('Experiment %s; Run %s' % (str(exper), str(run)))
            
        return s_return
    
    def scale(self, is_h2o2=False, time_add=27.5*60, shift_params=None):
        """
        Cut data and scale
        Scale H2O2 data
        Assumptions:
            1. Before the addition of H2O2, PH levels have reached a maximum defined by pip3_max
            2. Time: t[0] is 27.5 min. post CD1 addition
            3. Find scaling factor for 'PH' data
        """
        
        if is_h2o2:
            scale_fac = shift_params[0]
            pip3_ss_max = shift_params[1]
            pi3k_ss_max = shift_params[2]
        
        def cut(group):
            group = group.drop(group.index[range(self.stim_frame)])
            return group
        
        def adjust_time(group):
            group['time'] = group['time'] - group['time'][0]
            group['value'] = group['value'] - group['value'][0]
            return group
        
        def scale(data_in, data_min, scale_min, scale_max):
            data_out = (((scale_max-scale_min)*(data_in.astype(float)-data_min)) \
            /(data_in.max()-data_min)) + scale_min
            return data_out
            
        def init_values(group, t):
            init_value = group.loc[(group['time']==t), 'value']
            return init_value
    
        # cut and/or correct time
        if is_h2o2:
            self.data_zero = self.raw_long.copy()
            self.data_zero['time'] = self.data_zero['time'] - self.data_zero['time'][0]
            self.data_zero['time'] = self.data_zero['time'] + time_add
        else:
            self.data_cut = self.raw_long.groupby(['run','experiment','readout']).apply(cut)
            self.data_zero = self.data_cut.groupby(['run','experiment','readout']).apply(adjust_time)
    
        # scale data
        self.data_scaled = self.data_zero.copy()
        
        # keep original values
        self.data_scaled['value_raw'] = self.data_scaled['value']
        
        # scale
        if is_h2o2:
            self.data_scaled.loc[self.data_scaled['readout']=='PH','value'] = self.data_scaled.loc[self.data_scaled['readout']=='PH','value'] * scale_fac
            pip3_pre_shift = self.data_scaled.loc[(self.data_scaled['readout']=='PH')].groupby(['run','experiment']).apply(init_values, t=(1650,))
            self.pip3_pre_shift = pip3_pre_shift
            pip3_shift = np.max(pip3_pre_shift) - pip3_ss_max
            self.data_scaled.loc[self.data_scaled['readout']=='PH','value'] = self.data_scaled.loc[self.data_scaled['readout']=='PH','value'] - pip3_shift
            self.pip3_post_shift = self.data_scaled.loc[(self.data_scaled['readout']=='PH')].groupby(['run','experiment']).apply(init_values, t=(1650,))
            data_in = self.data_scaled.loc[self.data_scaled['readout']=='iSH2','value']
            self.data_scaled.loc[self.data_scaled['readout']=='iSH2','value'] = \
            scale(data_in, data_in.min(), 0, pi3k_ss_max)
            
        else:
            self.data_scaled.loc[self.data_scaled['readout']=='PH','value'] = \
            scale(self.data_scaled.loc[self.data_scaled['readout']=='PH','value'], 0, 0, self.pip3_max)
            self.data_scaled.loc[self.data_scaled['readout']=='iSH2','value'] = \
            scale(self.data_scaled.loc[self.data_scaled['readout']=='iSH2','value'], 0, 0, self.pi3k_max)
        
        def scale_inhib(group, post_inhib_idx):
            v = group['value']
            v[post_inhib_idx:len(v)] = scale(v[post_inhib_idx:len(v)], v.tail(10).mean(), 0, v.iloc[post_inhib_idx])
            group['value_post_inhib_scaled'] = v
            return group
            
        # the data should return to a baseline=0 (pi3k only if h2o2 data)
        if is_h2o2:
            post_inhib_idx = 30
            self.data_scaled['value_post_inhib_scaled'] =  self.data_scaled['value'] # TODO return to this to complete
            self.data_scaled.loc[(self.data_scaled['readout']=='iSH2')] = \
            self.data_scaled.loc[(self.data_scaled['readout']=='iSH2')].groupby(['run','experiment']).apply(scale_inhib, post_inhib_idx)
        else:
            post_inhib_idx = self.inhib_frame - self.stim_frame
            self.data_scaled = self.data_scaled.groupby(['run','experiment','readout']).apply(scale_inhib, post_inhib_idx)

