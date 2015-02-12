import numpy as np
import scipy
import pandas as pd
import os
from os.path import join as pj
import pip3bayes # used in get_data()
import pylab as pl
from pysb.integrate import Solver
from collections import namedtuple
import scipy.interpolate as inter
from scipy.integrate import odeint
import shutil

__all__ = ["ModelOptions", "Model", "OptimizeOptions", "get_data", "setup_cluster"]

"""
A set of classes and functions to prepare data and models for BayesSB
"""

def setup_cluster(d, iters, cluster_dir, run_script, path, home=False):
    """
    A function that writes the bash files to run 'run_script' on the cluster
    TODO if dir arguments don't have trailing /, add one
    """
    
    if os.path.exists(cluster_dir):
        shutil.rmtree(cluster_dir)
    os.makedirs(cluster_dir)
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    # the first line of each bash script
    line_1 = "#!/bin/bash\n"
    
    if home:
        file_name = '%shome_run.sh' % cluster_dir
        f = open(file_name, 'w')
        f.write(line_1)
        # jobs
        for key in d.information.keys():
            for j in range(len(d.information[key])):
                values = list(d.information[key])
                for k in iters:
                    line_2 = 'python %s %s %s %s %s\n' % (run_script, key, values[j], iters[k], path)
                    f.write(line_2)
        f.close()
        
    else: 
        count = 0
        file_names = []
  
        # the second line
        for key in d.information.keys():
            for j in range(len(d.information[key])):
                values = list(d.information[key])
                for k in iters:
                    line_2 = 'python %s %s %s %s %s' % (run_script, key, values[j], iters[k], path)
                    file_name = '%srun_%s.sh' % (cluster_dir, count)
                    f = open(file_name, 'w')
                    f.write(line_1 + line_2)
                    f.close()
                    file_names.append(file_name)
                    count += 1

        # master script
        line_2 = "echo 'Running single cell optimization'\n"
        for k in file_names:
            line_3 = 'bsub -o ./optOutSC.txt -R "select[gpfs]" %s\n' % k
            line_2 = line_2 + line_3
  
        file_name = '%soptRunSC.sh' % cluster_dir
        f = open(file_name, 'w')
        f.write(line_1 + line_2)
        f.close()
        os.system('chmod 755 %s*' % cluster_dir)
 
def get_data(filename='SingleCell_corrected.csv'):
    """
    Get the file path of the single cell data in the package
    """
    packagedir = pip3bayes.__path__[0]
    dirname = pj(os.path.dirname(packagedir), 'share', 'data')
    fullname = os.path.join(dirname, filename)
    return fullname

def extract_species_idx(model, name):
    """
    Returns the model index of the named species (used internally only in 'prep_data' module)
    """
    model_species_to_str = []
    for i in range(len(model.species)):
        model_species_to_str.append(str(model.species[i]))
    species_idx = model_species_to_str.index(name)
    return species_idx

class ModelOptions:
    """
    User-defined options to pass to Model class
        'mod_sim': This is equivalent to 'Events' in Copasi - a dictionary that contains perturbation information
        'gap': The gap between successive measurements (should be constant)
        'time_end': The time of the last measurement
        'name': The name of the model
    """
    def __init__(self, name):
        # set to None if there are no perturbations in the model
        self.mod_sim = {
            'time' : [2160], \
            'species' : [['sh2_p110(b=None)','sh2(b=None)']], \
            'value' : [[0,0]] \
            }
        self.gap = 15
        self.time_start = 0
        self.time_end = 2700
        self.name = name
        self.add_end_time()
        
    def add_end_time(self):
        self.mod_sim['time'].append(self.time_end+self.gap)
        

class Model():
    """
    Adds some additional functions and variables to the pysb model class. This takes different arguments depending on the simulation type ("pysb", "user", "logic")
    
    'type': The type of model:
        'pysb': A model as defined through pysb
        'user': A user-defined model (i.e. a set of ODEs)
        'logic': A CellNOpt-style model
    
    'pysb'
    
    'model': The pysb model
    'options': As defined in ModelOptions()
     
    'user'
    
    'user_params': A dictionary giving the set of parameter names and values
    'user_odes': A function of the type f(x, t, u) that can be passed to odeint()
    'init_cond': A list of initial conditions for the model species
    'obs_idx': A list of indices that matches 'observables' to their respective ODEs
    'observables': A list of names that correspond to those in 'data' i.e. the traces to fit
    """
    
    def __init__(self, options, model=None, type="pysb", user_params=None, user_odes=None, init_cond=None, data=None, obs_idx=None, observables=None):
        
        if type=="pysb":
            
            self.type = type
            self.model = model
            self.options = options
        
            # set up solver for each time interval
            times = []
            solvers = []
                
            if options.mod_sim is None:
                times.append(np.arange(options.time_start, options.time_end, options.gap))
                solvers.append(Solver(self.model, times[0]))
            else:
                times.append(np.arange(options.time_start, options.mod_sim['time'][0], options.gap))
                solvers.append(Solver(self.model, times[0]))         
                for i in range(len(options.mod_sim['time'])-1):
                    times.append(np.arange(options.mod_sim['time'][i], options.mod_sim['time'][i+1], options.gap))
                    solvers.append(Solver(self.model, times[i+1]))
        
            self.times = times
            self.solvers = solvers

        if type=="user":
            
            self.options = options
            self.times = np.arange(0, options.time_end+options.gap, options.gap)
            Parameter = namedtuple('Parameter', 'name value')
            Model = namedtuple('Model', 'parameters')
            params_list = set()
            for key, value in user_params.iteritems():
                params_list.add(Parameter(key, value))
            # order is list here, find the indices for 'user_params' order
            params_names = [x.name for x in params_list]
            self.params_idx = [params_names.index(x) if x in params_names else None for x in user_params.keys()]
            self.model = Model(params_list)
            spline = inter.InterpolatedUnivariateSpline(self.times, data['iSH2'])
            self.real_ish2 = spline.derivative()
            self.user_odes = user_odes
            self.init_cond = init_cond
            self.obs_idx = obs_idx
            self.observables = observables # TODO move this
            self.type = type
            
    def simulate(self, position=None, plot=False, data=None, exper=None, run=None, save=False):
        """
        Simulate the model
        TODO Expand this to include different simulation options
        i.e. (including H2O2 data?)
        """
        
        if position is None:
            position = [x.value for x in self.model.parameters] # TODO returning wrong order for user model
            if self.type=="user":
                position = [position[i] for i in self.params_idx]
        
        if self.type=="pysb":
        
            self.solvers[0].run(position)
            sim_all = [self.solvers[0].yobs]
            time = self.times[0]
    
            if self.options.mod_sim is not None: # stitch together simulations after perturbations
                time = np.hstack(self.times)
                for i in range(1,len(self.times)):
                    init_values = self.solvers[i-1].y[-1].copy()
                    for j in range(len(self.options.mod_sim['species'][i-1])):
                        j_idx = extract_species_idx(self.model, self.options.mod_sim['species'][i-1][j])
                        init_values[j_idx] = self.options.mod_sim['value'][i-1][j]
                    self.solvers[i].run(y0=init_values, param_values=position)
                    sim_all.append(self.solvers[i].yobs)
        
            y_sim = np.hstack(sim_all)
            
        if self.type=="user":
            
            time = self.times
            my_args = [self.real_ish2]
            for i in range(len(position)):
                my_args.append(position[i])
            my_args = tuple(my_args)
            soln = odeint(self.user_odes, self.init_cond, self.times, args=(my_args,))
            soln = soln[:,self.obs_idx] 
            formats = ['float64'] * len(self.obs_idx)
            y_sim = np.core.records.fromarrays(soln.transpose(), names=self.observables, formats=formats)
            
        if plot:
            if self.type=="pysb":
                obs = [i.name for i in self.model.observables]
            if self.type=="user":
                obs = self.observables
            pl.ion()
            pl.figure()
            for i in range(len(obs)):
                pl.plot(time, y_sim[obs[i]], label=obs[i])
                if data is not None:
                    if len(time)==len(data):
                        pl.plot(time, data[obs[i]], '--', label=obs[i])
                    else:
                        time_range = range(len(time)-len(data), len(time))
                        pl.plot(time[time_range], data[obs[i]], '--', label=obs[i])       
            pl.legend()
            pl.xlabel("Time (s)")
            pl.ylabel("Molecules/cell")
            pl.title('Experiment %s; Run %s' % (str(exper), str(run)))
            if save:
                pl.savefig('figures/%s_%s_%s.pdf' % (str(exper), str(run), self.options.name), format='pdf')
            else:
                pl.show()
        
        return(y_sim)

class OptimizeOptions:
    """
    This contains the (data) options for running MCMC
    """
    def __init__(self, exper, run, m, data):
        
        exper = str(exper)
        run = str(run)
        
        self.write = True # write mcmc results to file
        self.path = 'save/' # output directory for above
        self.root = '%s_%s_%s' % (exper, run, m.options.name) # output file root
        self.iter = 0
        self.fitted_obs = ['PH','iSH2'] # what species should be fitted
        self.weight_value = 2 # if not equal to one, weigh the fit to certain time points
        self.weight_range = range(144,154) # time index for weighting
        self.weights = np.ones((len(data), len(self.fitted_obs)))
        self.exp_var = 0.02
        self.nsteps = 20000
        self.n_hessian = 5
        self.params = None
        self.prior_mean = None
        self.prior_var = None
        self.use_pysb = True
