import bayessb
import numpy as np
import os
import pickle

__all__ = ["optimize"]

def optimize(m, data, o_opts, write=True, m_2=None, data_2=None, init_random=True):
    """
    The wrapper for optimization if argument 'm_2' is not none, 2 experiments are
    being optimized together
    """
    
    def extract_records(array, names):
        """Convert a structured array and list of names into an n-d array"""
        return np.vstack([array[name] for name in names]).T

    def likelihood(mcmc, position):
        
        if m_2 is None:
            if m.type=="pysb":
                position = mcmc.cur_params(position)
            if m.type=="user":
                position = np.power(10, position)
            y_sim = m.simulate(position=position, plot=False)
            y_sim_obs = extract_records(y_sim, o_opts.fitted_obs)
            # set nan (NAs) to 0
            y_sim_obs = np.nan_to_num(y_sim_obs)
            y_exper_obs = extract_records(data, o_opts.fitted_obs)
            return np.sum(((y_exper_obs - y_sim_obs) * o_opts.weights) ** 2 / (2 * o_opts.exp_var ** 2))
        
        else:
            position = mcmc.cur_params(position)
            
            # experiment 1
            y_sim_1 = m.simulate(position=position, plot=False)
            y_sim_obs_1 = extract_records(y_sim_1, o_opts.fitted_obs)
            y_sim_obs_1 = np.nan_to_num(y_sim_obs_1)
            y_exper_obs_1 = extract_records(data, o_opts.fitted_obs)

            # experiment 2
            y_sim_2 = m_2.simulate(position=position, plot=False)
            y_sim_obs_2 = extract_records(y_sim_2, o_opts.fitted_obs)
            y_sim_obs_2 = np.nan_to_num(y_sim_obs_2)
            y_exper_obs_2 = extract_records(data_2, o_opts.fitted_obs)
            
            # not fitting all data for experiment 2
            keep = range(110,191)
            y_exper_all = np.vstack((y_exper_obs_1, y_exper_obs_2))
            y_sim_all = np.vstack((y_sim_obs_1, y_sim_obs_2[keep,:]))
            return np.sum(((y_exper_all - y_sim_all) * o_opts.weights) ** 2 / (2 * o_opts.exp_var ** 2))
                    
    def prior(mcmc, position):
        """Distance to original parameter values"""
        return np.sum((position - o_opts.prior_mean) ** 2 / ( 2 * o_opts.prior_var))

    def step(mcmc):
        """Print out some statistics every 20 steps"""
        if mcmc.iter % 20 == 0:
            print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  prior=%g  post=%g' % \
            (mcmc.iter, mcmc.sig_value, mcmc.T, float(mcmc.acceptance)/(mcmc.iter+1), \
            mcmc.accept_likelihood, mcmc.accept_prior, mcmc.accept_posterior)

    opts = bayessb.MCMCOpts()
    opts.nsteps = o_opts.nsteps
    opts.likelihood_fn = likelihood
    opts.prior_fn = prior
    opts.step_fn = step
    opts.use_hessian = True
    opts.hessian_period = opts.nsteps / o_opts.n_hessian
    opts.model = m.model
    opts.estimate_params = o_opts.params
    opts.tspan = np.hstack(m.times)
    
    position = [x.value for x in o_opts.params]
    if m.type=="user":
        position = [position[i] for i in m.params_idx]
    
    if init_random:
        rand_position = []
        for i in position:
            r = np.random.normal(i, float(i)/2)
            r = np.max([r,0.01]) # cannot be < 0
            rand_position.append(r)
        position = rand_position
    
    opts.initial_values = position
    bayessb._use_pysb = o_opts.use_pysb
    mcmc = bayessb.MCMC(opts)
    mcmc.run()

    if(write):
        # save
        if not os.path.exists('%s%s/' % (o_opts.path, o_opts.root)):
            os.makedirs('%s%s/' % (o_opts.path, o_opts.root))

        # artefacts - causes problems pickling otherwise
        mcmc.solver = None 
        mcmc.options.prior_fn = None
        mcmc.options.likelihood_fn = None
        mcmc.options.step_fn = None
        
        if m.type=="user":
            mcmc.options.estimate_params = []
            mcmc.options.model = []
        
        outfile = open('%s%s/%s_iter_%s.pck' % (o_opts.path, o_opts.root, o_opts.root, o_opts.iter), 'w')
        pickle.dump(mcmc, outfile)
        outfile.close()
        
        outfile = open('%s%s/%s_iter_%s_positions.pck' % (o_opts.path, o_opts.root, o_opts.root, o_opts.iter), 'w')
        pickle.dump(mcmc.positions, outfile)
        outfile.close()       

        out_likelihood = open('%s%s/%s_iter_%s_likelihood.pck' % (o_opts.path, o_opts.root, o_opts.root, o_opts.iter), 'w')
        pickle.dump(mcmc.likelihoods[o_opts.nsteps-1], out_likelihood)
        out_likelihood.close()

    return(mcmc)
