import matplotlib.cm as cm
import matplotlib.gridspec as mgridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import itertools
import multiprocessing
import glob
import pickle
import re
import os

__all__ = ["scatter", "surf", "extract_names", "get_best_likelihood"]

def plot_likelihoods(my_path, count):

    """
    Plots the top n likelihoods per experiment
    """

    # my_path = '/Volumes/aidan/MCMCout/fret/12-06-2014/'
    # count = 10
    
    # list sub-directories (per experiment)
    exper_dirs = os.walk(my_path).next()[1]
    
    # get top-10 likelihoods per folder and store
    to_plot = np.empty((count, len(exper_dirs)))
    for i in range(len(exper_dirs)):
        my_files = glob.glob('%s%s/%s_best_likelihood_*' % (my_path, exper_dirs[i], exper_dirs[i]))
        likelihoods = np.empty(len(my_files))
        for j in range(len(my_files)):
            in_file = my_files[j]
            likelihoods[j] = pickle.load(open(in_file,'rb'))
        to_plot[:,i] = np.sort(likelihoods)[0:count]
    
    # plot
    plt.figure()
    plt.boxplot(to_plot)
    locs, labels = plt.xticks()
    plt.xticks(locs, exper_dirs, fontsize=10, rotation=45)
    plt.show()
    return to_plot

def get_best_likelihood(my_path, my_root):
    
    """
    This function collects the data from the server and
    returns the best result (by likelihood)
    """

    # get likelihood results from file
    my_files = glob.glob('%s/%s/%s_best_likelihood_*' % (my_path, my_root, my_root))
    likelihoods = np.empty((len(my_files),2))
    for a in range(len(my_files)):
        in_file = my_files[a]
        l = pickle.load(open(in_file,'rb'))
        i = re.search('([0-9]+)\.pck',my_files[a])  
        likelihoods[a,] = [int(i.group(1)), l]

    likelihoods = likelihoods[np.argsort(likelihoods[:,0])]

    # plot likelihoods by prior
    plt.plot(likelihoods[:,0],likelihoods[:,1])
    plt.show()

    # check out best solution
    best = likelihoods[np.argmin(likelihoods[:,1],0)][0]
    a = int(best)
    return a

def extract_names(param_object, name):
	"""Returns the index of the named parameter in the parameter list"""
	all_names = [y.name for x,y in enumerate(param_object)]
	return all_names.index(name)

def scatter(mcmc, params_to_display, mask=True):

    """
    Displays a grid of scatter plots for each 2-D projection of an MCMC walk.
    """

    # number of dimensions in position vector
    ndims = len(params_to_display)
    
    # vector of booleans indicating accepted MCMC moves
    accepts = mcmc.accepts.copy()
    
    # mask off the annealing (burn-in) phase, or up to a user-specified step
    if mask is True:
        mask = mcmc.options.anneal_length
    if mask is False:
        mask = 0
    
    accepts[0:mask] = 0
    
    # grab position vectors and posterior values from accepted moves
    positions = mcmc.positions[accepts][:,params_to_display]
    posteriors = mcmc.posteriors[accepts]
    params_values_orig = [y.value for x,y in enumerate(mcmc.options.estimate_params) \
    if x in params_to_display]
    params_names = [y.name for x,y in enumerate(mcmc.options.estimate_params) \
    if x in params_to_display]
    
    # calculate actual range of values on each dimension
    maxes = positions.max(0)
    mins = positions.min(0)
    ranges = abs(maxes - mins)
    
    # use 2% of the maximum range as a margin for all scatter plots
    margin = max(ranges) * 0.02
    
    # calculate upper and lower plot limits based on min/max plus the margin
    lims_top = maxes + margin
    lims_bottom = mins - margin
    
    # calculate new ranges based on limits
    lim_ranges = abs(lims_top - lims_bottom)

    plt.figure()
    
    # build a GridSpec which allocates space based on these ranges
    # include a minimum size
    gs = mgridspec.GridSpec(ndims, ndims) # width_ratios=lim_ranges, height_ratios=lim_ranges[-1::-1])
    
    # build an axis locator for each dimension
    locators = []
    for i, r in enumerate(lim_ranges):
        # place ticks on the integers, unless there is no integer within the
        # given dimension's calculated range
        nbins = np.ceil(r) + 1
        locators.append(mticker.MaxNLocator(nbins=nbins, steps=[2, 10]))

    fignum = 0
    
    for py in reversed(range(len(params_to_display))):
        for px in range(len(params_to_display)):
            ax = plt.subplot(gs[fignum])
            ax.tick_params(left=False, right=True, top=True, bottom=False,
                labelleft=False, labelright=False, labeltop=False,
                labelbottom=False, direction='in')
            ax.yaxis.set_major_locator(mticker.LinearLocator())
            ax.xaxis.set_major_locator(mticker.LinearLocator())
            if px == py:
                # consistent bin density
                bins = 200 * lim_ranges[px] / np.sum(lim_ranges)
                ax.hist(positions[:,px], bins=bins, histtype='stepfilled',
                    color='salmon', ec='tomato')
                ax.vlines(params_values_orig[px], *ax.get_ylim(),
                    color='red', linewidth=2)
                ax.set_xlim(lims_bottom[px], lims_top[px])
                ax.yaxis.set_major_locator(mticker.LinearLocator())
            else:
                # 2-D scatter plots off the diagonal
                ax.plot(positions[:,px], positions[:,py], color='darkblue',
                    alpha=0.2)
                ax.scatter(positions[:,px], positions[:,py], s=1, color='darkblue',
                    alpha=0.2)
                ax.set_xlim(lims_bottom[px], lims_top[px])
                ax.set_ylim(lims_bottom[py], lims_top[py])
            # parameter name labels along left and bottom edge of the grid
            if px == 0:
                ax.set_ylabel(params_names[py], weight='black', size='large',
                    labelpad=10, rotation='horizontal',
                    horizontalalignment='right')
            if py == 0:
                ax.set_xlabel(params_names[px], weight='black', size='large',
                    labelpad=10,)
            
            # tick labels along the right and top edge of the grid
            if True:
                ax.tick_params('y', labelright=True, labelsize=7)
            if py == ndims - 1:
                ax.tick_params('x', labeltop=True, labelsize=7)
            # move to next figure in the gridspec
            fignum += 1
            
def surf(mcmc, dim0, dim1, mask=True, step=1, parallelize=True, margin=0.1, walk=True, rejects=True, gridsize=20):


    """
    Displays the posterior of an MCMC walk on a 3-D surface.
    """ 

    # mask off the annealing (burn-in) phase, or up to a user-specified step
    if mask is True:
        mask = mcmc.options.anneal_length
    elif mask is False:
        mask = 0

    # create masked versions of a few vectors of interest
    display_slice = slice(mask, None, step)
    accept_mask = mcmc.accepts[display_slice]
    reject_mask = mcmc.rejects[display_slice]
    posteriors = mcmc.posteriors[display_slice]

    # filter out the position elements we aren't plotting
    positions = mcmc.positions[display_slice, (dim0,dim1)]

    # build grid of points for sampling the posterior surface
    pos_min = positions.min(0)
    pos_max = positions.max(0)

    # enforce square aspect ratio in x-y plane by recomputing pos_min/max
    pos_max_range = np.max(pos_max - pos_min)
    pos_mean = np.mean([pos_min, pos_max], 0)
    pos_min = pos_mean - pos_max_range / 2
    pos_max = pos_mean + pos_max_range / 2
    margin_offset = (pos_max - pos_min) * margin
    pos_min -= margin_offset
    pos_max += margin_offset

    p0_vals = np.linspace(pos_min[0], pos_max[0], gridsize)
    p1_vals = np.linspace(pos_min[1], pos_max[1], gridsize)
    p0_mesh, p1_mesh = np.meshgrid(p0_vals, p1_vals)

    # calculate posterior value at all gridsize*gridsize points
    posterior_mesh = np.empty_like(p0_mesh)
    position_base = np.median(mcmc.positions, axis=0)

    # use multiprocessing to make use of multiple cores
    idx_iter = itertools.product(range(gridsize), range(gridsize))
    inputs = ((p0_mesh[i0,i1], p1_mesh[i0,i1]) for i0, i1 in idx_iter)
    inputs = itertools.product([mcmc], [position_base], [dim0], [dim1], inputs)
    map_args = surf_calc_mesh_pos, inputs
    # if parallelize:
        # try:
            # pool = multiprocessing.Pool()
            # outputs = pool.map(*map_args)
            # pool.close()
        # except KeyboardInterrupt:
            # pool.terminate()
            # raise
    # else:
    outputs = map(*map_args)

    for i0 in range(gridsize):
        for i1 in range(gridsize):
            posterior_mesh[i0, i1] = outputs[i0 * gridsize + i1]

    posterior_mesh[np.isinf(posterior_mesh)] = 'nan'
    pmesh_min = np.nanmin(posterior_mesh)
    pmesh_max = np.nanmax(posterior_mesh)
    
    # plot 3-d surface
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    polys = ax.plot_surface(p0_mesh, p1_mesh, posterior_mesh,
        rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0.02, alpha=0.2,
        vmin=pmesh_min, vmax=pmesh_max
    )
    if walk:
        ax.plot(positions[accept_mask,0], positions[accept_mask,1],
        posteriors[accept_mask], c='k')
    if rejects:
        ax.scatter(positions[reject_mask,0], positions[reject_mask,1],
        posteriors[reject_mask], marker='.', c='k', alpha=0.3, s=1)

    # ax.set_xlabel('log10(%s)' % mcmc.options.estimate_params[dim0].name)
    # ax.set_ylabel('log10(%s)' % mcmc.options.estimate_params[dim1].name)
    ax.set_zlabel('-ln(posterior)')
    plt.show()

def surf_calc_mesh_pos(args):
    try:
        mcmc, position_base, dim0, dim1, param_vals = args
        p0_val, p1_val = param_vals
        position = position_base.copy()
        position[dim0] = p0_val
        position[dim1] = p1_val
        return mcmc.calculate_posterior(position)[0]
    except KeyboardInterrupt:
        raise RuntimeError()
