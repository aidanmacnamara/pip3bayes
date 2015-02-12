# define model: parameters, initial conditions, ODEs

pip3_max = 200
pi3k_max = 100

user_params = {'alpha':1, 'k_a':0.1, 'phi':0.005}
init_cond = [0,0,0]
obs_idx = [0,1]
observables = ['iSH2','PH']
def user_odes(x, t, u): 
    ish2_cd1_lck = x[0]
    ph = x[1]
    pip3 = x[2]
    ph_total = u[4]
    d_1 = u[0](t)
    d_2 = -u[1]*ph + u[2]*(ph_total-ph)*pip3
    d_3 = -u[3]*pip3 + ish2_cd1_lck
    return [d_1, d_2, d_3]