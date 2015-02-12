from collections import OrderedDict

# define model: parameters, initial conditions, ODEs

pip3_max = 1
pi3k_max = 1
stim_frame = 20
inhib_frame = 164

# sort alphabetically
user_params = OrderedDict([
    ('kcat_f',0.1),
    ('kcat_r',0.1),
    ('km_f',0.1), 
    ('km_r',0.1),
    ('pi_total',10),
    ('pten',1)
])

init_cond = [0,0]
obs_idx = [0,1]
observables = ['iSH2','PH']
def user_odes(x, t, u):
    sh2_total = x[0]
    pip3 = x[1]
    pi_total = u[5]
    pten = u[6]
    
    d_1 = u[0](t)
    d_2 = ((u[1]*(pi_total-pip3)*sh2_total)/(u[3] + (pi_total-pip3))) - ((u[2]*pip3*pten)/(u[4] + pip3))
    return [d_1, d_2]