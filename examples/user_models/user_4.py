from collections import OrderedDict

# define model: parameters, initial conditions, ODEs

pip3_max = 1
pi3k_max = 1
stim_frame = 20
inhib_frame = 164

# sort alphabetically
user_params = OrderedDict([
    ('k_active',0.1),
    ('k_f',0.1),
    ('k_r',0.1), 
    ('p110_total',1),
    ('pi_total',10),
    ('pten',1)
])

init_cond = [0,0,0] # sh2_total, pip_3, sh2_p110
obs_idx = [0,1]
observables = ['iSH2','PH']
def user_odes(x, t, u):
    sh2_total = x[0]
    pip3 = x[1]
    sh2_p110 = x[2]
    p110_total = u[4]
    pi_total = u[5]
    pten = u[6]

    d_1 = u[0](t)
    d_2 = u[2]*(pi_total-pip3)*sh2_p110 - u[3]*pten*pip3
    d_3 = u[1]*(sh2_total-sh2_p110)*(p110_total-sh2_p110)
    return [d_1, d_2, d_3]