from pysb import *

__all__ = ["model"]

def catalyze(enz, sub, prod, kf, kr, kc):
    """2-step catalytic process"""
    r1_name = 'bind_%s_%s' % (sub.name, enz.name)
    r2_name = 'produce_%s_via_%s' % (prod.name, enz.name)
    E = enz(b=None)
    S = sub(b=None)
    ES = enz(b=1) % sub(b=1)
    P = prod(b=None)
    Rule(r1_name, E + S <> ES, kf, kr)
    Rule(r2_name, ES >> E + P, kc)

# instantiate a model
Model()

# species
Monomer('sh2', ['b'])
Monomer('sh2_p110', ['b'])
Monomer('pip2', ['b'])
Monomer('pip3', ['b'])
Monomer('pten', ['b'])

# reactions
# 1. activation of sh2
Parameter('k_activate', 0.1)
Rule('sh2_p110_activation', sh2(b=None) >> sh2_p110(b=None), k_activate)

# 2. pip3 production
Parameter('kf_1', 0.1)
Rule('pip3_activation', pip2(b=None) + sh2_p110(b=None) >> pip3(b=None) + sh2_p110(b=None), kf_1)
# Parameter('kf_1', 0.1)
# Parameter('kr_1', 0.1)
# Parameter('kc_1', 0.1)
# catalyze(sh2_p110, pip2, pip3, kf_1, kr_1, kc_1)

# 3. pip3 degradation
Parameter('kf_2', 0.1)
Rule('pip3_degradation', pip3(b=None) + pten(b=None) >> pip2(b=None) + pten(b=None), kf_2)
# Parameter('kf_2', 0.1)
# Parameter('kr_2', 0.1)
# Parameter('kc_2', 0.1)
# catalyze(pten, pip3, pip2, kf_2, kr_2, kc_2)

# initial conditions
Parameter('pip2_0', 10)
Parameter('pten_0', 1)
Parameter('sh2_0', 1)

Initial(pip2(b=None), pip2_0)
Initial(pten(b=None), pten_0)
Initial(sh2(b=None), sh2_0)

# observables
Observable('PH', pip3(b=None))
Observable('iSH2', sh2_p110(b=None))
