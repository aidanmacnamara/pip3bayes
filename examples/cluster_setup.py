import pip3bayes as bayes

pip3_max = 0.5
pi3k_max = 0.08
stim_frame = 20
inhib_frame = 164

# all data
d = bayes.Data('test', bayes.get_data(), pip3_max, pi3k_max, stim_frame, inhib_frame)
d.scale()

bayes.setup_cluster(d, [0,1,2], './cluster/', './examples/run.py', './save/model_1_scale_free/', home=True)