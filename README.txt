* A model of PIP3 metabolism using Bayesian parameter estimation

* Test GIT

————— Part 1 —————

* Models
	* user_1: A single forward and reverse reaction with ‘PI3K’ defined as ‘iSH2’ i.e. a spline of the data
		max_pip3 = 1
		max_pi3k = 1
	* user_2: As Bandara et al. 2009 with ‘PI3K’ still defined by the data
		max_pip3 = 200
		max_pi3k = 100
		ph_total = max(ph_per_trace)
	* user_3: As user_1 but with Michaelis-Menten kinetics for PI3K and PTEN activity
	
	### user_4/5 models are flawed ### (if the observable consists of 2 species, it cannot be used)

	* user_4: ‘iSH2 (model)’ == ‘iSH2 (data)’ and a tunable p110_total parameter (the goal here is to include the ‘iSH2’ trace in the model but its effect is to saturate p110 at the PM)
	* user_5: As user_4 but with Michaelis-Menten kinetics for PI3K and PTEN activity


————— Part 2 —————	

* Modify user_4 for post-FK506 dynamics
	* user_6: TOADD


————— Preliminary —————

	* pip3bayes.model_1: A pysb model with ‘PI3K’ modeled as having a finite source 
		max_pip3 = 0.5
		max_pi3k = 0.08 (uM)
		initial conditions:
		pi3k_source = 0.08
		pten = 0.08
		pten = 10
	* pip3bayes.model_2: As model_1 but with a degradation term for the H2O2 inhibition of PIP3
	* model_1_scale_free: As model 1 but with:
		max_pip3 = 1
		max_pi3k = 1
		pten = 1
		pi3k_source = 1
		pip2 = 100
	* model_1_scale_free: as model 2 with model_1_scale_free parameters
	* user_1_scale_free: As user_1 but with ‘model_1_scale_free’ parameters (where applicable)


————— 02/02/2015 —————

* Make sure all files are correctly annotated (from paper >> R plots >> raw data)


————— 02/12/2014 —————

* Move the scaling data code out of ‘prepdata’ - getting very verbose

* Some hacks currently
	* keep = range(110,191) - line 48 optimize.py
	* post_inhib_idx = 30 - line 215 scale_data.py


————— 26/11/2014 —————

* Stitch the H2O2 data together

* Modify ‘model_2’ to incorporate H2O2 inhibition


————— 12/11/2014 —————

* Is the ‘PH’ variation explained by ‘iSH2’? Or is it noise?

* Is the ‘iSH2’ trace representative of active PI3K? - Test 2 alternative model sets to explore this

* Can these models predict the H2O2 data?

* What does the model say about PIP3 turnover at the PM?


————— TODO —————

* Be careful with references:
	* model = m.model
	* m_opts = m.options
	* mcmc.options = o_opts

* Check the likelihood weighting is working okay

