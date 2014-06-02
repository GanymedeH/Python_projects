ADP implementation for improving network durability readme file
Last updated May 10, 2013 for the final implementation corresponding to the thesis "Using Reinforcement Learning to Improve Network Durability" by Erik Hammel, submitted for publication May 2013. The code isn't copyrighted, anyone may alter it as s/he sees fit. The only stipulation is that Erik Hammel should be credited as the original author in any future adaptations.

The contents of this folder contains 3 types of files: python scripts, pickles, and text files. The scripts obviously run the code, while the pickle files are containers that contain variable data that the code uses to set up various configurations of the analysis. The details for each file are given below:

Scripts:
Note that in all files, a global function qpick is defined. This is merely a shortcut for the pickle I/O.
dataReports.py: contains quick-and-dirty scripts to plot various performance metrics of the analysis. Different functions correspond to different graphs.
		plot_saa_benchmarks: this is the most pertinent; these graphs are spread throughout the thesis. It plots algorithm performance on an LP as compared to the IP solutions at time 0 and time T, and does so over training iterations. It has the capability to plot multiple LPs on one graph, it just requires the user to input a list of the files.
		
		plot_time_benchmarks: this plots the runtime performance of the training code. As-is, this function is useless because the runtime measurements put in place in the main code are commented out.
		
		plot_vdata: plots the values of the VFA coefficients for each training iteration (on multiple plots). They can be viewed at intervals (modulo time t), they create separate directories for all the plots, and there is an option to sort the coefficient values to get a nice smooth ascending curve.
		
		
exeADP: this is file that actually executes the ADP algorithm itself (the main function, if you will). Running the main script at the bottom will run the training program and output a pickle of VFA coefficient values for all training iterations. It can also output the scenario/outage outcomes generated if the "savestate" option is set to a filename string. State information can likewise be loaded from an external source by specifying the appropriate filename string to "usestate". The file contains a copious amount of comments, but none describing overall purposes of each function, so they are briefly explained here in alphabetical order:
		adp_policy_search: executes the control step of the RL algorithm. Makes time-dependent modifications to the model and directly makes changes to the CPLEX API parameters. Randomizes the return policy if so specified.
		
		adp_true_val_fn_calc: calculates the target value used in the prediction step. 
		
		adp_vfa_update_*: these 4 functions execute the prediction step in the various ways described in the thesis: gradient descent, least-squares, the toggling method, and minimum norm method.
		
		adp_warm_start: initializes the VFA coefficients to values that are calculated by running an optimization problem. As-is, this is not used since we initialize the VFA coefficients to a constant value for all arcs.
		
		alt_*: these 3 functions correspond to the alternate ADP approach suggested by Professor Mitchell as outlined in notes 2012-09-27 Mitchell.pdf (provided). This method never really worked in research and can safely be ignored, but the code is there in case someone wishes to expand upon it.
		
		convert_cplex_structs: converts the CPLEX variable name designations to the integer indices of the internal data structures. For debugging purposes mostly.
		
		create_models_temp_framework: a very basic CPLEX object initialization that adds the installation decisions and any installed arcs' flow costs.
		
		execute_ADP: this is the main function that invokes all of the substeps of the RL procedure. The individual steps as well as the various input option contingencies are pretty well documented.
		
		execute_alt_method: this is the "execute_ADP" equivalent for the alternate method. It invokes all of the alt_* functions.
		
		find_constr_4_vars: takes in CPLEX variable name designations or the internal data structures' index for a set of variables and returns the constraints they are in as well as solution values (if a model is given). This is used for debugging purposes.
		
		generate_per_arc_random_policy: self-explanatory. Corresponds to the same-named step of the RL algorithm provided in Chapter 7 of the thesis.
		
		generate_random_policy: a cruder policy randomizer, this simply chooses the optimal policy with probability 1-eps or a completely random set of installations with probability eps.
		
		init_model_components: initializes the modelComponents variable (aka "mC" everywhere else). This dictionary contains all of the options available for running different configurations. The individual options are outlined below:
				"numscens": the number of outage scenarios incorporated into the policy search for all t < T
				"netobj": the data structure containing all information of the network, as contained in networkGenerator.py
				"mcsim": the data structure containing the simulation engine, as contained in simulationEngine.py
				"slackcap": capacity for slack variables
				"slackcost": per-unit-flow cost objective function coefficient for slack variables
				"pertinstbudget": arc installation constraint for an individual time step
				"numfeats": number of separable sets of VFA coefficients for every installation decision. Currently the code is hard-wired to only accept the "ctr_mass_quad_rel_frame" feature class which is hard-wired to 4 quadrants, so numfeats should always be set to 4.
				"explorepolicy": the value of eps when randomizing the policy from the control step.
				"polmethod": method for finding the optimal policy. This is a complex string consisting of several sub-options in the following format:
						<control method>_<pass term>, where <control method> can be "ud_vfa_combined", "ud_vfa_separate", or "node_constr", corresponding to policy search using the combined control method, disjoint control method, and node-constrained control method, respectively. The first 2 are clearly discussed in the thesis in Chapters 4 and 7, respectively. The node-constrained method is an alternate policy search suggested by Erik Hammel that associates VFA coefficients to nodes rather than installable arcs in order to incorporate a sense of flow into the basis function calculations. The theory is provided in notes 2012-10-12 Self.PDF (provided). This method was shown not to work, but is included in the software in case others may want to expand upon it in the future. <pass_term> can be "pass_udterm", "pass_vfaterm", or "pass_bothterms", corresponding to the portion of the control step's objective that gets used in the prediction step: the unmet demand cost, the VFA term, or both. These options apply only when "ud_vfa_combined" is used; it can be safely omitted otherwise.
				"updmethod": method for updating VFA coefficients in the prediction step. This is a complex string consisting of several sub-options in the following format:
						<update algorithm>_<MSE function>_<descent direction>_<step size type>, where <update algorithm> can be "grad", "least_squares", "lp_toggle", or "min_norm", corresponding to the respective methods in adp_vfa_update_*. <MSE function> can be "linear", "1norm_lin", "inverse", "root_inverse", or "neg_exp" when "grad" is used, "linear" when "least_squares" is used, and is not applicable for the other 2 update algorithms (and so can be safely omitted from the input string). The various options correspond to different VFA function forms: linear is the linear sum of basis function * VFA coefficient * decision variable as is the only method present in the thesis. The other forms of the VFA were shown not to work in practice, but they are included for future development. Their specific functional forms are provided in adp_policy_search, and the theory behind why these might work is provided in Chapter 2 of the thesis. <descent direction> can be "_direct_zeta", "_compl_zeta", or "_withvalfn_zeta". The first corresponds to the basic descent direction provided in Chapter 5 of the thesis, and the last explained in Chapter 7. "_compl_zeta" is a descent direction not discussed in the thesis (as it did not work in practice) that updates VFA coefficients corresponding to all arcs which were NOT installed. Finally, <step size type> can be "tdepss", "constss", "harmonicss", or "mcclainss", the last 3 of which are explained in Chapter 7 of the thesis. "tdepss" is simply a harmonic mean that multiplies the iteration count with the number of time steps (T-t).
				"truevalfn": method for calculating the target value in the prediction step. Valid values are "IPv1a", "IPv1b", "IPv2a", "IPv2b", and "LPv1". "IPv1" corresponds to a time-invariant target calculation using the network state at the last time step, and "IPv2" represents a time-dependent target calculation (it takes into account only partial arc install sets); versions "a" allow the IP calculation to only install "pertinstbudget" number of arcs; version "b" allows the IP to install the cumulative budget over all time steps. Finally, "LPv1" essentially just provides the unmet demand cost given all installations made in the control step at time T.
				"bftype": basis function type. This is a complex string consisting of several sub-options in the following format:
						<base calculation>_<sign>_<multiplicative factor>_<optional weight>_<optional scaling factor>_<persistence option>, where <base_calculation> can be "sp", "ones", or "node_constr". The first refers to the shortest path calculation used in the thesis. "ones" uses a descent direction that updates ALL VFA coefficients at every time step and every iteration, and "node_constr" is the corresponding basis functions to the policy method parameter above. Note that if either of these two options are used, every other component of the input string may be safely omitted, as they are not applicable. <sign> is either "pos" or "neg" and reflects that fact that in the standard formulations in the thesis, the VFA coefficients of the combined method are added to the unmet demand in the control step (thus using "pos"itive basis functions) whereas the disjoint method discounts the VFA calculation from unmet demand (thus using "neg"ative basis functions), although one may use either setting for either policy method. <multiplicative factor> is either "probabilistic" or "bottleneckcapratio", the former referring to the unembellished shortest path calculation and the latter applying the bottleneck capacity ratio factor as described in Chapter 7 of the thesis. <optional weight> is of the form "wt=#" where # is a multiplicative constant to the bottleneck capacity ratio. It may be safely omitted; doing so defaults the value to 1. <optional scaling factor> is of the form "sf=#" where # is a number between 0 and 1 as described by the sigma term in the appropriate section of Chapter 7 of the thesis. It too can be safely omitted; it defaults to 1 when "bottleneckcapratio" is used and 0 when "probabilistic" is used. Finally, <persistence option>, when present, must be "cum"; otherwise, it can be left out. Including this substring into the string causes the VFA coefficients corresponding to arcs that are installed in previous time steps to be updated. In other words, once an arc is chosen for installation and its VFA coefficient is altered at time t, it is also altered at every time thereafter. Omitting this substring causes only arcs installed in a particular time step to have the corresponding VFA coefficients altered for that time step. Note that this is independent of the "_withvalfn_zeta" substring of "updmethod"; that is, even though one might choose the persistence option, the target value installation decisions may still cancel out the update for a particular arc.
				"timeTSAAscens": as "numscens", except it specifies the number of scenarios incorporated into the control step at the last time step. The parameter is separate in case one is relying on "v1"-type "truevalfn"s.
				"Lscens": as "numscens", except it specifies the number of scenarios incorporated into the target calculation.
				The remaining parameters of modelComponents are internal CPLEX variables and should not be altered.
		
		trash_model: CPLEX garbage collection. Used because CPLEX has so many memory leaks.
		
		update_models_constraints: inputs the appropriate constraints corresponding to scenario calculations into the CPLEX data structure. It is cumulative; that is, if scenario information already exists within the model, calling this function with new scenario information will be appended to the existing constraints. To avoid this, one must trash the model and start over.
		
		update_models_orig_scens: inputs the original arcs' information into the CPLEX data structure, initializes the objective function and variable bounds, and optionally calls update_models_constraints. This function is cumulative by default, however can be made non-persistent by setting the reset function argument to True.
		
		
networkGenerator.py: this file serves 2 purposes: 1) generate random networks using a group of parameters and 2) store network data for existing random networks. Generating the random network is a separate exercise that, once a network is chosen, need not be run again for the purposes of training. Once the network is created it simply stores arc and node information to be accessed by exeADP.py and simulationEngine.py. The file contains a copious amount of comments for both the internal functions and the parameters involved in generating the random network. They should be sufficient to understand the code. generate_random_network is the main source for the generation algorithm; classify_nodes, place_edges, and place_nodes perform the task of establishing the physical connections and labelling them appropriately; balance_demand, ensure_connectedness, and ensure_feasibility are operations performed to make sure the network conforms to an operable flow network; the rest are utility functions: construct_network builds a network like generate_random_network, only it uses parsed data from parsers.py (to create the NHC and CLARC networks). create_temp_net constructs the set of installable arcs and their connections to the original network. specify_params initializes the parameters used in random network generation. Each parameter is well document as to how it's used. The rest should be self-explanatory.
		


parsers.py: this file reads in the tabular text files (which are created by exporting database queries) and parses them into the appropriate lists to send to networkGenerator for network creation. It outlines the columnar order for the text files if one wishes to import another network in the future.

 
setup.py: a simple script for compiling the source code into an executable.


simulationEngine.py: this file drives the Monte Carlo simulation process that feeds initialization data to exeADP for the purposes of establishing objective functions and constraints. It first simulates a forecast map of outage probabilities, then samples those probabilities to alter capacities at a given time step. The file contains a copious amount of comments, but none describing overall purposes of each function, so they are briefly explained here in alphabetical order:
		calc_feature_ctrmassquad: populates the feature list by determining which quadrant an arc is in by calculating its center of mass, where the weight is given by the outage %. It is a time-indexed list of dictionaries, which means that different dictionary keys or even different separable-type features can exist at different times to set the VFA coefficient values that are used in the adp_policy_search in exeADP.py. Currently the separable sets of VFA coefficients are just stored in a list, so that the quadrant number itself can be used to provide the index into which set of VFA coefficients to use. If one were interested in introducing a new means to partition VFA coefficients, its indexing scheme need not be done this way; any method will suffice so long as the hash method is consistent throughout the code.
		
		clear_scenarios: deletes all outage-altered capacities. If this is not invoked and there exists such information, new outage-altered capacities will just be appended to the network arc information. See sample_scenarios.
		
		eval_adp_basis_functions: as it says. This is where "bftype" (see init_model_components in exeADP.py) specifications are used. The basis functions are calculated to provide coefficients in the VFA term in the control step.
		
		extract_subnet: internal function that effectively determines which 'zone' an arc belongs to for the purposes of assigning it an outage probability within a particular subinterval of the full range of possible outage percentages. It's OK not to fully understand this function. ^_^
		
		generate_distributions: executes the parameter-specified Monte Carlo simulation to assign outage probabilities to each arc for each time step. The basic way it does this is by establishing a gradient of percentages across concentric circular regions.
		
		sample_scenarios: determines capacities of original arcs by random selection according to the outage probabilities of the arcs (i.e. a weighted coin flip). Note that sampling scenarios APPENDS the capacities as an entirely separate set within the data structure. If one wishes to use a new sample set, s/he should run clear_scenarios first.
		
		specify_params: initializes parameters used in forecast determination. Each parameter is well document as to how it's used.


testADP.py: this file serves as the procedure that performs the entire analysis discussed in the thesis. It initializes components using tuneable parameters, trains the model, generates the sets of forecast and outage data, executes IPs using scenario configurations at time 0 and time T, and finally runs policy performance with regards to LPs configured with these generated scenarios. Performance data is organized per training iteration per scenario per time step. Because of the large dimensionality of the performance calculations, performance for training iterations can be reduced to a regularly spaced set of iterations (modulo by some user-specified factor), and for the purposes of visualization the average performance across scenarios is observed, although other metrics can certainly be introduced with ease (an operation performed in dataReports.py). Finally, only the policy performance is observed at the last time step since the goal of the algorithm is to provide an entire temporal sequence of installations. The following functions exist within the class:
		generate_sample_data: produces the set of scenarios, each scenario containing a multiple groups of outage % and randomly determined capacity based on those percentages for all arcs. The set of scenarios can be saved or loaded; further, the scenario at time 0 can be fixed across all scenarios.
		
		run_IPs_per_forecast: arc installation decisions are determined at a given time step. Installation decisions are determined per scenario, each scenario's implementation being a sample average approximation over multiple outage samples (i.e. multiple sets of constraints corresponding to different combinations of outages); the solutions are then aggregated into a mean performance later.
		
		run_LPs_per_forecast: goes through training iterations, evaluating the VFA using the generated scenarios, then evaluates its performance when subjected to an LP. This is done per scenario, each scenario's implementation being a sample average approximation over multiple outage samples (i.e. multiple sets of constraints corresponding to different combinations of outages); the solutions are then aggregated into a mean performance later.
		
		run_one_big_IPSAA: determines one set of arc installations by running an IP using a sample average approximation in the control step over ALL outage samples for every sampled scenario that requires all constraints be met simultaneously.
		
		run_one_big_LPSAA: as run_LPs_per_forecast, but executes a single LP using a sample average approximation construction as used in run_one_big_IPSAA.
		
		



Pickle files:
an200_013_40_tests_60_samples.pickle: a saved set of scenarios generated by generate_sample_data which holds 40 forecast timelines, each timeline containing 60 sets of outage outcomes.

an200_013_network_data_tnet_arcwt0.pickle: network data structure corresponding to the 200 node random network described in the thesis.

an200_013_state_data_0p03eps_0p03mcclain_3T_301N.pickle: saved scenario and outage outcome data on the random network for 301 training iterations using "explorepolicy" of .03, "mcclainss" of .03, and a 3 time step horizon.

an200_013_test_samples.pickle: a saved set of scenarios for the random network, generated by generate_sample_data which holds 20 forecast timelines, each timeline containing 50 sets of outage outcomes.

NewC_network_data.pickle: network data structure corresponding to the CLARC network described in the thesis.

NHC_network_data.pickle: network data structure corresponding to the New Hanover County network described in the thesis.

NHC_state_data_0p03eps_0p03mcclain_3T_201N.pickle.pickle: saved scenario and outage outcome data on the New Hanover County network for 201 training iterations using "explorepolicy" of .03, "mcclainss" of .03, and a 3 time step horizon.

NHC_test_samples.pickle: a saved set of scenarios for the New Hanover County network, generated by generate_sample_data which holds 20 forecast timelines, each timeline containing 50 sets of outage outcomes.



Text files:
NewC_Master_Arc_Query.txt: extracted arc data from CLARC network's master database, provided by Ryan Loggins. Used in parsers.py and subsequently fed into a networkGenerator object.

NewC_Nodes_Basic_Attributes_Query.txt: extracted node data from CLARC network's master database, provided by Ryan Loggins. Used in parsers.py and subsequently fed into a networkGenerator object.

NHC_Master_Arc_Query.txt: extracted arc data from New Hanover County network's master database, provided by Burak Cavdoroglu and Ryan Loggins. Used in parsers.py and subsequently fed into a networkGenerator object.

NHC_Nodes_Basic_Attributes_Query.txt: extracted node data from New Hanover County network's master database, provided by Burak Cavdoroglu and Ryan Loggins. Used in parsers.py and subsequently fed into a networkGenerator object.

readme.txt: this document.



PDF files:
notes 2012-09-27 Mitchell.pdf: see alt_* functions in exeADP.py.

notes 2012-10-12 Self.PDF: see "node_constr" parameters of init_model_components function in exeADP.py.