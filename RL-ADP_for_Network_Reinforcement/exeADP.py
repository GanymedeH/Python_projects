import random
import sys
import time
import copy
import cplex
import networkx as nx
import networkGenerator as ngen
import simulationEngine as sim
import numpy as np



def qpick(fname, mode="r", data=None):
    import cPickle
    f=open(fname,mode)
    if mode == "r":
        data = cPickle.load(f)
    elif mode == "w":
        cPickle.dump(data, f)
    f.close()
    return data






def convert_cplex_structs(model=None, cpobj=None):
    if type(cpobj) == type(int()):
        cpobj = (model.variables.get_names())[cpobj]
    if "SlackForNode" in cpobj:
        cpobj = cpobj.lstrip("SlackForNode")
        # this is a hack to pare off the scenario label
        cpobj = cpobj.split("S")[0]
        return int(cpobj)
    else:
        fromnode = int(cpobj.split('(')[1].split(',')[0])
        tonode = int(cpobj.split(')')[0].split(',')[-1])
        return tuple([fromnode, tonode])








#TODO8: this is really slow.....

# this function should only be run once both temp and original arcs have been established in the model; otherwise,
#    errors can arise in the form of infeasibilities.
def update_models_constraints(mC, models, netobj=None, mrange=[], useSlack=True):
    import cplex

    if type(models) != type([0]):
        print "update_models_constraints: takes as an argument a list of models, not just a model."
        return
    if netobj == None:
        netobj = mC["netobj"]

    for m in mrange:
        # figure out how many scenarios we have and establish a set of constraints for each
        for k in xrange(mC["numscens"]):
            ovnSlice = mC["orig_var_names"][k*netobj.dinet.number_of_edges(): \
                                            (k+1)*netobj.dinet.number_of_edges() ]
            if useSlack:
                svnSlice = mC["slak_var_names"][k*len(netobj.demands): (k+1)*len(netobj.demands) ]
            for i,q in netobj.dinet.nodes_iter(data=True):
                # supply constraints
                if q["type"] == "supply":
                    supind = [s for s in ovnSlice if "OrigArc("+str(i)+"," in s] + \
                             [s for s in mC["temp_var_flows"] if "TempFlow("+str(i)+"," in s]
                    # recall that supply was stored as negative demand
                    models[m].linear_constraints.add(lin_expr = [cplex.SparsePair(ind = supind, val = [1]*len(supind)) ],\
                                             senses   = ["L"] , rhs = [-1*q["demand"] ], \
                                             names    = ["SConstr"+str(i)+"Scen"+str(k)] )
                # demand constraints
                elif q["type"] == "demand":
                    if useSlack:
                        demind = [d for d in ovnSlice if ","+str(i)+")Scen" in d] + \
                                 [d for d in mC["temp_var_flows"] if "TempFlow(" in d and ","+str(i)+")" in d] + \
                                 [d for d in svnSlice if "SlackForNode"+str(i)+"Scen" in d]
                    else:
                        demind = [d for d in ovnSlice if ","+str(i)+")Scen" in d] + \
                                 [d for d in mC["temp_var_flows"] if "TempFlow(" in d and ","+str(i)+")" in d]
                    models[m].linear_constraints.add(lin_expr = [cplex.SparsePair(ind = demind, val = [1]*len(demind)) ],\
                                             senses   = ["E"] , rhs = [q["demand"]], \
                                             names    = ["DConstr"+str(i)+"Scen"+str(k)] )
                # flow balance constraints
                elif q["type"] == "trans":
                    outflow = [f for f in ovnSlice if "OrigArc("+str(i)+"," in f] + \
                              [f for f in mC["temp_var_flows"] if "TempFlow("+str(i)+"," in f]
                    inflow  = [f for f in ovnSlice if ","+str(i)+")Scen" in f] + \
                              [f for f in mC["temp_var_flows"] if "TempFlow(" in f and ","+str(i)+")" in f]
                    models[m].linear_constraints.add(lin_expr = [cplex.SparsePair(ind = outflow + inflow, \
                                                                          val = [1]*len(outflow) + \
                                                                                [-1]*len(inflow)) ], \
                                               senses = ["E"], rhs = [0], \
                                               names  = ["FConstr"+str(i)+"Scen"+str(k)] )

        # limit any newly decided flow to its capacity if installed
        models[m].linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [mC["temp_var_insts"][i], \
                                                                     mC["temp_var_flows"][i]], \
                                                              val = [mC["temp_var_caps"][i], -1] ) for i in \
                                                                             xrange(len(mC["temp_var_insts"])) ], \
                                 senses   = ["G"]*len(mC["temp_var_insts"]), \
                                 rhs      = [0]*len(mC["temp_var_insts"]), \
                                 names    = ["TPConstr"+str(i) for i in xrange(len(mC["temp_var_insts"])) ] )
    #TODO9: can only install a certain number of arcs total (i.e. over all time steps)
    return







# only the temp variables have any meaning; the original and slack variables have upper bounds of 0 so that they
#    can be "earmarked" in the model, but require updating for the model to be universally feasible. the original
#    arcs and slack variables are typically updated in update_models_orig_scens.
def create_models_temp_framework(tsteps, mC):
    models = [cplex.Cplex() for i in xrange(tsteps)]
    for t in xrange(tsteps):
        models[t].set_results_stream(None)
        models[t].set_log_stream(None)              #open("cplex_warmstart_time"+str(t)+".log","w"))
        models[t].set_warning_stream(None)          #open("cplex_iter"+str(n)+"_time"+str(t)+"_warns.txt","w"))


        #NOTE: we are NOT including original arc or slack variables, as they are scenario dependent
        models[t].variables.add(obj   = mC["inst_var_costs"] + mC["flow_var_costs"],  \
                                lb    = [0]*len(mC["temp_var_insts"] + mC["temp_var_flows"]), \
                                ub    = [1]*len(mC["temp_var_insts"]) + mC["temp_var_caps"],  \
                                types = "B"*len(mC["temp_var_insts"]) + "C"*len(mC["temp_var_flows"]), \
                                names = mC["temp_var_insts"] + mC["temp_var_flows"])
    return models






def find_constr_4_vars(model, variables, vartype="edge"):
    allnames = model.variables.get_names()
    allconstr = model.linear_constraints.get_names()
    def append_constr(j, model, allnames, decisions):
        constrdecs = ""
        for k in xrange(len(model.linear_constraints.get_rows(j).ind)):
            constrdecs += str(model.linear_constraints.get_rows(j).val[k]) + "*[" + \
                allnames[model.linear_constraints.get_rows(j).ind[k]] + "=" + \
                str(decisions[model.linear_constraints.get_rows(j).ind[k]]) + "] + "
        constrdecs = constrdecs.rstrip(" +")
        consense = model.linear_constraints.get_senses(j)
        if consense == "G":
            constrdecs += " >= "
        elif consense == "L":
            constrdecs += " <= "
        elif consense == "E":
            constrdecs += " = "
        else:
            constrdecs += " R "
        constrdecs += str(model.linear_constraints.get_rhs(j))
        return constrdecs

    if model.solution.is_primal_feasible():
        decisions = model.solution.get_values()
    else:
        decisions = [None]*len(allnames)
    constrs = []
    for i in xrange(len(variables)):
        if type(variables[i]) == type(str()) and vartype == "edge":
            varidx = [j for j in xrange(len(allnames)) if variables[i] == allnames[j] ][0]
        elif type(variables[i]) == type(str()) and vartype == "constraint":
            varidx = [j for j in xrange(len(allconstr)) if variables[i] == allconstr[j] ][0]
        else:
            varidx = variables[i]
        if vartype == "edge":
            constri = []
            for j in xrange(model.linear_constraints.get_num()):
                if varidx in model.linear_constraints.get_rows(j).ind:
                    constri.append( append_constr(j, model, allnames, decisions) )
        elif vartype == "constraint":
            constri = append_constr(varidx, model, allnames, decisions)
        constrs.append((variables[i], constri))
    return constrs












def execute_ADP(tsteps, bigN, mC, features, recout="recvs", delout="delvs", savestate=None, usestate=None):
    # these are for statistical reporting purposes (change in v and historic v over iterations)
    delv = []
    recv = []
    if "sp_" in mC["bftype"]:
        basisname = "_sp_basis_val"
    elif "ones" in mC["bftype"]:
        basisname = "_1s_basis_val"
    elif "node_constr" in mC["bftype"]:
        basisname = "_nc_basis_val"
    else:
        print "execute_ADP: invalid basis function specification."
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    if savestate != None:
        iterstates = [ dict() for n in xrange(bigN)]
    if usestate != None:
        usestate = qpick(usestate)

    ## BEGIN STEP 0
    #NOTE: warm start isn't a good idea when you have no idea what to expect, so turn it off for now
##    adp_warm_start(mC, tsteps, features)
    ## END STEP 0

    ## BEGIN STEP 1
    # have initial v's, generate new scenario for each time step (from newly generated distributions!)
    #     and repeat process immediately above

    #DEBUG begin
    proctimedbg = open("cplex_solve_time_trials.txt","w")
    mla = 0
    #DEBUG end
    print "main algorithm begin"
    for n in xrange(bigN):
        #DEBUG begin
        mlt = time.time()
        oldv = copy.deepcopy(mC["theta"])
        #DEBUG end
        sys.stdout.write("Iteration "+str(n+1)+": \n")

        # create the "core" model: a general framework that is common to all cplex executions as they relate to ADP
        # since the objective depends on the policy method, this just adds the temp variables to the models
        models = create_models_temp_framework(tsteps, mC)
        if "ud_vfa_separate" in mC["polmethod"]:
            # there are 2 sets of models: 1 performs network flow and the other finds 'discounts' to the flow
            spdmodels = create_models_temp_framework(tsteps, mC)

        ## BEGIN STEP 1.1/1.2
        # create damage distributions for each time step. we change the distributions every iteration as well,
        #    because the v's must learn on all possible forecasts simultaneously to give fair representation
        #    over all features and to prevent bias occurring between algorithm executions: v's would learn based
        #    off of scenarios for only 1 distribution which specializes it to that distribution (no matter how
        #    it was generated). running the adp algorithm again would produce an entirely different policy no
        #    matter how the forecasts are simulated, and the goal here is to have a single policy whatever
        #    behavior we encounter.
        if usestate == None:
            mC["mcsim"].generate_distributions(mC["netobj"].dinet, tsteps)
            mC["mcsim"].eval_adp_basis_functions(mC["netobj"], tsteps, mC)
        else:
            for t in xrange(tsteps):
                nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(t)+"_dmg_pct", usestate[n]["t"+str(t)+"_dmg_pct"])
                nx.set_edge_attributes(mC["netobj"].tnet, "t"+str(t)+basisname, \
                                       usestate[n]["t"+str(t)+basisname])
        # determine which features are present in this iteration.
        #NOTE: we may have to move this depending on how features depend on scope
        if "ctr_mass_quad_rel_frame" in features:
            featurelist = mC["mcsim"].calc_feature_ctrmassquad(mC["netobj"], tsteps, features)
        ## END STEP 1.1/1.2

        ## BEGIN STEP 1.3
        #DEBUG begin
        procadpla = 0
        #DEBUG end

        # we're not performing SAA here, so make sure we only have 1 sample; however, even though there is only one
        #    sample, the model is now technically a SAA (where the sample size is 1).
        mC["mcsim"].clear_scenarios(mC["netobj"].dinet)
        if usestate == None:
            mC["numscens"] = mC["mcsim"].sample_scenarios(mC["netobj"].dinet, ksamples=1, trange=xrange(tsteps))
        else:
            for t in xrange(tsteps):
                nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(t)+"_capacity", \
                                       usestate[n]["t"+str(t)+"_capacity"])
        if savestate != None:
            for t in xrange(tsteps):
                iterstates[n]["t"+str(t)+"_dmg_pct"] = nx.get_edge_attributes(mC["netobj"].dinet, \
                                                                              "t"+str(t)+"_dmg_pct")
                iterstates[n]["t"+str(t)+basisname] = nx.get_edge_attributes(mC["netobj"].tnet, \
                                                                           "t"+str(t)+basisname)
                iterstates[n]["t"+str(t)+"_capacity"] = nx.get_edge_attributes(mC["netobj"].dinet, \
                                                                               "t"+str(t)+"_capacity")

        if "ud_vfa_combined" in mC["polmethod"] or "node_constr" in mC["polmethod"]:
            update_models_orig_scens(models, mC, mrange=xrange(len(models)-1), reset=True, updateConstr=True)
        # mC["polmethod"] ud_vfa_separate is a non-network oriented IP, thus doesn't use network constraints
        elif "ud_vfa_separate" in mC["polmethod"]:
            # there are 2 sets of models: 1 performs network flow and the other finds 'discounts' to the flow
            #    because of this, we need to update the model and run it before we update the network flow model
            # there are no constraints for the vfa optimization so multiple scenarios is moot, hence upd all tsteps
            update_models_orig_scens(models, mC, numscens=1, mrange=xrange(len(models)), reset=True, \
                                     updateConstr=False)
        else:
            print "execute_ADP: invalid policy method specified."
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        mC["instchoices"] = [ [0]*len(mC["temp_var_insts"]) for i in xrange(tsteps)]
        mC["previnsts"] = [[] for i in xrange(tsteps)]
        retpolicy = [ [0 for j in mC["temp_var_insts"]] for i in xrange(tsteps)]
        slacks = [None for i in xrange(tsteps)]
        baseobjvals = [None for i in xrange(tsteps)]
        vfavals = [None for i in xrange(tsteps)]
        thetaupdvfavals = [None for i in xrange(tsteps)]

        # t=tsteps-1 is special because it could be a standard flow model as used in the warm start; if so,
        #    the general policy_search can be invoked normally just by setting v's to 0
        for t in xrange(tsteps-1):
            #DEBUG begin
##                sys.stdout.write("Time "+str(t)+" ")
            #DEBUG end
            retpolicy[t],tmpflows,slacks[t],baseobjvals[t],vfavals[t] = \
                adp_policy_search(t, featurelist, models, mC, randpol=True, dbout=True, timedbg=proctimedbg, \
                                  adpla=procadpla)
            if "ud_vfa_combined" in mC["polmethod"] or "node_constr" in mC["polmethod"]:
                if "_pass_udterm" in mC["polmethod"]:
                    thetaupdvfavals[t] = baseobjvals[t]
                elif "_pass_vfaterm" in mC["polmethod"]:
                    thetaupdvfavals[t] = vfavals[t]
                elif "_pass_bothterms" in mC["polmethod"]:
                    thetaupdvfavals[t] = baseobjvals[t] + vfavals[t]
                else:
                    print "execute_ADP: invalid passed parameter specification in combined policy search method."
                    __dbg = raw_input("execution halted, press ENTER to exit")
                    sys.exit()

        # now perform for the last time step -- it's separate in case we wish to run an SAA model (which can only
        #    be used if we're not using a saved state)
        if usestate == None:
            # a complete hack... to allow for time T to be an SAA, need to modify numscens (it's used in objective)
            mC["numscens"] = mC["timeTSAAscens"]
            mC["mcsim"].clear_scenarios(mC["netobj"].dinet)
            mC["mcsim"].sample_scenarios(mC["netobj"].dinet, ksamples=mC["timeTSAAscens"], trange=[tsteps-1])

        if "ud_vfa_combined" in mC["polmethod"] or "node_constr" in mC["polmethod"]:
            update_models_orig_scens(models, mC, numscens=mC["numscens"], mrange=[len(models)-1], reset=True, \
                                     updateConstr=True)

        retpolicy[tsteps-1],tmpflows,slacks[tsteps-1],baseobjvals[tsteps-1],vfavals[tsteps-1] = \
            adp_policy_search(tsteps-1, featurelist, models, mC, randpol=True, dbout=True, timedbg=proctimedbg, \
                              adpla=procadpla)
        if "ud_vfa_combined" in mC["polmethod"] or "node_constr" in mC["polmethod"]:
            if "_pass_udterm" in mC["polmethod"]:
                thetaupdvfavals[tsteps-1] = baseobjvals[tsteps-1]
            elif "_pass_vfaterm" in mC["polmethod"]:
                thetaupdvfavals[tsteps-1] = vfavals[tsteps-1]
            elif "_pass_bothterms" in mC["polmethod"]:
                thetaupdvfavals[tsteps-1] = baseobjvals[tsteps-1] + vfavals[tsteps-1]
            else:
                print "execute_ADP: invalid updmethod subtype specified for given polmethod."
                __dbg = raw_input("execution halted, press ENTER to exit")
                sys.exit()
        elif "ud_vfa_separate" in mC["polmethod"]:
            # now run the network flow problem to obtain the 2nd component of the vfaval
            # another complete hack... to allow for time t to be an SAA, need to modify numscens (used in obj)
            #NOTE: we're using timeTSAAscens for all tsteps since we're looking for the most accurate flow cost
            baknumscens = mC["numscens"]
            mC["numscens"] = mC["timeTSAAscens"]
            mC["mcsim"].clear_scenarios(mC["netobj"].dinet)
            mC["numscens"] = mC["mcsim"].sample_scenarios(mC["netobj"].dinet, ksamples=mC["timeTSAAscens"], \
                                                          trange=xrange(tsteps))
            update_models_orig_scens(spdmodels, mC, numscens=mC["numscens"], mrange=xrange(len(spdmodels)), \
                                     reset=True, updateConstr=True)
            for t in xrange(tsteps):
                spdmodels[t].variables.set_types([ (mC["temp_var_insts"][j], spdmodels[t].variables.type.continuous) \
                                               for j in xrange(len(mC["temp_var_insts"])) ])
                spdmodels[t].variables.set_lower_bounds( \
                    [ (mC["temp_var_insts"][i], retpolicy[t-1][i]) for i in xrange(len(mC["temp_var_insts"])) ])
                spdmodels[t].variables.set_upper_bounds( \
                    [ (mC["temp_var_insts"][i], retpolicy[t-1][i]) for i in xrange(len(mC["temp_var_insts"])) ])
                spdmodels[t].set_problem_type(spdmodels[t].problem_type.LP)
                try:
                    spdmodels[t].solve()
                    if "grad_" in mC["updmethod"] or "least_squares_" in mC["updmethod"]:
                        # we add vfavals because it's already been negated within adp_policy_search. we force it to
                        #    be nonpositive however, since it makes little practical sense to allow for incurred cost
                        #    due to attempting to alleviate unmet demand. further, we force the discounted cost,
                        #    thetaupdvfavals, to be nonnegative since it doesn't make sense to receive money for
                        #    meeting demand.
                        thetaupdvfavals[t] = max(spdmodels[t].solution.get_objective_value() + min(vfavals[t], 0), 0)
                    elif "min_norm_" in mC["updmethod"]:
                        # this method doesn't use the vfa decisions directly
                        thetaupdvfavals[t] = spdmodels[t].solution.get_objective_value()
                    else:
                        print "execute_ADP: invalid updmethod specified for polmethod ud_vfa_separate."
                        __dbg = raw_input("execution halted, press ENTER to exit")
                        sys.exit()
                except:
                    print "execute_ADP step 1.3 using separated policy method: no solution exists in flow problem"
            mC["numscens"] = baknumscens
        ## END STEP 1.3

        ## BEGIN STEP 1.4
        trueudcost, trueinstdecs = adp_true_val_fn_calc(mC, tsteps, retpolicy)
        if "_withvalfn_zeta" in mC["updmethod"]:
            for t in xrange(tsteps):
                retpolicy[t] = [2*trueinstdecs[t][i]-retpolicy[t][i] for i in xrange(len(retpolicy[t])) ]
        ## END STEP 1.4

        ## BEGIN STEP 1.5
        if "lp_toggle" in mC["updmethod"]:
            pbla = adp_vfa_update_lp_toggle(tsteps, models, mC, baseobjvals, featurelist)
        elif "grad_" in mC["updmethod"]:
            # if we have node-based basis functions, we need to re-evaluate them based on the new decisions
            if "node_constr" in mC["polmethod"]:
                mC["mcsim"].eval_adp_basis_functions(mC["netobj"], tsteps, mC, constrdecs=retpolicy, slackdecs=slacks)
            adp_vfa_update_grad_descent(trueudcost, featurelist, tsteps, mC, vfaobj=thetaupdvfavals, decvars=retpolicy)
        elif "least_squares_" in mC["updmethod"]:
            adp_vfa_update_least_squares(trueudcost, featurelist, tsteps, mC, thetaupdvfavals, decvars=retpolicy)
        elif "min_norm_" in mC["updmethod"]:
            adp_vfa_update_min_norm(trueudcost, featurelist, tsteps, mC, thetaupdvfavals, decvars=trueinstdecs)
        else:
            print "execute_ADP: invalid updmethod specified."
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        #DEBUG begin
##        proctimedbg.write("       perturbation loop avg time: "+\
##                      str(pbla/float(len(mC["temp_var_insts"])))+" sec\n")
        #DEBUG end
        ## END STEP 1.5

        # garbage collection
        trash_model(models)
        if "ud_vfa_separate" in mC["polmethod"]:
            trash_model(spdmodels)

        #DEBUG begin
        mla += time.time() - mlt
        proctimedbg.write(" main loop avg time: "+str(mla/float(n+1))+" sec\n")
        #DEBUG end

        delv.append([ [ [mC["theta"][t][f][i]-oldv[t][f][i] for i in \
                         xrange(len(mC["theta"][t][f])) ] for f in xrange(len(mC["theta"][t])) ] \
                      for t in xrange(len(mC["theta"])) ])
        recv.append(copy.deepcopy(mC["theta"]))
    proctimedbg.flush()
    proctimedbg.close()

    qpick(delout+".pickle","w",delv)
    qpick(recout+".pickle","w",recv)
    if savestate != None:
        qpick(savestate+".pickle","w",iterstates)
    return









def adp_true_val_fn_calc(mC, tsteps, boundpolicy):
    trueudcost = [None for t in xrange(tsteps)]
    trueinstdecs = [None for t in xrange(tsteps)]
    # formulate IP with last time step's forecast to determine "true" unmet demand
    lmodels = create_models_temp_framework(tsteps, mC)
    mC["mcsim"].clear_scenarios(mC["netobj"].dinet)
    mC["mcsim"].sample_scenarios(mC["netobj"].dinet, ksamples=mC["Lscens"], trange=[tsteps-1])
    # the true value function uses time tsteps-1's forecast data for all time steps no matter what
    for t in xrange(tsteps-1):
        for af,at,ad in mC["netobj"].dinet.edges_iter(data=True):
            ad["t"+str(t)+"_capacity"] = ad["t"+str(tsteps-1)+"_capacity"]

    update_models_orig_scens(lmodels, mC, numscens=mC["Lscens"], mrange=xrange(len(lmodels)), reset=True, \
                             updateConstr=True)
    if mC["truevalfn"] == "IPv1a" or mC["truevalfn"] == "IPv1b":
        lmodels[-1].linear_constraints.add( \
                    lin_expr = [cplex.SparsePair(ind = mC["temp_var_insts"], val = [1]*len(mC["temp_var_insts"]) ) ], \
                                senses = ["L"], rhs = [mC["pertinstbudget"]*tsteps], names = ["PerTInstConstr"] )
        if mC["truevalfn"] == "IPv1b":
            # this is essentially lmodels[tsteps-1] for IPv2a/b, it just uses it for all time steps
            lmodels[-1].variables.set_lower_bounds(zip(mC["temp_var_insts"], boundpolicy[tsteps-2] ) )
    elif mC["truevalfn"] == "LPv1":
        lmodels[tsteps-1].variables.set_types([ (mC["temp_var_insts"][j], \
                                                lmodels[tsteps-1].variables.type.continuous) \
                                       for j in xrange(len(mC["temp_var_insts"])) ])
        lmodels[tsteps-1].variables.set_lower_bounds( zip(mC["temp_var_insts"], boundpolicy[tsteps-1]) )
        lmodels[tsteps-1].variables.set_upper_bounds( zip(mC["temp_var_insts"], boundpolicy[tsteps-1]) )
        lmodels[tsteps-1].set_problem_type(lmodels[tsteps-1].problem_type.LP)
    elif "IPv2" in mC["truevalfn"]:
        # 2 ways for time-dependent TVFs: allow installs of only a time step's budget of arcs, or the remainder of
        #    the total cumulative budget. note that at time T they produce the same objective value
        if mC["truevalfn"] == "IPv2a":
            for t in xrange(tsteps):
                lmodels[t].linear_constraints.add( \
                        lin_expr = [cplex.SparsePair(ind = mC["temp_var_insts"], \
                                                     val = [1]*len(mC["temp_var_insts"]) ) ], \
                                    senses = ["L"], rhs = [mC["pertinstbudget"]*(t+1)], names = ["PerTInstConstr"] )
        elif mC["truevalfn"] == "IPv2b":
            for t in xrange(tsteps):
                lmodels[t].linear_constraints.add( \
                        lin_expr = [cplex.SparsePair(ind = mC["temp_var_insts"], \
                                                     val = [1]*len(mC["temp_var_insts"]) ) ], \
                                    senses = ["L"], rhs = [mC["pertinstbudget"]*tsteps], names = ["PerTInstConstr"])
        for t in xrange(1,tsteps):
            # when accounting for previously installed arcs, time 0 has none
            lmodels[t].variables.set_lower_bounds(zip(mC["temp_var_insts"], boundpolicy[t-1] ) )
    else:
        "adp_true_val_fn_calc: invalid value function type specified."
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()
    try:
        if mC["truevalfn"] == "IPv1a" or mC["truevalfn"] == "IPv1b" or mC["truevalfn"] == "LPv1":
            lmodels[-1].solve()
            trueudcost = [lmodels[-1].solution.get_objective_value() for t in xrange(tsteps)]
            trueinstdecs = [ [int(round(i)) for i in lmodels[-1].solution.get_values(mC["temp_var_insts"])] \
                             for t in xrange(tsteps)]
        elif mC["truevalfn"] == "IPv2a" or mC["truevalfn"] == "IPv2b" or mC["truevalfn"] == "LPv2":
            for t in xrange(tsteps):
                lmodels[t].solve()
                trueudcost[t] = lmodels[t].solution.get_objective_value()
                trueinstdecs[t] = [int(round(i)) for i in lmodels[t].solution.get_values(mC["temp_var_insts"])]
    except:
        print "execute_ADP step 1.4: no solution to true value function exists"
    trash_model(lmodels)
    # every time we update_models_orig_scens it overwrites mC["numscens"] for its own purposes - reset to be safe
    mC["numscens"] = 1
    return trueudcost, trueinstdecs







def update_models_orig_scens(models, mC, netobj=None, numscens=1, mrange=[], reset=False, updateConstr=True):
    # account for resets: delete existing variables
    if netobj == None:
        netobj = mC["netobj"]

    if reset:
        for t in mrange:
            try:
                models[t].variables.delete(mC["orig_var_names"] + mC["slak_var_names"])
            except:
                pass
        mC["orig_var_names"] = []
        mC["orig_var_caps"]  = []
        mC["orig_var_costs"] = []
        mC["slak_var_names"] = []
        mC["slak_var_caps"]  = []
        mC["slak_var_costs"] = []
        mC["numscens"] = 0

    addOVNs = []
    # the indexing order is ovu[k][t][i] (scen#, time, arcID)
    addOVUs = [ [] for t in xrange(len(models))]
    addOVCs = []
    addSVNs = []
    addSVUs = []
    addSVCs = []

    for i in xrange(mC["numscens"],mC["numscens"]+numscens):
        #NOTE: can i add these to the model all at once? if so, use flattened model and 'extend', otherwise need to
        #      put the cplex calls inside this loop
        addOVNs.extend(["OrigArc("+str(j)+","+str(k)+")Scen"+str(i) for j,k in netobj.dinet.edges_iter() ])
        try:
            for t in mrange:
                addOVUs[t].extend([d["t"+str(t)+"_capacity"][i] for k,l,d in \
                                                                         netobj.dinet.edges_iter(data=True)])
        except:
            print "update_models_orig_scens: either invalid time or sample value given, exiting"
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        addOVCs.extend([d["weight"] for j,k,d in netobj.dinet.edges_iter(data=True) ])
        svnextend = ["SlackForNode"+str(j)+"Scen"+str(i) for j,k in netobj.dinet.nodes_iter(data=True) if \
                                                                  k["type"] == "demand"]
        addSVNs.extend(svnextend)
        #ASSUMPTION: slack capacities and costs are the same across sample scenarios
        addSVUs.extend([mC["slackcap"]]*len(svnextend))
        addSVCs.extend([mC["slackcost"]]*len(svnextend))

    mC["numscens"] += numscens
    # now the arduous task of updating the cplex side of the model
    for m in mrange:
        # step 1: update objective to reflect a new SAA cost coefficient in the objective (if original objective
        #         still has the variables)
        if not reset:
            objOVCcoeff = models[m].objective.get_linear(mC["orig_var_names"])
            objSVCcoeff = models[m].objective.get_linear(mC["slak_var_names"])
            models[m].objective.set_linear([ (mC["orig_var_names"][j], (1./mC["numscens"])*objOVCcoeff[j]) for j in\
                                             xrange(len(mC["orig_var_names"])) ] + \
                                           [ (mC["slak_var_names"][j], (1./mC["numscens"])*objSVCcoeff[j]) for j in\
                                             xrange(len(mC["slak_var_names"])) ] )

        # step 2: add the new variables using the updated cost in the objective for their contribution
        #NOTE: from this point on, indexing the solution vector by index number can have unintended effects!
        models[m].variables.add(obj   = [(1./mC["numscens"])*i for i in addOVCs ] + \
                                         [(1./mC["numscens"])*i for i in addSVCs],  \
                                lb    = [0]*len(addOVNs + addSVNs), \
                                ub    = addOVUs[m] + addSVUs, \
                                types = "C"*len(addOVNs + addSVNs), \
                                names = addOVNs + addSVNs)

    # step 3: update the components dictionary (mostly for bookkeeping purposes at this point)
    mC["orig_var_names"].extend(addOVNs)
    # just stick the values in as compiled: OVU[k][t][i], we'll see if i need them properly sorted later
    mC["orig_var_caps"].extend(addOVUs)
    mC["orig_var_costs"].extend(addOVCs)
    mC["slak_var_names"].extend(addSVNs)
    mC["slak_var_caps"].extend(addSVUs)
    mC["slak_var_costs"].extend(addSVCs)

    # step 4: update linear constraints to reflect new scenario restrictions
    if updateConstr:
        update_models_constraints(mC, models, netobj, mrange)
    return




#TODO8: this is really slow.....
def trash_model(models):
    for t in xrange(len(models)):
        models[t].linear_constraints.delete()
        models[t].variables.delete()
        models[t].cleanup(sys.maxint-1)
    del models
    return






def init_model_components(tsteps, slackcost, pertinstbudget, numfeats, poleps, netobj=None, nxnet=None, tnet=None, \
                          polmethod="ud_vfa_combined", updmethod="grad_linear", truevalfn="IPv1", \
                          bftype="sp_probabilistic", tfc=0, ttsaas=1, lscens=1, saveNet=True):
    modelComponents = {}
    # from scratch
    if nxnet == None and (netobj == None or not any(netobj.dinet)):
        print "init_model_components: Pass in a network first!"
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()
    elif netobj == None:
        netobj = ngen.networkGenerator(dn=nxnet, tn=tnet)
    elif netobj.dinet == None:
        netobj.dinet = nxnet

    #DEPRECATED
##    netobj.set_param("forcetermini",False)


    # 1 way to create temp arcs is to simply provide arcs between nodes that don't already exist; i.e., any 0
    # entry in the incidence matrix is a temp arc with capacity equal to the sum of all supply (or sum of all
    # demand if it's different and higher than the supply amount). temp arc cost ("weight") is defaulted to 1.
    ##    tnet = ngen.complement_arc_capacities(nxnet, [netobj.get_param("sub")]*netobj.get_param("snodes"), \
    ##                                            [netobj.get_param("dub")]*netobj.get_param("dnodes") )
    # a more realistic way: temp arcs only between hubs, a few relays, and from suppliers to nearby hubs/relays
    if tnet == None and netobj.tnet == None:
        print "Generating temp network... "
        netobj.tnet = netobj.create_temp_net()
    elif netobj.tnet == None:
        netobj.tnet = tnet

    if saveNet:
        netobj.save_network()

    #ASSUMPTION: initialization means you start from scratch, which means there aren't any scenarios!
    #NOTE: i have to use edges_iter for everything because a straight value dump from the dict doesn't sort indices
    modelComponents["numscens"] = 0
    modelComponents["netobj"]   = netobj
    modelComponents["mcsim"]    = sim.simulationEngine()

    modelComponents["slackcap"] = sys.maxint-1
    modelComponents["slackcost"] = slackcost
    modelComponents["pertinstbudget"] = pertinstbudget
    modelComponents["numfeats"] = numfeats
    modelComponents["explorepolicy"] = poleps
    modelComponents["polmethod"] = polmethod
    modelComponents["updmethod"] = updmethod
    modelComponents["truevalfn"] = truevalfn
    modelComponents["bftype"] = bftype
    #NOTE: if using a saved state, ttsaas must be 1!
    modelComponents["timeTSAAscens"] = ttsaas
    modelComponents["Lscens"] = lscens

    # there's so many vectorizations to do, it's actually quicker to append in a single loop
    modelComponents["temp_var_insts"] = []
    modelComponents["inst_var_costs"] = []
    modelComponents["orig_var_names"] = []
    modelComponents["orig_var_caps"]  = []
    modelComponents["orig_var_costs"] = []
    modelComponents["temp_var_flows"] = []
    modelComponents["temp_var_caps"]  = []
    modelComponents["slak_var_names"] = []

    for i,j,d in netobj.tnet.edges_iter(data=True):
        modelComponents["temp_var_insts"] += ["TempArc("+str(i)+","+str(j)+")"]
        modelComponents["inst_var_costs"] += [d["weight"]]
        modelComponents["temp_var_flows"] += ["TempFlow("+str(i)+","+str(j)+")"]
        modelComponents["temp_var_caps"]  += [d["capacity"]]

    for i,j,d in netobj.dinet.edges_iter(data=True):
        modelComponents["orig_var_names"] += ["OrigArc("+str(i)+","+str(j)+")Scen0"]
        # these will be overwritten once sampling is performed
        modelComponents["orig_var_caps"]  += [d["capacity"]]
        modelComponents["orig_var_costs"] += [d["weight"]]

    modelComponents["flow_var_costs"] = [tfc]*netobj.tnet.number_of_edges()

    modelComponents["slak_var_names"] = ["SlackForNode"+str(i)+"Scen0" for i,j in \
                                         netobj.dinet.nodes_iter(data=True) if j["type"] == "demand"]
    modelComponents["slak_var_caps"]  = [modelComponents["slackcap"]]*len(modelComponents["slak_var_names"])
    modelComponents["slak_var_costs"] = [slackcost]*len(modelComponents["slak_var_names"])

    # there are 2 major types of indexing scheme for VFA coefficients: one way provides a coefficient corresponding
    #    to each installable arc, while the other provides a coefficient corresponding to each node, which captures
    #    the behavior of node constraints in the model, rather than installation decisions (as in the former case)
    if "sp_" in bftype:
        #NOTE: v's indexing scheme is theta[t][f][i]
        #NOTE: 0 may be a controversial value for theta as it be be considered degenerate.. further consideration is
        #      needed. right now we're doing a very crude warm start that initializes theta values to be the cost of
        #      (half) total unmet demand spread out across the number of temp arcs. the last time step must have
        #      theta = 0 however (there is no future value to assess at the last stage).
        crudewarmv = float(sum(d["demand"] for n,d in modelComponents["netobj"].dinet.nodes_iter(data=True) if \
                               d["type"] == "demand")) #* slackcost <-- separate this & incorporate into basis func
        crudewarmv /= modelComponents["netobj"].tnet.number_of_edges() # * 2 <--- artifact from flow cost represen.
        #NOTE: rather, let's try total cost of unmet demand per allowable install for every arc. the rationale here
        #      is that the unmet demand flow will be spread amongst the budgeted installed arcs. also, rescale based
        #      off of expected basis value function.
##        crudewarmv /= (modelComponents["pertinstbudget"] * 0.2)
        modelComponents["theta"] = [ [ [crudewarmv for i in xrange(netobj.tnet.number_of_edges())] \
                                       for j in xrange(numfeats)] for k in xrange(tsteps)]
##        modelComponents["theta"][-1] = [ [0 for i in xrange(netobj.tnet.number_of_edges())] \
##                                         for j in xrange(numfeats)]
    elif "ones" in mC["bftype"]:
        modelComponents["theta"] = [ [ [0 for i in xrange(netobj.tnet.number_of_edges())] \
                                               for j in xrange(numfeats)] for k in xrange(tsteps)]
    elif "node_constr" in mC["bftype"]:
        #NOTE: when we have a theta for every node rather than every installable arc, we consider a warm start where
        #      thetas associated with demand nodes have the corresponding demand values evenly distributed among
        #      them, the supply thetas to have negative this even distribution, and the transshipment nodes to be 0
        #ASSUMPTION: to be continued....
        modelComponents["theta"] = [ [ [d["demand"] for n,d in \
                                        modelComponents["netobj"].dinet.nodes_iter(data=True)] \
                                       for j in xrange(numfeats)] for k in xrange(tsteps)]
    else:
        print "init_model_components: invalid basis function specification."
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    modelComponents["thetaupdcnt"] = [ [1 for i in xrange(numfeats)] for k in xrange(tsteps)]
    modelComponents["instchoices"] = [ [0]*len(modelComponents["temp_var_insts"]) for i in xrange(tsteps)]
    modelComponents["previnsts"] = [[] for i in xrange(tsteps)]

    return modelComponents






def adp_policy_search(t, featurelist, models, mC, randpol=True, dbout=False, timedbg=None, adpla=None):

    # reset from previous scenario (or initialize for n = 0)
    tmpflows = [0]*len(mC["temp_var_flows"])
    slacks = [0]*len(mC["slak_var_names"])
    #TODO8: inst var costs are really not needed anywhere at this point, consider making edge "weight" 0 in creation
    if "ud_vfa_combined" in mC["polmethod"] or "node_constr" in mC["polmethod"]:
        inst_var_costs = copy.copy(mC["inst_var_costs"])
    elif "ud_vfa_separate" in mC["polmethod"]:
        inst_var_costs = [0 for i in mC["inst_var_costs"]]
    else:
        print "adp_policy_search: invalid policy method specified."
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    objval = None
    unmetdemandcost = None

    #ASSUMPTION: costs (for everything but tmp installs & slacks) in the obj funct are 0 (see paper)
    #ASSUMPTION: model variables established in warm start
    #NOTE: to be clear, theta[t] has the data from the calculations done in step t+1
    if "sp_" in mC["bftype"]:
        basisname = "_sp_basis_val"
        edgecnt = 0
        # regarding the following for-block, it only makes sense (right?) to apply it to install-arc based basis
        #    functs if an arc has been previously installed, then with respect to the state, we must treat it as
        #    though it were an original arc. since no other original arcs are subject to the shortest path basis
        #    function (i.e. their phi value is zero) then so too is the previously installed arc's. original arcs'
        #    install costs are also zero.
        # the above rationale depends on how you want to represent the descent direction: sometimes travelling
        #    (aka modifying) the same direction (installed arc) is desirable, other times it is not, e.g. combined
        #    vs separated policy search objectives, respectively.
        for af, at, ad in mC["netobj"].tnet.edges_iter(data=True):
            if int(mC["instchoices"][t-1][edgecnt]):
                if "_cum" not in mC["bftype"]:
                    ad["t"+str(t)+basisname] = 0.
                inst_var_costs[edgecnt] = 0
            edgecnt += 1
    elif "ones" in mC["bftype"]:
        basisname = "_1s_basis_val"
    elif "node_constr" in mC["bftype"]:
        basisname = "_nc_basis_val"
    else:
        print "adp_policy_search: invalid basis function type specified."
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    # for certain setups, including unmet demand costs in the objective (at non-terminal time steps) is undesirable.
    #    we check for the flag, but leave it to the calling function to know when the terminal time step is.
    if "ud_vfa_combined" in mC["polmethod"]:
        slak_var_costs = copy.copy(mC["slak_var_costs"])
    elif "ud_vfa_separate" in mC["polmethod"]:
        slak_var_costs = [0 for i in xrange(len(mC["slak_var_costs"])) ]
    elif "node_constr" in mC["polmethod"]:
        # this is handled later, i just want to make sure a valid option has been passed in
        slak_var_costs = [None for i in xrange(len(mC["slak_var_costs"])) ]

    #ASSUMPTION: the only features being used are ctr_mass_quad_rel_frame. using any other
    #            will break the following code.
    # kind of a kluge.. the ctr_mass feature returns 1 thru 4 indicating cartesian quadrant
    #    so subtract 1 to provide the index to the list. to repeat: the only vfa param
    #    being updated is the one exhibiting the specific feature.
    # this if-block is to correctly set the objective function give the time and vfa form (sh_path or n_constr)
    if "ctr_mass_quad_rel_frame" in featurelist[t]:
        theta = mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1]
        # really confusing, so i'll try to notationally condense: part of the objective is
        #    pi[i](c[i] + phi[i]*vf_n[i] ) where phi is the basis function (sp_probabilistic, node_constr, etc),
        #    pi is the install decision (for i temp arcs), vf_n is the vfa for feature #, and c is install cost.
        basis_vals_t = [d["t"+str(t)+basisname] for i,j,d in mC["netobj"].tnet.edges_iter(data=True)]
        if "ud_vfa_combined" in mC["polmethod"] or "ud_vfa_separate" in mC["polmethod"]:
            if "_linear" in mC["updmethod"] or "lp_toggle" in mC["updmethod"]:
                vfa = [basis_vals_t[i]*theta[i] for i in xrange(len(mC["temp_var_insts"])) ]
            # the checks against div by 0 are only going to work when theta is identically zero, so overflow can occur
            #    this is only in place for the special case of the last time step, which sets all thetas to 0
            elif "_inverse" in mC["updmethod"]:
                vfa = []
                for i in xrange(len(mC["temp_var_insts"])):
                    if theta[i] != 0:
                        vfa.append(basis_vals_t[i]/theta[i])
                    else:
                        vfa.append(0)
            elif "_root_inverse" in mC["updmethod"]:
                # the inverse power term, meant to be between 1 and positive infinity, non-inclusive.
                invpowterm = 2.0
                # the following allows theta terms to be negative by taking the negative term out of the root calculation
                # note: this should be the same value as the neghack in adp_vfa_update_grad_descent
                neghack = True
                vfa = []
                for i in xrange(len(mC["temp_var_insts"])):
                    if theta[i] > 0:
                        vfa.append(basis_vals_t[i]/(theta[i]**(1/invpowterm)))
                    elif theta[i] < 0 and neghack:
                        vfa.append(-basis_vals_t[i]/((-theta[i])**(1/invpowterm)))
                    else:
                        vfa.append(0)
            elif "_neg_exp" in mC["updmethod"]:
                vfa = [basis_vals_t[i]*np.exp(-theta[i]) for i in xrange(len(mC["temp_var_insts"])) ]
            else:
                print "adp_policy_search: invalid update method specified."
                __dbg = raw_input("execution halted, press ENTER to exit")
                sys.exit()
        elif "node_constr" in mC["polmethod"]:
            # thetas correspond to nodes, so we need to convert to relate to arcs (+-theta_to - theta_from)
            #    there is only 1 situation in which the form is +theta_to-theta_from: when theta_to is a trans node
            arccoeffs = [0 for i in xrange(mC["netobj"].tnet.number_of_edges())]
            for taidx in xrange(len(mC["temp_var_insts"])):
                af,at = convert_cplex_structs(None, mC["temp_var_insts"][taidx])
                if mC["netobj"].tnet.node[at]["type"] == "trans":
                    arccoeffs[taidx] = theta[at] - theta[af]
                else:
                    arccoeffs[taidx] = -theta[at] - theta[af]
            if "_linear" in mC["updmethod"] or "lp_toggle" in mC["updmethod"]:
                vfa = [basis_vals_t[i]*arccoeffs[i] for i in xrange(len(mC["temp_var_insts"])) ]
            else:
                # implement the other methods later, as needed
                sys.stdout.write("incompatible updmethod chosen for adp_policy_search")
                __dbg = raw_input("execution halted, press ENTER to exit")
                sys.exit()

            #TODO7: this is awkward and cheating bigtime, need to find a better way to determine last time step w/o
            #       violating encapsulation and data obstruction
            if t == len(models)-1:
                slak_var_costs = copy.copy(mC["slak_var_costs"])
            else:
                # slack cost is now also replaced by a theta value
                for skidx in xrange(len(mC["slak_var_names"])):
                    slak_var_costs[skidx] = -theta[ convert_cplex_structs(None, mC["slak_var_names"][skidx]) ]
        else:
            print "in adp_policy_search, invalid policy method given"
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        # bring it all together into the model
        models[t].objective.set_linear( \
            [ (mC["temp_var_insts"][i], vfa[i] + inst_var_costs[i]) for i in xrange(len(mC["temp_var_insts"])) ] +\
            [ (i, (1./mC["numscens"])*0) for i in mC["orig_var_names"] ] + \
            [ (i, 0) for i in mC["temp_var_flows"] ] + \
            [ (mC["slak_var_names"][i], (1./mC["numscens"])*slak_var_costs[i]) for i in xrange(len(slak_var_costs))])
    else:
        print "in adp_policy_search, the objective hasn't been set because the appropriate feature isn't present"
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    # if an arc was chosen for installation previously, it is hereafter always installed: installing
    #    an arc @ t-1 means the decision choice @ time t must be to install it, so that when the
    #    decision choices are recorded @ time t, the lower bounds will include those of t-1
    models[t].variables.set_lower_bounds(zip(mC["temp_var_insts"], [int(i) for i in mC["instchoices"][t-1] ] ) )


    # can only install a certain number of arcs per time step
    # the 1st (obvious) condition in the if statement is to avoid an unnecessary error message from cplex
    if models[t].linear_constraints.get_num() > 0 and "PerTInstConstr" in models[t].linear_constraints.get_names():
        models[t].linear_constraints.delete("PerTInstConstr")
    models[t].linear_constraints.add( \
        lin_expr = [cplex.SparsePair(ind = [i for i in mC["temp_var_insts"] if i not in mC["previnsts"][t-1]], \
                                     val = [1]*(len(mC["temp_var_insts"])-len(mC["previnsts"][t-1])) ) ], \
                    senses = ["L"], rhs = [mC["pertinstbudget"]], names = ["PerTInstConstr"] )
    try:
        if dbout:
            #DEBUG begin
            adplt = time.time()
            #DEBUG end
        models[t].solve()
        if dbout:
            #DEBUG begin
            adpla += time.time() - adplt
            timedbg.write("   ADP loop avg time: "+str(adpla/float(t+1))+" sec\n")
            #DEBUG end

        optimalpolicy = [int(round(i)) for i in models[t].solution.get_values(mC["temp_var_insts"]) ]

        if randpol:
            # previously installed arcs musn't be considered to be randomly uninstalled, so just pass the newly
            #    installed list of arcs
            considerpolicy = []
            # they also mustn't be considered when a random arc is to be installed
            exclarcs = []
            for i in xrange(len(optimalpolicy)):
                if not optimalpolicy[i]:
                    considerpolicy.append(0)
                elif mC["instchoices"][t-1][i]:
                    considerpolicy.append(0)
                    exclarcs.append(i)
                else:
                    considerpolicy.append(1)
            mC["instchoices"][t] = generate_per_arc_random_policy(mC["pertinstbudget"],mC["explorepolicy"],\
                                                                  considerpolicy, banarcs=set(exclarcs))
            # now a new batch of arcs is being installed, include the previously installed arcs into the total set
            for i in xrange(len(mC["instchoices"][t-1])):
                if mC["instchoices"][t-1][i]:
                    mC["instchoices"][t][i] = 1
        else:
            mC["instchoices"][t] = optimalpolicy

##        if mC["explorepolicy"] > 0 and random.uniform(0,1) < mC["explorepolicy"]:
##            mC["instchoices"][t] = generate_random_policy(t,mC)
##        else:
##            mC["instchoices"][t] = optimalpolicy

        tmpflows = models[t].solution.get_values(mC["temp_var_flows"])
        slacks = models[t].solution.get_values(mC["slak_var_names"])
        objval = models[t].solution.get_objective_value()
        unmetdemandcost = (1./mC["numscens"])*sum( slak_var_costs[i]*slacks[i] for i in \
                                                   xrange(len(mC["slak_var_names"])) )
        vfacost = sum( vfa[i]*mC["instchoices"][t][i] for i in xrange(len(vfa)) )
        mC["previnsts"][t].extend(mC["temp_var_insts"][i] for i in xrange(len(mC["instchoices"][t])) if \
                                                                                           mC["instchoices"][t][i])
    except:
        sys.stderr.write("exeADP: Solve didn't work in adp_policy_search\n")
        models[t].conflict.refine(models[t].conflict.all_constraints())
        conflicts = models[t].conflict.get()
        conflicts = [i for i in xrange(len(conflicts)) if conflicts[i] != -1]
        cgs = models[t].conflict.get_groups(conflicts)
        ubcs=[j[0][1] for i,j in cgs if j[0][0] == 2]
        lbcs=[j[0][1] for i,j in cgs if j[0][0] == 1]
        lccs=[j[0][1] for i,j in cgs if j[0][0] == 3]
        constrByVar = find_constr_4_vars(models[t], models[t].variables.get_names())
        conflConstrs = find_constr_4_vars(models[t], models[t].linear_constraints.get_names(lccs), \
                                          vartype="constraint")
        pass

##    # reset costs
##    mC["inst_var_costs"] = []
##    for i,j,d in netobj.tnet.edges_iter(data=True):
##        mC["inst_var_costs"] += [d["weight"]]
    if dbout:
        #DEBUG start
##            sys.stdout.write("installing: ")
##            for i in xrange(len(mC["instchoices"][t])):
##                if mC["instchoices"][t][i] != 0:
##                    sys.stdout.write("["+mC["temp_var_insts"][i]+", "+str(tmpflows[i])+"];  ")
##            sys.stdout.write("\n")
##            firsttime = True
##            printnl = False
##            for i in xrange(len(slacks)):
##                if slacks[i] != 0:
##                    if firsttime:
##                        sys.stdout.write("time "+str(t)+" unmet demand: ")
##                        firsttime = False
##                        printnl = True
##                    sys.stdout.write("["+str(mC["slak_var_names"][i])+": "+str(slacks[i])+"];  ")
##            if printnl:
##                sys.stdout.write("\n")
        pass
        #DEBUG end
    return mC["instchoices"][t],tmpflows,slacks,unmetdemandcost,vfacost






def adp_vfa_update_grad_descent(trueobj, featurelist, tsteps, mC, vfaobj=None, decvars=None):
    # the following are constant terms involved in calculating step sizes
    # a linear scaling factor in the inverse-proportional and harmonic step sizes. a setting of 1 is the "standard"
    k = 40.0
    # a constant value representing the step size for all iterations (i realize this is a redundant variable)
    M = 0.90
    # a target constant step size for the McClain formula step size
    targetalpha = 0.03

    if "sp_" in mC["bftype"]:
        basisname = "_sp_basis_val"
    elif "ones" in mC["bftype"]:
        basisname = "_1s_basis_val"
    elif "node_constr" in mC["bftype"]:
        basisname = "_nc_basis_val"
    else:
        print "adp_vfa_update_grad_descent: invalid basis function type specified."
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    # the decision variables list must correspond to the number of installable arcs exactly
    # this vector of decision variables alters the descent direction to align with the currently determined minimum
    #    it's actually only needed for arc-based methods. the decisions have already been factored into
    #    node_constr's basis values.
    if decvars != None and len(decvars[0]) != mC["netobj"].tnet.number_of_edges():
        sys.stdout.write("Invalid decision variable vector length. Defaulting to all 1s\n")
        decvars = [[1 for i in xrange(mC["netobj"].tnet.number_of_edges())] for j in xrange(tsteps)]
    elif decvars == None:
        decvars = [[1 for i in xrange(mC["netobj"].tnet.number_of_edges())] for j in xrange(tsteps)]

    #ASSUMPTION: there are no fewer than 2 time steps (makes sense; 1 step is trivial)
    for t in xrange(tsteps):
        theta = mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1]
        # these are functions of temp arcs, node_constr's are functions of nodes
        if "ud_vfa_combined" in mC["polmethod"] or "ud_vfa_separate" in mC["polmethod"]:
            bvector = [ad["t"+str(t)+basisname] for af, at, ad in mC["netobj"].tnet.edges_iter(data=True)]
            if "_direct_zeta" in mC["updmethod"] or "_withvalfn_zeta" in mC["updmethod"]:
                bvector = [bvector[i]*decvars[t][i] for i in xrange(len(bvector))]
            elif "_compl_zeta" in mC["updmethod"]:
                bvector = [bvector[i]*(1-decvars[t][i]) for i in xrange(len(bvector))]
            else:
                print "adp_vfa_update_grad_descent: invalid zeta specification in updmethod"
                __dbg = raw_input("execution halted, press ENTER to exit")
                sys.exit()
        elif "node_constr" in mC["polmethod"]:
            # it does however require a re-evaluation of basis values *after* the policy has been found
            #ASSUMPTION: re-evaluation of basis values once a new policy is found is to be done before invoking
            #            this function
            bvector = [d["t"+str(t)+basisname] for n,d in mC["netobj"].tnet.nodes_iter(data=True)]
        else:
            print "adp_vfa_update_grad_descent: invalid policy method"
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        n = mC["thetaupdcnt"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1]

        # there are several ways to compute the step size; we consider 4 here. the 1st type is to make the step size
        #    inversely proportional to both iteration count and time. this is purely speculative as to its efficacy.
        #    the other 3 are conventional: 1) a constant step size, 2) harmonic step size (inversely proportional to
        #    iteration count), and 3) McClain's form (recursively defined from previous iteration values). rather
        #    than making it a runtime option, we just hardcode it here to find the best then leave it; however, in
        #    the future this may need to be tweaked to suit different conditions (e.g. differences in scenario MC
        #    simulations causing wide changes in outage distributions and time dependencies). all parameters involved
        #    in calculating the step size are given at the beginning of the function.

        # 0) step size alpha is proportional to time: let the last time step with nonzero theta values (time T-1)
        #    have an undisturbed coefficient (i.e. step size == 1/n), thus tsteps-1 - t
        if "_tdepss" in mC["updmethod"]:
            alpha = k / ( k + n*(tsteps-1-t) )
        # 1) constant step size
        elif "_constss" in mC["updmethod"]:
            alpha = M
        # 2) harmonic step size
        elif "_harmonicss" in mC["updmethod"]:
            alpha = k / (k+n-1)
        # 3) McClain step size. recursive formula requires either keeping track of previous iteration value or
        #    recalculating every time. this is incredibly wasteful, but storing the step size is a huge pain atm.
        elif "_mcclainss" in mC["updmethod"]:
            alpha = 1.0
            for i in xrange(1,n):
                alpha = alpha / (1+alpha - targetalpha)
        else:
            print "adp_vfa_update_grad_descent: invalid step size specification in updmethod"
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        # it's a little confusing which way theta values should trend using the discounted framework...
        if "_negalpha" in mC["updmethod"]:
            alpha *= -1

        if "_linear" in mC["updmethod"]:
            if vfaobj != None:
                vfa = vfaobj[t]
            else:
                vfa = sum(theta[i]*bvector[i] for i in xrange(len(theta)) )
##            # the max operator is used to try and prevent divergent oscillations in theta values over iterations
##            mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
##                [max(0, theta[i] + (k/(n*(k+tsteps-1-t)) )*(trueobj[t] - vfa)*bvector[i]) for i in xrange(len(bvector))]
##            mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
##                [theta[i] + (k/(n*(k+tsteps-1-t)))*(trueobj[t] - max(0,vfa))*bvector[i] for i in xrange(len(bvector))]
            mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
                [theta[i] + alpha*(trueobj[t] - vfa)*bvector[i] for i in xrange(len(bvector))]

            #DEBUG begin
##            sys.stdout.write("    gradient 1-norm @ time "+str(t)+": "+str(abs(trueobj[t] - vfa))+"\n")
            #DEBUG end
        elif "_1norm_lin" in mC["updmethod"]:
            # a modulator for the norm. setting to 0 gives descent of pure 1-norm, setting to 1 gives descent of 2-norm
            m = 0.0
            if vfaobj != None:
                vfa = vfaobj[t]
            else:
                vfa = sum(theta[i]*bvector[i] for i in xrange(len(theta)) )
            if trueobj[t] - vfa > 0:
                # don't bother multiplying by (1+m), we only care about direction. let the magnitude be from alpha
                mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
                    [theta[i] + alpha*((trueobj[t] - vfa)**m)*bvector[i] for i in xrange(len(bvector))]
            elif trueobj[t] - vfa < 0:
                mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
                    [theta[i] - alpha*((-1*trueobj[t] + vfa)**m)*bvector[i] for i in xrange(len(bvector))]
            else:
                # gradient when m = 0 is actually undefined, but the 0 vector is a subgradient
                mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = 0
        elif "_inverse" in mC["updmethod"]:
            if vfaobj != None:
                vfa = vfaobj[t]
            else:
                vfa = sum(bvector[i]/theta[i] for i in xrange(len(theta)) )
            try:
                mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
                    [theta[i] - alpha*(trueobj[t] - vfa)*(bvector[i]/theta[i]**2) for i in xrange(len(bvector))]
            except ZeroDivisionError:
                sys.stderr.write("Division by zero in vfa update method inverse, blanking theta vector\n")
                mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = [0 for i in xrange(len(bvector))]
        elif "_root_inverse" in mC["updmethod"]:
            # the inverse power term, meant to be between 1 and positive infinity, non-inclusive.
            m = 2.0
            # the following allows theta terms to be negative by taking the negative term out of the root calculation
            neghack = True
            if vfaobj != None:
                vfa = vfaobj[t]
            else:
                if neghack:
                    vfa =  sum(bvector[i]/(theta[i]**(1/m)) for i in xrange(len(theta)) if theta[i] > 0 )
                    vfa += sum(bvector[i]/(-(-theta[i])**(1/m)) for i in xrange(len(theta)) if theta[i] < 0 )
                else:
                    vfa = sum(bvector[i]/(theta[i]**(1/m)) for i in xrange(len(theta)) )
            try:
                if neghack:
                    for i in xrange(len(bvector)):
                        if theta[i] > 0:
                            mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1][i] = \
                                theta[i] - alpha*(trueobj[t] - vfa)*(bvector[i]/(m*theta[i]**(1./m+1)))
                        elif theta < 0:
                            mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1][i] = \
                                theta[i] - alpha*(trueobj[t] - vfa)*(bvector[i]/(-m*(-theta[i])**(1./m+1)))
                        else:
                            mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1][i] = 0
                else:
                    mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
                        [max(0.000001, theta[i] - alpha*(trueobj[t] - vfa)*(bvector[i]/(m*theta[i]**(1./m+1))) ) \
                            for i in xrange(len(bvector))]
            except ZeroDivisionError:
                sys.stderr.write("Division by zero in vfa update method root inverse, blanking theta vector\n")
                mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = [0 for i in xrange(len(bvector))]
        elif "_neg_exp" in mC["updmethod"]:
            if vfaobj != None:
                vfa = vfaobj[t]
            else:
                vfa = sum(bvector[i]*np.exp(-theta[i]) for i in xrange(len(theta)) )
            mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
                [theta[i] - alpha*(trueobj[t] - vfa)*(bvector[i]*np.exp(-1*theta[i])) for i in xrange(len(bvector))]

        else:
            sys.stderr.write("adp_vfa_update_grad_descent: invalid update method chosen.\n")
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        mC["thetaupdcnt"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] += 1

    return







def adp_vfa_update_least_squares(trueudcost, featurelist, tsteps, mC, thetaupdvfavals, decvars=None):
    # the least squares update method is only meant to be used for a linearly represented VFA
    if "_linear" not in mC["updmethod"]:
        print "adp_vfa_update_least_squares: invalid (nonlinear) update method specified."
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    # the following are constant terms involved in calculating step sizes
    # a linear scaling factor in the inverse-proportional and harmonic step sizes. a setting of 1 is the "standard"
    k = 40.0
    # a constant value representing the step size for all iterations (i realize this is a redundant variable)
    M = 0.90
    # a target constant step size for the McClain formula step size
    targetalpha = 0.03

    if "sp_" in mC["bftype"]:
        basisname = "_sp_basis_val"
    elif "ones" in mC["bftype"]:
        basisname = "_1s_basis_val"
    elif "node_constr" in mC["bftype"]:
        basisname = "_nc_basis_val"
    else:
        print "adp_vfa_update_least_squares: invalid basis function specification"
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    # the decision variables list must correspond to the number of installable arcs exactly
    # this vector of decision variables alters the descent direction to align with the currently determined minimum
    #    it's actually only needed for arc-based methods. the decisions have already been factored into
    #    node_constr's basis values.
    if decvars != None and len(decvars[0]) != mC["netobj"].tnet.number_of_edges():
        sys.stdout.write("Invalid decision variable vector length. Defaulting to all 1s\n")
        decvars = [[1 for i in xrange(mC["netobj"].tnet.number_of_edges())] for j in xrange(tsteps)]
    elif decvars == None:
        decvars = [[1 for i in xrange(mC["netobj"].tnet.number_of_edges())] for j in xrange(tsteps)]

    #ASSUMPTION: there are no fewer than 2 time steps (makes sense; 1 step is trivial)
    for t in xrange(tsteps):
        theta = np.array(mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1])
        # these are functions of temp arcs, node_constr's are functions of nodes
        if "ud_vfa_combined" in mC["polmethod"] or "ud_vfa_separate" in mC["polmethod"]:
            bvector = [ad["t"+str(t)+basisname] for af, at, ad in mC["netobj"].tnet.edges_iter(data=True)]
            if "_direct_zeta" in mC["updmethod"] or "_withvalfn_zeta" in mC["updmethod"]:
                bvector = [bvector[i]*decvars[t][i] for i in xrange(len(bvector))]
            elif "_compl_zeta" in mC["updmethod"]:
                bvector = [bvector[i]*(1-decvars[t][i]) for i in xrange(len(bvector))]
            else:
                print "adp_vfa_update_least_squares: invalid zeta specification in updmethod"
                __dbg = raw_input("execution halted, press ENTER to exit")
                sys.exit()
        elif "node_constr" in mC["polmethod"]:
            # it does however require a re-evaluation of basis values *after* the policy has been found
            #ASSUMPTION: re-evaluation of basis values once a new policy is found is to be done before invoking
            #            this function
            bvector = [d["t"+str(t)+basisname] for n,d in mC["netobj"].tnet.nodes_iter(data=True)]
        else:
            print "adp_vfa_update_least_squares: invalid policy method"
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        n = mC["thetaupdcnt"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1]

        # 0) step size alpha is proportional to time: let the last time step with nonzero theta values (time T-1)
        #    have an undisturbed coefficient (i.e. step size == 1/n), thus tsteps-1 - t
        if "_tdepss" in mC["updmethod"]:
            alpha1 = k / ( k + n*(tsteps-1-t) )
            alpha0 = k / ( k + (n-1)*(tsteps-1-t) )
        # 1) constant step size
        elif "_constss" in mC["updmethod"]:
            alpha1 = M
            alpha0 = M
        # 2) harmonic step size
        elif "_harmonicss" in mC["updmethod"]:
            alpha1 = k / (k+n-1)
            alpha0 = k / (k+n-2)
        # 3) McClain step size. recursive formula requires either keeping track of previous iteration value or
        #    recalculating every time. this is incredibly wasteful, but storing the step size is a huge pain atm.
        elif "_mcclainss" in mC["updmethod"]:
            alpha1 = 1.0
            alpha0 = 1.0
            for i in xrange(1,n):
                alpha0 = alpha1
                alpha1 /= (1+alpha1 - targetalpha)
        else:
            print "adp_vfa_update_least_squares: invalid step size specification in updmethod"
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        # this recursive formula was taken directly from Powell's book, pg 350-353
        eps = 1.0
        lamda = alpha0*( (1-alpha1) / alpha1)
        # this is a recursive variable, but because we would have to store multiple depending on the number of
        #    separable features and update them independently like thetaupdcnt, it's too much of a hassle
        B = eps*np.eye(len(bvector))
        bvector = np.array(bvector)
        gamma = lamda + np.dot( np.dot(np.transpose(bvector), B), bvector)
        for i in xrange(n-1):
            B = (1/lamda)*(B - (1/gamma)*( np.dot( np.dot( np.dot(B,bvector),np.transpose(bvector)),B) ) )
            gamma = lamda + np.dot( np.dot(np.transpose(bvector), B), bvector)

        H = (1/gamma)*B
        epshat = thetaupdvfavals[t] - trueudcost[t]
        mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = list(theta - np.dot(H,bvector)*epshat)

        mC["thetaupdcnt"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] += 1
    return








def adp_vfa_update_min_norm(trueudcost, featurelist, tsteps, mC, thetaupdvfavals, decvars=None):
    # the following are constant terms involved in calculating step sizes
    # a linear scaling factor in the inverse-proportional and harmonic step sizes. a setting of 1 is the "standard"
    k = 40.0
    # a constant value representing the step size for all iterations (i realize this is a redundant variable)
    M = 0.90
    # a target constant step size for the McClain formula step size
    targetalpha = 0.03

    if "sp_" in mC["bftype"]:
        basisname = "_sp_basis_val"
    elif "ones" in mC["bftype"]:
        basisname = "_1s_basis_val"
    elif "node_constr" in mC["bftype"]:
        basisname = "_nc_basis_val"
    else:
        print "adp_vfa_update_grad_descent: invalid basis function type specified."
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    # the decision variables list must correspond to the number of installable arcs exactly
    if decvars != None and len(decvars[0]) != mC["netobj"].tnet.number_of_edges():
        sys.stdout.write("Invalid decision variable vector length. Defaulting to all 1s\n")
        decvars = [[1 for i in xrange(mC["netobj"].tnet.number_of_edges())] for j in xrange(tsteps)]
    elif decvars == None:
        decvars = [[1 for i in xrange(mC["netobj"].tnet.number_of_edges())] for j in xrange(tsteps)]

    # the min norm minimizes over theta values, so all new variable names must be generated
    #ASSUMPTION: the order of theta variables in the list directly corresponds to the order in the temp_var_insts
    #            list -- no checking is done to see if they correspond to the same temp arc
    thetavars = ["Theta("+str(i)+","+str(j)+")" for i,j in mC["netobj"].tnet.edges_iter()]
    mnmodels = [cplex.Cplex() for t in xrange(tsteps)]

    for t in xrange(tsteps):
        n = mC["thetaupdcnt"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1]

        # there are several ways to compute the step size; we consider 4 here. the 1st type is to make the step size
        #    inversely proportional to both iteration count and time. this is purely speculative as to its efficacy.
        #    the other 3 are conventional: 1) a constant step size, 2) harmonic step size (inversely proportional to
        #    iteration count), and 3) McClain's form (recursively defined from previous iteration values). rather
        #    than making it a runtime option, we just hardcode it here to find the best then leave it; however, in
        #    the future this may need to be tweaked to suit different conditions (e.g. differences in scenario MC
        #    simulations causing wide changes in outage distributions and time dependencies). all parameters involved
        #    in calculating the step size are given at the beginning of the function.

        # 0) step size alpha is proportional to time: let the last time step with nonzero theta values (time T-1)
        #    have an undisturbed coefficient (i.e. step size == 1/n), thus tsteps-1 - t
        if "_tdepss" in mC["updmethod"]:
            alpha = k / ( k + n*(tsteps-1-t) )
        # 1) constant step size
        elif "_constss" in mC["updmethod"]:
            alpha = M
        # 2) harmonic step size
        elif "_harmonicss" in mC["updmethod"]:
            alpha = k / (k+n-1)
        # 3) McClain step size. recursive formula requires either keeping track of previous iteration value or
        #    recalculating every time. this is incredibly wasteful, but storing the step size is a huge pain atm.
        elif "_mcclainss" in mC["updmethod"]:
            alpha = 1.0
            for i in xrange(1,n):
                alpha = alpha / (1+alpha - targetalpha)
        else:
            print "adp_vfa_update_min_norm: invalid step size specification in updmethod"
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        mnmodels[t].set_results_stream(None)
        mnmodels[t].set_log_stream(None)
        mnmodels[t].set_warning_stream(None)

        theta = mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1]
        trueudthetaidxs = [i for i in xrange(len(decvars[t])) if decvars[t][i]]
        thetaweights = [1.0 for i in theta]
        for i in trueudthetaidxs:
            thetaweights[i] = (mC["netobj"].tnet.number_of_edges() - len(trueudthetaidxs)) / \
                               float(len(trueudthetaidxs))
        # these are functions of temp arcs, node_constr's are functions of nodes
        if "ud_vfa_combined" in mC["polmethod"] or "ud_vfa_separate" in mC["polmethod"]:
            bvector = [ad["t"+str(t)+basisname] for af, at, ad in mC["netobj"].tnet.edges_iter(data=True)]
        else:
            print "adp_vfa_update_min_norm: invalid policy method (node_constr is currently unsupported)."
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()

        mnmodels[t].variables.add(obj   = [-2*thetaweights[i]*theta[i] for i in xrange(len(theta))], \
                                  lb    = [-cplex.infinity]*len(theta), types = "C"*len(theta), names = thetavars)
        mnmodels[t].objective.set_quadratic_coefficients([(thetavars[i],thetavars[i],thetaweights[i]) \
                                                          for i in xrange(len(thetavars)) ])
        mnmodels[t].linear_constraints.add( \
            lin_expr = [cplex.SparsePair(ind = [thetavars[i] for i in trueudthetaidxs], \
                                         val = [bvector[i]   for i in trueudthetaidxs] ) ], \
                        senses = ["E"],  rhs = [trueudcost[t]-thetaupdvfavals[t]], names = ["DiscountMatching"] )
        mnmodels[t].linear_constraints.add( \
            lin_expr = [cplex.SparsePair(ind = [thetavars[i], thetavars[j]], \
                                         val = [bvector[i], -bvector[j]] ) for j in trueudthetaidxs \
                            for i in xrange(len(theta)) if i not in trueudthetaidxs], \
                        senses = ["G"]*len(trueudthetaidxs)*(len(theta)-len(trueudthetaidxs)), \
                        rhs    = [0]*len(trueudthetaidxs)*(len(theta)-len(trueudthetaidxs)), \
                        names  = ["Theta"+str(i)+"gteTheta"+str(j) for j in trueudthetaidxs \
                                     for i in xrange(len(theta)) if i not in trueudthetaidxs] )
        mnmodels[t].set_problem_type(cplex.Cplex.problem_type.QP)
        try:
            mnmodels[t].solve()
            mnvals = mnmodels[t].solution.get_values(thetavars)
            mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = [theta[i] + alpha*(mnvals[i]-theta[i]) \
                                                                           for i in xrange(len(theta))]
        except:
            sys.stderr.write("exeADP: Solve didn't work in adp_vfa_update_min_norm\n")
            mnmodels[t].conflict.refine(mnmodels[t].conflict.all_constraints())
            conflicts = mnmodels[t].conflict.get()
            conflicts = [i for i in xrange(len(conflicts)) if conflicts[i] != -1]
            cgs = mnmodels[t].conflict.get_groups(conflicts)
            ubcs=[j[0][1] for i,j in cgs if j[0][0] == 2]
            lbcs=[j[0][1] for i,j in cgs if j[0][0] == 1]
            lccs=[j[0][1] for i,j in cgs if j[0][0] == 3]
            constrByVar = find_constr_4_vars(mnmodels[t], mnmodels[t].variables.get_names())
            conflConstrs = find_constr_4_vars(mnmodels[t], mnmodels[t].linear_constraints.get_names(lccs), \
                                              vartype="constraint")

        mC["thetaupdcnt"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] += 1
    return








def adp_vfa_update_lp_toggle(tsteps, models, mC, baseobjvals, featurelist):
    #NOTE: what we're establishing here is v[t-1][f_n][i] = a*v[t-1][f_n][i] + b*(V[t]-V_i[t]); that is, the
    #      change to the value funct approx for time step t is an accumulating average of the existing value
    #      and a "gradient step" given by the difference in objectives between the adp solution and the
    #      perturbed model (same as adp, but toggle the install choice for the arc in question). f_n
    #      illustrates that only the vfa parameters corresponding to the features present in this iteration
    #      are being updated. a & b are some constants, typically related to stepsize or some convex combo.

    # the first time step is excluded because the last time step doesn't use val funct approx. (see below)
    for t in xrange(1,tsteps):
        # to start, consider just slack (to reduce uncertainty in analyzing results)
        models[t].objective.set_linear( \
            [ (mC["temp_var_insts"][j], 0) for j in xrange(len(mC["temp_var_insts"])) ]+ \
            [ (mC["orig_var_names"][j], (1./mC["numscens"])*0) for j in xrange(len(mC["orig_var_names"])) ]+ \
            [ (mC["temp_var_flows"][j], 0) for j in xrange(len(mC["temp_var_flows"])) ]+ \
            [ (mC["slak_var_names"][j], (1./mC["numscens"])*mC["slak_var_costs"][j]) \
                                           for j in xrange(len(mC["slak_var_names"])) ] )

        #DEBUG begin
        pbla = 0
        #DEBUG end

        for i in xrange(len(mC["temp_var_insts"])):
            #ASSUMPTION: if we're considering an arc that wasn't installed at a time where the objective cannot be
            #            improved upon (i.e. no unmet demand), then we shouldn't bother with updating -- doing so
            #            only serves to lower the coeff value toward 0 (due to the moving average formula).
            #NOTE: this assumption becomes invalid if the update takes into account the cost of installing the arc.
            if baseobjvals[t] == 0 and not mC["instchoices"][t-1][i]:
                continue

            models[t].variables.set_lower_bounds(mC["temp_var_insts"][i],1-int(mC["instchoices"][t-1][i]))
            # if an arc was chosen to be installed, make sure it stays uninstalled by lowering the ub
            if mC["instchoices"][t-1][i]:
                models[t].variables.set_upper_bounds(mC["temp_var_insts"][i],0)
            # if by perturbing the problem we install an arc, the install budget for this t is potentially violated
            #    so some means must reflect that in the perturbed model
            else:
                # one way: "borrow" an install slot from this time period by reducing the install budget by 1.
                #          warning: if the install budget is 1, this will consistently provide infeasible problems
                # another way is to just let it slide, which is what i'll do for now
##                models[t].linear_constraints.set_rhs("PerTInstConstr", mC["pertinstbudget"]-1)
                pass

            # new way: we wish to see how crucial an arc is and so we don't allow any other installs, making
            #    this an LP with much faster runtime. a downside of this is that you're not comparing apples to
            #    apples: fixing all arcs means you aren't installing anything in this time step, whereas the
            #    control model DOES get to install arcs in the current time step. fixing the installs on the
            #    vfa side is only ensuring that arcs will have been installed in the previous time step only.
            #NOTE: the easiest way i can ensure good performance is by changing the install variable type
            #      to continuous and set ub=lb, then switch them back afterwards. update: this is really slow
            models[t].variables.set_types([ (mC["temp_var_insts"][j], \
                                             models[t].variables.type.continuous) \
                                           for j in xrange(len(mC["temp_var_insts"])) ])
            models[t].variables.set_lower_bounds([ (mC["temp_var_insts"][j], \
                                                    mC["instchoices"][t-1][j]) \
                                                   for j in xrange(len(mC["temp_var_insts"])) \
                                                  if j != i] )
            models[t].variables.set_upper_bounds([ (mC["temp_var_insts"][j], \
                                                    mC["instchoices"][t-1][j]) \
                                                   for j in xrange(len(mC["temp_var_insts"])) \
                                                  if j != i] )
            models[t].set_problem_type(models[t].problem_type.LP)
            try:
                #DEBUG begin
                pblt = time.time()
                #DEBUG end
                models[t].solve()
                #DEBUG begin
                pbla += time.time() - pblt
                #DEBUG end

                # there are 2 cases: if we installed the arc (and have perturbed it to be uninstalled) and
                #    the perturbed solution is better than the pure solution, the value approx for that arc
                #    should increase; else, if the perturbed problem forces a previously uninstalled arc to
                #    be installed and the solution is better, it means we were better off installing it and
                #    its value approx should decrease.
                if mC["instchoices"][t-1][i]:
                    updobjval = baseobjvals[t] - models[t].solution.get_objective_value()
                else:
                    updobjval = models[t].solution.get_objective_value() - baseobjvals[t]

                #ASSUMPTION: the only features being used are ctr_mass_quad_rel_frame. using any other
                #            will break the following code.
                # kind of a kluge.. the ctr_mass feature returns 1 thru 4 indicating cartesian quadrant
                #    so subtract 1 to provide the index to the list. to repeat: the only vfa param
                #    being updated is the one exhibiting the specific feature at the time of decision making (t-1).
                if "ctr_mass_quad_rel_frame" in featurelist[t-1]:
                    #NOTE: negative values of v are degenerate: they will always be installed even if
                    #      not used; this should change upon incorporation of install & flow costs into
                    #      update step. I also think 0 is degenerate - nothing should be "free", but this needs to
                    #      considered further analytically i think....
                    n = mC["thetaupdcnt"][t-1][featurelist[t-1]["ctr_mass_quad_rel_frame"]-1]
                    mC["theta"][t-1][featurelist[t-1]["ctr_mass_quad_rel_frame"]-1][i] = \
                        max( (1 - (1./(n)))*mC["theta"][t-1][featurelist[t-1]["ctr_mass_quad_rel_frame"]-1][i] + \
                                  (1./(n))*updobjval, 0.0)
                    mC["thetaupdcnt"][t-1][featurelist[t-1]["ctr_mass_quad_rel_frame"]-1] += 1
            except:
                sys.stderr.write("exeADP: Solve in update step didn't work for arc "+mC["temp_var_insts"][i]+ \
                                 " in time "+str(t)+", maintaining original value\n")
            # restore the bound to what it was
            models[t].variables.set_lower_bounds(mC["temp_var_insts"][i], int(mC["instchoices"][t-1][i]))
            if mC["instchoices"][t-1][i]:
                models[t].variables.set_upper_bounds(mC["temp_var_insts"][i],1)
            else:
##                models[t].linear_constraints.set_rhs("PerTInstConstr", mC["pertinstbudget"])

                pass

            models[t].set_problem_type(models[t].problem_type.MILP)
            models[t].variables.set_lower_bounds([ (mC["temp_var_insts"][j], \
                                                    mC["instchoices"][t-1][j]) \
                                                   for j in xrange(len(mC["temp_var_insts"])) \
                                                  if j != i] )
            models[t].variables.set_upper_bounds([ (mC["temp_var_insts"][j], 1) \
                                                   for j in xrange(len(mC["temp_var_insts"])) \
                                                  if j != i] )
            models[t].variables.set_types([ (mC["temp_var_insts"][j], \
                                             models[t].variables.type.binary) \
                                           for j in xrange(len(mC["temp_var_insts"])) ])
    return pbla





def generate_random_policy(t,mC):
    # simple: pick up to the max # of installs this time step (weighted to install the max), and insert them
    newpolicy = mC["instchoices"][t-1]
    samplesize = int(round(random.triangular(low=0, high=mC["pertinstbudget"], mode=mC["pertinstbudget"])))
    sampleset = [iidx for iidx in xrange(len(newpolicy)) if not int(newpolicy[iidx]) ]
    polpicks = random.sample(sampleset, samplesize)
    for i in polpicks:
        newpolicy[i] = 1
    return newpolicy






def generate_per_arc_random_policy(numarcs, eps, optpolicy, banarcs=set([])):
    # each installation decision here is determined at random individually, including the possibility of no install
    # an arc has a probability of being selected (1-eps)% if it's installed in optpolicy or
    #    eps/len(mC["temp_var_insts]) if not -- this includes a null arc indicating that an arc is not installed
    #    repeat this selection procedure mC["pertinstbudget"] times
    narcs = len(optpolicy)
    newpolicy = [0 for i in xrange(narcs) ]
    optidxs = [i for i in xrange(len(optpolicy)) if optpolicy[i] ]
    # the number of arcs that were selected may be less than the total install budget
    noptinst = len(optidxs)
    repopt = []
    rndchoices = 0
    for a in xrange(numarcs):
        pick = random.uniform(0,1)
        try:
            repopt.append(optidxs.pop())
        except:
            # there was, in fact, no arc to be installed and so the optimal choice is to leave it this way
            if pick > eps:
                # don't break, because we may randomly choose to install on subsequent iterations
                continue
        if pick > eps:
            newpolicy[repopt[-1]] = 1
            # epsilon is divided into equal portions of either selecting no arc or selecting an arc from the set of
            #    arcs not including the original set that was chosen for installation nor the banned arcs
        elif pick <= float(eps)/(narcs-noptinst-len(banarcs)):
            # null arc selected => leave uninstalled
            continue
        else:
            # wait until all "rolls of the dice" have been made to prevent random selection of optimal arcs
            rndchoices += 1
    for r in xrange(rndchoices):
        # for suboptimal arcs, they're all equally likely to be chosen, so pick one uniformly
        while True:
            suboptarc = random.sample(set(xrange(narcs)) - set(repopt) - banarcs, 1)
            if newpolicy[suboptarc[0]]:
                continue
            else:
                newpolicy[suboptarc[0]] = 1
                break
    return newpolicy






def adp_warm_start(mC, tsteps, features):
    print "warm start begin"
    # create damage distributions for each time step. all parameters are in the function, tweak them there.
    mC["mcsim"].generate_distributions(mC["netobj"].dinet, tsteps)
    mC["mcsim"].eval_adp_basis_functions(mC["netobj"], tsteps, mC)
    featurelist = mC["mcsim"].calc_feature_ctrmassquad(mC["netobj"], tsteps, features)

    # we're not performing SAA here, so make sure we only have 1 sample; however, even though there is only one
    #    sample, the model is now technically a SAA (where the sample size is 1).
    mC["mcsim"].clear_scenarios(mC["netobj"].dinet)
    mC["numscens"] = mC["mcsim"].sample_scenarios(mC["netobj"].dinet)
    update_models_orig_scens(models, mC, mrange=xrange(len(models)), reset=True)

    #DEBUG start
##        sys.stdout.write("time ")
    #DEBUG end

    models, orig_var_caps = create_models_temp_framework(tsteps, mC)
    for t in xrange(tsteps):
        try:
            models[t].solve()
            #NOTE: since we seek to minimize cost, the values get flipped from what was chosen (i.e. we prefer
            #      they be chosen again)
            tmpflows = models[t].solution.get_values(mC["temp_var_flows"])
            maxtmpflows = max(tmpflows)+1
            # the v's represent value derived from the future time step, so the values found here are
            #    actually for the previous time step -- the last time step runs a flow problem as normal
            # warm start only the sampled feature. set all the others to maxtmpflows
            for f in xrange(mC["numfeats"]):
                # kind of a kluge.. the ctr_mass feature returns 1 thru 4 indicating cartesian quadrant
                #    so add 1 to the loop index to check for the key. to repeat: the only vfa param
                #    being updated is the one exhibiting the specific feature.
                if "ctr_mass_quad_rel_frame" in featurelist[t] and \
                           f+1 == featurelist[t]["ctr_mass_quad_rel_frame"]:
                    # vfa params will be derived from cost of unmet demand (slack flow cost), so give a relative
                    #    measure of value using the temp flows with slack costs
                    mC["theta"][t-1][f] = [mC["slackcost"]*(maxtmpflows-w) for w in tmpflows]
                else:
                    mC["theta"][t-1][f] = [mC["slackcost"]*maxtmpflows for w in tmpflows]
        except:
            sys.stderr.write("exeADP: cplex solve didn't work in adp_warm_start")

    #DEBUG start
##        sys.stdout.write("\n")
    mla = 0
    #DEBUG end

    # place holders for last time step since they're not used (it's a regular flow model)
    mC["theta"][-1] = [[0 for i in xrange(len(mC["temp_var_insts"]))] for f in xrange(mC["numfeats"])]

    trash_model(models)
    return









def execute_alt_method(tsteps, mC, bigN, features, recout="recthetas"):

    #TODO9: test this
    def control_exhaustive_search(mC, numforecasts, comboset):
##        totcombs = 0
##        for ccs in xrange(mC["pertinstbudget"]+1):
##            totcombs += comb(len(mC["temp_var_insts"]), ccs)
##        sys.stdout.write("   begin search of action space ("+str(int(totcombs))+" total)... ")
##        bests = []
##        #NOTE: this may not be needed. if phi and theta remain positive for all time and iterations, then it will
##        #      always behoove us to install the maximum number of arcs every time period
##        for combosize in xrange(mC["pertinstbudget"]+1):
##            for acombo in itertools.combinations(comboset, combosize):
##                vfa = [[0 for i in xrange(mC["numfeats"])] for nf in xrange(numforecasts) ]
##                for nf in xrange(numforecasts):
##                    edgecnt = 0
##                    nx.set_edge_attributes(mC["netobj"].tnet, "t"+str(t)+"_sp_basis_val", basisattr[nf])
##                    for af, at, ad in mC["netobj"].tnet.edges_iter(data=True):
##                        # the arcs in the combo list are the ones that we set phi = 0
##                        if "TempArc("+str(af)+","+str(at)+")" in acombo:
##                            edgecnt += 1
##                            continue
##                        else:
##                            for feat in xrange(mC["numfeats"]):
##                                #NOTE: we're calculating the recursive form of the VFA which requires info of time
##                                #      t+1; however, getting access to probabilities (phi's) of the next time step
##                                #      violates nonanticipativity (i think -- see notes 2012-10-01)
##                                vfa[nf][feat] += mC["theta"][t+1][feat][edgecnt]*ad["t"+str(t)+"_sp_basis_val"]
##                        edgecnt += 1
##                avgvfa = sum(vfa[i][j] for i in xrange(len(vfa)) for j in xrange(len(vfa[i])) ) /  \
##                          float(numforecasts*numfeats)
##                if avgvfa < vaval:
##                    bests = acombo
##                    vaval = avgvfa
##                for nf in xrange(numforecasts):
##                    for feat in xrange(mC["numfeats"]):
##                        if vfa[nf][feat] < timeTv.values()[0]:
##                            # the key is the fcast,feat pair for identification purposes
##                            timeTv = dict([ ( str(nf)+","+str(feat), vfa[nf][feat]),  ])
##        sys.stdout.write("done\n")
        return bests

    import testADP
    import itertools
##    from scipy.misc import comb

    rectheta = []
    numforecasts = 20
    samplesize = 50
    sys.stdout.write("Iteration ")

    # need to warm start thetas for the last time step with this method
    mC["theta"][tsteps-1] = copy.deepcopy(mC["theta"][0])

    for n in xrange(bigN):
        sys.stdout.write(str(n)+": \n")
        sys.stdout.write("   Backing up network... ")
        netobj = copy.deepcopy(mC["netobj"])
        sys.stdout.write("done.\n")

        sys.stdout.write("   Generating forecast timelines... ")
        # generate the candidate forecast & generate the training forecasts
        tfcastdmgattr, tfcastoutattr, cfcastdmgattr, cfcastoutattr = \
            testADP.generate_sample_data(mC, numforecasts, tsteps, samplesize, fixfc0=False)
        sys.stdout.write("done.\n")

        # for candidate, just pick the first forecast
        sys.stdout.write("   Calculating VFA terms... ")
        rankedarcs = alt_vfa_calculation(tsteps, mC, netobj, features, numforecasts, tfcastdmgattr)
        sys.stdout.write("done.\n")

        sys.stdout.write("   Executing forward step... ")
        policies = alt_control_step(tsteps, mC, netobj, samplesize, rankedarcs)
        sys.stdout.write("done.\n")

        sys.stdout.write("   Executing backward step... ")
        cfcastdmgattr = cfcastdmgattr[0]
        cfcastoutattr = cfcastoutattr[0]
        # run max flow problem on candidate forecast's state at time T =: vr (value of return)
        vrmodel = create_models_temp_framework(tsteps, mC)
        nx.set_edge_attributes(netobj.dinet, "t"+str(tsteps-1)+"_capacity", cfcastoutattr[tsteps-1])
        # use samplesize number of outages, obtained from generate_sample_data
        update_models_orig_scens(vrmodel, mC, netobj, numscens=len(cfcastoutattr[tsteps-1].values()[0]), \
                                        mrange=[tsteps-1], reset=True, updateConstr=True)

        vrmodel[tsteps-1].variables.set_types([ (mC["temp_var_insts"][j], \
                                                vrmodel[tsteps-1].variables.type.continuous) \
                                       for j in xrange(len(mC["temp_var_insts"])) ])
        # set lb of arcs to policies[tsteps-1]
        vrmodel[tsteps-1].variables.set_lower_bounds( \
            [ (mC["temp_var_insts"][i], policies[tsteps-1][i]) for i in xrange(len(mC["temp_var_insts"])) ])
        vrmodel[tsteps-1].variables.set_upper_bounds( \
            [ (mC["temp_var_insts"][i], policies[tsteps-1][i]) for i in xrange(len(mC["temp_var_insts"])) ])
        vrmodel[tsteps-1].set_problem_type(vrmodel[tsteps-1].problem_type.LP)
        try:
            vrmodel[tsteps-1].solve()
            vr = vrmodel[tsteps-1].solution.get_objective_value()
        except:
            print "alt method exeADP: test LP SAA model infeasible"

        # calculate candidate forecast timeline's sp_basis_vals: set_edge_attributes in netobj to candidate's var
        nx.set_edge_attributes(netobj.dinet, "t"+str(tsteps-1)+"_dmg_pct", cfcastdmgattr[tsteps-1])
        mC["mcsim"].eval_adp_basis_functions(netobj, tsteps, mC)
        featurelist = mC["mcsim"].calc_feature_ctrmassquad(netobj, tsteps, features)

##        alt_vfa_update_recursive(tsteps, mC, featurelist, netobj, vr)
        # slight difference for the alternate method: the thetas at time T are used, and it's implicitly not used
        #    in the adp update steps, so add 1 to tsteps to force it to be used
        adp_vfa_update_linear(vr, featurelist, tsteps+1, mC)
        sys.stdout.write("done.\n")
        rectheta.append(copy.deepcopy(mC["theta"]))

    qpick(recout+".pickle","w",rectheta)



    #TODO: try and salvage, but for now gonna scrap the whole thing and start over
    while code is not used:
##    sys.stdout.write("Begin alternate method\n")
##    smin = [[] for t in xrange(tsteps)]
##    comboset = copy.copy(mC["temp_var_insts"])
##    timeTv = dict([ ("-1,-1", sys.maxint), ])
##    for t in xrange(tsteps-1):
##
##        sys.stdout.write("Time step "+str(t)+": \n")
##
##        vaval = sys.maxint
##        basisattr = [dict() for nf in xrange(numforecasts)]
##
##        sys.stdout.write("   generating forecast data... ")
##
##        for nf in xrange(numforecasts):
##            mC["mcsim"].generate_distributions(mC["netobj"].dinet, tsteps)
##            mC["mcsim"].eval_adp_basis_functions(mC["netobj"], tsteps, mC)
##            basisattr[nf] = nx.get_edge_attributes(mC["netobj"].tnet, "t"+str(t)+"_sp_basis_val")
##        sys.stdout.write("done\n")
##
##        # one way: check every possible action (set of arcs to install) & return the best per t
##        bests = control_exhaustive_search(mC, numforecasts, comboset)
##
##        sys.stdout.write("   begin cleanup... ")
##        # convert arc text to numeric listing
##        numbs = [0 for i in xrange(len(mC["temp_var_insts"]))]
##        for an in bests:
##            numbs[mC["temp_var_insts"].index(an)] = 1
##        actuals = generate_per_arc_random_policy(mC["pertinstbudget"], mC["explorepolicy"], numbs)
##
##        #TODO: we have numeric arrays and we need to reconvert to text
##
##        comboset = list(set(comboset) - set(actuals))
##        # t+1 reflects the fact that the arcs aren't installed until after
##        smin[t+1] = smin[t] + actuals
##        sys.stdout.write(" done.\n")
##        # now to run max flows
        pass

    return








def alt_vfa_update_recursive(tsteps, mC, featurelist, netobj, vr):
    damper = 1.0

    # first calculate for time step T
    # need to copy & alter adp_vfa_update_linear as it doesn't evaluate for the last time step
    theta = mC["theta"][-1][featurelist[-1]["ctr_mass_quad_rel_frame"]-1]
    bvector = [ad["t"+str(tsteps-1)+"_sp_basis_val"] for af, at, ad in netobj.tnet.edges_iter(data=True)]
    vfa = sum(theta[i]*bvector[i] for i in xrange(len(theta)) )
    n = mC["thetaupdcnt"][-1][featurelist[-1]["ctr_mass_quad_rel_frame"]-1]

    mC["theta"][-1][featurelist[-1]["ctr_mass_quad_rel_frame"]-1] = \
        [theta[i] + (1./n)*(vr - vfa)*bvector[i] for i in xrange(len(bvector))]
    mC["thetaupdcnt"][-1][featurelist[-1]["ctr_mass_quad_rel_frame"]-1] += 1

    # for previous time steps, recalculate VFA with new thetas - just take a chunk of code out of adp_vfa_upd
    for t in xrange(tsteps-2,-1,-1):
        thetatp1 = mC["theta"][t+1][featurelist[t+1]["ctr_mass_quad_rel_frame"]-1]
        bvectortp1 = [ad["t"+str(t+1)+"_sp_basis_val"] for af, at, ad in netobj.tnet.edges_iter(data=True)]
        vrhat = sum(thetatp1[i]*bvectortp1[i] for i in xrange(len(thetatp1)) )

        theta = mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1]
        bvector = [ad["t"+str(t)+"_sp_basis_val"] for af, at, ad in netobj.tnet.edges_iter(data=True)]
        vfa = sum(theta[i]*bvector[i] for i in xrange(len(theta)) )
        n = mC["thetaupdcnt"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1]
        alpha = damper/(n*(damper+tsteps-2-t))

        mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] = \
            [theta[i] + alpha*(vrhat - vfa)*bvector[i] for i in xrange(len(bvector))]
        mC["thetaupdcnt"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1] += 1

    return







def alt_vfa_calculation(tsteps, mC, netobj, features, numforecasts, tfcastdmgattr):
    from operator import itemgetter

    cumbfd = [dict() for t in xrange(tsteps)]
    # for each training forecast, calculate basis functions (t+1 depends on t due to installed arcs)
    for nf in xrange(numforecasts):
        # basis function calculator requires a single sample, but will process for all time steps
        for t in xrange(tsteps):
            nx.set_edge_attributes(netobj.dinet, "t"+str(t)+"_dmg_pct", tfcastdmgattr[nf][t])
        mC["mcsim"].eval_adp_basis_functions(netobj, tsteps, mC)
        featurelist = mC["mcsim"].calc_feature_ctrmassquad(netobj, tsteps, features)
        # averaging theta*phi values
        for t in xrange(tsteps):
            for k,v in nx.get_edge_attributes(netobj.tnet, "t"+str(t)+"_sp_basis_val").iteritems():
                # i wish i knew how dictionaries are iterated...
                idx = mC["temp_var_insts"].index("TempArc"+str(k).replace(' ','') )
                if k not in cumbfd[t]:
                    cumbfd[t][k] = v*mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1][idx]
                else:
                    cumbfd[t][k] += v*mC["theta"][t][featurelist[t]["ctr_mass_quad_rel_frame"]-1][idx]
    rankedarcs = sorted(list(cumbfd[t].items()), key=itemgetter(1), reverse=True)
    return rankedarcs







def alt_control_step(tsteps, mC, netobj, samplesize, rankedarcs):
    banarcs = set([])
    policies = [ [0 for i in xrange(len(mC["temp_var_insts"])) ] for t in xrange(tsteps)]
    for t in xrange(tsteps):
        # select the pertinstbudget arcs corresponding to the largest average (over forecasts) theta*phi values
        considerpolicy = [0 for i in xrange(len(mC["temp_var_insts"] )) ]
        #ASSUMPTION: the max number (pertinstbudget) of arcs will be installed every time step. to do otherwise
        #            would cause unintended bugs in the following slice!
        toppicks = rankedarcs[t*mC["pertinstbudget"]:(t+1)*mC["pertinstbudget"]]
        for l,r in toppicks:
            # they're all being divided by the same thing, so i'm going to omit this for the sake of ease
##                cumbfd[t][k] /= float(numforecasts)
            # put the selections in the format of a binary vector (like mC["instchoices"])
            considerpolicy[mC["temp_var_insts"].index("TempArc"+str(l).replace(' ','') ) ] = 1
        # modify action via epsilon-greedy arc selection
        curpolicy = generate_per_arc_random_policy(mC["pertinstbudget"], mC["explorepolicy"], considerpolicy, \
                                                   banarcs)
        # apply selections to candidate forecast and record state -- transfer arc from tnet to dinet
        # extract arc attributes, modify to suit dinet, add to dinet, then remove from tnet
        for i in xrange(len(curpolicy)):
            if curpolicy[i]:
                arc = convert_cplex_structs(None, mC["temp_var_insts"][i])
                arcattr = netobj.tnet.get_edge_data(arc[0], arc[1])
                # there are extra temp entries but they're harmless (i.e. not used) in dinet
                arcattr["weight"] = 0
                arcattr["temp"] = False
                for tt in xrange(tsteps):
                    arcattr["t"+str(tt)+"_dmg_pct"] = 0.0
                    arcattr["t"+str(tt)+"_capacity"] = [arcattr["capacity"] for j in xrange(samplesize)]
                    arcattr["t"+str(tt)+"_sp_basis_val"] = 0.0
                # add arc to original infrastructure, but don't remove from temp: it causes unexpected errors,
                #   and changing its attributes this way essentially makes it invisible to temp net calculations
                #   EDIT: scratch that, it (may) cause errors to have it both in orig and temp, so leave as-is
##                    netobj.dinet.add_edge(arc[0], arc[1], arcattr)
##                    netobj.tnet.remove_edge(arc[0], arc[1])
                # now incorporate into state timeline
                policies[t][i] = 1
                banarcs.add(i)
                # remove arc from all future considerations
                for tt in xrange(t+1,tsteps):
                    policies[tt][i] = 1
    return policies













if __name__ == "__main__":
    import networkGenerator as ng
    #TODO8: refactoring should include creation of ADP param set
    runtime = time.time()
    tsteps = 3
    bigN = 300
    tmpflowcost = 0
    #NOTE: for now, slack cost is fine-tuned to be slightly above # temp arcs, but still enough such that 1 unit of
    #      slack flow is greater than installing all temp arcs (even though this is impossible for low pertinstbudget)
    slackcost = 10
    pertinstbudget = 7 #8
    # numfeats relates to cartesian quadrant of center of mass of dmg %s (NE[0], NW[1], SW[2], SE[3])
    features = ["ctr_mass_quad_rel_frame"]
    numfeats = 4
    epsgreedy = 0.03
    pmethod = "ud_vfa_separate"
    umethod = "grad_linear_direct_zeta_mcclainss"
    bft = "sp_bottleneckcapratio"
    tvfn = "IPv2a"
    timeTSAAscens = 10
    lAtTscens = 10
    savest = None#"an200_013_state_data_0p03eps_3T_300N"
    usest = None#"an200_013_state_data_0p03eps_3T_300N"#None
    nets = qpick("an200_013_network_data.pickle")
    nobj = ng.networkGenerator(dn=nets[0],tn=nets[1])
##    nobj.generate_random_network()
##    state = qpick("teststate.pickle")

    ## DEBUG stress testing via many loops of the main code
    numTrials = 1
    for dbg in xrange(numTrials):
        modelComponents = init_model_components(tsteps, slackcost, pertinstbudget, numfeats, epsgreedy, netobj=nobj,\
                                                polmethod=pmethod, updmethod=umethod, truevalfn=tvfn, bftype=bft, \
                                                tfc=tmpflowcost, ttsaas=timeTSAAscens, \
                                                lscens=lAtTscens, saveNet=False)

##        execute_alt_method(tsteps, modelComponents, bigN, features, recout="an200_013_alt_method_thetas")

        execute_ADP(tsteps, bigN, modelComponents, features, recout="an200_013_lin_0p03_upd_thetas", \
                    delout="an200_013_lin_0p03_upd_del_thetas", savestate=savest, usestate=usest)

    ##DEBUG end stress test loop

    print "done."
    print "total runtime: "+str(time.time() - runtime)+" seconds."
