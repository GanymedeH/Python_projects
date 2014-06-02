import exeADP
import simulationEngine
import networkGenerator
import dataReports
import cplex
import sys
import time
import copy
import networkx as nx
from numpy import abs

def qpick(fname, mode="r", data=None):
    import cPickle
    f=open(fname,mode)
    if mode == "r":
        data = cPickle.load(f)
    elif mode == "w":
        cPickle.dump(data, f)
    f.close()
    return data






def run_LPs_per_forecast(numforecasts, mC, tsteps, adpftdmgattr, adpoutageattr, features, lpoutageattr, \
                         poleval="adp"):
    policyEvalModels = [None for nf in xrange(numforecasts)]
    LPSAAModel       = [None for nf in xrange(numforecasts)]
    lpobjval         = [None for nf in xrange(numforecasts)]
    polpicks         = [None for nf in xrange(numforecasts)]

    # find multiple policy evaluations for robustness
    for nf in xrange(numforecasts):
        # step 1: evaluation of policy
        if poleval == "alt":
            netobj = copy.deepcopy(mC["netobj"])
            samplesize = len(lpoutageattr[0].values()[0])
            rankedarcs = exeADP.alt_vfa_calculation(tsteps, mC, netobj, features, numforecasts, adpftdmgattr)
            policies = exeADP.alt_control_step(tsteps, mC, netobj, samplesize, rankedarcs)
            polpicks[nf] = policies[-1]
        elif poleval == "adp":
            policyEvalModels[nf] = exeADP.create_models_temp_framework(tsteps, mC)
            # override netobj's info with the stored scenario data
            for t in xrange(tsteps):
                nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(t)+"_dmg_pct", adpftdmgattr[nf][t])
                nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(t)+"_capacity", adpoutageattr[nf][t])

            # solve policy search problem to find install arcs
            mC["mcsim"].eval_adp_basis_functions(mC["netobj"], tsteps, mC)
            featurelist = mC["mcsim"].calc_feature_ctrmassquad(mC["netobj"], tsteps, features)

            if "ud_vfa_combined" in mC["polmethod"] or "node_constr" in mC["polmethod"]:
                exeADP.update_models_orig_scens(policyEvalModels[nf], mC, mrange=xrange(len(policyEvalModels[nf])), \
                                                            reset=True, updateConstr=True)
            elif "ud_vfa_separate" in mC["polmethod"]:
                exeADP.update_models_orig_scens(policyEvalModels[nf], mC, mrange=xrange(len(policyEvalModels[nf])), \
                                                            reset=True, updateConstr=False)
            else:
                print "run_LPs_per_forecast: invalid policy method specified."
                __dbg = raw_input("execution halted, press ENTER to exit")
                sys.exit()

            mC["instchoices"] = [ [0]*len(mC["temp_var_insts"]) for i in xrange(tsteps)]
            mC["previnsts"] = [[] for i in xrange(tsteps)]
            for t in xrange(tsteps):
                optpolicy,tflows,slks,bov,vfav = exeADP.adp_policy_search(t, featurelist, policyEvalModels[nf], mC,\
                                                                          randpol=False)
            # this is all we care about -- the policy evaluation at the last time step (i.e. all installed arcs)
            polpicks[nf] = optpolicy

            #DEBUG begin
##            idx = [ i for i in xrange(len(optpolicy)) if int(round(optpolicy[i]))]
##            instnames = [mC["temp_var_insts"][i] for i in idx]
##            ltdict = dict()
##            for ins in instnames:
##                crds = exeADP.convert_cplex_structs(cpobj=ins)
##                if mC["netobj"].tnet.edge[crds[0]][crds[1]]["linetype"] not in ltdict:
##                    ltdict[mC["netobj"].tnet.edge[crds[0]][crds[1]]["linetype"]] = 1
##                else:
##                    ltdict[mC["netobj"].tnet.edge[crds[0]][crds[1]]["linetype"]] += 1
##            print "\nLinetype counts for forecast "+str(nf)
##            for k,v in ltdict.iteritems():
##                print k+": "+str(v)
##            print "policy eval slacks: "+str(len([i for i in slks if int(round(i,4)) != 0]))
            #DEBUG end



        # step 2: fix the temp arcs to the decisions and run an LP SAA over these samples at time T
        LPSAAModel[nf] = exeADP.create_models_temp_framework(tsteps, mC)
        nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(tsteps-1)+"_capacity", lpoutageattr[nf])
        exeADP.update_models_orig_scens(LPSAAModel[nf], mC, numscens=len(lpoutageattr[nf].values()[0]), \
                                        mrange=[tsteps-1], reset=True, updateConstr=True)

        LPSAAModel[nf][tsteps-1].variables.set_types([ (mC["temp_var_insts"][j], \
                                                LPSAAModel[nf][tsteps-1].variables.type.continuous) \
                                       for j in xrange(len(mC["temp_var_insts"])) ])
        LPSAAModel[nf][tsteps-1].variables.set_lower_bounds( \
            zip(mC["temp_var_insts"], polpicks[nf]) )
        LPSAAModel[nf][tsteps-1].variables.set_upper_bounds( \
            zip(mC["temp_var_insts"], polpicks[nf]) )
        LPSAAModel[nf][tsteps-1].set_problem_type(LPSAAModel[nf][tsteps-1].problem_type.LP)
        try:
            LPSAAModel[nf][tsteps-1].solve()
            lpobjval[nf] = LPSAAModel[nf][tsteps-1].solution.get_objective_value()

            #DEBUG begin
##            slaks = LPSAAModel[nf][tsteps-1].solution.get_values(mC["slak_var_names"])
##            print "avg # of slacks: "+str(len([i for i in slaks if int(round(i,4)) != 0]) / \
##                                          float(len(lpoutageattr[nf].values()[0])) )
            #DEBUG end

        except:
            print "testADP: test LP SAA model infeasible"
            LPSAAModel[nf][tsteps-1].conflict.refine(LPSAAModel[nf][tsteps-1].conflict.all_constraints())
            conflicts = LPSAAModel[nf][tsteps-1].conflict.get()
            conflicts = [i for i in xrange(len(conflicts)) if conflicts[i] != -1]
            cgs = LPSAAModel[nf][tsteps-1].conflict.get_groups(conflicts)
            ubcs=[j[0][1] for i,j in cgs if j[0][0] == 2]
            lbcs=[j[0][1] for i,j in cgs if j[0][0] == 1]
            lccs=[j[0][1] for i,j in cgs if j[0][0] == 3]
            constrByVar = exeADP.find_constr_4_vars(LPSAAModel[nf][tsteps-1], \
                                                    LPSAAModel[nf][tsteps-1].variables.get_names())
            conflConstrs = exeADP.find_constr_4_vars(LPSAAModel[nf][tsteps-1], \
                                                     LPSAAModel[nf][tsteps-1].linear_constraints.get_names(lccs), \
                                                     vartype="constraint")
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()
    # this is INCREDIBLY slow, and there's tsteps*numforecasts models, so we'll sacrifice memory for now
##    exeADP.trash_model([k for l in LPSAAModel for k in l])
##    exeADP.trash_model([k for l in policyEvalModels for k in l])
    return polpicks, lpobjval






def run_IPs_per_forecast(numforecasts, mC, tsteps, ipoutageattr, useThisT=None):
    if useThisT == None:
        useThisT = tsteps-1
    IPSAAModel       = [None for nf in xrange(numforecasts)]
    ipobjval         = [None for nf in xrange(numforecasts)]
    sys.stdout.write("Processing forecast ")
    for nf in xrange(numforecasts):
        sys.stdout.write(str(nf)+".. ")
        IPSAAModel[nf] = exeADP.create_models_temp_framework(tsteps, mC)
        nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(useThisT)+"_capacity", ipoutageattr[nf])
        exeADP.update_models_orig_scens(IPSAAModel[nf], mC, numscens=len(ipoutageattr[nf].values()[0]), \
                                        mrange=[useThisT], reset=True, updateConstr=True)
        # restricting the number of arcs to install as in the training procedure
        #ASSUMPTION: even if we observe a problem at a time step other than T, we will allow a full install budget
        IPSAAModel[nf][useThisT].linear_constraints.add( \
            lin_expr = [cplex.SparsePair(ind = mC["temp_var_insts"], \
                                         val = [1]*len(mC["temp_var_insts"]) ) ], \
                        senses = ["L"], rhs = [mC["pertinstbudget"]*tsteps], names = ["PerTInstConstr"] )

        try:
            IPSAAModel[nf][useThisT].solve()
            ipobjval[nf] = IPSAAModel[nf][useThisT].solution.get_objective_value()

            #DEBUG begin
##            slks = IPSAAModel[nf][useThisT].solution.get_values(mC["slak_var_names"])
##            insts = [int(round(i)) for i in IPSAAModel[nf][useThisT].solution.get_values(mC["temp_var_insts"]) ]
##            idx = [ i for i in xrange(len(insts)) if int(round(insts[i]))]
##            instnames = [mC["temp_var_insts"][i] for i in idx]
##            ltdict = dict()
##            for ins in instnames:
##                crds = exeADP.convert_cplex_structs(cpobj=ins)
##                if mC["netobj"].tnet.edge[crds[0]][crds[1]]["linetype"] not in ltdict:
##                    ltdict[mC["netobj"].tnet.edge[crds[0]][crds[1]]["linetype"]] = 1
##                else:
##                    ltdict[mC["netobj"].tnet.edge[crds[0]][crds[1]]["linetype"]] += 1
##            print "\nLinetype counts"
##            for k,v in ltdict.iteritems():
##                print k+": "+str(v)
##            # len of capacities => # of scenarios in the SAA
##            print "avg # of slacks: "+str(len([i for i in slks if int(round(i,4)) != 0]) / \
##                                          float(len(ipoutageattr[nf].values()[0])) )
            #DEBUG end

        except:
            print "testADP: test IP SAA model infeasible"
    # this is INCREDIBLY slow, and there's tsteps*numforecasts models, so we'll sacrifice memory for now
##    exeADP.trash_model([k for l in IPSAAModel for k in l])
    return ipobjval






#TODO8: tere's a LOT of unnecessary processing/redundancy with run_LPs_per_forecast in here. strip it down to that
#       like run_one_big_IPSAA
def run_one_big_LPSAA(numforecasts, mC, tsteps, adpftdmgattr, adpoutageattr, features, lpoutageattr):
    policyEvalModels = [None for nf in xrange(numforecasts)]
    polpicks         = [None for nf in xrange(numforecasts)]
    lpobjval         = [None for nf in xrange(numforecasts)]

    LPSAAModel = exeADP.create_models_temp_framework(tsteps, mC)
    nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(tsteps-1)+"_capacity", lpoutageattr)
    exeADP.update_models_orig_scens(LPSAAModel, mC, numscens=len(lpoutageattr.values()[0]), mrange=[tsteps-1], \
                                    reset=True, updateConstr=True)
    LPSAAModel[tsteps-1].variables.set_types([ (mC["temp_var_insts"][j], \
                                                LPSAAModel[tsteps-1].variables.type.continuous) \
                                               for j in xrange(len(mC["temp_var_insts"])) ])
    LPSAAModel[tsteps-1].set_problem_type(LPSAAModel[tsteps-1].problem_type.LP)

    # find multiple policy evaluations for robustness
    for nf in xrange(numforecasts):
        # step 5: evaluation of policy
        policyEvalModels[nf] = exeADP.create_models_temp_framework(tsteps, mC)
        # override netobj's info with the stored scenario data
        for t in xrange(tsteps):
            nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(t)+"_dmg_pct", adpftdmgattr[nf][t])
            nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(t)+"_capacity", adpoutageattr[nf][t])

        mC["mcsim"].eval_adp_basis_functions(mC["netobj"], tsteps, mC)
        featurelist = mC["mcsim"].calc_feature_ctrmassquad(mC["netobj"], tsteps, features)
        exeADP.update_models_orig_scens(policyEvalModels[nf], mC, mrange=xrange(len(policyEvalModels[nf])), \
                                        reset=True, updateConstr=True)

        # solve policy search problem to find install arcs
        mC["instchoices"] = [ [0]*len(mC["temp_var_insts"]) for i in xrange(tsteps)]
        mC["previnsts"] = [[] for i in xrange(tsteps)]
        for t in xrange(tsteps):
            #TODO6: this may be semantically out of date
            optpolicy,tflows,slks,bov,vfav = exeADP.adp_policy_search(t, featurelist, policyEvalModels[nf], mC, \
                                                                       randpol=False)
        # this is all we care about -- the policy evaluation at the last time step (i.e. all installed arcs)
        polpicks[nf] = optpolicy

        # step 6: fix the temp arcs to the decisions made in 5 and run an LP SAA over these samples at time T
        LPSAAModel[tsteps-1].variables.set_lower_bounds( \
            [ (mC["temp_var_insts"][i], polpicks[nf][i]) for i in xrange(len(mC["temp_var_insts"])) ])
        LPSAAModel[tsteps-1].variables.set_upper_bounds( \
            [ (mC["temp_var_insts"][i], polpicks[nf][i]) for i in xrange(len(mC["temp_var_insts"])) ])
        try:
            LPSAAModel[tsteps-1].solve()
            lpobjval[nf] = LPSAAModel[tsteps-1].solution.get_objective_value()
        except:
            print "testADP: test LP SAA model infeasible"
            LPSAAModel[tsteps-1].conflict.refine(LPSAAModel[tsteps-1].conflict.all_constraints())
            conflicts = LPSAAModel[tsteps-1].conflict.get()
            conflicts = [i for i in xrange(len(conflicts)) if conflicts[i] != -1]
            cgs = LPSAAModel[tsteps-1].conflict.get_groups(conflicts)
            ubcs=[j[0][1] for i,j in cgs if j[0][0] == 2]
            lbcs=[j[0][1] for i,j in cgs if j[0][0] == 1]
            lccs=[j[0][1] for i,j in cgs if j[0][0] == 3]
            constrByVar = exeADP.find_constr_4_vars(LPSAAModel[tsteps-1], \
                                                    LPSAAModel[tsteps-1].variables.get_names())
            conflConstrs = exeADP.find_constr_4_vars(LPSAAModel[tsteps-1], \
                                                     LPSAAModel[tsteps-1].linear_constraints.get_names(lccs), \
                                                     vartype="constraint")
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()
    exeADP.trash_model(LPSAAModel)
    exeADP.trash_model([k for l in policyEvalModels for k in l])

    return polpicks, lpobjval






def run_one_big_IPSAA(numforecasts, mC, tsteps, ipoutageattr, samplesize):
    IPSAAModel = exeADP.create_models_temp_framework(tsteps, mC)
    nx.set_edge_attributes(mC["netobj"].dinet, "t"+str(tsteps-1)+"_capacity", ipoutageattr)
    exeADP.update_models_orig_scens(IPSAAModel, mC, numscens=numforecasts*samplesize, mrange=[tsteps-1], reset=True,\
                                    updateConstr=True)
    IPSAAModel[tsteps-1].linear_constraints.add( \
        lin_expr = [cplex.SparsePair(ind = mC["temp_var_insts"], \
                                     val = [1]*len(mC["temp_var_insts"])) ], \
                    senses = ["L"], rhs = [tsteps*mC["pertinstbudget"]], names = ["TInstConstr"] )

    try:
        IPSAAModel[tsteps-1].solve()
        vcontr = IPSAAModel[tsteps-1].solution.get_objective_value()
    except:
        print "testADP: control IP SAA model infeasible"
        __dbg = raw_input("execution halted, press ENTER to exit")
        sys.exit()

    exeADP.trash_model(IPSAAModel)
    sys.stdout.write("Done.\n")
    return vcontr








def generate_sample_data(mC, numforecasts, tsteps, samplesize, fixfc0=True, useThisLPT=None, useThisIPT=None):

    gennet = copy.deepcopy(mC["netobj"].dinet)
    # for policy evaluation, there are tsteps dicts, one for every t#_capacity. each dict maps an edge to a
    #    singleton-sized list of capacities
    adpoutageattr = [[dict() for t in xrange(tsteps)] for nf in xrange(numforecasts)]
    adpftdmgattr = [[dict() for t in xrange(tsteps)] for nf in xrange(numforecasts)]
    # these attr dicts will have samplesize-sized capacity values
    mcoutageattr = [[dict() for t in xrange(tsteps)] for nf in xrange(numforecasts)]
    mcftdmgattr = [[dict() for t in xrange(tsteps)] for nf in xrange(numforecasts)]
##    # for benchmarking, there is only 1 dict, but it maps an edge to a samplesize list of capacities
##    lpoutageattr = [dict() for nf in xrange(numforecasts)]
##    ipoutageattr = [dict() for nf in xrange(numforecasts)]
##    # this gathers some stats on the samples generated, currently not in use
##    outagepcts = dict(zip(mC["temp_var_insts"],[[] for i in xrange(len(mC["temp_var_insts"] )) ] ))
##    # the SAA dicts may be at any time we want (typically either tsteps-1 or 0)
##    if useThisLPT == None:
##        useThisLPT = tsteps-1
##    if useThisIPT == None:
##        useThisIPT = tsteps-1

    # specifies if we want the same forecast at time 0 for every separate forecast timeline
    if fixfc0:
        fixedepictr = mC["mcsim"].generate_distributions(gennet, 1)[0]
        fixedattr = nx.get_edge_attributes(gennet, "t0_dmg_pct")

    for nf in xrange(numforecasts):
        mC["mcsim"].clear_scenarios(gennet)
        if fixfc0:
            mC["mcsim"].generate_distributions(gennet, tsteps, sfcastattr=fixedattr, sfcastepictr=fixedepictr)
        else:
            mC["mcsim"].generate_distributions(gennet, tsteps)
        # generate samplesize+1 to get policy eval set and the rest for SAA test
        mC["mcsim"].sample_scenarios(gennet, ksamples=1)
        for t in xrange(tsteps):
            adpftdmgattr[nf][t] = nx.get_edge_attributes(gennet, "t"+str(t)+"_dmg_pct")
            adpoutageattr[nf][t] = nx.get_edge_attributes(gennet, "t"+str(t)+"_capacity")
        mC["mcsim"].clear_scenarios(gennet)
        mC["mcsim"].sample_scenarios(gennet, ksamples=samplesize)
        for t in xrange(tsteps):
            mcftdmgattr[nf][t] = nx.get_edge_attributes(gennet, "t"+str(t)+"_dmg_pct")
            mcoutageattr[nf][t] = nx.get_edge_attributes(gennet, "t"+str(t)+"_capacity")
##        lpftdmgattr[nf] = nx.get_edge_attributes(gennet, "t"+str(useThisLPT)+"_dmg_pct")
##        lpoutageattr[nf] = nx.get_edge_attributes(gennet, "t"+str(useThisLPT)+"_capacity")
##        ipftdmgattr[nf] = nx.get_edge_attributes(gennet, "t"+str(useThisIPT)+"_dmg_pct")
##        ipoutageattr[nf] = nx.get_edge_attributes(gennet, "t"+str(useThisIPT)+"_capacity")

##        for k,v in lpoutageattr[nf].iteritems():
##            if k not in outagepcts:
##                outagepcts[k] = []
##            outagepcts[k].append(float(len([i for i in v if i == 0])) / len(v) )
# return remnants: lpftdmgattr, ipftdmgattr, adpoutageattr, lpoutageattr, ipoutageattr
    return adpftdmgattr, adpoutageattr, mcftdmgattr, mcoutageattr






if __name__=="__main__":
    from numpy import inner
    # when comparing objectives, one should normalize the difference between objectives with the optimal choice
    # (i.e. performance ~= |V_adp - V_opt| / V_opt)

    # test0: sample M outage scenarios at time T from L forecast timelines, and get an additional sample over all
    #        time periods for each timeline. using the latter, run the ADP policy search as a policy evaluation to
    #        obtain arc installation decisions for each forecast timeline. fix temp arc capacities using these
    #        install decisions and run an LP SAA on the M scenarios for each timeline. these represent the quantity
    #        to benchmark. finally, run an IP for each of the L timelines where install decisions are to be made
    #        over the M scenarios of each timeline. this can be done at time 0 to present an upper bound on
    #        performance (it should: we have access to future data in ADP), and/or at time tsteps-1 to provide a
    #        lower bound (the IP knows the "true" forecast and has all available resources at that time) this is the
    #        control to benchmark against. most likely reporting metrics will be to take the average over the
    #        objectives of the different timelines.
    #
    # test1: sample a timeline forecast, a set of outage outcomes for each forecast and an outage forecast for the
    #        the last time step representing the "actual" damage. the adp process is executed for the first set of
    #        outcomes. then a min cost flow without the vfa coeffs is performed using the actual outages outcome.
    #        to maintain relatively similar-in-magnitude objectives, vfa decisions will not be considered in V_adp.

    # test2: sample a timeline forecast and a post-event outage outcome (sampled from the final time step's outage
    #        distribution). find install decisions per time step using the forecasts (not the outage sample) with
    #        the trained vfa (no further updating). determine unmet demand by running a deterministic min cost flow
    #        problem with fixed original and temp arc capacities (LP), whose objective no longer contains the vfa
    #        coeffs. compare the objective value with the objective value obtained by running the MCF problem with
    #        fixed (using the sample outages) orig arc capacities and constrained temp arc install decisions.
    #    Discussion: (1) policy evaluation can happen many ways in test1 (the easiest way), a single outcome is
    #                    sampled and run against. one could do other things, such as possibly run each time step as
    #                    a saa, or simply use the last time step with the vfa coeffs.
    #                (2) providing the control MILP with a cumulative install budget (tsteps*pertinstbudget) will
    #                    almost certainly reduce unmet demand to zero; at the very least, it will undoubtedly always
    #                    have lower objective value due to temp arc install costs. this raises the question: what
    #                    (if any) changes in params to the model should we make for the control? the above test is a
    #                    worst-case scenario, since the req's and costs for installing post-event are no worse/
    #                    larger than pre-event (which is almost never the case). The worst-case testing, however,
    #                    provides a means to observe relative error changes over iterations.
    #                (3) test2's objective values are magnitude (of vfa coeffs) dependent, so may not work well.

    # test3: eval policy as above, but compare its objective value to the objective of a sample average stochastic
    #        program using the last time step's outage distribution as the sample space.
    #        --> represents the policy's efficacy relative to the "best" mitigation (pre-event) decision

    # by invoking exeADP, we obtain the orig & temp network structures, and a historic record of vfa coeffs
    #ADDFEATURE: be able to customize adp params (likely by making exeADP object), e.g.:
    #            params below, data output methods, and param dicts for initializing network & simulation environ

    ## exeADP parameter specifications
    tsteps = 3
    bigN = 201
    tmpflowcost = 0
    slackcost = 10#200
    pertinstbudget = 7#8#1
    features = ["ctr_mass_quad_rel_frame"]
    numfeats = 4
    epsgreedy = 0.03
    pmeth = "ud_vfa_separate"
    bft = "sp_neg_bottleneckcapratio_wt=1_sf=0"
    umeth = "least_squares_linear_withvalfn_zeta_mcclainss"
    tvfn = "IPv2a"
    timeTSAAscens = 10
    lAtTscens = 10
    fileprefix = "NHC_"
    savtrast = fileprefix+"state_data_0p03eps_0p03mcclain_3T_201N.pickle"#None
    usetrainst = None#fileprefix+"state_data_0p03eps_0p03mcclain_3T_301N.pickle"

    ## testADP parameter specifications
    # how many times to repeat the ENTIRE process... we'll stick w/ 1 for now and let it be assumed below
    testtrials = 1
    numforecasts = 20
    samplesize = 50
    viewmodulo = 5
    cutoff = 200
    # the following are meant to provide filenames to output/input (resp.) sample data for use across runs
    savesamples = fileprefix+"test_samples.pickle"
    usesamples = None#fileprefix+"test_samples.pickle"

    ## BEGIN TEST0

    sys.stdout.write("Test initializing... ")
    runtime = time.time()
    nets = qpick(fileprefix+"network_data.pickle")

    # step 1: initialize adp model with forecast timeline to be used by the control model
    mC = exeADP.init_model_components(tsteps, slackcost, pertinstbudget, numfeats, epsgreedy, nxnet=nets[0], \
                                      tnet=nets[1], polmethod=pmeth, updmethod=umeth, truevalfn=tvfn, bftype=bft, \
                                      tfc=tmpflowcost, ttsaas=timeTSAAscens, lscens=lAtTscens, saveNet=False)
    sys.stdout.write("Done in "+str(time.time() - runtime)+" seconds\n")

    # step 2: train the adp model - this needs to be done only once across many analyses
    sys.stdout.write("Training ADP agent... ")
    runtime = time.time()
    exeADP.execute_ADP(tsteps, bigN, mC, features, \
                       recout=fileprefix+"sh_path_0p03mcclain_ss_0p03eps_lin_upd_thetas", savestate=savtrast, \
                       usestate=usetrainst)
    # recall recvs structure: v[n][t][f][i]
    recvs = qpick(fileprefix+"sh_path_0p03mcclain_ss_0p03eps_lin_upd_thetas.pickle")

    # there is a different policy evaluation every time we have different v's OR timelines
    policyeval = [[None for nf in xrange(numforecasts)] for vtfi in recvs]
    # and for each policy evaluation, we assess its value using an LP SAA at time tsteps-1
    vtest   = [[None for nf in xrange(numforecasts)] for vtfi in recvs]
    vcontr  = [None for nf in xrange(numforecasts)]
    sys.stdout.write("Done in "+str(time.time() - runtime)+" seconds\n")

    # step 3: generate and keep all scenarios
    sys.stdout.write("Generating timelines & scenarios... ")
    runtime = time.time()
    if usesamples != None:
        adpftdmgattr, adpoutageattr, mcftdmgattr, mcoutageattr = qpick(usesamples)
    else:
        adpftdmgattr, adpoutageattr, mcftdmgattr, mcoutageattr = \
            generate_sample_data(mC, numforecasts, tsteps, samplesize, fixfc0=False)
    if savesamples != None:
        qpick(savesamples, "w", (adpftdmgattr, adpoutageattr, mcftdmgattr, mcoutageattr) )
    sys.stdout.write("Done in "+str(time.time() - runtime)+" seconds\n")


    # step 4: construct and run an IP w/ all arcs uninstalled at time T using every sample
##    # only 1 objective across all v values, forecast timelines, and outage samples
##    sys.stdout.write("Calculating "+str(samplesize*numforecasts)+"-sample single SAA baseline... ")
##    runtime = time.time()
##    #TODO9: ipoutageattr is no longer 1 big dict, need to do the forecast mashup here, vis-a-vis:
##    for k,v in lpoutageattr[nf].iteritems():
##        if k not in ipoutageattr:
##            ipoutageattr[k] = []
##        ipoutageattr[k].extend(v)
##    # end example code
##    vcontr = run_one_big_IPSAA(numforecasts, mC, tsteps, ipoutageattr, samplesize)
##    sys.stdout.write("Done in "+str(time.time() - runtime)+" seconds\n")

    # use time 0 to indicate performance of a competing 2-stage model
    sys.stdout.write("Calculating "+str(samplesize)+"-sample per forecast IP SAA baselines... ")
    runtime = time.time()
    #ASSUMPTION: the "useThis" time overrides MUST be the same between generating outcomes and running them!!
##    useThisIPT = 0
##    ipoutageattr = [mcoutageattr[nf][useThisIPT] for nf in xrange(numforecasts)]
##    vcontr = run_IPs_per_forecast(numforecasts, mC, tsteps, ipoutageattr, useThisT=useThisIPT)
##    qpick(fileprefix+"vcontrol_T0.pickle", mode="w", data=vcontr)
##    # run both 1st and last time periods to serve as upper and lower bounds on LP optimality, respectively
##    ipoutageattr = [mcoutageattr[nf][-1] for nf in xrange(numforecasts)]
##    vcontr = run_IPs_per_forecast(numforecasts, mC, tsteps, ipoutageattr, useThisT=tsteps-1)
##    qpick(fileprefix+"vcontrol_Tfinal.pickle", mode="w", data=vcontr)

    sys.stdout.write("Done in "+str(time.time() - runtime)+" seconds\n")


    # step 5: gather stats for each iteration
    vncnt = 0
    sys.stdout.write("Benchmarking ADP training @ iteration ")
    runtime = time.time()
    useThisLPT = tsteps-1
    lpoutageattr = [mcoutageattr[nf][useThisLPT] for nf in xrange(numforecasts)]
    for vtfi in recvs:
        if vncnt % viewmodulo != 0 or vncnt > cutoff:
            vncnt += 1
            continue
        sys.stdout.write(str(vncnt)+".. ")
        # we must maintain the current policy vfa coefficients
        mC["theta"] = vtfi

        policyeval[vncnt], vtest[vncnt] = run_LPs_per_forecast(\
            numforecasts, mC, tsteps, adpftdmgattr, adpoutageattr, features, lpoutageattr, poleval="adp")
##        policyeval[vncnt], vtest[vncnt] = run_one_big_LPSAA(numforecasts, mC, tsteps, adpftdmgattr, \
##                                                               adpoutageattr, features, ipoutageattr)
        vncnt += 1
    sys.stdout.write("\n Benchmarking done in "+str(time.time() - runtime)+" seconds\n")

    qpick(fileprefix+"vtest.pickle", mode="w", data=vtest)
    qpick(fileprefix+"policyEvalOutages.pickle", mode="w", data=adpoutageattr)
    qpick(fileprefix+"policyInstallPicks.pickle", mode="w", data=policyeval)
    print "Done."





























        # an attempt at test1, now broken due to refactoring in exeADP
##        print "Iteration: "+str(itercnt)
##        itercnt += 1
##
##        # step 3a: sample forecast timeline and hypothetical outcomes for each time step (test1 discussion (1)).
##        mcomp["mcsim"].generate_distributions(netobj.dinet, tsteps)
##        mcomp["mcsim"].eval_adp_basis_functions(netobj, tsteps)
##        # use only 1 sample
##        exeADP.update_models_orig_scens(mcomp, reset=True)
##        flist = mcomp["mcsim"].calc_feature_ctrmassquad(netobj, tsteps, features=["ctr_mass_quad_rel_frame"])
##        models, orig_var_caps = exeADP.create_models_temp_framework(tsteps, mcomp)
##
##        mcomp["instchoices"] = [ [0]*len(mcomp["temp_var_insts"]) for i in xrange(tsteps)]
##        mcomp["previnsts"] = [[] for i in xrange(tsteps)]
##        baseobjvals = [None for i in xrange(tsteps)]
##
##        mcomp["theta"] = vtfi
##
##        # step 3b: execute trained adp model in time order
##        #NOTE: currently using a saa model for the last time step (better means of approximating performance?)
##        for t in xrange(len(vtfi)-1):
##            tmpflows,slacks,baseobjvals[t] = exeADP.adp_policy_search(t, flist, models, mcomp)
##
##        #TODO8: SAA part goes right here
##
##        # extract the objective value omitting vfa coefficients (note, still no original flow costs)
##        #NOTE: the overwritten output of adp_policy_search is ok, since we only want the last t's values
##        objadp = inner(mcomp["inst_var_costs"], mcomp["instchoices"][-1]) + inner(mcomp["flow_var_costs"], tmpflows) + \
##                 inner(mcomp["slak_var_costs"], slacks)
##
##        # step 3c: determine (sample) the "actual" post-event damages from the last time step's forecast
##        exeADP.update_models_orig_scens(modelComponents, reset=True)
##
##        # step 3d: solve model for actual damages w/o vfa coeffs
##        exeADP.trash_model(models)
##        models, orig_var_caps = exeADP.create_models_temp_framework(tsteps, mcomp)
##        models[-1].linear_constraints.add(lin_expr = [cplex.SparsePair(ind = mcomp["temp_var_insts"], \
##                                          val = [1]*len(mcomp["temp_var_insts"]) )], senses = ["L"], \
##                                          rhs = [tsteps*pertinstbudget], names = ["PerTInstConstr"] )
##
##        # step 3e: solve & store the objective (as the control value)
##        try:
##            models[-1].solve()
##            objideal = models[-1].solution.get_objective_value()
##        except:
##            print "testADP: cplex didn't work for control model"
##            break
##
##        # for test 2:
####        # step 3f: freeze the temp arc capacity constraints with the install decisions made in step 3b and run max-
####        #          flow LP (i.e. min cost slack)
####        models[-1].linear_constraints.delete("PerTInstConstr")
####        models[-1].objective.set_linear( \
####            [ (mcomp["temp_var_insts"][i], 0) for i in xrange(len(mcomp["temp_var_insts"])) ] + \
####            [ (mcomp["orig_var_names"][i], (1./mcomp["numscens"])*0) for i in xrange(len(mcomp["orig_var_names"])) ] + \
####            [ (mcomp["temp_var_flows"][i], 0) for i in xrange(len(mcomp["temp_var_flows"])) ] + \
####          [(mcomp["slak_var_names"][i], (1./mcomp["numscens"])*mcomp["slak_var_costs"][i]) for i in xrange(len(mcomp["slak_var_names"]))])
####        if any(mcomp["instchoices"][-1]):
####            models[-1].variables.set_lower_bounds([ (mcomp["temp_var_insts"][i], 1) for i in \
####                                                    xrange(len(mcomp["temp_var_insts"])) if mcomp["instchoices"][-1][i]])
####        if not any(mcomp["instchoices"][-1]):
####            models[-1].variables.set_upper_bounds([ (mcomp["temp_var_insts"][i], 0) for i in \
####                                                    xrange(len(mcomp["temp_var_insts"])) if not mcomp["instchoices"][-1][i]])
####
####        try:
####            models[-1].solve()
####            objpolicy = models[-1].solution.get_objective_value()
####        except:
####            print "testADP: cplex didn't work for test2-type policy evaluation"
####            break
##
##        # step 3g: performance = |V_adp - V_opt| / V_adp; store for iterate i, i = i+1, rinse and repeat.
##        if objadp < 1e-4:
##            perfvsiter.append(0.0)
##        else:
##            perfvsiter.append(abs(objadp-objideal)/float(objadp))
