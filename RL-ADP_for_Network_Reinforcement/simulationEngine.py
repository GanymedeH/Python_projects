import exeADP
import networkx as nx
import random as rnd
from numpy import sqrt
from numpy import log
from numpy import prod

def qpick(fname, mode="r", data=None):
    import cPickle
    f=open(fname,mode)
    if mode == "r":
        data = cPickle.load(f)
    elif mode == "w":
        cPickle.dump(data, f)
    f.close()
    return data








class simulationEngine:
    def __init__(self):
        self.specify_params()
        return



    def specify_params(self):
        self.params = dict()
        # define distribution refining thus: start with an outage probability map, and apply a growth
        #    factor to a % of arcs and a shrink factor to the rest. rinse and repeat.

        # the origin epicenter is chosen at random from among the specified top percentile of nodes with the most
        #    neighbors. neighbors are counted within a "epi0srchrad" percentage of the map span.
        self.params["clustpctile"] = 0.55                                   #DEBUG 0.55/0.85
        self.params["epi0srchrad"] = 0.3
        # region size will be dictated by a % of an edge that lies within the radius of effect. the radius of
        #    effect is specified as a percentage of the map span-the longest distance between 2 nodes of the graph
        self.params["edgepct"] = 0.75
        self.params["radofeff"] = 0.7                                       #DEBUG 0.4/0.5
        # as the radius of effect, there is a radius of epicenter movement: in subsequent time steps, the epicenter
        #    can move a distance from the current one up to a certain distance away. the radius of shift is also
        #    given as a percentage of the map span.
        self.params["radofshift"] = 0.6                                     #DEBUG 0.3/0.3
        # to initialize, arcs will be given either outagechance1 (for likely targets) or outagechance0 (else)
        self.params["outagechance0"] = 0.02                                 #DEBUG 0.08/0.08
        self.params["outagechance1"] = 0.30                                 #DEBUG 0.16/0.16
        # in subsequent time steps, the bounds of uncertainty change by a multiplicative factor (positive means the
        #    uncertainty interval increases by that percentage), but the intervals are divided evenly among the
        #    ever-increasing number of prediction regions
        # TODO9: generalize this more... custom confidence interval changes, general gradient function, ...?
        self.params["unclbfactor"] = -1.0                                   #DEBUG -0.25/-0.25
        self.params["uncubfactor"] = -0.2                                   #DEBUG 0.16/0.16
##        # in subsequent time steps, outage chances are determined via an exponential growth/shrink by continually
##        #    multiplying by such a factor. growby is for arcs within the radius, shrink by is for all else
##        # all growth/shrink factors are randomly determined, varying between 0 and the value specified here. the
##        #    distribution is triangular with mode = the factor
##        self.params["growby"] = 0.5
##        self.params["shrinkby"] = 0.4
        #TODO9: potential for ros to exceed roe which doesn't always make sense, add flag for preference
        # in subsequent time steps, the radius of effect will decrease by yet another factor. these are stochastic:
        #    the actual shrink% is between 0 and this param, triangularly weighted to the full amount
        self.params["shrinkroeby"] = 0.7                                    #DEBUG 0.6/0.4
        # as we shrink the roe, we also shrink the radius for the maximum length the epicenter can shift over time
        self.params["shrinkrosby"] = 0.9                                    #DEBUG 0.8/0.27
        return






    def generate_distributions(self, nxnet, tsteps, sfcastattr=None, sfcastepictr=None):
        mapXlb = min(j['pos'][0] for i,j in nxnet.nodes_iter(data=True))
        mapXub = max(j['pos'][0] for i,j in nxnet.nodes_iter(data=True))
        mapYlb = min(j['pos'][1] for i,j in nxnet.nodes_iter(data=True))
        mapYub = max(j['pos'][1] for i,j in nxnet.nodes_iter(data=True))
        mapSpan = sqrt( (mapXub - mapXlb)**2 + (mapYub - mapYlb)**2 )
        roses = [self.params["radofshift"]]
        epiregions = [roses[0]*mapSpan]
        # to start, we include the default "unaffected" region to ensure all edges are handled
        roes = [self.params["radofeff"]]
        effregions = [mapSpan, roes[0]*mapSpan]

        # pick an epicenter that's in the thick of things; gives better chance of high arc exposure, more realistic
        busynodes = []
        for i,d in nxnet.nodes_iter(data=True):
            nodelist = []
            # start at 1 since the 1st dist measure is always to the node itself
            # since dists are sorted, once out of localarea, we're out for good
            for j,k in d["dists"][1:]:
                if k < self.params["epi0srchrad"]*mapSpan:
                    nodelist.append((j,k))
                else:
                    break
            busynodes.append((i, nodelist))

        # specifying the forecast % attributes at time 0 and the original epicenter are meant to be done together.
        #    you CAN specify either one or the other, but the resulting forecast timeline won't be "organic"
        if sfcastepictr != None:
            busynodes = dict(busynodes)
            epicenters = [sfcastepictr]
        else:
            # this ensures that only the top # of nodes (dictated by specified percentile) are considered "busy"
            busynodes = sorted(busynodes, key=lambda x: len(x[1]), reverse=True)
            busynodes = dict([ busynodes[i] for i in \
                               xrange(int(min(round((1-self.params["clustpctile"])*nxnet.number_of_nodes()), \
                                              len(busynodes)) )) ])
            # the epicenter can be in any of possibly several clusters, but once we pick one, the center stays in the
            #    general region for all time updates
            epicenters = [rnd.choice(busynodes.keys() )]

        # the "unaffected" region (region 0) can have any epicenter, do this for ease
        epicenters.append(epicenters[0])
        busynodes = busynodes[epicenters[0]]
        #DEBUG start
    ##    print "damage scenario at time 0 centered at node "+str(epicenter)
        #DEBUG end

        athome = False
        for i in xrange(1,tsteps):
            if not athome:
                # the next epicenter gets chosen from the nodes in the given grouping of the shrunken area
                busyholder = []
                for j in xrange(len(busynodes)):
                    nodelist = []
                    # start at 1 since the 1st dist measure is always to the node itself
                    # since dists are sorted, once out of localarea, we're out for good
                    for k in xrange(1,nxnet.number_of_nodes() ):
                        if nxnet.node[busynodes[j][0]]['dists'][k][1] < epiregions[-1]:
                            nodelist.append(nxnet.node[busynodes[j][0]]['dists'][k])
                        else:
                            break
                    busyholder.append((busynodes[j][0], nodelist))

                # randomly choose using a triangular distribution weighted to the full amount. if the region has
                #    shrunk to the point where there are no more neighbors close enough, keep the old epicenter
                if not any([l for k,l in busyholder if any(l)]):
                    #DEBUG start
##                    print "true epicenter found!"
                    #DEBUG end
                    athome = True
                    epicenters.append(epicenters[-1])
                else:
                    busynodes = sorted(busyholder, key=lambda x: len(x[1]), reverse=True)
                    if len(busynodes) == 1:
                        epicenters.append(busynodes[0][0])
                    else:
                        epicenters.append(busynodes[int(round(rnd.triangular(0,len(busynodes)-1,0) )) ][0] )
                    busynodes = dict(busynodes)[epicenters[-1]]

                    #DEBUG start
##                    print "damage scenario at time "+str(i)+" centered at node "+str(epicenters[-1])
                    #DEBUG end
            else:
                epicenters.append(epicenters[-1])

            # update rads of effect and shift
            roses.append(roses[-1]*(1-rnd.triangular(0,1,1)*self.params["shrinkrosby"]) )
            epiregions.append(roses[-1]*mapSpan)
            roes.append(roes[-1]*(1-rnd.triangular(0,1,1)*self.params["shrinkroeby"]) )
            effregions.append(roes[-1]*mapSpan)


            ##        for af, at, ad in nxnet.edges_iter(data=True):
            ##            if (af, at) in affectedarcs[0][0]:
            ##                ad["t0_dmg_pct"] = outageub
            ##            else:
            ##                ad["t0_dmg_pct"] = outagelb
##                        ad["t"+str(t)+"_dmg_pct"] = ad["t"+str(t-1)+"_dmg_pct"] *   \
##                                                                (1+rnd.triangular(0,1,1)*self.params["growby"])
        affectedarcs = self.extract_subnet(epicenters, effregions, nxnet)
        outageub = self.params["outagechance1"]
        outagelb = self.params["outagechance0"]

        if sfcastattr != None:
            nx.set_edge_attributes(nxnet, "t0_dmg_pct", sfcastattr)
            nx.set_edge_attributes(nxnet, "t0_capacity", dict(zip(nxnet.edges(), \
                                                                  [[] for i in xrange(nxnet.number_of_edges()) ] )))
            # rescale the outage chance interval to the next time step
            outagelb *= (1+self.params["unclbfactor"])
            outageub *= (1+self.params["uncubfactor"])
            startat = 1
        else:
            startat = 0

        for t in xrange(startat, tsteps):
            for af, at, ad in nxnet.edges_iter(data=True):
                # region 0 is the unaffected area - by definition, all edges are included in that region @"time -1"
                # the "+1" indicates that the first list is this degenerate list
                for regn in xrange(len(affectedarcs[t+1])):
                    if (af, at) in affectedarcs[t+1][regn]:
                        # calculate which interval to choose the random value
                        #NOTE: at time 0, we have 2 regions: outagechance0 (region 0) and outagechance1 (region 1)
                        unclb = outagelb + (outageub-outagelb)*regn/float(t+2)
                        uncub = unclb + (outageub-outagelb)/float(t+2)
                        ad["t"+str(t)+"_dmg_pct"] = rnd.uniform(unclb,uncub)
                        # provide initialized sample space
                        ad["t"+str(t)+"_capacity"] = []
                        # as per extract_subnet, an edge should be present in only 1 region
                        break
            # rescale the outage chance interval
            outagelb *= (1+self.params["unclbfactor"])
            outageub *= (1+self.params["uncubfactor"])
        return epicenters




    # netobj is specified separately from mC because we may wish to calculate basis functions on an arbitrary net
    def eval_adp_basis_functions(self, netobj, tsteps, mC, constrdecs=None, slackdecs=None):
        if "sp_" in mC["bftype"]:
##            cumdmgpct = 0.0
##            # it would be more appropriate to calculate the average slack cost among the arcs in the
##            #    shortest path only, but finding that relationship is nontrivial currently, and costs are
##            #    all the same (currently) anyway.
##            avgslakcost = sum(i for i in mC["slak_var_costs"]) / float( len(mC["slak_var_names"]) )
            #NOTE: while not directly needed, the default sp_ bftype is "sp_probabilistic"
            if "_probabilistic" in mC["bftype"]:
                pass
            elif "_bottleneckcapratio" in mC["bftype"]:
                maxflow = 0
                for af,at,ad in netobj.dinet.edges_iter(data=True):
                    maxflow = max(maxflow, ad["capacity"])
            else:
                print "eval_adp_basis_functions: invalid shortest path basis function sub-type."
                __dbg = raw_input("execution halted, press ENTER to exit")
                sys.exit()

            for t in xrange(tsteps):
                if not any(nx.get_edge_attributes(netobj.dinet, "t"+str(t)+"_dmg_pct") ):
                    continue

                for af,at,ad in netobj.dinet.edges_iter(data=True):
                    # find shortest path whose weights are the negative log of % an arc stays up: the chance a path
                    #    stays up is the product of chances for each arc, so taking the log lets us find the largest
                    #    path through adding edge weights. using the negative log of up% instead of the log of dmg%
                    #    lets us avoid negative edge weights, but more importantly quantifies the likelihood that at
                    #    least one arc in such a path is damaged (former) vs the likelihood that *every* arc in the
                    #    path is damaged (latter). we only need 1 arc in an existing path to be damaged to motivate
                    #    the installation of an arc, so it makes more sense to use the former.
                    netobj.dinet[af][at]["t"+str(t)+"_sp_wt"] = -1*log(1-netobj.dinet[af][at]["t"+str(t)+"_dmg_pct"])
                # find shortest path between nodes of every temp arc
                #ASSUMPTION: all temp nodes exist in the original network
                for af,at,ad in netobj.tnet.edges_iter(data=True):
                    # path might not exist
                    try:
                        path = nx.dijkstra_path(netobj.dinet, af, at, weight="t"+str(t)+"_sp_wt")
                        ad["t"+str(t)+"_sp_basis_val"] = 1 - \
                            prod([1-netobj.dinet[path[i]][path[i+1]]["t"+str(t)+"_dmg_pct"] \
                                  for i in xrange(len(path)-1)] )
                        if "_bottleneckcapratio" in mC["bftype"]:
                            # find bottleneck flow
                            if "wt=" in mC["bftype"]:
                                wt = float(mC["bftype"].split("wt=")[-1].split("_")[0])
                            else:
                                wt = 1
                            if "sf=" in mC["bftype"]:
                                sf = float(mC["bftype"].split("sf=")[-1].split("_")[0])
                            else:
                                sf = 0
                            minflow = min([netobj.dinet[path[i]][path[i+1]]["capacity"] for i in xrange(len(path)-1)])
                            bncapval = ((float(minflow) / maxflow) * wt)
                            ad["t"+str(t)+"_sp_basis_val"] *= (bncapval + sf*(1-bncapval))
                        # this is most commonly used for a separated policy search objective
                        if "_neg" in mC["bftype"]:
                            ad["t"+str(t)+"_sp_basis_val"] *= -1
                    except nx.NetworkXNoPath:
                        # if no path exists, we assume the probability of installing the temp arc (relative to
                        #    others) is 100%
                        ad["t"+str(t)+"_sp_basis_val"] = 1.
                        if "_bottleneckcapratio" in mC["bftype"]:
                            ad["t"+str(t)+"_sp_basis_val"] *= 100.
                        if "_neg" in mC["bftype"]:
                            ad["t"+str(t)+"_sp_basis_val"] *= -1
                        continue
                    except:
                        print "simulationEngine: misc. error in finding shortest path basis value"
                        break
##                    cumdmgpct += ad["t"+str(t)+"_sp_basis_val"]
##                # adjust order of magnitude to be representative of flow cost (to pair with theta that has
##                #    an order of magnitude comparitive to flow amount)
##                # k := scaling factor -- depends on magnitude of probabilities (below)
##                k = round(1.0 / (cumdmgpct/netobj.tnet.number_of_edges()), 1)
##                for af,at,ad in netobj.tnet.edges_iter(data=True):
##                    ad["t"+str(t)+"_sp_basis_val"] *= k*avgslakcost
        elif "ones" in mC["bftype"]:
            # simplest of them all: everything has a weight of 1
            for t in xrange(tsteps):
                for af,at,ad in netobj.tnet.edges_iter(data=True):
                    if "_neg" in mC["bftype"]:
                        ad["t"+str(t)+"_1s_basis_val"] = -1
                    else:
                        ad["t"+str(t)+"_1s_basis_val"] = 1
        elif "node_constr" in mC["bftype"]:
            # we actually just use the capacity. this is not technically the basis function (the whole node constr
            #    eval is), but when rearranging terms to consolidate arc decisions, this ends up being the case.
            for t in xrange(tsteps):
                for af,at,ad in netobj.tnet.edges_iter(data=True):
                    ad["t"+str(t)+"_nc_basis_val"] = ad["capacity"]
                # we also check to see if node constraint basis functions need updating
                if constrdecs != None and slackdecs != None:
                    for n,d in netobj.tnet.nodes_iter(data=True):
                        # search in-arcs and out-arcs. capacities for both types are negated, unless the destination
                        #    node is a trans node. i.e. add in-arcs for trans nodes
                        inarcs = netobj.tnet.in_edges(n)
                        outarcs = netobj.tnet.out_edges(n)
                        # find indices into constrdecs using mC["temp_var_insts"] as a map key
                        inidxs = []
                        outidxs = []
                        for dkidx in xrange(len(mC["temp_var_insts"])):
                            if exeADP.convert_cplex_structs(None, mC["temp_var_insts"][dkidx]) in inarcs:
                                inidxs.append(dkidx)
                            if exeADP.convert_cplex_structs(None, mC["temp_var_insts"][dkidx]) in outarcs:
                                outidxs.append(dkidx)
                        # if supply node, negate "demand" value (to convert it to a positive "supply" value)
                        if d["type"] == "supply":
                            d["t"+str(t)+"_nc_basis_val"] = -d["demand"]
                        else:
                            d["t"+str(t)+"_nc_basis_val"] = d["demand"]
                        iicnt = 0
                        for iaf,iat in inarcs:
                            # if trans node, add in-arc capacities
                            if d["type"] == "trans":
                                d["t"+str(t)+"_nc_basis_val"] += netobj.tnet.edge[iaf][iat]["capacity"] * \
                                                                  constrdecs[t][inidxs[iicnt]]
                            else:
                                d["t"+str(t)+"_nc_basis_val"] -= netobj.tnet.edge[iaf][iat]["capacity"] * \
                                                                  constrdecs[t][inidxs[iicnt]]
                            iicnt += 1
                        oocnt = 0
                        for oaf,oat in outarcs:
                            # always substract out-arc capacities
                            d["t"+str(t)+"_nc_basis_val"] -= netobj.tnet.edge[oaf][oat]["capacity"] * \
                                                                constrdecs[t][outidxs[oocnt]]
                            oocnt += 1
                        # if demand node, subtract slack contribution
                        if d["type"] == "demand":
                            didx = [i for i in xrange(len(mC["slak_var_names"])) \
                                    if "SlackForNode"+str(n)+"Scen" in mC["slak_var_names"][i]][0]
                            d["t"+str(t)+"_nc_basis_val"] -= slackdecs[t][didx]
        else:
            print "eval_adp_basis_functions: invalid basis function specification."
            __dbg = raw_input("execution halted, press ENTER to exit")
            sys.exit()
        return







    def calc_feature_ctrmassquad(self, netobj, tsteps, features=["ctr_mass_quad_rel_frame"]):
        featurelist = [{} for t in xrange(tsteps)]
        for t in xrange(tsteps):
            for f in features:
                xabslb = netobj.get_param("mapXlb")
                xabsub = netobj.get_param("mapXub")
                yabslb = netobj.get_param("mapYlb")
                yabsub = netobj.get_param("mapYub")
                if f == "ctr_mass_quad_rel_frame":
                    # calculate centroid weighted by dmg %
                    xctr = 0
                    yctr = 0
                    totdmgpct = 0
                    xrellb = xabsub
                    xrelub = xabslb
                    yrellb = yabsub
                    yrelub = yabslb
                    for af,at,ad in netobj.dinet.edges_iter(data=True):
                        # half the dmg % on each end node is equivalent to using full dmg % on the midpt
                        xctr += (netobj.dinet.node[af]["pos"][0] + netobj.dinet.node[at]["pos"][0]) * \
                                (ad["t"+str(t)+"_dmg_pct"]/2.)
                        yctr += (netobj.dinet.node[af]["pos"][1] + netobj.dinet.node[at]["pos"][1]) * \
                                (ad["t"+str(t)+"_dmg_pct"]/2.)
                        totdmgpct += ad["t"+str(t)+"_dmg_pct"]
                        xrellb = min(netobj.dinet.node[af]["pos"][0], netobj.dinet.node[at]["pos"][0], xrellb)
                        xrelub = max(netobj.dinet.node[af]["pos"][0], netobj.dinet.node[at]["pos"][0], xrelub)
                        yrellb = min(netobj.dinet.node[af]["pos"][1], netobj.dinet.node[at]["pos"][1], yrellb)
                        yrelub = max(netobj.dinet.node[af]["pos"][1], netobj.dinet.node[at]["pos"][1], yrelub)
                    xctr /= totdmgpct
                    yctr /= totdmgpct
                    # rotating from pos x-axis, each quadrant will contain its "starting" ray (eg -x axis in quad3)
                    if xctr > (xrellb+xrelub)/2. and yctr >= (yrellb+yrelub)/2.:
                        featurelist[t][f] = 1
                    elif xctr <= (xrellb+xrelub)/2. and yctr > (yrellb+yrelub)/2.:
                        featurelist[t][f] = 2
                    elif xctr < (xrellb+xrelub)/2. and yctr <= (yrellb+yrelub)/2.:
                        featurelist[t][f] = 3
                    elif xctr >= (xrellb+xrelub)/2. and yctr < (yrellb+yrelub)/2.:
                        featurelist[t][f] = 4
        return featurelist










    def sample_scenarios(self, nxnet, ksamples=1, trange=[]):
        # now that we have all the uncertainties determined, get a sample. i put this separate from the preamble
        #    and within the preceding loop a) for readability and b) in the hopes that this will be extracted
        #    to an external function

        # as iterations increase, storing every realization separately will make the data structure huge. this
        #    will be mitigated by having a clear_scenarios function, which must be used diligently!
        # for an arc's capacity value is keyed by t#_capacity and maps to a list indexed by scen #
        #ASSUMPTION: the number of scenarios is the same for all time steps
        for af, at, ad in nxnet.edges_iter(data=True):
            if not any(trange):
                i = 0
                # sampling over all t is a convenience, but this may cause unnecessary overhead when you just want
                #    samples for a single time step
                while "t"+str(i)+"_dmg_pct" in ad:
                    if "t"+str(i)+"_capacity" not in ad:
                        nxnet[af][at]["t"+str(i)+"_capacity"] = []
                    for k in xrange(ksamples):
                        if rnd.random() < ad["t"+str(i)+"_dmg_pct"]:
                            nxnet[af][at]["t"+str(i)+"_capacity"].append(0)
                        else:
                            nxnet[af][at]["t"+str(i)+"_capacity"].append(nxnet[af][at]["capacity"])
                    i += 1
            else:
                for t in trange:
                    if "t"+str(t)+"_capacity" not in ad:
                        nxnet[af][at]["t"+str(t)+"_capacity"] = []
                    for k in xrange(ksamples):
                        if rnd.random() < ad["t"+str(t)+"_dmg_pct"]:
                            nxnet[af][at]["t"+str(t)+"_capacity"].append(0)
                        else:
                            nxnet[af][at]["t"+str(t)+"_capacity"].append(nxnet[af][at]["capacity"])
        return len(nxnet[af][at]["t"+str(0)+"_capacity"])






    def clear_scenarios(self, nxnet, kvals=None):
        # default kvals wipes out all scenarios
        if kvals == None:
            for af, at, ad in nxnet.edges_iter(data=True):
                i = 0
                while "t"+str(i)+"_capacity" in ad:
                    nxnet[af][at]["t"+str(i)+"_capacity"] = []
                    i += 1
        #TODO9: clearing scenarios currently not implemented for anything else but all
    ##    else:
    ##        for af, at, ad in nxnet.edges_iter(data=True):
    ##            i = 0
    ##            while i < len(kvals) and "t"+str(i)+"_capacity" in ad:
    ##                nxnet[af][at]["t"+str(kvals[i])+"_capacity"]
    ##                i += 1
        return







    def extract_subnet(self, epicenters, effregions, nxnet):
        from numpy import sqrt
        from math import copysign
        if len(epicenters) != len(effregions):
            print "invalid data entry: number of regions don't match number of epicenters"
            return

        # every time step has arcs that are affected in regions of the same size as all previous time steps, plus
        #    a new one (a cumulative number of regions). this reflects an increase in granularity of the gradient
        #    of the outage distributions
        affectedarcs = [ [ [] for j in xrange(i+1)] for i in xrange(len(epicenters))]
        for tidx in xrange(len(affectedarcs)):
            for i,j in nxnet.edges_iter():
                edgedist = dict(nxnet.node[i]["dists"])[j]
                # the following equations assume the circle's center is at the origin, so shift everything
                ox,oy = nxnet.node[epicenters[tidx]]["pos"]
                ix,iy = nxnet.node[i]["pos"]
                ix -= ox
                iy -= oy
                jx,jy = nxnet.node[j]["pos"]
                jx -= ox
                jy -= oy
                dx = jx - ix
                dy = jy - iy
                detij = ix*jy - iy*jx
                # unfortunately we must repeat the process for previous regions due to time-changing epicenters
                for ridx in xrange(len(affectedarcs[tidx])-1, -1, -1):
                    # now find the edge length that coincides with the epiregion
                    discr = (effregions[ridx]**2)*(edgedist**2) - detij**2
                    if  discr <= 0:
                        # discriminant <= 0 => no overlap with the effected region
                        continue
                    # the values for x & y intersections aren't aligned, have to figure out how to pair them
                    # cross product == 0 => colinearity
                    rx1 = (detij*dy + copysign(dx*sqrt(discr),dy)) / (edgedist**2)
                    ry1 = (-1*detij*dx + abs(dy)*sqrt(discr)) / (edgedist**2)
                    rx2 = (detij*dy - copysign(dx*sqrt(discr),dy)) / (edgedist**2)
                    ry2 = (-1*detij*dx - abs(dy)*sqrt(discr)) / (edgedist**2)
                    if abs((jx-ix)*(ry1-iy) - (rx1-ix)*(jy-iy)) > 5e-8:
                        tmpy = ry1
                        ry1 = ry2
                        ry2 = tmpy
                    if dict(nxnet.node[epicenters[tidx]]["dists"])[i] <= effregions[ridx]:
                        if dict(nxnet.node[epicenters[tidx]]["dists"])[j] <= effregions[ridx]:
                            rx1 = ix
                            ry1 = iy
                            rx2 = jx
                            ry2 = jy
                        elif rx1 >= min(ix,jx) and rx1 <= max(ix,jx) and ry1 >= min(iy,jy) and ry1 <= max(iy,jy):
                            rx2 = ix
                            ry2 = iy
                        else:
                            rx1 = ix
                            ry1 = iy
                    elif dict(nxnet.node[epicenters[tidx]]["dists"])[j] <= effregions[ridx]:
                        if rx1 >= min(ix,jx) and rx1 <= max(ix,jx) and ry1 >= min(iy,jy) and ry1 <= max(iy,jy):
                            rx2 = jx
                            ry2 = jy
                        else:
                            rx1 = jx
                            ry1 = jy
                    effregdist = sqrt( (rx2-rx1)**2 + (ry2-ry1)**2 )
                    # if a certain % of the edge is contained in the region of effect, it is itself affected
                    # be careful: this finds points of intersections for the infinite line defined by the 2 points;
                    #    to determine if the edge is actually contained in the radius, the points of intersection
                    #    must either be a convex combination of the original 2 points (if at least 1 endpoint is
                    #    outside the radius) or both endpoints must be inside the radius
                    if effregdist >= self.params["edgepct"]*edgedist and rx1 >= min(ix,jx) and rx1 <= max(ix,jx) \
                       and ry1 >= min(iy,jy) and ry1 <= max(iy,jy) and rx2 >= min(ix,jx) and rx2 <= max(ix,jx) and\
                       ry2 >= min(iy,jy) and ry2 <= max(iy,jy):
                        affectedarcs[tidx][ridx].append((i,j))
                        # we're going from smallest region to largest: if an arc is in a smaller region, it is
                        #    subject to the most recent effect and so should be excluded from other regions
                        break
        return affectedarcs







## deprecated code selects in-region arcs by adjacencies, regardless of physical orientations
    ##        affectedarcs = [[] for i in xrange(len(localareas))]
    ##        usedupnodes = []
    ##        dbgcount = 0
    ##        # this extracts the distances between nodes for each node
    ##        alldists = dict([ (n, ad['dists']) for n,ad in unnet.nodes_iter(data=True)])
    ##
    ##
    ##        while len(affectedarcs) < round(currpctedges*nxnet.number_of_edges()):
    ##            #DEBUG start
    ##            if dbgcount > 100000 or len(usedupnodes) >= nxnet.number_of_nodes():
    ##                print "something's wrong in extract_subnet"
    ##                break
    ##            #DEBUG end
    ##
    ##            for an in unnet.edges(nodenow):
    ##                if len(affectedarcs) >= round(currpctedges*nxnet.number_of_edges()):
    ##                    break
    ##                # make sure we include the appropriate directed arcs
    ##                if an in nxnet.edges() and an not in affectedarcs:
    ##                    affectedarcs.append(an)
    ##                if (an[1], an[0]) in nxnet.edges() and (an[1], an[0]) not in affectedarcs and \
    ##                   len(affectedarcs) < round(currpctedges*nxnet.number_of_edges()):
    ##                    affectedarcs.append((an[1], an[0]) )
    ##
    ##            usedupnodes.append(nodenow)
    ##
    ##            # this picks the next nodenow
    ##            # grab all the nodes the are within the localarea and pick one randomly (that hasn't already been
    ##            #    chosen). If all within the localarea have been chosen, take the next closest node outside it.
    ##            localnodes = set([i for i,j in alldists[nodenow] if j <= localareas])
    ##            localnodes = localnodes - set(usedupnodes)
    ##            if localnodes:
    ##                nodenow = rnd.choice(list(localnodes))
    ##            else:
    ##                for nn,nd in alldists[nodenow]:
    ##                    # list is sorted, just grab the 1st node past localarea that hasn't been picked before
    ##                    if nd > localareas and nn not in usedupnodes:
    ##                        nodenow = nn
    ##                        break
    ##
    ##            dbgcount += 1






if __name__=="__main__":
    import networkGenerator as ngen
    testnet = ngen.networkGenerator(qpick("dinet_run_1332283427_orig.pickle")[0])
    seng = simulationEngine()
    seng.generate_distributions(testnet.dinet, 6)
    pass


