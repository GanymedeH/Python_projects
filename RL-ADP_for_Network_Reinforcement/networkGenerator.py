import networkx as nx
import random as rnd
import time
from math import floor
from math import ceil



def qpick(fname, mode="r", data=None):
    import cPickle
    f=open(fname,mode)
    if mode == "r":
        data = cPickle.load(f)
    elif mode == "w":
        cPickle.dump(data, f)
    f.close()
    return data







class networkGenerator:

    def __init__(self, dn=None, tn=None, h=[], re=[], su=[], d=[], p={}):
        if type(dn) == type(nx.DiGraph()):
            self.dinet=dn
            #NOTE: this performs a graph copy, which is very computationally expensive!
##            self.unnet=dn.to_undirected()
            self.tnet=tn
            self.hubs = [i for i,j in dn.nodes_iter(data=True) if j["connection"] == "hub" and j["type"] != "supply"]
            self.relays = [i for i,j in dn.nodes_iter(data=True) if j["connection"] == "relay"]
            self.suppliers = [i for i,j in dn.nodes_iter(data=True) if j["type"] == "supply"]
            self.demands = [i for i,j in dn.nodes_iter(data=True) if j["type"] == "demand"]
        else:
            self.hubs = h
            self.relays = re
            self.suppliers = su
            self.demands = d
            self.phonebook = p

            self.dinet = nx.DiGraph()
##            self.unnet = nx.Graph()
            self.tnet = nx.DiGraph()

        self.specify_params()

        # this is a kluge for switchpath creation in lieu of generically generating it. it conflicts with random
        #    demand node specification
        self.switchpaths = []
        return




    def export_netobj(self, fname):
        qpick(fname, mode="w", data=self)
        return




    def specify_params(self):
        self.params = dict()

        ###########################################
        ## Network node and general map parameters
        # this very basic test set assumes no intermediate temp arcs; i.e., no temp nodes
        self.params["nodes"] = 200#500#10                                          #DEBUG 500
        self.params["snodes"] = 5#10#1                                          #DEBUG 10
        #DEPRECATED: currently not in use; demand is chosen by process of elimination
##        self.params["dnodes"] = 7                                          #DEBUG 7
        # fix supply and demand values for now, but arcs can have different capacities
        self.params["slb"] = self.params["sub"] = 1000                     #DEBUG 1000/1000
        # the following forces demand to balance: demand/node = [total supply]/dnodes. if True, ignores dlb/dub
        self.params["demreltosup"] = True                                  #DEBUG True
        #DEPRECATED: same reason as dnodes
##        self.params["dlb"] = self.params["dub"] = 400                      #DEBUG 400/400
        # define limits for potential distance between nodes
        self.params["mapXlb"] = self.params["mapYlb"] = 0                  #DEBUG 0
        self.params["mapXub"] = self.params["mapYub"] = 1000#100                #DEBUG 1000

        ###########################################
        ## Original network edge placement parameters
        # arc capacities should scale with supply/demand requirements
        self.params["caplb"] = 400#400#200                                         #DEBUG 400
        self.params["capub"] = 1200#1200#800                                         #DEBUG 1200
        # often a randomly generated network will require augmentation to be feasible, resulting in a graph that
        #    provides the bare minimum in flow capacities along transhipment arcs. if an arc is found to require
        #    its maximum capacity at standard operating procedure, it will be increased by "cappad" %.
        self.params["cappad"] = 2.5#2.5#0.0                                        #DEBUG 2.5
        # hubs are chosen by how many neighbors it has which depends on the search radius. the number of neighbors
        #    within the radius is calculated for all nodes, then the list is traversed from largest to smallest,
        #    adding hubs until the number numhubs is reached. nodes that are within a certain distance from the
        #    newly chosen hubs are skipped to prevent multiple hubs connecting to many same stems. this hub number
        #    refers to the number of TRANS hubs only.
        self.params["neighborhood"] = 150#100#4                                    #DEBUG 100
        # this number does NOT take into account the number of suppliers that exist in the model. i.e., the total
        #    number of actual hubs is numhubs+snodes
        self.params["numhubs"] = 6#15#1                                         #DEBUG 15
        self.params["minhubdist"] = 180#150#5                                      #DEBUG 150
        # the number of stems need not be just those within the neighborhood; if not, manually specify. either
        #    method doesn't ensure that all nodes will be connected to the network; these nodes will just be
        #    connected to hubs closest to it.
        self.params["stemsinball"] = True#True#False                                  #DEBUG True
        # if above is false, the number of stems for a hub is determined stochastically uniformly between bounds.
        self.params["nslb"] = 8                                            #DEBUG 8
        self.params["nsub"] = 16                                           #DEBUG 16
        # if a specific number of stems is specified, we can still prevent far away nodes by disconnecting if a
        #    node is a max dist from the hub. set to 0 if there is no maximum.
        self.params["distmax"] = 30                                        #DEBUG 30
        # nodes also have the potential of being stems to multiple hubs. this can force them to be a stem of a
        #    fixed max number of stems (the closest hubs are chosen).
        self.params["multistemmax"] = 3#3#2                                    #DEBUG 3
        # in regards to hubs connecting to other hubs, choose bounds for the number of outgoing connections there
        #    are for each hub.
        self.params["outtohublb"] = 0                                      #DEBUG 0
        self.params["outtohubub"] = 1#1#0                                      #DEBUG 1
        # relays are nodes that connect a terminus of one hub to another hub. the switch path can be of variable
        #    length; i.e., a switch path of length 3 connects a terminus to another, a 3rd, and then back to a hub,
        #    making the actual number of edges 4 in the path (the code measures its dist by node for simplicity).
        #    individual lines of a switch path can be limited by geospatial length via "splinemaxdist".
        # specifying the number of switchpaths and the length of them provides a total number of relays picked; e.g.
        #    deterministically setting path number lower and upper bounds to 3 and their lengths to 4 will tag 12
        #    nodes as relays. the only requirement that is enforced above this is that there will remain at least 1
        #    demand node, even if it means shorting the number of paths or path lengths. due to this and some other
        #    switchpath options below, you may have more or less than the number and/or length you actually specify.
        self.params["numsplb"] = 10#16#2                                         #DEBUG 16
        self.params["numspub"] = 14#20#2                                         #DEBUG 20
        #TODO7: create a parameter that specifies that an "anchor" relay has to be within a certain dist to a hub
        # if too many paths and/or too long of a path is designated and there
        #    not enough remaining nodes, "strictsplb" doesn't allow for truncated paths (due to getting down to the
        #    last remaining node). instead, it will remove the truncated path and allow those nodes to be demands.
        self.params["splenlb"] = 6#12#2                                         #DEBUG 12
        self.params["strictsplb"] = False#False#True                                   #DEBUG False
        self.params["splenub"] = 16#24#3                                         #DEBUG 24
        self.params["splinemaxdist"] = 60#60#15                                  #DEBUG 60
        # if True, this lets you travel back to the same hub from which you started
        self.params["relbacktrack"] = True                                 #DEBUG True
        # the following designates whether relays are allowed in more than one switch path. setting this to True
        #    may have the unintended consequence of making more switchpaths (through the creation of cycles) than
        #    specified above. it also creates the potential for the same arcs to be a part of 2 different switch-
        #    paths (same arc internally though). Forbidding this behavior will have a slightly different effect.
        self.params["allowxpaths"] = True                                  #DEBUG True
        self.params["reusespaths"] = True#True#False                                  #DEBUG True
        # these specify the number of connections a supply/demand node may have. it will randomly select up to the
        #    appropriate number within the search radius. if override is true, then if there are not enough nodes
        #    to connect to within the search radius, connections will be made outside, nearest first. setting the
        #    search radius to 0 randomly selects any node. by definition, supply nodes are hubs. setting the lower
        #    bound to anything less than 1 is irrelevant, as all nodes must have at least 1 connection.
        self.params["shublb"] = 1
        self.params["shubub"] = 3
        # it's a good idea to make the number of supplies proportional to other non-demand nodes
        #DEBUG        self.params["shublb"] = floor((self.params["numhubs"]+self.params["numsplb"]) / 4.0)
        #DEBUG        self.params["shubub"] = ceil( ((self.params["numhubs"]+self.params["numsplb"]) / 3.0) )
        self.params["supradius"] = 80#80#3                                       #DEBUG 80
        self.params["overridesrad"] = False#False#True                                 #DEBUG False
        self.params["dstemlb"] = 1                                         #DEBUG 1
        self.params["dstemub"] = 2#2#1                                         #DEBUG 2
        self.params["demradius"] = 50#50#5                                       #DEBUG 50
        self.params["overridedrad"] = False                                #DEBUG False
        #DEPRECATED: currently demand nodes are defined as whatever is left from creating hubs & switchpaths
        #            (i.e., forcetermini is always True)
        ##self.params["forcetermini"] = True

        ###########################################
        ## Temp network edge placement parameters
        # temp arcs are typically built to provide as much flow as needed
        self.params["tlb"] = self.params["tub"] = self.params["snodes"]*self.params["sub"]
        # when creating temp arcs, specify a max distance from supply nodes and from relay nodes to hubs/relays
        #    to connect to
        self.params["slinkrad"] = 90#90#50                                       #DEBUG 90
        self.params["rlinkrad"] = 90#90#50                                       #DEBUG 90
        # what % of temp arcs to remove: reflects an imperfect support structure that can engender unmet demand
        #    unmet demand can also be produced by altering the pertinstbudget parameter in exeADP
        self.params["deltmppct"] = 92.0#96.0#0.0                                     #DEBUG 96.0

        return self.params




    def set_param(self, param, value):
        self.params[param] = value
        return




    def get_param(self, param):
        if param in self.params:
            return self.params[param]
        else:
            print "invalid parameter"
            return




    # This function builds the network that will eventually be damaged then repaired repeatedly as part of the
    #    algorithm. It's built in such a way that arcs are inclined to connect spatially closer cities than those that
    #    are far from each. Further, nodes that have many spatially close neighbors will be chosen as "hubs" that
    #    connect to its neighbors and to some number of other hubs. Note that by this definition, supply nodes must be
    #    hubs and demand nodes must be stems.
    #ADDFEATURE: allow parametrized a star's stem length (currently it's 1, meaning every leaf directly connects to
    #    the hub. also, allow minimum distances between supply and demand nodes (currently there is none).
    def generate_random_network(self):
        print "Generating random network... "
        self.place_nodes()
        self.classify_nodes()
        self.place_edges()
        self.ensure_connectedness()
        self.ensure_feasibility()
        return self.dinet#, self.unnet








    def place_nodes(self):
        from numpy import sqrt
        from operator import itemgetter

        #TODO9: parametrize geospatial placement of nodes
        # generate coordinates first, then calculate distances in the main loop
        nodelocs = [(rnd.uniform(self.params["mapXlb"],self.params["mapXub"]), \
                     rnd.uniform(self.params["mapYlb"],self.params["mapYub"]) ) for i in \
                                                                                xrange(self.params["nodes"]) ]
        for i in xrange(self.params["nodes"]):
            # distances will help determine which arcs to extract later
            # weird to have to recalculate distances, but it's slower to have 2 separate loops
            # do the round() so that we can ensure identical recalculations of distances
            gendists = [(j, round(sqrt( (nodelocs[i][0] - nodelocs[j][0])**2  +  \
                                  (nodelocs[i][1] - nodelocs[j][1])**2), 8) ) for j in \
                                                                           xrange(self.params["nodes"]) ]
            gendists = sorted(gendists, key=itemgetter(1))
            self.dinet.add_node(i, demand=0, type="", dists=gendists, pos=nodelocs[i], connection="")

##        self.unnet = self.dinet.to_undirected()
        return






    # this function takes placed nodes and designates them as supply/demand/hub/stem. it should be run after
    #    place_nodes() but before place_edges()
    def classify_nodes(self):
        from operator import itemgetter
        if self.dinet.number_of_nodes() < 1:
            print "run place_nodes first!"
            return

        self.hubs = []
        self.relays = []
        self.suppliers = []
        self.demands = []
        self.phonebook = {}

        # the following steps have been put into internal functions for readability - they are not meant to be
        #    invoked as a standalone function!
        ## step 2
        def count_neighbors(self):
            numneighbors = []
            for i,d in self.dinet.nodes_iter(data=True):
                if d["type"] != "supply" and d["type"] != "demand":
                    ncnt = []
                    for j,k in d["dists"][1:]:
                        if k <= self.params["neighborhood"]:
                            ncnt.append(j)
                        # distances are sorted, so once we find a node outside, the rest are
                        else:
                            break
                    self.phonebook[i] = ncnt
                    self.dinet.node[i]["neighbors"] = ncnt
                    numneighbors.append((i,len(ncnt)))
            return numneighbors

        ## step 3
        def designate_hubs(self, numneighbors):
            numneighbors = sorted(numneighbors, key=itemgetter(1), reverse=True)
            self.dinet.node[numneighbors[0][0]]["connection"] = "hub"
            self.dinet.node[numneighbors[0][0]]["demand"] = 0
            self.dinet.node[numneighbors[0][0]]["type"] = "trans"
            self.hubs.append(numneighbors[0][0])
            neighboridx = 1
            validhub = True
            while len(self.hubs) < self.params["numhubs"]:
                #NOTE: this may have the effect that there may be fewer than the desired # of hubs (if minhubdist is
                #       too large, for example). the only exception to this is that supply hubs can be close to trans
                if neighboridx >= len(numneighbors):
                    print "not enough nodes fit the criteria to be hubs. number of hubs in use: "+str(len(self.hubs))
                    break
                hubdists = dict(self.dinet.node[numneighbors[neighboridx][0]]["dists"])
                for j in self.hubs:
                    if hubdists[j] < self.params["minhubdist"]:
                        validhub = False
                        break
                if validhub:
                    self.dinet.node[numneighbors[neighboridx][0]]["connection"] = "hub"
                    self.dinet.node[numneighbors[neighboridx][0]]["demand"] = 0
                    self.dinet.node[numneighbors[neighboridx][0]]["type"] = "trans"
                    self.hubs.append(numneighbors[neighboridx][0])
                else:
                    validhub = True
                neighboridx += 1
            return

        ## step 4
        def designate_trans(self):
            availrelays = [j for j,k in self.dinet.nodes_iter(data=True) if k["connection"] != "hub"]
            # bundle the trans nodes so switchpaths will look reasonable
            for i in xrange(rnd.randint(self.params["numsplb"], self.params["numspub"]) ):
                nextrelay = rnd.choice(availrelays)
                availrelays.remove(nextrelay)
                # this is logging a list of switchpath lengths
                self.switchpaths.append(1)
                if nextrelay not in self.relays:
                    self.relays.append(nextrelay)
                    self.dinet.node[nextrelay]["type"] = "trans"
                    self.dinet.node[nextrelay]["demand"] = 0
                    self.dinet.node[nextrelay]["connection"] = "relay"
                # -2 because a switchpath of length n contains n-1 relays and we just picked one above
                for j in xrange(rnd.randint(self.params["splenlb"], self.params["splenub"])-1):
                    # look in specified radius first; failing that, use the relay's original neighborhood; if all
                    #    else fails, we still need to satisfy the switchpath length, so find nearest neighbor
                    relneighbors1 = [k for k,l in self.dinet.node[nextrelay]["dists"] if \
                                    l <= self.params["splinemaxdist"] and k in availrelays]
                    relneighbors2 = [k for k,l in self.dinet.node[nextrelay]["dists"] if \
                                    l <= self.params["neighborhood"] and k in availrelays]
                    relneighbors3 = [k for k,l in self.dinet.node[nextrelay]["dists"] if k in availrelays]
                    if any(relneighbors1):
                        nextrelay = rnd.choice(relneighbors1)
                    elif any(relneighbors2):
                        nextrelay = rnd.choice(relneighbors2)
                    else:
                        #TODO9: perhaps parametrize this 1:3
                        nextrelay = rnd.choice(relneighbors3[1:3])
                    availrelays.remove(nextrelay)
                    self.switchpaths[-1] += 1
                    if nextrelay not in self.relays:
                        self.relays.append(nextrelay)
                        self.dinet.node[nextrelay]["type"] = "trans"
                        self.dinet.node[nextrelay]["demand"] = 0
                        self.dinet.node[nextrelay]["connection"] = "relay"
                    # we need to ensure there is at least 1 demand; or if we prefer, make sure that the smallest
                    #    length switchpath is no lower than that which was specified.
                    if len(availrelays) < 2:
                        break
                if len(availrelays) < 2 or \
                   (self.params["strictsplb"] and len(availrelays) <= self.params["splenlb"]):
                    break
            return


        ## step 5
        def designate_demand(self, specnodes):
            #DEPRECATED: it doesn't make sense to manually specify how many relays and how many demands, since
            #            specifying one defaults the other to the rest. currently, demand is placed under default.
##            # demand nodes are fixed to a specified number. the others are forced to be relays
##            if not self.params["forcetermini"]:
##                # demand value depends on how much supply is available (to ensure we have enough for feasibility)
##                # we have to do it for each case of forcetermini because the number of demands are different
##                if self.params["demreltosup"]:
##                    totsup = 0
##                    for i in self.suppliers:
##                        totsup -= self.dinet.node[i]["demand"]
##                    self.params["dlb"] = floor(float(totsup*self.params["snodes"]) / len(self.params["dnodes"]))
##                    self.params["dub"] = floor(float(totsup*self.params["snodes"]) / len(self.params["dnodes"]))
##                for i in xrange(self.params["snodes"], len(specnodes) ):
##                    self.dinet.node[specnodes[i]]["type"] = "demand"
##                    self.dinet.node[specnodes[i]]["demand"] = rnd.randint(self.params["dlb"],self.params["dub"])
##                    self.dinet.node[specnodes[i]]["connection"] = "terminus"
##                    self.demands.append(specnodes[i])
##            else:

            # demand nodes are specified by whatever are termini, after supply/hubs/trans are specified
            self.demands = [i for i in self.dinet.nodes_iter() if i not in self.hubs and i not in self.relays \
                            and i not in self.suppliers]
            # demand value depends on how much supply is available (to ensure we have enough for feasibility)
            # we have to do it for each case of forcetermini because the number of demands are different
            if self.params["demreltosup"]:
                totsup = -1*sum(j["demand"] for i,j in self.dinet.nodes_iter(data=True) if i in self.suppliers)
                self.params["dlb"] = int(floor(float(totsup)/len(self.demands)))
                self.params["dub"] = int(floor(float(totsup)/len(self.demands)))
            for i in self.demands:
                self.dinet.node[i]["type"] = "demand"
                self.dinet.node[i]["demand"] = rnd.randint(self.params["dlb"],self.params["dub"])
                self.dinet.node[i]["connection"] = "terminus"
                #DEPRECATED: see above in the same function
##            # no matter the demand specification, trans is by definition whatever is left
##            for i in self.dinet.nodes_iter():
##                if i not in self.hubs and i not in self.demands and i not in self.suppliers:
##                    self.dinet.node[i]["type"] = "trans"
##                    self.dinet.node[i]["demand"] = 0
##                    self.dinet.node[i]["connection"] = "relay"
##                    self.relays.append(i)
            return



        ## begin classify_nodes main code
        # step 1: choose random nodes as supply.
        specnodes = rnd.sample(range(self.dinet.number_of_nodes()), self.params["snodes"])
        for i in xrange(self.params["snodes"]):
            self.dinet.node[specnodes[i]]["type"] = "supply"
            self.dinet.node[specnodes[i]]["demand"] = -1*rnd.randint(self.params["slb"],self.params["sub"])
            self.dinet.node[specnodes[i]]["connection"] = "hub"
            self.suppliers.append(specnodes[i])
        # step 2: count the number of neighbors of each trans node (supply nodes are by default hubs) but are NOT
        #         counted toward the "numhubs" (standard hubs) parameter
        numneighbors = count_neighbors(self)
        # step 3: designate hubs from distance/neighborhood criteria
        designate_hubs(self, numneighbors)
        # step 4: all other nodes are demands/trans
        designate_trans(self)
        # step 5: mark the rest demands
        designate_demand(self, specnodes)
        # step 6: make sure the sum of supplies equals the sum of demands: just distribute deficit evenly
        self.balance_demand()

##        self.unnet = self.dinet.to_undirected()
        return




    # step 6 of classify_nodes must be a class function so that step 6 of place_edges can use it
    ## step 6
    def balance_demand(self):
        totsup = sum(j["demand"] for i,j in self.dinet.nodes_iter(data=True) if j["demand"] < 0)
        totdem = sum(j["demand"] for i,j in self.dinet.nodes_iter(data=True) if j["demand"] > 0)
        if totsup != totdem:
            deficit = totsup+totdem
            idx = 0
            while deficit != 0:
                if idx == len(self.demands):
                    idx = 0
                # this means there is more supply than demand
                if deficit < 0:
                    self.dinet.node[self.demands[idx]]["demand"] += 1
                    deficit += 1
                # more demand than supply
                elif deficit > 0:
                    self.dinet.node[self.demands[idx]]["demand"] -= 1
                    deficit -= 1
                idx += 1
        return







    # once the nodes are in place and classified, this attempts to build a network structure akin to an expected
    #    network (i.e., one that is more likely to have arcs between geospatially close nodes, and that has few
    #    cycles) as much as possible, yet generically and randomly enough to be empirically sound.
    def place_edges(self):

        from operator import itemgetter
        if self.dinet.number_of_nodes() < 1:
            print "run place_nodes first!"
            return
        elif self.hubs == [] or self.suppliers == []:
            print "run classify_nodes first!"
            return


        # the following steps have been put into internal functions for readability - they are not meant to be
        #    invoked as a standalone function!
        ## step 2
        def connect_hubs_to_termini(self, i):
            if self.params["stemsinball"]:
                for j,k in self.dinet.node[i]["dists"][1:]:
                    if k > self.params["neighborhood"]:
                        break
                    if self.dinet.in_degree(j) <= self.params["multistemmax"] and \
                       self.dinet.node[j]["connection"] == "terminus":

                        self.dinet.add_edge(i,j, capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
                                            weight=0, temp=False, linetype="delivery")
            else:
                numstems = rnd.randint(self.params["nslb"], self.params["nsub"])
                j = 0
                stemcnt = 0
                while self.dinet.node[i]["dists"][j][1] <= self.params["distmax"] and self.params["distmax"] > 0 \
                      and stemcnt < numstems:
                    j += 1
                    if self.dinet.in_degree(self.dinet.node[i]["dists"][j][0]) <= self.params["multistemmax"] and \
                       self.dinet.node[self.dinet.node[i]["dists"][j][0] ]["connection"] == "terminus":

                        self.dinet.add_edge(i,self.dinet.node[i]["dists"][j][0], \
                                            capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
                                            weight=0, temp=False, linetype="delivery")
                        stemcnt += 1
            return


        ## step 3
        def create_switch_paths(self):

            def hookup_switchpaths(self, usedrelays, availrelays, allowxpaths=False, reusespaths=False, \
                                   relbacktrack=False):
                for sp in self.switchpaths:
                    currspath = []
                    # this finds the smallest distance of (hub,unused relay) arcs
                    closestrelay = (self.hubs[0], availrelays[0])
                    hubrelaydist = dict(self.dinet.node[self.hubs[0]]["dists"])[availrelays[0]]
                    for i in self.hubs:
                        hubdists = dict(self.dinet.node[i]["dists"])
                        for j in availrelays:
                            if hubdists[j] < hubrelaydist:
                                hubrelaydist = hubdists[j]
                                closestrelay = (i,j)
                    # if we wish, we can go back to the hub we started at (you need to if there is only 1 hub)
                    orighub = closestrelay[0]
                    # edge might not exist
                    if (closestrelay[0],closestrelay[1]) not in self.dinet:
                        self.dinet.add_edge(closestrelay[0],closestrelay[1], \
                                            capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
                                            weight=0, temp=False, linetype="switch")
                    # if we allow cross paths, keep track of what we added to this path so it won't loop back on itself
                    currspath.append(closestrelay[1])
                    # even if we allow cross paths, remove this "anchor" relay so that it can't be the start of the
                    #    next switchpath
                    #TODO8: this prevents the relay from being used ever again, even in a non-anchor position; fix
##                if not allowxpaths:
                    availrelays.remove(closestrelay[1])
                    # with the anchor relay hooked up, get the rest as counted in the node designation.
                    for i in xrange(1,sp):
                        for j,k in self.dinet.node[closestrelay[1]]["dists"][1:]:
                            # with cross paths, it may happen that we come upon a node with a pre-existing switchpath
                            if (closestrelay[1],j) in self.dinet and not reusespaths:
                                continue
                            if j in availrelays and j not in currspath:
                                self.dinet.add_edge(closestrelay[1],j, \
                                                    capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
                                                    weight=0, temp=False, linetype="switch")
                                closestrelay = (closestrelay[1], j)
                                currspath.append(j)
                                if not allowxpaths:
                                    availrelays.remove(j)
                                break
                        if not any(availrelays):
                            break
                    usedrelays.update(currspath)
                    # now hook it to the nearest hub
                    # if we want the possibility of coming back to the starting hub, it would by default be
                    #    chosen if the switchpath only has 1 relay, so give it a fair chance by randomly deciding
                    if relbacktrack and random.sample(self.hubs,1)[0] == orighub:
                        self.dinet.add_edge(closestrelay[1],orighub, \
                                            capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
                                            weight=0, temp=False, linetype="switch")
                    else:
                        for i,j in self.dinet.node[closestrelay[1]]["dists"][1:]:
                            if i in self.hubs and i != orighub:
                                self.dinet.add_edge(closestrelay[1],i, \
                                                    capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
                                                    weight=0, temp=False, linetype="switch")
                                break
                    if not any(availrelays):
                        break

            import copy
            import random
            #NOTE: for now, make this really simple: relay nodes are already labelled and there are no cross-paths,
            #      so just record the number and length of relays and connect them via nearest neighbor
            availrelays = copy.copy(self.relays)
            usedrelays = set([])
            while any(availrelays):
                hookup_switchpaths(self, usedrelays, availrelays, allowxpaths=self.params["allowxpaths"], \
                                   reusespaths=self.params["reusespaths"], relbacktrack=self.params["relbacktrack"])
                availrelays = list( set(self.relays) - usedrelays)

            # if we allowed cross paths, there may be some "dangling" relays that never got chosen. just hook them
            #    into the nearest 2 relays and/or hub(s)
##            if set(self.relays) != usedrelays:
##                relaystotag = list( set(self.relays) - usedrelays)

##            if self.params["allowxpaths"] and set(self.relays) != usedrelays:
##                relaystotag = list( set(self.relays) - usedrelays)
##                for rtt in relaystotag:
##                    # we still make the check in case the act of hooking 1 dangling relay partially connects another
##                    if self.dinet.in_degree(rtt) == 0:
##                        indegset = False
##                    else:
##                        indegset = True
##                    if self.dinet.out_degree(rtt) == 0:
##                        outdegset = False
##                    else:
##                        outdegset = True
##                    for j,k in self.dinet.node[rtt]["dists"][1:]:
##                        if j in self.hubs or j in self.relays:
##                            # avoid inadvertent dead ends as well
##                            if not outdegset and (j,rtt) not in self.dinet:
##                                self.dinet.add_edge(rtt,j, \
##                                                    capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
##                                                    weight=0, temp=False, linetype="switch")
##                                outdegset = True
##                            elif not indegset and (rtt,j) not in self.dinet:
##                                self.dinet.add_edge(j,rtt, \
##                                                    capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
##                                                    weight=0, temp=False, linetype="switch")
##                                indegset = True
##                            if outdegset and indegset:
##                                break
            return


        ## step 4 & step 5
        def connect_end_nodes(self, nodetype):

            def override_connections(self, connectto, onbors, numxions):
                xionsmade = 0
                for j in onbors:
                    if xionsmade >= numxions:
                        break
                    connectto.append(j)
                    xionsmade += 1
                return
            if nodetype == "supply":
                endnodes = self.suppliers
                lowbnd = self.params["shublb"]
                upbnd = self.params["shubub"]
                radius = self.params["supradius"]
                override = self.params["overridesrad"]
            else:
                endnodes = self.demands
                lowbnd = self.params["dstemlb"]
                upbnd = self.params["dstemub"]
                radius = self.params["demradius"]
                override = self.params["overridedrad"]

            for i in endnodes:
                # this is only creating connections
                numxions = rnd.randint(lowbnd, upbnd) - self.dinet.in_degree(i) - self.dinet.out_degree(i)
                if numxions < 1:
                    continue
                nbors = [(j,k) for j,k in self.dinet.node[i]["dists"][1:] if \
                               (i,j) not in self.dinet and (j,i) not in self.dinet and \
                               self.dinet.node[j]["connection"] != "terminus" and \
                               self.dinet.node[j]["type"] != "supply" ]
                nnbors = [j for j,k in nbors if k <= radius]
                nbors = [j for j,k in nbors if j not in nnbors]
                if radius:
                    if numxions <= len(nnbors):
                        connectto = rnd.sample(nnbors, numxions)
                    else:
                        connectto = nnbors
                    if  override:
                        override_connections(self, connectto, nbors, numxions-len(nnbors) )
                    # supply or demand, we have to connect to something no matter what
                    elif not any(nnbors):
                        override_connections(self, connectto, nbors, 1)
                else:
                    # don't care what we connect to so long as it's a relay or hub
                    if numxions <= len(nnbors+nbors):
                        connectto = rnd.sample(nnbors+nbors, numxions)
                    else:
                        connectto = nnbors+nbors

                for j in connectto:
                    #ASSUMPTION: supply arcs should be able to transmit full flow down any of its arcs, while feed
                    #    arcs need only hold as much as the demand requires
                    if nodetype == "supply":
                        self.dinet.add_edge(i,j, capacity=-1*self.dinet.node[i]["demand"], weight=0, temp=False, \
                                            linetype="supply")
                    else:
                        self.dinet.add_edge(j,i, capacity=self.dinet.node[i]["demand"], weight=0, temp=False, \
                                            linetype="feed")
            return


        ## step 6
        def connect_isolates(self):
            # if there are any isolates, they are hubs or relays; connect_end_nodes assures no supply/demands
            singles = nx.isolates(self.dinet)
            for s in singles:
                # again, distances are in increasing order, so the first hub/relay we find is the closest
                for i,j in self.dinet.node[s]["dists"][1:]:
                    if i in self.hubs or i in self.relays:
                        self.dinet.add_edge(i,s,capacity=rnd.randint(self.params["caplb"],self.params["capub"]),\
                                            weight=0, temp=False, linetype="feed")
                        # since this is the only edge, the isolate becomes a demand node
                        self.dinet.node[s]["type"] = "demand"
                        self.dinet.node[s]["demand"] = rnd.randint(self.params["dlb"],self.params["dub"])
                        self.dinet.node[s]["connection"] = "terminus"
                        self.demands.append(s)
                        self.balance_demand()
                        break
            return


        ## begin place_edges main code
        for i in self.hubs:

            # step 1: connect trans hubs (supply nodes are handled later)
            #         min is used in case someone specifies more connections than there are hubs
            tohubs = rnd.sample([j for j in self.hubs if j != i], min(rnd.randint(self.params["outtohublb"],
                                                                                 self.params["outtohubub"]), \
                                                                     len(self.hubs)-1) )
            for j in tohubs:
                self.dinet.add_edge(i,j, capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
                                    weight=0, temp=False, linetype="main")
            # step 2: connect hubs to termini
            connect_hubs_to_termini(self, i)
        # step 3: create switch paths
        # relays are nodes that connect a terminus of one hub to another hub. the switch path can be of
        #    variable length; i.e., a switch path of length 3 connects a terminus to another, then another, and
        #    finally to another hub. individual lines of a switch path can be limited by geospatial length. once a
        #    terminus is part of a switch path, it's connection type becomes "relay".
        create_switch_paths(self)
        # step 4: connect supply nodes
        connect_end_nodes(self, "supply")
        # step 5: connect demand nodes
        connect_end_nodes(self, "demand")
        #TODO9: make sure every switchpath connects to at least 1 demand node
        # step 6: ensure no isolates -- this shouldn't even be necessary, it's just a failsafe
        connect_isolates(self)

##        self.unnet = self.dinet.to_undirected()
        return






    def construct_network(self, nIDs, nConnect, nDemand, nDists, nPos, nType, aIDs, aCap, aLnType, aTmp, aWt, \
                          giveStructure=False):
        # place nodes
        nDicts = []
        for n in xrange(len(nIDs)):
            nd = {"demand":nDemand[n]}
            nd["dists"] = nDists[n]
            nd["pos"] = nPos[n]
            if not giveStructure:
                nd["connection"] = nConnect[n]
                nd["type"] = nType[n]
                if nType == "supply":
                    self.suppliers.append(nIDs[n])
                elif nType == "demand":
                    self.demands.append(nIDs[n])
                elif nConnect == "hub":
                    self.hubs.append(nIDs[n])
                elif nConnect == "relay":
                    self.relays.append(nIDs[n])
            nDicts.append(nd)
        self.dinet.add_nodes_from( [ (nIDs[i], nDicts[i]) for i in xrange(len(nIDs))] )

        # place edges
        aDicts = []
        for a in xrange(len(aIDs)):
            ad = {"capacity":aCap[a]}
            ad["temp"] = aTmp[a]
            ad["weight"] = aWt[a]
            if not giveStructure:
                ad["linetype"] = aLnType[a]
            aDicts.append(ad)
        self.dinet.add_edges_from( [ (aIDs[j][0], aIDs[j][1], aDicts[j]) for j in xrange(len(aIDs))] )

        # giving structure to the network involves setting node connection and type, and arc linetype
        if giveStructure:
            # as with random network generation, start with classifying nodes
            for n,d in self.dinet.nodes_iter(data=True):
                ideg = self.dinet.in_degree(n)
                odeg = self.dinet.out_degree(n)
                if (ideg == 0 or odeg == 0) and d["demand"] == 0:
                    print "supply/demand node has no flow! deleting node"
                    self.dinet.remove_node(n)
                    # need to remove it from the dists lists too
                    for nn,dd in self.dinet.nodes_iter(data=True):
                        for didx in xrange(len(dd["dists"])):
                            if n == dd["dists"][didx][0]:
                                break
                        dd["dists"] = dd["dists"][:didx] + dd["dists"][didx+1:]
##                elif ideg != 0 and odeg != 0 and d["demand"] != 0:
##                    print "trans node has demand value! structural parameters unset"
                elif ideg == 0 and d["demand"] != 0:
                    # suppliers have negative demand in the standard model
                    d["demand"] = -1*abs(d["demand"])
                    d["type"] = "supply"
                    d["connection"] = "hub"
                    self.suppliers.append(n)
                elif odeg == 0 and d["demand"] != 0:
                    d["demand"] = abs(d["demand"])
                    d["type"] = "demand"
                    d["connection"] = "terminus"
                    self.demands.append(n)
                #ASSUMPTION: a node is a hub if the ?three-fourths? majority of its connections are outgoing
                # force demand to be zero
                elif odeg > 9*ideg and ideg != 0:
                    d["demand"] = 0
                    d["type"] = "trans"
                    d["connection"] = "hub"
                    self.hubs.append(n)
                elif odeg != 0 and ideg != 0:
                    d["demand"] = 0
                    d["type"] = "trans"
                    d["connection"] = "relay"
                    self.relays.append(n)
                elif odeg == 0 and ideg == 0:
                    # remove from the network
                    self.dinet.remove_node(n)
                else:
                    print "typical node behavior not found! structural parameters unset"
            for i,j,d in self.dinet.edges_iter(data=True):
                if self.dinet.node[i]["connection"] == "hub" and self.dinet.node[j]["connection"] == "hub":
                    d["linetype"] = "main"
                elif self.dinet.node[i]["connection"] == "supply":
                    d["linetype"] = "supply"
                elif self.dinet.node[i]["connection"] == "hub" and self.dinet.node[j]["connection"] == "terminus":
                    d["linetype"] = "delivery"
                elif self.dinet.node[j]["connection"] == "terminus":
                    d["linetype"] = "feed"
                # this is meant to imply that the arc connects a hub or relay node to or from a relay node
                else:
                    d["linetype"] = "switch"
        # this shouldn't be needed, but just in case
        self.ensure_connectedness()
        self.ensure_feasibility()
        return








    def ensure_connectedness(self):
        # there are 2 types of connectedness we need: to make sure every demand has access to a supply, and to make
        #    sure the graph is weakly connected
        # we just want to make sure we can reach every demand from at least 1 supply, so create a super-source
##        self.dinet.add_edges_from([(self.params["nodes"], i) for i in self.suppliers])
        for i in self.hubs:
            if not any([nx.has_path(self.dinet, j, i) for j in self.suppliers]):
##            try:
##                sp = nx.shortest_path(self.dinet, source=self.params["nodes"],target=i)
##            except:
                # if not, hook the hub into a supply, or a hub or relay, but check for supply path to those cases
                for j,k in self.dinet.node[i]["dists"][1:]:
                    if self.dinet.node[j]["type"] == "supply":
                        self.dinet.add_edge(j,i, capacity=rnd.randint(self.params["caplb"],self.params["capub"]), \
                                            weight=0, temp=False, linetype="supply")
                        break
                    elif self.dinet.node[j]["connection"] == "relay" or self.dinet.node[j]["connection"] == "hub":
##                        try:
                        if any([nx.has_path(self.dinet, source=k ,target=j) for k in self.suppliers]):
                            # if no exception is thrown there, we're good
                            if self.dinet.node[j]["connection"] == "relay":
                                self.dinet.add_edge(j,i, \
                                                    capacity=rnd.randint(self.params["caplb"],self.params["capub"]),\
                                                    weight=0, temp=False, linetype="switch")
                            else:
                                self.dinet.add_edge(j,i, \
                                                    capacity=rnd.randint(self.params["caplb"],self.params["capub"]),\
                                                    weight=0, temp=False, linetype="main")
                            break
                        else:
                            continue
##                        except:
##                            continue
        # second step: ensure weak connectedness
        #TODO?: connectedness... pain in the butt and i might not need to do it.
        if not nx.is_weakly_connected(self.dinet):
            print "graph is disconnected!"
##            wcomps = nx.weakly_connected_component_subgraphs(self.dinet)
##            for wc in wcomps:
##                # connect the closest non-terminal nodes
##        self.dinet.remove_node(self.params["nodes"])
        return






    def ensure_feasibility(self):
        import cplex
        import exeADP

        # if we have to add temp arcs to have the original undamaged network be feasible, we need to give them
        #    proper bounds on flow
        origtemplb = self.params["tlb"]
        origtempub = self.params["tub"]
        self.params["tlb"] = self.params["caplb"]
        self.params["tub"] = self.params["capub"]
        # for now, we need a "complete" temp net
        savedelparam = self.params["deltmppct"]
        self.params["deltmppct"] = 0.0
        trialarcs = self.create_temp_net()
        feasible = False
        while not feasible:
            trialrun = cplex.Cplex()
            try:
                trialModel = exeADP.init_model_components(1,0,0,0,0,nxnet=self.dinet, tnet=trialarcs, saveNet=False)
                orig_var_caps  = [d["capacity"] for i,j,d in self.dinet.edges_iter(data=True)]
                ivc,ovc,fvc,tvi,tvc,ovn,tvf = trialModel["inst_var_costs"],trialModel["orig_var_costs"],\
                                              trialModel["flow_var_costs"],trialModel["temp_var_insts"],\
                                              trialModel["temp_var_caps"],trialModel["orig_var_names"],\
                                              trialModel["temp_var_flows"]

                trialrun.set_results_stream(None)
                trialrun.set_log_stream(None)
                trialrun.set_warning_stream(None)
                trialrun.variables.add(obj    = ivc + ovc + fvc,
                                        lb    = [0]*(len(tvi) + len(ovn) + len(tvf)), \
                                        ub    = [1]*len(tvi) + orig_var_caps + tvc, \
                                        types = "B"*(len(tvi))+"C"*(len(ovc) + len(fvc)),\
                                        names = tvi + ovn + tvf)
                # our "scenario" has no outages, but force it to 1 so that the constraints can be applied
                trialModel["numscens"] = 1
                exeADP.update_models_constraints(trialModel, [trialrun], mrange=[0], useSlack=False)

                trialrun.solve()
                decisions = trialrun.solution.get_values()
                #DEBUG begin
##                print "Generated network feasible?: "+str(trialrun.solution.is_primal_feasible())
                #DEBUG end
                feasible = True
            except:
##                trialrun.feasopt(trialrun.feasopt.all_constraints())
##                decisions = trialrun.solution.get_values()
                #DEBUG begin
##                print "Generated network feasible?: "+str(trialrun.solution.is_primal_feasible())
                #DEBUG end
                self.params["slinkrad"] += 1
                self.params["rlinkrad"] += 1
                trialarcs = self.create_temp_net()
                exeADP.trash_model([trialrun])

##                trialrun.conflict.refine(trialrun.conflict.all_constraints())
##                conflicts = trialrun.conflict.get()
##                conflicts = [i for i in xrange(len(conflicts)) if conflicts[i] != -1]
##                cgs = trialrun.conflict.get_groups(conflicts)
##                ubcs=[j[0][1] for i,j in cgs if j[0][0] == 2]
##                lbcs=[j[0][1] for i,j in cgs if j[0][0] == 1]
##                lccs=[j[0][1] for i,j in cgs if j[0][0] == 3]
##                constrByVar = exeADP.find_constr_4_vars(trialrun, trialrun.variables.get_names())
##                conflConstrs = exeADP.find_constr_4_vars(trialrun, trialrun.linear_constraints.get_names(lccs), \
##                                                  vartype="constraint")
                pass

        #DEBUG begin
##        print "Final radii for connecting temp arcs for suppliers: "+str(self.params["slinkrad"])
##        print "Final radii for connecting temp arcs for relays: "+str(self.params["rlinkrad"])
        #DEBUG end

        for i in xrange(len(tvi)):
            # if there's flow on it and it's a temp node, we've exposed a point of infeasibility
            if int(round(decisions[i])) > 0 or decisions[len(tvi+ovn)+i] > 0:
                if decisions[len(tvi+ovn)+i] <= 0:
                    continue

                # add it to the main network
                addedge = exeADP.convert_cplex_structs(trialrun, tvi[i])
                # temp arcs are redundant to original; if flow was put on it, it's in addition to the original flow
                if addedge in self.dinet:
                    #DEBUG begin
##                    print "augmenting "+tvi[i]+" with capacity "+ \
##                          str(int(ceil(decisions[len(tvi+ovn)+i])) )+" to the network"
                    #DEBUG end
                    self.dinet.edge[addedge[0]][addedge[1]]["capacity"] += int(ceil(decisions[len(tvi+ovn)+i]))
                    self.dinet.edge[addedge[0]][addedge[1]]["capacity"] = int( (1+self.params["cappad"]) * \
                                                        self.dinet.edge[addedge[0]][addedge[1]]["capacity"] )
                else:
                    #DEBUG begin
##                    print "adding "+tvi[i]+" with capacity "+ \
##                          str(int(ceil(decisions[len(tvi+ovn)+i])) )+" to the network"
                    #DEBUG end
                    self.dinet.add_edge(addedge[0],addedge[1], capacity=int( (1+self.params["cappad"]) * \
                                                                             ceil(decisions[len(tvi+ovn)+i])), \
                                        weight=0, temp=False, \
                                        linetype=trialarcs.edge[addedge[0]][addedge[1]]["linetype"])
        # put the flow bounds and del % back where they were
        self.params["tlb"] = origtemplb
        self.params["tub"] = origtempub
        self.params["deltmppct"] = savedelparam
##        self.unnet = self.dinet.to_undirected()
        return








    def complement_arc_capacities(netobj, supply, demand):
        import networkx as nx
        capval = max(sum(supply), sum(demand) )

        # this is very powerful, now all we have to do is add attributes and remove arcs that don't make sense
        netobj.tnet = nx.complement(netobj.dinet)
        for i in netobj.tnet.nodes_iter():
            netobj.tnet.node[i] = netobj.dinet.node[i]

        # make a copy of edges, otherwise the iteration breaks down (we're deleting dynamically)
        tempnet = netobj.tnet.copy()
        for i,j in tempnet.edges_iter():
            if netobj.tnet.node[i]["type"] == "trans" and netobj.tnet.node[j]["type"] == "supply":
                netobj.tnet.remove_edge(i,j)
            elif netobj.tnet.node[i]["type"] == "demand" and netobj.tnet.node[j]["type"] == "trans":
                netobj.tnet.remove_edge(i,j)
            elif netobj.tnet.node[i]["type"] == "demand" and netobj.tnet.node[j]["type"] == "supply":
                netobj.tnet.remove_edge(i,j)
            elif netobj.tnet.node[i]["type"] == "demand" and netobj.tnet.node[j]["type"] == "demand":
                netobj.tnet.remove_edge(i,j)
            elif netobj.tnet.node[i]["type"] == "supply" and netobj.tnet.node[j]["type"] == "demand":
                netobj.tnet.remove_edge(i,j)
            elif netobj.tnet.node[i]["type"] == "supply" and netobj.tnet.node[j]["type"] == "supply":
                netobj.tnet.remove_edge(i,j)
            elif i == j:
                netobj.tnet.remove_edge(i,j)
            else:
                netobj.tnet.edge[i][j] = {"capacity":capval, "weight":1, "temp":True}
        return netobj.tnet







    #TODO9: make a custom temp arc creator (pick specific connections, capacities, weights, costs, etc)
    #NOTE: temp arcs are designated as temp for all time, no matter what may have been chosen at a previous time.
    #       this works to our advantage: (1) there's no overhead in creating new "original" arcs; (2) the temp arcs
    #       with flow are ensured not to be wiped out; (3) all we have to do is increase the lower bound.
    #       We can worry about workgroup assignment later...
    def create_temp_net(self):
        tempnetwk = nx.DiGraph()
        tempnetwk.add_nodes_from(self.dinet.nodes(data=True))
        # 1st, add main lines
        for i in self.hubs:
            for j in self.hubs:
                if i != j:
                    tempnetwk.add_edge(i,j, capacity=rnd.randint(self.params["tlb"], self.params["tub"]), weight=1, \
                                  temp=True, linetype="main")
        # 2nd, connect nearby hubs/relays to relays
        for i in self.relays:
            for j in self.hubs+self.relays:
                if i != j and tempnetwk.node[j]["type"] != "supply" and \
                   dict(tempnetwk.node[i]["dists"])[j] <= self.params["rlinkrad"]:
                    tempnetwk.add_edge(j,i, capacity=rnd.randint(self.params["tlb"], self.params["tub"]), weight=1, \
                                  temp=True, linetype="switch")
        # 3rd, remove % of these edges to reflect an imperfect support network (making unmet demand possible)
        delthesearcs = rnd.sample(tempnetwk.edges(), \
                                  int(round(0.01*self.params["deltmppct"]*tempnetwk.number_of_edges() )) )
        tempnetwk.remove_edges_from(delthesearcs)
        # 4th, connect supplies to nearby hubs/relays
        for i in self.suppliers:
            for j in self.hubs+self.relays:
                if tempnetwk.node[j]["type"] != "supply" and dict(tempnetwk.node[i]["dists"])[j] <= self.params["slinkrad"]:
                    tempnetwk.add_edge(i,j, capacity=rnd.randint(self.params["tlb"], self.params["tub"]), weight=1, \
                                  temp=True, linetype="supply")
        # 4b: every supply must be connected to something, so make this happen
        for i in self.suppliers:
            if tempnetwk.out_degree(i) == 0:
                for j,k in tempnetwk.node[i]["dists"][1:]:
                    if (tempnetwk.node[j]["connection"] == "hub" and tempnetwk.node[j]["type"] != "supply"):
                        tempnetwk.add_edge(i,j, capacity=rnd.randint(self.params["tlb"], self.params["tub"]), weight=1,\
                                      temp=True, linetype="supply")
                        break
                    elif tempnetwk.node[j]["connection"] == "relay":
                        # this is valid, but we have to make sure it has a path to at least 1 hub
                        if any([nx.has_path(tempnetwk, source=j, target=l) for l in self.hubs]):
                            tempnetwk.add_edge(i,j, capacity=rnd.randint(self.params["tlb"], self.params["tub"]), \
                                          weight=1, temp=True, linetype="supply")
                            break
        # 5th, add a demand's connections as redundancies (in case those arcs go down)
        for i in self.demands:
            # just add 1, no redundancies needed
            die = self.dinet.in_edges(i)[0]
            tempnetwk.add_edge(die[0],die[1], capacity=rnd.randint(self.params["tlb"], self.params["tub"]), weight=1, \
                          temp=True, linetype=self.dinet.edge[die[0]][die[1]]["linetype"])
        return tempnetwk






    def plot_network(self, net="orig", savefig=False, figfile=None, toscreen=False):
        # show what the network looks like
        import matplotlib.pyplot as mpl
        if net == "orig":
            neto = self.dinet
        elif net == "temp":
            neto = self.tnet

        nposns = dict( [ (i, j['pos']) for i,j in neto.nodes_iter(data=True)] )
        ncolors = []
        ecolors = []
        for i,j in neto.nodes_iter(data=True):
            if j["type"] == "demand":
                ncolors += 'r'
            elif j["type"] == "supply":
                ncolors += 'g'
            elif j["connection"] == "terminus":
                ncolors += 'm'
            elif j["connection"] == "hub":
                ncolors += 'c'
            elif j["type"] == "trans":
                ncolors += 'w'
        if len(ncolors) != neto.number_of_nodes():
            print "you have stray node labels!"
        for i,j,k in neto.edges_iter(data=True):
            if k["linetype"] == "supply":
                ecolors += 'g'
            elif k["linetype"] == "feed":
                ecolors += 'r'
            elif k["linetype"] == "main":
                ecolors += 'c'
            else:
                ecolors += 'k'
        if len(ecolors) != neto.number_of_edges():
            print "you have stray edge labels!"

        mpl.figure()
        nx.draw_networkx(neto, pos=nposns, node_color=ncolors, edge_color=ecolors)
        mpl.axis([self.params["mapXlb"],self.params["mapXub"],self.params["mapYlb"],self.params["mapYub"]])
        mpl.title(net+" network")
        mpl.xlabel("x coordinate")
        mpl.ylabel("y coordinate")
        if savefig:
            if figfile:
                mpl.savefig(figfile, format="png")
            else:
                mpl.savefig(net+"_net_"+str(int(time.time()))+".png", format="png")
        if toscreen:
            mpl.show(block=False)
        return





    def save_network(self, filename=None, whichnet="all"):
        if not filename:
            import time
            filename = "dinet_tnet_run_"+str(time.time())+".pickle"
        if whichnet == "all":
            qpick(filename,"w",(self.dinet,self.tnet))
##        elif whichnet == "both_orig":
##            qpick(filename,"w",(self.dinet,self.unnet))
        elif whichnet == "orig":
            qpick(filename,"w",self.dinet)
        elif whichnet == "temp":
            qpick(filename,"w",self.tnet)
##        elif whichnet == "undirected":
##            qpick(filename,"w",self.unnet)
        return







if __name__ == "__main__":
##    nets = qpick("analysis_networks_3.pickle")
##    test = networkGenerator(dn=nets[0],tn=nets[1])
##    test.ensure_feasibility()

    for i in xrange(20):
        print "\ngenerating network pair %(cnt)03d" % {"cnt": i}
        test = networkGenerator()
        test.generate_random_network()

        runID = str(int(time.time()))
        test.plot_network(savefig=True, figfile="orig_net_%(cnt)03d.png" % {"cnt": i} )

        test.tnet = test.create_temp_net()
        test.plot_network(net="temp", savefig=True, figfile="temp_net_%(cnt)03d.png" % {"cnt": i} )

        test.save_network(filename="analysis_networks_%(cnt)03d.pickle" % {"cnt": i} )

    print "Done."