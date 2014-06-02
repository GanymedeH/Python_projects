import sys
from numpy import sqrt
from operator import itemgetter
import networkGenerator as ngen

if __name__ == "__main__":
    arcf = open("NewC_Master_Arc_Query.txt")
    nodf = open("NewC_Nodes_Basic_Attributes_Query.txt")

    #ASSUMPTION: arccolns is ID_ARCs_Master, From_ID_Master_Node, To_ID_Master_Node, Infra_From, Infra_To,
    #                            Power_Capacity, Definition
    arccolns = arcf.readline().split()
    #ASSUMPTION: nodcolns is ID_Master_Node, Definition, Infras_Type, Power_Use, Node_Type, MinGX, MinGY
    nodcolns = nodf.readline().split()

    arcdata = dict([(int(ad.split()[0]), ad.split()[1:]) for ad in arcf.readlines()])
    noddata = dict([(int(nd.split()[0]), nd.split()[-5:]) for nd in nodf.readlines()])

    arcf.close()
    nodf.close()

    # prune any arcs that don't have any endpoints in the list of (approved) nodes
    prunekeys = [k for k,v in arcdata.iteritems() if int(v[0]) not in noddata or int(v[1]) not in noddata]
    for k in prunekeys:
        arcdata.pop(k)
    # prune any singleton nodes
    tonodes = [int(v[0]) for k,v in arcdata.iteritems()]
    fromnodes = [int(v[1]) for k,v in arcdata.iteritems()]
    prunekeys = [k for k,v in noddata.iteritems() if k not in tonodes and k not in fromnodes]
    for k in prunekeys:
        noddata.pop(k)

    # make lists for network generation

    # node keys: connection (hub/relay/terminus), demand (-supply,+demand), dists (list of node,dist tuples),
    #            pos, type (supply/trans/demand)
    #            connection: hub -> mostly outgoing connects to demands;
    #                        relay -> several incoming arcs but at least 1 outgoing; terminus -> all incoming arcs
    nodeids = []
    nodecxn = []
    nodedem = []
    nodedst = []
    nodepos = []
    nodetyp = []
    # shift the coordinates to "reasonable" numbers
    posxmin = sys.maxint
    posymin = sys.maxint
    for v in noddata.itervalues():
        if int(v[3]) < posxmin:
            posxmin = int(v[3])
        if int(v[4]) < posymin:
            posymin = int(v[4])

    for k,v in noddata.iteritems():
        nodeids.append(int(k))
        # determining connection is hard. wait until the data is in the networkx object to find it
        nodedem.append(int(v[1]))
        # distances can wait to be calculated until all nodes' positions are set
        nodepos.append( ( (int(v[3])-posxmin)/100, (int(v[4])-posymin)/100 ) )
        # type requires knowing if the node has incoming or outgoing arcs, unless demand is 0. need to wait tho

    # now calc dists
    for i in xrange(len(nodepos)):
        # do the round() so that we can ensure identical recalculations of distances
        gendists = [(nodeids[j], round(sqrt( (nodepos[i][0] - nodepos[j][0])**2  +  \
                                             (nodepos[i][1] - nodepos[j][1])**2), 8)) for j in xrange(len(nodepos))]
        nodedst.append(sorted(gendists, key=itemgetter(1)))


    # arc keys: capacity, linetype (main/supply/switch/feed/delivery), temp (T/F), weight (1/0)
    #           linetype: main -> hub to hub; supply -> supply to any; switch -> relay connects;
    #                     feed -> else to terminus; delivery -> hub to terminus;
    arcids = []
    arccap = []
    arclts = []
    arctmp = []
    arcwts = []
    for k,v in arcdata.iteritems():
        arcids.append( (int(v[0]),int(v[1])) )
        arccap.append(int(v[4]))
        # linetype needs to know which nodes its connecting, which isn't done until arcs are placed
        arctmp.append(False)
        arcwts.append(0)


    netobj = ngen.networkGenerator()
    netobj.construct_network(nodeids, nodecxn, nodedem, nodedst, nodepos, nodetyp, arcids, arccap, arclts, arctmp,\
                             arcwts, giveStructure=True)

    # the temp network is the exact original network, just set the temp and weight values accordingly
    # the copy operation is expensive, but it only happens once and i'm lazy
    netobj.tnet = netobj.dinet.copy()
    for i,j,d in netobj.tnet.edges_iter(data=True):
        d["temp"] = True
##        d["weight"] = 1

    netobj.save_network(filename="NewC_network_data.pickle")

    print "done."

