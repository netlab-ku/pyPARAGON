
from Paragon.GraphletFrequency import GraphletSelect
from Paragon.Graphlets_lite import GraphletSearch 
from Paragon.Flux import PageRankFlux   
from Paragon.NetworkAnalysis import Community_Analysis

def GraphletFrequency(network, Node_list):
    return GraphletSelect(network, Node_list)


def GraphletGuidance(network,node_list=None):
    return GraphletSearch(network,node_list=node_list)


def NetworkInference(network,guide_network=None,edge_attribute=None):
    if guide_network==None:
        return PageRankFlux(network=network,motif_network=network,edge_attribute=edge_attribute)
    else:
        return PageRankFlux(network,motif_network=guide_network,edge_attribute=edge_attribute)
def CommunityAnalysis(subnetwork):
    return Community_Analysis(subnetwork)
