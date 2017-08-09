
def renameGraph(graphName):
    if graphName == "RHG-100000000-nodes-2000000000-edges":
        return "rhg-100M-2G"
    if graphName == "rgg_n26_s0":
        return "rgg26"
    if graphName == "delaunay_n24":
        return "del24"
    return graphName
