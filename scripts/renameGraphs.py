
def renameGraph(graphName):
    if graphName == "RHG-100000000-nodes-2000000000-edges":
        return "rhg"
    if graphName == "rgg_n26_s0":
        return "rgg26"
    if graphName == "delaunay_n24":
        return "del24"
    return graphName
