import numpy as np
import networkx as nx



def fcnetwork_avgCC(weight):
    """
        calculate fc network avg CC

        Arg:
            weight: weight matrix (nchns * nchns)
    """
    nchns = weight.shape[0]
    G = nx.Graph()
    G.add_nodes_from(np.arange(1, nchns + 1))
    for i in range(0, nchns - 1):
        for j in range(i+1, nchns):
            if abs(weight[i,j]) > 0:
                G.add_edge(i, j, weight= weight[i, j])

    avg_CC = nx.average_clustering(G)
    del G

    return avg_CC