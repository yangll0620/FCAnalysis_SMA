import os
import sys
import scipy.io as sio
import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle 
import math
import pandas as pd
import cv2
from statsmodels.stats.multitest import multipletests
import networkx as nx

folder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/pipeline/NHPs/Pinky/1_dataAnaly/FCAnaly/Rest/areabased/py/m1_rest_Pinky_FCVisual'
fc = pd.read_pickle(os.path.join(folder, 'fc_rest_freq26_28.pickle'))


lowweight = 0.1

ciCOH = fc['ciCOH']['normal']

weight = abs(ciCOH)
weight[np.where(weight<lowweight)] = 0 # weight > lowweight



G = nx.Graph()
G.add_nodes_from(np.arange(1, 17))
for i in range(0, weight.shape[0] - 1):
    for j in range(i+1, weight.shape[0]):
        if weight[i,j] > 0:
            G.add_edge(i, j, weight= weight[i, j])

print(list(G.edges))
print(nx.average_clustering(G))
print(nx.degree_centrality(G))

nx.degree_centrality


