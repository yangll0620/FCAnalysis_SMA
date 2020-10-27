import numpy as np
import pandas as pd
import sys
from statsmodels.stats.multitest import multipletests


codefolder = '/home/lingling/Desktop/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code'
sys.path.append(codefolder)

from m1_rest_calcciCOH import fc_select
from connAnalyTool.fc_visual_time import threshold_fc_overtime, ciCOH_visual_save, pvals_fc_overtime

file = '/home/lingling/Desktop/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/pipeline/NHPs/Pinky/1_dataAnaly/FCAnaly/Rest/areabased/py/m1_rest_calcciCOH/fc_rest_freq13_15.pickle'
freq = [13, 15]

fc = pd.read_pickle(file)


fcGraph_prefix = 'visualFC_rest' + '_freq' + str(freq[0]) + '_' + str(freq[1])
fc_lThaCortex = fc_select(fc, ['lVA', 'lVLo/VPLo', 'lSMA', 'rSMA','M1'])

ciCOH = fc_lThaCortex['ciCOH']['normal']
ntemp, ntrials, fs = fc_lThaCortex['setup']['ntemp_normal'], fc_lThaCortex['setup']['ntrials_normal'], fc_lThaCortex['setup']['fs']
f = (freq[0] + freq[1])/2
t = ntemp / fs

pvals = pvals_fc_overtime(ciCOH, ntrials, ntemp, f, t)


pvals_vec = []
nchns = pvals.shape[0]
for chi in range(nchns-1):
    for chj in range(chi + 1, nchns):
        pvals_vec.append(pvals[chi, chj])
pvals_vec = np.asarray(pvals_vec)

rejs, pvals_corr, alphacSidak, alphacBonf = multipletests(pvals_vec, alpha = 0.05, method = 'bonferroni')
rejs, pvals_corr, alphacSidak, alphacBonf = multipletests(pvals_vec, alpha = 0.05, method = 'fdr_bh')