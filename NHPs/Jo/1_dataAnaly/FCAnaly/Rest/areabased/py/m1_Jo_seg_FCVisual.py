import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import re
import glob
import numpy as np
from scipy.stats import norm
from mne.stats import fdr_correction


codefolder = re.match('.*exp.*code', os.getcwd()).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder
from connAnalyTool.fc_visual_time import weight_visual_save


from lfp_ciCOH_extract import lfp_extract, calc_ciCOHs_rest
from statistical import pval_permciCOH_rest
from chnArea_set import assign_coord2chnArea



animal =  re.search('NHPs/[a-zA-Z]*/', os.getcwd()).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
corresfolder, correparentfolder = code_corresfolder(__file__)


freq = [24, 26]
inputfolder = os.path.join(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'Rest', 'm4_restData_filtered' + str(freq[0]) + '_' + str(freq[1]) + '_eqLen')


savefolder = corresfolder
area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')

def segfc_visual(onefile):
    
    # lfpdata: nchns * ntemp * nsegs
    lfpdata, chnAreas, fs = lfp_extract([onefile])

    nchns, _, nsegs = lfpdata.shape
    seg_ciCOHs = np.zeros(shape = (nchns, nchns, nsegs))
    mus = np.zeros(shape = (nchns, nchns, nsegs))
    stds = np.zeros(shape = (nchns, nchns, nsegs))
    shuffleN = 1000
    for segi in range(nsegs):
        
        seglfp = lfpdata[:, :, segi]
        seg_ciCOHs[:, :, segi] = calc_ciCOHs_rest(np.expand_dims(seglfp, axis = 2))

        
        # Fit a normal distribution to the data using permutation test
    for segi in range(nsegs):
        print("segi = " + str(segi))
        seglfp = lfpdata[:, :, segi]
        permlfp = seglfp.copy()
        permciCOHs = np.zeros(shape=(nchns, nchns, shuffleN))
        for shi in range(shuffleN):
            permlfp = np.transpose(permlfp, axes=(1, 0))
            np.random.shuffle(permlfp)
            permlfp = np.transpose(permlfp, axes=(1, 0))

            permciCOHs[:, :, shi] = calc_ciCOHs_rest(np.expand_dims(permlfp, axis = 2))
        
        for i in range(nchns -1):
            for j in range(i + 1, nchns):
                mus[i, j, segi], stds[i, j, segi] = norm.fit(permciCOHs[i, j, :])



    # permutation test: use the lfp data whose ciCOHs are the largest to get  distribution
    [i, j] = np.unravel_index(np.argmax(ciCOHs), shape = ciCOHs.shape)
    lfp1, lfp2 = lfpdata[i, :, :], lfpdata[j, :, :]
    _, mu, std = pval_permciCOH_rest(lfp1, lfp2, ciCOHs[i, j], shuffleN = 1000)
    pvals = norm.sf(abs(ciCOHs), loc = mu, scale = std) * 2


    # multiple comparison correction, get weights
    reject, pval_corr = fdr_correction(pvals, alpha = 0.05, method='indep')
    [rows, cols]= np.where(reject == True)
    weight = np.zeros(ciCOHs.shape)
    if len(rows) > 0:
        weight[rows, cols] = ciCOHs[rows, cols]


    # visual and save
    filename = os.path.basename(onefile)
    datestr = re.search('[0-9]{8}', filename).group()
    cond = re.search('_[a-z]*_[0-9]{8}', filename).group()[1:-9]

    save_prefix = 'all'
    saveFCGraph = os.path.join(savefolder, cond + '_' + save_prefix + '_' + datestr + '.png')
    weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                        savefile = saveFCGraph, texts = None, threds_edge = None)


def main():

    onefile = glob.glob(os.path.join(inputfolder, '*_mild_20151111_*'))[0]
    segfc_visual(onefile)


    # files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
    # files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))
    # files_moderate = glob.glob(os.path.join(inputfolder, '*_moderate_*'))


    # files_normal.sort()
    # dailyfc_visual(files_normal)

    # files_mild.sort()
    # dailyfc_visual(files_mild)


    # files_moderate.sort()
    # dailyfc_visual(files_moderate)


if __name__ == '__main__':
    main()