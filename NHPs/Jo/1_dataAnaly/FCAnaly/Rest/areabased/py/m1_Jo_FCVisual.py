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

def fc_combdays_visual(files):
    
    lfpdata, chnAreas, fs = lfp_extract(files)


    ciCOHs = calc_ciCOHs_rest(lfpdata)


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
    filename = files[0]
    cond = re.search('_[a-z]*_[0-9]{8}', filename).group()[1:-9]
    freqstr = 'freq' + re.search('_filtered[0-9]*_[0-9]*', filename).group()[len('_filtered'):]

    save_prefix = 'all'
    texts = dict()
    texts[cond] = [-80, 50, 15]
    saveFCGraph = os.path.join(savefolder, freqstr + '_' + cond + '_' + save_prefix + '.png')
    weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                        savefile = saveFCGraph, texts = None, threds_edge = None)


def main():
    files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
    files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))
    files_moderate = glob.glob(os.path.join(inputfolder, '*_moderate_*'))


    fc_combdays_visual(files_normal)

    fc_combdays_visual(files_mild)

    fc_combdays_visual(files_moderate)


if __name__ == '__main__':
    main()