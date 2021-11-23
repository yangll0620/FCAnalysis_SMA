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
from subArea import ciCOH_select



animal =  re.search('NHPs/[a-zA-Z]*/', os.getcwd()).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
corresfolder, correparentfolder = code_corresfolder(__file__)


freq = [26, 28]
inputfolder = os.path.join(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'Rest', 'm4_restData_filtered' + str(freq[0]) + '_' + str(freq[1]) + '_eqLen')


savefolder = corresfolder
area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')

def subArea_dailyfc_visual(files):
    
    for onefile in files:
        lfpdata, chnAreas, fs = lfp_extract([onefile])

        if lfpdata.shape[2] < 80:
            continue


        print(onefile)
        ciCOHs = calc_ciCOHs_rest(lfpdata)




        # permutation test: use the lfp data whose ciCOHs are the largest to get  distribution
        [i, j] = np.unravel_index(np.argmax(ciCOHs), shape = ciCOHs.shape)
        lfp1, lfp2 = lfpdata[i, :, :], lfpdata[j, :, :]
        _, mu, std = pval_permciCOH_rest(lfp1, lfp2, ciCOHs[i, j], shuffleN = 1000)


        cond = re.search('_[a-z]*_[0-9]{8}', files[0]).group()[1:-9]
        datestr = re.search('[0-9]{8}', os.path.basename(onefile)).group()


        ### left thalamus and SMA/M1 ###
        save_prefix = 'leftThaCor_' 
        areas_used = ['lVA', 'lVLo/VPLo', 'lSMA', 'rSMA','M1']

        # subareas selection
        ciCOH_new, chnAreas_new = ciCOH_select(ciCOHs, chnAreas, areas_used)
        
        
        # multiple comparison correction, get weight matrix
        pvals = norm.sf(abs(ciCOH_new), loc = mu, scale = std) * 2
        reject, pval_corr = fdr_correction(pvals, alpha = 0.05, method='indep')
        [rows, cols]= np.where(reject == True)
        weight = np.zeros(ciCOH_new.shape)
        if len(rows) > 0:
            weight[rows, cols] = ciCOH_new[rows, cols]

        # visual and save
        saveFCGraph = os.path.join(savefolder, cond + '_' + save_prefix + '_' + datestr + '.png')
        texts = dict()
        texts[datestr] = [80, 50, 15]
        weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file, chnAreas_new), 
                            savefile = saveFCGraph, texts = None, threds_edge = None)
        del ciCOH_new, chnAreas_new, save_prefix, areas_used
        del saveFCGraph, weight



def main():
    files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
    files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))
    files_moderate = glob.glob(os.path.join(inputfolder, '*_moderate_*'))


    files_normal.sort()
    subArea_dailyfc_visual(files_normal)

    files_mild.sort()
    subArea_dailyfc_visual(files_mild)


    files_moderate.sort()
    subArea_dailyfc_visual(files_moderate)


if __name__ == '__main__':
    main()