import os, sys
import scipy.io as sio
import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle 
import math
import pandas as pd
import cv2

import re
codefolder = re.match('.*exp.*code', __file__).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder

from lfpextract import lfp_align2_reachonset, lfp_align2_returnonset, lfp_align2_targetonset
from testfunc import assign_coord2chnArea, calcciCOH_from_lfptrials, fc4drawing
from visualFC import ciCOH_visual_save
from simulated.python.threshold_ciCOH import threshold_ciCOH_sin, corr_threshold_ciCOH_sin_BH

_, _, pipelinefolder, _= exp_subfolders()
corresfolder, correparentfolder = code_corresfolder(__file__)


def phaseFC_extract(phase):

    fcfile = os.path.join(savefolder, saveFCfilename_prefix + phase + '.pickle')
    if not os.path.exists(fcfile):

        files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
        files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))
        files_moderate = glob.glob(os.path.join(inputfolder, '*_moderate_*'))

        if phase == 'reach':
            tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate = [0, 0.5], [0, 0.7], [0, 0.7]

            lfptrials_normal, chnAreas, fs = lfp_align2_reachonset(files_normal, tdur_trial = tdur_trial_normal, tmin_reach = 0.5, tmax_reach = 1.5)
            lfptrials_mild, _, _ = lfp_align2_reachonset(files_mild, tdur_trial = tdur_trial_mild, tmin_reach = 0.5, tmax_reach = 1.5)
            lfptrials_moderate, _, _ = lfp_align2_reachonset(files_moderate, tdur_trial = tdur_trial_moderate, tmin_reach = 0.5, tmax_reach = 1.5)

        if phase == 'return':
            tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate = [0, 0.5], [0, 0.9], [0, 1]

            lfptrials_normal, chnAreas, fs = lfp_align2_returnonset(files_normal, tdur_trial = tdur_trial_normal, tmin_return = 0.5, tmax_return = 1.5)
            lfptrials_mild, _, _ = lfp_align2_returnonset(files_mild, tdur_trial = tdur_trial_mild, tmin_return = 0.5, tmax_return = 1.5)
            lfptrials_moderate, _, _ = lfp_align2_returnonset(files_moderate, tdur_trial = tdur_trial_moderate, tmin_return = 0.5, tmax_return = 1.5)


        if phase == 'base':
            tdur_trial = [-0.5, 0]

            lfptrials_normal, chnAreas, fs = lfp_align2_targetonset(files_normal, tdur_trial = tdur_trial, tmin_return = 0.5, tmax_return = 1)
            lfptrials_mild, _, _ = lfp_align2_targetonset(files_mild, tdur_trial = tdur_trial, tmin_return = 0.5, tmax_return = 1)
            lfptrials_moderate, _, _ = lfp_align2_targetonset(files_moderate, tdur_trial = tdur_trial, tmin_return = 0.5, tmax_return = 1)
            
        fc =  fc4drawing(chnAreas, fs, chnInf_file, lfptrials_normal = lfptrials_normal, lfptrials_mild = lfptrials_mild, lfptrials_moderate = lfptrials_moderate)

        with open(fcfile, 'wb') as fp:
            pickle.dump(fc, fp, protocol=pickle.HIGHEST_PROTOCOL)

    else:
        fc = pickle.load(open(fcfile, 'rb'))

    return fc


def threshold_fc(ciCOH, ntrials, ntemp, f, t):
    """
        identify threshold for fc
    """

    df = pd.DataFrame(columns = ['ntemp', 'ntrials', 'f','threshold', 'mu', 'std'])

    record_exist = False
    if(os.path.exists(savefile_threshold)): # file_threshold exists
        
        df = pd.read_pickle(savefile_threshold)
        
        # check if threshold under the ntemp and f exist 
        mask = (df['ntemp'] == ntemp) & (df['ntrials'] == ntrials) & (df['f'] == f)
           
        if(df.loc[mask].shape[0] == 1): # record exist
            
            record_exist = True
            
            df1 = df.loc[mask]

            threshold = df1['threshold'].to_numpy()
            mu, std = df1['mu'].to_numpy(), df1['std'].to_numpy()

            del df1

    if not record_exist: # run threshold_ciCOH_sin if record not exist
        threshold, mu, std = threshold_ciCOH_sin(ntimes = 500, ntrials = ntrials, ntemp = ntemp, 
                                                 f = f, t = t,ploton= False)


        # store the new record
        thred_sets = dict()
        thred_sets['ntemp'],  thred_sets['ntrials'], thred_sets['f'] = ntemp, ntrials, f
        thred_sets['mu'], thred_sets['std'] = mu, std
        thred_sets['threshold'] = threshold

        df = df.append(thred_sets, ignore_index = True)

        # write to file_threshold
        df.to_pickle(savefile_threshold)




    nchns = ciCOH.shape[0]
    ciCOHs_vec = list()
    for chni in range(nchns -1):
        for chnj in range(chni+1, nchns):
            ciCOHs_vec.append(ciCOH[chni][chnj])
    ciCOHs_vec = np.asarray(ciCOHs_vec)

    corr_threshold, mu, std = corr_threshold_ciCOH_sin_BH(ciCOHs_actual = ciCOHs_vec, ntimes = 500, 
                                ntrials = ntrials, ntemp = ntemp, f = f, t = t, false_rate = 0.01, mu = mu, std = std)
    

    threshold = np.around(threshold, decimals=2)
    corr_threshold = np.around(corr_threshold, decimals=2)
    return threshold, corr_threshold

def fc_visual_save(fc, savefile_prefix):

    ### text setup for brain areas ###
    pos_text_lefttop1 = [-80, 50, 30]
    pos_text_lefttop2 = [-80, 70, 10]
    pos_text_Down1 = [0, 500, 15]
    pos_text_Down2 = [0, 550, 15]
    pos_text_Down3 = [0, 570, 15]
    
    texts_org = dict()
    texts_org['STN'] = [200, 120, 20]
    texts_org['GP'] = [430, 210, 20]
    text_meta = animal + ": " + phase + ", [" + str(freq[0])  + " "+ str(freq[1]) + "] Hz"
    texts_org[text_meta] = pos_text_Down1

      

    for cond in fc['ciCOH'].keys():
        ciCOH = fc['ciCOH'][cond]
        ntrials, ntemp = fc['setup']['ntrials_' + cond], fc['setup']['ntemp_' + cond]
        t = ntemp/fc['setup']['fs']

        threshold, corr_threshold = threshold_fc(ciCOH = ciCOH, ntrials = ntrials, ntemp = ntemp, 
                                    f = (freq[0] + freq[1])//2, t = t)
        
        texts = texts_org.copy()
        lowweight = threshold
        texts[cond] = pos_text_lefttop1
        text_thred = 'thred = ' + str(lowweight)
        texts[text_thred] = pos_text_lefttop2
        text_ntrials = 'ntrials = ' + str(ntrials)
        texts[text_ntrials] = pos_text_Down2
        text_temp = 't = ' + str(t)
        texts[text_temp] = pos_text_Down3

        saveFCGraph = os.path.join(savefolder, savefile_prefix + '_lowweight' + str(lowweight) + '_'  + cond + '.png')

        igplot = ciCOH_visual_save(ciCOH = ciCOH, chnInf = fc['chnInf'], lowweight = lowweight, 
                                savefile = saveFCGraph, texts = texts, threds_edge = None)

        del texts[cond], texts[text_ntrials]

    return igplot

def combine_imgs():
    phases = ['base', 'reach', 'return']
    conds = ['normal', 'mild', 'moderate']

    # find all freqstr
    files = glob.glob(os.path.join(savefolder, 'ciCOH_*_reach_*.png'))
    freqstrs = list()
    for f in files:
        fname = os.path.split(f)[1] 
        
        freqstr = re.search('freq[0-9]*_[0-9]*', fname).group()
        freqstrs.append(freqstr)

    ufreqstr = set(freqstrs)

    print(ufreqstr)

    for freqstr in ufreqstr:

        for ci, cond in enumerate(conds):
            for pi, phase in enumerate(phases):
                files = glob.glob(os.path.join(savefolder, '*' + freqstr + '*_' + phase + '_*_' + cond + '.png'))
                file_fc = files[0]
                img = cv2.imread(file_fc)

                if pi == 0:
                    img_1cond = img
                    
                else:
                    img_1cond = np.concatenate((img_1cond, np.zeros((600, 2, 3)), img), axis = 1)

                del files, file_fc, img

            if ci == 0:
                imgs = img_1cond
                
            else:
                imgs = np.concatenate((imgs, np.zeros((2, img_1cond.shape[1], 3)), img_1cond), axis = 0)


        cv2.imwrite(os.path.join(savefolder, 'combined_' + freqstr +  '.png'), imgs)
        del imgs


animal =  re.search('NHPs/[a-zA-Z]*/', __file__).group()[len('NHPs/'):-1]
chnInf_folder = correparentfolder
chnInf_file = os.path.join(chnInf_folder, 'chn_brainArea_simCoord_BrainArea.csv')
savefolder = corresfolder

savefile_threshold = os.path.join(savefolder, 'threshold.pickle')

freq_opts = [[24, 26],[9, 11]]

for freq in freq_opts:
    inputfolder = os.path.join(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'm3_STKData_narrowfiltered' + str(freq[0]) + '_' + str(freq[1]))

    saveFCfilename_prefix = 'ciCOH_SKT_freq' + str(freq[0]) + '_' + str(freq[1])
    savefile_fcgraph_prefix = 'ciCOH_SKT_freq' + str(freq[0]) + '_' + str(freq[1])

    phase = 'base'
    fc = phaseFC_extract(phase)
    fc_visual_save(fc = fc, savefile_prefix = savefile_fcgraph_prefix + '_' + phase)

    phase = 'reach'
    fc = phaseFC_extract(phase)
    fc_visual_save(fc = fc, savefile_prefix = savefile_fcgraph_prefix + '_' + phase)

    phase = 'return'
    fc = phaseFC_extract(phase)
    fc_visual_save(fc = fc, savefile_prefix = savefile_fcgraph_prefix + '_' + phase)


combine_imgs()





