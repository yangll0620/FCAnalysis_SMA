import os, sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import re
import shutil
import pickle 
import ruptures as rpt
import pandas as pd
import ruptures as rpt


codefolder = re.match('.*exp.*code', os.getcwd()).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder
from connAnalyTool.network_visual import assign_coord2chnArea
from connAnalyTool.fc_visual_time import weight_visual_save


animal =  re.search('NHPs/[a-zA-Z]*/', os.getcwd()).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
correparentfolder = os.getcwd().replace('code', 'pipeline')
codefilename = os.path.split(__file__.replace('.py', ''))[-1]
corresfolder = os.path.join(correparentfolder, codefilename)


savefolder = corresfolder
if os.path.isdir(savefolder):
    shutil.rmtree(savefolder)
os.mkdir(savefolder)

inputfolder = os.path.join(correparentfolder, 'm1_SKT_dycFC')
area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')

aligns = ['targetonset', 'reachonset']
if animal == 'Jo':
    conds = ['normal', 'mild', 'moderate']
elif animal == 'Kitty':
    conds = ['normal', 'moderate']

t_AOI = [-200, 200]


for subfreqfolder in os.listdir(inputfolder):

    inputsubfolder = os.path.join(inputfolder, subfreqfolder)
    
    for event_name in aligns:

        dynfc = pd.read_pickle(os.path.join(inputsubfolder, animal + '_' + event_name + '_dynfc.pickle'))
        savename_prefix = animal + '_changepoint_' + event_name



        chnAreas, fs = dynfc['chnAreas'], dynfc['fs']

        """ detection using dig_dynfc """
        ts_all = dynfc['trun_dynfc']['ts'] * 1000
        idxs_AOI = np.logical_and(ts_all>= t_AOI[0], ts_all<= t_AOI[1])
        ts = ts_all[idxs_AOI]
        cpointsDic = dict()
        for cond in conds:
            trun_dynfc_all = dynfc['trun_dynfc'][cond]
            trun_dynfc = trun_dynfc_all[:, :, idxs_AOI]

            dig_dynfc = trun_dynfc.copy()
            dig_dynfc[dig_dynfc > 0] = 1

            used_dynfc = dig_dynfc

            # algin  used_dynfc (nchns * nchns * ntemp) into signal (ntemp * (nchns * (nchns -1)/2))　
            nchns, _, ntemp = used_dynfc.shape
            signal = np.zeros((ntemp, round(nchns * (nchns -1) / 2)))
            areaPairs = []
            for ni in range(ntemp):
                i = 0
                for ci in range(nchns -1):
                    for cj in range(ci + 1, nchns):
                        signal[ni, i] = np.reshape(used_dynfc[ci, cj, ni], (1, 1))
                        i = i + 1

                        if ni == 0:
                            areaPairs.append(chnAreas[ci] +  ' <-> ' + chnAreas[cj])

            # detection　　
            algo = rpt.Pelt(model="rbf").fit(signal)
            cpoints = algo.predict(pen=15)
            #t_changepoint = np.array(cpoints) / fs * 1000 + ts[0]

            cpointsDic[cond] = cpoints 

            del algo, cpoints
            del nchns, ntemp, signal, areaPairs
        
        cpointSet = dict()
        cpointSet['fs'], cpointSet['chnAreas'] = fs, chnAreas
        cpointSet['t_AOI'] = t_AOI
        cpointSet['cpoints'] = cpointsDic


        savename_prefix = animal + '_' + subfreqfolder + '_' + event_name
        savefile_dynfc_cpoint = os.path.join(savefolder, savename_prefix + '_dynfc_cpoint' + '.pickle')
        with open(savefile_dynfc_cpoint, 'wb') as fp:
            pickle.dump(cpointSet, fp, protocol= pickle.HIGHEST_PROTOCOL)

        del cpointSet, savefile_dynfc_cpoint