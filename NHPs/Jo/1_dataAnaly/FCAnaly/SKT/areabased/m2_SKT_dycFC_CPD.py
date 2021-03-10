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


freq = [24, 26]
inputfolder = os.path.join(correparentfolder, 'm1_SKT_dycFC', 'freq' + str(freq[0]) + '_' + str(freq[1]))
savefolder = corresfolder

area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')


aligns = ['targetonset', 'reachonset']
conds = ['normal', 'mild', 'moderate']
t_AOI = [-200, 200]
for event_name in aligns:


    dynfc = pd.read_pickle(os.path.join(inputfolder, animal + '_' + event_name + '_dynfc.pickle'))
    savename_prefix = animal + '_changepoint_' + event_name



    chnAreas, fs = dynfc['chnAreas'], dynfc['fs']

    """ detection using trun_dync """
    ts_all = dynfc['trun_dynfc']['ts'] * 1000
    idxs_AOI = np.logical_and(ts_all>= t_AOI[0], ts_all<= t_AOI[1])
    ts = ts_all[idxs_AOI]
    for cond in conds:
        trun_dynfc_all = dynfc['trun_dynfc'][cond]
        trun_dynfc = trun_dynfc_all[:, :, idxs_AOI]
        #trun_dynfc[trun_dynfc > 0] = 1

        # algin  trun_dynfc (nchns * nchns * ntemp) into signal (ntemp * (nchns * (nchns -1)/2))　
        nchns, _, ntemp = trun_dynfc.shape
        signal = np.zeros((ntemp, round(nchns * (nchns -1) / 2)))
        areaPairs = []
        for ni in range(ntemp):
            i = 0
            for ci in range(nchns -1):
                for cj in range(ci + 1, nchns):
                    signal[ni, i] = np.reshape(trun_dynfc[ci, cj, ni], (1, 1))
                    i = i + 1

                    if ni == 0:
                        areaPairs.append(chnAreas[ci] +  ' <-> ' + chnAreas[cj])

        # detection　　
        algo = rpt.Pelt(model="rbf").fit(signal)
        cpoints = algo.predict(pen=15)
        t_changepoint = np.array(cpoints) / fs * 1000 + ts[0]



        # only show connection with M1
        plt.subplot(2,1, 1)
        idxs = [i for i in range(len(areaPairs)) if 'M1' in areaPairs[i] and np.any(signal[:, i])]    
        for i in range(len(idxs)):
            sig = signal[:, idxs[i]]
            plt.plot(ts, sig, label = areaPairs[idxs[i]])
            del sig

        # plot the change points
        for i in range(len(t_changepoint)):
            plt.plot([t_changepoint[i], t_changepoint[i]], plt.ylim(), '--')

        plt.title(event_name + ': ' + cond)
        plt.ylabel('M1-STN/GP') 
   

        plt.subplot(2,1, 2)
        idxs = [i for i in range(len(areaPairs)) if 'M1' not in areaPairs[i] and np.any(signal[:, i])]    
        for i in range(len(idxs)):
            sig = signal[:, idxs[i]]
            plt.plot(ts, sig / len(idxs) * (i + 1), label = areaPairs[idxs[i]])
            del sig

        # plot the change points
        for i in range(len(t_changepoint)):
            plt.plot([t_changepoint[i], t_changepoint[i]], plt.ylim(), '--')
        plt.ylabel('STN-GP')
        plt.xlabel('t/ms')

        # save figure
        plt.savefig(os.path.join(savefolder, savename_prefix  + cond + '.png'))
        plt.clf()



        """ Summary FC during segments """
        dig_fc = trun_dynfc.copy()
        dig_fc[dig_fc > 0] = 1
        for ci in range(len(cpoints) - 1):
            idx_str, idx_end  = cpoints[ci], cpoints[ci+1]
            trange = [t_changepoint[ci], t_changepoint[ci + 1]]
            
            sum1FC = np.sum(dig_fc[:, :, idx_str:idx_end], axis = 2)
            sum01 = idx_end - idx_str
            summaryFC = sum1FC / sum01
            summaryFC[summaryFC < 0.8] = 0

            
            saveFCGraph = os.path.join(savefolder, animal + '_summaryFC_' + event_name + cond + '_' + str(trange[0]) + '_' + str(trange[1]) + 'ms.png')
            texts = dict()
            texts['[' + str(trange[0]) + ' ' + str(trange[1]) + ' ]ms'] = [-80, 40, 15]
            texts[animal + ': ' + cond + ',' + event_name] = [300, 40, 15]
            weight_visual_save(summaryFC, chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                                savefile = saveFCGraph, texts = texts, threds_edge = None)
        
        ci = np.where(t_changepoint < 0)[0][-1]
        idx_str, idx_end  = cpoints[ci],  np.where(ts == 0)[0][0]
        trange = [t_changepoint[ci], 0]
        sum1FC = np.sum(dig_fc[:, :, idx_str:idx_end], axis = 2)
        sum01 = idx_end - idx_str
        summaryFC = sum1FC / sum01
        summaryFC[summaryFC < 0.8] = 0

        
        saveFCGraph = os.path.join(savefolder, animal + '_premov_' + event_name + cond + '_' + str(trange[0]) + '_' + str(trange[1]) + 'ms.png')
        texts = dict()
        texts['[' + str(trange[0]) + ' ' + str(trange[1]) + ' ]ms'] = [-80, 40, 15]
        texts[animal + ': ' + cond + ',' + event_name] = [300, 40, 15]
        weight_visual_save(summaryFC, chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                            savefile = saveFCGraph, texts = texts, threds_edge = None)


