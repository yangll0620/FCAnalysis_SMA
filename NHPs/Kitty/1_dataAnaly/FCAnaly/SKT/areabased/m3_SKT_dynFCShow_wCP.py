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
from connAnalyTool.fc_visual_time import weight_visual_save
from connAnalyTool.network_visual import assign_coord2chnArea, generate_video



animal =  re.search('NHPs/[a-zA-Z]*/', os.getcwd()).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
correparentfolder = os.getcwd().replace('code', 'pipeline')
codefilename = os.path.split(__file__.replace('.py', ''))[-1]
corresfolder = os.path.join(correparentfolder, codefilename)

freq = [21, 23]
inputfolder1 = os.path.join(correparentfolder, 'm1_SKT_dycFC')
inputfolder2 = os.path.join(correparentfolder, 'm2_SKT_dycFC_CPD')


savefolder = os.path.join(corresfolder)
if os.path.isdir(savefolder):
    shutil.rmtree(savefolder)
os.mkdir(savefolder)


area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')


aligns = ['targetonset', 'reachonset']
if animal == 'Jo':
    conds = ['normal', 'mild', 'moderate']
elif animal == 'Kitty':
    conds = ['normal', 'moderate']


for subfreqfolder in os.listdir(inputfolder1):

    inputsubfolder = os.path.join(inputfolder1, subfreqfolder)
    
    savesubfolder = os.path.join(savefolder, subfreqfolder)
    if os.path.isdir(savesubfolder):
        shutil.rmtree(savesubfolder)
    os.mkdir(savesubfolder)

    for event_name in aligns:

        savename_prefix = animal + '_' + event_name

        # load dynfc and cpointSet
        dynfc = pd.read_pickle(os.path.join(inputsubfolder, animal + '_' + event_name + '_dynfc.pickle'))
        cpointSet = pd.read_pickle(os.path.join(inputfolder2, animal + '_' + subfreqfolder + '_' + event_name + '_dynfc_cpoint.pickle'))


        t_AOI = cpointSet['t_AOI']
        
        fs = round(dynfc['fs'])
        chnAreas = dynfc['chnAreas']
        
        fps = round(fs / 10) # fps for generated video
        ts_all = np.round(dynfc['trun_dynfc']['ts'] * 1000, decimals=0)
        tmpfolder = os.path.join(savesubfolder, 'temp')
        
        # extract only the t_AOI
        idxs_tAOI = np.logical_and(ts_all >= t_AOI[0], ts_all <= t_AOI[1])
        ts = ts_all[idxs_tAOI]


        
        for cond in conds:
            trun_dynfc_all = dynfc['trun_dynfc'][cond]
            trun_dynfc = trun_dynfc_all[:,:, idxs_tAOI]
            
            cpoints = cpointSet['cpoints'][cond]
            t_changepoint = (np.array(cpoints) -1 ) / fs * 1000 + ts[0]

            dig_fc = trun_dynfc.copy()
            dig_fc[dig_fc > 0] = 1


            """ plot dyn_FC video in t_AOI with change point """ 
            os.mkdir(tmpfolder)
            images = []
            segi = 0
            for ni in range(ts.shape[0]):
                saveFCGraph = os.path.join(tmpfolder, str(ni) + '.png')
                texts = dict()
                texts['t = ' + str(round(ts[ni],3)) + 'ms'] = [-80, 40, 15]
                texts[animal + ': ' + cond + ',' + event_name] = [300, 40, 15]

                texts['Segi = ' + str(segi)] = [-80, 240, 15]

                # add change point marker if it is
                if ni + 1 in cpoints and (ni + 1) != ts.shape[0]:
                    segi = segi + 1

                weight_visual_save(trun_dynfc[:, :, ni], chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = dynfc['chnAreas']), 
                                        savefile = saveFCGraph, texts = texts, threds_edge = None)

                images.append(saveFCGraph)

            generate_video(genvideofile = os.path.join(savesubfolder, savename_prefix + '_wCP_' + cond + '.avi'), fps = fps, images = images)
            shutil.rmtree(tmpfolder)
            del images, segi
            """ end plot dyn_FC video """ 



            """ Summary FC during segments """
            for ci in range(len(cpoints)):

                if ci == 0 : 
                    if cpoints[ci] != 0: # [0 cpoints[0]]
                        idx_str, idx_end = 0, cpoints[ci]
                        trange = [ts[0], t_changepoint[ci]]
                    else:
                        continue

                else:
                    idx_str, idx_end = cpoints[ci-1], cpoints[ci]
                    trange = [t_changepoint[ci -1], t_changepoint[ci]]

                sum1FC = np.sum(dig_fc[:, :, idx_str:idx_end], axis = 2)
                sum01 = idx_end - idx_str
                summaryFC = sum1FC / sum01
                summaryFC[summaryFC < 0.8] = 0

                
                saveFCGraph = os.path.join(savesubfolder, savename_prefix + '_summaryFC_'  + cond + '_' + str(ci) +  '.png')
                texts = dict()
                texts['[' + str(trange[0]) + ' ' + str(trange[1]) + ' ]ms'] = [-80, 40, 15]
                texts[animal + ': ' + cond + ',' + event_name] = [300, 40, 15]
                weight_visual_save(summaryFC, chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                                    savefile = saveFCGraph, texts = texts, threds_edge = None)
            """ End summary FC """

            

            """ Show Change Points on dig_fc """
            areas = ['M1', 'stn', 'gp']
            colors = ['r', 'g','b','c','m','y']

            nsub = 5
            subi = 0

            ### M1-stn, M1-gp, stn-gp
            for ai in range(len(areas)):
                area = areas[ai]
                idxs1 = [i for i in range(len(chnAreas)) if area in chnAreas[i]]
                
                for aj in range(ai + 1, len(areas)):
                    area2 = areas[aj]

                    subi = subi + 1
                    idxs2 = [i for i in range(len(chnAreas)) if area2 in chnAreas[i]]
                    
                    plt.subplot(nsub, 1, subi)
                    
                    for idx1 in idxs1:
                        site1 = chnAreas[idx1]
                        for idx2 in idxs2:
                            site2 = chnAreas[idx2]
                            sig = dig_fc[idx1, idx2, :]

                            if np.any(sig) and np.any(sig - 1): # only show sig with change
                                plt.plot(ts, sig, label = site1 + '<>' + site2)

                            del sig
                    
                    # plot the change points
                    for i in range(len(t_changepoint)):
                        if t_changepoint[i] != ts[-1]:
                            plt.plot([t_changepoint[i], t_changepoint[i]], plt.ylim(), colors[i] + '--')

                    plt.xlim([ts[0], ts[-1]])
                    plt.ylabel(area + '-' + area2)

                    if subi < nsub:
                        ticks, labels = plt.xticks()
                        plt.xticks(ticks = ticks, labels = [])

                    if subi == 1:
                        plt.title(savename_prefix + ': ' + cond)
                    
                    del area2, idxs2
                
                del area, idxs1

            
            ### stn-stn, gp-gp
            for ai in range(len(areas)):
                area = areas[ai]
                if area == 'M1':
                    continue

                idxs = [i for i in range(len(chnAreas)) if area in chnAreas[i]]                                    
                subi = subi + 1
                
                plt.subplot(nsub, 1, subi)
                
                for idx1 in idxs:
                    site1 = chnAreas[idx1]
                    for idx2 in idxs:
                        if idx1 == idx2:
                            continue

                        site2 = chnAreas[idx2]
                        sig = dig_fc[idx1, idx2, :]

                        if np.any(sig) and np.any(sig - 1): # only show sig with change
                            plt.plot(ts, sig, label = site1 + '<>' + site2)

                        del sig
                
                # plot the change points
                for i in range(len(t_changepoint)):
                    if t_changepoint[i] != ts[-1]:
                        plt.plot([t_changepoint[i], t_changepoint[i]], plt.ylim(), colors[i] + '--')

                plt.xlim([ts[0], ts[-1]])
                plt.ylabel(area + '-' + area)

                if subi < nsub:
                    ticks, labels = plt.xticks()
                    plt.xticks(ticks = ticks, labels = [])

                if subi == 1:
                    plt.title(savename_prefix + ': ' + cond)

                
                del area, idxs

            # save figure
            plt.savefig(os.path.join(savesubfolder, savename_prefix + '_digFC_CP_' + cond + '.png'))
            plt.clf()
            
            
            
            del trun_dynfc_all, trun_dynfc, cpoints
            
        


        
        
        
        
        
        del savename_prefix, dynfc, cpointSet, t_AOI, fps, ts_all, tmpfolder, idxs_tAOI, ts
