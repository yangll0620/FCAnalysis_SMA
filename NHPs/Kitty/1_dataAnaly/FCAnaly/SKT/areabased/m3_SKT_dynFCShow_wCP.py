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

        
        ### plot dyn_FC video in t_AOI with change point ###
        fps = round(dynfc['fs'] / 10)
        ts_all = dynfc['trun_dynfc']['ts'] * 1000
        tmpfolder = os.path.join(savesubfolder, 'temp')
        
        # extract only the t_AOI
        idxs_tAOI = np.logical_and(ts_all >= t_AOI[0], ts_all <= t_AOI[1])
        ts = ts_all[idxs_tAOI]

        
        for cond in conds:
            trun_dynfc_all = dynfc['trun_dynfc'][cond]
            trun_dynfc = trun_dynfc_all[:,:, idxs_tAOI]
            cpoints = cpointSet['cpoints'][cond]
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
            
            del trun_dynfc_all, trun_dynfc, cpoints, images, segi

        
        del savename_prefix, dynfc, cpointSet, t_AOI, fps, ts_all, tmpfolder, idxs_tAOI, ts