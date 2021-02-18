import os, sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import re
import shutil
import pickle 

codefolder = re.match('.*exp.*code', os.getcwd()).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder
from connAnalyTool.fc_visual_time import weight_visual_save


from connAnalyTool.lfpextract import lfp_align2, Event 
from connAnalyTool.fc_extract import dynciCOH_from_lfptrials
from connAnalyTool.fc_extract import pval_perm_dynciCOH_SKT, truncate_dynfc, digitized_dynfc
from connAnalyTool.network_metrics import fcnetwork_avgCC
from connAnalyTool.matrix_operations import cosSimilarity_SVDComps
from connAnalyTool.network_visual import assign_coord2chnArea, generate_video



animal =  re.search('NHPs/[a-zA-Z]*/', os.getcwd()).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
corresfolder, correparentfolder = code_corresfolder(__file__)


freq = [26, 28]
inputfolder = os.path.join(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'm2_STKData_narrowfiltered' + str(freq[0]) + '_' + str(freq[1]))
savefolder = os.path.join(corresfolder, 'freq' + str(freq[0]) + '_' + str(freq[1]))
if os.path.isdir(savefolder):
    shutil.rmtree(savefolder)
os.mkdir(savefolder)

area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')



files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))

areas_used = ['M1', 'stn', 'gp']


t_minmax_reach_normal = [0.6, 1]
t_minmax_reach_mild = [0.6, 1]
t_minmax_return_normal = [0.5, 1]
t_minmax_return_mild = [0.6, 1.3]


tdur_targetonset = [-0.5, 0.6]
tdur_reachonset = [-0.5, 0.6]
tdur_returnonset = [-0.2, 0.6]


aligns = dict()
event = dict()
event['tdur_trial'] = tdur_targetonset
event['event_name'] = 'targetonset'
aligns[Event.TARGETONSET] = event
del event

event = dict()
event['tdur_trial'] = tdur_reachonset
event['event_name'] = 'reachonset'
aligns[Event.REACHONSET] = event
del event

event = dict()
event['tdur_trial'] = tdur_returnonset
event['event_name'] = 'returnonset'
aligns[Event.RETURNONSET] = event
del event



for align2 in aligns.keys():
    tdur_trial = aligns[align2]['tdur_trial']
    event_name = aligns[align2]['event_name']
    

    savename_prefix = animal + '_' + event_name
    

    lfptrials_normal, chnAreas, fs = lfp_align2(files_normal, align2 = align2, tdur_trial = tdur_trial, t_minmax_reach = t_minmax_reach_normal, t_minmax_return = t_minmax_return_normal)
    lfptrials_mild, _, _ = lfp_align2(files_mild, align2 = align2, tdur_trial = tdur_trial, t_minmax_reach = t_minmax_reach_mild, t_minmax_return = t_minmax_return_mild)


    ###  generate the used chnAreas and lfptrials ###
    idx_del = []
    for i, area in enumerate(chnAreas):
        if 'stn' in area or 'gp' in area:
            continue
        if area not in areas_used:
            idx_del.append(i)
    # generate the used chnAreas
    idx_del.reverse()
    for i in idx_del:
        del chnAreas[i]
    # generate used lfptrials
    lfptrials_normal = np.delete(lfptrials_normal, idx_del, axis = 0)
    lfptrials_mild = np.delete(lfptrials_mild, idx_del, axis = 0)
    del idx_del

     

    # balance trial numbers for normal and mild conditions
    ntrials_normal = lfptrials_normal.shape[2]
    ntrials_mild = lfptrials_mild.shape[2]
    ntrials = min(ntrials_normal, ntrials_mild)
    lfptrials_normal = lfptrials_normal[:, :, 0:ntrials]
    lfptrials_mild = lfptrials_mild[:, :, 0:ntrials]
    print('ntrials_normal: '+ str(ntrials_normal) + ', ntrials_mild: '+ str(ntrials_mild))



    # calculate the dynamic ciCOH
    dynciCOH_normal = dynciCOH_from_lfptrials(lfptrials_normal)
    dynciCOH_mild = dynciCOH_from_lfptrials(lfptrials_mild)



    # digitized dynamic fc based on permutation test
    pvals_normal = pval_perm_dynciCOH_SKT(dynciCOH = dynciCOH_normal, lfptrials = lfptrials_normal)
    trun_dynfc_normal = truncate_dynfc(dynciCOH = dynciCOH_normal, pvals = pvals_normal)
    digi_dynfc_normal = digitized_dynfc(dynciCOH = dynciCOH_normal, pvals = pvals_normal)

    pvals_mild = pval_perm_dynciCOH_SKT(dynciCOH = dynciCOH_mild, lfptrials = lfptrials_mild)
    trun_dynfc_mild = truncate_dynfc(dynciCOH = dynciCOH_mild, pvals = pvals_mild)
    digi_dynfc_mild = digitized_dynfc(dynciCOH = dynciCOH_mild, pvals = pvals_mild)



    # calc  dynamic avg CC
    dynfc = digi_dynfc_normal
    ntemp = dynfc.shape[2]
    dyn_avgCC = np.zeros(shape = (ntemp, ))
    for ti in range(ntemp):
        dyn_avgCC[ti] = fcnetwork_avgCC(dynfc[:, :, ti])
    dyn_avgCC_normal = dyn_avgCC
    del dyn_avgCC, dynfc

    dynfc = digi_dynfc_mild
    ntemp = dynfc.shape[2]
    dyn_avgCC = np.zeros(shape = (ntemp, ))
    for ti in range(ntemp):
        dyn_avgCC[ti] = fcnetwork_avgCC(dynfc[:, :, ti])
    dyn_avgCC_mild = dyn_avgCC
    del dyn_avgCC, dynfc


    # plot dyn_avg_CC
    t_dynavgcc = np.array([*range(0, ntemp, 1)]) / fs + tdur_trial[0]


    plt.plot(t_dynavgcc, dyn_avgCC_normal, 'b', label = 'normal')
    plt.plot(t_dynavgcc, dyn_avgCC_mild, 'r', label = 'mild')
    plt.xticks(ticks= [-0.4, -0.2, 0, 0.2, 0.4, 0.6], labels = ['-0.4', '-0.2', event_name, '0.2', '0.4', '0.6'])
    plt.legend()
    plt.xlabel('time/s')
    plt.ylabel('avg_CC')
    plt.title(animal)

    plt.savefig(os.path.join(savefolder, savename_prefix + '_dyn_avg_CC' + '.png'))
    plt.clf()
    print('save' +  savename_prefix + '_dyn_avg_CC' + ' at ' + os.path.join(savefolder))


    # cos diff along time
    cosdiff_normal = cosSimilarity_SVDComps(trun_dynfc_normal)
    cosdiff_mild = cosSimilarity_SVDComps(trun_dynfc_mild)


    # plot dyn_avg_CC
    t_cosdiff = np.array([*range(0, ntemp-1, 1)]) / fs + tdur_trial[0]


    plt.subplot(2, 1, 1)
    plt.plot(t_cosdiff, cosdiff_normal, 'b', label = 'normal')
    plt.legend()
    plt.title(animal)
    plt.xticks([])
    plt.ylabel('cos diff')
    plt.subplot(2, 1, 2)
    plt.plot(t_cosdiff, cosdiff_mild, 'r', label = 'mild')
    plt.xticks(ticks= [-0.4, -0.2, 0, 0.2, 0.4, 0.6], labels = ['-0.4', '-0.2', event_name, '0.2', '0.4', '0.6'])
    plt.legend()
    plt.xticks([])
    plt.ylabel('cos diff')


    plt.savefig(os.path.join(savefolder, savename_prefix + '_cos diff' + '.png'))
    plt.clf()
    print('save' +  savename_prefix + '_cos diff' + '.png' + ' at ' + os.path.join(savefolder))


    # plot dyn_FC video
    fps = round(fs / 10)
    tmpfolder = os.path.join(savefolder, 'temp')

    trun_dynfc = trun_dynfc_normal
    cond = 'normal'
    os.mkdir(tmpfolder)
    images = []
    for ni in range(ntemp):
        saveFCGraph = os.path.join(tmpfolder, str(ni) + '.png')
        texts = dict()
        texts['t = ' + str((round(ni/fs + tdur_trial[0], 3)) * 1000) + 'ms'] = [-80, 40, 15]
        texts[animal + ': ' + cond + ',' + event_name] = [300, 40, 15]
        weight_visual_save(trun_dynfc[:, :, ni], chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                                savefile = saveFCGraph, texts = texts, threds_edge = None)

        images.append(saveFCGraph)

    generate_video(genvideofile = os.path.join(savefolder, savename_prefix + '_' + cond + '.avi'), fps = fps, images = images)
    shutil.rmtree(tmpfolder)
    del trun_dynfc, cond



    trun_dynfc = trun_dynfc_mild
    cond = 'mild'
    os.mkdir(tmpfolder)
    images = []
    for ni in range(ntemp):
        saveFCGraph = os.path.join(tmpfolder, str(ni) + '.png')
        texts = dict()
        texts['t = ' + str((round(ni/fs + tdur_trial[0], 3)) * 1000) + 'ms'] = [-80, 40, 15]
        texts[animal + ': ' + cond + ',' + event_name] = [300, 40, 15]
        weight_visual_save(trun_dynfc[:, :, ni], chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                                savefile = saveFCGraph, texts = texts, threds_edge = None)

        images.append(saveFCGraph)

    generate_video(genvideofile = os.path.join(savefolder,  savename_prefix + '_' + cond + '.avi'), fps = fps, images = images)
    shutil.rmtree(tmpfolder)
    del trun_dynfc, cond


    
    t_dynfc = np.array([*range(0, ntemp, 1)]) / fs + tdur_trial[0]


    # save all the data into .pickle
    trun_dynfc = dict()
    trun_dynfc['normal'] = trun_dynfc_normal
    trun_dynfc['mild'] = trun_dynfc_mild
    trun_dynfc['ts'] = t_dynfc


    cosdiff = dict()
    cosdiff['normal'] = cosdiff_normal
    cosdiff['mild'] = cosdiff_mild
    cosdiff['ts'] = t_cosdiff


    dynfc = dict()
    dynfc['trun_dynfc'] = trun_dynfc
    dynfc['cosdiff'] = cosdiff
    dynfc['chnAreas'] = chnAreas
    dynfc['fs'] = fs 

    savefile_dynfc = os.path.join(savefolder, savename_prefix + '_dynfc' + '.pickle')
    with open(savefile_dynfc, 'wb') as fp:
        pickle.dump(dynfc, fp, protocol=pickle.HIGHEST_PROTOCOL)
    del trun_dynfc, cosdiff, dynfc, savefile_dynfc

    del tdur_trial, event_name, savename_prefix
