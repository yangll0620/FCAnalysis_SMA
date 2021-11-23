import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import re
import glob
import numpy as np
from scipy.stats import norm
from mne.stats import fdr_correction
from mne.connectivity import spectral_connectivity
import pickle
import networkx as nx


codefolder = re.match('.*exp.*code', os.getcwd()).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder
from connAnalyTool.fc_visual_time import weight_visual_save


from lfp_ciCOH_extract import lfp_extract, calc_ciCOHs_rest
from statistical import pval_perm_imcohMNE_rest
from chnArea_set import assign_coord2chnArea


freqs = [26, 28]

animal =  re.search('NHPs/[a-zA-Z]*/', os.getcwd()).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
corresfolder, correparentfolder = code_corresfolder(__file__)



inputfolder = os.path.join(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'Rest', 'm4_restData_eqLen')


savefolder = corresfolder
area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')


def imcohs_combdays_calc(files):

    filename = os.path.basename(files[0])
    cond = re.search('_[a-z]*_[0-9]{8}', filename).group()[1:-9]
    freqstr = 'freq' + str(freqs[0]) + '_' + str(freqs[1])
    fcfile_pickle =  os.path.join(savefolder, freqstr + '_' + cond  + '.pickle')

    if os.path.exists(fcfile_pickle):
        return fcfile_pickle
        
    
    lfpdata, chnAreas, fs = lfp_extract(files)


    [imcohs, _, _, _, _]= spectral_connectivity(data = np.transpose(lfpdata, axes = (2, 0, 1)), 
                                            method= 'imcoh', sfreq = fs,fmin= 26, fmax = 28, faverage = True)
    imcohs = np.squeeze(imcohs)
    imcohs = imcohs + np.transpose(imcohs, axes = (1, 0))

    # permutation test: use the lfp data whose ciCOHs are the largest to get  distribution
    [i, j] = np.unravel_index(np.argmax(imcohs), shape = imcohs.shape)
    lfp1, lfp2 = lfpdata[i, :, :], lfpdata[j, :, :]
    mu, std = pval_perm_imcohMNE_rest(lfp1, lfp2, fs = fs, fmin = freqs[0], fmax = freqs[1], shuffleN = 300)
    pvals = norm.sf(abs(imcohs), loc = mu, scale = std) * 2
    del lfp1, lfp2


    fc = dict()
    fc['imcohs'] = imcohs
    fc['pvals'] = pvals
    fc['chnAreas'] = chnAreas

    # save
    with open(fcfile_pickle, 'wb') as f:
        pickle.dump(fc, f)

    return fcfile_pickle


def fc_visual(fcfile_pickle):

    with open(fcfile_pickle, 'rb') as handle:
        fc = pickle.load(handle)


    imcohs = fc['imcohs']
    pvals = fc['pvals']
    chnAreas = fc['chnAreas']


    # multiple comparison correction, get weights
    reject, pval_corr = fdr_correction(pvals, alpha = 0.05, method='indep')
    [rows, cols]= np.where(reject == True)
    weight = np.zeros(imcohs.shape)
    if len(rows) > 0:
        weight[rows, cols] = imcohs[rows, cols]


    for co in ['normal', 'mild', 'moderate']:
        if co in fcfile_pickle:
            cond = co 

    save_prefix = 'all'
    folder, filename = os.path.split(fcfile_pickle)[0], os.path.split(fcfile_pickle)[1]
    saveFCGraph = os.path.join(folder,  'visual_'  + filename[:-len('.pickle')] + '_' + save_prefix + '.png')
    texts = dict()
    texts[cond] = [-80, 40, 15]
    texts[animal] = [80, 20, 20]
    weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                        savefile = saveFCGraph, texts = texts, threds_edge = None)



def fc_visual_subAreas(fcfile_pickle, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS'):

    with open(fcfile_pickle, 'rb') as handle:
        fc = pickle.load(handle)


    imcohs = fc['imcohs']
    pvals = fc['pvals']
    chnAreas = fc['chnAreas']


    idxs_remain = []
    chnAreas_new = []
    for ci, carea in enumerate(chnAreas):
        for sarea in subareas:
            if sarea.lower() in carea.lower():
                idxs_remain.append(ci)
                chnAreas_new.append(carea)

    idxs_remain = np.array(idxs_remain)

    tmp = imcohs[idxs_remain, :]
    tmp = tmp[:, idxs_remain]
    imcohs = tmp
    tmp = pvals[idxs_remain, :]
    tmp = tmp[:, idxs_remain]
    pvals = tmp

    chnAreas = chnAreas_new


    # multiple comparison correction, get weights
    reject, pval_corr = fdr_correction(pvals, alpha = 0.05, method='indep')
    [rows, cols]= np.where(reject == True)
    weight = np.zeros(imcohs.shape)
    if len(rows) > 0:
        weight[rows, cols] = imcohs[rows, cols]


    for co in ['normal', 'mild', 'moderate']:
        if co in fcfile_pickle:
            cond = co 


    folder, filename = os.path.split(fcfile_pickle)[0], os.path.split(fcfile_pickle)[1]
    saveFCGraph = os.path.join(folder,  'visual_'  + filename[:-len('.pickle')] + '_' + subtitle + '.png')
    texts = dict()
    texts[cond] = [-80, 40, 15]
    texts[animal] = [80, 20, 20]
    weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                        savefile = saveFCGraph, texts = texts, threds_edge = None)

def fc_metrics(fcfile_pickle):
    """
        cc: average Clustering Coefficient

        nbc: Node Betweenness centrality ()
    """

    with open(fcfile_pickle, 'rb') as handle:
        fc = pickle.load(handle)


    imcohs = fc['imcohs']
    pvals = fc['pvals']


    # multiple comparison correction, get weights
    reject, pval_corr = fdr_correction(pvals, alpha = 0.05, method='indep')
    [rows, cols]= np.where(reject == True)
    weight = np.zeros(imcohs.shape)
    if len(rows) > 0:
        weight[rows, cols] = imcohs[rows, cols]



    weight = abs(weight)

    G = nx.Graph()
    G.add_nodes_from(np.arange(0, weight.shape[0]))

    for i in range(0, weight.shape[0] - 1):
        for j in range(i+1, weight.shape[0]):
            if weight[i,j] > 0:
                G.add_edge(i, j, weight= weight[i, j])


    cc = nx.average_clustering(G)
    nbcs = nx.degree_centrality(G)


    folder, filename = os.path.split(fcfile_pickle)[0], os.path.split(fcfile_pickle)[1]
    metricfile = os.path.join(folder,  'metric_'  + filename)

    metrics = dict()
    metrics['cc'] = cc
    metrics['nbcs'] = nbcs
    metrics['chnAreas'] =  fc['chnAreas']

    with open(metricfile, 'wb') as f:
        pickle.dump(metrics, f)
        

    return metricfile

def fc_metrics_subareas(fcfile_pickle, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS'):
    """
        cc: average Clustering Coefficient

        nbc: Node Betweenness centrality ()
    """

    with open(fcfile_pickle, 'rb') as handle:
        fc = pickle.load(handle)


    imcohs = fc['imcohs']
    pvals = fc['pvals']
    chnAreas = fc['chnAreas']

    idxs_remain = []
    chnAreas_new = []
    for ci, carea in enumerate(chnAreas):
        for sarea in subareas:
            if sarea.lower() in carea.lower():
                idxs_remain.append(ci)
                chnAreas_new.append(carea)

    idxs_remain = np.array(idxs_remain)

    tmp = imcohs[idxs_remain, :]
    tmp = tmp[:, idxs_remain]
    imcohs = tmp

    tmp = pvals[idxs_remain, :]
    tmp = tmp[:, idxs_remain]
    pvals = tmp


    chnAreas = chnAreas_new
    
    # multiple comparison correction, get weights
    reject, pval_corr = fdr_correction(pvals, alpha = 0.05, method='indep')
    [rows, cols]= np.where(reject == True)
    weight = np.zeros(imcohs.shape)
    if len(rows) > 0:
        weight[rows, cols] = imcohs[rows, cols]



    weight = abs(weight)

    G = nx.Graph()
    G.add_nodes_from(np.arange(0, weight.shape[0]))

    for i in range(0, weight.shape[0] - 1):
        for j in range(i+1, weight.shape[0]):
            if weight[i,j] > 0:
                G.add_edge(i, j, weight= weight[i, j])


    cc = nx.average_clustering(G)
    nbcs = nx.degree_centrality(G)


    folder, filename = os.path.split(fcfile_pickle)[0], os.path.split(fcfile_pickle)[1]
    metricfile = os.path.join(folder,  'metric_' + subtitle + '_'+ filename)

    metrics = dict()
    metrics['cc'] = cc
    metrics['nbcs'] = nbcs
    metrics['chnAreas'] =  fc['chnAreas']

    with open(metricfile, 'wb') as f:
        pickle.dump(metrics, f)

    return metricfile

def main():
    files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
    files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))
    files_moderate = glob.glob(os.path.join(inputfolder, '*_moderate_*'))

    # calculate imcohs, use permutation test (run long)
    fcfile_pickle_normal = imcohs_combdays_calc(files_normal)
    fcfile_pickle_mild = imcohs_combdays_calc(files_mild)
    fcfile_pickle_moderate =imcohs_combdays_calc(files_moderate)


    # visual all the fc
    fc_visual(fcfile_pickle_normal)
    fc_visual(fcfile_pickle_mild)
    fc_visual(fcfile_pickle_moderate)

    # calculate fc metrics
    metricfile_normal= fc_metrics_subareas(fcfile_pickle_normal, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS')
    metricfile_mild = fc_metrics_subareas(fcfile_pickle_mild, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS')
    metricfile_moderate = fc_metrics_subareas(fcfile_pickle_moderate, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS')


    with open(metricfile_normal, 'rb') as handle:
        metric = pickle.load(handle)
        cc_normal = metric['cc']
        nbcs_normal = metric['nbcs']
        chnAreas = metric['chnAreas']


    with open(metricfile_mild, 'rb') as handle:
        metric = pickle.load(handle)
        cc_mild = metric['cc']
        nbcs_mild = metric['nbcs']
        


    with open(metricfile_moderate, 'rb') as handle:
        metric = pickle.load(handle)
        cc_moderate = metric['cc']
        nbcs_moderate = metric['nbcs']

    print([cc_normal, cc_mild, cc_moderate])
    print([nbcs_normal, nbcs_mild, nbcs_moderate])

    
    # # visual subAreas fc
    # fc_visual_subAreas(fcfile_pickle_normal, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS')
    # fc_visual_subAreas(fcfile_pickle_mild, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS')
    # fc_visual_subAreas(fcfile_pickle_moderate, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS')
    
if __name__ == '__main__':
    main()