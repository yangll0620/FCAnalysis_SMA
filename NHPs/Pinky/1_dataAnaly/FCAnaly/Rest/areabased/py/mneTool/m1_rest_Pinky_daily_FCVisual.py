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


def imcohs_daily_calc(onefile):

        
    lfpdata, chnAreas, fs = lfp_extract([onefile])

    if lfpdata.shape[2] < 80:
        return None


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
    filename = os.path.basename(onefile)
    datestr = re.search('[0-9]{8}', filename).group()
    cond = re.search('_[a-z]*_[0-9]{8}', filename).group()[1:-9]
    freqstr = 'freq' + str(freqs[0]) + '_' + str(freqs[1])
    fcfile_pickle =  os.path.join(savefolder, freqstr + '_' + cond + '_' + datestr + '.pickle')
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


    folder, filename = os.path.split(fcfile_pickle)[0], os.path.split(fcfile_pickle)[1]
    

    save_prefix = 'all'
    saveFCGraph = os.path.join(folder,  'visual_'  + filename[:-len('.pickle')] + '_' + save_prefix + '.png')
    texts = dict()
    texts[cond] = [-80, 40, 15]
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

def fcmetrics_alldates(files):
    

    nfiles = len(files)
    ccs = np.zeros(shape=(nfiles, ))
    with open(files[0], 'rb') as handle:
        metrics = pickle.load(handle)
    nbcs = np.zeros(shape = (nfiles, len(metrics['nbcs'])))    
    for fi, file in enumerate(files):

        with open(file, 'rb') as handle:
            metrics = pickle.load(handle)

        ccs[fi] = metrics['cc']

        for ki, key in enumerate(metrics['nbcs'].keys()):
            nbcs[fi, ki] = metrics['nbcs'][key]

    return ccs, nbcs


        
def main():
    # files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
    # files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))
    # files_moderate = glob.glob(os.path.join(inputfolder, '*_moderate_*'))


    # for onefile in files_normal:
    #     fcfile_pickle = imcohs_daily_calc(onefile)

    #     if fcfile_pickle is not None:
    #         fc_visual(fcfile_pickle)
    #         fc_metrics(fcfile_pickle)


    # for onefile in files_mild:
    #     fcfile_pickle = imcohs_daily_calc(onefile)

    #     if fcfile_pickle is not None:
    #         fc_visual(fcfile_pickle)
    #         fc_metrics(fcfile_pickle)


    # for onefile in files_moderate:
    #     fcfile_pickle = imcohs_daily_calc(onefile)

    #     if fcfile_pickle is not None:
    #         fc_visual(fcfile_pickle)
    #         fc_metrics(fcfile_pickle)

    fcfiles = glob.glob(os.path.join(savefolder, 'freq*.pickle'))
    for fcfile_pickle in fcfiles:
        fc_metrics_subareas(fcfile_pickle, subareas = ['M1', 'STN', 'GP'], subtitle = 'M1DBS')

    cond = 'normal'
    metricfiles_normal = glob.glob(os.path.join(savefolder, 'metric_M1DBS_*_' + cond + '_*.pickle'))
    ccs_normal, nbcs = fcmetrics_alldates(metricfiles_normal)


    cond = 'mild'
    metricfiles_normal = glob.glob(os.path.join(savefolder, 'metric_M1DBS_*_' + cond + '_*.pickle'))
    ccs_mild, nbcs = fcmetrics_alldates(metricfiles_normal)



    cond = 'moderate'
    metricfiles_normal = glob.glob(os.path.join(savefolder, 'metric_M1DBS_*_' + cond + '_*.pickle'))
    ccs_moderate, nbcs = fcmetrics_alldates(metricfiles_normal)


    plt.plot(np.zeros(shape = ccs_normal.shape), ccs_normal, 'r*', 
             np.zeros(shape = ccs_mild.shape) + 1, ccs_mild, 'bs',
             np.zeros(shape = ccs_moderate.shape) + 2, ccs_moderate, 'g^')

    plt.xticks(ticks = [0,1,2], labels = ['normal', 'mild', 'moderate'])
    plt.ylabel('daily Clustering Coefficient')
    plt.title(animal + ' M1DBS')

    plt.savefig(os.path.join(savefolder, 'visualcc_M1DBS.png'))


if __name__ == '__main__':
    main()