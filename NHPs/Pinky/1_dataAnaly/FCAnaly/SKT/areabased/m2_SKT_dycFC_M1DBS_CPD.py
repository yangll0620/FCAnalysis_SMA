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


animal =  re.search('NHPs/[a-zA-Z]*/', os.getcwd()).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
correparentfolder = os.getcwd().replace('code', 'pipeline')
codefilename = os.path.split(__file__.replace('.py', ''))[-1]
corresfolder = os.path.join(correparentfolder, codefilename)


freq = [26, 28]
inputfolder = os.path.join(correparentfolder, 'm1_SKT_dycFC_M1DBS', 'freq' + str(freq[0]) + '_' + str(freq[1]))
savefolder = corresfolder
if os.path.isdir(savefolder):
    shutil.rmtree(savefolder)
os.mkdir(savefolder)


aligns = ['targetonset', 'reachonset']
conds = ['normal', 'mild']
t_AOI = [-300, 300]
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
        trun_dynfc[trun_dynfc > 0] = 1

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
            plt.plot(ts, sig / len(idxs) * (i + 1), label = areaPairs[idxs[i]])
            del sig

        # plot the change points
        for i in range(len(t_changepoint)):
            plt.plot([t_changepoint[i], t_changepoint[i]], [0, 1], '--')

        # save figure
        plt.savefig(os.path.join(savefolder, savename_prefix  + cond + '.png'))
        plt.clf()





        # rbf = np.zeros((ntemp-1, ))
        # for ni in range(ntemp - 1):
        #     rbf[ni] = np.exp(-np.sum(np.power(signal[ni + 1, :] - signal[ni, :], 2)))


        # plt.plot(ts[0:-1], rbf)
        # plt.title(animal  + ' : rbf in ' + cond)
        # plt.savefig(os.path.join(savefolder, savename_prefix + '_rbf_' + cond + '.png'))
        # plt.clf()



""" detection using cosdiff """
# cosdiff_normal = dynfc['cosdiff']['normal']

# cosdiff = cosdiff_normal
# ntemp = cosdiff.shape[0]
# signal = np.reshape(cosdiff, (ntemp, 1))

# # detection
# algo = rpt.Pelt(model="rbf").fit(signal)
# result = algo.predict(pen=10)

# print(result)

# # display
# rpt.display(signal, result)
# plt.show()

