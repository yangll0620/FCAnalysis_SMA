import os, sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import re

codefolder = re.match('.*exp.*code', os.getcwd()).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder
from connAnalyTool.fc_visual_time import weight_visual_save


from lfpextract import lfp_align2_reachonset
from fc_extract import dynciCOH_from_lfptrials
from fc_extract import pval_perm_dynciCOH_SKT, truncate_dynfc, digitized_dynfc
from network_metrics import fcnetwork_avgCC
from matrix_operations import cosSimilarity_SVDComps
from network_visual import assign_coord2chnArea, generate_video



animal =  re.search('NHPs/[a-zA-Z]*/', os.getcwd()).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
corresfolder, correparentfolder = code_corresfolder(__file__)


freq = [26, 28]
inputfolder = os.path.join(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'm2_STKData_narrowfiltered' + str(freq[0]) + '_' + str(freq[1]))
savefolder = corresfolder
area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')



files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))



tdur_trial = [-0.5, 0.6]
lfptrials_normal, chnAreas, fs = lfp_align2_reachonset(files_normal, tdur_trial = tdur_trial, tmin_reach = tdur_trial[1], tmax_reach = 1)
lfptrials_mild, _, _ = lfp_align2_reachonset(files_mild, tdur_trial = tdur_trial, tmin_reach = tdur_trial[1], tmax_reach = 1)


# balance trial numbers for normal and mild conditions
ntrials_normal = lfptrials_normal.shape[2]
ntrials_mild = lfptrials_mild.shape[2]
ntrials = min(ntrials_normal, ntrials_mild)
lfptrials_normal = lfptrials_normal[:, :, 0:ntrials]
lfptrials_mild = lfptrials_mild[:, :, 0:ntrials]


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
t = np.array([*range(0, ntemp, 1)]) / fs + tdur_trial[0]

prefix = 'dyn_avg_CC'
plt.plot(t, dyn_avgCC_normal, 'b', label = 'normal')
plt.plot(t, dyn_avgCC_mild, 'r', label = 'mild')
plt.xticks(ticks= [-0.4, -0.2, 0, 0.2, 0.4, 0.6], labels = ['-0.4', '-0.2', 'reach onset', '0.2', '0.4', '0.6'])
plt.legend()
plt.xlabel('time/s')
plt.ylabel('avg_CC')
plt.title(animal)

plt.savefig(os.path.join(savefolder, prefix + '.png'))
plt.clf()
print('save' +  prefix + ' at ' + os.path.join(savefolder))


# cos diff along time
cosdiff_normal = cosSimilarity_SVDComps(trun_dynfc_normal)
cosdiff_mild = cosSimilarity_SVDComps(trun_dynfc_mild)

# plot dyn_avg_CC
t = np.array([*range(0, ntemp-1, 1)]) / fs + tdur_trial[0]

prefix = 'cos diff'
plt.subplot(2, 1, 1)
plt.plot(t, cosdiff_normal, 'b', label = 'normal')
plt.legend()
plt.title(animal)
plt.xticks([])
plt.subplot(2, 1, 2)
plt.plot(t, cosdiff_mild, 'r', label = 'mild')
plt.xticks(ticks= [-0.4, -0.2, 0, 0.2, 0.4, 0.6], labels = ['-0.4', '-0.2', 'reach onset', '0.2', '0.4', '0.6'])
plt.legend()
plt.xlabel('time/s')
plt.ylabel('cos diff')


plt.savefig(os.path.join(savefolder, prefix + '.png'))
plt.clf()
print('save' +  prefix + ' at ' + os.path.join(savefolder))


# plot dyn_FC video
images = []
fps = fs
for ni in range(ntemp):
    saveFCGraph = os.path.join(savefolder, str(ni) + '.png')
    texts = dict()
    texts['t = ' + str((round(ni/fs + tdur_trial[0], 3)) * 1000) + 'ms'] = [-80, 40, 15]
    weight_visual_save(trun_dynfc_normal[:, :, ni], chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                            savefile = saveFCGraph, texts = texts, threds_edge = None)

    images.append(saveFCGraph)

generate_video(genvideofile = os.path.join(savefolder, animal + '_normal.mp4'), fps = fps, images = images)
for ni in range(ntemp):
    os.remove(os.path.join(savefolder, str(ni) + '.png'))



images = []
for ni in range(ntemp):
    saveFCGraph = os.path.join(savefolder, str(ni) + '.png')
    texts = dict()
    texts['t = ' + str((round(ni/fs + tdur_trial[0], 3)) * 1000) + 'ms'] = [-80, 40, 15]
    weight_visual_save(trun_dynfc_mild[:, :, ni], chnInf = assign_coord2chnArea(area_coord_file = area_coord_file, chnAreas = chnAreas), 
                            savefile = saveFCGraph, texts = texts, threds_edge = None)

    images.append(saveFCGraph)

generate_video(genvideofile = os.path.join(savefolder, animal + '_mild.mp4'), fps = fps, images = images)
for ni in range(ntemp):
    os.remove(os.path.join(savefolder, str(ni) + '.png'))
