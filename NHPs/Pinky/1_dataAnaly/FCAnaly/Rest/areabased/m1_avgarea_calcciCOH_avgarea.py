# # Functionality of this notebook: 
# 
# * calculate the ciCOH for normal and mild LFP data in rest

import os, sys
import scipy.io as sio
import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle 
import math
import pandas as pd


# extract the exp folder path
currfolder = os.getcwd()
codefolder = currfolder[0 : currfolder.find('code')+len('code')]

# add path the exp folder
sys.path.append(codefolder)

# import_nbmodule used for import package in .ipynb
import import_nbmodule

# import util/folder_extract.pynb 
from util.folder_extract import exp_subfolders, code_corresfolder

# import ciCoherence_overtime in connAnalyTool/synchronization_indices.ipynb
from connAnalyTool.synchronization_indices import ciCoherence_acrosstrials
from connAnalyTool.synchronization_indices import ciCoherence_overtime

# %% [markdown]
# ## exp subfolders & code_corresfolder

# %%
_, _, pipelinefolder, _= exp_subfolders()


# %%
get_ipython().run_cell_magic('javascript', '', 'IPython.notebook.kernel.execute(\'nb_name = "\' + IPython.notebook.notebook_name + \'"\')')


# %%
nb_name = nb_name[0: nb_name.find('.ipynb')]

# corresfolder
corresfolder, correparentfolder = code_corresfolder(os.getcwd(), nb_name)

# %% [markdown]
# ## global parameters

# %%
freq = [26, 28]

sameDates4SKTRest = True

# %% [markdown]
# ## input setup

# %%
inputfolder = os.path.join(pipelinefolder, 'NHP_Pinky', '0_dataPrep', 'Rest', 'm6_restData_averagedAcrossOneArea')


# %%
variablesinLoadfile = ['lfpsegs', 'fs', 'chnAreas']


# %%
chnInf_folder = os.path.join(pipelinefolder, 'NHP_Pinky', '1_dataAnaly', 'FCAnaly','Rest')
chnInf_file = os.path.join(chnInf_folder, 'chn_brainArea_simCoord_BrainArea.csv')


# %%
if sameDates4SKTRest:
    
    sameDatesInfFile = os.path.join(pipelinefolder, 'NHP_Pinky', '0_dataPrep', 'Pinky_sameDatesUsedforSTKRest.csv')

# %% [markdown]
# ## Save Setup

# %%
savefolder = corresfolder
savefilename =  'ciCOH_rest' + '_freq' + str(freq[0]) + '_' + str(freq[1])

if sameDates4SKTRest:
    savefilename = savefilename + '_samedays'

# %% [markdown]
# ## extract lfp

# %%
def lfpallfiles_extract(files):
    if 'lfpdata' in locals():
        del lfpdata
    
    for i, file in enumerate(files):
        
        ### load data
        matdat = sio.loadmat(file, variable_names = variablesinLoadfile, 
                             struct_as_record = False, squeeze_me = True) 
        
        
        ### extract the noused channels, only calculate once
        if i == 0:
            
            # chnAreas
            chnAreas = matdat['chnAreas'].tolist()
            
            # fs: sample rate
            fs = matdat['fs'] 
             
        

        ### dealing lfp data
        
        # lfp (np.ndarray): ntemporal * nchns * ntrials
        lfpdata_1file = matdat['lfpsegs']
        
        # concatenate to lfpdata for all files
        if 'lfpdata' not in locals():
            lfpdata = lfpdata_1file
        else:
            lfpdata = np.concatenate((lfpdata, lfpdata_1file), axis = 2)
          
    
    return lfpdata, chnAreas


# %%
files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))


if sameDates4SKTRest: ## extract only the file with the dates used in datestrings_rest
    
    df = pd.read_csv(sameDatesInfFile)
    
    datestrings_rest = df['datestrings_rest']
    
    
    
    files = files_normal
    files_used = []
    for file in files:

        # extract the date string for the file, format '_20170915_'
        idx = file.find('_tdt')
        datestring = '_' + file[idx-8:idx] + '_'

        # if the date in datestrings_rest
        if datestring in set(datestrings_rest):
            files_used.append(file)

    files_normal = files_used
    del files_used


    files = files_mild
    files_used = []
    for file in files:

        # extract the date string for the file, format '_20170915_'
        idx = file.find('_tdt')
        datestring = '_' + file[idx-8:idx] + '_'

        # if the date in datestrings_rest
        if datestring in set(datestrings_rest):
            files_used.append(file)

    files_mild = files_used
    del files_used


# %%
print(len(files_normal))
print(len(files_mild))


# %%
lfpwins_normal, chnAreas = lfpallfiles_extract(files_normal)
lfpwins_mild, _ = lfpallfiles_extract(files_mild)


# %%
print(lfpwins_normal.shape)
print(lfpwins_mild.shape)

# %% [markdown]
# ## Assign the xy coord of each chnArea

# %%
# load channel coord from chnInf_file
df = pd.read_csv(chnInf_file, header = 0)

# fill in the x,y coordinates of each area in chnAreas based on the values in df_chninf
coord_x, coord_y = np.zeros(shape = [len(chnAreas), ]), np.zeros(shape = [len(chnAreas), ])
for chnArea in chnAreas:
    
    mask_area = (chnArea == df['brainarea'])

    x, y = df['simulated_x'][mask_area].to_numpy(), df['simulated_y'][mask_area].to_numpy()

    coord_x[mask_area], coord_y[mask_area] = x, y
    
    del mask_area, x, y

df_chninf = pd.DataFrame(data = {'chnAreas': chnAreas, 'coord_x': coord_x, 'coord_y': coord_y})
    
del df, coord_x, coord_y

# %% [markdown]
# ## Calculate ciCOH
# %% [markdown]
# ### balance mild and normal trials

# %%
# select the smaller trial number
ntrials_normal, ntrials_mild = lfpwins_normal.shape[2], lfpwins_mild.shape[2]

ntrials = min([ntrials_normal, ntrials_mild])

# balance trials by randomly selecting ntrials
idx_ntrials = np.random.randint(ntrials_normal, size = ntrials)
lfpwins_normal = lfpwins_normal[:,:,idx_ntrials]

idx_ntrials = np.random.randint(ntrials_mild, size = ntrials)
lfpwins_mild = lfpwins_mild[:,:,idx_ntrials]

# %% [markdown]
# ###  normal ciCOH

# %%
lfpwins_allfiles = lfpwins_normal

### calculate ciCOH
ntempo, nchns, nwins = lfpwins_allfiles.shape
ciCOH_allWins = np.zeros((nchns, nchns, nwins))
for wini in range(nwins):
    
    if wini % 100 == 0:
        print("wini = " + str(wini) + "/" + str(nwins))
    
    for chni in range(nchns -1):
        signal1 = lfpwins_allfiles[:, chni, wini]
        
        for chnj in range(chni+1, nchns):
            signal2 = lfpwins_allfiles[:, chnj, wini]
            
            # ciCOH_allWins assignment
            ciCOH_allWins[chni, chnj, wini] = ciCoherence_overtime(signal1, signal2)
            
            # symmetrical
            ciCOH_allWins[chnj, chni, wini] = ciCOH_allWins[chni, chnj, wini]
            
            del signal2
        del signal1
        
print("ciCOH calculated!!")


# %%
ciCOH = np.mean(ciCOH_allWins, axis = 2)
ciCOH_normal = ciCOH
del ciCOH, ciCOH_allWins

# %% [markdown]
# ### mild ciCOH

# %%
lfpwins_allfiles = lfpwins_mild

### calculate ciCOH
ntempo, nchns, nwins = lfpwins_allfiles.shape
ciCOH_allWins = np.zeros((nchns, nchns, nwins))
for wini in range(nwins):
    
    if wini%100 ==0:
        print("wini = " + str(wini) + "/" + str(nwins))
        
    for chni in range(nchns -1):
        signal1 = lfpwins_allfiles[:, chni, wini]
        
        for chnj in range(chni+1, nchns):
            signal2 = lfpwins_allfiles[:, chnj, wini]
            
            # ciCOH_allWins assignment
            ciCOH_allWins[chni, chnj, wini] = ciCoherence_overtime(signal1, signal2)
            
            # symmetrical
            ciCOH_allWins[chnj, chni, wini] = ciCOH_allWins[chni, chnj, wini]
            
            del signal2
        del signal1
print("ciCOH calculated!!")


# %%
ciCOH = np.mean(ciCOH_allWins, axis = 2)
ciCOH_mild = ciCOH
del ciCOH, ciCOH_allWins

# %% [markdown]
# ## save ciCOH

# %%
ciCOH = dict()
ciCOH['normal'], ciCOH['mild']  = ciCOH_normal, ciCOH_mild, 


fc = dict()
fc['ciCOH'] = ciCOH
fc['chnInf'] = df_chninf
fc['ntrials'] = ntrials


# %%
with open(os.path.join(savefolder, savefilename + '.pickle'), 'wb') as fp:
    pickle.dump(fc, fp, protocol=pickle.HIGHEST_PROTOCOL)


# %%
os.path.join(savefolder, savefilename)


