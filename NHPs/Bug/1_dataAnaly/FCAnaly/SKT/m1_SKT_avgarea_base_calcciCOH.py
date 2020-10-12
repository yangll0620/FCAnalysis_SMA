import os, sys
import numpy as np
import scipy.io as sio

# add exp code folder to path
currfolder = os.getcwd()
codefolder = currfolder[0: currfolder.find('code') + len('code')]
sys.path.append(codefolder)

from util.folder_extract import code_corresfolder
from connAnalyTool.synchronization_indices import ciCoherence_acrosstrials

folder_code_corres, folder_codecorresparent = code_corresfolder(__file__)

savefolder = code_corres_folder


def lfpallfiles_extract(files):
    if 'lfpdatas' in locals():
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
        lfpdata_1file = matdat['lfpdata']
        
        # concatenate to lfpdata for all files
        if 'lfpdatas' not in locals():
            lfpdatas = lfpdata_1file
        else:
            lfpdatas = np.concatenate((lfpdatas, lfpdata_1file), axis = 2)
          
    
    return lfpdatas, chnAreas



def ciCOHs_alltrials(lfptrials):
    """
        calculate ciCOH from all trials

        arg:
            lfptrials: (ntemp, nchns, ntrials)


        return:
            ciCOH: ciCOH values (nchns, nchns) 

    """ 

    nchns = lfptrials.shape[1]

    ciCOH = np.zeros((nchns, nchns))
    for chni in range(nchns -1):
        
        sig1 = lfptrials[:, chni, :] # sig1: ntemp * ntrials
        sig1 = np.transpose(sig1) # sig1: ntrials * ntemp
        
        for chnj in range(chni+1, nchns):
            
            sig2 = lfptrials[:, chnj, :]
            sig2 = np.transpose(sig2) # sig2: ntrials * ntemp
            
            ciCOHs = ciCoherence_acrosstrials(sig1, sig2) # ciCOHs: 1 * ntemp
            
            # average across time
            ciCOH[chni, chnj] = np.mean(ciCOHs)
            
            # symmetrical
            ciCOH[chnj, chni] = ciCOH[chni, chnj]
            
            del sig2
        del sig1
        

    return ciCOH