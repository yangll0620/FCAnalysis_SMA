import os
import sys
import scipy.io as sio
import numpy as np
import re

codefolder = re.match('.*exp.*code', os.getcwd()).group()
sys.path.append(codefolder)
from connAnalyTool.synchronization_indices import ciCoherence_overtime

def lfp_extract(files):
    """
        extract all rest lfp from files

        Args:
            files: extract using glob.glob

        Returns:
            lfpdata: nareas * ntemp * nsegs

            chnAreas: list for area name

    """
        
    if 'lfpdata' in locals():
        del lfpdata
    
    for i, file in enumerate(files):
        
        ### load data
        matdat = sio.loadmat(file, variable_names = ['lfpsegs', 'lfpdata', 'fs', 'chnAreas'], 
                             struct_as_record = False, squeeze_me = True) 
        
        
        
        ### extract the noused channels, only calculate once
        if i == 0:
            
            # chnAreas
            chnAreas = matdat['chnAreas'].tolist()
            
            # fs: sample rate
            fs = matdat['fs'] 
             
        

        ### dealing lfp data
        
        # lfp (np.ndarray): nareas * ntemp * ntrials or ntemp * nareas * ntrials
        if 'lfpdata' in matdat.keys():
            lfpdata_1file = matdat['lfpdata']
        elif 'lfpsegs' in matdat.keys():
            lfpdata_1file = matdat['lfpsegs']

        n1, n2, n3 = lfpdata_1file.shape
        if n1 > n2:  # ntemp * nareas * ntrials
            lfpdata_1file = np.transpose(lfpdata_1file, (1, 0, 2))
        
        # concatenate to lfpdata for all files
        if 'lfpdata' not in locals():
            lfpdata = lfpdata_1file
        else:
            lfpdata = np.concatenate((lfpdata, lfpdata_1file), axis = 2)
          
    
    return lfpdata, chnAreas, fs



def calc_ciCOHs_rest(lfpdata):
    """
        calculate the ciCOHs for rest using ciCoherence_overtime

        Arg:
            lfpdata: nchns * ntemp * nsegs

        Return:
            ciCOH: nchns * nchns
    """

    nchns, _, nsegs = lfpdata.shape
    ciCOHs = np.zeros((nchns, nchns, nsegs))
    for segi in range(nsegs):
        
        for chni in range(nchns -1):
            signal1 = lfpdata[chni, :, segi]
            
            for chnj in range(chni+1, nchns):
                signal2 = lfpdata[chnj, :, segi]
                
                # ciCOHs assignment
                ciCOHs[chni, chnj, segi] = ciCoherence_overtime(signal1, signal2)
                
                # symmetrical
                ciCOHs[chnj, chni, segi] = ciCOHs[chni, chnj, segi]
                
                del signal2
            del signal1
            
    ciCOH = np.mean(ciCOHs, axis = 2)

    return ciCOH