import os
import sys
import scipy.io as sio
import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle 
import math
import pandas as pd
import cv2


import re
codefolder = re.match('.*exp.*code', __file__).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder



def lfpallfiles_extract(files):
    """
        extract all rest segments from files

        Args:
            files: extract using glob.glob

        Returns:
            lfpdata:

            chnAreas

    """
        
    if 'lfpdata' in locals():
        del lfpdata
    
    for i, file in enumerate(files):
        
        ### load data
        matdat = sio.loadmat(file, variable_names = ['lfpsegs', 'fs', 'chnAreas'], 
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


def main():
    _, _, pipelinefolder, _= exp_subfolders()
    corresfolder, correparentfolder = code_corresfolder(__file__)

    animal =  re.search('NHPs/[a-zA-Z]*/', __file__).group()[len('NHPs/'):-1]


if __name__ == "__main__":
     main()
