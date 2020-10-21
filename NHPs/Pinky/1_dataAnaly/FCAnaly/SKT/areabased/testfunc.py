import pandas as pd
import numpy as np
import os, sys
import re
codefolder = re.match('.*exp.*code', __file__).group()
sys.path.append(codefolder)
from connAnalyTool.synchronization_indices import ciCoherence_acrosstrials

def assign_coord2chnArea(chnInf_file, chnAreas):
    # load channel coord from chnInf_file
    df = pd.read_csv(chnInf_file, header = 0)

    # fill in the x,y coordinates of each area in chnAreas based on the values in df_chninf
    coord_x, coord_y = np.zeros(shape = [len(chnAreas), ]), np.zeros(shape = [len(chnAreas), ])
    for i, chnArea in enumerate(chnAreas):
        
        mask_area = (df['brainarea']==chnArea)

        x, y = df['simulated_x'][mask_area].to_numpy(), df['simulated_y'][mask_area].to_numpy()

        coord_x[i], coord_y[i] = x, y
        
        del mask_area, x, y

    df_chninf = pd.DataFrame(data = {'chnAreas': chnAreas, 'coord_x': coord_x, 'coord_y': coord_y})
        
    return df_chninf


def calcciCOH_from_lfptrials(lfptrials):
    """
    calculate ciCOH from lfptrials

    args:
        lfptrials:  nchns * ntemp * ntrials

    return:
        ciCOH: nchns * nchns
    """
    
    nchns = lfptrials.shape[0]

    ciCOH = np.zeros((nchns, nchns))
    for chni in range(nchns -1):
        
        sig1 = lfptrials[chni, :, :] # sig1: ntemp * ntrials
        sig1 = np.transpose(sig1) # sig1: ntrials * ntemp
        
        for chnj in range(chni+1, nchns):
            
            sig2 = lfptrials[chnj, :, :]
            sig2 = np.transpose(sig2) # sig2: ntrials * ntemp
            
            ciCOHs = ciCoherence_acrosstrials(sig1, sig2) # ciCOHs: 1 * ntemp
            
            # average across time
            ciCOH[chni, chnj] = np.mean(ciCOHs)
            
            # symmetrical
            ciCOH[chnj, chni] = ciCOH[chni, chnj]
            
            del sig2
        del sig1
        
    return ciCOH


def dyciCOH_from_lfptrials(lfptrials):
    """
    calculate ciCOH from lfptrials

    args:
        lfptrials:  nchns * ntemp * ntrials

    return:
        ciCOH: nchns * nchns
    """
    
    nchns, ntemp, _ = lfptrials.shape

    ciCOH = np.zeros((nchns, nchns, ntemp))
    for chni in range(nchns -1):
        
        sig1 = lfptrials[chni, :, :] # sig1: ntemp * ntrials
        sig1 = np.transpose(sig1) # sig1: ntrials * ntemp
        
        for chnj in range(chni+1, nchns):
            
            sig2 = lfptrials[chnj, :, :]
            sig2 = np.transpose(sig2) # sig2: ntrials * ntemp
            
            ciCOHs = ciCoherence_acrosstrials(sig1, sig2) # ciCOHs: 1 * ntemp
            
            # average across time
            ciCOH[chni, chnj, :] = ciCOHs
            
            # symmetrical
            ciCOH[chnj, chni, :] = ciCOH[chni, chnj]
            
            del sig2
        del sig1
        
    return ciCOH



def dyfc4drawing(chnAreas, fs, chnInf_file, lfptrials_normal = None, lfptrials_mild = None, lfptrials_moderate = None):
    """
        args:
            lfptrials_normal: nchns * ntemp * ntrials

        return:
            fc: dictionary for drawing functional connectivity 
    """

    if lfptrials_normal is None and lfptrials_mild is None and lfptrials_moderate is None:
        print('no lfptrials is given!')
        return None

    ciCOH = dict()
    setup = dict()
    setup['fs'] = fs


    if lfptrials_normal is not None and lfptrials_mild is not None:
        ntrials_normal = lfptrials_normal.shape[2]
        ntrials_mild = lfptrials_mild.shape[2]

        ntrials = np.min([ntrials_normal, ntrials_mild])

        lfptrials_normal = lfptrials_normal[:, :, 0:ntrials]
        lfptrials_mild = lfptrials_mild[:, :, 0:ntrials]


    if lfptrials_normal is not None:
        ciCOH_normal = dyciCOH_from_lfptrials(lfptrials_normal)
        ciCOH['normal']  = ciCOH_normal
        
        setup['ntemp_normal'] = lfptrials_normal.shape[1]
        setup['ntrials_normal'] = lfptrials_normal.shape[2]
    
    if lfptrials_mild is not None:
        ciCOH_mild = dyciCOH_from_lfptrials(lfptrials_mild)
        ciCOH['mild'] = ciCOH_mild

        setup['ntemp_mild'] = lfptrials_mild.shape[1]
        setup['ntrials_mild'] = lfptrials_mild.shape[2]

       
    fc = dict()
    fc['ciCOH'] = ciCOH
    fc['chnInf'] = assign_coord2chnArea(chnInf_file, chnAreas)
    fc['setup']= setup

    return fc



def fc4drawing(chnAreas, fs, chnInf_file, lfptrials_normal = None, lfptrials_mild = None, lfptrials_moderate = None):
    """
        args:
            lfptrials_normal: nchns * ntemp * ntrials

        return:
            fc: dictionary for drawing functional connectivity 
    """

    if lfptrials_normal is None and lfptrials_mild is None and lfptrials_moderate is None:
        print('no lfptrials is given!')
        return None

    ciCOH = dict()
    setup = dict()
    setup['fs'] = fs


    if lfptrials_normal is not None and lfptrials_mild is not None:
        ntrials_normal = lfptrials_normal.shape[2]
        ntrials_mild = lfptrials_mild.shape[2]

        ntrials = np.min([ntrials_normal, ntrials_mild])

        lfptrials_normal = lfptrials_normal[:, :, 0:ntrials]
        lfptrials_mild = lfptrials_mild[:, :, 0:ntrials]


    if lfptrials_normal is not None:
        ciCOH_normal = calcciCOH_from_lfptrials(lfptrials_normal)
        ciCOH['normal']  = ciCOH_normal
        
        setup['ntemp_normal'] = lfptrials_normal.shape[1]
        setup['ntrials_normal'] = lfptrials_normal.shape[2]
    
    if lfptrials_mild is not None:
        ciCOH_mild = calcciCOH_from_lfptrials(lfptrials_mild)
        ciCOH['mild'] = ciCOH_mild

        setup['ntemp_mild'] = lfptrials_mild.shape[1]
        setup['ntrials_mild'] = lfptrials_mild.shape[2]

       
    fc = dict()
    fc['ciCOH'] = ciCOH
    fc['chnInf'] = assign_coord2chnArea(chnInf_file, chnAreas)
    fc['setup']= setup

    return fc