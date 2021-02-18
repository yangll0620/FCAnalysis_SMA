import os, sys
import numpy as np
from scipy.stats import norm
from mne.stats import fdr_correction

import re
codefolder = re.match('.*exp.*code', os.getcwd()).group()
sys.path.append(codefolder)
from connAnalyTool.synchronization_indices import ciCoherence_acrosstrials

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


def dynciCOH_from_lfptrials(lfptrials):
    """
    calculate dynamic ciCOH from lfptrials

    args:
        lfptrials:  nchns * ntemp * ntrials

    return:
        ciCOH: nchns * nchns *  ntemp
    """
    
    nchns, ntemp, _ = lfptrials.shape

    dynciCOH = np.zeros((nchns, nchns, ntemp))
    for chni in range(nchns -1):
        
        sig1 = lfptrials[chni, :, :] # sig1: ntemp * ntrials
        sig1 = np.transpose(sig1) # sig1: ntrials * ntemp
        
        for chnj in range(chni+1, nchns):
            
            sig2 = lfptrials[chnj, :, :]
            sig2 = np.transpose(sig2) # sig2: ntrials * ntemp
            
            dynciCOH[chni, chnj, :] = ciCoherence_acrosstrials(sig1, sig2) # ciCOHs: 1 * ntemp
            
            
            # symmetrical
            dynciCOH[chnj, chni, :] = dynciCOH[chni, chnj, :]
            
            del sig2
        del sig1
        
    return dynciCOH




def permdist_dynciCOH_SKT(lfp1, lfp2, shuffleN = 100):
    """
        permutation test distribution

        Arg:
            lfp1, lfp2:  ntemp * nsegs(ntrials)

            shuffleN: the total shuffle times


        Return:
            mu, std: the mu and std of the fitted normal distribution
    """

    permlfp1, permlfp2 = lfp1.copy(), lfp2.copy()
    ntemp = permlfp1.shape[0]
    perm_dynciCOHs = np.zeros(shape=(shuffleN,ntemp))
    for i in range(shuffleN):

        # shuffle permlfp2 along ntrials
        permlfp2 = np.transpose(a = permlfp2, axes = (1, 0))
        np.random.shuffle(permlfp2)
        permlfp2 = np.transpose(a = permlfp2, axes = (1, 0))

        permlfp = np.concatenate((np.expand_dims(permlfp1, axis = 0), np.expand_dims(permlfp2, axis = 0)), axis = 0)

        perm_dynciCOHs[i, :] = dynciCOH_from_lfptrials(permlfp)[0, 1, :]
        

        del permlfp


    # Fit a normal distribution to the data:
    perm_vec = np.reshape(a = perm_dynciCOHs, newshape = (shuffleN * ntemp, ))
    mu, std = norm.fit(perm_vec)

    return mu, std

    


def pval_perm_dynciCOH_SKT(dynciCOH, lfptrials):
    """
        pvalues using permutation test for dynamic ciCOHs

        Arg:
            dynciCOH: dynamic ciCOHs [nchns * nchns * ntemp]

            lfptrials: the lfp trial data calculatingthe dynciCOH

        Return:
            pvals: p-value for each value in dynciCOH, shape = dynciCOH.shape
    """
    
    # 
    [i, j, _] = np.unravel_index(np.argmax(dynciCOH), shape = dynciCOH.shape)
    lfp1, lfp2 = lfptrials[i, :, :], lfptrials[j, :, :]
    mu, std = permdist_dynciCOH_SKT(lfp1, lfp2, shuffleN = 100)
    pvals = norm.sf(abs(dynciCOH), loc = mu, scale = std) * 2
    
    
    return pvals



def truncate_dynfc(dynciCOH, pvals):
    """
        truncate fc to be 0 if not significant



        Arg:
            dynciCOH: dynamic ciCOHs [nchns * nchns * ntemp]

            pvals: p-value for each value in dynciCOH, shape = dynciCOH.shape


        Return:

            trunc_dynfc: truncated dynamic fc (value is 0 or 1)

    """

    # multiple comparison correction, get truncate dynfc
    reject, _ = fdr_correction(pvals, alpha = 0.05, method='indep')
    [rows, cols, ts]= np.where(reject == True)
    trunc_dynfc = np.zeros(dynciCOH.shape)
    if len(rows) > 0:
        trunc_dynfc[rows, cols, ts] = abs(dynciCOH[rows, cols, ts])


    return trunc_dynfc


def digitized_dynfc(dynciCOH, pvals):
    """
        digitized fc to be 1 or 0



        Arg:
            dynciCOH: dynamic ciCOHs [nchns * nchns * ntemp]

            pvals: p-value for each value in dynciCOH, shape = dynciCOH.shape


        Return:

            digi_dynfc: digitized dynamic fc (value is 0 or 1)

    """

    # multiple comparison correction, get digitized dynfc
    reject, _ = fdr_correction(pvals, alpha = 0.05, method='indep')
    [rows, cols, ts]= np.where(reject == True)
    digi_dynfc = np.zeros(dynciCOH.shape)
    if len(rows) > 0:
        digi_dynfc[rows, cols, ts] = 1


    return digi_dynfc

