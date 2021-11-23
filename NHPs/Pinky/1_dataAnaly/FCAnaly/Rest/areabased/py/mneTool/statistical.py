import numpy as np
from scipy.stats import norm
from lfp_ciCOH_extract import calc_ciCOHs_rest
from mne.connectivity import spectral_connectivity

def pval_permciCOH_rest(lfp1, lfp2, actciCOH, shuffleN = 1000):
    """
        

        Arg:
            lfp1, lfp2:  ntemp * nsegs(ntrials)

            shuffleN: the total shuffle times

            actciCOH: an actual ciCOH value (positive or negative)

        Return:
            pval: the p-value base on permutation test 
            mu, std: the mu and std of the fitted normal distribution
    """

    permlfp1, permlfp2 = lfp1.copy(), lfp2.copy()
    permciCOHs = np.zeros(shape=(shuffleN,))
    for i in range(shuffleN):

        # shuffle permlfp2
        permlfp2 = np.transpose(permlfp2, axes=(1, 0))
        np.random.shuffle(permlfp2)
        permlfp2 = np.transpose(permlfp2, axes=(1, 0))

        permlfp = np.concatenate((np.expand_dims(permlfp1, axis = 0), np.expand_dims(permlfp2, axis = 0)), axis = 0)
        
        ciCOHM = calc_ciCOHs_rest(permlfp)

        permciCOHs[i] = ciCOHM[0, 1]

        del ciCOHM, permlfp


    # Fit a normal distribution to the data:
    mu, std = norm.fit(permciCOHs)

    pval = norm.sf(abs(actciCOH), loc = mu, scale = std) 

    return pval, mu, std


def pval_perm_imcohMNE_rest(lfp1, lfp2, fs, fmin, fmax, shuffleN = 1000):
    """
        

        Arg:
            lfp1, lfp2:  ntemp * nsegs(ntrials)

            shuffleN: the total shuffle times

            fs: sampling rate

            fmin, fmax: min and max interested frequency

        Return:
            mu, std: the mu and std of the fitted normal distribution
    """

    permlfp1, permlfp2 = lfp1.copy(), lfp2.copy()
    perm_imcohs = np.zeros(shape=(shuffleN,))
    for i in range(shuffleN):

        # shuffle permlfp2 along ntemp
        np.random.shuffle(permlfp2)

        permlfp = np.concatenate((np.expand_dims(permlfp1, axis = 0), np.expand_dims(permlfp2, axis = 0)), axis = 0)

        [imcohs, _, _, _, _]= spectral_connectivity(data = np.transpose(permlfp, axes = (2, 0, 1)), 
                                            method= 'imcoh', sfreq = fs, fmin= fmin, fmax = fmax, faverage = True, n_jobs = 100)
        

        perm_imcohs[i] = imcohs[1, 0, 0]

        del imcohs, permlfp


    # Fit a normal distribution to the data:
    mu, std = norm.fit(perm_imcohs)

    return mu, std

