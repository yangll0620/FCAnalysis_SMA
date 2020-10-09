import numpy as np

from addnoise import add_gaussiannoise_SNR

def gen_series_noconnection_sin(ntrials, ntemp, f, t, Desired_SNR_dB = 20):
    """
        generate two time series data using sinusodal function without phase locking 
        (i.e. the phase diff varies in [0 2*pi}])
        
        inputs:
            
            ntrials: the trials number
            
            ntemp: the total temporal number
            
            f: the frequency of the two time series (Hz)
            
            t: the total time duration for the time series (default 1s)
            
            Desired_SNR_dB: the add desired SNR (dB) gaussian noise to signal
        
        return:
        
            sig1, sig2: the generated no phase locking two sinusodal time series (ntrials * ntemp)
    """
    
    ts = np.linspace(0, t, ntemp)
    
    # generated sigs1, s2 signals:  ntrials * ntemp
    sig1, sig2 = np.empty(shape=[0,ntemp]), np.empty(shape=[0,ntemp])
    for triali in range(ntrials):

        # random phase diff in range [0 2*pi)
        phi = 2 * np.pi * np.random.rand(1) 

        # generate sin time series s1 and s2 (phase diff is phi)
        s1 = np.sin(2 * np.pi * f * ts)
        s2 = np.sin(2 * np.pi * f * ts + phi)


        # add normal distribution noise
        s1 = add_gaussiannoise_SNR(s1, Desired_SNR_dB = Desired_SNR_dB)
        s2 = add_gaussiannoise_SNR(s2, Desired_SNR_dB = Desired_SNR_dB)



         # append the time serie of the new trial
        sig1 = np.append(sig1, np.expand_dims(s1, axis = 0), axis=0)
        sig2 = np.append(sig2, np.expand_dims(s2, axis = 0), axis=0)

        del phi, s1, s2

    return sig1, sig2