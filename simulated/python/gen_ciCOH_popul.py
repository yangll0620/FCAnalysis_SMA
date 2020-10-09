import numpy as np
from scipy.signal import hilbert

from simulated_timeseries import gen_series_noconnection_sin



def ciCoherence_acrosstrials(signal1, signal2):
    """
        corrected imaginary coherency function across trials
    
        ref:
            C. J. Stam, G. Nolte, and A. Daffertshofer, “Phase lag index: assessment of functional connectivity 
            from multi channel EEG and MEG with diminished bias from common sources,” 
            Human brain mapping, vol. 28, no. 11, pp. 1178–1193, 2007.

        Args: 
            signal1, signal2: n_epochs * n_times

        Kwargs:
            none

        Returns:
            iCOH: corrected imaginary coherence value for signal1 and signal 2 (1 * n_times) 
    """
    
    # analytic_s1: : n_epochs * n_times
    analytic_s1 = hilbert(signal1, axis = 1)
    
    # analytic_s2: : n_epochs * n_times
    analytic_s2 = hilbert(signal2, axis = 1) 

    
    # phase1_instant: : n_epochs * n_times
    phase1_instant = np.angle(analytic_s1)
    # phase2_instant: : n_epochs * n_times
    phase2_instant = np.angle(analytic_s2)
    # delta_phase2: : n_epochs * n_times
    delta_phase = phase1_instant - phase2_instant 
    A1 = np.abs(analytic_s1) # A1: : n_epochs * n_times
    A2 = np.abs(analytic_s2) # A2: : n_epochs * n_times

    exp_phase = np.exp(1j * delta_phase)
    G_12 = np.mean(np.multiply(np.multiply(A1,A2), exp_phase),axis =0) # numerator : 1 * n_times
    G_11 = np.mean(np.square(A1), axis = 0) # expect_A1Square : 1 * n_times
    G_22 = np.mean(np.square(A2), axis = 0) # expect_A2Square : 1 * n_times

    C = np.divide(G_12, np.sqrt(np.multiply(G_11, G_22)))
    ciCOH = np.divide(np.imag(C), np.sqrt(1-np.square(np.real(C))))

    return ciCOH


def gen_ciCOH_population_sin(ntimes, ntrials, ntemp, f, t, Desired_SNR_dB = 20):
    
        """
        using sinc function to generate the ciCOH population of no connections using sinc function
        
        @paras:
            ntimes: the repeated time (can be set nchns * nchns)
            
            ntrials: the number of trials 
            
            ntemp: the length of the temporal data
            
            f: the frequency of the two time series (Hz)
            
            t: the total time duration for the time series (default 1s)
        
        @return:
            ciCOHs_popul: the generated ciCOH population (ntimes, )
            
        """
        
        
        
        ciCOHs_popul = np.zeros((ntimes))
        for timei in range(ntimes):

            if timei  % 100 ==0:
                print("run the sinc simulation at timei = " + str(timei) + "/" + str(ntimes))

            # generate the two time series for one time sig1, sig2: ntrials * ntemp
            sig1, sig2 = gen_series_noconnection_sin(ntrials = ntrials, ntemp = ntemp, f = f, t = t, 
                                                 Desired_SNR_dB = Desired_SNR_dB)

            # calculate the ciCOH for one time
            ciCOHs = ciCoherence_acrosstrials(sig1, sig2)


            ciCOHs_popul[timei] = np.mean(ciCOHs)

            del ciCOHs
            
        return ciCOHs_popul