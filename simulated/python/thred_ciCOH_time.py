from scipy.stats import norm
import numpy as np
from scipy.signal import hilbert
from scipy import stats


def add_gaussiannoise_SNR(signal, Desired_SNR_dB):
    """
        add desired SNR (dB) gaussian noise to signal
        SNR(dB) = 10 * log10(power_signal/power_noise)

        @ parameter:
            signal: (n_times,)
            Desired_SNR_dB: desired SNR in dB
        
        @ return 
            signal_noisy: (n_times,)
    """
    
    n_times = signal.shape[0]
    noise = np.random.normal(loc=0.0, scale=1.0, size=(n_times,))

    power_signal = np.dot(abs(signal), abs(signal))/n_times
    power_noise = np.dot(abs(noise), abs(noise))/n_times 

    k = (power_signal * pow(10,(-Desired_SNR_dB/10)))/power_noise # scale factor
    noise_new = np.sqrt(k) * noise

    signal_noisy = signal + noise_new

    return signal_noisy


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


def ciCoherence_overtime(signal1, signal2):
    """
        corrected imaginary coherency function over time

        ref:
            R. Bruña, F. Maestú, and E. Pereda, “Phase Locking Value revisited: teaching new tricks to an old dog,” 
            Journal of neural engineering, vol. 15, no. 5, p. 056011, 2018.
            
        Args: 
            signal1, signal2: 1 * n_times

        Kwargs:
            none

        Returns:
            ciCOH: corrected imaginary coherence value for signal1 and signal 2 (scalar) 
    """
    
    # analytic_s1: : 1 * n_times
    analytic_s1 = hilbert(signal1)

    # analytic_s2: : 1 * n_times
    analytic_s2 = hilbert(signal2) 



    # phase1_instant: : 1 * n_times
    phase1_instant = np.angle(analytic_s1)
    
    # phase2_instant: : 1 * n_times
    phase2_instant = np.angle(analytic_s2) 
    
    # delta_phase2: : 1 * n_times
    delta_phase = phase1_instant - phase2_instant 
    
    # exponent of delta_phase
    exp_phase = np.exp(1j * delta_phase)


    
    # A1: Amplitude of analytic_s1 (1 * n_times)
    A1 = np.abs(analytic_s1)
    # A2: Amplitude of analytic_s1 (1 * n_times)
    A2 = np.abs(analytic_s2) 



    # complex cross-spectrum (1 * times) of two signals
    G_12 = np.multiply(np.multiply(A1,A2), exp_phase)

    # auto-spectrum of signal1 and signal2 individually (1 * times)
    G_11, G_22 = np.square(A1), np.square(A2)



    # coherence (1 * 1)
    C = np.mean(np.divide(G_12, np.sqrt(np.multiply(G_11, G_22))))

    # corrected imaginary coherence
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
            ciCOHs = np.zeros((ntrials))
            for triali in range(ntrials):

                s1, s2 = sig1[triali, :], sig2[triali, :]
                ciCOHs[triali] = ciCoherence_overtime(s1, s2)

                del s1, s2


            ciCOHs_popul[timei] = np.mean(ciCOHs)

            del ciCOHs
            
        return ciCOHs_popul


def threshold_ciCOH_sin(ntimes, ntrials, ntemp, f, t, alpha = 0.05, ploton = True):
    
        """
        using sinc function to simulated the no connections and identify the threshold for has connection
        
        @paras:
            ntimes: the repeated time (can be set nchns * nchns)
            
            ntrials: the number of trials 
            
            ntemp: the length of the temporal data
            
            f: the frequency of the two time series (Hz)
            
            t: the total time duration for the time series (default 1s)
            
            alpha: the critirial (5%, or 1%)
            
            ploton: show a figure if True
        
        @return:
            threshold: the ciCOH threshold value (a scalar)
            
            
        Example Usage:
        
            thres = threshold_ciCOH_sin(ntimes = 1000, ntrials = 100, ntemp = 500, f = 30, t = 1, alpha = 0.01, ploton = True)
            
        """
        
        print("identifying the ciCOH threshold using sinc.....")
        
        # generate the ciCOH population 
        ciCOHs = gen_ciCOH_population_sin(ntimes = ntimes, ntrials = ntrials, ntemp = ntemp, f = f, t = t, Desired_SNR_dB = 20)
        
        
        
        ### Identify the threshold
        
        mu, std = norm.fit(ciCOHs) # Fit a normal distribution to the data
        if alpha == 0.05:
            threshold = np.around(mu + 2*std, decimals=2)
        
        elif alpha == 0.01:
            threshold = np.around(mu + 3*std, decimals=2)

        
        
        
        ### plot the ciCOH distribution of ciCOHs of the simulated data 
        if ploton:
            text_title = "ntrials = " +  str(ntrials) + ", ntemp = " + str(ntemp) +\
                        ",f = "  + str(f) + "Hz, ntimes = " + str(ntimes)
            
            _plot_ciCOH_distributions(ciCOHs, text_title = text_title)
            
            
    
        print("threshold = " + str(threshold) + ", mu = " + str(mu) + ", std = " + str(std))
            
            
        return threshold, mu, std




## Multiple comparison correction
def corr_threshold_ciCOH_sin_BH(ciCOHs_actual, ntimes, ntrials, ntemp, f, t, false_rate = 0.25, mu = None, std = None):
    """
    
        control the false discovery reate using Benjamini-Hochberg procedure
        
        
        @paras:
            
            ciCOHs_actual:the ciCOH values (nvalues, )
            
            ntimes: the repeated time (can be set nchns * nchns)
            
            ntrials: the number of trials 
            
            ntemp: the length of the temporal data
            
            mu, std: the mean and std values of the null hypothesis no connection probability distribution
            
            f: the frequency of the two time series (Hz)
            
            t: the total time duration for the time series (default 1s)
            
            
            
        ref:
            http://www.biostathandbook.com/multiplecomparisons.html
        
    """
    
            
    print("identifying the ciCOH corrected threshold using sinc and Benjamini-Hochberg procedure....")

    if mu is None or std is None:
    
        # generate the ciCOH population 
        ciCOHs_simulated = gen_ciCOH_population_sin(ntimes = ntimes, ntrials = ntrials, ntemp= ntemp, f = f, t = 1, Desired_SNR_dB = 20)

        # evaluate the null hypothesis distribution mean and std using Normal distribution
        mu, std = norm.fit(ciCOHs_simulated)

    
    
    ### BH corrected
    
    ciCOHs_actual = abs(ciCOHs_actual)
    
    # calculate the pvalue for each abs(ciCOH) value
    pvalues = stats.norm.sf(ciCOHs_actual, loc = mu, scale = std) * 2
    
    # sort the pvalues
    ind = np.argsort(pvalues)
    pvalues_sorted = pvalues[ind]
    ciCOHs_sorted = ciCOHs_actual[ind]
    
    # generate the Benjamini-Hochberg critical values
    criticalvalue_bh = np.array([*range(1, len(pvalues) + 1, 1)])/len(pvalues) * false_rate

    # find the index which has the largest P value what pvalue < criticalvalue_bh
    for i in range(len(pvalues_sorted)):
        if(pvalues_sorted[i] >= criticalvalue_bh[i]):
            break
    
    
    
    corrected_threshold = ciCOHs_sorted[i-1]
    
    print("corrected threshold = " + str(corrected_threshold))
    
    return corrected_threshold, mu, std