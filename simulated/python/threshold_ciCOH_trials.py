from scipy.stats import norm
import numpy as np


from gen_ciCOH_popul import gen_ciCOH_population_sin

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