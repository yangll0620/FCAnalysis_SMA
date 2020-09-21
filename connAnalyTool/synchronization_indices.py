#!/usr/bin/env python
# coding: utf-8

# # Functionality: synchronization index calculation


from scipy.signal import hilbert
import numpy as np


# ## corrected imaginary coherence (ciCoherence)

# ### 1. ciCoherence_acrosstrials
# 
# signals connected  in a given  time  should have a stable phase-difference along trials.
#             
# $$
# COH_{a,b}(t) = \frac{1}{N}\sum_{n=1}^{N}\frac{G_{ab}(t,n)}{\sqrt{G_{aa}(t,n)G_{bb}(t,n)}} 
# $$
# 
# $$
# iCOH_{a,b}(t) = \Im(COH_{a,b}(t))
# $$
# 
# $$
# ciCOH_{a,b}(t) = \frac{iCOH_{a,b}(t)}{\sqrt{1-(\Re(COH_{a,b}(t)))^2}}
# $$
# 
# ref:
#     
#     C. J. Stam, G. Nolte, and A. Daffertshofer, “Phase lag index: assessment of functional connectivity 
#     from multi channel EEG and MEG with diminished bias from common sources,” 
#     Human brain mapping, vol. 28, no. 11, pp. 1178–1193, 2007.
#     
#     Meini  Tang,  Yao  Lu  and  Lingling  Yang. Temporal–Spatial Patterns in Dynamic Functional Brain Network for Self-Paced Hand Movement. IEEE Transactions on Neural Systems and Rehabilitation Engineering, 2019


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


# ### 2. ciCoherence_overtime
# 
# assessing ciCoherence as a stable phase-difference over time.
#             
# $$
# COH_{a,b} = \frac{1}{T}\sum_{t=1}^{T}\frac{G_{ab}(t)}{\sqrt{G_{aa}(t)G_{bb}(t)}} 
# $$
# 
# 
# $$
# ciCOH_{a,b}(t) = \frac{\Im(COH_{a,b})}{\sqrt{1-(\Re(COH_{a,b}))^2}}
# $$
# 
# ref:
#     R. Bruña, F. Maestú, and E. Pereda, “Phase Locking Value revisited: teaching new tricks to an old dog,” Journal of neural engineering, vol. 15, no. 5, p. 056011, 2018.

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
            iCOH: corrected imaginary coherence value for signal1 and signal 2 (scalar) 
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

