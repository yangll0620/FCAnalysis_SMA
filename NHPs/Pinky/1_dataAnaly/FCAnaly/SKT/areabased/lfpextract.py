import os, sys
import scipy.io as sio
import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle 
import math
import pandas as pd


from enum import Enum
class Event(Enum):
    TARGETONSET = 0
    REACHONSET = 1
    REACH = 2
    RETURNONSET = 3
    MOUTH = 4


def lfp_align2_targetonset(files, tdur_trial, tmin_return, tmax_return):
    """
        extract lfp data respect to targetonset

        return:
            lfptrials: nchns * ntemp * ntrials
    """

    variablesinLoadfile = ['lfpdata', 'fs', 'chnAreas', 'idxevent']

    coli_targetonset = 0
    
    tdur_trial = np.array(tdur_trial)

    if 'lfptrials' in locals():
        del lfptrials
    
    for i, file in enumerate(files):
        
        ### load data
        matdat = sio.loadmat(file, variable_names = variablesinLoadfile, 
                             struct_as_record = False, squeeze_me = True) 
        
        
        ### extract chnAreas
        if i == 0:  
            chnAreas = matdat['chnAreas'].tolist()
            fs = matdat['fs']
        

        ### dealing lfp data
        
        # lfp (np.ndarray): ntemporal * nchns * ntrials
        lfpdata_1file = matdat['lfpdata']
        idxevent = matdat['idxevent']
        

        ntrials = lfpdata_1file.shape[2]
        for tri in range(ntrials):

            idxdur = np.round(tdur_trial * fs).astype(int) + idxevent[tri, coli_targetonset]
            lfp_phase_1trial = lfpdata_1file[:,idxdur[0]:idxdur[1], tri]
            lfp_phase_1trial = np.expand_dims(lfp_phase_1trial, axis = 2)

            if 'lfptrials' not in locals():
                lfptrials = lfp_phase_1trial
            else:
                lfptrials = np.concatenate((lfptrials, lfp_phase_1trial), axis = 2)

    
    return lfptrials, chnAreas, fs

def lfp_align2_returnonset(files, tdur_trial = [0, 0.5], tmin_return = 0.5, tmax_return = 1):
    """
        extract lfp data respect to returnonsest

        return:
            lfptrials: nchns * ntemp * ntrials
    """

    variablesinLoadfile = ['lfpdata', 'fs', 'chnAreas', 'idxevent']

    coli_returnonset = 3
    coli_mouth = 4
    
    tdur_trial = np.array(tdur_trial)

    if 'lfptrials' in locals():
        del lfptrials
    
    for i, file in enumerate(files):
        
        ### load data
        matdat = sio.loadmat(file, variable_names = variablesinLoadfile, 
                             struct_as_record = False, squeeze_me = True) 
        
        
        ### extract chnAreas
        if i == 0:  
            chnAreas = matdat['chnAreas'].tolist()
            fs = matdat['fs']
        

        ### dealing lfp data
        
        # lfp (np.ndarray): ntemporal * nchns * ntrials
        lfpdata_1file = matdat['lfpdata']
        idxevent = matdat['idxevent']
        

        ntrials = lfpdata_1file.shape[2]
        for tri in range(ntrials):

            t_return = (idxevent[tri, coli_mouth] - idxevent[tri, coli_returnonset]) / fs
            if t_return < tmin_return or t_return > tmax_return:
                continue

            idxdur = np.round(tdur_trial * fs).astype(int) + idxevent[tri, coli_returnonset]
            lfp_phase_1trial = lfpdata_1file[:,idxdur[0]:idxdur[1], tri]
            lfp_phase_1trial = np.expand_dims(lfp_phase_1trial, axis = 2)

            if 'lfptrials' not in locals():
                lfptrials = lfp_phase_1trial
            else:
                lfptrials = np.concatenate((lfptrials, lfp_phase_1trial), axis = 2)

    
    return lfptrials, chnAreas, fs



def lfp_align2_reachonset(files, tdur_trial = [0, 0.5], tmin_reach = 0.5, tmax_reach = 1):
    """
        extract lfp data respect to reachonset

        return:
            lfptrials: nchns * ntemp * ntrials
    """

    variablesinLoadfile = ['lfpdata', 'fs', 'chnAreas', 'idxevent']

    coli_reachonset = 1
    coli_reach = 2
    
    tdur_trial = np.array(tdur_trial)

    if 'lfptrials' in locals():
        del lfptrials
    
    for i, file in enumerate(files):
        
        ### load data
        matdat = sio.loadmat(file, variable_names = variablesinLoadfile, 
                             struct_as_record = False, squeeze_me = True) 
        
        
        ### extract chnAreas
        if i == 0:  
            chnAreas = matdat['chnAreas'].tolist()
            fs = matdat['fs']
        

        ### dealing lfp data
        
        # lfp (np.ndarray): ntemporal * nchns * ntrials
        lfpdata_1file = matdat['lfpdata']
        idxevent = matdat['idxevent']
        

        ntrials = lfpdata_1file.shape[2]
        for tri in range(ntrials):

            t_reach = (idxevent[tri, coli_reach] - idxevent[tri, coli_reachonset]) / fs
            if t_reach < tmin_reach or t_reach > tmax_reach:
                continue

            idxdur = np.round(tdur_trial * fs).astype(int) + idxevent[tri, coli_reachonset]
            lfp_phase_1trial = lfpdata_1file[:,idxdur[0]:idxdur[1], tri]
            lfp_phase_1trial = np.expand_dims(lfp_phase_1trial, axis = 2)

            if 'lfptrials' not in locals():
                lfptrials = lfp_phase_1trial
            else:
                lfptrials = np.concatenate((lfptrials, lfp_phase_1trial), axis = 2)

    
    return lfptrials, chnAreas, fs



def lfp_align2(files, align2 = Event.REACHONSET, tdur_trial = [0, 0.5], t_minmax_reach = None, t_minmax_return = None):
    """
        extract lfp data respect to targetonset, reachonset, reach and returnonset separately


        Args:
            align2: the event to be aligned (e.g Event.REACHONSET)

            tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5, 0.6])
            
            t_minmax_reach, t_minmax_return : min and max reach/return (s) for selecting trials (e.g [0.5, 1])

        return:
            lfptrials: nchns * ntemp * ntrials

            chnAreas:

            fs:
    """

    variablesinLoadfile = ['lfpdata', 'fs', 'chnAreas', 'idxevent']

    coli_targetonset, coli_reachonset,  coli_reach, coli_returnonset, coli_mouth = 0, 1, 2, 3, 4
    
    tdur_trial = np.array(tdur_trial)

    if align2 == Event.TARGETONSET:
        coli_align2 = coli_targetonset
    elif align2 == Event.REACHONSET:
        coli_align2 = coli_reachonset
    elif align2 == Event.REACH:
        coli_align2 = coli_reach
    elif align2 == Event.RETURNONSET:
        coli_align2 = coli_returnonset
    else:
        print("align2 Event is not valid. ")
        lfptrials, chnAreas, fs = [], [],  []
        return lfptrials, chnAreas, fs

    if 'lfptrials' in locals():
        del lfptrials
    
    for i, file in enumerate(files):
        
        ### load data
        matdat = sio.loadmat(file, variable_names = variablesinLoadfile, 
                             struct_as_record = False, squeeze_me = True) 
        
        
        ### extract chnAreas
        if i == 0:  
            chnAreas = matdat['chnAreas'].tolist()
            fs = matdat['fs']
        

        ### dealing lfp data
        
        # lfp (np.ndarray): ntemporal * nchns * ntrials
        lfpdata_1file = matdat['lfpdata']
        idxevent = matdat['idxevent']
        

        ntrials = lfpdata_1file.shape[2]
        for tri in range(ntrials):

            t_reach = (idxevent[tri, coli_reach] - idxevent[tri, coli_reachonset]) / fs
            t_return = (idxevent[tri, coli_mouth] - idxevent[tri, coli_returnonset]) / fs
            if t_minmax_reach is not None and (t_reach < t_minmax_reach[0] or t_reach > t_minmax_reach[1]):
                continue
            if t_minmax_return is not None and (t_return < t_minmax_return[0] or t_return >  t_minmax_return[1]):
                continue

            idxdur = np.round(tdur_trial * fs).astype(int) + idxevent[tri, coli_align2]
            lfp_phase_1trial = lfpdata_1file[:,idxdur[0]:idxdur[1], tri]
            lfp_phase_1trial = np.expand_dims(lfp_phase_1trial, axis = 2)

            if 'lfptrials' not in locals():
                lfptrials = lfp_phase_1trial
            else:
                lfptrials = np.concatenate((lfptrials, lfp_phase_1trial), axis = 2)

    
    return lfptrials, chnAreas, fs
    