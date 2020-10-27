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
from statsmodels.stats.multitest import multipletests


import re
codefolder = re.match('.*exp.*code', __file__).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder
from connAnalyTool.synchronization_indices import ciCoherence_overtime
from connAnalyTool.fc_visual_time import threshold_fc_overtime, ciCOH_visual_save, pvals_fc_overtime


def lfp_extract(files):
    """
        extract all rest lfp from files

        Args:
            files: extract using glob.glob

        Returns:
            lfpdata: nareas * ntemp * nsegs

            chnAreas: list for area name

    """
        
    if 'lfpdata' in locals():
        del lfpdata
    
    for i, file in enumerate(files):
        
        ### load data
        matdat = sio.loadmat(file, variable_names = ['lfpsegs', 'lfpdata', 'fs', 'chnAreas'], 
                             struct_as_record = False, squeeze_me = True) 
        
        
        
        ### extract the noused channels, only calculate once
        if i == 0:
            
            # chnAreas
            chnAreas = matdat['chnAreas'].tolist()
            
            # fs: sample rate
            fs = matdat['fs'] 
             
        

        ### dealing lfp data
        
        # lfp (np.ndarray): nareas * ntemp * ntrials or ntemp * nareas * ntrials
        if 'lfpdata' in matdat.keys():
            lfpdata_1file = matdat['lfpdata']
        elif 'lfpsegs' in matdat.keys():
            lfpdata_1file = matdat['lfpsegs']

        n1, n2, n3 = lfpdata_1file.shape
        if n1 > n2:  # ntemp * nareas * ntrials
            lfpdata_1file = np.transpose(lfpdata_1file, (1, 0, 2))
        
        # concatenate to lfpdata for all files
        if 'lfpdata' not in locals():
            lfpdata = lfpdata_1file
        else:
            lfpdata = np.concatenate((lfpdata, lfpdata_1file), axis = 2)
          
    
    return lfpdata, chnAreas, fs





def assign_coord2chnArea(area_coord_file, chnAreas):
    """
        assign the xy coord of each area


        Args:
            area_coord_file: file containing the x y coord for each area (normally predefined)


            chnAreas: a list of areas representing the corresponding area for each channel


        Return:
            df_chninf: DataFrame with chnAreas, coord_x and coord_y

    """
    # load channel coord from area_coord_file
    df = pd.read_csv(area_coord_file, header = 0)

    # fill in the x,y coordinates of each area in chnAreas based on the values in df_chninf
    coord_x, coord_y = np.zeros(shape = [len(chnAreas), ]), np.zeros(shape = [len(chnAreas), ])
    for i, chnArea in enumerate(chnAreas):
        
        mask_area = (df['brainarea'] == chnArea)
        
        if len(df['brainarea'][mask_area].index) == 0:
            continue

        x, y = df['simulated_x'][mask_area].to_numpy(), df['simulated_y'][mask_area].to_numpy()

        coord_x[i], coord_y[i] = x, y
        
        
        del mask_area, x, y

    df_chninf = pd.DataFrame(data = {'chnAreas': chnAreas, 'coord_x': coord_x, 'coord_y': coord_y})
        
    return  df_chninf


def calc_ciCOHs_rest(lfpdata):
    """
        calculate the ciCOHs for rest using ciCoherence_overtime

        Arg:
            lfpdata: nchns * ntemp * nsegs

        Return:
            ciCOH: nchns * nchns
    """

    nchns, _, nsegs = lfpdata.shape
    ciCOHs = np.zeros((nchns, nchns, nsegs))
    for segi in range(nsegs):
        
        if segi % 100 == 0:
            print("segi = " + str(segi) + "/" + str(nsegs))
        
        for chni in range(nchns -1):
            signal1 = lfpdata[chni, :, segi]
            
            for chnj in range(chni+1, nchns):
                signal2 = lfpdata[chnj, :, segi]
                
                # ciCOHs assignment
                ciCOHs[chni, chnj, segi] = ciCoherence_overtime(signal1, signal2)
                
                # symmetrical
                ciCOHs[chnj, chni, segi] = ciCOHs[chni, chnj, segi]
                
                del signal2
            del signal1
            
    ciCOH = np.mean(ciCOHs, axis = 2)

    return ciCOH

def fc_extract(savefilename):
    """
        dict fc extraction

        Return:
            fc
        Output:
            written fc to savefolder/savefilename 
    """

    ### lfpdata extract ###
    files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
    files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))
    files_moderate = glob.glob(os.path.join(inputfolder, '*_moderate_*'))


    lfpdata_normal, chnAreas, fs = lfp_extract(files_normal)
    lfpdata_mild, _, _ = lfp_extract(files_mild)
    lfpdata_moderate, _ , _= lfp_extract(files_moderate)



    ### balance mild, normal and moderate trials ###
    ntrials_normal, ntrials_mild = lfpdata_normal.shape[2], lfpdata_mild.shape[2]
    ntrials_moderate = lfpdata_moderate.shape[2]
    ntrials = min([ntrials_normal, ntrials_mild, ntrials_moderate])

    # balance trials by randomly selecting ntrials
    idx_ntrials = np.random.randint(ntrials_normal, size = ntrials)
    lfpdata_normal = lfpdata_normal[:,:,idx_ntrials]

    idx_ntrials = np.random.randint(ntrials_mild, size = ntrials)
    lfpdata_mild = lfpdata_mild[:,:,idx_ntrials]

    idx_ntrials = np.random.randint(ntrials_moderate, size = ntrials)
    lfpdata_moderate = lfpdata_moderate[:,:,idx_ntrials]


    ### calc ciCOH for each cond ###
    ciCOH_normal = calc_ciCOHs_rest(lfpdata_normal)
    ciCOH_mild = calc_ciCOHs_rest(lfpdata_mild)
    ciCOH_moderate = calc_ciCOHs_rest(lfpdata_moderate)


    ### dict fc generation ###
    ciCOH = dict()
    ciCOH['normal']  = ciCOH_normal
    ciCOH['mild'] = ciCOH_mild
    ciCOH['moderate'] = ciCOH_moderate


    setup = dict()
    setup['fs'] = fs
    setup['ntemp_normal'] = lfpdata_normal.shape[1]
    setup['ntrials_normal'] = lfpdata_normal.shape[2]
    setup['ntemp_mild'] = lfpdata_mild.shape[1]
    setup['ntrials_mild'] = lfpdata_mild.shape[2]
    setup['ntemp_moderate'] = lfpdata_moderate.shape[1]
    setup['ntrials_moderate'] = lfpdata_moderate.shape[2]



    fc = dict()
    fc['ciCOH'] = ciCOH
    fc['chnAreas'] = chnAreas
    fc['setup']= setup



    ### save ####
    with open(os.path.join(savefolder, savefilename + '.pickle'), 'wb') as fp:
        pickle.dump(fc, fp, protocol=pickle.HIGHEST_PROTOCOL)

    
    return fc



def fc_visual_save(fc, lowweight, savenamefile_prefix):
    """
        visual fc and save the figure

        Arg:
            fc: the dict fc
    """


    ### text setup for brain areas ###
    pos_text_lefttop1 = [-80, 50, 30]
    pos_text_righttop1 = [120, 50, 30]
    pos_text_lefttop2 = [-80, 70, 10]
    pos_text_Down1 = [0, 500, 15]
    pos_text_Down2 = [0, 550, 15]
    
    texts_org = dict()


    # plot
    df_chninf = assign_coord2chnArea(area_coord_file, fc['chnAreas'])
    print(df_chninf)
    for cond in fc['ciCOH'].keys():
        ciCOH = fc['ciCOH'][cond]
        ntrials, ntemp = fc['setup']['ntrials_' + cond], fc['setup']['ntemp_' + cond]


        texts = texts_org.copy()
        texts[cond] = pos_text_righttop1
        text_thred = 'thred = ' + str(np.round(lowweight, decimals = 2))
        texts[text_thred] = pos_text_lefttop2
        texts[text_task] = pos_text_lefttop1
        text_ntrials = 'ntrials = ' + str(ntrials)
        texts[text_ntrials] = pos_text_Down2
        texts[text_ntrials] = pos_text_Down2
        

        saveFCGraph = os.path.join(savefolder, savenamefile_prefix + '_lw' + str(np.round(lowweight, decimals = 2)) + '_'  + cond + '.png')

        igplot = ciCOH_visual_save(ciCOH = ciCOH, chnInf = df_chninf, lowweight = lowweight, 
                                savefile = saveFCGraph, texts = texts, threds_edge = None)

        del texts[cond], texts[text_ntrials]

    del df_chninf


def fc_select(fc, areas_used):
    """
        select fc only contains the areas in areas_used
        
        Arg:
            fc: the fc dict
            areas_used: list containing the areas used
    """



    chnAreas = fc['chnAreas'].copy()


    ### extract idx_del ###
    if 'GP' in areas_used:
        areas_used.remove('GP')
        areas_used = areas_used + ['gp0-1', 'gp1-2', 'gp2-3', 'gp3-4', 'gp4-5', 'gp5-6', 'gp6-7']

    if 'STN' in areas_used:
        areas_used.remove('STN')
        areas_used = areas_used + ['stn0-1', 'stn1-2', 'stn2-3', 'stn3-4', 'stn4-5', 'stn5-6', 'stn6-7']

    idx_del = []
    for i, area in enumerate(chnAreas):
        if area not in areas_used:
            idx_del.append(i)
    idx_del.reverse()

    ###  generate the used df_chninf ###
    for i in idx_del:
        del chnAreas[i]


    ### generate used ciCOH ###
    ciCOH = dict()
    for cond in fc['ciCOH'].keys():
        ciCOH_1cond = fc['ciCOH'][cond].copy()
        ciCOH_1cond = np.delete(ciCOH_1cond, idx_del, axis = 0)
        ciCOH_1cond = np.delete(ciCOH_1cond, idx_del, axis = 1)

        ciCOH[cond] = ciCOH_1cond


    fc_new = dict()
    fc_new['chnAreas'] = chnAreas
    fc_new['ciCOH'] = ciCOH
    fc_new['setup'] = fc['setup'].copy()


    return fc_new


def main():

    ### fc extract ###
    savename_fc = 'fc_rest' + '_freq' + str(freq[0]) + '_' + str(freq[1])
    if os.path.exists(os.path.join(savefolder, savename_fc + '.pickle')):
        print("reading the exist fc values")
        fc = pd.read_pickle(os.path.join(savefolder, savename_fc + '.pickle'))
    else:
        fc = fc_extract(savename_fc)




    # find lowweight
    pvals_vec = []
    ciCOH_vec = []
    for cond in fc['ciCOH'].keys():
        ciCOH = fc['ciCOH'][cond]
        ntrials, ntemp = fc['setup']['ntrials_' + cond], fc['setup']['ntemp_' + cond]

        pvals = pvals_fc_overtime(ciCOH = ciCOH, ntrials = ntrials, ntemp = ntemp, f = (freq[0] + freq[1])//2, t = ntemp/fc['setup']['fs'])
        
        nchns = pvals.shape[0]
        for chi in range(nchns-1):
            for chj in range(chi + 1, nchns):
                pvals_vec.append(pvals[chi, chj])
                ciCOH_vec.append(ciCOH[chi, chj])
        
    pvals_vec = np.asarray(pvals_vec)
    ciCOH_vec = np.asarray(ciCOH_vec)
    rejs, _, _, _ = multipletests(pvals_vec, alpha = 0.05, method = 'fdr_bh')
    lowweight = min(ciCOH_vec[rejs])
    print(lowweight)


    ### fc visual and save ##
    fcGraph_prefix = 'vFC_rest' + '_f' + str(freq[0]) + '_' + str(freq[1])


    # All
    fc_visual_save(fc,  lowweight, fcGraph_prefix + 'ALL')



    # right thalamus and SMA/M1
    fc_new = fc_select(fc, ['VA', 'VLo/VPLo', 'SMA', 'Sensory Motor'])
    fc_visual_save(fc_new, lowweight, fcGraph_prefix + 'CorRTha')
    del fc_new



    # right thalamus and GP
    fc_new = fc_select(fc, ['VA', 'VLo/VPLo', 'GP'])
    fc_visual_save(fc_new, lowweight, fcGraph_prefix + 'GPRTha')
    del fc_new





animal =  re.search('NHPs/[a-zA-Z]*/', __file__).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
corresfolder, correparentfolder = code_corresfolder(__file__)


freq = [13, 15]
inputfolder = os.path.join(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'Rest', 'm4_restData_avgAreafiltered11_13_noM1PMC')
area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')
savefolder = corresfolder
text_task = 'Rest'

if __name__ == '__main__':
    main()


