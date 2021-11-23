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
codefolder = re.match('.*exp.*code', os.getcwd()).group()
sys.path.append(codefolder)
from util.folder_extract import exp_subfolders, code_corresfolder
from connAnalyTool.synchronization_indices import ciCoherence_overtime
from connAnalyTool.fc_visual_time import threshold_fc_overtime, ciCOH_visual_save, pvals_fc_overtime, weight_visual_save
from mne.stats import fdr_correction

import networkx as nx



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
            print("calc ciCOH segi = " + str(segi) + "/" + str(nsegs))
        
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

def dailyfc_extract(lfpdata_1file):
    """
        dict fc extraction

        Return:
            fc
        Output:
            
    """

    ### lfpdata extract ###


    lfpdata, chnAreas, fs = lfp_extract(lfpdata_1file)


    ### calc ciCOH for each cond ###
    ciCOH = calc_ciCOHs_rest(lfpdata)


    ### dict fc generation ###

    setup = dict()
    setup['fs'] = fs
    setup['ntemp'] = lfpdata.shape[1]
    setup['ntrials'] = lfpdata.shape[2]


    fc = dict()
    fc['ciCOH'] = ciCOH
    fc['chnAreas'] = chnAreas
    fc['setup']= setup

    
    return fc



def fc_visual_save(fc, lowweight, savenamefile_prefix):
    """
        visual fc and save the figure

        Arg:
            fc: the dict fc
    """


    ### text setup for brain areas ###
    pos_text_lefttop1 = [-80, 50, 30]
    pos_text_middletop1 = [120, 50, 30]
    pos_text_lefttop2 = [-80, 70, 10]
    pos_text_leftDown1 = [-80, 550, 30]
    pos_text_leftDown2 = [-80, 570, 10]
    pos_text_leftDown3 = [-80, 580, 10]
    
    texts_org = dict()

    lowweight = np.round(lowweight, decimals = 2) 

    # plot
    df_chninf = assign_coord2chnArea(area_coord_file, fc['chnAreas'])
    for ci, cond in enumerate(fc['ciCOH'].keys()):
        ciCOH = fc['ciCOH'][cond]
        ntrials, ntemp = fc['setup']['ntrials_' + cond], fc['setup']['ntemp_' + cond]


        texts = texts_org.copy()
        
        text_thred = 'thred = ' + str(np.round(lowweight, decimals = 2))
        text_ntrials = 'ntrials = ' + str(ntrials)

        texts[cond] = pos_text_middletop1
        texts[text_task] = pos_text_leftDown1
        texts[text_ntrials] = pos_text_leftDown2
        texts[text_thred] = pos_text_leftDown3
        

        saveFCGraph = os.path.join(savefolder, savenamefile_prefix + '_lw' + str(np.round(lowweight, decimals = 2)) + '_'  + cond + '.png')

        igplot = ciCOH_visual_save(ciCOH = ciCOH, chnInf = df_chninf, lowweight = lowweight, 
                                savefile = saveFCGraph, texts = texts, threds_edge = None)

        del texts[cond], texts[text_ntrials]

        img = cv2.imread(saveFCGraph)
        if ci == 0:
            imgs = img
        else:
            imgs = np.concatenate((imgs, np.zeros((img.shape[0], 5, 3)),img), axis = 1)

        os.remove(saveFCGraph)

    # combine all conditions
    print(imgs.shape)
    saveFCGraph_comb = os.path.join(savefolder, 'comb_' + savenamefile_prefix + '_lw' + str(np.round(lowweight, decimals = 2))  + '.png')
    cv2.imwrite(saveFCGraph_comb, imgs)

def ciCOH_select(ciCOH, chnAreas, areas_used):
    """
        select ciCOH only contains the areas in areas_used
        
        Arg:
            ciCOH: the ciCOH nchns * nchns
            chnAreas: area list
            areas_used: list containing the areas used
    """

    chnAreas_new = chnAreas.copy()
    ciCOH_new = ciCOH.copy()

    ### extract idx_del ###
    if 'GP' in areas_used:
        areas_used.remove('GP')
        areas_used = areas_used + [area for area in chnAreas_new if 'gp' in area]

    if 'STN' in areas_used:
        areas_used.remove('STN')
        areas_used = areas_used + [area for area in chnAreas_new if 'stn' in area]

    idx_del = []
    for i, area in enumerate(chnAreas_new):
        if area not in areas_used:
            idx_del.append(i)
    idx_del.reverse()

    ###  generate the used df_chninf ###
    for i in idx_del:
        del chnAreas_new[i]


    ### generate used ciCOH ###
    ciCOH_new = np.delete(ciCOH_new, idx_del, axis = 0)
    ciCOH_new = np.delete(ciCOH_new, idx_del, axis = 1)


    return ciCOH_new, chnAreas_new

def comb_fc(filepatt):
    """
        combine all fc figures belong to same 
    """


    files = glob.glob(os.path.join(savefolder, filepatt))
    print(filepatt)
    
    if files == []:
        imgs = []
        print('No files found for ' + filepatt)
        return


    imgs = np.empty((600, 600, 3))
    for fi, file in enumerate(files):
        img = cv2.imread(file)
        
        if fi == 0:
            imgs = img
        else:
            imgs = np.concatenate((imgs, img), axis = 2)

    idx = filepatt.find('freq')
    comb_fcGraph = os.path.join(savefolder, 'comb_' + filepatt[idx: -len('.mat')])
    cv2.imwrite(comb_fcGraph, imgs)
    print(comb_fcGraph)

    # find lowweight
    pvals_vec = []
    ciCOH_vec = []
    pvals = pvals_fc_overtime(ciCOH = ciCOH, ntrials = ntrials, ntemp = ntemp, f = f, t = t)

    reject, pval_corr = fdr_correction(pvals, alpha=0.05, method='indep')
    

    
    lowweight = min(ciCOH_vec[rejs])

    return lowweight



def graph_metrics(weight):
    nchns = weight.shape[0]
    G = nx.Graph()
    G.add_nodes_from(np.arange(1, nchns + 1))
    for i in range(0, nchns - 1):
        for j in range(i+1, nchns):
            if weight[i,j] > 0:
                G.add_edge(i, j, weight= weight[i, j])

    avg_CC = nx.average_clustering(G)

    del G
    return avg_CC
    
def dailyfc_visual(files):

    ### fc extract ###
    for onefile in files:
        
        filename = os.path.basename(onefile)
        datestr = re.search('[0-9]{8}', filename).group()

        lfpdata, chnAreas, fs = lfp_extract([onefile])

        if 'cond' not in locals():
            cond = re.search('_[a-z]*_[0-9]{8}', filename).group()[1:-9]
        
        if 'lfpdatas' not in locals():
            lfpdatas = lfpdata
            datestrs = datestr
        else:
            lfpdatas = np.concatenate((lfpdatas, lfpdata), axis = 2)
            datestrs = datestrs + '_' + datestr 


        del lfpdata, datestr
        

        ### if enough lfpdatas
        if lfpdatas.shape[2] >= 500 / 5:

            lfp1, lfp2 = lfpdatas[:, 0:500, :], lfpdatas[:, 125:625, :]
            lfp3, lfp4 = lfpdatas[:, 250:750, :], lfpdatas[:, 375:875, :]
            lfp5= lfpdatas[:, 500:, :]
            lfpdatas = np.concatenate((lfp1, lfp2, lfp3, lfp4, lfp5), axis=2)

            idx_ntrials = np.random.randint(lfpdatas.shape[2], size = 500)
            lfpdatas = lfpdatas[:, :, idx_ntrials]
            nchns, ntemp, ntrials = lfpdatas.shape
            
            
            ### calc ciCOH for each cond ###
            ciCOH = calc_ciCOHs_rest(lfpdatas)
            ciCOH = abs(ciCOH)

            
            ### all ##
            save_prefix = 'all'
            
            # get weight matrix
            pvals = pvals_fc_overtime(ciCOH = ciCOH, ntrials = ntrials, ntemp = ntemp, f = (freq[0] + freq[1])/2, t = ntemp/fs)
            reject, pval_corr = fdr_correction(pvals, alpha = 0.1, method='indep')
            [rows, cols]= np.where(reject)
            weight = np.zeros(ciCOH.shape)
            if len(rows) > 0:
                weight[rows, cols] = ciCOH[rows, cols]

            # visual and save
            saveFCname = cond + '_'  + save_prefix + '_' + datestrs + '.png'
            saveFCGraph = os.path.join(savefolder, saveFCname)
            weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file, chnAreas), 
                                savefile = saveFCGraph, texts = None, threds_edge = None)


            # network metric
            avg_CC = graph_metrics(weight)
            d = {saveFCname: avg_CC}
            with open(os.path.join(savefolder, 'avgCC.csv'), 'a+') as f:
                for key in d.keys():
                    f.write("%s,%s\n"%(key,d[key]))
            
            del avg_CC, d
            del pvals, reject, pval_corr, rows, cols
            del saveFCGraph, weight, save_prefix, saveFCname

            
            
            ### left thalamus and SMA/M1 ###
            save_prefix = 'leftThaCor_' 
            areas_used = ['lVA', 'lVLo/VPLo', 'lSMA', 'rSMA','M1']

            # subareas selection
            ciCOH_new, chnAreas_new = ciCOH_select(ciCOH, chnAreas, areas_used)
            
            
            # get weight matrix
            pvals = pvals_fc_overtime(ciCOH = ciCOH_new, ntrials = ntrials, ntemp = ntemp, f = (freq[0] + freq[1])/2, t = ntemp/fs)
            reject, pval_corr = fdr_correction(pvals, alpha = 0.1, method='indep')
            [rows, cols]= np.where(reject)
            weight = np.zeros(ciCOH.shape)
            if len(rows) > 0:
                weight[rows, cols] = ciCOH[rows, cols]

            # visual and save
            saveFCGraph = os.path.join(savefolder, cond + '_' + save_prefix + '_' + datestrs + '.png')
            weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file, chnAreas_new), 
                                savefile = saveFCGraph, texts = None, threds_edge = None)
            del ciCOH_new, chnAreas_new, save_prefix, areas_used
            del saveFCGraph, weight




            ### right thalamus and SMA/M1 ###
            save_prefix = 'rightThaCor'
            areas_used = ['rVA', 'rVLo/VPLo', 'lSMA', 'rSMA','M1']
            
            # subareas selection
            ciCOH_new, chnAreas_new = ciCOH_select(ciCOH, chnAreas, areas_used)

            # get weight matrix
            pvals = pvals_fc_overtime(ciCOH = ciCOH_new, ntrials = ntrials, ntemp = ntemp, f = (freq[0] + freq[1])/2, t = ntemp/fs)
            reject, pval_corr = fdr_correction(pvals, alpha = 0.1, method='indep')
            [rows, cols]= np.where(reject)
            weight = np.zeros(ciCOH.shape)
            if len(rows) > 0:
                weight[rows, cols] = ciCOH[rows, cols]

            # visual and save
            saveFCGraph = os.path.join(savefolder, cond + '_' + save_prefix + '_' + datestrs + '.png')
            weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file, chnAreas_new), 
                                savefile = saveFCGraph, texts = None, threds_edge = None)
            del ciCOH_new, chnAreas_new, save_prefix, areas_used
            del saveFCGraph, weight
           

            
            ### right thalamus and GP ###
            save_prefix = 'gpRightTha'
            areas_used = ['rVA', 'rVLo/VPLo', 'GP']
            
            # subareas selection
            ciCOH_new, chnAreas_new = ciCOH_select(ciCOH, chnAreas, areas_used)

            # get weight matrix
            pvals = pvals_fc_overtime(ciCOH = ciCOH_new, ntrials = ntrials, ntemp = ntemp, f = (freq[0] + freq[1])/2, t = ntemp/fs)
            reject, pval_corr = fdr_correction(pvals, alpha = 0.1, method='indep')
            [rows, cols]= np.where(reject)
            weight = np.zeros(ciCOH.shape)
            if len(rows) > 0:
                weight[rows, cols] = ciCOH[rows, cols]

            # visual and save
            saveFCGraph = os.path.join(savefolder, cond + '_' + save_prefix + '_' + datestrs + '.png')
            weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file, chnAreas_new), 
                                savefile = saveFCGraph, texts = None, threds_edge = None)
            del ciCOH_new, chnAreas_new, save_prefix, areas_used
            del saveFCGraph, weight



            ### left thalamus and GP ###
            save_prefix = 'gpLeftTha'
            areas_used = ['lVA', 'lVLo/VPLo', 'GP']
            
            # subareas selection
            ciCOH_new, chnAreas_new = ciCOH_select(ciCOH, chnAreas, areas_used)

            # get weight matrix
            pvals = pvals_fc_overtime(ciCOH = ciCOH_new, ntrials = ntrials, ntemp = ntemp, f = (freq[0] + freq[1])/2, t = ntemp/fs)
            reject, pval_corr = fdr_correction(pvals, alpha = 0.1, method='indep')
            [rows, cols]= np.where(reject)
            weight = np.zeros(ciCOH.shape)
            if len(rows) > 0:
                weight[rows, cols] = ciCOH[rows, cols]

            # visual and save
            saveFCGraph = os.path.join(savefolder, cond + '_' + save_prefix + '_' + datestrs + '.png')
            weight_visual_save(weight, chnInf = assign_coord2chnArea(area_coord_file, chnAreas_new), 
                                savefile = saveFCGraph, texts = None, threds_edge = None)
            del ciCOH_new, chnAreas_new, save_prefix, areas_used
            del saveFCGraph, weight



            del lfpdatas, idx_ntrials, datestrs
            del ciCOH
            




def generate_video(genvideofile, images): 
    """
    
    """
      

    frame = cv2.imread( images[0])
  
    # setting the frame width, height width 
    # the width, height of first image 
    height, width, layers = frame.shape   
  
    video = cv2.VideoWriter(genvideofile, 0, 1, (width, height))  
  
    # Appending the images to the video one by one 
    for image in images:  
        video.write(cv2.imread(image))  
      
    # Deallocating memories taken for window creation 
    cv2.destroyAllWindows()  
    video.release()  # releasing the video generated 
    print("Generated video")
            


def main():

    files_normal = glob.glob(os.path.join(inputfolder, '*_normal_*'))
    files_mild = glob.glob(os.path.join(inputfolder, '*_mild_*'))
    files_moderate = glob.glob(os.path.join(inputfolder, '*_moderate_*'))

    files_normal.sort()
    dailyfc_visual(files_normal)

    files_mild.sort()
    dailyfc_visual(files_mild)


    files_moderate.sort()
    dailyfc_visual(files_moderate)




animal =  re.search('NHPs/[a-zA-Z]*/', __file__).group()[len('NHPs/'):-1]
_, _, pipelinefolder, _= exp_subfolders()
corresfolder, correparentfolder = code_corresfolder(__file__)


freq = [26, 28]
inputfolder = os.path.join(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'Rest', 'm4_restData_filtered' + str(freq[0]) + '_' + str(freq[1]) + '_eqLen')
area_coord_file = os.path.join(correparentfolder, 'chn_brainArea_simCoord_BrainArea.csv')
savefolder = corresfolder
text_task = 'Rest'

if __name__ == '__main__':
    # main()

    cond = 'normal'
    save_prefix = 'all'
    images =  glob.glob(os.path.join(savefolder, cond + '_' + save_prefix + '_[0-9]*'))
    images.sort()
    genvideofile = os.path.join(savefolder, cond + '_' + save_prefix + '.avi')
    generate_video(genvideofile, images)


    cond = 'mild'
    save_prefix = 'all'
    images =  glob.glob(os.path.join(savefolder, cond + '_' + save_prefix + '_[0-9]*'))
    images.sort()
    genvideofile = os.path.join(savefolder, cond + '_' + save_prefix + '.avi')
    generate_video(genvideofile, images)


    cond = 'moderate'
    save_prefix = 'all'
    images =  glob.glob(os.path.join(savefolder, cond + '_' + save_prefix + '_[0-9]*'))
    images.sort()
    genvideofile = os.path.join(savefolder, cond + '_' + save_prefix + '.avi')
    generate_video(genvideofile, images)

