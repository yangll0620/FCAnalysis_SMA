import pandas as pd
import os
import numpy as np
import sys
from scipy.stats import norm


import re
codefolder = re.match('.*exp.*code', __file__).group()
sys.path.append(codefolder)

from simulated.python.thred_ciCOH_time import threshold_ciCOH_sin, corr_threshold_ciCOH_sin_BH
from util.folder_extract import exp_subfolders
from connAnalyTool import graphVisual

_, _, pipelinefolder, _= exp_subfolders()
savefile_threshold = os.path.join(pipelinefolder, 'fc_thred_time.pickle')


def threshold_fc_overtime(ciCOH, ntrials, ntemp, f, t):
    """
        identify threshold for fc

        Args:
            ciCOH: nchns * nchns

        Returns:
            threshold:

            corr_threshold:

        Output:
            save the threshold into savefile_threshold
    """

    df = pd.DataFrame(columns = ['ntemp', 'ntrials', 'f','threshold', 'mu', 'std'])

    record_exist = False
    if(os.path.exists(savefile_threshold)): # file_threshold exists
        
        df = pd.read_pickle(savefile_threshold)
        
        # check if threshold under the ntemp and f exist 
        mask = (df['ntemp'] == ntemp) & (df['ntrials'] == ntrials) & (df['f'] == f)
           
        if(df.loc[mask].shape[0] == 1): # record exist
            
            record_exist = True
            print("read threshold from file")

            
            df1 = df.loc[mask]

            threshold = df1['threshold'].to_numpy()
            mu, std = df1['mu'].to_numpy(), df1['std'].to_numpy()

            del df1

    if not record_exist: # run threshold_ciCOH_sin if record not exist
        threshold, mu, std = threshold_ciCOH_sin(ntimes = 500, ntrials = ntrials, ntemp = ntemp, 
                                                 f = f, t = t,ploton= False)


        # store the new record
        thred_sets = dict()
        thred_sets['ntemp'],  thred_sets['ntrials'], thred_sets['f'] = ntemp, ntrials, f
        thred_sets['mu'], thred_sets['std'] = mu, std
        thred_sets['threshold'] = threshold

        df = df.append(thred_sets, ignore_index = True)

        # write to file_threshold
        df.to_pickle(savefile_threshold)




    nchns = ciCOH.shape[0]
    ciCOHs_vec = list()
    for chni in range(nchns -1):
        for chnj in range(chni+1, nchns):
            ciCOHs_vec.append(ciCOH[chni][chnj])
    ciCOHs_vec = np.asarray(ciCOHs_vec)

    corr_threshold, mu, std = corr_threshold_ciCOH_sin_BH(ciCOHs_actual = ciCOHs_vec, ntimes = 500, 
                                ntrials = ntrials, ntemp = ntemp, f = f, t = t, false_rate = 0.01, mu = mu, std = std)
    

    threshold = np.around(threshold, decimals=2)
    corr_threshold = np.around(corr_threshold, decimals=2)
    
    
    return threshold, corr_threshold



def pvals_fc_overtime(ciCOH, ntrials, ntemp, f, t):
    """
        return pvalus for each value of ciCOH

        Args:
            ciCOH: nchns * nchns

        Returns:
            pvals:

        Output:
            save the threshold into savefile_threshold
    """

    df = pd.DataFrame(columns = ['ntemp', 'ntrials', 'f','threshold', 'mu', 'std'])

    record_exist = False
    if(os.path.exists(savefile_threshold)): # file_threshold exists
        
        df = pd.read_pickle(savefile_threshold)
        
        # check if threshold under the ntemp and f exist 
        mask = (df['ntemp'] == ntemp) & (df['ntrials'] == ntrials) & (df['f'] == f)
           
        if(df.loc[mask].shape[0] == 1): # record exist

            print("read threshold from file")
            
            record_exist = True
            
            df1 = df.loc[mask]

            threshold = df1['threshold'].to_numpy()
            mu, std = df1['mu'].to_numpy(), df1['std'].to_numpy()

            del df1

    if not record_exist: # run threshold_ciCOH_sin if record not exist
        threshold, mu, std = threshold_ciCOH_sin(ntimes = 500, ntrials = ntrials, ntemp = ntemp, 
                                                 f = f, t = t,ploton= False)


        # store the new record
        thred_sets = dict()
        thred_sets['ntemp'],  thred_sets['ntrials'], thred_sets['f'] = ntemp, ntrials, f
        thred_sets['mu'], thred_sets['std'] = mu, std
        thred_sets['threshold'] = threshold

        df = df.append(thred_sets, ignore_index = True)

        # write to file_threshold
        df.to_pickle(savefile_threshold)


    pvals = norm.sf(ciCOH, loc = mu, scale = std)
    
    return pvals

def ciCOH_visual_save(ciCOH, chnInf, lowweight, savefile, texts = None, threds_edge = None):
    """
        
        Args:
            
            ciCOH (np.sdarray): ciCOH matrix (nchns, nchns)
            
            chnInf (dict): dictionary for channels, ('chnsAreas', coord_x, coord_y)
            
            
            lowweight: the threshold lowweight, only weight>lowweight is treated as connection
            
            savefile: file to save the visualized figure
            
            texts:

            threds_edge: threds for edges
            
        Output:
            the visualizaton of ciCOH is saved in tobesavedfile
            
    """
    
    weight = abs(ciCOH)

    # weight > lowweight
    weight[np.where(weight<lowweight)] = 0


    igplot = weight_visual_save(weight, chnInf, savefile, texts = None, threds_edge = None)

    return igplot

def weight_visual_save(weight, chnInf, savefile, texts = None, threds_edge = None):
    """
        
        Args:
            
            weight (np.sdarray):  weight matrix (nchns, nchns)
            
            chnInf (dict): dictionary for channels, ('chnsAreas', coord_x, coord_y)
            
            
            savefile: file to save the visualized figure
            
            texts:

            threds_edge: threds for edges
            
        Output:
            the visualizaton of weight is saved in tobesavedfile
            
    """

    weight = abs(weight)

    # chnInf_coord
    x, y = chnInf['coord_x'].to_numpy(), chnInf['coord_y'].to_numpy()
    x, y = np.expand_dims(x, axis = 1), np.expand_dims(y, axis = 1)
    chn_coord = np.concatenate((x, y), axis = 1)


    # chnAreas format: gp1-2, stn1-2
    vsarea = chnInf['chnAreas'].to_list()
    vsname = chnInf['chnAreas'].to_list()
    if 'stn1-2' in vsarea or  'gp1-2' in vsarea:
        for i, barea in enumerate(vsarea):

            if re.match('stn[0-6]-[1-7]', barea):
                vsarea[i] = 'STN'
                vsname[i] = barea
            
            if re.match('gp[0-6]-[1-7]', barea):
                vsarea[i] = 'GP'
                vsname[i] = barea

    # create new graph
    graph = graphVisual.graph_create(weight, vsarea = vsarea, vscoord = chn_coord, vsname = vsname)

    
    ### set graph visualization style ###
    visual_style = dict()

    # set the edge color base one thred
    if threds_edge is not None:
        colors_opt = ['gray', 'blue' ,'green','red']
        
        colors = ['black'] * len(graph.es['weight'])
        for i, w in enumerate(graph.es['weight']):
            if w > threds_edge[2]: 
                colors[i] = colors_opt[3]
            elif w > threds_edge[1]:
                colors[i] = colors_opt[2]
            elif w > threds_edge[0]:
                colors[i] = colors_opt[1]           
        visual_style['edge_color'] = colors


    # set the vertex color
    vertexcolors_opt = ['red', 'blue', 'green' ,'yellow', 'purple', 'gray', 'magenta', 'cyan', 'black', 'Azure', 'violet']
    vcolors = ['black'] * len(graph.vs['area'])
    uiqareas = list(set(graph.vs["area"]))
    uiqareas.sort()

    print(uiqareas)
    print('uniqareas = ' + str(len(uiqareas)))
    for i, uiarea in enumerate(uiqareas):

        # set the chn belong to the same area to be the same color 
        for vi, area in enumerate(graph.vs["area"]):
            if uiarea == area:
                vcolors[vi] = vertexcolors_opt[i] 
    visual_style['vertex_color'] = vcolors


    # apply the set visual style  
    visual_style = graphVisual.graph_style(graph, visual_style)
    

     # plot graph
    igplot = graphVisual.graph_plot(graph, visual_style, texts = texts)
    
    # save graph
    igplot.save(savefile)
    print("Figure saved to " + savefile)
    
    return igplot