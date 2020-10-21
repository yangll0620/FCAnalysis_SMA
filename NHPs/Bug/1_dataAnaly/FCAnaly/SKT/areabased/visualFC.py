import os, sys
import numpy as np
import pandas as pd
import pickle
from igraph.drawing.text import TextDrawer
import cairo

import re
codefolder = re.match('.*exp.*code', __file__).group()
sys.path.append(codefolder)
from connAnalyTool import graphVisual


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
    weight[np.where(weight<lowweight)] =0

    
    # chnInf_coord
    x, y = chnInf['coord_x'].to_numpy(), chnInf['coord_y'].to_numpy()
    x, y = np.expand_dims(x, axis = 1), np.expand_dims(y, axis = 1)
    chn_coord = np.concatenate((x, y), axis = 1)


    # chnAreas format: gp1-2, stn1-2
    vsarea = chnInf['chnAreas'].to_list()
    vsname = chnInf['chnAreas'].to_list()
    if 'stn1-2' in vsarea:
        for i, barea in enumerate(vsarea):

            if re.match('stn[0-6]-[1-7]', barea):
                vsarea[i] = 'STN'
                vsname[i] = barea[3:]
            
            if re.match('gp[0-6]-[1-7]', barea):
                vsarea[i] = 'GP'
                vsname[i] = barea[2:]

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
    vertexcolors_opt = ['red', 'blue', 'green' ,'yellow', 'purple', 'gray', 'magenta', 'cyan', 'black']
    vcolors = ['black'] * len(graph.vs['area'])
    uiqareas = list(set(graph.vs["area"]))
    uiqareas.sort()

    print(uiqareas)
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