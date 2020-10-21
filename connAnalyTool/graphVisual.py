#!/usr/bin/env python
# coding: utf-8

# # Functionality:  visualization of graph

# ## Import Packages

# In[1]:


import igraph as ig
import igraph.drawing.colors as igcolors
from igraph.drawing.text import TextDrawer
import cairo

import numpy as np


# ## Random Matrix

# In[3]:


def randweight_sysmetric(nchns=64):
    w = np.random.rand(nchns, nchns)
    
    weight = (w + w.T)/2
    
    return weight


# ## Create Graph

# ### generate a graph from its adjacency matrix

# In[4]:


def graph_create(weight, tagdirected = False, vsarea = None, vscoord = None, vsname = None, edge_color = None):
    """
        Create graph based on weight matrix
        
        Args:
            
            weight: weight matrix (should be symmetrics if indirected graph)
            
            tagdirected: tag for if directed (True) or indirected (False, default)
            
            vsname: a list of string representing vertex name
            
            vschni: chni number for each channel
            
        return:
            graph: the generated graph
            
    """

    # change nan to 0
    ijs_nan = np.where(np.isnan(weight))

    for k in range(len(ijs_nan[0])):
        i_nan, j_nan = ijs_nan[0][k], ijs_nan[1][k]
        weight[i_nan][j_nan] = 0
    
    #  adjacency matrix
    adjmatrix = (weight>0).tolist()
    

    # create graph from adjacency matrix
    if not tagdirected:
        """ undirected graph """

        # check the symmetric of weights matrix
        if not np.allclose(weight, weight.T, rtol=1e-05, atol=1e-08):
            print("weight matrix is not symmetric")
            return []
        
        
        # generate a undirected graph from its adjacency matrix (weights>0)
        graph = ig.Graph.Adjacency(adjmatrix, "UNDIRECTED")
        


        # weights for the undirected graph, lower triangular 
        wght = np.triu(weight)

    else:
        """ directed graph """

        # generate a undirected graph from its adjacency matrix (weights>0)
        graph = ig.Graph.Adjacency(adjmatrix, "DIRECTED")

        # weights for the directed graph
        wght = weight
        
    
    
    #### vertex 
    
    # set the vs attributes area
    if vsarea is not None:
        graph.vs["area"] = vsarea

    
    # set the vs position
    if vscoord is not None:
        graph.vs["coord"] = vscoord
        
    
    # set the vertex color
    if vsname is not None:
        graph.vs['name'] = vsname
    
    
    
    
    ###e edge
    
    # set the attribute of weight for graph element es
    graph.es["weight"] = wght[wght.nonzero()]
    
    
    # set the vertex color
    if vsname is not None:
        graph.es['color'] = edge_color
    
    return graph


# ## Graph Visulization

# ### Set visulization style

# In[17]:


def graph_style(graph, visual_style = None):
    """
    
        args:
            graph: graph
            
        return:
            visual_style: a dictionary containing the visualization style
    """

    if visual_style is None:
        visual_style = dict()
        print("visual_style is None")

    # vertex color
#     if not 'vertex_color' in visual_style.keys():
    
#         try:
#             "set vertex  color based on vs attribute area"

#             # unique areas
#             uniqarea = list(set(graph.vs["area"]))

#             # set palette
#             palette = igcolors.ClusterColoringPalette(len(uniqarea))

#             # set vertex color
#             visual_style['vertex_color'] = [palette[uniqarea.index(area)] for area in graph.vs["area"]]

#         except KeyError:
#             visual_style['vertex_color'] =  'black'
        
    
    # vertex size
    outdegree = graph.outdegree()
    if max(outdegree) > 0:
        visual_style['vertex_size'] = [x/max(outdegree)*10+5 for x in outdegree]
    
    # vertex label size
    visual_style['vertex_label_size'] = 2

    # vertex label distance
    visual_style['vertex_label_dist'] = 2

    # vertex label color
    visual_style['vertex_label_color'] = 'black'
    
    # vertex label name
    try:
        visual_style["vertex_label"] = graph.vs['name']
        
    except KeyError:
        print("No vertex name")
    
    # vertext label size
    visual_style["vertex_label_size"] = 10
    
    visual_style["vertex_label_angle"] = 0;

    # the outdegree for each vertex
    outdegree = graph.outdegree()
    
    
    visual_style["margin"] = [10,10,30,10];
    

    # set edge width if not exist
    if not 'edge_width' in visual_style.keys():
        visual_style['edge_width'] = graph.es['weight']
        
    
    # set edge color if not exist
    if not 'edge_color' in visual_style.keys():
        print("Use default edge color")

        

    # layout
    try:
        "set vertex  layout based on vs attribute coordinate coord"
        visual_style['layout'] = ig.Layout(graph.vs["coord"])
    
    except KeyError:
        print("Does not have a user defined layout, use the default.")
    
    return visual_style


# ### Actual Plot

# In[6]:


def graph_plot(graph, visual_style, texts = None):
    """
        @param graph:
        
        @param visual_style:
        
        @param texts: a dictionary storing the text to be printed with text[key]: [x,y, fontsize]
        e.g texts["M1"] = [0, 180, 30]
        
        @return igplot: an igraph.drawing.Plot object
    """
    
    # an igraph.drawing.Plot object
    igplot = ig.plot(graph, **visual_style)
    
    # plot text
    if texts is not None:
        
        igplot.redraw()
            
        # Context object
        ctx = cairo.Context(igplot.surface)
        
        
        # draw each text in ig.plot
        for text in texts:
            val = texts[text]
            
            # fontsize
            try:
                fontsize = val[2]
            
            except IndexError:
                fontsize = 20
            
            
            # set the font size
            ctx.set_font_size(fontsize)
            
            # TextDrawer Object
            drawer = TextDrawer(ctx, text, halign=TextDrawer.CENTER)
            
            # draw the text at the coordinates coord
            drawer.draw_at(x=val[0], y=val[1], width=300)
            
    return igplot


# ## Test Section

# In[18]:


if False:
    weight = randweight_sysmetric(nchns=4)
    print(weight)

    graph = graph_create(weight)


    colors = ['black'] * len(graph.es['weight'])
    for i, w in enumerate(graph.es['weight']):
        if w>0.5:
            colors[i] = 'red'

    print(colors)
    visual_style = dict()     
    visual_style['edge_color'] = colors

    visual_style = graph_style(graph, visual_style)
    igplot = graph_plot(graph, visual_style)
    igplot.show()

