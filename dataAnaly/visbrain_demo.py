#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np

from visbrain.objects import ConnectObj, SceneObj, SourceObj, BrainObj
from visbrain.io import download_file


# In[2]:


# Download data
arch = np.load(download_file('phase_sync_delta.npz', astype='example_data'))
nodes, edges = arch['nodes'], arch['edges']


# In[12]:


# Create the scene with a black background
sc = SceneObj(size=(1500, 600))


# Coloring method
color_by = 'strength'
# Because we don't want to plot every connections, we only keep connections
# above .7
select = edges > .7
# Define the connectivity object
c_default = ConnectObj('default', nodes, edges, select=select, line_width=2.,
                       cmap='Spectral_r', color_by=color_by)
# Then, we define the sources
s_obj = SourceObj('sources', nodes, color='#ab4642', radius_min=15.)
sc.add_to_subplot(c_default, title='Color by connectivity strength')
# And add connect, source and brain objects to the scene
sc.add_to_subplot(s_obj)
sc.add_to_subplot(BrainObj('B3'), use_this_cam=True)


# In[13]:


# Finally, display the scene
sc.preview()

