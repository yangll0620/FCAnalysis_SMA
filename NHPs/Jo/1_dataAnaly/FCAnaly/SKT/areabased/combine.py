import os, sys
import numpy as np
import cv2 



imgs = np.empty(shape=[600,0,3])

cond = 'normal'
lowweight = 0.07
file_fc = os.path.join(savefolder, 
                       savefile_fcgraph_prefix + '_lowweight' + str(lowweight) + '_' + cond + '.png') 
img = cv2.imread(file_fc)
imgs = np.concatenate((imgs, img), axis = 1)

cond = 'mild'
lowweight = 0.09
file_fc = os.path.join(savefolder, 
                       savefile_fcgraph_prefix + '_lowweight' + str(lowweight) + '_' + cond + '.png') 
img = cv2.imread(file_fc)
imgs = np.concatenate((imgs, img), axis = 1)

cond = 'moderate'
lowweight = 0.09
file_fc = os.path.join(savefolder, 
                       savefile_fcgraph_prefix + '_lowweight' + str(lowweight)+ '_' + cond + '.png') 
img = cv2.imread(file_fc)
imgs = np.concatenate((imgs, img), axis = 1)



# write the combined img back
cv2.imwrite(os.path.join(savefolder, 'combined_' + phase + '_' + strfreq  + '.png'), imgs)