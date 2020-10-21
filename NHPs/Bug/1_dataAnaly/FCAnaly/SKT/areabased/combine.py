import os, sys
import numpy as np
import cv2 
import glob
import re

savefolder = '/home/lingling/Insync/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/pipeline/NHPs/Bug/1_dataAnaly/FCAnaly/SKT/areabased/m1_SKT'
phases = ['base', 'reach', 'return']
conds = ['normal', 'mild']

# find all freqstr
files = glob.glob(os.path.join(savefolder, 'ciCOH_*_reach_*.png'))
freqstrs = list()
for f in files:
    fname = os.path.split(f)[1] 
    
    freqstr = re.search('freq[0-9]*_[0-9]*', fname).group()
    freqstrs.append(freqstr)

ufreqstr = set(freqstrs)

print(ufreqstr)
 
img_white = np.zeros((600, 2, 3)) 
for freqstr in ufreqstr:

    for ci, cond in enumerate(conds):
        for pi, phase in enumerate(phases):
            files = glob.glob(os.path.join(savefolder, '*_' + phase + '_*_' + cond + '.png'))
            file_fc = files[0]
            img = cv2.imread(file_fc)

            if pi == 0:
                img_1cond = img
                
            else:
                img_1cond = np.concatenate((img_1cond, img_white, img), axis = 1)

            del files, file_fc, img

        if ci == 0:
            imgs = img_1cond
            
        else:
            img_black = np.zeros((2, img_1cond.shape[1], 3))
            imgs = np.concatenate((imgs, img_black, img_1cond), axis = 0)


    cv2.imwrite(os.path.join(savefolder, 'combined_' + freqstr +  '.png'), imgs)
    del imgs
