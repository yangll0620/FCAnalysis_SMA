# # folder extraction return absoluate folder path
# 
# * exp_subfolders(): return datafolder, codefolder, pipelinefolder, outputfolder
# 
# * code_corresfolder(): return the corresponding folder for codefolder/codefile,  create the corresponding folder and/or the'out' and 'store' subfolders


import os


def exp_subfolders(): 
    """ 
        exp_subfolders return datafolder, codefolder, pipelinefolder, outputfolder
        
    Usage:
        
        datafolder, codefolder, pipelinefolder, outputfolder = exp_subfolders()

    Returns:
        datafolder: data folder 
        
        codefolder: code folder 
        
        pipelinefolder: pipeline folder
        
        outputfolder: output folder

    """
    
    currfolder = os.getcwd()
    
    # the exp folder path
    expfolder = currfolder[0 : currfolder.find('exp')+len('exp')]
    
    # data folder
    datafolder = os.path.join(expfolder, 'data')
    
    # code folder
    codefolder = os.path.join(expfolder, 'code')
    
    # pipeline folder
    pipelinefolder = os.path.join(expfolder, 'pipeline')
    
    # output folder
    outputfolder = os.path.join(expfolder, 'output')
    
    
    return datafolder, codefolder, pipelinefolder, outputfolder





def code_corresfolder(codefolder, codefile, makefolder = True, makesubfolder = False):
    """
        code_corresfolder 
            return the corresponding folder for codefolder/codefile
            
            create the corresponding folder if not exist (makefolder = True)
                   ,and the out and store subfolders of the corresponding folder (makesubfolder = True) 
            
        args:
        
            codefolder: the folder for the codefile
            
            codefile: code file name without suffix (i.e codefile = examplecode for examplecode.ipynb)
            
            makefolder: tag for creating codefolder (default True)
            
            makesubfolder: tag for creating the 'out' and 'store' subfolders of codefolder (default False)
            
        
        return:
            codecorresfolder: the corresponding folder for codefolder/codefile
            
            codecorresparentfolder: the corresponding parent folder for codefolder/codefile
    """
    

    # extract the substring of codefolder after 'code'
    subfolder = codefolder[codefolder.find('code') + len('code'):]
    
    # delete the first character if the first one is '/'
    if subfolder[0] == '/':
        subfolder = subfolder[1:]

        
    _, _, pipelinefolder, _ = exp_subfolders()
    
    # return the corresponding folder in pipeline folder
    codecorresfolder = os.path.join(pipelinefolder, subfolder, codefile)
    codecorresparentfolder = os.path.join(pipelinefolder, subfolder)
    
    # make folder if needed
    if makefolder == True and os.path.isdir(codecorresfolder) == False:
        os.makedirs(codecorresfolder)
        
    # make folder if needed
    if makesubfolder == True: 
        
        storefolder = os.path.join(codecorresfolder, 'store') 
        if os.path.isdir(storefolder) == False:
            os.makedirs(storefolder)

            
        outfolder = os.path.join(codecorresfolder, 'out') 
        if os.path.isdir(outfolder) == False:
            os.makedirs(outfolder)
        
    return codecorresfolder, codecorresparentfolder