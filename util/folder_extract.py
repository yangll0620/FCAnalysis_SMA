# # folder extraction return absoluate folder path
# 
# * exp_subfolders(): return datafolder, codefolder, pipelinefolder, outputfolder
# 
# * code_corresfolder(): return the corresponding folder for codefolder/codefilepath,  create the corresponding folder and/or the'out' and 'store' subfolders


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





def code_corresfolder(codefilepath, makefolder = True):
    """
        code_corresfolder 
            return the corresponding folder for codefilepath
            create the corresponding folder if not exist (makefolder = True)
            
        args:
            
            codefilepath: code file path with suffix (i.e path/to/code/examplecode.py)  
            makefolder: tag for creating the corresponding codefolder (default True)
            
                 
        return:
            code_corres_folder: the corresponding folder for codefilepath
            codecorresparentfolder: the corresponding parent folder
    """
    


    path_no_suffix = codefilepath.replace('.py', '')
    code_corres_folder = path_no_suffix.replace('code', 'pipeline')
    code_corresparent_folder = os.path.split(code_corres_folder)[0]
    
    # make folder if needed
    if makefolder == True and os.path.isdir(code_corres_folder) == False:
        os.makedirs(code_corres_folder)
        
        
    return code_corres_folder, code_corresparent_folder