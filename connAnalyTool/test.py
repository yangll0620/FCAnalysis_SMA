import os, sys

# extract the exp folder path
currfolder = os.getcwd()
codefolder = currfolder[0 : currfolder.find('code')+len('code')]

# add path the exp folder
sys.path.insert(1, codefolder)
print(codefolder)



# # import util/folder_extract
# from util.folder_extract import exp_subfolders, code_corresfolder

# py_name = os.path.basename(__file__)[0:-3]

# # corresfolder
# corresfolder, correparentfolder = code_corresfolder(os.getcwd(), py_name)

# print(corresfolder)