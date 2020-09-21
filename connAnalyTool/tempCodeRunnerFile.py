# import util/folder_extract
from util.folder_extract import exp_subfolders, code_corresfolder

py_name = os.path.basename(__file__)[0:-3]

# corresfolder
corresfolder, correparentfolder = code_corresfolder(os.getcwd(), py_name)

print(corresfolder)