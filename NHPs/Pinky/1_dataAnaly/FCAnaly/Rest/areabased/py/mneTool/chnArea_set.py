import pandas as pd
import numpy as np

def assign_coord2chnArea(area_coord_file, chnAreas):
    """
        assign the xy coord of each area


        Args:
            area_coord_file: file containing the x y coord for each area (normally predefined)


            chnAreas: a list of areas representing the corresponding area for each channel


        Return:
            df_chninf: DataFrame with chnAreas, coord_x and coord_y

    """
    # load channel coord from area_coord_file
    df = pd.read_csv(area_coord_file, header = 0)

    # fill in the x,y coordinates of each area in chnAreas based on the values in df_chninf
    coord_x, coord_y = np.zeros(shape = [len(chnAreas), ]), np.zeros(shape = [len(chnAreas), ])
    for i, chnArea in enumerate(chnAreas):
        
        mask_area = (df['brainarea'] == chnArea)
        
        if len(df['brainarea'][mask_area].index) == 0:
            continue

        x, y = df['simulated_x'][mask_area].to_numpy(), df['simulated_y'][mask_area].to_numpy()

        coord_x[i], coord_y[i] = x, y
        
        
        del mask_area, x, y

    df_chninf = pd.DataFrame(data = {'chnAreas': chnAreas, 'coord_x': coord_x, 'coord_y': coord_y})
        
    return  df_chninf
