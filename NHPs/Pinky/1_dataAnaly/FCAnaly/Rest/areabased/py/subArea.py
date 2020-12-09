import numpy as np


def ciCOH_select(ciCOH, chnAreas, areas_used):
    """
        select ciCOH only contains the areas in areas_used
        
        Arg:
            ciCOH: the ciCOH nchns * nchns
            chnAreas: area list
            areas_used: list containing the areas used
    """

    chnAreas_new = chnAreas.copy()
    ciCOH_new = ciCOH.copy()

    ### extract idx_del ###
    if 'GP' in areas_used:
        areas_used.remove('GP')
        areas_used = areas_used + [area for area in chnAreas_new if 'gp' in area]

    if 'STN' in areas_used:
        areas_used.remove('STN')
        areas_used = areas_used + [area for area in chnAreas_new if 'stn' in area]

    idx_del = []
    for i, area in enumerate(chnAreas_new):
        if area not in areas_used:
            idx_del.append(i)
    idx_del.reverse()

    ###  generate the used df_chninf ###
    for i in idx_del:
        del chnAreas_new[i]


    ### generate used ciCOH ###
    ciCOH_new = np.delete(ciCOH_new, idx_del, axis = 0)
    ciCOH_new = np.delete(ciCOH_new, idx_del, axis = 1)


    return ciCOH_new, chnAreas_new
