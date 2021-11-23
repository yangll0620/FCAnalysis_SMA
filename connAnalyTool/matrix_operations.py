import numpy as np
from numpy.linalg import svd
from scipy import stats

def similarity_cosin(u1, u2):
    """
        calulate the cosine similarity of two vectors u1 and u2
    """

    similarity_cos = np.dot(u1,u2)/(np.linalg.norm(u1)* np.linalg.norm(u2))
    
    return similarity_cos




def cosSimilarity_SVDComps(matrix_Con):
    """
        detect the change points using connectivity matrix based on cosine similarity of SVD components
		each matrix (n_chns * n_chns) is first decomposed using SVD
		
        Then, cosSimilarity is calcuated using weighted cosin similarity
		
        @ parameter 
            matrix_Con: connectivity matrix, n_chns * n_chns * n_times
		
        return 
            diff: cos similarity differences
    """

    n_chns, n_chns, n_times = matrix_Con.shape

    # SVD decomposition for the connectivity matrix in each time i_time
    Us = np.zeros((n_chns, n_chns, n_times))
    Ss = np.zeros((n_chns, n_times))
    for i_time in range(n_times):
        matrix_conEach = matrix_Con[:,:,i_time]
        u,s,vh = svd(matrix_conEach, full_matrices = True)
        Us[:,:,i_time] = u
        Ss[:,i_time] = s # s is a vector
        del matrix_conEach, u, s, vh


    # calcuate different matrix based on cosine similarity 
    diff = np.zeros((n_times-1,))
    for i_time in range(1, n_times):
        U_current, S_current = Us[:,:,i_time], Ss[:,i_time]
        U_previous, S_previous = Us[:,:,i_time - 1],Ss[:,i_time -1] 
        n_comp = U_current.shape[1]
        similarity_cos = np.zeros((n_comp,))
        for i_comp in range(U_current.shape[1]):
            u1, u2 = U_current[:,i_comp], U_previous[:,i_comp]
            similarity_cos[i_comp] = similarity_cosin(u1, u2)
            del u1, u2
        S_mean = (S_current + S_previous)/2
        weight = S_mean /np.sum(S_mean)
        diff[i_time-1] = np.sqrt(np.sum(weight * (similarity_cos*similarity_cos)))
        del weight, S_mean
        del U_current, U_previous, S_previous, S_current
        del similarity_cos
        del n_comp, i_comp


    return diff