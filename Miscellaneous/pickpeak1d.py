from __future__ import print_function

def peaks1d(fid, threshold=0.1):
    '''
    Extract peaks from 1d from FID
    '''
    listpk = np.where(((fid > threshold*np.ones(fid.shape))&# thresholding
                    (fid > np.roll(fid,  1, 0)) &     
                    (fid > np.roll(fid, -1, 0)) )) # roll 1 and -1 on axis 0
    return listpk[0]