import numpy as np
from cpython cimport array

def distFun(d, distFunParams):
    """
    Select the distance function parameters

    Parameters
    ----------
    disFunParams: dic
    d: series
        x["distSign"] column

    Returns
    -------
    out : list
    """


    cdef int[:] obs_max = distFunParams["obs_max"]
    cdef int[:] obs_min = distFunParams["obs_min"]
    cdef int[:] head_coef = distFunParams["head_coef"]
    cdef int[:] tail_coef = distFunParams["tail_coef"]
    cdef int[:] fit = distFunParams["cubicFit"]

    ##Put everything together to get the final function
    d = np.log(d)

    out = []

    for dist in d:
        if dist > obs_max:
            out.append(tail_coef[0] + dist*tail_coef[1])
        elif dist < obs_min:
            out.append(head_coef[0] + dist*head_coef[1])
        else:
            out.append(fit[0] + fit[1]*dist + fit[2]*(dist**2) + fit[3]*(dist**3))

    return np.exp(out)