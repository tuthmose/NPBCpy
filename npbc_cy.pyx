#cython: language_level=3

# G Mancini July 2021
#
# Middle level python and cython functions to analyze Link402
# MDMC simulations

from math import acos

import numpy as np
import numpy.linalg as LA
import scipy as sp

from cython import cdivision, boundscheck, wraparound
from cython.parallel import prange
    
from libcpp cimport bool    
from libc.math cimport fmin as cmin
from libc.math cimport fmax as cmax
from libc.math cimport acos as cacos
from libc.math cimport exp  as cexp
from libc.math cimport sqrt as csqrt

### python functions
def compute_neighbours(nsolute, nsolvent, natoms_solvent, nall, top, X):
    """
    return the coordinates of the first nsolvent molecules (each with 
    natoms_solvent atoms) out of nall
    """
    w = np.asarray([a.element.mass for a in top.atoms])
    com = np.average(X[:nsolute], axis=0, weights=w[:nsolute])
    D = list()
    for s in range(nall):
        nstart = nsolute + natoms_solvent*s
        com1 = np.average(X[nstart:nstart+natoms_solvent], axis=0, weights=w[nstart:nstart+natoms_solvent])
        dist = np.linalg.norm(com-com1)
        D.append(dist)
    order = np.argsort(D)
    coords = np.empty((nsolute + nsolvent*natoms_solvent, 3))
    coords[:nsolute] = X[:nsolute]
    for i in range(nsolvent):
        res = order[i]
        for j in range(natoms_solvent):
            k = i*natoms_solvent + nsolute + j
            atom = res*natoms_solvent + nsolute + j
            coords[k] = X[atom]
    return coords

def findvec(coords, normV=False):
    """
    find vector normal to the plane spanned by the atom triplet
    or vector that bisects 1 0 2 angle
    coords is a 3x3 array
    """    
    a = coords[1] - coords[0]
    b = coords[2] - coords[1]
    if normV:
        return np.cross(b,a)
    else:
        return a+b

def dotprod(v1,v2):
    """
    orientation of selected vector wrt radial dir
    """
    V1 = v1/LA.norm(v1)
    V2 = v2/LA.norm(v2)
    dot = np.dot(V1,V2)
    return dot

def calc_angle(coords, a0, a1, a2):
    """
    select atoms a0, a1, a2 and calculate
    a2_a0 /\ a2_a1 angle
    """
    v1 = coords[a2] - coords[a0]
    v2 = coords[a2] - coords[a1]
    cos = dotprod(v1,v2)
    return acos(cos)
    
def triple_product(v0, v1, v2):
    """
    Return the triple triple_product i.e.
    the volume spanned by v0, v1, v2
    """
    v3 = np.cross(v0, v1)
    V  = np.dot(v2, v3)
    return V

### cython functions
@wraparound(False)  
@boundscheck(False)
def calc_fhb(coords, what, rdf, adf, bmax, dmax, hmax, H, D, A, radius):
    """
    calculate ADF (what=1) or FHB (what=2) function on a given frame
    """
    # definitions
    cdef int i, j, k, hi, aj, dk, l
    cdef int NH, ND, NA
    cdef double b0, b1, d0, d1, h0, h1
    cdef double peak_rdf, width_rdf, peak_adf, width_adf
    cdef double const, n1, n2, n3, theta, out, Ar, Bt, R, Rmin, Rmax
    cdef double [:,:] cX = coords
    cdef int [:] cH = np.asarray(H,dtype=np.intc)
    cdef int [:] cD = np.asarray(D,dtype=np.intc)
    cdef int [:] cA = np.asarray(A,dtype=np.intc)    
    cdef double v1[3]
    cdef double v2[3]
    cdef double v3[3]
    # check what to do
    Rmin = radius[0]
    Rmax = radius[1]
    if what != 1 and what != 2:
        print("what you provided to calc_fhb " + str(what))
        raise ValueError("Unknown value")
    const = 180.0/np.pi
    b0 = bmax[0]
    b1 = bmax[1]
    d0 = dmax[0]
    d1 = dmax[1]
    h0 = hmax[0]
    h1 = hmax[1]    
    if what == 2:
        peak_rdf  = rdf[0]
        width_rdf = rdf[1]
        peak_adf  = adf[0] / const
        width_adf = adf[1] / const
    elif what == 1:
        peak_rdf  = 0.
        width_rdf = 0.
        peak_adf  = 0.
        width_adf = 0.
    NH = len(H)
    ND = len(D)
    NA = len(A)
    values = list()
    # start a lot of nested loops (ugly!)
    for j in range(NA):
        out = 0.
        R = csqrt(cX[cA[j],0]*cX[cA[j],0]+cX[cA[j],1]*cX[cA[j],1]+cX[cA[j],2]*cX[cA[j],2])
        if R < Rmin or R > Rmax:
            continue
        for i in range(NH):
            if cA[j] != cH[i]:
                n1 = 0.
                for l in range(3):
                    v1[l] = cX[cA[j],l] - cX[cH[i],l]
                    n1 = n1 + v1[l]*v1[l]
                n1 = csqrt(n1)
                if n1 >= h0 and n1 <= h1:
                    for k in range(ND):
                        if cA[j] != cD[k] and cD[k] != cH[i]:
                            n2 = 0.
                            n3 = 0.
                            for l in range(3):
                                v2[l] = cX[cA[j],l] - cX[cD[k],l]
                                n2 = n2 + v2[l]*v2[l]
                                v3[l] = cX[cD[k],l] - cX[cH[i],l]
                                n3 = n3 + v3[l]*v3[l]
                            n2 = csqrt(n2)
                            n3 = csqrt(n3)
                            if n2>=d0 and n2<=d1 and n3>=b0 and n3<=b1:
                                theta = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(n1*n2)
                                if what==1:
                                    out = const*cacos(theta)
                                    values.append(out)
                                elif what==2:
                                    theta = cacos(theta)
                                    if (peak_rdf - n1) < 0.:
                                        Ar = cexp( -(peak_rdf - n1)*(peak_rdf - n1)/( 2.*width_rdf*width_rdf ))
                                    else:
                                        Ar = 1.
                                    if (peak_adf - theta) < 0.:
                                        Bt = cexp( -(peak_adf - theta)*(peak_adf - theta)/( 2.*width_adf*width_adf ))
                                    else:
                                        Bt = 1.
                                    out = out + Ar * Bt
        if what == 2:
            values.append(out)                                    
    return values                                
