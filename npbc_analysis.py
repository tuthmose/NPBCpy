# G Mancini Sept 2021

import math
import mdtraj as md
import numpy as np
import scipy as sp

import npbc_cy

# FUNCTIONS    

## hbonds

def cont_hbonds(coords, nbins, rdf, adf, bmax, hmax, dmax, H, D, A, norm):
    """
    calculate F_HB function; if A contains more than one atom,
    return an average value and a histogram
    """
    values = npbc_cy.calc_fhb(coords, 2, rdf, adf, bmax, dmax, hmax, H, D, A)
    if len(values) > 0:
        his, rsp = np.histogram(values, bins=nbins, range=(0., 3.))
        if norm:
            return np.average(values), his/len(values)
        else:
            return np.average(values), his
    else:
        return None, None

def calc_hbonds(first_frame, last_frame, shift, traj, rdf, adf, bmax, hmax, dmax, nbins, \
    norm, H, D, A):
    """
    read frames from xtcfile, then loop over particles and distances  
    and calculate histogram for g(r); return numpy arrays 
    """
    timeF = list()
    histH = np.zeros(nbins,dtype=np.float64)
    if first_frame == -1: 
        first_frame = 0
    frame = first_frame
    #loop over all frames
    for frame in range(first_frame, last_frame):
        #calculate rdf for this frame
            X = traj.xyz[frame]
            fhb, hist = cont_hbonds(X+shift, nbins, rdf, adf, bmax, hmax, dmax, H, D, A, norm)
            if fhb is None:
                histH = histH + np.zeros(nbins)
                timeF.append(-1)
            else: 
                histH = histH + hist
                timeF.append(fhb)
    print("--- Read ",frame," frames")
    return frame, histH, timeF

## density layers

def collect_dens(atoms, ntot, Coords, natoms, radii2, MM):
    """
    loop over molecules
    """
    COM2 = list()
    Coords = MM*Coords[atoms]
    for i in range(0, ntot, natoms):
        #squared center of mass of each molecule        
        COM2.append(np.sum((np.sum(Coords[i:i+natoms,:],axis=0))**2))
    rho,edges = np.histogram(COM2,bins=radii2)
    return rho      
    
def calc_density(first_frame, last_frame, shift, vol, from_wall, traj, atoms, natoms, nmol, \
    radii, rmax, M):
    """
    read frames from xtcfile, then    
    loop over particles; returns numpy arrays 
    """   
    nbins = len(radii)-1
    atoms = tuple(atoms)
    RHO  = np.zeros(nbins, dtype=np.float64)
    RHO2 = np.zeros(nbins, dtype=np.float64)    
    ntot = len(atoms)
    Mtot = np.sum(M[:natoms])
    M = M / Mtot
    # element by element np array mult
    M = np.vstack((M, M, M)).T
    radii2 = radii**2
    print("--- Reading frames")
    frame = first_frame
    while True:
        if frame >= last_frame:
            break
        else:
            X = traj.xyz[frame]
            rho = collect_dens(atoms, ntot, X+shift, natoms, radii2 ,M)
            frame += 1
            RHO  = RHO  + rho
            RHO2 = RHO2 + rho*rho
    print("--- Done reading frames")
    print("--- Read ", frame," frames")
    frame += 1
    nm3_l = 10**(-24)
    norm = (vol*sp.constants.N_A*nm3_l)
    RHO2 = np.sqrt((RHO2/frame - (RHO/frame)**2))
    halfpoints = [radii[i-1] + (radii[i]-radii[i-1])/2.0 for i in range(1,nbins+1)]
    if from_wall is False:
        RHO  = np.vstack((halfpoints,RHO/(frame*norm),RHO2/norm))
    else:
        x_from_boundary = [ (rmax - i) for i in halfpoints]
        RHO  = np.vstack((x_from_boundary,RHO/(frame*norm),RHO2/norm))
    return RHO

## nearest neighbour

def calculate_distance(coords, nbins, nneigh, groupA, groupB):
    """
    calculate distance, sort and return nneigh distances
    """
    A = coords[groupA]
    B = coords[groupB]
    D = sp.spatial.distance.cdist(A, B, metric=Myarg.metric)
    DS = np.sort(np.ravel(D))
    return 10.*DS[nneigh]

def calc_nearest_dist(first_frame, last_frame, traj, nneigh, groupA, groupB):
    """
    read frames from xtcfile, then loop over particles and distances  
    and calculate histogram for g(r); return numpy arrays 
    """
    timeD = list()
    if first_frame == -1: 
        first_frame = 0
    frame = first_frame
    #loop over all frames
    for frame in range(first_frame, last_frame):
        #calculate for this frame
            X = traj.xyz[frame]
            dist = calculate_distance(X+shift, nbins, nneigh, groupA, groupB)
            timeD.append(dist)
    print("--- Read ",frame," frames")
    return frame, timeD

## nmol

#FUNCTIONS
def find_coms(natoms, coords, weights):
    """
    find center of mass of molecules of natoms 
    """
    nmol = coords.shape[0] // natoms
    com = np.empty((nmol, 3))
    for n in range(0, nmol, natoms):
        com[n] = np.average(coords[i:i+natoms], weights=weights[i:i+natoms], axis=0)
    return com

def calculate_coord_number(coords, weights, cutoff, natoms, nearest, groupA, groupB, groupC):
    """
    calculate number of molecules within cutoff
    """
    if natoms[0] > 1: 
        A = find_coms(natoms[0], coords[groupA], weights[groupA], com=[0])
    else:
        A = coords[groupA]
    if natoms[1] > 1: 
        B = find_coms(natoms[1], coords[groupB], weights[groupB], com[1])
    else:
        B = coords[groupB]
    nB = np.empty((A.shape[0]), dtype='int')
    DB = sp.spatial.distance.cdist(A, B, metric=Myarg.metric)    
    if groupC == False:
        for a in range(A.shape[0]):
            dd1 = DB[a] >= cutoff[0]
            dd2 = DB[a] <= cutoff[1]
            dd3 = dd1 & dd2
            nB[a] = np.count_nonzero(dd3)
    else:
        if natoms[2]>1 and not nearest: 
            C = find_coms(natoms[2], coords[groupC], weights[groupC])
        else:
            C = coords[groupC]
        DC = sp.spatial.distance.cdist(A, C, metric=Myarg.metric)
        if natoms[2]>1 and nearest:
            DC2 = np.empty((len(A), len(C)//natoms[2]))
            for a in range(A.shape[0]):
                for c in range(C.shape[0]//natoms[2]):
                    DC2[a] = np.min(DC[a,c:c+natoms[2]])
            print(DC2[a])
            DC = DC2
            print(DC.shape, C.shape, natoms[2])
        for a in range(A.shape[0]):
            dd1 = DB[a] >= cutoff[0]
            dd2 = DB[a] <= cutoff[1]
            dd3 = dd1 & dd2
            dd1 = DC[a] >= cutoff[2]
            dd2 = DC[a] <= cutoff[3]
            dd4 = dd1 & dd2
            dd3 = dd3 & dd4
            nB[a] = np.count_nonzero(dd3)        
    return nB

def calc_nmol(first_frame, last_frame, traj, natoms, cutoff, nearest, groupA, groupB, groupC):
    """
    read frames from xtcfile, then loop over particles and distances  
    and calculate histogram for g(r); return numpy arrays 
    """
    timeN = list()
    if first_frame == -1: 
        first_frame = 0
    frame = first_frame
    top = traj.topology
    weights = np.asarray([a.element.mass for a in top.atoms])
    #loop over all frames
    for frame in range(first_frame, last_frame):
        #calculate for this frame
            X = traj.xyz[frame]
            nmol = calculate_coord_number(X+shift, weights, cutoff, natoms, nearest, groupA, groupB, groupC)
            timeN.append(nmol)
    print("--- Read ",frame," frames")
    timeN = np.asarray(timeN)
    return frame, timeN

## orientation

def collect(normV, Coords, atoms, W, hfalp, versor, axis, nbins, natoms):
    """
    loop over molecules
    """
    Theta = np.zeros(nbins)
    cosines = np.zeros(nbins)
    sines   = np.zeros(nbins)
    if versor:
        for i in atoms:
            com = np.average(Coords[i:i+natoms,:], weights=W[i:i+natoms], axis=0)
            z = com[axis]
            mybin = (np.abs(hfp-z)).argmin()
            vnormal = npbc_cy.findvec(Coords[i:i+natoms,:], normV)
            cosangle = dotprod(versor, vnormal)
            cosines[mybin] += cosangle
            sines[mybin]   += sqrt(1.-cosangle**2)
    else:
        for a in range(0, len(atoms), natoms):
            i = atoms[a]
            vnormal = npbc_cy.findvec(Coords[i:i+natoms,:], normV)   
            com = np.average(Coords[i:i+natoms,:], weights=W[i:i+natoms],axis=0)
            r = np.linalg.norm(com)
            mybin = (np.abs(hfalp - r)).argmin()
            cosangle = npbc_cy.dotprod(com, vnormal)
            cosines[mybin] += cosangle
            sines[mybin]   += math.sqrt(1.-cosangle**2)
    #circ mean
    Theta = np.arctan2(sines, cosines)
    return Theta, Theta**2

def calc_orient(normV, first_frame, last_frame, traj, shift, top, group, axis, vol, \
    radii, weights, natoms):
    """
    read frames from xtcfile, then loop over particles and distances  
    and calculate histogram for g(r); return numpy arrays 
    """
    rad2deg = 180.0/np.pi
    nbins = len(radii)-1
    THETA  = np.zeros(nbins,dtype=np.float64)
    THETA2 = np.zeros(nbins,dtype=np.float64)    
    if axis == False:
        halfpoints = np.asarray([radii[i-1] + (radii[i]-radii[i-1])/2.0 for i in range(1,nbins+1)])
    else:
        raise ValueError('Cylindrical version nyi')
    #
    if first_frame == -1: 
        first_frame = 0
    versor = False
    frame = first_frame
    tf = list()
    #loop over all frames
    for frame in range(first_frame, last_frame):
        #calculate for this frame
        X = traj.xyz[frame]
        theta, theta2 = collect(normV, X+shift, group, weights, halfpoints, versor, axis, nbins, natoms)
        THETA  = THETA  + theta
        THETA2 = THETA2 + theta2
        tf.append((frame, theta))
    print("--- Read ",frame," frames")
    THETA  = theta/frame
    THETA2 = np.sqrt(THETA2/frame - THETA**2)
    THETA  = np.vstack((halfpoints,-(THETA*rad2deg-90.0),THETA2*rad2deg)).T
    tf = np.asarray(tf)
    return tf, THETA

## rdf

#FUNCTIONS    
def calculate_histogram(coords, nbins, rmax, dmax, ref, target, vol_dmax):
    """
    calculate g(r) for current frame
    """
    rho_d = .0    
    rdf = np.zeros(nbins,dtype=np.float64) 
    cn  = np.zeros(nbins,dtype=np.float64) 
    REF = coords[ref]
    TARG = coords[target]
    dREF = np.sqrt(np.sum(REF**2,axis=-1))
    mask = np.logical_and(dREF>=rmax[0], dREF<=rmax[1])
    REF = REF[mask]
    if len(np.ravel(REF))==0:
        return 0.0, 0.0, np.zeros(nbins), np.zeros(nbins)
    na = len(REF)
    nb = len(TARG)
    if na==0 or nb==0:
        return None, None, None, None
    for rp in range(na):
        D = np.sqrt(np.sum((TARG-REF[rp])**2,axis=-1))
        D = D[(D<dmax) & (D>0.0)]
        his,rsp = np.histogram(D,bins=nbins, range=(0.0,dmax))
        rho_tmp = len(D)/vol_dmax
        cn = cn + his
        rdf = rdf + his/rho_tmp
        rho_d = rho_d + rho_tmp
    cn = cn/na
    rdf = rdf/na
    rho_d = rho_d/na
    return na, rho_d, cn, rdf

def normalize(calc_cn, donorm, smooth, rmax, dmax, nbins, cn, gofr, frame):
    """
    normalize g(r) and calculate CN(r)
    """
    x = np.linspace(0.0,dmax,nbins+1)
    xp = x[1:]
    xm = x[:-1]
    xh = xp-0.5*(xp-xm)
    deltaV = (4.0*np.pi/3.0)*(xp**3-xm**3)
    if calc_cn:
       CN = np.cumsum(cn)/frame
    else:
       CN = None
    norm = frame*deltaV
    if donorm:
        gofr = gofr/norm
    if smooth:
        xs = xh[np.where(xh>smooth)]
        ns = len(xs)
        nS = len(xh)-ns
        print("--- Smoothing over",ns," points out of",nbins)
        gofr = np.piecewise(gofr,[xh<=smooth,xh>smooth],\
        [gofr[:nS],(1.+(gofr[nS:]-1.)*np.exp((xs/smooth-1)**2))])
    return gofr, xh, CN

def calc_rdf(first_frame, last_frame, nbins, calc_cn, smooth, norm, \
             RSphere, rmax, dmax, shift, traj, Ref, Target):
    """
    read frames from xtcfile, then loop over particles and distances  
    and calculate histogram for g(r); return numpy arrays 
    """
    RDF = np.zeros(nbins,dtype=np.float64)
    CN  = np.zeros(nbins,dtype=np.float64)
    frame = first_frame
    ## density
    vol = (4.0*np.pi/3.0)*(RSphere**3)
    rho = float(len(Target))/vol    
    print("--- Number density in system is ",rho)
    vol_dmax = (4.0*np.pi/3.0)*(dmax**3)
    rho  = .0
    Nref = 0
    #loop over all frames
    while True:
        if frame >= last_frame:
            break
        else:
        #calculate rdf for this frame
            X = traj.xyz[frame]
            na,rho_d,cn,rdf = calculate_histogram(X+shift, nbins, rmax, dmax, Ref, Target, vol_dmax)
            if rdf is None:
                continue
            RDF += rdf
            rho += rho_d
            Nref  += na
            CN    += cn
            frame += 1
    frame += 1
    Nref = Nref/frame
    print("--- Read ",frame," frames")
    print("--- Average number of reference molecules ",Nref)
    print("--- Average number density in dmax ",rho/(frame))
    gofr,xbins,coord_num = normalize(calc_cn, norm, smooth, rmax, dmax, nbins, CN, RDF, frame)
    return 10.0*xbins, gofr, coord_num

