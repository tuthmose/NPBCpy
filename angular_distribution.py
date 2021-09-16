#!/usr/bin/env python3
# G Mancini Jul 2021

from scipy.spatial.distance import cdist,pdist,squareform
from sys  import argv, exit
import argparse as argp
import mdtraj as md
import numpy as np
import scipy as sp

import myparse
import calclow

Parse = argp.ArgumentParser(description='Angular distribution function for spherical boxes between  Hydrogen (H)\,
                            donor (D) and acceptor (A) atoms')
# template and gaussian settings
Parse.add_argument('-i','--input',help='input XTC trajectory',default=False,action='store')
Parse.add_argument('-a','--adf',help='ASCII angular d. f. in [0,pi]',default=False,action='store')
Parse.add_argument('-A','--aver',help='ASCII Average angle for each frame',default=False,action='store')
Parse.add_argument('-t','--topology',help='topology PBD file',default=False,action='store')
Parse.add_argument('-b','--begin',help='first frame to use',default=-1,type=int)
Parse.add_argument('-e','--end',help='last frame to use',default=-1,type=int)
Parse.add_argument('-C','--do_cython',help='use Cython version with len(A)=1',default=False,\
    action='store_true')
Parse.add_argument('-I','--index',help='index file name (GMX like)',default=False,action='store')
Parse.add_argument("-m","--molecules",dest="select",default=False,action="store",\
    nargs=3,help="select atom 1 (hydrogen), 2 (donor) and 3 (acceptor); com nyi; use group numbers")
Parse.add_argument("-n","--nbins",action="store",default=500,\
    help="number of bins to use")
Parse.add_argument("-B","--bmax",action="store",default=False,nargs=2,\
    help="minimum and maximum distance between atom1 and atom2")
Parse.add_argument("-H","--hmax",action="store",default=False,nargs=2,\
    help="minimum and maximum distance between atom1 and atom3")
Parse.add_argument("-D","--dmax",action="store",default=False,nargs=2,\
    help="minimum and maximum distance between atom2 and atom3")
Parse.add_argument("-R","--rsphere",default=False,action="store",\
    help="radius for spherical boxes (angs)")
Parse.add_argument("-N","--norm",default=True,action="store_false",\
    help="normalization; default is to normalize")
Parse.add_argument("-x","--shift",action="store",default=False,nargs=3,\
    help="shift coordinates by x,y,z")
Myarg = Parse.parse_args()
print(Myarg)

# SETUP
print ("WARNING: internally working in nm, I/O is angstroms")

if not Myarg.input:
    raise ValueError("Missing input file name")

if not Myarg.topology:
    raise ValueError("Missing topology file name")

if not Myarg.adf:
    raise ValueError("Missing ADF file name")

try:
    RSphere = float(Myarg.rsphere)/10.
except:
    raise ValueError("ERROR: sphere radius not set")

if Myarg.bmax: 
    bmax = np.array(list(map(float,Myarg.bmax)))/10.
else:
    raise ValueError("ERROR: bmax not set")
    
if Myarg.hmax: 
    hmax = np.array(list(map(float,Myarg.hmax)))/10.
else:
    raise ValueError("ERROR: hmax not set")

if Myarg.dmax: 
    dmax = np.array(list(map(float,Myarg.dmax)))/10.
else:
    raise ValueError("ERROR: dmax not set")

if Myarg.shift:
    shift = np.array(map(float,Myarg.shift))/10.
else:
    shift = np.zeros(3)
    
try:
    nbins = int(Myarg.nbins)
except:
    raise ValueError("ERROR: wrong number of bins")

if not Myarg.index:
    raise ValueError("ERROR: missing index file to use with -m")
if not Myarg.select:
    raise ValueError("ERROR: -m not selected")
else: 
    select = list(map(int, Myarg.select))
    
#FUNCTIONS    

def calculate_histogram(coords, nbins, bmax, hmax, dmax, H, D, A):
    """
    calculate average angle and assign histogram bin)
    """
    # sort groups to loop from small to big
    if not Myarg.do_cython:
        distB = cdist(coords[H],coords[D])
        distH = cdist(coords[H],coords[A])
        distD = cdist(coords[D],coords[A])
        m1 = distB >= bmax[0] 
        m2 = distB <= bmax[1] 
        maskB = m1 & m2
        m1 = distH >= hmax[0] 
        m2 = distH <= hmax[1] 
        maskH = m1 & m2
        m1 = distD >= dmax[0] 
        m2 = distD <= dmax[1] 
        maskD = m1 & m2
        values = list()
        for i, hi in enumerate(H):
            for j, aj in enumerate(A):
                if hi != aj and maskH[i,j]:
                    for k, dk in enumerate(D):
                        if hi != dk and aj != dk and maskD[k,j] and maskB[i,k]:
                            angle = (180./np.pi)*calclow.calc_angle(coords, hi, dk, aj)
                            values.append(angle)
    else:
        #calc_fhb(coords, what, rdf, adf, bmax, dmax, hmax, H, D, A)
        values = calclow.calc_fhb(coords, 1, None, None, bmax, dmax, hmax, H, D, A)        
    if len(values) > 0:
        his,rsp = np.histogram(values, bins=nbins, range=(0.0,180.))
        if Myarg.norm:
            return np.average(values), his/len(values)
        else:
            return np.average(values), his
    else:
        return None, None

def read_trajectory(first_frame, last_frame, nbins, bmax, hmax, dmax, traj, H, D, A):
    """
    read frames from xtcfile, then loop over particles and distances  
    and calculate histogram for g(r); return numpy arrays 
    """
    ADF = np.zeros(nbins,dtype=np.float64)
    timeA = list()
    if first_frame == -1: 
        first_frame = 0
    frame = first_frame
    ## density
    vol = (4.0*np.pi/3.0)*(RSphere**3)
    print("--- Number density for atom 1 is ",float(len(H))/vol)
    print("--- Number density for atom 2 is ",float(len(D))/vol)
    print("--- Number density for atom 3 is ",float(len(A))/vol)
    #loop over all frames
    for frame in range(first_frame, last_frame):
        #calculate rdf for this frame
            X = traj.xyz[frame]
            angle, adf = calculate_histogram(X+shift, nbins, bmax, hmax, dmax, H, D, A)
            if angle is None:
               timeA.append(-1)
               ADF = ADF + np.zeros(nbins)
            else: 
                ADF = ADF + adf
                timeA.append(angle)
    print("--- Read ",frame," frames")
    if Myarg.norm:
        ADF = ADF/frame
    return frame, ADF, timeA

# RUN
groups = myparse.parse_index(Myarg.index, select)
traj, first_frame, last_frame = myparse.loadtrj(Myarg.begin, Myarg.end, Myarg.input, top=Myarg.topology)
frame, ADF, timeA = read_trajectory(first_frame, last_frame, nbins, bmax, hmax, dmax, \
    traj, groups[0], groups[1], groups[2])

x = np.linspace(0., 180., nbins+1)
out_adf = np.vstack(([(x[i]+x[i+1])/2. for i in range(nbins)], ADF)).transpose()
np.savetxt(Myarg.adf+".dat", out_adf, fmt="%9.6f")

if Myarg.aver:
    out_time = np.vstack((np.linspace(Myarg.begin, frame, frame-Myarg.begin), timeA)).transpose()
    np.savetxt(Myarg.aver+".dat", out_time, fmt="%9.6f")
    
quit()    
