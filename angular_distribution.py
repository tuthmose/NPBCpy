#!/usr/bin/env python3
# G Mancini Jul 2021

from scipy.spatial.distance import cdist,pdist,squareform
from sys  import argv, exit
import argparse as argp
import mdtraj as md
import numpy as np
import scipy as sp

import npbc_analysis
import npbc_io

Parse = argp.ArgumentParser(description='Angular distribution function for spherical boxes between  Hydrogen (H), donor (D) and acceptor (A) atoms')
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
Parse.add_argument("-R","--radius",default=False,action="store",nargs=2,\
    help="include molecules between Rmin and Rmax")
Parse.add_argument("--rsphere",default=False,action='store',help='radius of spherical box')
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

if not Myarg.rsphere:
    raise ValueError("ERROR: missing spherical box radius")
else:
    rsphere = float(Myarg.rsphere)/10.

if not Myarg.radius:
    radius = (0., rsphere)
else:
    radius = list(map(float, Myarg.radius))

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
    
# RUN
groups = npbc_io.parse_index(Myarg.index, select)
traj, first_frame, last_frame = npbc_io.loadtrj(Myarg.begin, Myarg.end, Myarg.input, top=Myarg.topology)

adf, timeA = npbc_analysis.calc_adf(first_frame, last_frame, shift, Myarg.do_cython, nbins, \
    bmax, hmax, dmax, traj, Myarg.norm, radius, groups[0], groups[1], groups[2])

np.savetxt(Myarg.adf+".dat", adf, fmt="%9.6f")

if Myarg.aver:
    np.savetxt(Myarg.aver+".dat", timeA, fmt="%9.6f")
    
quit()    
