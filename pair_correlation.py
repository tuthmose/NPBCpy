#!/usr/bin/env python3
# G Mancini Jul 2021

from sys  import argv, exit
import argparse as argp
import mdtraj as md
import numpy as np
import scipy as sp

import libl402
import myparse

Parse = argp.ArgumentParser(description='Pair correlation function (g(r)) for spherical boxes based on atom groups')
# template and gaussian settings
Parse.add_argument('-i','--input',help='input XTC trajectory',default=False,action='store')
Parse.add_argument('-o','--output',help='output filename for ASCII g(r)',default=False,action='store')
Parse.add_argument('-t','--topology',help='topology PBD file',default=False,action='store')
Parse.add_argument('-I','--index',help='index file name (GMX like)',default=False,action='store')
Parse.add_argument('-b','--begin',help='first frame to use',default=-1,type=int)
Parse.add_argument('-e','--end',help='last frame to use',default=-1,type=int)
Parse.add_argument("-m","--molecules",dest="select",default=False,action="store",\
    nargs=2,help="select reference \and target atoms; com nyi; use group numbers")
Parse.add_argument("-n","--nbins",action="store",default=500,\
    help="number of bins to use; default=500")
Parse.add_argument("-c","--cn",action="store_true",default=False,\
    help="calculate g(r) integral (coordination number)")
Parse.add_argument("-r","--rmax",action="store",default=False,nargs=2,\
    help="minimum and maximum distance (angs) from the center of the sphere for reference atoms")
Parse.add_argument("-d","--dmax",action="store",default=False,
    help="maximum distance from the center of the sphere for target atoms  (angs)")
Parse.add_argument("-R","--rsphere",default=False,action="store",\
    help="radius for spherical boxes (angs)")
Parse.add_argument("-N","--norm",default=True,action="store_false",\
    help="normalization; default is to normalize")
Parse.add_argument("-s","--smooth",action="store",default=0.0,\
    help="for r>smooth, transform g(r) as 1+[g(r)-1]exp(-r/smooth-1)**2; default is 0, i.e. do not smooth")
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

if not Myarg.output:
    raise ValueError("Missing output file name")

try:
    RSphere = float(Myarg.rsphere)/10.
except:
    raise ValueError("ERROR: sphere radius not set")
    
if Myarg.rmax: 
    rmax = np.array(list(map(float,Myarg.rmax)))/10.
else:
    raise ValueError("ERROR: rmax not set")

if Myarg.dmax: 
    dmax = float(Myarg.dmax)/10.
else:
    raise ValueError("ERROR: dmax not set")

smooth =float(Myarg.smooth)

if Myarg.shift:
    shift = np.array(map(float,Myarg.shift))/10.
else:
    shift = np.zeros(3)
    
try:
    nbins = int(Myarg.nbins)
    thick = dmax/(float(nbins))
    print("--- Using ",nbins," bins ",thick," nm wide")
except:
    raise ValueError("ERROR: wrong number of bins")

if not Myarg.index:
    raise ValueError("ERROR: missing index file to use with -m")
if not Myarg.select:
    raise ValueError("ERROR: -m not selected")
else: 
    select = list(map(int, Myarg.select))
    
# RUN
groups = myparse.parse_index(Myarg.index, select)
traj, first_frame, last_frame = myparse.loadtrj(Myarg.begin, Myarg.end, Myarg.input, top=Myarg.topology)
X,Rdf,CN = libl402.calc_rdf(first_frame, last_frame, nbins, Myarg.cn, Myarg.smooth, \
    Myarg.norm, RSphere, rmax, dmax, shift, traj, groups[0], groups[1])

if Myarg.cn:
    Out = np.vstack((X, Rdf, CN))
else:
    Out = np.vstack((X, Rdf))
Out = Out.transpose() 
np.savetxt(Myarg.output+".dat",Out,fmt="%9.6f")
    
quit()    
