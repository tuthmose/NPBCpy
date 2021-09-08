# G Mancini Jul 2021

from sys  import argv, exit
import argparse as argp
import mdtraj as md
import numpy as np
import scipy as sp

import npbc_cy
import npbc_analysis
import npbc_io

Parse = argp.ArgumentParser(description='Continuous F function for hydrogen bonds based on radial distribution\
                             functions between H (hydrogen) and A (acceptor) and D (donor) and A \
                             and angular distribution functions for H-D-A')
# template and gaussian settings
Parse.add_argument('-i','--input',help='input XTC trajectory',default=False,action='store')
Parse.add_argument('-t','--topology',help='topology PBD file',default=False,action='store')
Parse.add_argument('-o','--output',help='output file',default=False,action='store')
Parse.add_argument('-O','--histogram',help='histogram output file',default=False,action='store')
Parse.add_argument('-b','--begin',help='first frame to use',default=-1,type=int)
Parse.add_argument('-e','--end',help='last frame to use',default=-1,type=int)
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
Parse.add_argument("-r","--rdf",action="store",default=False,nargs=2,\
    help="rdf first peak position and half width")
Parse.add_argument("-a","--adf",action="store",default=False,nargs=2,\
    help="adf first peak position and half width")
Parse.add_argument("-x","--shift",action="store",default=False,nargs=3,\
    help="shift coordinates by x,y,z")
Parse.add_argument("-N","--norm",default=True,action="store_false",\
    help="normalization; default is to normalize")
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

if not Myarg.index:
    raise ValueError("ERROR: missing index file to use with -m")

if not Myarg.select:
    raise ValueError("ERROR: -m not selected")
else: 
    select = list(map(int, Myarg.select))
    
if not Myarg.rdf:
    raise ValueError("ERROR: -r not selected")
else: 
    rdf = np.array(list(map(float, Myarg.rdf)))/10.

if not Myarg.adf:
    raise ValueError("ERROR: -a not selected")
else: 
    adf = list(map(float, Myarg.adf))
    
try:
    nbins = int(Myarg.nbins)
except:
    raise ValueError("ERROR: wrong number of bins")

# RUN
groups = npbc_io.parse_index(Myarg.index, select)
traj, first_frame, last_frame = npbc_io.loadtrj(Myarg.begin, Myarg.end, Myarg.input, top=Myarg.topology)
frame, histH, timeF = \
    npbc_analysis.calc_hbonds(first_frame, last_frame, shift, traj, rdf, adf, bmax, hmax, dmax, \
        nbins, Myarg.norm, groups[0], groups[1], groups[2])

out = np.vstack((np.linspace(Myarg.begin, frame, frame-Myarg.begin), timeF)).transpose()
np.savetxt(Myarg.output+"_fhb.dat", out, fmt="%9.6f")

if Myarg.histogram:
    x = np.linspace(0., 3., nbins+1)
    out_hist = np.vstack(([(x[i]+x[i+1])/2. for i in range(nbins)], histH)).transpose()
    np.savetxt(Myarg.histogram+"_hist_fhb.dat", out_hist, fmt="%9.6f")
   
quit()    
