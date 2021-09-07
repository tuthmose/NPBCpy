# G Mancini Aug 2021

from sys  import argv, exit
import argparse as argp
import mdtraj as md
import numpy as np
import scipy as sp

import calclow
import ibl402
import myparse

Parse = argp.ArgumentParser(description='n-neighbour distance calculation')
# template and gaussian settings
Parse.add_argument('-i','--input',help='input XTC trajectory',default=False,action='store')
Parse.add_argument('-t','--topology',help='topology PBD file',default=False,action='store')
Parse.add_argument('-o','--output',help='output file',default=False,action='store')
Parse.add_argument('-O','--histogram',help='histogram output file',default=False,action='store')
Parse.add_argument('-b','--begin',help='first frame to use',default=-1,type=int)
Parse.add_argument('-e','--end',help='last frame to use',default=-1,type=int)
Parse.add_argument('-I','--index',help='index file name (GMX like)',default=False,action='store')
Parse.add_argument("-m","--molecules",dest="select",default=False,action="store",\
    nargs=2,help="distance between this atoms; com nyi; use group numbers")
Parse.add_argument("-M","--metric",action="store",default="euclidean",help="distance type")
Parse.add_argument("-n","--nbins",action="store",default=200,help="number of bins to use")
Parse.add_argument("-N","--neigh",action="store",default=1,help="print nth neighbour distances")
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

if Myarg.shift:
    shift = np.array(map(float,Myarg.shift))/10.
else:
    shift = np.zeros(3)

if not Myarg.index:
    raise ValueError("ERROR: missing index file to use with -m")

if not Myarg.select:
    raise ValueError("ERROR: --molecules not selected")
else: 
    select = list(map(int, Myarg.select))
   
try:
    nbins = int(Myarg.nbins)
except:
    raise ValueError("ERROR: wrong number of bins")

nneigh = int(Myarg.neigh)

# RUN
groups = myparse.parse_index(Myarg.index, select)
traj, first_frame, last_frame = myparse.loadtrj(Myarg.begin, Myarg.end, Myarg.input, top=Myarg.topology)
frame, timeD = libl402.calc_nearest_dist(first_frame, last_frame, traj, nneigh, groups[0], groups[1])

out = np.vstack((np.linspace(Myarg.begin, frame, frame-Myarg.begin), timeD)).transpose()
np.savetxt(Myarg.output+"_dist_n"+str(nneigh)+".dat", out, fmt="%9.6f")

if Myarg.histogram:
    hist, edges = np.histogram(timeD, bins=Myarg.nbins)
    hist = np.append(hist, (-1))
    out = np.vstack((edges, hist)).transpose()
    np.savetxt(Myarg.histogram+"_hist_n"+str(nneigh)+".dat", out, fmt="%9.6f")
   
quit()    
