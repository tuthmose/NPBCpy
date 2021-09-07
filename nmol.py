# G Mancini Aug 2021

from sys  import argv, exit
import argparse as argp
import mdtraj as md
import numpy as np
import scipy as sp

import calclow
import myparse

Parse = argp.ArgumentParser(description='print the number of molecules from groups B and \
            (optionally) C between rmin and rmax of the COM of group A')
# template and gaussian settings
Parse.add_argument('-i','--input',help='input XTC trajectory',default=False,action='store')
Parse.add_argument('-t','--topology',help='topology PBD file',default=False,action='store')
Parse.add_argument('-o','--output',help='output file',default=False,action='store')
Parse.add_argument('-O','--histogram',help='histogram output file',default=False,action='store')
Parse.add_argument('-b','--begin',help='first frame to use',default=-1,type=int)
Parse.add_argument('-e','--end',help='last frame to use',default=-1,type=int)
Parse.add_argument('-I','--index',help='index file name (GMX like)',default=False,action='store')
Parse.add_argument("-m","--molecules",dest="select",default=False,action="store",\
    nargs=3,help="number of groups 2 and 3 within interval of group 1; use group numbers")
Parse.add_argument("-c","--cutoff",action="store",default=False,nargs=4,help="minimum and maximum distance")
Parse.add_argument("--nearest",action="store_true",default=False,help="for group C use nearest atom")
Parse.add_argument("-M","--metric",action="store",default="euclidean",help="distance type")
Parse.add_argument("-n","--nbins",action="store",default=200,help="number of bins to use")
Parse.add_argument("-N","--natoms",action="store",default=False,nargs=3,\
                   help="number of atoms in each molecule in group 1 and 2")
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

if not Myarg.select:
    raise ValueError("ERROR: --cutoff not selected")
else: 
    cutoff = np.asarray(list(map(float, Myarg.cutoff)))/10.
    assert cutoff[1] > cutoff[0]
    if select[2] != -1:
        assert cutoff[3] > cutoff[2]

try:
    nbins = int(Myarg.nbins)
except:
    raise ValueError("ERROR: wrong number of bins")

natoms = tuple(map(int,Myarg.natoms))
assert natoms[0] >= 0
assert natoms[1] >= 1
if select[2] != -1:
    assert natoms[2] != -1

# RUN
if select[2] != -1:
    groups = myparse.parse_index(Myarg.index, select)
else:
    groups = myparse.parse_index(Myarg.index, select[:2])
    groups.append(False)
traj, first_frame, last_frame = myparse.loadtrj(Myarg.begin, Myarg.end, Myarg.input, top=Myarg.topology)
if natoms[0] > 1:
    if traj.xyz[0][groups[0]].shape[0] % natoms[0] != 0:
        raise ValueError("number of molecules in group A is not integer")
if natoms[1] > 1:
    if traj.xyz[0][groups[1]].shape[0] % natoms[1] != 0:
        raise ValueError("number of molecules in group B is not integer")
if natoms[2] > 1:
    if traj.xyz[0][groups[2]].shape[0] % natoms[1] != 0:
        raise ValueError("number of molecules in group C is not integer")        
    
frame, timeN = libl402.calc_nmol(first_frame, last_frame, traj, natoms, cutoff, Myarg.nearest, groups[0], groups[1], groups[2])

out = np.column_stack((np.linspace(Myarg.begin, frame, frame-Myarg.begin), timeN))
np.savetxt(Myarg.output+"_dist_nmol.dat", out, fmt="%9.6f")

if Myarg.histogram:
    hist, edges = np.histogram(timeD, bins=Myarg.nbins)
    hist = np.append(hist, (-1))
    out = np.vstack((edges, hist)).transpose()
    np.savetxt(Myarg.histogram+"_hist_nmol.dat", out, fmt="%9.6f")
   
quit()    
