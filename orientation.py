#!/usr/bin/env python3
# G Mancini Aug 2021

from math import sqrt
from sys  import argv, exit

import argparse as argp
import mdtraj as md
import numpy as np
import scipy as sp

import npbc_io
import npbc_analysis
import npbc_cy

Parse = argp.ArgumentParser(description='orientation of a vector defined by a atom triplet wrt to the sphere axis')
# template and gaussian settings
Parse.add_argument('-i','--input',help='input XTC trajectory',default=False,action='store')
Parse.add_argument('-t','--topology',help='topology PBD file',default=False,action='store')
Parse.add_argument('-o','--output',help='output file',default=False,action='store')
Parse.add_argument('-b','--begin',help='first frame to use',default=-1,type=int)
Parse.add_argument('-e','--end',help='last frame to use',default=-1,type=int)
Parse.add_argument('-s','--skip',help='skip every nth frame',default=-1,type=int)
Parse.add_argument('-I','--index',help='index file name (GMX like)',default=False,action='store')
Parse.add_argument('-l','--layers',help='number of layers to use',default=10,type=int)
Parse.add_argument("-m","--molecules",dest="select",default=False,action="store",\
    nargs=3,help="group of atoms defining triplets; com nyi; use group numbers")
Parse.add_argument("-x","--shift",action="store",default=False,nargs=3,\
    help="shift coordinates by x,y,z")
Parse.add_argument("-R","--rad-sphere",dest="sphere",default=False,action="store",\
    help="radius for spherical boxes (nm)")
Parse.add_argument("-n","--nbins",action="store",type=int,default=10,\
    help="number of bins to use; default=10")
Parse.add_argument("-r","--rmax",action="store",default=False,nargs=2,\
    help="minimum and maximum distance (angs) from the center of the sphere for reference atoms")
Parse.add_argument("-a","--axis",dest="axis",default=False,action="store",\
    help="calculate distribution along an axis wrt to a triclinic box")
Parse.add_argument("-c","--cutoff",dest="cutoff",default=False,action="store",nargs=2)
Parse.add_argument("-V","--volume",dest="volume",default=False,action="store_true",\
    help="Use constant volume layers instead of constant radius")
Parse.add_argument("-v","--normV",dest="normV",default=False,action="store_true",\
    help="Use vector orthogonal to plane spanned by triplet instead of bisecting one")
Parse.add_argument("-N","--normV",dest="normV",default=False,action="store_true",\
    help="Use vector orthogonal to plane spanned by triplet instead of bisecting one")
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
    raise ValueError("ERROR: -m not selected")
else: 
    select = list(map(int, Myarg.select))

try:
    nlayers = int(Myarg.layers)
except:
    raise ValueError("ERROR: wrong number of layers")
        
if Myarg.cutoff:
    rmin = 0.1*float(Myarg.cutoff[0])
    rmax = 0.1*float(Myarg.cutoff[1])
    assert rmax>0 and rmax <= RSphere
    assert rmin>0 and rmin<rmax

if Myarg.axis:
    axis = int(Myarg.axis)
    raise ValueError("Polyhedral bins nyi")
else:
    axis = False

if Myarg.rmax: 
    rmax = np.array(list(map(float,Myarg.rmax)))/10.
else:
    raise ValueError("ERROR: rmax not set")

if not Myarg.natoms:
    print("Perception of residues from topology nyi")
    raise ValueError("Missing number of atoms in target molecule type")
else:
    natoms = int(Myarg.natoms)

# RUN
group = myparse.parse_index(Myarg.index, select)
if len(group) % 3 != 0:
    raise ValueError("ERROR: selected group should contain a multiple of 3 atoms")

traj, first_frame, last_frame = myparse.loadtrj(Myarg.begin, Myarg.end, Myarg.input, top=Myarg.topology)
top = traj.topology()
W = np.array([a.element.mass for a in atoms[target]],dtype=np.float32)

theta = npbc_analysis.calc_orient(Myarg.normV, first_frame, last_frame, traj, \
    shift, top, group, axis, vol, radii, W, natoms)

np.savetxt(Myarg.outname+".dat",theta.T,fmt="%16.9f")
 
quit()    
