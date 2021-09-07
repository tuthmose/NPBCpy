# G Mancini Jul 2021

from collections import defaultdict
from sys import argv, exit
import numpy as np
import scipy as sp

import myparse

Parse = argp.ArgumentParser(description='Calculate density mean and variance  for a molecular group \
                    in concentric layers of constant volume or radius')
Parse.add_argument('-i','--input',help='input XTC trajectory',default=False,action='store')
Parse.add_argument('-o','--output',help='output filename for ASCII g(r)',default=False,action='store')
Parse.add_argument('-t','--topology',help='topology PBD file',default=False,action='store')
Parse.add_argument("-I","--index",dest="index",default=False,action="store",
    help="GROMACS index file to use with --molecules")
Parse.add_argument('-b','--begin',help='first frame to use',default=-1,type=int)
Parse.add_argument('-e','--end',help='last frame to use',default=-1,type=int)
Parse.add_argument("-m","--molecules",default=False,action="store",\
Parse.add_argument("-N","--natoms",default=False,action="store",help="number of atoms
    in select molecules")
help="select target molecules; use group numbers")
Parse.add_argument("-R","--rsphere",default=False,action="store",\
    help="radius for spherical boxes (angs)")
Parse.add_argument("-n","--nbins",action="store",type=int,default=500,\
    help="number of bins to use; default=500")
Parse.add_argument("-V","--volume",default=False,action="store_true",\
    help="Use constant volume layers instead of constant radius")
Parse.add_argument("-w","--from_wall",default=False,action="store_true",\
    help="plot distance from sphere boundary")
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

if not Mayrg.natoms:
    print("Perception of residues from topology nyi")
    raise ValueError("Missing number of atoms in target molecule type")
else:
    natoms = int(Myarg.natoms)
    
# CONSTANTS
DIM = 3
NA  = sp.constants.N_A
nm3_l = 10**(-24)
csph = (4.0/3.0)*np.pi
    
# RUN
target = myparse.parse_index(Myarg.index, select)

traj, first_frame, last_frame = myparse.loadtrj(Myarg.begin, Myarg.end, Myarg.input, top=Myarg.topology)
top = traj.topology
W = np.array([a.element.mass for a in list(top.atoms)[target]],dtype=np.float32)
nmol = len(target)/natoms

vol,Radii = myparse.gen_bin_radii(target, Nat, nbins, RSphere)
rho = libl402.calc_density(first_frame, last_frame, vol, traj, target, natoms, nmol, Radii, W)
np.savetxt(options.outname+".dat",RHO.T,fmt="%15.6f")
        
quit()        
