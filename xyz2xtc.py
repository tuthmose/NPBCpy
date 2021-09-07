# G Mancini Jul 2021

from sys  import argv, exit
import argparse as argp
import mdtraj as md
import numpy as np
import scipy as sp

import myparse

Parse = argp.ArgumentParser(description='Convert xyz to XTC using a PDB topology')
# template and gaussian settings
Parse.add_argument('-i','--input',help='input xyz trajectory',default=False,action='store')
Parse.add_argument('-t','--topology',help='topology PDB file',default=False,action='store')
Parse.add_argument('-o','--out',help='output (.xtc)',default=False,action='store')
Parse.add_argument('-b','--begin',help='first frame to use',default=-1,type=int)
Parse.add_argument('-e','--end',help='last frame to use',default=-1,type=int)
Parse.add_argument('-d','--deltat',help='timestep between frames',default=False,action='store',type=float)
Parse.add_argument('-I','--index',help='index file name (GMX like)',default=False,action='store')
Parse.add_argument("-m","--molecules",dest="select",default=0,action="store",\
    help="molecules to be written (default all); com nyi; use group numbers")
Parse.add_argument("-x","--shift",action="store",default=False,nargs=3,\
    help="shift coordinates by x,y,z")
Myarg = Parse.parse_args()
print(Myarg)

#SETUP
print ("WARNING: internally working in nm, I/O is angstroms")

if not Myarg.input:
    raise ValueError("Missing input file name")

if not Myarg.topology:
    raise ValueError("Missing topology file name")

if not Myarg.deltat:
    raise ValueError("Missing timestep")

if not Myarg.out:
    raise ValueError("Missing output file name")
else:
    out = Myarg.out + ".xtc"

if Myarg.shift:
    shift = np.array(map(float,Myarg.shift))/10.
else:
    shift = np.zeros(3)
    
#RUN
top = md.load(Myarg.topology)

if Myarg.select != 0:
    atoms = myparse.parse_index(Myarg.index, select)
else:
    atoms = np.arange(top.n_atoms)

natoms = len(atoms)
if natoms != top.n_atoms:
    newtrj = top.atom_slice(atoms)
else:
    newtrj = top
    
myxyz, natoms, nframes = myparse.loadxyz(Myarg.input, atoms)
time = np.asarray([i*Myarg.deltat for i in range(nframes)])

with md.formats.XTCTrajectoryFile(Myarg.out, 'w') as f:
    f.write(myxyz)
 
quit()
