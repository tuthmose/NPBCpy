# G Mancini July 2021

#
# collection of functions for parsing (e.g. Gromacs index files)
# and basic I/O

from collections import defaultdict
from itertools import islice

import mdtraj as md
import numpy as np
import scipy as sp

def create_hole(solute, solvent, rsphere, radii, elec, tol, outname):
    """
    create hole using covalent radii and same eq as proxima with
    role of tolerance reversed
    """
    solute_top = solute.topology
    solvent_top = solvent.topology
    solute_atoms = [a.element.symbol for a in solute_top.atoms]
    solvent_atoms = np.asarray([a.element.symbol for a in solvent_top.atoms])
    solvent_natoms = len(solvent_top.select("resid 1"))
    residues = list(range(solvent_top.n_residues))
    remove = list()
    rsphere = rsphere/10.
    solvent.xyz[0] = solvent.xyz[0] - rsphere
    for res in range(len(residues)):
        ratoms = solvent_top.select("resid " + str(res))
        if len(ratoms) == 0: 
            continue
        for jatom in range(solute.n_atoms):
            jelem = solute_atoms[jatom]
            if jelem == "VS": continue
            D = 10.*np.linalg.norm(solvent.xyz[0][ratoms] - solute.xyz[0][jatom], axis=1)
            C = list()
            for i in solvent_atoms[ratoms]:
                if i=="VS" or i=="X" or i=="LP":
                    continue
                c = radii[i] + radii[jelem]-0.07*(elec[i]-elec[jelem])**2 + tol
                C.append(c)
            if np.min(D) <= np.max(np.asarray(C)):
                remove.append(res)
            #break
        #break
    print(remove)
    okres = list(set(residues).difference(remove))
    rlist = "resid " + ' '.join(list(map(str,okres)))
    solvent.restrict_atoms(solvent_top.select(rlist))
    solvent.save(outname)
    return okres

def sphere_radii(atoms, nbins, const_vol, rmin, rmax):
    """
    calculate inner and outer radius of each concentric shell
    """
    NA  = sp.constants.N_A
    nm3_l = 10**(-24)
    radius = (rmax - rmin)
    Vtot = (4.0/3.0)*np.pi*(radius**3)
    if const_vol:
        vol  = Vtot/nbins    
        print("--- Using concentrinc shells of volume: ", vol, " angstroem^3")
        # v_layer_i = 4/3pi(r_i^3-r_i-1^3) = 4/3pi r^3/n
        # r0=0 -> r1^3 = r^3/n; r2^3 = r1^3 + r^3/n=2r^3/n ...
        const = (radius**3)/nbins
        c = 1.0/3.0
        r = [rmin]
        for i in range(1, nbins):
            r.append((i*const)**(c))
            print("--- Layer ",i," radii= ",r[i-1],"->",r[i]," nm, Vol,",csph*(r[i]**3-r[i-1]**3))
        r.append(rmax)
        print("--- Layer ",nbins," radii= ",r[-2],"->",r[-1]," nm, Vol,",csph*(r[-1]**3-r[-2]**3))
    else:
        rlayer = radius/nbins
        csph = 4.0*np.pi/3.0
        r   = [rmin]
        vol = [0.]
        print("--- Using concentric shells of radius :", rlayer, " angstroem")
        for i in range(1, nbins+1):
            r.append(i*rlayer)
            vol.append(csph*(r[i]**3 - r[i-1]**3))
            print("--- Layer ",i," radii= ",r[i-1],"->",r[i]," nm, Vol,",vol[i])
        vol = np.asarray(vol)
        vol = vol[1:]
    r = np.asarray(r)
    ntot = len(atoms)
    nmol = ntot/3
    rho_mean = nmol/Vtot
    rho_moll = rho_mean/(NA*nm3_l)
    print("--- # Atoms and triplets ",ntot, nmol)
    print("--- Average density (for triplets), is: ",rho_mean," molecules/nm3, i.e. ",rho_moll," mol/l")
    return vol, r

def parse_index(ndxfile, mols):
    """
    parse index file; return dict with
    group_name -> atomic indexes in it
    """
    print("--- Parsing index file")
    system = defaultdict(list)
    keys = list()
    ndx = open(ndxfile,"r")
    for line in ndx.readlines():
        rec = line.split()
        if rec[0] == "[" and len(rec)>2:
            k = ' '.join(rec[1:-1])
            keys.append(k)
            current_key = k
        elif rec[0] == "[" and rec[1] == "]":
            keys.append("no group")
            current_key = "no group"
        else:
            atoms = list(map(int,rec))
            system[current_key] = system[current_key]+atoms
    ndx.close()
    st = "--- Found " + str(len(keys)) + " atom groups"
    print(st)
    groups = list()
    atom_groups = list()
    for m, mol in enumerate(mols):
        groups.append(keys[mol])
        st = "--- Selected groups " + str(mols[m]) + ":" + str(groups[m])
        print(st)
        atom_groups.append([int(i)-1 for i in system[groups[m]]])
    return atom_groups

def loadtrj(begin, end, trjname, top):
    """
    load XTC files and checks frame boundaries
    """
    first_frame = int(begin)
    last_frame  = int(end)
    if first_frame!=-1 and last_frame!=-1:
        if last_frame<=first_frame and last_frame!=0:
            print("ERROR: first frame greater or equal than last frame")
            raise ValueError
    traj = md.load(trjname, top=top)
    if last_frame == -1:
        last_frame = traj.n_frames
    elif last_frame > traj.n_frames:
        raise ValueError("ERROR: not enough frames")
    return traj, first_frame, last_frame

def loadxyz(xyzfile, atoms):
    """
    open and read xyz trajectory
    """
    traj = open(xyzfile,"r")
    natoms = int((traj.readline()).split()[0])
    if len(atoms) != natoms:
        select = True
    else:
        select = False
    traj.close()
    nframes = 0
    frames = list()
    traj = open(xyzfile,"r")
    while True:
        lines_gen = islice(traj, natoms+2)
        lines = [i for i in lines_gen][2:]
        if len(lines) != natoms: 
            break
        F = list()
        for nL, L in enumerate(lines):
            record = list(map(float, L.split()[1:]))
            if not select or nL in atoms:
                F.append(record)
        frames.append(F)
        nframes += 1
    traj.close()
    frames = 0.1*np.asarray(frames)
    frames.shape = (nframes, natoms, 3)
    return frames, natoms, nframes
