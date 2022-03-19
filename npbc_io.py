# G Mancini July 2021

#
# collection of functions for parsing (e.g. Gromacs index files)
# and basic I/O

# TODO eliminate atom counting and use residues and mdtraj.topology

from collections import defaultdict
from itertools import islice

import mdtraj as md
import numpy as np
import re
import scipy as sp

def reglog():
    """
    define regula expressions to parse input file
    gdv/j19 Link402 ver March 22
    """
    REGX = dict()
    REGX['get_data'] = \
    re.compile(r'Statistics\swill\sbe(\s|\w)+steps\.\n(.*)Simulation\sconcluded',re.DOTALL)
    REGX['get_mean_Ep'] = \
    re.compile(r'Average\spotential\senergy\s:\s+(.*)\s\+\/\-\s+(.*)\n')
    REGX['get_mean_T']  = \
    re.compile(r'Average\stemperature\s+:\s+([0-9\.]+)\s\+\/\-\s+([0-9\.]+)')
    REGX['get_acc_ratio'] = re.compile(r'Total\sacceptance\sratio\s+:\s+(.*)\s\(')
    REGX['get_atomic_ratio'] = re.compile(r'Atomic\stranslation\s.*\s+:\s+(.*)\s\(')
    REGX['get_transl_ratio'] = re.compile(r'Fragment\stranslation.*\s+:\s+(.*)\s\(')
    REGX['get_rota_ratio'] = re.compile(r'Fragment\srotation.*\s+:\s+(.*)\s\(')
    return REGX

def getmean(regx, log):
    """
    using the regular expressions defined in reglog, return averages
    of energies and acceptance ratios
    log is a read gdv log file
    """
    data = dict()
    try:
        data['mean_temp'] = float(regx['get_mean_T'].search(log).group(1))
        data['sd_temp'] =  float(regx['get_mean_T'].search(log).group(2))
    except:
        data['mean_temp'] = None
        data['sd_temp'] =  None
    try:
        data['mean_Ep'] = float((regx['get_mean_Ep'].search(log).group(1)).replace('D','E'))
        data['sd_Ep'] =  float((regx['get_mean_Ep'].search(log).group(2)).replace('D','E'))
    except:
        data['mean_Ep'] = None
        data['sd_Ep'] =  None
    try:
        data['mean_ratio'] = float(regx['get_acc_ratio'].search(log).group(1))
    except:
        data['mean_ratio'] = None
    try:
        data['mean_at_ratio'] = float(regx['get_atomic_ratio'].search(log).group(1))
    except:
        data['mean_at_ratio'] = None
    try:
        data['mean_tr_ratio'] = float(regx['get_transl_ratio'].search(log).group(1))
    except:
        data['mean_tr_ratio'] = None
    try:
        data['mean_rt_ratio'] = float(regx['get_rota_ratio'].search(log).group(1)) 
    except:
        data['mean_rt_ratio'] = None
    return data

def getsimul(regx, averages, log, use_ts=False, SepMD=False, SepMC=False):
    """
    get data point from simulations; average values are repeated for 
    fast plot with xmgrace/GNUplot
    when printing all data in one array missing data in MC lines is
    filled with 0s
    """
    ts_fields = list(range(1,7))
    sp_fields = [0] + list(range(2,7))
    toskip = ('Trajectory', 'MMDT2A', 'MMDT2F', 'Step', '(#)')
    data = regx['get_data'].search(log).group(2)
    lines = data.split("\n")
    if SepMD:
        mddata = list()
    else:
        mddata = None
    if SepMC:
        mcdata = list()
    else:
        mcdata = None
    alldata = list()
    for line in lines:
        record = line.split()
        if len(record)==0 or record[0] in toskip:
            continue
        if SepMD and len(record)==7:
            if use_ts:
                dataline = [float(record[f]) for f in ts_fields[:-1]] + \
                    [float(record[ts_fields[-1]].replace('D','E'))] + \
                    [averages['mean_temp'], averages['mean_Ep']]
            else:
                dataline = [float(record[f]) for f in sp_fields[:-1]] + \
                    [float(record[sp_fields[-1]].replace('D','E'))] + \
                    [averages['mean_temp'], averages['mean_Ep']]
            mddata.append(dataline)
        if SepMC and len(record) == 5:
            acc = lambda x: 0 if x=='F' else 1
            dataline = [float(record[0]),float(record[-1]),averages['mean_Ep'],acc(record[3]),\
                 averages['mean_ratio']]
            mcdata.append(dataline)
        if len(record) == 5:
            alldata.append([float(record[0]),float(record[-1]),averages['mean_Ep'],0.0])
        elif len(record) == 7:
            alldata.append([float(record[0]),float(record[4]), \
                float(record[-1].replace('D','E')),averages['mean_Ep']])
    return np.asarray(alldata), np.asarray(mddata), np.asarray(mcdata)
    
    

def sortmol(reference, target, natom1, natom2, coords, nearest=True, metric="euclidean"):
    """
    sort atoms in target with respect to their distance to
    reference and rewrite target and natom consecutive atoms
    in that order
    all target atoms must be already ordered and begin after
    reference+natom1
    """
    R = coords[reference][np.newaxis]
    D = sp.spatial.distance.cdist(R, coords[target], metric=metric)[0]
    order = np.argsort(D)
    #if nearest == False:
    #    order = order[::-1]
    print(D[30])
    D = D[order]
    neworder = list(range(natom1))
    for o in order:
        neworder = neworder + list(range((o+1)*3,(o+1)*3+natom2))
    return neworder, order, D

def create_hole(solute, solvent, rsphere_solute, rsphere_solvent, radii, elec, tol, outname):
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
    rsphere_solute = rsphere_solute/10.
    rsphere_solvent = rsphere_solvent/10.
    solute.xyz[0] = solute.xyz[0] - rsphere_solute
    solvent.xyz[0] = solvent.xyz[0] - rsphere_solvent
    for res in range(len(residues)):
        ratoms = solvent_top.select("resid " + str(res))
        if len(ratoms) == 0: 
            continue
        for jatom in range(solute.n_atoms):
            jelem = solute_atoms[jatom]
            if jelem == "VS" or jelem == "LP":
                continue
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
    okres = list(set(residues).difference(remove))
    rlist = "resid " + ' '.join(list(map(str,okres)))
    solvent.restrict_atoms(solvent_top.select(rlist))
    solvent.save(outname)
    return okres

def sphere_radii(atoms, natoms, nbins, const_vol, rmin, rmax):
    """
    calculate inner and outer radius of each concentric shell
    """
    NA  = sp.constants.N_A
    nm3_l = 10**(-24)
    radius = (rmax - rmin)
    csph = 4.0*np.pi/3.0
    Vtot = csph*(radius**3)
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
    nmol = ntot/natoms
    rho_mean = nmol/Vtot
    rho_moll = rho_mean/(NA*nm3_l)
    print("--- # Atoms and molecules",ntot, nmol)
    print("--- Average density (for molecules), is: ",rho_mean," molecules/nm3, i.e. ",rho_moll," mol/l")
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
