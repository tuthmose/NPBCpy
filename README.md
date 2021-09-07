# Collection of scripts for the analysis of MD trajectories under Non Periodic (spherical) Boundary Conditions (NPBC)

G Mancini aug 2021

Atom / molecule selections is done using GROMACS index files. Trajectory I/O is in xtc.
Dimensions are specified by each script but usually are angstroems, kJ/mol (even if
internally distances are calculated in nm).
A PDB topology file is needed for most calculations

## dependencies
- numpy
- scipy
- mdtraj
- cython

## files

setup: to compile Cython functions

### libraries
- npbc\_io: includes functions for parsing index files, generating data for concentric layers
- npbc\_analysis: includes functions to read trajectories and calculate the needed quantities frame per frame
- npbc\_cy: includes some basic math/numpy based functions and Cython versions of slow analysis (e. g. ADFs)

### scripts
- angular\_distribution: calculate ADF for three atom groups, in given distance interval
- calc\_glob: (awk script) calculate mean field VdW potential from data formatted by get\_glob.awk (see DOI: 10.1021/acs.jctc.0c00454) 
- continuos\_hbonds: calculate HB strenght as a continuous function in [0,1] based on distance and angular cutoff (see DOI: 10.1002/jcc.24683)
- density\_layer: density in concentric layers of constant radius or volume
- distance: distance between a atom group and its nth nearest neighbour from another group as a function of time
- gen\_connM: (awk script) generate connectivity matrix for G16 input file for a molecule type with a central atom
- get\_data\_l402: parse Gaussian16 log file to extract data from a classical MDMC run in Link402 (see: 10.1021/acs.jctc.0c00454) 
- get\_glob: see calc\_glob
- nmol: get the number of atoms/molecules within a a given distance interval from another one
- orientation: calculate the angle of a vector calculated from an atom triplet wrt to the radius of the spherical box
- pair\_correlation: calculate the radial distribution function (g(r)) normalizing correctly for spherical NPBC
- pbd2com: (awk script) generate a G16 input file from a PDB and other template files
- sphere.vmd: example of selection of a sphere in VMD

### notebooks
- Solvate.ipynb: example of carving a hole in a spherical box using covalent radii adjusted as in (DOI: acs.jcim.0c00076)
