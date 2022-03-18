#!/usr/bin/env python3
# G Mancini June 2020

import argparse as argp
import math
import re
import sys

import npbc_io

Parse = argp.ArgumentParser(description='Get trajectory information (energy, temperature, etc) from Link402.j12 MD/MC log file')
# default: print only number of steps and potential energy
Parse.add_argument('--MC',help='print montecarlo data in separate file (<log_file_prefix.mc.dat>',\
    default=False,action='store_true')
Parse.add_argument('--MD',help='print montecarlo data in separate file (<log_file_prefix.MD.dat>',\
    default=False,action='store_true')
Parse.add_argument('--ts',help='label MD points with time step',default=False,action='store_true')
Parse.add_argument('--base_line',help='print base line with average energies',\
    default=False,action='store_true')
Parse.add_argument('--prefix',help='output files name prefix (default: same as input)',\
    default=None,action='store')
Parse.add_argument('--log',help='G16 log file to parse',\
    default=None,action='store')
Myarg = Parse.parse_args()
print("+++ Options set")
print(Myarg)

# set file name and open output files
prf = re.compile('(.*)\.log')
if Myarg.prefix is not None:
    prefix = Myarg.prefix
else:
    try:
        prefix = prf.search(Myarg.log).group(1)
    except:
        raise ValueError('Please use .log suffix for gaussian output files')
if Myarg.MC:
    mcname = prefix + ".mc.dat"
    mcfile = open(mcname, 'w')
if Myarg.MD:
    mdname = prefix + ".MD.dat"
    mdfile = open(mdname, 'w')
outname = prefix + ".dat"
outfile = open(outname, "w")

log_regex = npbc_io.getlog()

# open and parse file
print("+++ Parsing file ", Myarg.log)
log_file = open(Myarg.log,'r')
log = log_file.read()

if Myarg.base_line:
    average_data = npbc_io.getmean(log_regx, log)

alldata, mddata, mcdata = npbc_io.get_simul(log_regex, average_data, log, Myarg.ts, Myarg.MD,\
    Myarg.MC)

        outline = "".join([d + " " for d in dataline]) +"\n"
log_file.close()
outfile.close()
if Myarg.MC:
    mcfile.close()
if Myarg.MD:
    mdfile.close()

quit()
