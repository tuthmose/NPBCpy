#!/usr/bin/env python3
# G Mancini June 2020

import argparse as argp
import numpy as np
import gzip
import os
import re
import shutil
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
Parse.add_argument('--zip',help='log file is gzipped G16 log file to parse',\
    default=False,action='store_true')
Myarg = Parse.parse_args()
print("+++ Options set")
print(Myarg)

#
def read_uncompressed(compressed_log):
    with gzip.open(compressed_log, 'rb') as f_in:
        with open('file.txt', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    f_in = open('file.txt', "r")
    log = f_in.read()   
    f_in.close()
    os.remove("file.txt")
    return log

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
if Myarg.MD:
    mdname = prefix + ".MD.dat"
outname = prefix + ".dat"

log_regex = npbc_io.reglog()

# open and parse file
if Myarg.zip:
    log = read_uncompressed(Myarg.log)
else:
    print("+++ Parsing file ", Myarg.log)
    log_file = open(Myarg.log,'r')
    log = log_file.read()

if Myarg.base_line:
    average_data = npbc_io.getmean(log_regex, log)

alldata, mddata, mcdata = npbc_io.getsimul(log_regex, average_data, log, Myarg.ts, Myarg.MD,\
    Myarg.MC)

if Myarg.MC:
    np.savetxt(mcname, mcdata)
if Myarg.MD:
    np.savetxt(mdname, mddata)
np.savetxt(outname, alldata)

quit()
