#!/usr/bin/env python3
# G Mancini June 2020


import argparse as argp
import math
import re
import sys

Parse = argp.ArgumentParser(description='Get trajectory information (energy, temperature, etc) from Link402.j12 MD/MC log file')
# default: print only number of steps and potential energy
Parse.add_argument('--MC',help='print montecarlo data in separate file (<log_file_prefix.mc.dat>',\
    default=False,action='store_true')
Parse.add_argument('--MD',help='print montecarlo data in separate file (<log_file_prefix.MD.dat>',\
    default=False,action='store_true')
Parse.add_argument('--ts',help='label MD points with time step',default=False,action='store_true')
Parse.add_argument('--base_line',help='print base line with average energies',\
    default=False,action='store_true')
Parse.add_argument('--prefix',help='output files name prefix (defaulti: same as input)',\
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

# G16 regular expressions
get_data = re.compile(r'Statistics\swill\sbe(\s|\w)+steps\.\n(.*)Simulation\sconcluded',re.DOTALL)
get_mean_Ep = re.compile(r'Average\spotential\senergy\s:\s+(.*)\s\+\/\-\s+(.*)\n')
get_mean_T  = re.compile(r'Average\stemperature\s+:\s+([0-9\.]+)\s\+\/\-\s+([0-9\.]+)')
get_acc_ratio = re.compile(r'Total\sacceptance\sratio\s+:\s+(.*)\s\(')
get_atomic_ratio = re.compile(r'Atomic\stranslation\s.*\s+:\s+(.*)\s\(')
get_transl_ratio = re.compile(r'Fragment\stranslation.*\s+:\s+(.*)\s\(')
get_rota_ratio = re.compile(r'Fragment\srotation.*\s+:\s+(.*)\s\(')

# open and parse file
print("+++ Parsing file ", Myarg.log)
log_file = open(Myarg.log,'r')
log = log_file.read()

# base line
if Myarg.base_line:
    try:
        mean_temp = get_mean_T.search(log).group(1)
        sd_temp =  get_mean_T.search(log).group(2)
    except:
        print("No average temperature read")
        mean_temp = ""
        sd_temp =  ""
    try:
        mean_Ep = get_mean_Ep.search(log).group(1)
        sd_Ep =  get_mean_Ep.search(log).group(2)
    except:
        print("No average potential read")
        mean_Ep = ""
        sd_Ep =  ""
    try:
        mean_ratio = get_acc_ratio.search(log).group(1)
    except:
        print("No acceptance ratio read")
        mean_ratio = ""
    try:
        mean_at_ratio = get_atomic_ratio.search(log).group(1)
    except:
        print("No atomic acceptance ratio read")
        mean_at_ratio = ""
    try:
        mean_tr_ratio = get_transl_ratio.search(log).group(1)
    except:
        print("No translation acceptance ratio read")
        mean_tr_ratio = ""
    try:
        mean_rt_ratio = get_rota_ratio.search(log).group(1) 
    except:
        print("No rotation acceptance ratio read")
        mean_rt_ratio = ""

ts_fields = list(range(1,7))
sp_fields = [0] + list(range(2,7))
toskip = ('Trajectory', 'MMDT2A', 'MMDT2F', 'Step', '(#)')

data = get_data.search(log).group(2)
lines = data.split("\n")
for line in lines:
    record = line.split()
    if len(record)==0 or record[0] in toskip:
        continue
    if Myarg.MD and len(record) == 7:
        if Myarg.ts:
            dataline = [record[f] for f in ts_fields] + [mean_temp, mean_Ep]
        else:
            dataline = [record[f] for f in sp_fields] + [mean_temp, mean_Ep]
        outline = "".join([d + " " for d in dataline]) +"\n"
        mdfile.write(outline)
    if Myarg.MC and len(record) == 5:
        acc = lambda x: '0' if x=='F' else '1'
        dataline = record[0] + " " + record[-1] + " " + mean_Ep\
            + " " + acc(record[3]) + " " + mean_ratio +"\n"
        mcfile.write(dataline)
    if len(record) == 5:
        outfile.write(record[0] + " " + record[-1] + " " + mean_Ep + "\n")
    elif len(record) == 7:
        outfile.write(record[0] + " " + record[4] + " " + mean_Ep + "\n")

log_file.close()
outfile.close()
if Myarg.MC:
    mcfile.close()
if Myarg.MD:
    mdfile.close()

quit()
