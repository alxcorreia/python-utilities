# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

#Usage: 
#  $ python ./compile-data6.py /home/acorreia/datafolder /home/acorreia/outputfolder    outputfilename

import os
import sys
import glob
import pandas as pd
import numpy as np


datadir=sys.argv[1]
curdir=os.getcwd()
outdir=curdir
if len(sys.argv)>2 :
    outdir=sys.argv[2]
os.chdir(outdir)
outfile=sys.argv[3]


if os.path.isfile(outfile):
    outfn = open(outfile, 'a')
else:
    header='sat,site,season,year,julian,boxsize,type,N,min_temp(K),max_temp(K),median_temp(K),mean_temp(K),stdev_temp(K)\n'
    outfn = open(outfile, 'w')
    outfn.writelines(header)    

flist = sorted(glob.glob(datadir+'/*-temperatures.csv'))

count=0
while count < len(flist):
    filename = flist[count]
    filename = filename.split('/')[len(filename.split('/'))-1]
    sat      = filename.split('.')[0]
    year     = filename.split('.')[1]
    julian   = filename.split('.')[2]
    last     = filename.split('.')[3]
    boxsize  = last.split('-')[0]
    season   = last.split('-')[1]
    site     = last.split('-')[2]
    ttype    = last.split('-')[4]
    data     = pd.read_csv(flist[count],header=None,names=['t'])
    temp = data.t
    n=str(len(temp))
    tmin=str(min(temp))
    tmax=str(max(temp))
    tmedian=str(np.median(temp))
    tavg = str(np.average(temp))
    tstd = str(np.std(temp))
    outstr = sat+','+site+','+season+','+year+','+julian+','+boxsize+','+ttype+','+n+','+tmin+','+tmax+','+tmedian+','+tavg+','+tstd+'\n'
    outfn.writelines(outstr)
    count = count + 1
outfn.close()
os.chdir(curdir)
exit()


