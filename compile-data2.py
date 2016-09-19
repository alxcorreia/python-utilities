# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

#Usage: 
#  $ python ./compile-data.py /home/acorreia/datafolder /home/acorreia/outputfolder

import os
import sys
import glob

datadir=sys.argv[1]
curdir=os.getcwd()
outdir=curdir
if len(sys.argv)>2 :
    outdir=sys.argv[2]
os.chdir(outdir)

flist=sorted(glob.glob(datadir+'/*alldata.csv'))

filename0=flist[0]
filename0=filename0.split('/')[len(filename0.split('/'))-1]
outfile0=filename0[:-24]+'dailydata.csv'

sat0    = filename0.split('.')[0]
year0   = filename0.split('.')[1]
julian0 = filename0.split('.')[2]

fh0 = open (flist[0],'r')
data0=fh0.read()
fh0.close()

fhout = open(outfile0, 'w')
fhout.write(data0)
fhout.close()

        
data=''
count=1
while count < len(flist):
    filename=flist[count]
    filename=filename.split('/')[len(filename.split('/'))-1]
    outfile=filename[:-24]+'dailydata.csv'

    sat    = filename.split('.')[0]
    year   = filename.split('.')[1]
    julian = filename.split('.')[2]

    fh =  open (flist[count],'r')
    data=fh.read()
    fh.close()
    
    if ((sat == sat0) & (year == year0) & (julian == julian0)):
        fhout = open(outfile, 'a')
        fhout.write(data)
        fhout.close()
    else :
        fhout = open(outfile, 'w')
        fhout.write(data)
        fhout.close()
        sat0,year0,julian0=sat,year,julian
    count = count + 1
os.chdir(curdir)
exit()


