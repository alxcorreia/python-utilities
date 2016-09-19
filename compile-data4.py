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

outfile='goesdata.csv'


fhout = open(outfile, 'a')
#header='site,year,julian,hour,min,sec,sat,RefVIS,RefIR,Temp\n'
#fhout.write(header)

#site='AH'
#site='AF'
site='JP'
        
data=''
count=0
while count < len(flist):
    filename=flist[count]
    filename=filename.split('/')[len(filename.split('/'))-1]
    #outfile=filename[:-24]+'dailydata.csv'

    sat    = filename.split('.')[0]
    year   = filename.split('.')[1]
    julian = filename.split('.')[2]
    hhmmss = filename.split('.')[3]
    hhmmss = hhmmss.split('_')[0]
    hh,mm,ss=hhmmss[:-4],hhmmss[2:-2],hhmmss[4:]


    with open(flist[count],'r') as fh:
        for line in fh:
            RefVIS,RefIR,Temp=line.split(',')[0],line.split(',')[1],line.split(',')[2]
            outstr=site+','+year+','+julian+','+hh+','+mm+','+ss+','+sat+','+RefVIS+','+RefIR+','+Temp
            fhout.write(outstr)
    count = count + 1
fhout.close()
os.chdir(curdir)
exit()


