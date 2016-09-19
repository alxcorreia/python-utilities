# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:15:03 2015

@author: acorreia
"""

# Feb19/2015 This script reads CSV files in a given folder and calculates statistics about them. The input files are supposed to have 3 numeric columns 
# containing real numbers. The output files will contain all the original data and statistics about each column.  
#
# Usage: 
#  $ python ./compile-stats2.py /home/acorreia/datafolder /home/acorreia/outputfolder DRY AH
#

import os
import sys
import glob
import pandas as pd
import numpy as np


datadir=sys.argv[1]
curdir=os.getcwd()
outdir=sys.argv[2]
os.chdir(outdir)
season=sys.argv[3]
site=sys.argv[4]


# Define VIS/IR thresholds to detect water, ice or mixed phase pixels
cloudVISthreshold  = 0.20   # VIS reflectance above this limit indicates cloud pixel 
cloudTempthreshold = 283.0  # Temp in Kelvin below this limit indicates cloud pixel
icedry   = 17.2399  # VIS/IR above this limit indicates ice pixel
icewet   = 20.2301  # VIS/IR above this limit indicates ice píxel
waterdry =  4.6801  # VIS/IR below this limit indicates water pixel
waterwet =  7.3172  # VIS/IR below this limit indicates water pixel

# Define percentage of pixels to be considered as 'hot ice', e.g. 0.10 means 10% hottest ice pixels
hoticepercent = 0.10

# Define percentage of pixels to be considered as 'cold water', e.g. 0.10 means 10% coldest water pixels
coldwaterpercent = 0.10



flist = sorted(glob.glob(datadir+'/*.csv'))


count=0
while count < len(flist):
    filename = flist[count]
    filename = filename.split('/')[len(filename.split('/'))-1]
    outfile  = filename[:-4]+'_stat.csv'
    sat      = filename.split('.')[0]
    year     = filename.split('.')[1]
    julian   = filename.split('.')[2]
    last     = filename.split('.')[3]
    hhmmss   = last.split('_')[0]
    boxsize  = last.split('_')[1]
    data     = pd.read_csv(flist[count],header=None,names=['vis','ir','t'])
    rvis,rir,temp=data.vis,data.ir,data.t
    
    nonzero=(rir>0.0)
    ratio=rvis[nonzero]/rir[nonzero]
    icelim = icewet
    waterlim = waterwet
    if season == 'DRY':
        icelim = icedry
        waterlim = waterdry
   
    
    nv=str(len(rvis))
    vismin=str(min(rvis))
    vismax=str(max(rvis))
    vismed=str(np.median(rvis))
    visavg=str(np.average(rvis))
    visstd=str(np.std(rvis))

    ni=str(len(rir))
    irmin=str(min(rir))
    irmax=str(max(rir))
    irmed=str(np.median(rir))
    iravg=str(np.average(rir))
    irstd=str(np.std(rir))    
    
    nt=str(len(temp))
    tmin=str(min(temp))
    tmax=str(max(temp))
    tmed=str(np.median(temp))
    tavg=str(np.average(temp))
    tstd=str(np.std(temp))
    
    minratio,minphase='---','---'
    if (min(rir)>0.0):
        minphase='mix'
        minratio=float(vismin)/float(irmin)
        if ( (minratio>icelim) and (float(vismin)>cloudVISthreshold) and (float(tmin)<cloudTempthreshold) ):
            minphase='ice'
        if ( (minratio<waterlim) and (float(vismin)>cloudVISthreshold) and (float(tmin)<cloudTempthreshold) ):
            minphase='water'
    minratio=str(minratio)

    maxratio,maxphase='---','---'
    if (max(rir)>0.0):
        maxphase='mix'
        maxratio=float(vismax)/float(irmax)
        if ( (maxratio>icelim) and (float(vismax)>cloudVISthreshold) and (float(tmax)<cloudTempthreshold) ):
            maxphase='ice'
        if ( (maxratio<waterlim) and (float(vismax)>cloudVISthreshold) and (float(tmax)<cloudTempthreshold) ):
            maxphase='water'
    maxratio=str(maxratio)

    medratio,medphase='---','---'
    if (np.median(rir)>0.0):
        medphase='mix'
        medratio=float(vismed)/float(irmed)
        if ( (medratio>icelim) and (float(vismed)>cloudVISthreshold) and (float(tmed)<cloudTempthreshold) ):
            medphase='ice'
        if ( (medratio<waterlim) and (float(vismed)>cloudVISthreshold) and (float(tmed)<cloudTempthreshold) ):
            medphase='water'
    medratio=str(medratio)    
    
    avgratio,avgphase='---','---'
    if (np.average(rir)>0.0):
        avgphase='mix'
        avgratio=float(visavg)/float(iravg)
        if ( (avgratio>icelim) and (float(visavg)>cloudVISthreshold) and (float(tavg)<cloudTempthreshold) ):
            avgphase='ice'
        if ( (avgratio<waterlim) and (float(visavg)>cloudVISthreshold) and (float(tavg)<cloudTempthreshold) ):
            avgphase='water'
    avgratio=str(avgratio)
    
    
    if os.path.isfile(outfile):
        outfn = open(outfile, 'a')
    else:
        header='sat,site,season,year,julian,hhmmss,boxsize,r_vis,r_ir,temp(K),VIS/IR_ratio,phase\n'
        outfn = open(outfile, 'w')
        outfn.writelines(header)    

    
    with open(flist[count],'r') as fh0:
        for line in fh0:
            rv,ri,tp=line.split(',')[0],line.split(',')[1],line.split(',')[2]
            tp=float(tp)
            ratio,phase='---','---'
            if (float(ri)>0.0):
                phase='mix'
                ratio=float(rv)/float(ri)
                if ( (ratio>icelim) and (float(rv)>cloudVISthreshold) and (float(tp)<cloudTempthreshold) ):
                    phase='ice'
                if ( (ratio<waterlim) and (float(rv)>cloudVISthreshold) and (float(tp)<cloudTempthreshold) ):
                    phase='water'
            ratio=str(ratio)
            tp=str(tp)
            outstr=sat+','+site+','+season+','+year+','+julian+','+hhmmss+','+boxsize+','+rv+','+ri+','+tp+','+ratio+','+phase+'\n'
                        
            outfn.writelines(outstr)


    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',N,'+nv+','+ni+','+nt+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',MIN,'+vismin+','+irmin+','+tmin+','+minratio+','+minphase+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',MAX,'+vismax+','+irmax+','+tmax+','+maxratio+','+maxphase+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',AVG,'+visavg+','+iravg+','+tavg+','+avgratio+','+avgphase+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',MED,'+vismed+','+irmed+','+tmed+','+medratio+','+medphase+'\n'
    outfn.writelines(outstr)
    outstr = sat+','+site+','+season+','+year+','+julian+','+hhmmss+',STD,'+visstd+','+irstd+','+tstd+'\n'
    outfn.writelines(outstr+'\n')

    count = count + 1
outfn.close()
os.chdir(curdir)
exit()

