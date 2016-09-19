# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:15:03 2015

@author: acorreia
"""

# 
# Usage: 
#  $ python ./coldwater-hotice.py /home/acorreia/datafolder /home/acorreia/outputfolder outputfilename DRY AH
#

import os
import sys
import glob
import pandas as pd
import numpy as np

datadir=sys.argv[1]
curdir=os.getcwd()
outdir=sys.argv[2]
outfile=sys.argv[3]
season=sys.argv[4]
site=sys.argv[5]
os.chdir(outdir)

header='sat,site,season,year,julian,hhmmss,boxsize,N_cw,min_temp_cw,max_temp_cw,median_temp_cw,mean_temp_cw,stdev_temp_cw,N_hi,min_temp_hi,max_temp_hi,median_temp_hi,mean_temp_hi,stdev_temp_hi\n'
outfn = open(outfile, 'w')
outfn.writelines(header) 

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

flist = sorted(glob.glob(datadir+'/*_stat.csv'))

count=0
while count < len(flist):
    filename = flist[count]
    filename = filename.split('/')[len(filename.split('/'))-1]
    sat      = filename.split('.')[0]
    
    year     = filename.split('.')[1]
    
    julian   = filename.split('.')[2]
    last     = filename.split('.')[3]
    hhmmss   = last.split('_')[0]
    boxsize  = last.split('_')[1]
    
    data     = pd.read_csv(flist[count],header=1,names=['asat','asite','aseason','ayear','ajulian','ahhmmss','aboxsize','avis','air','atemp','aratio','aphase'])
   
    watertemp,icetemp = '-9999','-9999'
        
    idx=(data.aboxsize=='box12') & (data.aphase=='water')
    watertemp=data[idx][['atemp']].values
    watertemp=watertemp.flatten()
    
    idx=(data.aboxsize=='box12') & (data.aphase=='ice')
    icetemp=data[idx][['atemp']].values
    icetemp=icetemp.flatten()
    
    nwater=np.shape(watertemp)[0]
    nice=np.shape(icetemp)[0]
    
    coldwater,coldwatermed,coldwatersd,coldwatern,coldwatermin,coldwatermax='-9999','-9999','-9999','-9999','-9999','-9999'
    hotice,hoticemed,hoticesd,hoticen,hoticemin,hoticemax='-9999','-9999','-9999','-9999','-9999','-9999'
    
    if (nwater>0):
        sortedwater = np.sort(watertemp,axis=None)
        waterset = sortedwater[0:int(coldwaterpercent*len(sortedwater))+1]
        if (np.shape(waterset)[0]==0):
            waterset=sortedwater
        coldwater=str(np.average(waterset))
        coldwatermed=str(np.median(waterset))
        coldwatermin=str(min(waterset))
        coldwatermax=str(max(waterset))
        coldwatersd=str(np.std(waterset))
        coldwatern=str(len(waterset))
        
    if (nice>0):
        sortedice   = np.sort(icetemp,axis=None)
        iceset   = sortedice[int((1-hoticepercent)*len(sortedice))-1:len(sortedice)-1]
        if (np.shape(iceset)[0]==0):
            iceset=sortedice
        hotice=str(np.average(iceset))
        hoticemed=str(np.median(iceset))
        hoticemin=str(min(iceset))
        hoticemax=str(max(iceset))
        hoticesd=str(np.std(iceset))
        hoticen=str(len(iceset))
        
    outstring = sat+','+site+','+season+','+str(year)+','+str(julian)+','+str(hhmmss)+','+boxsize+','
    outstring = outstring + str(coldwatern)+','+coldwatermin+','+coldwatermax+','+coldwatermed+','+coldwater+','+coldwatersd+','
    outstring = outstring + hoticen+','+hoticemin+','+hoticemax+','+hoticemed+','+hotice+','+hoticesd+'\n'
    
    outfn.writelines(outstring)

    count = count + 1
outfn.close()
os.chdir(curdir)
exit()



        
