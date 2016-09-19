# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/acorreia/.spyder2/.temp.py
"""



import os
import sys
import glob
import calendar

startyear=2000
endyear=2014
datadir='/home/acorreia/Desktop/phase/data/18z/CG/001'
outdir='/home/acorreia/Desktop/phase'



os.chdir(datadir)
f1list = sorted(glob.glob(datadir+'/*BAND_01.nc'))
f2list = sorted(glob.glob(datadir+'/*BAND_02.nc'))
f4list = sorted(glob.glob(datadir+'/*BAND_04.nc'))


print 'Number of CH1 files: ',len(f1list)
print 'Number of CH2 files: ',len(f2list)
print 'Number of CH4 files: ',len(f4list)

#exit()

year=startyear
while year <= endyear:
    leap=0
    if (calendar.isleap(year)):
        leap=1
    for i in range (1,366+leap):
        if i<10:
            julian='00'+str(i)
        elif i<100:
            julian='0'+str(i)
        else:
            julian=str(i)
        
        print 'Year and julian day: ',str(year),julian
        
        fn1='/goes??.'+str(year)+'.'+julian+'.??????.BAND_01.nc'
        fn2='/goes??.'+str(year)+'.'+julian+'.??????.BAND_02.nc'
        fn4='/goes??.'+str(year)+'.'+julian+'.??????.BAND_04.nc'
        
        
        f1list = sorted(glob.glob(datadir+fn1))
        f2list = sorted(glob.glob(datadir+fn2))
        f4list = sorted(glob.glob(datadir+fn4))
        
        numok=0
        if (len(f1list)==len(f2list)) and (len(f1list)==len(f4list)):
            numok=1

        dataok=1
        if ((len(f1list)==0) or (len(f2list)==0) or (len(f4list)==0)):
            dataok=0

        if numok==0:
            print 'Number of CH1, CH2, CH4 files is different: ',len(f1list),len(f2list),len(f4list)
            print
        if dataok==0:
            print 'No data for CH1, CH2 or CH4: ',len(f1list),len(f2list),len(f4list)
            print
        if ((numok==1)and(dataok==1)):
            print '.'
    year=year+1
        


exit()

        
