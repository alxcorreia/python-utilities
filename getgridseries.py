# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 10:57:31 2015

@author: acorreia
"""

#
# Extracts time series at one or more pressure levels from a .nc (netCDF) file at a single grid point
#
# Usage: 
# $ python ./getgridseries.py /home/acorreia/datafolder/ncfile.nc outputfile level=750
# $ python ./getgridseries.py /home/acorreia/datafolder/ncfile.nc outputfile level=750,300
# $ python ./getgridseries.py /home/acorreia/datafolder/ncfile.nc outputfile level=ALL
#

import sys
from netCDF4 import Dataset
import numpy as np
import csv
import pandas as pd

file    = sys.argv[1]
outfile = sys.argv[2]
mylevel = sys.argv[3]
mylevel = mylevel.split('=')[1]

print ' '
print 'Reading file: ', file
print '-------------------------'
print ' '

# Open datafile for read and get time in a vector
fh = Dataset(file, mode='r')
ref=pd.datetime(1900,1,1,0,0,0,0)
time=np.int64(fh.variables['time'][:])
td=pd.TimedeltaIndex(time,'h')
timearray=td+ref

#Get lat lon
mylon=float(fh.variables['longitude'][:][0])
mylat=float(fh.variables['latitude'][:][0])
if (mylon>180):
    mylon=mylon-360
print 'Latitude:  ',str(mylat)+' N'
print 'Longitude: ',str(mylon)+' E'

# Get levels in the input file and read variable names
level=fh.variables['level'][:]
myvars = fh.variables.keys()

# Set a vector with the required levels to be processed
mylevel=mylevel.split(',')
if (mylevel[0]=='ALL'):
    mylevel=level.astype(str)
#i=0
#while i<len(mylevel):
#    if i==0:
#        levelvec=int(mylevel[i])*np.ones(len(time))
#    levelvec=np.concatenate((levelvec,int(mylevel[i])*np.ones(len(time))))
    
# Prepare a header for the output file
header0='lat,lon,year,month,day,hour,level'
year,month,day,hour=timearray.year,timearray.month,timearray.day,timearray.hour
#latvec,lonvec=int(mylat)*np.ones(len(time)*len(mylevel)),int(mylon)*np.ones(len(time)*len(mylevel))
#output0=np.column_stack((latvec,lonvec,year,month,day,hour,levelvec))    
latvec,lonvec=int(mylat)*np.ones(len(time)),int(mylon)*np.ones(len(time))
output0=np.column_stack((latvec,lonvec,year,month,day,hour))    
    
    



# Iterate through levels and variables to extract each time series



# Iterate through levels and variables to extract each time series
j=0
while (j<len(mylevel)):
    print 'Level:     ',mylevel[j]+' hPa'
    if (j==0):
        header=header0
    i=0
    while (i < len(myvars)):        
        if (i==0):
            #print i,j
            #levelvec=int(mylevel[j])*np.ones(len(time))
            #header=header0+',level'
            levelvec=int(mylevel[j])*np.ones(len(time))
            output=np.column_stack((output0,levelvec))
            
        if ((myvars[i]!='longitude') and (myvars[i]!='latitude') and (myvars[i]!='time') and (myvars[i]!='level')): 
            #print i,j
            print 'Extracting: ',myvars[i]
            serieslevels=fh.variables[myvars[i]][:]
            
            id=np.where(level==np.int(mylevel[j]))[0]
            series=serieslevels[:,id,:,:]
            myseries=series.squeeze()            
            output=np.column_stack((output,myseries))
            if (j==0):
                header=header+','+myvars[i]
        i=i+1
    if (j==0):
        outlist=output
    else:
        outlist=np.row_stack((outlist,output))
    print '-------------------------'
    j=j+1
fh.close()

header=header.split(',')

csv_out = open(outfile, 'wb')
dw = csv.DictWriter(csv_out, delimiter=',', fieldnames=header)
dw.writeheader()
mywriter = csv.writer(csv_out)
mywriter.writerows(outlist)
csv_out.close()
print 'Saved file:',outfile



exit()
  

