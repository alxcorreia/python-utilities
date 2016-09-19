# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 10:57:31 2015

@author: acorreia
"""

#
# Extracts time series at one or more pressure levels from a .nc (netCDF) file 
#
# Usage: 
# $ python ./getncseries.py /home/acorreia/datafolder/ncfile.nc outputfile lat=-11 lon=-62 level=750
# $ python ./getncseries.py /home/acorreia/datafolder/ncfile.nc outputfile lat=-11 lon=-62 level=750,300
# $ python ./getncseries.py /home/acorreia/datafolder/ncfile.nc outputfile lat=-11 lon=-62 level=ALL
#

import sys
from netCDF4 import Dataset
import numpy as np
import csv
import pandas as pd

file     = sys.argv[1]
outfile0 = sys.argv[2]
mylat    = sys.argv[3]
mylat    = mylat.split('=')[1]
mylon    = sys.argv[4]
mylon    = mylon.split('=')[1]
mylevel  = sys.argv[5]
mylevel  = mylevel.split('=')[1]

print ' '
print 'Reading file: ', file
print '-------------------------'
print ' '
print 'Latitude:  ',mylat+' N'
print 'Longitude: ',mylon+' E'

# Get x and y coordinates from the lat, lon
x=int(90-int(mylat))
y=int(mylon)
if (y < 0):
    y=360+y

# Open datafile for read and get time in a vector
fh = Dataset(file, mode='r')
ref=pd.datetime(1900,1,1,0,0,0,0)
time=np.int64(fh.variables['time'][:])
td=pd.TimedeltaIndex(time,'h')
timearray=td+ref

# Prepare a header for the output file
header0='lat,lon,year,month,day,hour'
year,month,day,hour=timearray.year,timearray.month,timearray.day,timearray.hour
latvec,lonvec=int(mylat)*np.ones(len(time)),int(mylon)*np.ones(len(time))
output0=np.column_stack((latvec,lonvec,year,month,day,hour))

# Get levels in the input file and read variable names
level=fh.variables['level'][:]
myvars = fh.variables.keys()

# Set a vector with the required levels to be processed
mylevel=mylevel.split(',')
if (mylevel[0]=='ALL'):
    mylevel=level.astype(str)


# Iterate through levels and variables to extract each time series
j=0
while (j<len(mylevel)):
    print 'Level:     ',mylevel[j]+' hPa'
    i=0
    while (i < len(myvars)):        
        if (i==0):
            levelvec=int(mylevel[j])*np.ones(len(time))
            header=header0+',level'
            output=np.column_stack((output0,levelvec))
        if ((myvars[i]!='longitude') and (myvars[i]!='latitude') and (myvars[i]!='time') and (myvars[i]!='level')):       
            print 'Extracting: ',myvars[i]
            serieslevels=fh.variables[myvars[i]][:]
            id=np.where(level==np.int(mylevel[j]))[0]
            series=serieslevels[:,id,:,:]
            s=series.squeeze()
            myseries = s[:,x,y]
            output=np.column_stack((output,myseries))
            header=header+','+myvars[i]
        i=i+1
    outfile = outfile0[:-4]+'-'+mylevel[j]+'.csv'
    header=header.split(',')
    csv_out = open(outfile, 'wb')
    dw = csv.DictWriter(csv_out, delimiter=',', fieldnames=header)
    dw.writeheader()
    mywriter = csv.writer(csv_out)
    mywriter.writerows(output)
    csv_out.close()
    print 'Saved file:',outfile
    print '-------------------------'
    j=j+1
fh.close()

exit()
    