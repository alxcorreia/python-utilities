# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 10:57:31 2015

@author: acorreia
"""

#
# Extracts time series at a single pressure level from a .nc (netCDF) file 
#
# Usage: 
# $ python ./getncseries-singlelevel.py /home/acorreia/datafolder/ncfile.nc outputfile lat=-11 lon=-62 level=750
#

import sys
from netCDF4 import Dataset
import numpy as np
import csv
import pandas as pd

file    = sys.argv[1]
outfile = sys.argv[2]
mylat   = sys.argv[3]
mylat   = mylat.split('=')[1]
mylon   = sys.argv[4]
mylon   = mylon.split('=')[1]
mylevel = sys.argv[5]
mylevel = mylevel.split('=')[1]
outfile = outfile[:-4]+'-'+mylevel+'.csv'

print ' '
print 'Reading file: ', file
print '-------------------------'
print ' '
print 'Latitude:  ',mylat+' N'
print 'Longitude: ',mylon+' E'
print 'Level:     ',mylevel+' hPa'

# Get x and y coordinates from the lat, lon
x=int(90-int(mylat))
y=int(mylon)
if (y < 0):
    y=360+y

# Open datafile for read
fh = Dataset(file, mode='r')

# Get time in a vector
ref=pd.datetime(1900,1,1,0,0,0,0)
time=np.int64(fh.variables['time'][:])
td=pd.TimedeltaIndex(time,'h')
timearray=td+ref

# Prepare a header for the output file
header='year,month,day,hour'
year,month,day,hour=timearray.year,timearray.month,timearray.day,timearray.hour
output=np.column_stack((year,month,day,hour))

# Find what levels are in the input file and locate the one we want
level=fh.variables['level'][:]
id=np.where(level==np.int(mylevel))[0]

# Find all the variables in the input file
myvars = fh.variables.keys()

# Iterate through the variable list to extract their time series
i=0
while (i < len(myvars)):    
    if ((myvars[i]!='longitude') and (myvars[i]!='latitude') and (myvars[i]!='time') and (myvars[i]!='level')):
        print 'Extracting: ',myvars[i]
        serieslevels=fh.variables[myvars[i]][:]
        series=serieslevels[:,id,:,:]
        s=series.squeeze()
        myseries = s[:,x,y]
        output=np.column_stack((output,myseries))
        header=header+','+myvars[i]
    i=i+1

# Prepare the final header 
header=header.split(',')

# Write the output in a CSV file
csv_out = open(outfile, 'wb')
dw = csv.DictWriter(csv_out, delimiter=',', fieldnames=header)
dw.writeheader()
mywriter = csv.writer(csv_out)
mywriter.writerows(output)
csv_out.close()
fh.close()

print 'Saved file:',outfile
print '-------------------------'
exit()
    