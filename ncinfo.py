# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 10:57:31 2015

@author: acorreia
"""

#
# Shows information about all variables in a .nc (netCDF) file 
#
# Usage: 
# $ python ./ncinfo.py /home/acorreia/datafolder/ncfile.nc 
#

import sys
from netCDF4 import Dataset
file = sys.argv[1]
print ' '
print 'Showing info on: ', file
print '-------------------------'
print ' '
fh = Dataset(file, mode='r')
for obj in fh.variables.values():
    print obj
    print ' '
print '-------------------------'
exit()
    