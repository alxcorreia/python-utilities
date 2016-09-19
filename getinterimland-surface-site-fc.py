#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 08:22:31 2015

@author: acorreia
"""


# Get 1x1 deg ERA-Interim-Land surface data for a given site and time span

# Usage:   python ./getinterimland-surface-site outputdir XX YYYYMMDD YYYYMMDD
# Example: python ./getinterimland-surface-site /home/acorreia/outputdir AH 20000101 20141231

import os
import sys
import numpy as np
import time
from ecmwfapi import ECMWFDataServer
import ecmwfapi

def getsite (lesite):
#mylat,mylon,tz=-3.0,-60.0,-4 #Manaus
#mylat,mylon,tz= -9.87083333,-56.10388889, -4   #Alta Floresta
#mylat,mylon,tz=-10.76000000,-62.35777778, -4   #Abracos Hill
#mylat,mylon,tz=-10.93388889,-61.85194444, -4   #Ji-Parana
#T0e: ACONVEX (EMBRAPA):  2o53'39.27''S, 59o58'18.46''W [-2.894242, -59.971794]
#T2  : TIWA:              3o8'21.12''S,  60o7'53.52''W  [-3.1392,   -60.131533]
#T3  : ARM (Manacapuru):  3o12'47.82''S, 60o35'55.32''W [-3.213283, -60.5987  ]
#Cuiaba-Miranda  CB  : 15o43'44''S, 56o01'15''W [-15.728889, -56.020833]
#Rio Branco      RB  : 09o57'25''S, 67o52'08''W [-09.956944, -67.868889]
#Campo Grande    CG  : 20o26'16''S, 54o32'16''W [-20.437778, -54.537778]
#Belterra        BT  : 02o38'52''S, 54o57'07''W [-02.647778, -54.951944]

    sitelist = np.array(['AF','AH','JP','MN','CB','RB','CG','BT'])
    failoc=len(np.where(sitelist == lesite)[0])
    if (failoc == 0):
        print 'Failed to recognize the site ',lesite
        print 'Possible sites are: ',sitelist
        print 'Exiting now.'
        sys.exit()

    # Define the site location and time zone
    site,mylat,mylon,tz= 'AF',-10,-56, -4   #Alta Floresta
    if lesite == 'AH':
        site,mylat,mylon,tz='AH',-11,-62, -4   #Abracos Hill
    if lesite == 'JP':
        site,mylat,mylon,tz='JP',-11,-62, -4   #Ji-Parana    
    if lesite == 'MN':
        site,mylat,mylon,tz='MN',-3,-60,-4   #Manaus
    if lesite == 'CB':
        site,mylat,mylon,tz='CB',-16, -56,-4  #Cuiaba
    if lesite == 'RB':
        site,mylat,mylon,tz='RB',-10, -68,-4  #Rio Branco
    if lesite == 'CG':
        site,mylat,mylon,tz='CG',-20, -55,-4  #Campo Grande
    if lesite == 'BT':
        site,mylat,mylon,tz='BT',-3, -55,-4  #Belterra
    return site,mylat,mylon,tz

def retrieval (mydates,myloc,myoutfile,mytype,mytimes,mystep):
    server = ECMWFDataServer()
    server.retrieve({
        'repres'  : "SH",
        'dataset' : "interim_land",     
        'stream'  : "OPER",
        'expver'  : "0002",
        'class'   : "EI",
        'type'    : mytype,
        'levtype' : "SFC",
        'resol'   : "AV", 
        'param'   : "45.128/144.128/146.128/147.128/176.128/177.128/182.128/228.128/243.128",
        'time'    : mytimes,    
        'step'    : mystep,   
        'date'    : mydates,
        'area'    : myloc,
        'grid'    : "1/1",
        'format'  : "netcdf",
        'target'  : myoutfile
        })
    return


    
#  MAIN LOOP


# Usage:   python ./getinterimland-surface-site.py /h/home/acorreia/outputdir XX YYYYMMDD YYYYMMDD
# Example:   python ./getinterimland-surface-site.py /home/acorreia/Desktop/phase/code AH 20000101 20141231

# Read input parameters
outdir    = sys.argv[1]
lesite    = sys.argv[2]
inidate   = sys.argv[3]
enddate   = sys.argv[4]

iniyear = int(inidate[:-4])
endyear = int(enddate[:-4])
max_retries = 3 #number of retries to get data 
os.chdir(outdir)

#mytype="AN"
#mytimes="00/06/12/18"
#mystep="00"

#mytype="AN" or "FC" (Analysis or Forecast)
#mytimes="00/06/12/18"
#mystep="00" or "06" (eg: 00 if analysis, 06 if forecast)
# For instance, many variables like latent heat at the surface are not available in the 18Z analysis.
# In order to get that info we can order 12Z analysis with 06h forecast:
mytype="FC"
mytimes="00"
mystep="00/06/12/18"


# Get site location
[site,mylat,mylon,tz]=getsite(lesite)
print 'Processing site: ',site
print 'Lat/Lon: ',mylat,mylon
print 'Time zone: ',tz

for yy in range(iniyear, endyear+1):
    myloc=str(mylat)+'/'+str(mylon)+'/'+str(mylat)+'/'+str(mylon)
    
    if (yy == iniyear):
        retrievalstartdate=str(yy)+str(inidate[-4:])
    else:
        retrievalstartdate=str(yy)+'0101'
    if (yy == endyear):
        retrievalenddate=str(yy)+str(enddate[-4:])
    else:
        retrievalenddate=str(yy)+'1231'

    myoutfile ='erainterimland-sfc-'+site+'-'+retrievalstartdate+'-'+retrievalenddate+'.nc'
    mydates=retrievalstartdate+'/to/'+retrievalenddate

    i=1
    while (i<=max_retries):
        try:
            retrieval(mydates,myloc,myoutfile,mytype,mytimes,mystep)
            break
        except (RuntimeError, NameError, TypeError, Warning, ValueError, ecmwfapi.api.APIException ):
            print 'Unexpected error:', sys.exc_info()[0]
            print 'This was failed try #',i
            print 'Max number of retries is',max_retries
            print
            if (i<max_retries):
                print 'Waiting 60 seconds for next retry'
                time.sleep(60)
            i=i+1
exit()

