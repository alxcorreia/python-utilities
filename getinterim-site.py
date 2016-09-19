#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 08:22:31 2015

@author: acorreia
"""


# Get 1x1 deg ERA-Interim data for a given site and time span

# Usage:   python ./getinterim-site outputdir XX YYYYMMDD YYYYMMDD
# Example: python ./getinterim-site /home/acorreia/outputdir AH 20000101 20141231

import os
import sys
import numpy as np
from ecmwfapi import ECMWFDataServer

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

def retrieval (mydates,myloc,myoutfile):
    server = ECMWFDataServer()
    server.retrieve({
        'repres'  : "SH",
        'dataset' : "interim",
        'expver'  : "0001",
        'stream'  : "OPER",
        'class'   : "EI",
        'type'    : "AN",
        'levtype' : "PL",
        'resol'   : "AV",
        'levelist': "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
        'param'   : "z/q/t/w/d/vo/o3/u/v/pv/ciwc/cc/r/clwc",
        'time'    : "00/06/12/18",    
        'step'    : "0",   
        'date'    : mydates,
        'area'    : myloc,
        'grid'    : "1/1",
        'format'  : "netcdf",
        'target'  : myoutfile
        })
    return

#  MAIN LOOP


# Usage:   python ./getinterim-site.py /h/home/acorreia/outputdir XX YYYYMMDD YYYYMMDD
# Example:   python ./getinterim-site.py /home/acorreia/Desktop/phase/code AH 20000101 20141231

# Read input parameters
outdir    = sys.argv[1]
lesite    = sys.argv[2]
inidate   = sys.argv[3]
enddate   = sys.argv[4]

iniyear = int(inidate[:-4])
endyear = int(enddate[:-4])

os.chdir(outdir)

# Get site location
[site,mylat,mylon,tz]=getsite(lesite)
print 'Processing site: ',site
print 'Lat/Lon: ',mylat,mylon
print 'Time zone: ',tz

for yy in range(iniyear, endyear+1):
    myloc=str(mylat)+'/'+str(mylon)+'/'+str(mylat)+'/'+str(mylon)
    myoutfile ='erainterim-'+site+'-'+str(yy)+'.nc'
    if (yy == iniyear):
        retrievalstartdate=str(yy)+str(inidate[-4:])
    else:
        retrievalstartdate=str(yy)+'0101'
    if (yy == endyear):
        retrievalenddate=str(yy)+str(enddate[-4:])
    else:
        retrievalenddate=str(yy)+'1231'

    mydates=retrievalstartdate+'/to/'+retrievalenddate
    print mydates,myloc,myoutfile
    retrieval(mydates,myloc,myoutfile)
exit()

