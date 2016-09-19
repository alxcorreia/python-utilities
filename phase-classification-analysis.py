# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:15:03 2015

@author: acorreia
"""


# Usage: 
#  $ python ./phase-classification-analysis.py /home/acorreia/datafolder /home/acorreia/outputfolder DRY AH
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
outfile=site+'-'+season+'-joined.csv'

flist = sorted(glob.glob(datadir+'/*alldata.csv'))
frame = pd.DataFrame()
list = []
for file in flist:
    str=''
    df = pd.read_csv(file,index_col=None, header=None)
    list.append(df)
frame = pd.concat(list)
csv_out = open(outfile, 'wb')
frame.to_csv(csv_out,index=False,header=None)
csv_out.close()


data     = pd.read_csv(outfile,header=None,names=['vis','ir','t','sza','vza'])
data.fillna(-9999)
rvis,rir,temp=data.vis,data.ir,data.t
rvis,rir,temp=np.float64(rvis),np.float64(rir),np.float64(temp)
rvis,rir,temp=np.array(rvis),np.array(rir),np.array(temp)

iceonly=temp<231
wateronly=temp>274

waterfile=site+'-'+season+'-water.csv'
icefile=site+'-'+season+'-ice.csv'

header='vis,ir,t\n'
outfnw = open(waterfile, 'w')
outfni = open(icefile, 'w')
outfnw.writelines(header) 
outfni.writelines(header) 
outfnw.close()
outfni.close()


wrvis=rvis[wateronly]
wrir=rir[wateronly]
wtemp=temp[wateronly]

irvis=rvis[iceonly]
irir=rir[iceonly]
itemp=temp[iceonly]


wrvis=wrvis.reshape((-1,1))
wrir=wrir.reshape((-1,1))
wtemp=wtemp.reshape((-1,1))


irvis=irvis.reshape((-1,1))
irir=irir.reshape((-1,1))
itemp=itemp.reshape((-1,1))

water=np.column_stack((wrvis,wrir,wtemp))
ice=np.column_stack((irvis,irir,itemp))


frame = pd.DataFrame(data=water)
csv_out = open(waterfile, 'a+b')
frame.to_csv(csv_out,index=False,header=None)
csv_out.close()

frame = pd.DataFrame(data=ice)
csv_out = open(icefile, 'a+b')
frame.to_csv(csv_out,index=False,header=None)
csv_out.close()

os.chdir(curdir)
exit()

