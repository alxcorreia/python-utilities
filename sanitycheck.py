# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 10:24:41 2015

@author: acorreia
"""
import os
import glob
import numpy as np

def sanity(datadir):
    flag=False
    missinglist=np.array([])
    os.chdir(datadir)
    listch1=sorted(glob.glob('*BAND_01.nc')) #BAND1 in 4km resolution
    listch2=sorted(glob.glob('*BAND_02.nc'))
    listch4=sorted(glob.glob('*BAND_04.nc'))
    l1,l2,l4=len(listch1),len(listch2),len(listch4)
    if ( l1 != l2 ) or ( l1 != l4 ):
        flag=True
    maxlist=listch1
    dlist1=listch2
    dlist2=listch4
    if np.argmax([l1,l2,l4])==1:
        maxlist=listch2
        dlist1=listch1
        dlist2=listch4
    if np.argmax([l1,l2,l4])==2:
        maxlist=listch4
        dlist1=listch1
        dlist2=listch2
    i,j,k=0,0,0
    while i<len(maxlist):
        sat,year,julian,hhmmss=maxlist[i].split('.')[0],maxlist[i].split('.')[1],maxlist[i].split('.')[2],maxlist[i].split('.')[3]
        if (j<len(dlist1)) and (j>=0):
            sat1,year1,julian1,hhmmss1=dlist1[j].split('.')[0],dlist1[j].split('.')[1],dlist1[j].split('.')[2],dlist1[j].split('.')[3]
            if (sat1!=sat) or (year1!=year) or (julian1!=julian) or (hhmmss1!=hhmmss):
                flag=True
                prefix=maxlist[i]
                prefix=prefix.split('BAND')[0]
                posfix=dlist1[j]
                posfix=posfix.split('BAND')[1]
                missinglist=np.append(missinglist,prefix+'BAND'+posfix)
                j=j-1            
        else:
            flag=True
            missinglist=np.append(missinglist,maxlist[i])
            
        
        if (k<len(dlist2)) and (k>=0):
            sat2,year2,julian2,hhmmss2=dlist2[k].split('.')[0],dlist2[k].split('.')[1],dlist2[k].split('.')[2],dlist2[k].split('.')[3]
            if (sat2!=sat) or (year2!=year) or (julian2!=julian) or (hhmmss2!=hhmmss):
                flag=True
                prefix=maxlist[i]
                prefix=prefix.split('BAND')[0]
                posfix=dlist2[j]
                posfix=posfix.split('BAND')[1]
                missinglist=np.append(missinglist,prefix+'BAND'+posfix)
                k=k-1
            
        else:
            flag=True
            missinglist=np.append(missinglist,maxlist[i])
               
        
        
        i=i+1
        j=j+1
        k=k+1
        

    return flag,missinglist
    
    
datadir='/home/acorreia/Desktop/phase/data/alldata2005-2007'

[myflag,mylist]=sanity(datadir)

print myflag
print mylist
    
