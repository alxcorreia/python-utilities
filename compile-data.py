# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

#Usage: 
#  $ python ./compile-data.py /home/acorreia/datafolder /home/acorreia/outputfolder

import os
import sys
import glob

datadir=sys.argv[1]
curdir=os.getcwd()
outdir=curdir
if len(sys.argv)>2 :
    outdir=sys.argv[2]
os.chdir(outdir)

listvis=sorted(glob.glob(datadir+'/*cloudreflectance063.csv'))
listird=sorted(glob.glob(datadir+'/*cloudreflectance390.csv'))
listtmp=sorted(glob.glob(datadir+'/*cloudtemperature.csv'))

count=0
while count < len(listvis):
    fh1, fh2, fh3 = open(listvis[count], 'r'), open(listird[count], 'r'), open(listtmp[count], 'r')
    outfile = listvis[count]
    name=listvis[count]
    outfile=name.split('/')[len(name.split('/'))-1]
    outfile=outfile[:-23]+'alldata.csv'
    fh4 = open(outfile, 'w')
    for line in fh1:
        d1=line[:-2]
        d2=fh2.readline()[:-2]
        d3=fh3.readline()[:-2]
        outstr=d1+','+d2+','+d3+'\n'
        fh4.writelines(outstr)
    fh1.close()
    fh2.close()
    fh3.close()
    fh4.close()
    count = count + 1
os.chdir(curdir)
exit()


