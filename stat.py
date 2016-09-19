# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""
#def flatten(lista):
#    flat = [val for sublist in lista for val in sublist]
#    return flat
#            


def goesstat ( file, type ):
    from netCDF4 import Dataset
    import numpy as np
    #import pylab as P
    import matplotlib.pyplot as plt
    #from mpl_toolkits.basemap import Basemap
    import datetime
    #plt.close("all")
    #VIS DEFS
    mVIS=0.6118208
    kVIS=0.001160
    X0=29.0
    C=1.255  # For FEB2013
    RadVISUnits='W/(m2 sr um)'
    RefVISUnits='Unitless'
    #IR DEFS
    m2=227.3889
    m3=38.8383
    m4=5.2285
    m6=5.5297
    bb2=68.2167
    bb3=29.1287
    bb4=15.6854
    bb6=16.5892
    RadIRUnits='mW/(m2 sr cm-1)'
    TUnits='K'
    C1=1.191066e-5   #  mW/(m2 sr cm-4)
    C2=1.438833      #  K/(cm-1)
    n2=2561.7421
    n3=1522.5182
    n4=937.23449
    n6=749.82589
    a2=-1.4755462
    a3=-4.1556932
    a4=-0.52227011
    a6=-0.16089410
    b2=1.0028656
    b3=1.0142082
    b4=1.0023802
    b6=1.0006896
    g2=-5.8203946e-7
    g3=-8.0255086e-6
    g4=-2.0798856e-6
    g6=-3.9853774e-7
    band=file.split('.')[4]
    procmode='IR'
    if band == 'BAND_01':
        procmode = 'VIS'
        #if time < 09 UTC, nao processar (madrugada)
    print 'Reading: ', file
    print 'Mode: ', procmode
    fh = Dataset(file, mode='r')
    #lons = fh.variables['lon'][:]
    #lats = fh.variables['lat'][:]
    datas = fh.variables['data'][:]
    fh.close()




    # Get some parameters for the Stereographic Projection    
    #f=plt.figure(facecolor='Black')
    #f.set_size_inches(7,6,forward='True')
    #ax = f.add_axes([0,0,1,1])
    #m = Basemap(resolution='h',projection='geos', lon_0=-75,lat_0=0, llcrnrlat=-5.98,llcrnrlon=-63.0,urcrnrlon=-57.0,urcrnrlat=0.0)
    #lon, lat = lons, lats
    #xi, yi = m(lon, lat)


    #Convert 16bit to 10bit
    Raw=datas/32.
    RawUnits='Unitless'
    # Data processing for VIS
    if procmode == 'VIS':
        if type == 'RAW':
            #cs = m.pcolormesh(xi,yi,np.squeeze(Raw),cmap='gray',vmin=-150, vmax=850)
            ds=Raw
            myunits=RawUnits
        if type == 'RAD':
            RadPre=mVIS*(Raw-X0)
            RadVIS=RadPre * C
            ds=RadVIS
            #cs = m.pcolormesh(xi,yi,np.squeeze(RadVIS),cmap='gray')
            myunits=RadVISUnits
        if type == 'ALB':
            RefPre=kVIS*(Raw-X0)
            RefVIS=RefPre * C
            ds=RefVIS
            #cs = m.pcolormesh(xi,yi,np.squeeze(RefVIS),cmap='gray')
            myunits=RefVISUnits
    # Data processing for IR
    if procmode == 'IR':
        if band == 'BAND_02':
            if type == 'RAD':
                RadIR=(Raw - bb2)/m2
                ds=RadIR
            if type == 'TMP':
                RadIR=(Raw - bb2)/m2
                Teff2=(C2 * n2)/(np.log(1.+(C1 * np.power(n2,3) /RadIR)))
                TIR=a2+b2*Teff2+g2*( np.power(Teff2,2) )
                ds=TIR
        elif band == 'BAND_03':
            if type == 'RAD':
                RadIR=(Raw - bb3)/m3
                ds=RadIR
            if type == 'TMP':
                RadIR=(Raw - bb3)/m3
                Teff3=(C2 * n3)/(np.log(1.+(C1 * np.power(n3,3) /RadIR)))
                TIR=a3+b3*Teff3+g3*( np.power(Teff3,2) )
                ds=TIR
        elif band == 'BAND_04':
            if type == 'RAD':
                RadIR=(Raw - bb4)/m4
                ds=RadIR
            if type == 'TMP':
                RadIR=(Raw - bb4)/m4
                Teff4=(C2 * n4)/(np.log(1.+(C1 * np.power(n4,3) /RadIR)))
                TIR=a4+b4*Teff4+g4*( np.power(Teff4,2) )
                ds=TIR
        elif band == 'BAND_06':
            if type == 'RAD':
                RadIR=(Raw - bb6)/m6
                ds=RadIR
            if type == 'TMP':
                RadIR=(Raw - bb6)/m6
                Teff6=(C2 * n6)/(np.log(1.+(C1 * np.power(n6,3) /RadIR)))
                TIR=a6+b6*Teff6+g6*( np.power(Teff6,2) )
                ds=TIR
        if type == 'RAW':
            #cs = m.pcolormesh(xi,yi,np.squeeze(Raw),cmap='gray_r',vmin=-80, vmax=850)
            ds=Raw
            myunits=RawUnits
        if type == 'RAD':
            #cs = m.pcolormesh(xi,yi,np.squeeze(RadIR),cmap='gray_r')
            myunits=RadIRUnits
        if type == 'TMP':
            #cs = m.pcolormesh(xi,yi,np.squeeze(TIR),cmap='gray_r')
            myunits=TUnits
    # Add Coastlines, States, and Country Boundaries
    #m.drawcoastlines(linewidth=0.5)
    #m.drawstates(linewidth=1.0)
    #m.drawcountries(linewidth=1.0)
    
    # Add legend
    
    #plt.text(0.5, 0.5,'matplotlib', horizontalalignment='center', verticalalignment='center', bbox=dict(facecolor='black' ))
    # Add Colorbar
    #cbar = m.colorbar(cs, location='bottom', pad="10%")
    #cbar.set_label(myunits)
    fn=file.split('/')
    fn=fn[len(fn)-1][:-2]
    year=fn.split('.')[1]
    julian=fn.split('.')[2]
    hhmm=fn.split('.')[3]
    h=datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(julian) - 1)
    bandnum=band.split('_')[1]
    month=h.month
    if procmode == 'IR':
        procmode = '  '+procmode
    if month < 10:
        month='0'+str(month)
    #label1='     GOES-13 Band '+bandnum+'  '+procmode+'               '
    #label2='          '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
    #plt.text(0.0, 0.018,label1,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
    #plt.text(0.237, 0.018,label2,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
    #' GOES-13 Band 1 VIS       25-FEB-2013 161514 UTC ' -60 -3
    print "output: ", fn
    print "type: ", type
    print "units: ", myunits
    print "min: ", np.min(ds)
    print "mean: ", np.mean(ds)
    print "median: ", np.median(ds)
    print "max: ", np.max(ds)
    print "N:  ",len(np.ndarray.flatten(ds))
    
    minmin=np.min(ds)
    maxmax=np.max(ds)
    #minmin, maxmax = 0, 1.500 

    
    hhh, bins = np.histogram(ds, bins=100,range=(minmin,maxmax))
    
    width = bins[1] - bins[0]
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hhh,align='center', width=width)
    plt.show()

    #plt.show()
    #plt.savefig(fn, dpi=200.,bbox_inches='tight',pad_inches=0.0)
    #plt.close("all")
    print "=========================================="
    
    
    
    return
    



import os
import sys
import glob


#Arguments: [1]folderpath with nc files; [2]RAW/RAD/ALB/TMP [3]:output folder

#Usage: 
#  $ python ./stat.py /home/acorreia/datafolder ALB /home/acorreia/outputdir

datadir=sys.argv[1]
mymode=sys.argv[2]
curdir=os.getcwd()
outdir=curdir
if len(sys.argv)>3 :
    outdir=sys.argv[3]

os.chdir(outdir)

list=glob.glob(datadir+'/*.nc')

count=0
while count < len(list):
    goesstat(list[count],mymode)
    count = count + 1
os.chdir(curdir)
exit()


