# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

def sunzen ( date, time, tz, lat, lon ):
    #function to calculate the sun zenith angle and the Sun-Earth distance
    #
    #input: date object such as 2013-02-25, time object such as 14:00:00
    #input: time zone such as -4 for Manaus
    #input: lat, lon such as -3.0, -60.0 for Manaus
    #output: sun zenith angle in degrees
    #output: Sun-Earth distance in astronomical units

    import numpy as np
    import datetime
    
    days=date-datetime.date(1900,1,1)
    days=days.days+2
    time=(time.hour/24.)+(time.minute/60./24.)+(time.second/60./60./24.)
    
    # Calculate Julian Day and Julian Century
    JD = days + 2415018.5 + time - tz/24.
    JC = (JD - 2451545)/36525
    
    # Calculate Geometric Mean Anomaly Sun (deg) and Geometric Mean Lon Sun (deg)
    GMAS = 357.52911 + JC*(35999.05029 - (0.0001537)*JC)
    GMLS = np.mod(280.46646 + JC*(36000.76983 + JC*0.0003032),360)
    
    # Calculate Sun Eq of Ctr, Sun True Long (deg),  Sun True Anom (deg) and Sun App Long (deg)
    SEC = np.sin(np.deg2rad(GMAS))*(1.914602-JC*(0.004817+0.000014*JC)) + np.sin(np.deg2rad(2*GMAS))*(0.019993-0.000101*JC) + np.sin(np.deg2rad(3*GMAS))*0.000289
    STL = GMLS + SEC
    STA = GMAS + SEC
    SAL = STL - 0.00569 - 0.00478 * np.sin(np.deg2rad(125.04 - 1934.136*JC))
    
    # Calculate Mean Obliq Ecliptic (deg), Obliq Correction (deg), Eccent Earth Orbit
    MOE = 23 + (26 + ((21.448 - JC*(46.815 + JC*(0.00059 - JC*0.001813))))/60)/60
    OC = MOE + 0.00256*np.cos(np.deg2rad(125.04 - 1934.136*JC))
    EEO = 0.016708634 - JC *(0.000042037 + 0.0000001267*JC)
        
    #Calculate the Equation of Time (min) and the True Solar Time (min)
    vary = np.tan(np.deg2rad(OC/2))*np.tan(np.deg2rad(OC/2))
    EOT = 4*np.rad2deg(vary*np.sin(2*np.deg2rad(GMLS)) - 2*EEO*np.sin(np.deg2rad(GMAS)) + 4*EEO*vary*np.sin(np.deg2rad(GMAS)) * np.cos(2*np.deg2rad(GMLS)) - 0.5*vary*vary*np.sin(4*np.deg2rad(GMLS)) - 1.25*EEO*EEO*np.sin(2*np.deg2rad(GMAS))    )
    TST = np.mod(time*1440 + EOT + 4*lon -60*tz, 1440)
    
    # Calculate the Hour Angle (deg) and the Sun Declination (deg)
    HA = TST/4 - 180
    id=TST<0
    HA[id]=TST[id]/4 + 180
    SD = np.rad2deg(np.arcsin(np.sin(np.deg2rad(OC))*np.sin(np.deg2rad(SAL))))
    
    #Calculate Sun-Earth distance (AUs)
    SED = (1.000001018 * (1 - EEO * EEO))/(1 + EEO*np.cos(np.deg2rad(STA)))    
    
    #Calculate the Sun Zenith Angle (deg)
    SZA = np.rad2deg(np.arccos(np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(SD)) + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(SD)) *np.cos(np.deg2rad(HA))  ))
    
    return SZA, SED
    

def goesimg ( file, type ):
    from netCDF4 import Dataset
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import datetime
    plt.close("all")
    #VIS DEFS
    mVIS=0.6118208
    kVIS=0.001160
    X0=29.0
    C=1.255  # For FEB2013
    RadVISUnits='Radiance  [W/(m2 sr um)]'
    RefVISUnits='Reflectance'
    #IR DEFS
    m2=227.3889
    m3=38.8383
    m4=5.2285
    m6=5.5297
    bb2=68.2167
    bb3=29.1287
    bb4=15.6854
    bb6=16.5892
    RadIRUnits='Radiance  [mW/(m2 sr cm-1)]'
    TUnits='Brightness Temperature  [K]'
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
        #if time < 09 UTC, do not process night time in Manaus - TBD
    print 'Reading: ', file
    print 'Mode: ', procmode
    fh = Dataset(file, mode='r')
    lons = fh.variables['lon'][:]
    lats = fh.variables['lat'][:]
    datas = fh.variables['data'][:]
    fh.close()
    
    print "lon", np.shape(lons)
    print "lat", np.shape(lats)
    print "data", np.shape(datas)
    
    fn=file.split('/')
    fn=fn[len(fn)-1][:-3]+'-'+type+'.png'
    year=fn.split('.')[1]
    julian=fn.split('.')[2]
    hhmm=fn.split('.')[3]
    h=datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(julian) - 1)    
     
    # Get some parameters for the Geostationary map projection    
    f=plt.figure(facecolor='Black')
    f.set_size_inches(7,6,forward='True')
    ax = f.add_axes([0,0,1,1])
    m = Basemap(resolution='h',projection='geos', lon_0=-75,lat_0=0, llcrnrlat=-5.98,llcrnrlon=-63.0,urcrnrlon=-57.0,urcrnrlat=0.0)
    lon, lat = lons, lats
    xi, yi = m(lon, lat)
    
    #Convert 16bit to 10bit
    Raw=datas/32.
    RawUnits='Counts'
    
    # Data processing for VIS channel
    if procmode == 'VIS':
        if type == 'RAW':
            cs = m.pcolormesh(xi,yi,np.squeeze(Raw),cmap='gray',vmin=-50, vmax=1035)
            print np.shape(cs)
            myunits=RawUnits
            myticks=[0,200, 400,600,800, 1000]
            mylabels=["0","200","400","600","800", "1000"]
        if type == 'RAD':
            RadPre=mVIS*(Raw-X0)
            RadVIS=RadPre * C
            cs = m.pcolormesh(xi,yi,np.squeeze(RadVIS),cmap='gray', vmin =-20, vmax=520)
            myunits=RadVISUnits
            myticks=[0,100, 200,300,400, 500] 
            mylabels=["0","100","200","300","400","500"] 
        if type == 'ALB':
            date=datetime.date(h.year,h.month,h.day)
            time=datetime.time(int(hhmm[:-4]),int(hhmm[2:-2]),int(hhmm[4:]))
            tz=-4  #For Manaus
            RefPre=kVIS*(Raw-X0)
            RefVIS=RefPre * C
            SZA, SED = sunzen(date, time, tz, lat, lon)
            RefVIS=RefVIS*(SED*SED/np.cos(np.deg2rad(SZA)))
            cs = m.pcolormesh(xi,yi,np.squeeze(RefVIS),cmap='gray',vmin=-0.05,vmax=1.15)
            myunits=RefVISUnits
            myticks=[0,0.2, 0.4,0.6,0.8, 1.0]
            mylabels=["0.0","0.2","0.4","0.6","0.8", "1.0"]
            
    # Data processing for IR
    if procmode == 'IR':
        if band == 'BAND_02':
            if type == 'RAW':
                cs = m.pcolormesh(xi,yi,np.squeeze(Raw),cmap='gray_r',vmin=-5, vmax=370)
                print np.shape(cs)
                myunits=RawUnits
                myticks=[0,50,100,150,200,250,300,350]
                mylabels=["0","50","100","150","200","250","300", "350"]
            if type == 'RAD':
                RadIR=(Raw - bb2)/m2
                cs = m.pcolormesh(xi,yi,np.squeeze(RadIR),cmap='gray_r', vmin=-0.05, vmax=1.45)
                myunits=RadIRUnits
                myticks=[0.0,0.2,0.4,0.6,0.8,1.0,1.2, 1.4]
                mylabels=["0.0","0.2","0.4","0.6","0.8","1.0","1.2","1.4"]
            if type == 'TMP':
                RadIR=(Raw - bb2)/m2
                Teff2=(C2 * n2)/(np.log(1.+(C1 * np.power(n2,3) /RadIR)))
                TIR=a2+b2*Teff2+g2*( np.power(Teff2,2) )
                cs = m.pcolormesh(xi,yi,np.squeeze(TIR),cmap='gray_r', vmin=245, vmax=305)
                myunits=TUnits
                myticks=[300,290,280,270, 260, 250]
                mylabels=["300","290","280","270","260","250"]
        elif band == 'BAND_03':
            if type == 'RAW':
                cs = m.pcolormesh(xi,yi,np.squeeze(Raw),cmap='gray_r',vmin=-5, vmax=275)
                myunits=RawUnits
                myticks=[0,50,100,150, 200, 250]
                mylabels=["0","50","100","150","200","250"]
            if type == 'RAD':
                RadIR=(Raw - bb3)/m3
                cs = m.pcolormesh(xi,yi,np.squeeze(RadIR),cmap='gray_r', vmin=-0.5,vmax=6.5)
                myunits=RadIRUnits
                myticks=[0,1,2,3,4,5,6]
                mylabels=["0.0","1.0","2.0","3.0","4.0","5.0", "6.0"]
            if type == 'TMP':
                RadIR=(Raw - bb3)/m3
                Teff3=(C2 * n3)/(np.log(1.+(C1 * np.power(n3,3) /RadIR)))
                TIR=a3+b3*Teff3+g3*( np.power(Teff3,2) )
                cs = m.pcolormesh(xi,yi,np.squeeze(TIR),cmap='gray_r', vmin=205, vmax=255)
                myunits=TUnits
                myticks=[250,240, 230, 220, 210, 200]
                mylabels=["250","240","230","220","210", "200"]
        elif band == 'BAND_04':
            if type == 'RAW':
                cs = m.pcolormesh(xi,yi,np.squeeze(Raw),cmap='gray_r',vmin=-10, vmax=620)
                print np.shape(cs)
                myunits=RawUnits
                myticks=[0,100,200, 300,400,500,600]
                mylabels=["0","100","200","300","400","500","600"]
            if type == 'RAD':
                RadIR=(Raw - bb4)/m4
                cs = m.pcolormesh(xi,yi,np.squeeze(RadIR),cmap='gray_r',vmin=-5,vmax=125)
                myunits=RadIRUnits
                myticks=[0,20, 40,60,80, 100,120]
                mylabels=["0","20","40","60","80","100","120"]
            if type == 'TMP':
                RadIR=(Raw - bb4)/m4
                Teff4=(C2 * n4)/(np.log(1.+(C1 * np.power(n4,3) /RadIR)))
                TIR=a4+b4*Teff4+g4*( np.power(Teff4,2) )
                cs = m.pcolormesh(xi,yi,np.squeeze(TIR),cmap='gray_r', vmin=190, vmax=310)
                myunits=TUnits
                myticks=[300,280, 260, 240, 220, 200]
                mylabels=["300","280","260","240","220", "200"]
        elif band == 'BAND_06':
            if type == 'RAW':
                cs = m.pcolormesh(xi,yi,np.squeeze(Raw),cmap='gray_r',vmin=-10, vmax=620)
                myunits=RawUnits
                myticks=[0,100,200, 300,400, 500, 600]
                mylabels=["0","100","200","300","400", "500","600"]
            if type == 'RAD':
                RadIR=(Raw - bb6)/m6
                cs = m.pcolormesh(xi,yi,np.squeeze(RadIR),cmap='gray_r', vmin=-5, vmax=125)
                myunits=RadIRUnits
                myticks=[0,20, 40,60,80, 100,120]
                mylabels=["0","20","40","60","80","100","120"]
            if type == 'TMP':
                RadIR=(Raw - bb6)/m6
                Teff6=(C2 * n6)/(np.log(1.+(C1 * np.power(n6,3) /RadIR)))
                TIR=a6+b6*Teff6+g6*( np.power(Teff6,2) )
                cs = m.pcolormesh(xi,yi,np.squeeze(TIR),cmap='gray_r', vmin=190, vmax=310)
                myunits=TUnits
                myticks=[300,280, 260, 240, 220, 200]
                mylabels=["300","280","260","240","220", "200"]

    # Add Coastlines, States, and Country Boundaries
    m.drawcoastlines(linewidth=0.5)
    m.drawstates(linewidth=1.0)
    m.drawcountries(linewidth=1.0)
    
    # Add legend
    cbar = plt.colorbar(orientation='horizontal', aspect=47.5,fraction=0.02,pad=0.0)
    cbar.set_ticks(myticks)
    cbar.set_ticklabels(mylabels)    
    cbar.set_label(myunits, size=10)
    bandnum=band.split('_')[1]
    month=h.month
    if procmode == 'IR':
        procmode = '  '+procmode
    if month < 10:
        month='0'+str(month)
    label1=' GOES-13 Band '+bandnum+'  '+procmode+'           '
    label2='     '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
    plt.text(0.005, 0.01,label1,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
    plt.text(0.22, 0.01,label2,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
 
    # Output filename 
    print "output: ", fn
    plt.savefig(fn, dpi=200.,bbox_inches='tight',pad_inches=0.0)
    plt.close("all")
    print "done"
    return
 

import os
import sys
import glob


#Arguments: [1]folderpath with nc files; [2]RAW/RAD/ALB/TMP [3]:output folder

#Usage: 
#  $ python ./batch.py /home/acorreia/datafolder RAW /home/acorreia/outputdir

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
    goesimg(list[count],mymode)
    count = count + 1
os.chdir(curdir)
exit()


