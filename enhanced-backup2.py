# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

# This script reads NetCDF GOES-13 files from the input folder, calculates reflectance in the VIS 0.67um (channel 1), reflectance in the IR 3.9um
# (channel 2), and brightness temperature in the 11um (channel 4). The data is combined in a RGB picture as per Rosenfeld and Lensky 1998. 

#Arguments: [1]folderpath with nc files; [2]:output folder

#Usage: 
#  $ python ./enhanced.py /home/acorreia/datafolder /home/acorreia/outputfolder

import os
import sys
import glob
import numpy as np
import datetime
from netCDF4 import Dataset
#import matplotlib as mpl
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
import pylab as py

# Global definitions
# VIS channel definitions
mVIS=0.6118208
kVIS=0.001160
X0=29.0
C=1.255  # For FEB2013
tz=-4  #For Manaus
RadVISUnits='Radiance  [W/(m2 sr um)]'
RefVISUnits='Reflectance'


# For Planck's function:
h = 6.62606896e-34 #   Planck's constant in  J.s
ls = 299792458     #   light speed in  m/s
k = 1.3806504e-23  #   Boltzmann constant in J/K

t0 = 0.75      # Bidirectional transmission function (at 3.9um) between GOES - cloud tops - GOES, by Kaufman and Nakajima 1993 
F0 = 9.68465   # Solar descending irradiance at TOA for 3.9um, by Platnick and Fontenla 2008, in W/[m2 um]

# IR channel definitions
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

def clear_frame(ax=None): 
    if ax is None: 
        ax = plt.gca() 
    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False) 
    for spine in ax.spines.itervalues(): 
        spine.set_visible(False) 
        
        
def linear(inputArray, scale_min=None, scale_max=None):
	"""Performs linear scaling of the input numpy array.

	@type inputArray: numpy array
	@param inputArray: image data array
	@type scale_min: float
	@param scale_min: minimum data value
	@type scale_max: float
	@param scale_max: maximum data value
	@rtype: numpy array
	@return: image data array
	
	"""		
	print "img_scale : linear"
	imageData=np.array(inputArray, copy=True)
	
	if scale_min == None:
		scale_min = imageData.min()
	if scale_max == None:
		scale_max = imageData.max()

	imageData = imageData.clip(min=scale_min, max=scale_max)
	imageData = (imageData -scale_min) / (scale_max - scale_min)
	indices = np.where(imageData < 0)
	imageData[indices] = 0.0
	indices = np.where(imageData > 1)
	imageData[indices] = 1.0
	
	return imageData


def getlatlon ( file ):
    fh = Dataset(file, mode='r')
    lon = fh.variables['lon'][:]
    lat = fh.variables['lat'][:]
    fh.close()
    
    return lat, lon


def getCorners(centers):
    one = centers[:-1,:]
    two = centers[1:,:]
    d1 = (two - one) / 2.
    one = one - d1
    two = two + d1
    stepOne = np.zeros((centers.shape[0] + 1,centers.shape[1]))
    stepOne[:-2,:] = one
    stepOne[-2:,:] = two[-2:,:]
    one = stepOne[:,:-1]
    two = stepOne[:,1:]
    d2 = (two - one) / 2.
    one = one - d2
    two = two + d2
    stepTwo = np.zeros((centers.shape[0] + 1,centers.shape[1] + 1))
    stepTwo[:,:-2] = one
    stepTwo[:,-2:] = two[:,-2:]
    return stepTwo


def Planck ( T, wavelength ):
    # T given in K, wavelength given in um
    #
    # Planck's black body emission: emitted power per unit area per unit solid angle per unit wavelength
    # B=2hc^2/lambda^5 * 1/(exp(hc/(lambda*k*T))-1)  
    lmbda = wavelength * 1e-6   # Converting um to meters
    B = ((2*h*ls*ls/(lmbda**5)) * (1/(np.exp(h*ls/(lmbda*k*T))-1)))/1e6  # Units of W/[m2 sr um]
    return B

def sunzen ( date, time, tz, lat, lon ):
    #function to calculate the sun zenith angle and the Sun-Earth distance
    #
    #input: date object such as 2013-02-25, time object such as 14:00:00
    #input: time zone such as -4 for Manaus
    #input: lat, lon such as -3.0, -60.0 for Manaus
    #output: sun zenith angle in degrees
    #output: Sun-Earth distance in astronomical units

    
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
    

def reflectance ( file, Temp='None' ):
    band=file.split('.')[4]
    procmode='IR'
    if band == 'BAND_01':
        procmode = 'VIS'
    print 'Reading: ', file
    print 'Mode: ', procmode
    fh = Dataset(file, mode='r')
    lon = fh.variables['lon'][:]
    lat = fh.variables['lat'][:]
    data = fh.variables['data'][:]
    fh.close()
    Raw=data/32.   #Convert 16bit to 10bit
    fn=file.split('/')
    
    year=file.split('.')[1]
    julian=file.split('.')[2]
    hhmm=file.split('.')[3]
    h=datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(julian) - 1)
    date=datetime.date(h.year,h.month,h.day)
    time=datetime.time(int(hhmm[:-4]),int(hhmm[2:-2]),int(hhmm[4:]))
    SZA, SED = sunzen(date, time, tz, lat, lon)
    mu0 = np.cos(np.deg2rad(SZA))     
    if procmode == 'VIS':
        RefPre=kVIS*(Raw-X0)
        RefVIS=RefPre * C
        RefVIS=RefVIS*(SED*SED/mu0)  # NOAA eq for reflectance, considering negligible path radiance: can be wrong if aerosols above clouds
        cs = RefVIS
    if procmode == 'IR':
        RadIR=(Raw - bb2)/m2
        RadIR=RadIR*235/358.6362665  # Converting mW/[m2 sr cm-1] to W/[m2 sr um] considering GOES-13 Channel 2 spectral response function
        Emission = Planck(T=Temp,wavelength=3.90)
        RefIR = (RadIR - Emission)/((t0 * F0 * mu0 / (np.pi*SED*SED) - Emission))
        cs = RefIR
    return cs
     
    
def temperature ( file ):
    band=file.split('.')[4]
    procmode='IR'
    print 'Reading: ', file
    print 'Mode: ', procmode
    fh = Dataset(file, mode='r')
    data = fh.variables['data'][:]
    fh.close()
    Raw=data/32.   #Convert 16bit to 10bit
    RadIR=(Raw - bb4)/m4
    Teff4=(C2 * n4)/(np.log(1.+(C1 * np.power(n4,3) /RadIR)))
    TIR=a4+b4*Teff4+g4*( np.power(Teff4,2) )
    
    return TIR

def makeRGB (R, G, B, name='my_rgb_image.png'):
    plt.close("all")
    f=plt.figure(facecolor='Black')
    ax = f.add_axes([0.1,0.1,0.9,0.9])
    img = np.zeros((R.shape[0], R.shape[1], 3), dtype=float) 
    img[:,:,0] = linear(R, scale_min=0)     
    img[:,:,1] = linear(G, scale_min=0, scale_max=np.max(G)/4)
    img[:,:,2] = linear(B, scale_min=248, scale_max=298)   
    py.clf()
    rr=py.imshow(img, aspect=1.8)
    
    fn=name.split('/')[6]
    fn='GOES13-enhanced-'+fn[7:-10]+'png'
    print fn
    
    year=name.split('.')[1]
    julian=name.split('.')[2]
    hhmm=name.split('.')[3]

    h=datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(julian) - 1)    
    month=h.month
    if month < 10:
        month='0'+str(month)
    
    
    label1=' GOES-13      '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
    plt.text(0.145, 0.013,label1,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
    plt.text(0.46, 0.395,'*',transform=ax.transAxes, fontsize=18,color='white')
    plt.text(0.48, 0.395,'Manaus',transform=ax.transAxes, fontsize=8,color='white')
    clear_frame()

    
    py.savefig(fn,dpi=200.,bbox_inches='tight',pad_inches=0.0,transparent=True)
    
    #plt.close("all")
#    f=plt.figure(facecolor='Black')
#    f.set_size_inches(7,6,forward='True')
#    
#    colorTuple=tuple(np.array([img[:,:,0].flatten(),img[:,:,1].flatten(),img[:,:,2].flatten()]).transpose().tolist())
#    m = Basemap(resolution='h',projection='cyl', lon_0=-75,lat_0=0, llcrnrlat=-5.98,llcrnrlon=-63.0,urcrnrlon=-57.0,urcrnrlat=0.0)
#    xi, yi = lon, lat
#    xib=np.hstack((xi, np.zeros((xi.shape[0], 1))))
#    xic=np.vstack((np.zeros((1,xib.shape[1])),xib))
#    yib=np.hstack((yi, np.zeros((yi.shape[0], 1))))
#    yic=np.vstack((np.zeros((1,yib.shape[1])),yib))
#    
#    
#    print "xi ",np.shape(xi)
#    print "xib ", np.shape(xib)
#    print "xic ", np.shape(xic)
#    cornerLats=getCorners(lat);cornerLons=getCorners(lon)
#    xCorners,yCorners=m(cornerLats,cornerLons,inverse=True)
#    m.pcolormesh(xic,yic,np.squeeze(img[:,:,0]),latlon=True,color=colorTuple)
#    m.drawcoastlines(linewidth=0.5)
#    m.drawstates(linewidth=1.0)
#    m.drawcountries(linewidth=1.0)    
#    plt.savefig('teste.png', dpi=200.,bbox_inches='tight',pad_inches=0.0)
#    plt.close("all")
#    
    
    return rr
    

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def rebinned ( X, factor ):
    X=np.squeeze(X)
    newx=int(np.trunc(np.shape(X)[0]/factor))
    newy=int(np.trunc(np.shape(X)[1]/factor))
    oldx=int(newx*factor)
    oldy=int(newy*factor)
    X=X[0:oldx,0:oldy]
    shape = [newx,newy]
    sh = shape[0],X.shape[0]//shape[0],shape[1],X.shape[1]//shape[1]
    r=X.reshape(sh).mean(-1).mean(1)
    return r

        
     
#    # Get some parameters for the Geostationary map projection    
#    f=plt.figure(facecolor='Black')
#    f.set_size_inches(7,6,forward='True')
#    ax = f.add_axes([0,0,1,1])
#    m = Basemap(resolution='h',projection='geos', lon_0=-75,lat_0=0, llcrnrlat=-5.98,llcrnrlon=-63.0,urcrnrlon=-57.0,urcrnrlat=0.0)
#    lon, lat = lons, lats
#    xi, yi = m(lon, lat)
#
#
#    # Add Coastlines, States, and Country Boundaries
#    m.drawcoastlines(linewidth=0.5)
#    m.drawstates(linewidth=1.0)
#    m.drawcountries(linewidth=1.0)
#    
#    # Add legend
#    cbar = plt.colorbar(orientation='horizontal', aspect=47.5,fraction=0.02,pad=0.0)
#    cbar.set_ticks(myticks)
#    cbar.set_ticklabels(mylabels)    
#    cbar.set_label(myunits, size=10)
#    bandnum=band.split('_')[1]
#    month=h.month
#    if procmode == 'IR':
#        procmode = '  '+procmode
#    if month < 10:
#        month='0'+str(month)
#    label1=' GOES-13 Band '+bandnum+'  '+procmode+'           '
#  #  label2='     '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
#    plt.text(0.005, 0.01,label1,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
#    plt.text(0.22, 0.01,label2,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
# 
#    # Output filename 
#   # print "output: ", fn
#   # plt.savefig(fn, dpi=200.,bbox_inches='tight',pad_inches=0.0)
#    plt.close("all")
#    print "done"
#    return
 


#  MAIN LOOP

datadir=sys.argv[1]
curdir=os.getcwd()
outdir=curdir
if len(sys.argv)>2 :
    outdir=sys.argv[2]
os.chdir(outdir)

listch1=glob.glob(datadir+'/*BAND_01.nc')
listch2=glob.glob(datadir+'/*BAND_02.nc')
listch4=glob.glob(datadir+'/*BAND_04.nc')

count=0
while count < len(listch1):
    r063=reflectance(listch1[count])
    r063=rebinned(r063,4)
    temp=temperature(listch4[count])
    r390=reflectance(listch2[count],temp)
    fig=makeRGB(r063,r390,temp,name=listch1[count])
    #label1=' GOES-13      '
    #f=fig.gcf()
    #f.text(100, 100,label1,horizontalalignment='left',verticalalignment='center', fontsize=8,color='white', bbox=dict(facecolor='black' ))
  
    
    #plt.imshow(fig)
    count = count + 1
os.chdir(curdir)
exit()


