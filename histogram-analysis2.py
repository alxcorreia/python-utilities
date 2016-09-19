# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:22:06 2014

@author: acorreia
"""

# This script reads NetCDF GOES-10, -12 or -13 files from the input folder, calculates reflectance in the VIS 0.67um (channel 1), reflectance in the IR 3.9um
# (channel 2), and brightness temperature in the 11um (channel 4). The data is combined in a RGB picture as per Rosenfeld and Lensky 1998. 

#Arguments: [1]folderpath with nc files; [2]:output folder

#Usage: 
#  $ python ./histogram-analysis.py /home/acorreia/datafolder /home/acorreia/outputfolder

import os
import sys
import glob
import numpy as np
import datetime
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pylab as py
import datetime as dt

# Global definitions

# For Planck's function:
h = 6.62606896e-34 #   Planck's constant in  J.s
ls = 299792458     #   light speed in  m/s
k = 1.3806504e-23  #   Boltzmann constant in J/K

t0 = 1  #0.75      # Bidirectional transmission function (at 3.9um) between GOES - cloud tops - GOES, t0=0.75 by Kaufman and Nakajima 1993 
F0 = 9.68465   # Solar descending irradiance at TOA for 3.9um, by Platnick and Fontenla 2008, in W/[m2 um]

# VIS channel definitions
RadVISUnits='Radiance  [W/(m2 sr um)]'
RefVISUnits='Reflectance'
X0=29.0

mVIS10=0.5582154
kVIS10=0.001110
C10=1.711 #Valid for GOES-10 from Jun2005-Dec2007

mVIS12=0.5771
kVIS12=0.001141
C12=1.0 #Valid for GOES-12 until Aug2007

mVIS13=0.6118208
kVIS13=0.001160
C13=1.255  # Valid for GOES-13 for FEB2013



# IR channel definitions
RadIRUnits='Radiance  [mW/(m2 sr cm-1)]'
TUnits='Brightness Temperature  [K]'
C1=1.191066e-5   #  mW/(m2 sr cm-4)
C2=1.438833      #  K/(cm-1)


m210=227.3889
m310=38.8383
m410=5.2285
m610=5.0273
bb210=68.2167
bb310=29.1287
bb410=15.6854
bb610=15.3332
n210=2552.9845
n310=1486.2212
n410=936.10260
n610=830.88473
a210=-0.63343694
a310=-0.66500842
a410=-0.36939128
a610=-0.32763317
b210=1.0013206
b310=1.0017857
b410=1.0017466
b610=1.0014057
g210=-4.2038547e-7
g310=-7.3885254e-7
g410=-1.4981835e-6
g610=-9.5563444e-7
factor10=164.6/252.267139805


m212=227.3889
m312=38.8383
m412=5.2285
m612=5.5297
bb212=68.2167
bb312=29.1287
bb412=15.6854
bb612=16.5892
n212=2562.45
n312=1536.43
n412=933.21
n612=751.91
a212=-0.727744
a312=-5.278975
a412=-0.534982
a612=-0.177244
b212=1.002131
b312=1.016476
b412=1.002693
b612=1.000138
g212=-1.173898e-6
g312=-7.754348e-6
g412=-2.667092e-6
g612=1.163496e-6
factor12=174.6/264.8402490292



m213=227.3889
m313=38.8383
m413=5.2285
m613=5.5297
bb213=68.2167
bb313=29.1287
bb413=15.6854
bb613=16.5892
n213=2561.7421
n313=1522.5182
n413=937.23449
n613=749.82589
a213=-1.4755462
a313=-4.1556932
a413=-0.52227011
a613=-0.16089410
b213=1.0028656
b313=1.0142082
b413=1.0023802
b613=1.0006896
g213=-5.8203946e-7
g313=-8.0255086e-6
g413=-2.0798856e-6
g613=-3.9853774e-7
factor13=232.4/354.1531901968
# RAD [W/m2/sr/um] = factor * RAD [mW/m2/sr/cm-1]




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
	#print "img_scale : linear"
	imageData=np.array(inputArray, copy=True)
	
	if scale_min == None:
		scale_min = imageData.min()
	if scale_max == None:
		scale_max = imageData.max()

	imageData = imageData.clip(min=scale_min, max=scale_max)
	factor=(scale_max - scale_min)
     	if factor == 0:
		factor=1
     
	imageData = (imageData -scale_min) / factor
	indices = np.where(imageData < 0)
	imageData[indices] = 0.0
	indices = np.where(imageData > 1)
	imageData[indices] = 1.0
	
	return imageData

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
    
def savereflectanceimg (data,file):
    procmode, vmin, vmax, myticks, mylabels = ' IR', -0.03, 1.03, [0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], ["0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"]
    bandnum=file.split('_')[1].split('.')[0]
    data=np.squeeze(data)
    img = np.zeros((data.shape[0], data.shape[1], 2), dtype=float) 
    img = linear(data, scale_min=0, scale_max=1) 
    if bandnum == '01':
        procmode, vmin, vmax, myticks, mylabels = 'VIS', -0.03 , 1.03, [0,0.1,0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], ["0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"]
    f=plt.figure(facecolor='Black')
    f.set_size_inches(7,6,forward='True')
    ax = f.add_axes([0.02,0.1,0.9,0.9])
    py.clf()
    py.imshow(np.squeeze(img), aspect=1.8, vmin=vmin, vmax=vmax,cmap='spectral')
    cbar = plt.colorbar(orientation='horizontal', aspect=46.5,fraction=0.02,pad=0.0)
    cbar.set_ticks(myticks)
    cbar.set_ticklabels(mylabels)    
    cbar.set_label('Reflectance', size=10)
    fn=file.split('/')
    fn=fn[len(fn)-1][:-3]+'-ALB.png'
    year=fn.split('.')[1]
    julian=fn.split('.')[2]
    hhmm=fn.split('.')[3]
    h=datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(julian) - 1)
    month=h.month
    if month < 10:
        month='0'+str(month)
    label1=' GOES-13 Band '+bandnum+'  '+procmode+'     '
    label2='     '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
    plt.text(0.198, 0.028,label1+label2,transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
    plt.text(0.55, 0.41,'*',transform=ax.transAxes, fontsize=18,color='white')
    plt.text(0.57, 0.41,'Manaus',transform=ax.transAxes, fontsize=8,color='white')
    clear_frame()
    print "output: ", fn
    plt.savefig(fn, dpi=200.,bbox_inches='tight',pad_inches=0.0)
    plt.close("all")
    return


def reflectance ( file, Temp='None',save='No' ):
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
        #RefVIS=RefVIS*(SED*SED/mu0)  # This line was commented out due to unrealistic values for reflectance in the afternoon. NOAA eq for reflectance, considering negligible path radiance: can be wrong if aerosols above clouds
        cs = RefVIS
    if procmode == 'IR':
        RadIR=(Raw - bb2)/m2
        RadIR=RadIR*factor  # Converting mW/[m2 sr cm-1] to W/[m2 sr um] considering GOES-13 Channel 2 spectral response function
        Emission = Planck(T=Temp,wavelength=3.90)
        #RefIR = (RadIR - Emission)/((t0 * F0 * mu0 / (np.pi*SED*SED) - Emission))
        RefIR = (RadIR - Emission)/((t0 * F0 * 1.0 / (np.pi*1.0) - Emission))
        cs = RefIR
    if save == 'Yes':
        savereflectanceimg(cs,file)
    return cs
     
    
def temperature ( file ):
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

def makeRGB (R, G, B,lat,lon,myx,myy, site,name='//////my_ratio_rgb_image.png'):
    plt.close("all")
    f=plt.figure(facecolor='Black')
    ax = f.add_axes([0.0,0.0,1.,1.],projection='rectilinear',position=[0.125,0.199,0.772,0.595])
    Rmax =np.median(R)+2*np.std(R)
    Rmin =np.median(R)-2*np.std(R)
    Rmin=0
    Gmax =np.median(G)+2*np.std(G)
    Gmin =np.median(G)-2*np.std(G)
    Gmin=0
    
    xlabel=float(myy)/float(np.shape(R)[1])
    ylabel=1.0-(float(myx)/float(np.shape(R)[0]))
    #print np.shape(R), myx, myy,xlabel,ylabel
    
    img = np.zeros((R.shape[0], R.shape[1], 3), dtype=float) 
    img[:,:,0] = linear(R, scale_min=Rmin, scale_max=Rmax)     
    img[:,:,1] = linear(G, scale_min=Gmin, scale_max=Gmax)  #np.max(G)/4)
    img[:,:,2] = linear(B, scale_min=250, scale_max=300)   #248-298  233-283
    
    #makelabels(img,lat,lon,ax)
    
    img[myx-3:myx+3,myy-4:myy+4,0],img[myx-3:myx+3,myy-4:myy+4,1],img[myx-3:myx+3,myy-4:myy+4,2]=0,0,0
    py.clf()
    rr=py.imshow(img, aspect=1.5)
    
    fn=name.split('/')[len(name.split('/'))-1]
    sat0    = fn.split('.')[0]
    year0   = fn.split('.')[1]
    julian0 = fn.split('.')[2]
    hhmm0   = fn.split('.')[3]
    
    fn=sat0+'-enhanced-'+year0+'.'+julian0+'.'+hhmm0+'.png'
    
    print 'MakeRGB fn:',fn
    #year=name.split('.')[1]
    #julian=name.split('.')[2]
    #hhmm=name.split('.')[3]
    #h=datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(julian) - 1)    
    #month=h.month
    #if month < 10:
    #    month='0'+str(month)
    #label1=' GOES-13      '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
    #plt.text(0.145, 0.013,label1,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
    
#    plt.text(0.5, 0.0,'0',transform=ax.transAxes, fontsize=6,color='white')
#    plt.text(0.5, 0.2,'2',transform=ax.transAxes, fontsize=6,color='white')
#    plt.text(0.5, 0.4,'4',transform=ax.transAxes, fontsize=6,color='white')
#    plt.text(0.5, 0.6,'6',transform=ax.transAxes, fontsize=6,color='white')
#    plt.text(0.5, 0.8,'8',transform=ax.transAxes, fontsize=6,color='white')
#    plt.text(0.5, 1.0,'1',transform=ax.transAxes, fontsize=6,color='white')
#    plt.text(0.0, 0.5,'0',transform=ax.transAxes, fontsize=6,color='black')
#    plt.text(0.2, 0.5,'2',transform=ax.transAxes, fontsize=6,color='black')
#    plt.text(0.4, 0.5,'4',transform=ax.transAxes, fontsize=6,color='black')
#    plt.text(0.6, 0.5,'6',transform=ax.transAxes, fontsize=6,color='black')
#    plt.text(0.8, 0.5,'8',transform=ax.transAxes, fontsize=6,color='black')
#    plt.text(1.0, 0.5,'1',transform=ax.transAxes, fontsize=6,color='black')
    
    plt.text(xlabel-0.02,ylabel+0.02,site,transform=ax.transAxes, fontsize=7,color='black')
    #plt.text(xlabel, ylabel,'0',transform=ax.transAxes, fontsize=8,color='white')
    clear_frame()
    py.savefig(fn,dpi=200.,bbox_inches='tight',pad_inches=0.0,transparent=True)
    plt.close("all")
    return rr

def makeratioRGB (R, G, B, name='//////my_ratio_rgb_image.png'):
    plt.close("all")
    f=plt.figure(facecolor='Black')
    ax = f.add_axes([0.1,0.1,0.9,0.9])
    Rmax =np.median(R)+3*np.std(R)
    #Gmax =np.median(G)+3*np.std(G)
    img = np.zeros((R.shape[0], R.shape[1], 3), dtype=float) 
    img[:,:,0] = linear(R, scale_min=0, scale_max=Rmax)     
    img[:,:,1] = linear(G, scale_min=248, scale_max=298)  #np.max(G)/4)
    img[:,:,2] = linear(B, scale_min=248, scale_max=298)   #248-298
    py.clf()
    rr=py.imshow(img, aspect=1.8)
    fn=name.split('/')[len(name.split('/'))-1]
    
    fn='GOES13-enhanced-'+fn[7:-10]+'ratio.png'
    print 'Ratio fn:',fn
    year=name.split('.')[1]
    julian=name.split('.')[2]
    hhmm=name.split('.')[3]
    h=datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(julian) - 1)    
    month=h.month
    if month < 10:
        month='0'+str(month)
    label1=' GOES-13      '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
    plt.text(0.145, 0.013,label1,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
    plt.text(0.459, 0.403,'*',transform=ax.transAxes, fontsize=18,color='white')
    plt.text(0.48, 0.405,'Manaus',transform=ax.transAxes, fontsize=8,color='white')
    clear_frame()
    py.savefig(fn,dpi=200.,bbox_inches='tight',pad_inches=0.0,transparent=True)
    plt.close("all")
    return rr


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


def phaseid (r063,r390,temp,name='//////my_ratio_rgb_image.png'):

    r063,r390,temp=np.squeeze(r063),np.squeeze(r390),np.squeeze(temp)
    ratio=r063/r390

    ice = (ratio>26.5) & (r063>0.15)     #ice= (temp<=233) & (r063>=0.4) 
    water= (ratio<13.2) & (r063>0.15)     #water = (temp >= 273) & (r063 >= 0.4)
    mix= (ratio >=13.2) & (ratio<=26.5) & (r063>0.15)

    R,G,B=r063,r063,r063
    R,G,B=R-r063,R-r063,R-r063
    R[water]=r063[water]
    G[mix]=r063[mix]
    B[ice]=r063[ice]

    ice= (temp<=233) & (r063>=0.4) 
    water = (temp >= 273) & (r063 >= 0.4)


    fh = open('waterhist063', 'w')
    str1 = ''.join(str(e)+'\n' for e in r063[water])
    fh.writelines(str1)
    fh.close()

    fh = open('waterhist390', 'w')
    str1 = ''.join(str(e)+'\n' for e in r390[water])
    fh.writelines(str1)
    fh.close()

    fh = open('icehist063', 'w')
    str1 = ''.join(str(e)+'\n' for e in r063[ice])
    fh.writelines(str1)
    fh.close()
 
    fh = open('icehist390', 'w')
    str1 = ''.join(str(e)+'\n' for e in r390[ice])
    fh.writelines(str1)
    fh.close()

    ice = (ratio>26.5) & (r063>0.2)     #ice= (temp<=233) & (r063>=0.4) 
    water= (ratio<13.2) & (r063>0.2)     #water = (temp >= 273) & (r063 >= 0.4)
    mix= (ratio >=13.2) & (ratio<=26.5) & (r063>0.2)

    cloud=r063>0.3
    fh = open('cloud063', 'w')
    str1 = ''.join(str(e)+'\n' for e in r063[cloud])
    fh.writelines(str1)
    fh.close()
    fh = open('cloud390', 'w')
    str1 = ''.join(str(e)+'\n' for e in r390[cloud])
    fh.writelines(str1)
    fh.close()
    fh = open('cloudtemp', 'w')
    str1 = ''.join(str(e)+'\n' for e in temp[cloud])
    fh.writelines(str1)
    fh.close()

    


    #makeratioRGB(ratio,temp,temp,name=listch1[count])

    #plt.close("all")
    #f=plt.figure(facecolor='Black')
    #ax = f.add_axes([0.1,0.1,0.9,0.9])
    #img = np.zeros((r063.shape[0], r063.shape[1], 3), dtype=float)
    #img[:,:,0],img[:,:,1],img[:,:,2]=r063,r063,r063
    #img[:,:,0] = linear(R)  #, scale_min=0, scale_max=Rmax)     
    #img[:,:,1] = linear(G)  #, scale_min=248, scale_max=298)  #np.max(G)/4)
    #img[:,:,2] = linear(B)  #, scale_min=248, scale_max=298)   #248-298
    #py.clf()
    #py.imshow(img, aspect=1.8)
    
    #fn=name.split('/')[len(name.split('/'))-1]
    #fn='GOES13-enhanced-'+fn[7:-10]+'phase.png'
    #print 'Phase fn:',fn
    #year=name.split('.')[1]
    #julian=name.split('.')[2]
    #hhmm=name.split('.')[3]
    #h=datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(julian) - 1)    
    #month=h.month
    #if month < 10:
    #    month='0'+str(month)
    #label1=' GOES-13      '+str(h.day)+'-'+str(month)+'-'+str(year)+'   '+str(hhmm)+' UTC          '
    #plt.text(0.147, 0.01,label1,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes, fontsize=8,color='white', bbox=dict(facecolor='black' ))
    #plt.text(0.459, 0.403,'*',transform=ax.transAxes, fontsize=18,color='white')
    #plt.text(0.48, 0.405,'Manaus',transform=ax.transAxes, fontsize=8,color='white')
    #clear_frame()
    #py.savefig(fn,dpi=200.,bbox_inches='tight',pad_inches=0.0,transparent=True)
    #plt.close("all")
    
    return

def getlatlon (file):
    fh = Dataset(file, mode='r')
    lon = fh.variables['lon'][:]
    lat = fh.variables['lat'][:]
    fh.close()
    return lat,lon

def makelabels (img,lat,lon,ax):
    
    R=img[:,:,0]
    

    site=['AF']
    mylat=[-9.87083333]
    mylon=[-56.10388889]
    
    site=['MN','AF','AH','JP']
    mylat=[-3.00, -9.87083333,-10.76000000,-10.93388889]
    mylon=[-60.0,-56.10388889,-62.35777778,-61.85194444]
        
    
    
    mylen,myang=mylat,mylat
    
   
    rlen=np.sqrt(lat*lat/(90*90) + lon*lon/(180*180))
    rang=np.arctan(lat/lon)
    print np.min(rlen), np.max(rlen), np.min(rang),np.max(rang)
    print np.shape(rlen),np.shape(rang)
    count=0
    print mylat,mylat[0]
    while count < len(site):
        print mylat
        mla= mylat[count]
        mlo=mylon[count]
        print mla
        if (mla>np.min(lat)) & (mla<np.max(lat)) & (mlo>np.min(lon)) & (mlo<np.max(lon)) :
            print mla
            mylen[count]=np.sqrt(mla*mla/(90*90) + mlo*mlo/(180*180))
            print 'jjj',mylen[count],np.shape(mylen[count])
            print mla
            print mylat[count],mylon[count],mylat[count]/mylon[count]
            myang[count]=np.arctan(mla/mlo)
            print 'myang', myang[count]
            difo=np.absolute((rlen-mylen[count])*(rang-myang[count]))
            loco=np.where(difo==np.min(difo))
            myx,myy=loco[0],loco[1]
            xlabel=float(myy)/float(np.shape(R)[1])
            ylabel=1.0-(float(myx)/float(np.shape(R)[0]))
            print xlabel, ylabel,myx,myy
            plt.text(xlabel-0.02,ylabel+0.02,site[count],transform=ax.transAxes, fontsize=7,color='black')
            img[myx-3:myx+3,myy-4:myy+4,0],img[myx-3:myx+3,myy-4:myy+4,1],img[myx-3:myx+3,myy-4:myy+4,2]=0,0,0
        count=count+1
    return
     
    
    
        
def getpostlaunchconstant(name):
    fn=name.split('/')[len(name.split('/'))-1]
    mysat    = fn.split('.')[0]
    myyear   = int(fn.split('.')[1])
    myjulian = int(fn.split('.')[2])
    refdate=dt.datetime(myyear-1,12,31)+dt.timedelta(myjulian)-dt.timedelta(days=30) #Subtract 30 days according to NOAA website
    initdate=dt.datetime(2003,4,1)
    mya=1.0848
    myb=0.0490
    if mysat == 'goes13':
        initdate=dt.datetime(2010,4,14)
        mya=1.1220 #For data acquired on Feb2013
        myb=0.0390 #For data acquired on Feb2013
    mydelta=refdate-initdate
    mydelta=mydelta.total_seconds()/(365.25*24*60*60)
    myc=1.0
    
    if ((myyear >= 2007) & (myjulian >= 213)) :
        myc=mya*np.exp(myb*mydelta)

    #print myc

    return myc    
    
    
    
    
    


def cloudprop (r063,r390,temp,mylat,mylon,site,name='my_ratio_rgb_image.png'):
    r063,r390,temp=np.squeeze(r063),np.squeeze(r390),np.squeeze(temp)
    lat, lon = getlatlon(name)
    rlen=np.sqrt(lat*lat/(90*90) + lon*lon/(180*180))
    rang=np.arctan(lat/lon)
    mylen=np.sqrt(mylat*mylat/(90*90) + mylon*mylon/(180*180))
    myang=np.arctan(mylat/mylon)
    #print mylen,myang
    dif=np.absolute((rlen-mylen)*(rang-myang))
    loc=np.where(dif==np.min(dif))
    myx,myy=loc[0],loc[1]
    #print myx,myy
    boxsize=[12] #[7,12,20]  Box size around the site, in pixels. Boxsize=7 means [x-7,y-7] to [x+7,y+7], or a square with 15 pixels of side.
    for b in boxsize:
        minx,miny,maxx,maxy=myx-b,myy-b,myx+b,myy+b
        r063b=r063[minx:maxx,miny:maxy]
        r390b=r390[minx:maxx,miny:maxy]
        tempb=temp[minx:maxx,miny:maxy]
        makeRGB(r063,r390,temp,lat,lon,myx,myy,site,name)
        cloud=(r063b>0.2)&(tempb<283.0)
        fh = open(name[:-11]+'_box'+str(b)+'_cloudreflectance063.csv', 'w')
        str1 = ''.join(str(e)+'\n' for e in r063b[cloud])
        fh.writelines(str1)
        fh.close()
        fh = open(name[:-11]+'_box'+str(b)+'_cloudreflectance390.csv', 'w')
        str1 = ''.join(str(e)+'\n' for e in r390b[cloud])
        fh.writelines(str1)
        fh.close()
        fh = open(name[:-11]+'_box'+str(b)+'_cloudtemperature.csv', 'w')
        str1 = ''.join(str(e)+'\n' for e in tempb[cloud])
        fh.writelines(str1)
        fh.close()
    return



#  MAIN LOOP

datadir=sys.argv[1]
curdir=os.getcwd()
outdir=curdir
if len(sys.argv)>2 :
    outdir=sys.argv[2]
os.chdir(outdir)

listch1=sorted(glob.glob(datadir+'/*BAND_01.nc'))
listch2=sorted(glob.glob(datadir+'/*BAND_02.nc'))
listch4=sorted(glob.glob(datadir+'/*BAND_04.nc'))

#site,mylat,mylon,tz='MN',-3.0,-60.0,-4 #Manaus
site,mylat,mylon,tz= 'AF',-9.87083333,-56.10388889, -4   #Alta Floresta
#site,mylat,mylon,tz='AH',-10.76000000,-62.35777778, -4   #Abracos Hill
#site,mylat,mylon,tz='JP',-10.93388889,-61.85194444, -4   #Ji-Parana

count=0
while count < len(listch1):
    mVIS,kVIS,C,m2,m3,m4,m6,bb2,bb3,bb4,bb6,factor=mVIS13,kVIS13,C13,m213,m313,m413,m613,bb213,bb313,bb413,bb613,factor13
    n2,n3,n4,n6,a2,a3,a4,a6,b2,b3,b4,b6,g2,g3,g4,g6=n213,n313,n413,n613,a213,a313,a413,a613,b213,b313,b413,b613,g213,g313,g413,g613
    SAT=listch1[count].split('/')[len(listch1[count].split('/'))-1].split('.')[0]
    if SAT=='goes12':
        mVIS,kVIS,C,m2,m3,m4,m6,bb2,bb3,bb4,bb6,factor=mVIS12,kVIS12,C12,m212,m312,m412,m612,bb212,bb312,bb412,bb612,factor12
        n2,n3,n4,n6,a2,a3,a4,a6,b2,b3,b4,b6,g2,g3,g4,g6=n212,n312,n412,n612,a212,a312,a412,a612,b212,b312,b412,b612,g212,g312,g412,g612
    if SAT=='goes10':
        mVIS,kVIS,C,m2,m3,m4,m6,bb2,bb3,bb4,bb6,factor=mVIS10,kVIS10,C10,m210,m310,m410,m610,bb210,bb310,bb410,bb610,factor10
        n2,n3,n4,n6,a2,a3,a4,a6,b2,b3,b4,b6,g2,g3,g4,g6=n210,n310,n410,n610,a210,a310,a410,a610,b210,b310,b410,b610,g210,g310,g410,g610
    C=getpostlaunchconstant(listch1[count])
    r063=reflectance(listch1[count])
    #r063=rebinned(r063,4)
    temp=temperature(listch4[count])
    r390=reflectance(listch2[count],temp)
    cloudprop(r063,r390,temp,mylat,mylon,site,name=listch2[count])
    count = count + 1
os.chdir(curdir)
exit()


