# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:57:15 2014

@author: acorreia
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter#, MaxNLocator, LinearLocator, FormatStrFormatter
from numpy import linspace
#from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D
#import time



import pandas as pd


myfile=sys.argv[1]
x=np.array([0.])
y,z=x,x

data=pd.read_csv(myfile,header=None,names=['x','y','z'])
RefVIS=data.x
RefIR=data.y
Temp=data.z
Temp=Temp-np.float64(273.15)
N=len(RefVIS)
print "Total number of points read: ",N 

fn=myfile.split('/')[len(myfile.split('/'))-1]
fn=fn[:-3]

ice   = (RefVIS > 0.4) & (Temp < -40) & (RefIR > 0)
water = (RefVIS > 0.4) & (Temp >  0 ) & (RefIR > 0)
iceratio=(RefVIS[ice]/RefIR[ice])
waterratio=(RefVIS[water]/RefIR[water])
print "Water: N,min,max,median,stdev:", np.shape(waterratio)[0], np.min(waterratio), np.max(waterratio), np.median(waterratio),np.std(waterratio)
print "Ice  : N,min,max,median,stdev:", np.shape(iceratio)[0], np.min(iceratio), np.max(iceratio), np.median(iceratio),np.std(iceratio)

icedata   = pd.DataFrame({'RefVIS': RefVIS[ice]  ,'RefIR':RefIR[ice]  ,'Temp':Temp[ice]})
waterdata = pd.DataFrame({'RefVIS': RefVIS[water],'RefIR':RefIR[water],'Temp':Temp[water]})

icedata.to_csv(fn+'icedata.csv')
waterdata.to_csv(fn+'waterdata.csv')



#ICE HISTOGRAM
hist,bin_edges=np.histogram(iceratio,bins=100)
hist=100.0*hist/(1.0*np.sum(hist))
plt.plot(bin_edges[0:-1],hist,drawstyle='steps',fillstyle='full',linewidth=2.5)
plt.xlabel('VIS/IR ratio')
plt.ylabel('$\mathrm{Relative\\ frequency\\ ( \% )}$')
plt.title(r'Ice particles')
plt.axvline(linewidth=1.5, color='r',x=np.median(iceratio))
plt.axvline(linewidth=1, color='r',x=(np.median(iceratio)+np.std(iceratio)),linestyle='--')
plt.axvline(linewidth=1, color='r',x=(np.median(iceratio)-np.std(iceratio)),linestyle='--')
plt.show()

#WATER HISTOGRAM
histo,bin_edgeso=np.histogram(waterratio,bins=100,range=(0.,15.))
histo=100.0*histo/(1.0*np.sum(histo))
plt.plot(bin_edgeso[0:-1],histo,drawstyle='steps',fillstyle='full',linewidth=2.5)
plt.xlabel('VIS/IR ratio')
plt.ylabel('$\mathrm{Relative\\ frequency\\ ( \% )}$')
plt.title(r'Water particles')
plt.axvline(linewidth=1.5, color='r',x=np.median(waterratio))
plt.axvline(linewidth=1, color='r',x=(np.median(waterratio)+np.std(waterratio)),linestyle='--')
plt.axvline(linewidth=1, color='r',x=(np.median(waterratio)-np.std(waterratio)),linestyle='--')
plt.show()


exit()








plt.jet()
#cmap = plt.colors.ListedColormap(['b','g','y','r'])
# Set up default x and y limits
#xlims = [min(x)-0.05,max(x)+0.05]
#ylims = [min(y)-0.05,max(y)+0.05]
xlims = [0.195,0.605]    #VIS  limits
ylims = [0.0,0.205]      #IR   limits
zlims = [-60.0,10.0]     #Temp limits

# Find the min/max of the data
xmin = min(xlims)
xmax = max(xlims)
ymin = min(ylims)
ymax = max(ylims)
zmin = min(zlims)
zmax = max(zlims) 
 
# Set up x and y labels
xlabel = '$\mathrm{Reflectance\\ at\\ 0.63\\ \mu m}$'
ylabel = '$\mathrm{Reflectance\\ at\\ 3.90\\ \mu m}$'
zlabel = '$\mathrm{Temperature\\ (\degree C)}$'
hlabel = '$\mathrm{Relative\\ frequency\\ ( \% )}$'
 
# Define the locations for the axes
left, width = 0.1, 0.85
bottom, height = 0.1, 0.85
bottom_h = left_h = left+width+0.02
 
# Set up the geometry of the three plots
rect_temperature = [left, bottom, width, height] # dimensions of temp plot
rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
 
# Make the 'main' temperature plot
# Define the number of bins
nxbins = 181
nybins = 90
nzbins = 90
nbins =  90
 
xbins = linspace(start = xmin, stop = xmax, num = nxbins)
ybins = linspace(start = ymin, stop = ymax, num = nybins)
zbins = linspace(start = zmin, stop = zmax, num = nzbins)

#xcenter = (xbins[0:-1]+xbins[1:])/2.0
#ycenter = (ybins[0:-1]+ybins[1:])/2.0
#zcenter = (zbins[0:-1]+zbins[1:])/2.0

#aspectratio   = 1.0*(xmax - 0)/(1.0*ymax - 0)
#aspectratioXZ = 1.0*(xmax - 0)/(1.0*np.abs(zmax) - 0)
#aspectratioYZ = 1.0*(ymax - 0)/(1.0*np.abs(zmax) - 0)

aspectratio    = (1.0*xmax-1.0*xmin)/(1.0*ymax-1.0*ymin)
aspectratioXZ  = (1.0*xmax-1.0*xmin)/(1.0*zmax-1.0*zmin)
aspectratioYZ  = (1.0*ymax-1.0*ymin)/(1.0*zmax-1.0*zmin)

#print "aspectratio: ",aspectratio,aspectratioXZ,aspectratioYZ
 
H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins),normed=True)
#print "shapes H,xedges,yedges: ", np.shape(H),np.shape(xedges),np.shape(yedges), np.min(H),np.max(H),np.sum(H)

IRedges,VISedges=xedges, yedges
mt=np.zeros((np.shape(xedges)[0],np.shape(yedges)[0]))+100.0

IRe1=IRedges[0:-1]
IRe2=IRedges[1:]
VIe1=VISedges[0:-1]
VIe2=VISedges[1:]



H=100.0*H/(1.0*np.sum(H))
H=H.clip(min=0,max=0.028)

HXZ, xxzedges,yxzedges = np.histogram2d(z,x,bins=(zbins,xbins),normed=True)
#print "shapes HXZ,xxzedges,yxzedges: ", np.shape(HXZ),np.shape(xxzedges),np.shape(yxzedges)

HXZ=100.0*HXZ/(1.0*np.sum(HXZ))


HYZ, xyzedges,yyzedges = np.histogram2d(z,y,bins=(zbins,ybins),normed=True)
#print "shapes HYZ,xyzedges,yyzedges: ", np.shape(HYZ),np.shape(xyzedges),np.shape(yyzedges)

HYZ=100.0*HYZ/(1.0*np.sum(HYZ))




# Set up the size of the figure
fig = plt.figure(1, figsize=(10.,8.))
 
# Make the three plots
axTemperature = plt.axes(rect_temperature) # temperature plot
 
# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()

# Plot the temperature data
cax = (axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax],
       interpolation='None', origin='lower',aspect=aspectratio))

cbaxes = fig.add_axes([.9, 0.1, 0.03, 0.85])
cbar=plt.colorbar(cax,cax=cbaxes)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_title(hlabel,rotation='vertical',position=(-0.8,0.7),fontsize=14)

cbar.set_clim(0.000, 0.028)
cbar.set_ticks([0.000,0.002,0.004,0.006,0.008,0.010, 0.012,0.014,0.016,0.018,0.020, 0.022,0.024, 0.026])#,0.035,0.040,0.045,0.050])

#Plot the axes labels
axTemperature.set_xlabel(xlabel,fontsize=16)
axTemperature.set_ylabel(ylabel,fontsize=16)

axTemperature.xaxis.set_ticks([0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60])
axTemperature.yaxis.set_ticks([0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20])

 
#Make the tickmarks pretty
ticklabels = axTemperature.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(14)
    label.set_family('serif')
 
ticklabels = axTemperature.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(14)
    label.set_family('serif')
 
#Set up the plot limits
axTemperature.set_xlim(xlims)
axTemperature.set_ylim(ylims)
 
#Set up the histogram bins
#xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
#ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)
 
#Show the plot
plt.draw()
 
# Save to a File
filename = myfile[:-3]+'IRxVIS'
plt.savefig(filename + '.png',format = 'png', transparent=False)






plt.jet()
# Set up the size of the figure
fig = plt.figure(2, figsize=(10.,8.0))
 
# Make the three plots
axTemperature = plt.axes(rect_temperature) # temperature plot
#axHistx = plt.axes(rect_histx) # x histogram
#axHisty = plt.axes(rect_histy) # y histogram
 
# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()
#axHistx.xaxis.set_major_formatter(nullfmt)
#axHisty.yaxis.set_major_formatter(nullfmt)

# Plot the temperature data
cax = (axTemperature.imshow(HXZ, extent=[xmin,xmax,zmin,zmax],
       interpolation='None', origin='lower',aspect=aspectratioXZ))


cbaxes = fig.add_axes([.9, 0.1, 0.03, 0.85])
cbar=plt.colorbar(cax,cax=cbaxes)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_title(hlabel,rotation='vertical',position=(-0.8,0.7),fontsize=14)



#Plot the axes labels
axTemperature.set_xlabel(xlabel,fontsize=16)
axTemperature.set_ylabel(zlabel,fontsize=16)

#axTemperature.xaxis.set_ticks([0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60])
#axTemperature.yaxis.set_ticks([0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20])

 
#Make the tickmarks pretty
ticklabels = axTemperature.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(14)
    label.set_family('serif')
 
ticklabels = axTemperature.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(14)
    label.set_family('serif')
 
#Set up the plot limits
axTemperature.set_xlim(xlims)
axTemperature.set_ylim(zlims)
 
#Set up the histogram bins
xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
zbins = np.arange(zmin, zmax, (zmax-zmin)/nbins)
 
#Plot the histograms
#axHistx.hist(x, bins=xbins, color = 'blue')
#axHisty.hist(z, bins=zbins, orientation='horizontal', color = 'red')
 
#Set up the histogram limits
#axHistx.set_xlim( min(x), max(x) )
#axHisty.set_ylim( min(z), max(z) )
 
#Make the tickmarks pretty
#ticklabels = axHistx.get_yticklabels()
#for label in ticklabels:
#    label.set_fontsize(8)
#    label.set_family('serif')
 
#Make the tickmarks pretty
#ticklabels = axHisty.get_xticklabels()
#for label in ticklabels:
#    label.set_fontsize(8)
#    label.set_family('serif')
 
#Cool trick that changes the number of tickmarks for the histogram axes
#axHisty.xaxis.set_major_locator(MaxNLocator(4))
#axHistx.yaxis.set_major_locator(MaxNLocator(4))
 
#Show the plot
plt.draw()
 
# Save to a File
filename = myfile[:-3]+'TempxVIS'
plt.savefig(filename + '.png',format = 'png', transparent=False)














plt.jet()
# Set up the size of the figure
fig = plt.figure(3, figsize=(10.0,8.0))
 
# Make the three plots
axTemperature = plt.axes(rect_temperature) # temperature plot
#axHistx = plt.axes(rect_histx) # x histogram
#axHisty = plt.axes(rect_histy) # y histogram
 
# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()
#axHistx.xaxis.set_major_formatter(nullfmt)
#axHisty.yaxis.set_major_formatter(nullfmt)

# Plot the temperature data
cax = (axTemperature.imshow(HYZ, extent=[ymin,ymax,zmin,zmax],
       interpolation='None', origin='lower',aspect=aspectratioYZ))


cbaxes = fig.add_axes([.9, 0.1, 0.03, 0.85])
cbar=plt.colorbar(cax,cax=cbaxes)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_title(hlabel,rotation='vertical',position=(-0.8,0.7),fontsize=14)



#Plot the axes labels
axTemperature.set_xlabel(ylabel,fontsize=16)
axTemperature.set_ylabel(zlabel,fontsize=16)

#axTemperature.xaxis.set_ticks([0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60])
#axTemperature.yaxis.set_ticks([0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20])

 
#Make the tickmarks pretty
ticklabels = axTemperature.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(14)
    label.set_family('serif')
 
ticklabels = axTemperature.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(14)
    label.set_family('serif')
 
#Set up the plot limits
axTemperature.set_xlim(ylims)
axTemperature.set_ylim(zlims)
 
#Set up the histogram bins
ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)
zbins = np.arange(zmin, zmax, (zmax-zmin)/nbins)
 
#Plot the histograms
#axHistx.hist(y, bins=ybins, color = 'blue')
#axHisty.hist(z, bins=zbins, orientation='horizontal', color = 'red')
 
#Set up the histogram limits
#axHistx.set_xlim( min(y), max(y) )
#axHisty.set_ylim( min(z), max(z) )
 
#Make the tickmarks pretty
#ticklabels = axHistx.get_yticklabels()
#for label in ticklabels:
#    label.set_fontsize(8)
#    label.set_family('serif')
 
#Make the tickmarks pretty
#ticklabels = axHisty.get_xticklabels()
#for label in ticklabels:
#    label.set_fontsize(8)
#    label.set_family('serif')
 
#Cool trick that changes the number of tickmarks for the histogram axes
#axHisty.xaxis.set_major_locator(MaxNLocator(4))
#axHistx.yaxis.set_major_locator(MaxNLocator(4))
 
#Show the plot
plt.draw()
 
# Save to a File
filename = myfile[:-3]+'TempxIR'
plt.savefig(filename + '.png',format = 'png', transparent=False)










i=1
while i < len (VISedges):
    j=1
    while j < len (IRedges) :
        loc = ((x>=VISedges[i-1])&(x<VISedges[i])&(y>=IRedges[j-1])&(y<IRedges[j]))
        if len(z[loc])>0 :
            mt[j-1,i-1]=np.median(z[loc])
        j=j+1
    i=i+1




# Set up the size of the figure
fig = plt.figure(4, figsize=(10.,8.))

 
# Make the three plots
axTemperature = plt.axes(rect_temperature,aspect=aspectratio) # temperature plot
#axHistx = plt.axes(rect_histx) # x histogram
#axHisty = plt.axes(rect_histy) # y histogram
 
# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()
#axHistx.xaxis.set_major_formatter(nullfmt)
#axHisty.yaxis.set_major_formatter(nullfmt)





#Plot the axes labels
axTemperature.set_xlabel(xlabel,fontsize=16)
axTemperature.set_ylabel(ylabel,fontsize=16)
 
axTemperature.xaxis.set_ticks([0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60])
axTemperature.yaxis.set_ticks([0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20])
 
#Make the tickmarks pretty
ticklabels = axTemperature.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(14)
    label.set_family('serif')
 
ticklabels = axTemperature.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(14)
    label.set_family('serif')
 
#Set up the plot limits
axTemperature.set_xlim(xlims)
axTemperature.set_ylim(ylims)




X,Y = np.meshgrid(VISedges,IRedges)
mm=plt.pcolormesh(X,Y,mt,vmin=-40.,vmax=10.)



cbaxes = fig.add_axes([.91, 0.1, 0.03, 0.85])
cbar=plt.colorbar(mm,cax=cbaxes)

cbar.ax.tick_params(labelsize=10)
cbar.ax.set_title(zlabel,rotation='vertical',position=(-0.8,0.6),fontsize=14)

cbar.set_clim(-40.0, 20.0)
plt.draw()
 
# Save to a File
filename = myfile[:-3]+'TempxIRxVIS'
plt.savefig(filename + '.png',format = 'png', transparent=False)

sys.exit()



