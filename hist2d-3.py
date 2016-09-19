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
import time
import pandas as pd
#import matplotlib as mpl


myfile=sys.argv[1]

x=np.array([0.])

y,z=x,x
start=time.clock()

data=pd.read_csv(myfile,header=None,names=['x','y','z'])
x=data.x   #VIS 0.63um
y=data.y   #IR  3.90um
z=data.z   #Temp K

end=time.clock()
print "elapsed time: ",end-start
print "Total number of points read: ",len(x) 
z=z-np.float64(273.15)
ny=1-y
reff=ny


a,b,c,d,e=-14.5168,6.68151,25.558,-23.8399,27.6392
reff=a+b*np.exp(c*ny+d)+e*ny


#plt.jet()

#xlims = [0.195,0.605]    #VIS  limits
#ylims = [0.0,0.205]      #IR   limits
#zlims = [-60.0,10.0]     #Temp limits
#nylims = [0.795,1.0]  # (1-refl 3.90) limits
nzlims = [10.0, -80.0] # inverted temperature limits
rlims = [0.0,45.0]  #reff um limits


# Find the min/max of the data
#xmin = min(xlims)
#xmax = max(xlims)
#ymin = min(ylims)
#ymax = max(ylims)
#zmin = min(zlims)
#zmax = max(zlims) 

nzmax = max(nzlims)
nzmin = min(nzlims)
#nymax = max(nylims)
#nymin = min(nylims)
rmax  = max (rlims)
rmin  = min (rlims)

 
# Set up x and y labels
#xlabel = '$\mathrm{Reflectance\\ at\\ 0.63\\ \mu m}$'
#ylabel = '$\mathrm{Reflectance\\ at\\ 3.90\\ \mu m}$'
nzlabel = '$\mathrm{Temperature\\ (\degree C)}$'
hlabel = '$\mathrm{Relative\\ frequency\\ ( \% )}$'

#nzlabel = zlabel
#nylabel = '$\mathrm{1\\ -\\ Reflectance\\ at\\ 3.90\\ \mu m}$'
rlabel = '$\mathrm{Effective\\ radius\\ (\mu m)}$'


 
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
#nxbins = 181
nybins = 50
nzbins = 50
#nbins =  90
 
#xbins = linspace(start = xmin, stop = xmax, num = nxbins)
#ybins = linspace(start = ymin, stop = ymax, num = nybins)
zbins = linspace(start = nzmin, stop = nzmax, num = nzbins)

#nnzbins = linspace(start = nzmax, stop = nzmin, num = nzbins)
#nnybins = linspace(start = nymin, stop = nymax, num = nybins)
rbins = linspace(start = rmin, stop = rmax, num = nybins)


aspectratioYZ  = (1.0*rmax-1.0*rmin)/(1.0*nzmax-1.0*nzmin)


nHYZ, nxyzedges,nyyzedges = np.histogram2d(z,reff,bins=(zbins,rbins),normed=True)
nHYZ=100.0*nHYZ/(1.0*np.sum(np.abs(nHYZ)))

print np.sum(nHYZ)

#Setting the limit for the colorbar
plt.jet()
colorbarmin=0.0
colorbarmax=1.
colorbarticknum=10
nHYZ = nHYZ.clip(min=colorbarmin,max=colorbarmax)
#mymin = (nHYZ == np.min(nHYZ))
#mymax = (nHYZ == np.max(nHYZ))
nHYZ[0,1] = colorbarmin
nHYZ[0,0] = colorbarmax


# Set up the size of the figure
fig = plt.figure(3, figsize=(10.0,8.0))
 
# Make the three plots
axTemperature = plt.axes(rect_temperature) # temperature plot
 
# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()

# Plot the temperature data
cax = (axTemperature.imshow(nHYZ, extent=[rmin,rmax,nzmax,nzmin],
       interpolation='None', origin='upper',aspect=aspectratioYZ))
#print np.shape(cax)
cbaxes = fig.add_axes([.9, 0.1, 0.03, 0.85])
cbar=plt.colorbar(cax,cax=cbaxes)
#                  norm=mpl.colors.Normalize(vmin=colorbarmin, vmax=colorbarmax))
#cbar=plt.colorbar(CS3)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_title(hlabel,rotation='vertical',position=(-0.8,0.7),fontsize=14)

cbar.set_clim(colorbarmin,colorbarmax)
cbar.set_ticks(np.linspace(colorbarmin,colorbarmax,num=colorbarticknum+1))
#[0.00,0.05,0.10,0.15,0.20,0.25])

#Plot the axes labels
axTemperature.set_xlabel(rlabel,fontsize=16)
axTemperature.set_ylabel(nzlabel,fontsize=16)
 
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
axTemperature.set_xlim(rlims)
axTemperature.set_ylim(nzlims)
 
#Show the plot
plt.draw()
 
# Save to a File
filename = myfile[:-3]+'TempxReff'
plt.savefig(filename + '.png',format = 'png', transparent=False)





sys.exit()
