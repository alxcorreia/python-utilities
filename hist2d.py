# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:57:15 2014

@author: acorreia
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator, LinearLocator, FormatStrFormatter
from numpy import linspace
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import time


# Define a function to make the ellipses
def ellipse(ra,rb,ang,x0,y0,Nb=100):
    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang
    co,si=np.cos(an),np.sin(an)
    the=linspace(0,2*np.pi,Nb)
    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
    return X,Y

myfile=sys.argv[1]

x=np.array([0.])
y,z=x,x
start=time.clock()
# Define the x and y data 
with open (myfile, 'r') as mf:
    print "Reading: ",myfile
    for line in mf:
        if line.startswith('site'):
            continue
        RefVIS,RefIR,Temp=float(line.split(',')[0]),float(line.split(',')[1]),float(line.split(',')[2])
        if len(x) == 1:
            x,y,z=np.array([RefVIS]),np.array([RefIR]),np.array([Temp])
        x,y,z = np.append(x,RefVIS),np.append(y,RefIR),np.append(z,Temp)

end=time.clock()
print "elapsed time: ",end-start
print "Total number of points read: ",len(x) 
z=z-273.15
print "z shape, min, max median: ",np.shape(z),np.min(z),np.max(z),np.median(z)

plt.jet()
# Set up default x and y limits
#xlims = [min(x)-0.05,max(x)+0.05]
#ylims = [min(y)-0.05,max(y)+0.05]
xlims = [0.195,0.605]    #VIS  limits
ylims = [0.0,0.205]      #IR   limits
zlims = [-80.0,10.0]     #Temp limits

# Find the min/max of the data
xmin = min(xlims)
xmax = max(xlims)
ymin = min(ylims)
ymax = max(ylims)
zmin = min(zlims)
zmax = max(zlims) 
 
# Set up x and y labels
xlabel = '$\mathrm{Reflectance\\ 0.63um}$'
ylabel = '$\mathrm{Reflectance\\ 3.90um}$'
zlabel = '$\mathrm{Temperature\\ (deg\\ C)}$'
 
# Define the locations for the axes
left, width = 0.12, 0.55
bottom, height = 0.12, 0.55
bottom_h = left_h = left+width+0.02
 
# Set up the geometry of the three plots
rect_temperature = [left, bottom, width, height] # dimensions of temp plot
rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
 
# Make the 'main' temperature plot
# Define the number of bins
nxbins = 97
nybins = 97
nzbins = 97
nbins = 50
 
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
print "shapes H,xedges,yedges: ", np.shape(H),np.shape(xedges),np.shape(yedges)


HXZ, xxzedges,yxzedges = np.histogram2d(z,x,bins=(zbins,xbins),normed=True)
print "shapes HXZ,xxzedges,yxzedges: ", np.shape(HXZ),np.shape(xxzedges),np.shape(yxzedges)


HYZ, xyzedges,yyzedges = np.histogram2d(z,y,bins=(zbins,ybins),normed=True)
print "shapes HYZ,xyzedges,yyzedges: ", np.shape(HYZ),np.shape(xyzedges),np.shape(yyzedges)





# Set up the size of the figure
fig = plt.figure(1, figsize=(9.5,9.5))
 
# Make the three plots
axTemperature = plt.axes(rect_temperature) # temperature plot
axHistx = plt.axes(rect_histx) # x histogram
axHisty = plt.axes(rect_histy) # y histogram
 
# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# Plot the temperature data
cax = (axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax],
       interpolation='nearest', origin='lower',aspect=aspectratio))

#Plot the axes labels
axTemperature.set_xlabel(xlabel,fontsize=18)
axTemperature.set_ylabel(ylabel,fontsize=18)
 
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
xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)
 
#Plot the histograms
axHistx.hist(x, bins=xbins, color = 'blue')
axHisty.hist(y, bins=ybins, orientation='horizontal', color = 'red')
 
#Set up the histogram limits
axHistx.set_xlim( min(x), max(x) )
axHisty.set_ylim( min(y), max(y) )
 
#Make the tickmarks pretty
ticklabels = axHistx.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(8)
    label.set_family('serif')
 
#Make the tickmarks pretty
ticklabels = axHisty.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(8)
    label.set_family('serif')
 
#Cool trick that changes the number of tickmarks for the histogram axes
axHisty.xaxis.set_major_locator(MaxNLocator(4))
axHistx.yaxis.set_major_locator(MaxNLocator(4))
 
#Show the plot
plt.draw()
 
# Save to a File
filename = myfile[:-3]+'IRxVIS'
plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)







# Set up the size of the figure
fig = plt.figure(2, figsize=(9.5,9.5))
 
# Make the three plots
axTemperature = plt.axes(rect_temperature) # temperature plot
axHistx = plt.axes(rect_histx) # x histogram
axHisty = plt.axes(rect_histy) # y histogram
 
# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# Plot the temperature data
cax = (axTemperature.imshow(HXZ, extent=[xmin,xmax,zmin,zmax],
       interpolation='nearest', origin='lower',aspect=aspectratioXZ))

#Plot the axes labels
axTemperature.set_xlabel(xlabel,fontsize=18)
axTemperature.set_ylabel(zlabel,fontsize=18)
 
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
axHistx.hist(x, bins=xbins, color = 'blue')
axHisty.hist(z, bins=zbins, orientation='horizontal', color = 'red')
 
#Set up the histogram limits
axHistx.set_xlim( min(x), max(x) )
axHisty.set_ylim( min(z), max(z) )
 
#Make the tickmarks pretty
ticklabels = axHistx.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(8)
    label.set_family('serif')
 
#Make the tickmarks pretty
ticklabels = axHisty.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(8)
    label.set_family('serif')
 
#Cool trick that changes the number of tickmarks for the histogram axes
axHisty.xaxis.set_major_locator(MaxNLocator(4))
axHistx.yaxis.set_major_locator(MaxNLocator(4))
 
#Show the plot
plt.draw()
 
# Save to a File
filename = myfile[:-3]+'TempxVIS'
plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)















# Set up the size of the figure
fig = plt.figure(3, figsize=(9.5,9.5))
 
# Make the three plots
axTemperature = plt.axes(rect_temperature) # temperature plot
axHistx = plt.axes(rect_histx) # x histogram
axHisty = plt.axes(rect_histy) # y histogram
 
# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# Plot the temperature data
cax = (axTemperature.imshow(HYZ, extent=[ymin,ymax,zmin,zmax],
       interpolation='nearest', origin='lower',aspect=aspectratioYZ))

#Plot the axes labels
axTemperature.set_xlabel(ylabel,fontsize=18)
axTemperature.set_ylabel(zlabel,fontsize=18)
 
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
axHistx.hist(y, bins=ybins, color = 'blue')
axHisty.hist(z, bins=zbins, orientation='horizontal', color = 'red')
 
#Set up the histogram limits
axHistx.set_xlim( min(y), max(y) )
axHisty.set_ylim( min(z), max(z) )
 
#Make the tickmarks pretty
ticklabels = axHistx.get_yticklabels()
for label in ticklabels:
    label.set_fontsize(8)
    label.set_family('serif')
 
#Make the tickmarks pretty
ticklabels = axHisty.get_xticklabels()
for label in ticklabels:
    label.set_fontsize(8)
    label.set_family('serif')
 
#Cool trick that changes the number of tickmarks for the histogram axes
axHisty.xaxis.set_major_locator(MaxNLocator(4))
axHistx.yaxis.set_major_locator(MaxNLocator(4))
 
#Show the plot
plt.draw()
 
# Save to a File
filename = myfile[:-3]+'TempxIR'
plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)











sys.exit()










#X = xcenter
#Y = ycenter
#Z = H
 
# Plot the temperature plot contours
#contourcolor = 'white'
#xcenter = np.mean(x)
#ycenter = np.mean(y)
#ra = np.std(x)
#rb = np.std(y)
#ang = 0
# 
#X,Y=ellipse(ra,rb,ang,xcenter,ycenter)
#axTemperature.plot(X,Y,"k:",ms=1,linewidth=2.0)
#axTemperature.annotate('$1\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
#                       textcoords='offset points', horizontalalignment='right',
#                       verticalalignment='bottom',fontsize=25)
# 
#X,Y=ellipse(2*ra,2*rb,ang,xcenter,ycenter)
#axTemperature.plot(X,Y,"k:",color = contourcolor,ms=1,linewidth=2.0)
#axTemperature.annotate('$2\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
#                       textcoords='offset points',horizontalalignment='right',
#                       verticalalignment='bottom',fontsize=25, color = contourcolor)
# 
#X,Y=ellipse(3*ra,3*rb,ang,xcenter,ycenter)
#axTemperature.plot(X,Y,"k:",color = contourcolor, ms=1,linewidth=2.0)
#axTemperature.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
#                       textcoords='offset points',horizontalalignment='right',
#                       verticalalignment='bottom',fontsize=25, color = contourcolor)
 




fig=plt.figure(4, figsize=(9.5,9.5))
ax=fig.gca(projection='3d')
x,y=np.meshgrid(xedges,yedges)
surf=ax.plot_surface(x,y,z,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
ax.set_zlim(np.median(z)-3.0*np.std(z),np.median(z)+3.0*np.std(z))
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf,shrink=0.5, aspect=5)
plt.show()
filename = myfile+'myplot3D'
plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)






