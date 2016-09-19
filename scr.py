from netCDF4 import Dataset
import numpy as np

#file = '/home/acorreia/Desktop/f.nc'
file='/home/acorreia/Downloads/goes13.2013.032.154518.BAND_04.nc'
#file='/home/acorreia/Downloads/hhh/goes13.2013.032.161518.BAND_01-1km.nc'



fh = Dataset(file, mode='r')

lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
datas = fh.variables['data'][:]

#data_units = fh.variables['data'].units
data_units = 'B4 4km'

fh.close()


import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


# Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()


f=plt.figure(num=15, figsize=(100, 40), dpi=180)




# When we give this Basemap instance our coordinate variables, it returns our plotting coordinates. This is how basemap knows where to put our projected data on the map.
#m = Basemap(width=5000000,height=3500000, resolution='l',projection='stere', lat_ts=0,lat_0=lat_0,lon_0=lon_0)
#m = Basemap(width=2500000,height=2000000, resolution='l',projection='stere', lon_0=lon_0,lat_0=lat_0)
m = Basemap(width=1500000,height=1000000, resolution='h',projection='stere', lon_0=-60,lat_0=-3)


# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
#lon, lat = np.meshgrid(lons, lats)

lon, lat = lons, lats
xi, yi = m(lon, lat)




#Colormap info
#The base colormaps are derived from those of the same name provided with Matlab:
#
#Colormap	Description
#autumn	sequential linearly-increasing shades of red-orange-yellow
#bone	sequential increasing black-white color map with a tinge of blue, to emulate X-ray film
#cool	linearly-decreasing shades of cyan-magenta
#copper	sequential increasing shades of black-copper
#flag	repetitive red-white-blue-black pattern (not cyclic at endpoints)
#gray	sequential linearly-increasing black-to-white grayscale
#hot	sequential black-red-yellow-white, to emulate blackbody radiation from an object at increasing temperatures
#hsv	cyclic red-yellow-green-cyan-blue-magenta-red, formed by changing the hue component in the HSV color space
#jet	a spectral map with dark endpoints, blue-cyan-yellow-red; based on a fluid-jet simulation by NCSA [1]
#pink	sequential increasing pastel black-pink-white, meant for sepia tone colorization of photographs
#prism	repetitive red-yellow-green-blue-purple-...-green pattern (not cyclic at endpoints)
#spring	linearly-increasing shades of magenta-yellow
#summer	sequential linearly-increasing shades of green-yellow
#winter	linearly-increasing shades of blue-green
#
#gist_earth mapmakers colors from dark blue deep ocean to green lowlands to brown highlands to white mountains
#gist_heat	sequential increasing black-red-orange-white, to emulate blackbody radiation from an iron bar as it grows hotter
#gist_ncar	pseudo-spectral black-blue-green-yellow-red-purple-white colormap from National Center for Atmospheric Research [2]
#gist_rainbow	runs through the colors in spectral order from red to violet at full saturation (like hsv but not cyclic)
#gist_stern	 Stern special color table from Interactive Data Language software
#The following colormaps are based on the ColorBrewer color specifications and designs developed by Cynthia Brewer:
#
#ColorBrewer Diverging (luminance is highest at the midpoint, and decreases towards differently-colored endpoints):
#
#Colormap	Description
#BrBG	brown, white, blue-green
#PiYG	pink, white, yellow-green
#PRGn	purple, white, green
#PuOr	orange, white, purple
#RdBu	red, white, blue
#RdGy	red, white, gray
#RdYlBu	red, yellow, blue
#RdYlGn	red, yellow, green
#Spectral	red, orange, yellow, green, blue
#ColorBrewer Sequential (luminance decreases monotonically):
#
#Colormap	Description
#Blues	white to dark blue
#BuGn	white, light blue, dark green
#BuPu	white, light blue, dark purple
#GnBu	white, light green, dark blue
#Greens	white to dark green
#Greys	white to black (not linear)
#Oranges	white, orange, dark brown
#OrRd	white, orange, dark red
#PuBu	white, light purple, dark blue
#PuBuGn	white, light purple, dark green
#PuRd	white, light purple, dark red
#Purples	white to dark purple
#RdPu	white, pink, dark purple
#Reds	white to dark red
#YlGn	light yellow, dark green
#YlGnBu	light yellow, light green, dark blue
#YlOrBr	light yellow, orange, dark brown
#YlOrRd	light yellow, orange, dark red
#ColorBrewer Qualitative:
#
#(For plotting nominal data, ListedColormap should be used, not LinearSegmentedColormap. Different sets of colors are recommended for different numbers of categories. These continuous versions of the qualitative schemes may be removed or converted in the future.)
#
#Accent
#Dark2
#Paired
#Pastel1
#Pastel2
#Set1
#Set2
#Set3
#Other miscellaneous schemes:
#
#Colormap	Description
#afmhot	sequential black-orange-yellow-white blackbody spectrum, commonly used in atomic force microscopy
#brg	blue-red-green
#bwr	diverging blue-white-red
#coolwarm	diverging blue-gray-red, meant to avoid issues with 3D shading, color blindness, and ordering of colors [3]
#CMRmap	Default colormaps on color images often reproduce to confusing grayscale images. The proposed colormap maintains an aesthetically pleasing color image that automatically reproduces to a monotonic grayscale with discrete, quantifiable saturation levels. [4]
#cubehelix	Unlike most other color schemes cubehelix was designed by D.A. Green to be monotonically increasing in terms of perceived brightness. Also, when printed on a black and white postscript printer, the scheme results in a greyscale with monotonically increasing brightness. This color scheme is named cubehelix because the r,g,b values produced can be visualised as a squashed helix around the diagonal in the r,g,b color cube.
#gnuplot	gnuplots traditional pm3d scheme (black-blue-red-yellow)
#gnuplot2	sequential color printable as gray (black-blue-violet-yellow-white)
#ocean	green-blue-white
#rainbow	spectral purple-blue-green-yellow-orange-red colormap with diverging luminance
#seismic	diverging blue-white-red
#nipy_spectral	black-purple-blue-green-yellow-red-white spectrum, originally from the Neuroimaging in Python project
#terrain	mapmakers colors, blue-green-yellow-brown-white, originally from IGOR Pro
#The following colormaps are redundant and may be removed in future versions. Its recommended to use the names in the descriptions instead, which produce identical output
#
#Colormap	Description
#gist_gray	identical to gray
#gist_yarg	identical to gray_r
#binary	identical to gray_r
#spectral	identical to nipy_spectral [5]










#Now, we can plot the data using one of the available plot types (pcolor, pcolormesh, contour, contourf, scatter, etc.). Here we use pcolor. Gridlines, colorbars, and axis labels can also be added at this point.

# Plot Data
#cs = m.pcolormesh(xi,yi,np.squeeze(datas),cmap='gray_r',vmin=5000,vmax=15000)
cs = m.pcolormesh(xi,yi,np.squeeze(datas),cmap='gray_r')


# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(data_units)




#plt.figure(figsize=(8, 6), dpi=80)

#f=plt.figure(num=15, figsize=(16, 12), dpi=200)
# Add Title
plt.title('01FEB2013 GOES-13 IR Channel 4')

#f=plt.gcf()
DPI = f.get_dpi()
print "DPI:", DPI
DefaultSize = f.get_size_inches()
print "Default size in Inches", DefaultSize
#print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

#f.set_size_inches( 16.,12. )
#f.set_dpi(200)
#f.set_size_inches(16.,12.)

#DPI = f.get_dpi()
#print "DPI:", DPI
#NewSize = f.get_size_inches()
#print "New size in Inches", NewSize

plt.show()


