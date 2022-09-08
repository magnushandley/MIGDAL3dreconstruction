import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from PIL import Image
from PIL.TiffTags import TAGS
from scipy import interpolate
from ITOstripcorews import *
from RFsingletiffws import *
from collections import Counter
import time

pathname = "/Users/magnus/Documents/MIGDAL/deconvolution/betatoreconstruct.tif"
pathnameraw = "/Users/magnus/Documents/MIGDAL/deconvolution/betatoreconstructnodark.tif"

SIGMA = 2.5 #sigma for derivative determination ~> Related to track width
lthresh = 10 #tracks with a response lower than this are rejected (0 accepts all)
uthresh = 0 #tracks with a response higher than this are rejected (0 accepts all)
minlen = 6 #minimum track length accepted
linkthresh = 1 #maximum distance to be linked
logim = False


braggimg = io.imread(pathnameraw)/20

#strip offset between image x coordinate and strips, in strips
offset = 3
#velocity in cm/mus
driftvel = 13

daqfilepath = "/Users/magnus/Documents/MIGDAL/datawithito/3D/MIG_Fe55_data_220803T152221.CAL/daq.txt"

event = 3054
tic = time.time()
x,y = returnlines(pathname,
    SIGMA,
    lthresh,
    uthresh,
    minlen,
    linkthresh,
    logim)
toc = time.time()
print(toc-tic)

#taking every ::nth value
x = x[::5]
y = y[::5]

print(len(x),len(y))
#deposition = np.log(RF.simplebragg(y,x,braggimg))



trackstrip = pixeltostrip(x,offset)


strips = ITO(daqfilepath,1000,event+1)

#strips.readevent(event)
#print(strips.data[event])
#np.save("CF4_2221_3055",strips.data[event])

eventdata = np.load("/Users/magnus/Documents/workshop_code/CF4_2221_3055.npy")
strips.data[event] = eventdata

peaks = np.zeros(120)
halfpeaks = strips.peaktime(event)
peaks[:60] = halfpeaks
peaks[60:] = halfpeaks

stripdep = np.zeros(120)
halfstripdep = strips.stripintegrals(event)
stripdep[:60] = halfstripdep
stripdep[60:] = halfstripdep
stripdepinterp = interpolate.interp1d(np.arange(120),stripdep)
deposition = stripdepinterp(trackstrip)
        
interpolatedpeaks = interpolate.interp1d(np.arange(120), peaks)
#

deltaz = np.zeros(len(x))

for i in range(len(trackstrip)):
    deltaz[i] = interpolatedpeaks(trackstrip[i])*driftvel

"""
BREAK FIRST DEGENERACY
"""
#identify any duplicate x coordinates, corresponding to the track turning back on itself
#print(x)
#d = Counter(trackstrip)
#duplicates = list([item for item in d if d[item]>1])
#print(duplicates)
#
#for i in range(len(trackstrip)):
#
#    if (trackstrip[i] not in duplicates):
#        deltaz[i] = peaks[int(trackstrip[i])]*driftvel
#
#    elif (trackstrip[i] in duplicates):
#        if (d[trackstrip[i]] == 2):
#            indicies = np.where(trackstrip[i] == trackstrip)[0]
#            ind_a = i
#            ind_b = indicies[0] if indicies[0] != i else indicies[1]
#            lum_a = integratecircle(x[ind_a],y[ind_a],10,braggimg)
#            lum_b = integratecircle(x[ind_b],y[ind_b],10,braggimg)
#
#            if (lum_a > lum_b):
#                deltaz[i] = peaks[int(trackstrip[i])]*driftvel
#            else:
#                if (secondpeaks[int(trackstrip[i])] != 0):
#                    print("SECONDARY PEAK USED")
#                    deltaz[i] = secondpeaks[int(trackstrip[i])]*driftvel
#                else:
#                    deltaz[i] = peaks[int(trackstrip[i])]*driftvel
#        elif (d[trackstrip[i]] > 2):
#            deltaz[i] = peaks[int(trackstrip[i])]*driftvel
#            print("CONFUSED, MULTIPLE CROSSINGS DETECTED")
        

#print(deltaz)
fig = plt.figure()
ax = plt.axes(projection='3d')
scatter_plot = ax.scatter3D(x,y,deltaz,c=deposition,s=deposition*150,cmap='turbo')

output = np.zeros((4,len(x)))

output[0] = x
output[1] = y
output[2] = deltaz
output[3] = deposition
np.savetxt("3Dtrack.txt",output)

ax.set_xlabel("x [pix]")
ax.set_ylabel("y [pix]")
ax.set_zlabel("z [cm]")
#ax.set_zlim(-0.8,-0.4)
#ax.set_xlim(220,320)
#ax.set_ylim(540,600)
plt.colorbar(scatter_plot, label="Charge deposited to strip [Au]")
plt.show()

