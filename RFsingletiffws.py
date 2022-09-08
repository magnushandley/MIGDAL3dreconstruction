import RFFunctions as RF
from astropy.io import fits
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
#from nearlocus import *

pathname = "/Users/magnus/Documents/MIGDAL/deconvolution/betatoreconstruct.tif"
pathnameraw = "/Users/magnus/Documents/MIGDAL/deconvolution/betatoreconstructnodark.tif"

SIGMA = 2.5 #sigma for derivative determination ~> Related to track width
lthresh = 10 #tracks with a response lower than this are rejected (0 accepts all)
uthresh = 0 #tracks with a response higher than this are rejected (0 accepts all)
minlen = 6 #minimum track length accepted
linkthresh = 1 #maximum distance to be linked
logim = False

def singletiff(pathname,sigma, lt, ut, minlen, linkthresh, logim = False):
    '''
    Create a figure containing the image of the particle track and its
    ridgeline. Two plots are produced in this figure. One is in linear space
    and the other is in log space. On both, the unlinked and linked ridges are
    plotted. The linked are plotted as a thin white line overlayed on the
    colored points of the unlinked lines.

    This can be iterated through.

    Parameters
    ----------
    pathname : str
        Location of the image files to be analysed.
    sigma : float
        Sigma for derivative determination (somehow relate to track width).
    lt : float
        Lower threshold for the ridgefinding algorithm.
        This excludes tracks whose hessian eigenvalues fall below lthresh.
    ut : float
        Upper threshold for the ridgefinding algorithm.
        This excludes tracks whose hessian eigenvalues exceed uthresh.
    minlen : int
        Minimum track length accepted.
    linkthresh : int
        The maximum distance between endpoints allowed for linking ridges.
    logim : BOOL, optional
        Do you want the RidgeFinder to operate on the image in log space?
        The default is False.

    Returns
    -------
    None.

    '''
    ##This will create a plot of tracks and ridge for a single image.
    ##It's meant to be used in cycling through plots but can be run alone.

    ##Get Raw images##
      
    img = io.imread(pathname)
    imnumb = os.path.basename(pathname) ## Get image name
    img = img/20 ## Scale the image


    ################################################################
    ##Initialize matplotlib
    plt.close()
    plt.rc('font', size=10)



    ################################################################
    #These are the Parameters to change to work with the RF        #
    ################################################################

    SIGMA = sigma #sigma for derivative determination ~> Related to track width
    lthresh = lt #tracks with a response lower than this are rejected (0 accepts all)
    uthresh = ut #tracks with a response higher than this are rejected (0 accepts all)
    minlen = minlen #minimum track length accepted
    linkthresh = linkthresh #maximum distance to be linked
    # thresh = 0.0001

    ## Logarithmically scales the image
    if logim:
        c = 255/np.log(1+np.max(img))
        img2 = c *(np.log(img+1))
    else:
        img2 = img

    c = 255/np.log(1+np.max(img))
    img3 = c *(np.log(img+1))


    ##Run Ridgefinder
    px, py, nx, ny, eigvals, valid = points_out = RF.find_points(img2, sigma=SIGMA, l_thresh = lthresh, u_thresh=uthresh)
    lines_before, junctions = RF.compose_lines_from_points(points_out)
    ##Link the Ridges
    nlines = RF.linklines(lines_before,minlen,linkthresh)
    lines=lines_before

    ##Set some plot properties
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,8), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(str(imnumb))

    ax1.imshow(img*100, cmap="magma")
    ax1.set_title('Linear Scaling')
#    ax2.imshow(img3, cmap="magma")
#    ax2.set_title('Logarithmic Scaling')
#    lim =RF.get_lines_bounding_box(lines)
    ax2.set_title('Bragg Curve')
    lim =RF.get_lines_bounding_box(lines)


    ##Run through all ridges found in image

    ##Create and plot the splinefit for all unlinked ridgepoints
#    for i, line in enumerate(lines):
#        if len(line[1]) > minlen:
#
#           ax1.set_xlim(lim)
#           ax1.set_ylim(lim)
##          ax2.set_xlim(lim)
##          ax2.set_ylim(lim)
#
#           x = px[line[1], line[0]]
#           y = py[line[1], line[0]]
#
#           ##Get the splinefit for the image
#           try:
#               new_points,der_points = RF.getspline(x,y,ss=1)
#               ax1.plot(new_points[1],new_points[0],'.')
##              ax2.plot(new_points[1],new_points[0],'.')
#           except Exception as e: print(e)

    ##Create and plot the splinefit for all linked ridgepoints
#    for i, line in enumerate(nlines):
#        if len(line[1]) > minlen:
#
#                x = px[line[1], line[0]]
#                y = py[line[1], line[0]]
#
#                ##Get the splinefit for the image
#                try:
#                    new_points1,der_points1 = RF.getspline(x,y,ss=2)
#                    ax1.plot(new_points1[1],new_points1[0],'-',color='white')
#                    ax2.plot(new_points1[1],new_points1[0],'-',color='white')
#                except Exception as e: print(e)

    for i, line in enumerate(nlines):
        if len(line[1]) > minlen:
            
            x = px[line[1], line[0]]
            y = py[line[1], line[0]]
                
                
                
                
            ##Get the splinefit for the image
            print(x,y)
            new_points1,der_points1 = RF.getspline(x,y,ss=2)
            Bragg = RF.simplebragg(new_points1[0],new_points1[1],img)
            
            roughsum = np.sum(Bragg)
            peak = np.where(Bragg == max(Bragg))[0][0]
            if peak < len(Bragg)/2:
            #Init Dir is at the end of the line list, so we reverse it
                new_points1[1]=new_points1[1][::-1]
                new_points1[0]=new_points1[0][::-1]
            else:
                new_points1[1]=new_points1[1]
                new_points1[0]=new_points1[0]
                
            ax1.plot(new_points1[1],new_points1[0],'-',color='white')
            # ax2.plot(new_points1[1],new_points1[0],'-',color='white')
            Bragg = RF.simplebragg(new_points1[0],new_points1[1],img)

            indir,(x0,y0) = RF.initdir_simdat(new_points1[1],new_points1[0],pix=1)
            slope = abs(np.tan(indir))
            leng = 10
            if indir >= -np.pi and indir < -np.pi/2:
                dy = -leng/(np.sqrt(1+slope**2))
                dx = slope*dy
            elif indir >= -np.pi/2 and indir < 0:
                dy = leng/(np.sqrt(1+slope**2))
                dx = -slope*dy
            elif indir >= 0 and indir < np.pi/2:
                dy = leng/(np.sqrt(1+slope**2))
                dx = slope*dy
            else:
                dy = -leng/(np.sqrt(1+slope**2))
                dx = -slope*dy
            colors = ["red","orange","green","blue","purple","yellow"]
            ax1.arrow(x0,y0,dx,dy,head_width=3, head_length=3,label=str(i))
            ax1.annotate(str(i),(x0,y0),color = "white")
            ax2.plot(Bragg,label=str(i))
            ax2.set_xlabel("Length along Track")
            ax2.set_ylabel("Pixel Intensity")
            # indir = indir*180/np.pi
            ax1.legend()
            ax2.legend()
            print(indir)

                    
    plt.show()
    plt.pause(0.1) ##This command allows the cycling to actually happen
    return()
    
def energydeposition(pathname,sigma, lt, ut, minlen, linkthresh,offsetdistance, logim = False):
    '''
    Create a figure containing the image of the particle track and its
    ridgeline. Two plots are produced in this figure. One is in linear space
    and the other is in log space. On both, the unlinked and linked ridges are
    plotted. The linked are plotted as a thin white line overlayed on the
    colored points of the unlinked lines.

    This can be iterated through.

    Parameters
    ----------
    pathname : str
        Location of the image files to be analysed.
    sigma : float
        Sigma for derivative determination (somehow relate to track width).
    lt : float
        Lower threshold for the ridgefinding algorithm.
        This excludes tracks whose hessian eigenvalues fall below lthresh.
    ut : float
        Upper threshold for the ridgefinding algorithm.
        This excludes tracks whose hessian eigenvalues exceed uthresh.
    minlen : int
        Minimum track length accepted.
    linkthresh : int
        The maximum distance between endpoints allowed for linking ridges.
    logim : BOOL, optional
        Do you want the RidgeFinder to operate on the image in log space?
        The default is False.

    Returns
    -------
    None.

    '''
    ##This will create a plot of tracks and ridge for a single image.
    ##It's meant to be used in cycling through plots but can be run alone.

    ##Get Raw images##
      
    img = io.imread(pathname)
    imgraw = io.imread(pathnameraw)
    imnumb = os.path.basename(pathname) ## Get image name
    img = img/20 ## Scale the image


    ################################################################
    ##Initialize matplotlib
    plt.close()
    plt.rc('font', size=10)



    ################################################################
    #These are the Parameters to change to work with the RF        #
    ################################################################

    SIGMA = sigma #sigma for derivative determination ~> Related to track width
    lthresh = lt #tracks with a response lower than this are rejected (0 accepts all)
    uthresh = ut #tracks with a response higher than this are rejected (0 accepts all)
    minlen = minlen #minimum track length accepted
    linkthresh = linkthresh #maximum distance to be linked
    # thresh = 0.0001

    ## Logarithmically scales the image
    if logim:
        c = 255/np.log(1+np.max(img))
        img2 = c *(np.log(img+1))
    else:
        img2 = img

    c = 255/np.log(1+np.max(img))
    img3 = c *(np.log(img+1))


    ##Run Ridgefinder
    px, py, nx, ny, eigvals, valid = points_out = RF.find_points(img2, sigma=SIGMA, l_thresh = lthresh, u_thresh=uthresh)
    lines_before, junctions = RF.compose_lines_from_points(points_out)
    ##Link the Ridges
    nlines = RF.linklines(lines_before,minlen,linkthresh)
    lines=lines_before

    ##Set some plot properties
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,8), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(str(imnumb))

    ax1.imshow(img, cmap="magma")
    ax1.set_title('Linear Scaling')
#    ax2.imshow(img3, cmap="magma")
#    ax2.set_title('Logarithmic Scaling')
#    lim =RF.get_lines_bounding_box(lines)
    ax2.set_title('Bragg Curve')
    lim =RF.get_lines_bounding_box(lines)

    for i, line in enumerate(nlines):
        if len(line[1]) > minlen:
            
            x = px[line[1], line[0]]
            y = py[line[1], line[0]]
                
                
                
                
            ##Get the splinefit for the image
            new_points1,der_points1 = RF.getspline(x,y,ss=2)
            Bragg = RF.simplebragg(new_points1[0],new_points1[1],img)
#            segmentedbragg(new_points1[0],new_points1[1],4,20,1.0,img)
            
            peak = np.where(Bragg == max(Bragg))[0][0]
            if peak < len(Bragg)/2:
            #Init Dir is at the end of the line list, so we reverse it
                new_points1[1]=new_points1[1][::-1]
                new_points1[0]=new_points1[0][::-1]
            else:
                new_points1[1]=new_points1[1]
                new_points1[0]=new_points1[0]
                
            ax1.plot(new_points1[1],new_points1[0],'-',color='white')
            # ax2.plot(new_points1[1],new_points1[0],'-',color='white')
            Bragg = RF.simplebragg(new_points1[0],new_points1[1],img)
            locus = createlocus(new_points1,offsetdistance)
            
            """
            IndexError occurs here due to range of image, added fudge, will change later
            """
            
            try:
                energy = sumnearby(locus[0],locus[1],imgraw)
            except IndexError:
                energy = 0
            
            print(energy)
            
            indir,(x0,y0) = RF.initdir_simdat(new_points1[1],new_points1[0],pix=1)
            slope = abs(np.tan(indir))
            leng = 10
            
            if indir >= -np.pi and indir < -np.pi/2:
                dy = -leng/(np.sqrt(1+slope**2))
                dx = slope*dy
            elif indir >= -np.pi/2 and indir < 0:
                dy = leng/(np.sqrt(1+slope**2))
                dx = -slope*dy
            elif indir >= 0 and indir < np.pi/2:
                dy = leng/(np.sqrt(1+slope**2))
                dx = slope*dy
            else:
                dy = -leng/(np.sqrt(1+slope**2))
                dx = -slope*dy
            colors = ["red","orange","green","blue","purple","yellow"]
            ax1.arrow(x0,y0,dx,dy,head_width=3, head_length=3,label=str(i))
            #,ec=colors[i],fc=colors[i]
            ax1.annotate(str(i)+' sum='+str(energy),(x0,y0),color = "white")
            ax2.plot(Bragg,label=str(i))
            #color=colors[i]
            ax2.set_xlabel("Length along Track")
            ax2.set_ylabel("Pixel Intensity")
            # indir = indir*180/np.pi
            ax1.legend()
            ax2.legend()
            

                    
    plt.show()
    plt.pause(0.1) ##This command allows the cycling to actually happen
    return()
 
def returnlines(pathname,sigma, lt, ut, minlen, linkthresh, logim = False, trackid = "All"):
    img = io.imread(pathname)
    imnumb = os.path.basename(pathname) ## Get image name
    img = img/20 ## Scale the image


    ################################################################
    ##Initialize matplotlib
    plt.close()
    plt.rc('font', size=10)



    ################################################################
    #These are the Parameters to change to work with the RF        #
    ################################################################

    SIGMA = sigma #sigma for derivative determination ~> Related to track width
    lthresh = lt #tracks with a response lower than this are rejected (0 accepts all)
    uthresh = ut #tracks with a response higher than this are rejected (0 accepts all)
    minlen = minlen #minimum track length accepted
    linkthresh = linkthresh #maximum distance to be linked
    # thresh = 0.0001

    ## Logarithmically scales the image
    if logim:
        c = 255/np.log(1+np.max(img))
        img2 = c *(np.log(img+1))
    else:
        img2 = img

    c = 255/np.log(1+np.max(img))
    
    ##Run Ridgefinder
    px, py, nx, ny, eigvals, valid = points_out = RF.find_points(img2, sigma=SIGMA, l_thresh = lthresh, u_thresh=uthresh)
    lines_before, junctions = RF.compose_lines_from_points(points_out)
    ##Link the Ridges
    nlines = RF.linklines(lines_before,minlen,linkthresh)
    x = []
    y = []
    
    for i, line in enumerate(nlines):
        if (len(line[1]) > minlen):
            
            newx = px[line[1], line[0]]
            newy = py[line[1], line[0]]
            new_points1,der_points1 = RF.getspline(newx,newy,ss=2)
            x.extend(new_points1[1])
            y.extend(new_points1[0])
            
        
    return(np.array(x),np.array(y))

#singletiff(pathname,
#    SIGMA,
#    lthresh,
#    uthresh,
#    minlen,
#    linkthresh,
#    logim)

#x,y = returnlines(pathname,
#    SIGMA,
#    lthresh,
#    uthresh,
#    minlen,
#    linkthresh,
#    logim)

