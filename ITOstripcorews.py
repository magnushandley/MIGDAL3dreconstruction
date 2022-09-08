import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from PIL import Image
from PIL.TiffTags import TAGS
from scipy import interpolate
from scipy.signal import find_peaks

class ITO:

    def __init__(self, filepath, maxsamples, eventstoconsider):
        self.numevents = eventstoconsider
        self.maxsamples = maxsamples
        self.filepath = filepath
        events = eventstoconsider + 1
        num_channels = 68
        self.channels = num_channels
        
        self.data = np.zeros((int(eventstoconsider),int(num_channels),2,maxsamples))
  
  
    def fullread(self):
        """
        Returns all data from the strips and PMT's as a numpy array, up until the event 'eventstoconsider'
        """
        deltas = []
        daqtimestamps = []
        prevtimestamp = 0
        
        maxsamples = self.maxsamples
        eventstoconsider = self.numevents
        filepath = self.filepath
        text_file = open(filepath, "r")
        
        filename = text_file.readline()
        events = int(text_file.readline())
        num_channels = int(text_file.readline())
        
        data = np.zeros((int(eventstoconsider),int(num_channels),2,maxsamples))

        for k in range(eventstoconsider):
            eventnum = text_file.readline()
            try:
                #timestamp given in seconds in file, converting to ms
                timestamp = float(text_file.readline())*1e3
                daqtimestamps.append(timestamp)
            except Exception as e:
                print(e)
                print(text_file.readline())
                break
                    
            if (k != 0):
                dt = timestamp - prevtimestamp
                deltas.append(dt)
                    
            prevtimestamp = timestamp
            for l in range(num_channels):
                channel = text_file.readline()
                samples = text_file.readline()
                    
                timerow = text_file.readline().split()
                    
                if (len(timerow) != maxsamples):
                    timerow = np.pad(timerow, (0, maxsamples - len(timerow)), 'constant')
                    
                data[k][l][0] = timerow

                amprow = text_file.readline().split()
                
                if (len(amprow) != maxsamples):
                    amprow = np.pad(amprow, (0, maxsamples - len(amprow)), 'constant')
    
                data[k][l][1] = amprow
                
        return(data)
 
    def readnextevent(self):
        """
        returns the data from the next event (from the line previously read) as a numpy array 'eventdata'
        """
        num_channels = self.channels
        text_file = self.file
        maxsamples = self.maxsamples
        eventnum = text_file.readline()
        eventdata = np.zeros((int(num_channels),2,maxsamples))
        
        try:
            #timestamp given in seconds in file, converting to ms
            timestamp = float(text_file.readline())*1e3
        except Exception as e:
            print(e)
            print(text_file.readline())
                    
        prevtimestamp = timestamp
        
        for l in range(num_channels):
            channel = text_file.readline()
            samples = text_file.readline()
                    
            timerow = text_file.readline().split()
                    
            if (len(timerow) != maxsamples):
                timerow = np.pad(timerow, (0, maxsamples - len(timerow)), 'constant')
                    
            eventdata[l][0] = timerow

            amprow = text_file.readline().split()
                
            if (len(amprow) != maxsamples):
                amprow = np.pad(amprow, (0, maxsamples - len(amprow)), 'constant')
    
            eventdata[l][1] = amprow
            
        return(eventdata)
        
    def lookuptable(self,offset,imagedir):
        """
        Generates a lookup table (by printing) that relates daq events to the point of the rolling shutter within an image at the moment of the daq trigger, based on a known offset.
        Inputs
        --------
        Offset: The offset between the daq timestamp and the camera timestamp, in ms (~18 for most of the data looked at)
        imagedir: The directory containing the imagefiles
        Outputs: Will print lines consisitng of a string with the format daq event number + , + rolling shutter position.
        """
        
        daqtimestamps = []
        imgtimestamps = []
        
        maxsamples = self.maxsamples
        eventstoconsider = self.numevents
        filepath = self.filepath
        
        text_file = open(filepath, "r")
        
        filename = text_file.readline()
        events = int(text_file.readline())
        num_channels = int(text_file.readline())
        
        for k in range(eventstoconsider):
            eventnum = text_file.readline()
            try:
                #timestamp given in seconds in file, converting to ms
                timestamp = float(text_file.readline())
                daqtimestamps.append(timestamp)
                    
            except Exception as e:
                print(e)
                print(text_file.readline())
                break

            for l in range(num_channels):
                channel = text_file.readline()
                samples = text_file.readline()
                timerow = text_file.readline().split()
                amprow = text_file.readline().split()
          
        imgfiles = sorted(glob.glob(imagedir+'*.TIFF'))


        for i in range(len(imgfiles)):
            with Image.open(imgfiles[i]) as img:
                meta_dict = {TAGS[key] : img.tag[key] for key in img.tag_v2}
                imgtime = float(meta_dict['ImageDescription'][0])-offset
                imgtimestamps.append(imgtime)

        timedifferences = np.subtract(imgtimestamps[1:],imgtimestamps[:-1])
        imgnumbers = np.arange(len(imgfiles))
        print(len(imgtimestamps),len(imgnumbers))
        interpolatednumbers = interpolate.interp1d(imgtimestamps, imgnumbers)

        minimgtime = np.min(imgtimestamps)
        maximgtime = np.max(imgtimestamps)

        for i in range(eventstoconsider):
            ts = daqtimestamps[i]
            if ((ts < minimgtime) or (ts > maximgtime)):
                if (ts < minimgtime):
                    delim = "<"
                if (ts > maximgtime):
                    delim = ">"
#                print("Out of range: "+str(minimgtime)+" "+str(ts)+" "+str(maximgtime)+delim)
            else:
                #outputs EVENT NUMBER, not index - subtract 1 when indexing later.
                print(i+1,',',interpolatednumbers(ts)+1)
                
        return(True)
    
    def readevent(self,event):
        """
        Reads a single event to the overall numpy array
        Inputs:
        -------
        event: index of the event to read
        Outputs:
        -------
        Writes the event data into data[event]
        """
        maxsamples = self.maxsamples
        eventstoconsider = self.numevents
        filepath = self.filepath
        text_file = open(filepath, "r")
        
        filename = text_file.readline()
        events = int(text_file.readline())
        num_channels = int(text_file.readline())
        
        data = np.zeros((int(eventstoconsider),int(num_channels),2,maxsamples))
        for k in range(event+1):
                eventnum = text_file.readline()
                try:
                    #timestamp given in seconds in file, converting to ms
                    timestamp = float(text_file.readline())*1e3
                except Exception as e:
                    print(e)
                    print(text_file.readline())
                    break
                    
                for l in range(num_channels):
                    channel = text_file.readline()
                    samples = text_file.readline()
                    
                    timerow = text_file.readline().split()
                    
                    if (len(timerow) != maxsamples):
                        timerow = np.pad(timerow, (0, maxsamples - len(timerow)), 'constant')

                    if (k==event):
                        data[k][l][0] = timerow

                    amprow = text_file.readline().split()
                    
                    if (len(amprow) != maxsamples):
                        amprow = np.pad(amprow, (0, maxsamples - len(amprow)), 'constant')

                    if (k==event):
                        data[k][l][1] = amprow
                        
                self.data = data
                
        return(data)

    def displayfullevent(self,event):
        """
        Displays the ITO strip data for an event as a heatmap, requires readevent to have been run before, and plt.plot() to be run after.
        """
        data = self.data
        eventdata = np.zeros((60,1000))
        for i in range(60):
            eventdata[i] = data[event][i][1]
            offsetmean = np.mean(eventdata[i][:200])
            eventdata[i] = np.subtract(eventdata[i], offsetmean)
        eventdata = np.transpose(eventdata)
        plt.pcolormesh(eventdata, cmap='magma')
        plt.ylabel("Time [samples]")
        plt.xlabel("Strip number")
        plt.title("ITO strip readout")

    def displaygivenevent(self,data):
        """
        Displays an event passed to this function, rather than from the main 'data' array. Used if you aren't using the class in quite the normal way.
        """
        eventdata = np.zeros((60,1000))
        for i in range(60):
            eventdata[i] = data[i][1]
        plt.pcolormesh(eventdata, cmap='inferno')
        plt.xlabel("Time [samples]")
        plt.ylabel("Strip number")
        plt.title("Dataset ...2221 event 3055")
        
    def peaktime(self,event):
        """
        Returns the sample number at which each ITO strip waveform is at a maximum - used for the simple 3D reconstruction.
        Inputs:
        --------
        event: Index of event to look at
        
        Outputs:
        --------
        maxtimes: a 1d 60 element numpy array containing the time of the peak for each strip, in microseconds.
        """
        data = self.data
        maximums = np.zeros(60)
        maxtimes = np.zeros(60)
        
        for i in range(60):
            max = np.amax(data[event][i][1])
            maximums[i] = max
            index = int(np.where(data[event][i][1] == max)[0][0])
            maxtimes[i] = data[event][i][0][index]
            
        return(maxtimes)
        
    def multipeaktime(self,event):
        """
        Returns the sample number(s) at which each ITO strip is at a maximum
        """
        data = self.data
        maxtimes = np.zeros(60)
        secondarypeaktimes = np.zeros(60)
        
        for i in range(60):
            rowdata = data[event][i][1]
            std = np.std(rowdata[:200])
            peaks = find_peaks(data[event][i][1],10*std,prominence = 10*std)
            heights = peaks[1]['peak_heights']
            prominences = peaks[1]['prominences']
            indices = peaks[0]
            
            if (len(heights) > 0):
                maxtimes[i] = data[event][i][0][indices[np.where(heights == np.amax(heights))[0][0]]]
                if (len(heights) > 1):
                    #sort peaks by prominence, and assume that the secondary peak is either the most prominent or second most prominent.
                    sortedprominences = np.sort(prominences)
                    index1 = np.where(sortedprominences[0] == prominences)[0][0]
                    index2 = np.where(sortedprominences[1] == prominences)[0][0]
                    if (heights[index1] == np.amax(heights)):
                        secondarypeaktimes[i] = data[event][i][0][indices[np.where(heights == heights[index2])[0][0]]]
                    else:
                        secondarypeaktimes[i] = data[event][i][0][indices[np.where(heights == heights[index1])[0][0]]]
                    
            
        return(maxtimes, secondarypeaktimes)
            
        
        
    def stripintegrals(self,event):
        """
        Returns the energy deposition on each strip, as a numpy array
        """
        data = self.data
        rowintegrals = np.zeros(60)
        
        for i in range(60):
            rowdata = data[event][i][1]
            offsetmean = np.mean(rowdata[:200])
            std = np.std(rowdata[:200])
            rowdata = np.subtract(rowdata, offsetmean)
            max = np.amax(rowdata)
            maxindex = int(np.where(rowdata == max)[0][0])
            
            #finds the points at which the signal drops below
            startfound = False
            j = maxindex
            
            while ((startfound == False) and (j >= 0)):
                if rowdata[j] <= 0:
                    startpoint = j
                    startfound = True
                else:
                    j -= 1
              
            endfound = False
            k = maxindex
            while ((endfound == False) and (k<1000)):
                if rowdata[k] <= 0:
                    endpoint = k
                    endfound = True
                else:
                    k += 1
              
            if (max > 5*std):
                integral = np.sum(rowdata[startpoint:endpoint])
                rowintegrals[i] = integral
            else:
                rowintegrals[i] = 0
            
        return(rowintegrals)
        
def pixeltostrip(px, offset):
    """
    Vectorised operation on lists here, not arranged for wraparoud yet
    """
    strip = np.add(0.084375*px, (10.89996 - offset))
    return(strip)
    
def integratecircle(x,y,r,img):
    """
    Takes x,y,r as integers, returns the sum over the value of all coordinates withing r of x,y in an image.
    """
    sum = 0
    xmin = x - r
    xmax = x + r
    ymin = y - r
    ymax = y + r
    xcoords = np.arange(xmin,xmax+1)
    ycoords = np.arange(ymin,ymax+1)

    for i in range(len(ycoords)):
        for j in range(len(xcoords)):
            if (((xcoords[j]-x)**2 + (ycoords[i]-y)**2) <= r**2):
                sum += img[i][j]
            
    return(sum)
    

            
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
