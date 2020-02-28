

'''
My class based on lvisClass.py

Changes are made in readLVIS and writeLVIS

#readLVIS:give a return value to indentify the flag 'p'.
p = Flase: There is no data contained in the choosen region
p = True : There are some data and have been loaded.

#writeLvis:Change the method of reprojection 

'''

from pyproj import Proj, transform # package for reprojecting data
from osgeo import gdal             # pacage for handling geotiff data
from osgeo import osr              # pacage for handling projection information
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d 
import h5py

class tiffhandle():
    def __init__(self,filename,setElev=False,minX=-1000000,maxX=1000000,minY=-10000000,maxY=1000000,res='None',Outname='None',onlyBounds=False):
        '''
        Class initialiser. Calls a function
        to read LVIS data within bounds
        minX,minY and maxX,maxY
        setElev=1 converts LVIS's stop and start
        elevations to arrays of elevation.
        onlyBounds sets "bounds" to the corner of the area of interest
        '''

        
    def readLVIS(self,filename,minX,minY,maxX,maxY):
        '''
        Read LVIS data from file
        '''
        # open file for reading
        f=h5py.File(filename,'r')
        # determine how many bins
        self.nBins=f['RXWAVE'].shape[1]
        # read coordinates for subsetting
        lon0=np.array(f['LON0'])       # longitude of waveform top
        lat0=np.array(f['LAT0'])       # lattitude of waveform top
        lonN=np.array(f['LON'+str(self.nBins-1)]) # longitude of waveform bottom
        latN=np.array(f['LAT'+str(self.nBins-1)]) # lattitude of waveform bottom
        # find a single coordinate per footprint
        tempLon=(lon0+lonN)/2.0
        tempLat=(lat0+latN)/2.0
        
        # dertermine which are in region of interest
        useInd=np.where((tempLon>=minX)&(tempLon<maxX)&(tempLat>=minY)&(tempLat<maxY))
        if(len(useInd)>0):
          useInd=useInd[0]
    
        if(len(useInd)==0):
          print("No data contained in that region")
          return False
    
        # save the subset of all data
        self.nWaves=len(useInd)
        self.lon=tempLon[useInd]
        self.lat=tempLat[useInd]
        
        # load sliced arrays, to save RAM
        self.lfid=np.array(f['LFID'])[useInd]          # LVIS flight ID number
        self.lShot=np.array(f['SHOTNUMBER'])[useInd]   # the LVIS shot number, a label
        self.waves=np.array(f['RXWAVE'])[useInd]       # the recieved waveforms. The data
        self.nBins=self.waves.shape[1]
        # these variables will be converted to easier variables
        self.lZN=np.array(f['Z'+str(self.nBins-1)])[useInd]       # The elevation of the waveform bottom
        self.lZ0=np.array(f['Z0'])[useInd]          # The elevation of the waveform top
        # close file
        f.close()
        # return to initialiser
        return True
      
    #####################################
    
    def writeTiff(self,data,x,y,res,filename="output.tif"):
      '''
      Make a geotiff from an array of points
      '''
      #print(x)
      
      # determine bounds
      minX=np.min(x)
      maxX=np.max(x)
      minY=np.min(y)
      maxY=np.max(y)
      
      # determine image size
      nX=int((maxX-minX)/res+1)
      nY=int((maxY-minY)/res+1)
    
      # pack in to array
      imageArr=np.full((nY,nX),-999.0)        # make an array of missing data flags
      xInds=np.array((x-minX)/res,dtype=int)  # determine which pixels the data lies in
      yInds=np.array((maxY-y)/res,dtype=int)  # determine which pixels the data lies in
      # this is a simple pack which will assign a single footprint to each pixel
      imageArr[yInds,xInds]=data
    
      # set geolocation information (note geotiffs count down from top edge in Y)
      geotransform = (minX, res, 0, maxY, 0, -res)
    
      # load data in to geotiff object
      dst_ds = gdal.GetDriverByName('GTiff').Create(filename, nX, nY, 1, gdal.GDT_Float32)
    
      dst_ds.SetGeoTransform(geotransform)    # specify coords
      srs = osr.SpatialReference()            # establish encoding
      srs.ImportFromEPSG(3031)              # epsg: 3031 to fit the Antarctica
      dst_ds.SetProjection(srs.ExportToWkt()) # export coords to file
      dst_ds.GetRasterBand(1).WriteArray(imageArr)  # write image to the raster
      dst_ds.GetRasterBand(1).SetNoDataValue(-999)  # set no data value
      dst_ds.FlushCache()                     # write to disk
      dst_ds = None
    
      print("Image written to",filename)
      return
        
    def estimateGround(self,sigThresh=5,statsLen=10,minWidth=3,sWidth=0.5):
      self.setElevations()
        # find noise statistics
      self.findStats(statsLen=statsLen)

    # set threshold
      threshold=self.setThreshold(sigThresh)

    # remove background
      self.denoise(threshold,minWidth=minWidth,sWidth=sWidth)

    # find centre of gravity of remaining signal
      self.CofG()


  #######################################################

    def setThreshold(self,sigThresh):

      threshold=self.meanNoise+sigThresh*self.stdevNoise
      return(threshold)


  #######################################################

    def CofG(self):
    # allocate space and put no data flags
      self.zG=np.full((self.nWaves),-999.0)

    # loop over waveforms
      for i in range(0,self.nWaves):
        if(np.sum(self.denoised[i])>0.0):   # avoid empty waveforms (clouds etc)
          self.zG[i]=np.average(self.z[i],weights=self.denoised[i])


  #######################################################

    def reproject(self,inEPSG,outEPSG):
      # set projections
      inProj=Proj(init="epsg:"+str(inEPSG))
      outProj=Proj(init="epsg:"+str(outEPSG))
    # reproject data
      x,y=transform(inProj,outProj,self.lon,self.lat)
      self.lon=x
      self.lat=y


  ##############################################

    def findStats(self,statsLen=10):
    # make empty arrays
      self.meanNoise=np.empty(self.nWaves)
      self.stdevNoise=np.empty(self.nWaves)

    # determine number of bins to calculate stats over
      res=(self.z[0,0]-self.z[0,-1])/self.nBins    # range resolution
      noiseBins=int(statsLen/res)   # number of bins within "statsLen"

    # loop over waveforms
      for i in range(0,self.nWaves):
        self.meanNoise[i]=np.mean(self.waves[i,0:noiseBins])
        self.stdevNoise[i]=np.std(self.waves[i,0:noiseBins])


  ##############################################

    def denoise(self,threshold,sWidth=0.5,minWidth=3):
    # find resolution
      res=(self.z[0,0]-self.z[0,-1])/self.nBins    # range resolution

    # make array for output
      self.denoised=np.full((self.nWaves,self.nBins),0)

    # loop over waves
      for i in range(0,self.nWaves):
        print("Denoising wave",i+1,"of",self.nWaves)

      # subtract mean background noise
        self.denoised[i]=self.waves[i]-self.meanNoise[i]

      # set all values less than threshold to zero
        self.denoised[i,self.denoised[i]<threshold[i]]=0.0

      # minimum acceptable width
        binList=np.where(self.denoised[i]>0.0)[0]
        for j in range(0,binList.shape[0]):       # loop over waveforms
          if((j>0)&(j<(binList.shape[0]-1))):    # are we in the middle of the array?
            if((binList[j]!=binList[j-1]+1)|(binList[j]!=binList[j+1]-1)):  # are the bins consecutive?
              self.denoised[i,binList[j]]=0.0   # if not, set to zero

      # smooth
        self.denoised[i]=gaussian_filter1d(self.denoised[i],sWidth/res)
          
  ##############################################
    
    def setElevations(self):
      self.z=np.empty((self.nWaves,self.nBins))
      for i in range(0,self.nWaves):    # loop over waves
        res=(self.lZ0[i]-self.lZN[i])/self.nBins
        self.z[i]=np.arange(self.lZ0[i],self.lZN[i],-1.0*res)   # returns an array of floats
