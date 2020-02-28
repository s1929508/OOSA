from lvisClass import lvisData
from myhandletif import tiffhandle 

import argparse
import numpy as np
import gdal,osr

def ToTiff(data,x,y,res,filename="lvis_image.tif"):
  '''
  Make a geotiff from an array of points
  '''
  #print(x)
  # determine bounds
  minX=np.min(x)
  maxX=np.max(x)
  minY=np.min(y)
  maxY=np.max(y)
  #print(maxX,minX)
  #print(maxX-minX)
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
  srs.ImportFromEPSG(3031)              # WGS84 lat/long
  dst_ds.SetProjection(srs.ExportToWkt()) # export coords to file
  dst_ds.GetRasterBand(1).WriteArray(imageArr)  # write image to the raster
  dst_ds.GetRasterBand(1).SetNoDataValue(-999)  # set no data value
  dst_ds.FlushCache()                     # write to disk
  dst_ds = None
  gdal.SetConfigOption('HFA_USE_RRD', 'YES')#Set Config Option
  print("Image written to",filename)
  return

def readCommands():
  '''
  Read commandline arguments
  '''
  p = argparse.ArgumentParser(description=("Handle a set of points in geopandas"))
  p.add_argument("--input", dest ="filename", type=str, default='/geos/netdata/avtrain/data/3d/oosa/assignment/lvis/2015/ILVIS1B_AQ2015_1017_R1605_043439.h5', help=("Input filename"))
  p.add_argument("--output", dest ="Outname", type=str, default='ILVIS1B_AQ2015_1017_R1605_043439.h5', help=("Output filename"))
  p.add_argument("--res", dest ="res", type=int, default=10, help=("The resolution"))
  cmdargs = p.parse_args()
  return cmdargs

      
if __name__=="__main__":
  '''Main block'''
  cmd=readCommands()
  filename=cmd.filename
  Outname=cmd.Outname
  res=cmd.res
  zG=[]
  lon=[]
  lat=[]
  # find bounds
  b=lvisData(filename,onlyBounds=True)
  #Set the bounds of this pics
  gap = nX = nY = 4
  for i in range(nX):
      for j in range(nY):
          # set some bounds
          x0=(b.bounds[2]-b.bounds[0])*i/gap+b.bounds[0]
          y0=(b.bounds[3]-b.bounds[1])*j/gap+b.bounds[1]
          x1=(b.bounds[2]-b.bounds[0])*(i+1)/gap+b.bounds[0]
          y1=(b.bounds[3]-b.bounds[1])*(j+1)/gap+b.bounds[1]
          lvis=tiffhandle(filename,minX=x0,minY=y0,maxX=x1,maxY=y1,res=res,Outname=Outname,onlyBounds=False)
          #Initialize the lvis file
          p = lvis.readLVIS(filename,minX=x0,minY=y0,maxX=x1,maxY=y1)
          if(p):
              print("There are some data in",[i,j],"bounds!")
              lvis.reproject(4326,3031)
              lvis.estimateGround()
              zG = np.append(zG,lvis.zG)
              lat = np.append(lat,lvis.lat)
              lon = np.append(lon,lvis.lon)
          else:
              print("There is no data in",[i,j],"bounds!")
  #Give all the data of one h5 file and put them into tiff.
  ToTiff(zG,lon,lat,res,filename=Outname+".tif")
