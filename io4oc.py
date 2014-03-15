#-*- coding: utf-8 -*-
from __future__ import division  
import arcpy
import numpy as np

def input(shpfile, field,cellsize,IniRasFile,bndFile=None):
 if IniRasFile is not None:
  F0, cellsize, extent, nrows,ncols=getRasterProperties(IniRasFile)  
  data,rcIx,env2, nrows2, ncols2= setupEnv(shpfile, field,cellsize,extent)  
 else:
  data,rcIx,extent, nrows, ncols= setupEnv(shpfile, field,cellsize,None,bndFile)
  F0=np.zeros((nrows,ncols))+np.mean(data[:,2])    
 return np.float64(F0), data[:,2],rcIx,extent, nrows, ncols

 
def setupEnv(sFullFilename, sFieldname,cellsize,extent=None,bndFile=None):
 xyzs=arcpy.da.FeatureClassToNumPyArray(sFullFilename, ("SHAPE@X","SHAPE@Y",sFieldname)) #read the shapefile 
 data=np.array([xyzs["SHAPE@X"],xyzs["SHAPE@Y"], xyzs[sFieldname]],dtype=np.float)  #3*pts 
 data=data.T  #pts*3
 
 if extent is None:  
  if bndFile is not None:
   arcpy.env.extent = arcpy.Describe(bndFile).Extent
   tempx=np.arange(arcpy.env.extent.XMin, arcpy.env.extent.XMax,cellsize); 
   tempy=np.arange(arcpy.env.extent.YMin, arcpy.env.extent.YMax,cellsize); 
   arcpy.env.extent = arcpy.Extent(arcpy.env.extent.XMin-0.5*cellsize, arcpy.env.extent.YMin-0.5*cellsize,
                   tempx.max()+cellsize,tempy.max()+cellsize) #set the environment extent 
  else:
   tempx=np.arange(data[:,0].min()-0.5*cellsize, data[:,0].max()+0.5*cellsize,cellsize); 
   tempy=np.arange(data[:,1].min()-0.5*cellsize, data[:,1].max()+0.5*cellsize,cellsize); 
   arcpy.env.extent = arcpy.Extent(data[:,0].min()-0.5*cellsize, data[:,1].min()-0.5*cellsize,
                   tempx.max()+cellsize,tempy.max()+cellsize) #set the environment extent 
  ncols=tempx.shape[0]; nrows=tempy.shape[0]    
 else:
  ncols=np.int(np.round((extent.XMax-extent.XMin)/cellsize))
  nrows=np.int(np.round((extent.YMax-extent.YMin)/cellsize))
  arcpy.env.extent = arcpy.Extent(extent.XMin, extent.YMin, extent.XMax,extent.YMax) #set the environment extent 
  #it is necessary to delete sample point lie out of the scope of raster
    
 extent=arcpy.env.extent
 arcpy.env.overwriteOutput=True    
 #mapping the sampling to the corresponding grid,0,1,nrows-1; 0,1,...,ncols-1 
 rowcol=np.array([np.ceil((data[:,1]-extent.lowerLeft.Y)/cellsize)-1,
                   np.ceil((data[:,0]-extent.lowerLeft.X)/cellsize)-1], dtype=np.int32)
 rowcol=rowcol.T  #pts*2
 index1=(rowcol[:,1]>=0) & (rowcol[:,1]<ncols); 
 rc1=rowcol[index1]; data1=data[index1]
 
 index2=(rc1[:,0]>=0) & (rc1[:,0]<nrows)
 rc2=rc1[index2]; data2=data1[index2] 
 return data2,rc2,extent, nrows, ncols


def getRasterProperties(sIniRasFile):
 iniM=arcpy.RasterToNumPyArray(sIniRasFile) # send all these arcpy to one util.py file 
 ras=arcpy.Raster(sIniRasFile)
 cellsize=float(np.double(arcpy.GetRasterProperties_management(sIniRasFile,"CELLSIZEX")))
 arcpy.env.extent =ras.extent # arcpy.Extent(XMin, YMin, XMax,YMax)#(-3, -3, 3, 3.0) 
 extent=arcpy.env.extent
 return iniM,cellsize,extent,ras.height,ras.width

def Output2raster(ArrayResult, lowerLeft,cellsize, sOutputFilename):
 #output the array with lowerleft to raster(.img) file
 if sOutputFilename is not None:  
  RasterResult = arcpy.NumPyArrayToRaster(np.flipud(ArrayResult), lowerLeft,cellsize, cellsize)  
  RasterResult.save(sOutputFilename)
 else:
  print sOutputFilename + ' not exist' 