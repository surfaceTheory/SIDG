'''
Created on 2013-12-20

@author: SIST
'''
import oc
import io4oc
data='gauss50unif20';   solver='knitro'
shpfile='c:\\HASMCode\\gauss\\'+data+'.shp';    field='GRID_CODE'; cellsize=1.0/50 
out_raster='C:\\HASMCode\\gauss.img'

## example 1
# F1, extent=oc.hasm(shpfile, field, cellsize,IniRasFile=None, reTol=1e-5,
#          options='hasm6',solver=solver,out_raster=out_raster)

##example 2, do the same thing as in the example 1
F1, extent=oc.hasm(shpfile, field, cellsize,IniRasFile=None, reTol=1e-5,
        options='hasm6',solver=solver)
io4oc.Output2raster(F1, extent.lowerLeft,cellsize, out_raster)  