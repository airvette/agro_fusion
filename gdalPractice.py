# -*- coding: utf-8 -*-
"""
Created on Sat Jan 09 08:31:43 2016

@author: Ruthie
"""

# gdal practice
from osgeo import gdal, ogr, osr
gdal.UseExceptions()

# Open the raster
dataset = gdal.Open('..\LC80420362015194LGN00_B4.TIF', gdal.GA_ReadOnly)
# Retrieve the spatial reference of the raster dataset
dstProj = osr.SpatialReference(wkt=dataset.GetProjection())

shapefile = ogr.Open('Melville\Melville.shp')
layer = shapefile.GetLayer(0)
feature = layer.GetFeature(0)
geometry = feature.GetGeometryRef()
#print geometry.ExportToWkt()
srcProj = osr.SpatialReference()
srcProj.SetWellKnownGeogCS('WGS84')
#print "Shapefile datum is ", layer.GetSpatialRef()

transform = osr.CoordinateTransformation(srcProj, dstProj)
geometry.Transform(transform)
print geometry.ExportToWkt()
minX,maxX,minY,maxY = geometry.GetEnvelope()
print minX,maxX,minY,maxY 
