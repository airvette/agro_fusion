# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 06:01:28 2015

@author: Ruthie
"""

# import needed packages
import math
from osgeo import osr, ogr, gdal
import numpy as np
from matplotlib import pyplot as plt
import cv2

# Enable GDAL/OGR exceptions
gdal.UseExceptions()

src = '..\LC80420362015194LGN00_B4.TIF'
src2 = '..\LC80420362015194LGN00_B3.TIF'
mask = 0
boundary = np.array([(34.671431, -120.331162),
                    (34.671426, -120.327284),
                    (34.670870, -120.327238),
                    (34.667834, -120.328600),
                    (34.665236, -120.329596),
                    (34.664424, -120.330089),
                    (34.664159, -120.330776),
                    (34.664282, -120.331291),
                    (34.671431, -120.331162)])
                    
# Set up the shapefile driver
driver = ogr.GetDriverByName("ESRI Shapefile")

# Create data source
dataSource = driver.CreateDataSource("Melville_Vineyard.shp")
                    
# Initialize spatial reference (SR) system so the lat long inputs can be 
# converted to machine readable coordinates
inputEPSG = 4326 # wgs84 lat long
outputEPSG = 32611 # wgs84 / UTM zone 11N
# Create spatial reference coordinate transform
inSpatialRef = osr.SpatialReference() # init input SR object
inSpatialRef.ImportFromEPSG(inputEPSG) # set SR object to specific coord sys
outSpatialRef = osr.SpatialReference() # init output SR object
outSpatialRef.ImportFromEPSG(outputEPSG) # set SR object to specific coord sys
# Create coordinate transform using input and output references
coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

# Create layer
layer = dataSource.CreateLayer("Vineyard", srs = outSpatialRef, geom_type = ogr.wkbPolygon)

# Create feature
vineRows = ogr.Feature(layer.GetLayerDefn())
vineRows.SetField("Name", "Vineyard Rows")

# Create the geometry associated with the property, create a linear ring
#  that outlines the property in the raster image.  This results in quicker
#  processing than searching the entire image for the property shape
propRing = ogr.Geometry(ogr.wkbLinearRing)
# to be implemented later:
# it = np.nditer(boundary, flags = ['multi-index'])
# while not it.finished:        

propRing.AddPoint(boundary[0,1],boundary[0,0])
propRing.AddPoint(boundary[1,1],boundary[1,0])
propRing.AddPoint(boundary[2,1],boundary[2,0])
propRing.AddPoint(boundary[3,1],boundary[3,0])
propRing.AddPoint(boundary[4,1],boundary[4,0])
propRing.AddPoint(boundary[5,1],boundary[5,0])
propRing.AddPoint(boundary[6,1],boundary[6,0])
propRing.AddPoint(boundary[7,1],boundary[7,0])
propRing.AddPoint(boundary[8,1],boundary[8,0])

propPoly = ogr.Geometry(ogr.wkbPolygon)
propPoly.AddGeometry(propRing)
#print propPoly.ExportToWkt() # print polygon values before transform

# NOTE: the order of lng and lat in the function are important

propPoly.Transform(coordTransform)  

#print propPoly.ExportToWkt() # print polygon values after transform

vineRows.SetGeometry(propPoly)