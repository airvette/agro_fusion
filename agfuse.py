# -*- coding: utf-8 -*-
# agfuse.py: Agro Fusion module
# This module features functions that are to be used for Agro Fusion purposes
"""
Created on Fri Jan 01 16:27:01 2016

@author: Jeff
"""

from osgeo import ogr, osr, gdal
import os, os.path, shutil
import numpy as np
gdal.UseExceptions()

def createProperty (name, use, boundaryPoints):
#This function creats a shapefile defining the name and coordinates of a 
# vineyard of interest.  The inputs accept a name string, a numpy array 
# containing coordinates and a destination file path.  If the file path does 
# not exist it is created and if files already exist in the folder they are
# deleted.  The intent of this function is to create the required shapfiles
# for use in vineyard analysis and for identifying a vineyard's boundary in a
# raster image
# Inputs: ---------------------------------------------------------------------
# - name = a string containing the name of the vineyard.  Will be added to the
#   appropriate name field in metadata
# - use = the type of use the property has (i.e."Agriculture", "Energy", etc)
# - boundaryPoints = a numpy array.  The expected input format for N points is:
#   boundaryPoints =   [[point 0 lat, point 0 long],
#                       [point 1 lat, point 1 long],
#                       [point 2 lat, point 2 long],
#                       ...
#                       [point N lat, point N long],
#                       [point 0 lat, point 0 long]] 
#   Where the last point closes the ring.  Lat and Long are using the WGS84
#   datum
# Output: ---------------------------------------------------------------------
# - no variable is returned, but the shape file is created in a folder with the
#   value the 'name' input variable
# Jeff Guido, Dec 2015

    # check for a folder that matches 'name'
    # if the directory exists
    if os.path.exists(name):
        # if the files exist delete them
        shutil.rmtree(name)
    os.mkdir(name)
    
    # Create shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile") # use the ESRI format
    
    # Create shapefile datasource
    dstFile = driver.CreateDataSource(name)
    
    # Create spatial reference and shapefile layer
    spatialRef = osr.SpatialReference()
    spatialRef.SetWellKnownGeogCS('WGS84')  
    dstLayer = dstFile.CreateLayer(name, spatialRef)
    
    # Set the metadata field types for features
    # Create Property Name field    
    fieldDef = ogr.FieldDefn("NAME", ogr.OFTString)
    fieldDef.SetWidth(50)
    dstLayer.CreateField(fieldDef)
    
    # Create Property use field (i.e. "Agriculture", "Energy", etc)    
    fieldDef = ogr.FieldDefn("USE", ogr.OFTString)
    fieldDef.SetWidth(50)
    dstLayer.CreateField(fieldDef)
    
    # Create boundary geometry
    # Create a linear ring
    boundary = ogr.Geometry(ogr.wkbLinearRing)
    # for each point in the boundary points variable
    for row in boundaryPoints:
        # create and add the point to the linear ring
        boundary.AddPoint(row[1],row[0]) #NOTE: indexes 0 and 1 are flipped (long comes first)
    # convert ring to polygon
    prop = ogr.Geometry(ogr.wkbPolygon)
    prop.AddGeometry(boundary)
    
    # Set feature parameters
    feature = ogr.Feature(dstLayer.GetLayerDefn()) # get layer defn from above
    feature.SetGeometry(prop) # add the property to the feature
    # Add vineyard name and type attributes from input variables
    feature.SetField("NAME", name) 
    feature.SetField("USE", use)
    # Create feature by setting it to this layer
    dstLayer.CreateFeature(feature)
       
    # close out function by destroying the data source and the feature
    feature.Destroy()
    dstFile.Destroy()

    return None
    
def isoRasterProp(raster, shapefile):
#This function takes a Landsat raster image and uses a shapefile to return a 
# numpy array that only includes the raster pixel values that are inside of the 
# property as identified by the shapefile.  The output is intended for 
# processing by other functions in this module
# Inputs: ---------------------------------------------------------------------
# - raster = the file path and file name of the raster image to be processed
# - shapefile = the file path and file name of the shapefile to be used for
#   this function.  The .shp file input needs to have been created by the 
#   createProperty() function that is part of this module (agfuse.py)
# Output: ---------------------------------------------------------------------
# - propPix = the raster pixels representing the property.  This will be a 2D
#   numpy array.  Since it is unlikely that the property will be rectangular, 
#   all values that are not inside the property will have a negative value.
#   Positive values are to be interpreted as valid pixels that are inside the
#   property
# Jeff Guido, Jan 2016

    # Utility functions that check if the files actually exist
    # Check the filepath of the raster
    # if the filepath does not exist
        # exit the program and throw and exception
    # Check the filepath of the shapefile
    # if the filepath does not exist
        # exit the program and throw and exception
    # Check to make sure that the shapefile coordinates are within the raster
    #  boundaries
        # exit program and throw execption if the shapefile is outside of the
        #  raster 
    
    # Open the raster
    # Retrieve the spatial reference of the raster
    
    # Open the shapefile
    # Clone the polygon geometry object that represents the property
    # Project geometry points to the raster's spatial reference
    # Get max and min x and y values of the transformed geometry
    # Determine the dimensions and offsets of a bounding box that will contain
    #  the entire property
    
    # Get the affine transform from the raster
    # Using the bounding box, read a portion of the raster image to an array
    #  where each array cell corresponds to a raster pixel
    # Create new affine transform using the bounding box offsets
    
    # Set all array cells that are outside of the property to a negative number
    
    # Perform cleanup actions
    # Destroy the geometry, band or whatever to free up memory    
    
    # Return the array
    return propPix