# -*- coding: utf-8 -*-
# agfuse.py: Agro Fusion module
# This module features functions that are to be used for Agro Fusion purposes
"""
Created on Fri Jan 01 16:27:01 2016

@author: Jeff
"""

from osgeo import ogr, osr, gdal
import os, os.path, shutil, mmap
import numpy as np
from matplotlib import pyplot as plt
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
    
def createFileName(platform, sensor, msnNo, orbPath, orbRow, yr, dy, gndID, ver, fileType):
# This function accepts scene identification data for a particular raster file.
#  The function will return a string that is the identifying file name for a
#  raster.  Currently this function is built for LANDSAT 8, but other naming
#  formats can be used in the future.
# Inputs: ---------------------------------------------------------------------
# - platform = the mission name (eg. Landsat)
# - sensor = the sensor configuration used (eg. OLI, TRS, Combined)
# - msnNo = mission number (for Landsat 8, the input would be 8)
# - orbPath = the orbital path that the desired scene was taken from (eg. 042)
# - orbRow = the orbital row of the desired scene (eg. 036)
# - yr = year scene was taken, must be 4 digits (eg. 2015)
# - dy = Julian day scene was taken (eg. 194)
# - gndID = Ground Station Identifier (eg. LGN)
# - ver = scene version (eg. 00)
# - fileType = type of file of the raster image.  This value also determines
#   the extension. For example, a Band 5 file would yield outputs of "_B5" 
#   indicating Landsat 8 Band 5 raster image and a file extension of ".TIF"
# Output: ---------------------------------------------------------------------
# - fileName = the name of the desired file, comprised of all relevant inputs
#   in one string
    
    if platform == 'Landsat' & msnNo == 8: # If a L8 raster
        # Follow L8 file-naming convention
        P = 'L' # 'L' for Landsat
        
        # Determine the sensor abbreviation 
        if sensor == 'OLI': # if only OLI sensor was used
            s = 'O'
        elif sensor == 'TIRS': # if only TIRS sensor was used
            s = 'T'
        else: # if both were used
            s = 'C'
        
        # Determine file type and extension abbreviations
        if fileType == 'Band3': # Green
            ft = '_B3.TIF'
        if fileType == 'Band4': # Red
            ft = '_B4.TIF'
        if fileType == 'Band5': # Near Infrared
            ft = '_B5.TIF'
        if fileType == 'MetaD': # metadata file
            ft = 'MTL.txt'
        if fileType == 'BandQ': # quality band file
            ft = 'BQA.TIF'
        if fileType == 'CkSum': # checksum file
            ft = 'MD5.txt'
        else:
            ft = 'xxx.TIF'            
            # Throw an error
        
        filename = '%s%s%s%s%s%s%s%s%s%s%s' % (P, s, )
            
        


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
    dataset = gdal.Open(raster, gdal.GA_ReadOnly)
    # Retrieve the spatial reference of the raster dataset
    
    # Open the shapefile
    shapesource = ogr.Open(shapefile)
    # Clone the polygon geometry object that represents the property
    workingLayer = shapesource.GetLayer(0)
    workingFeature = workingLayer.GetFeature(0)
    workingGeometry = workingFeature.GetGeometryRef()
    
    # Create the spatial transform
    # Create the source spatial reference and set to the shapefile datum
    srcProj = osr.SpatialReference()
    srcProj.SetWellKnownGeogCS('WGS84') # it is assumed that the shapefile inputs will always be the WGS84 datum
    # Create the destination spatial reference and set to the raster projection
    dstProj = osr.SpatialReference(wkt=dataset.GetProjection())
    # Create the transform object and perform the transform on the shapefile
    transform = osr.CoordinateTransformation(srcProj, dstProj)
    workingGeometry.Transform(transform) # Project geometry points to the raster's spatial reference
    
    # Get max and min x and y values of the transformed geometry
    minGeomX,maxGeomX,minGeomY,maxGeomY = workingGeometry.GetEnvelope()
    # Determine the dimensions and offsets of a bounding box that will contain the entire property
    bufferSize = 5 * 30 # buffer size is the number of pixles multiplied by the raster resolution (30 meters)      
    
    # Get the affine transform from the raster and define the raster bounding box in pixles
    geotransform = dataset.GetGeoTransform()
    minBoxX = int(((minGeomX - bufferSize) - geotransform[0])/geotransform[1])
    maxBoxX = int(((maxGeomX + bufferSize) - geotransform[0])/geotransform[1])
    minBoxY = int(((minGeomY - bufferSize) - geotransform[3])/geotransform[5])
    maxBoxY = int(((maxGeomY + bufferSize) - geotransform[3])/geotransform[5])
    # Using the bounding box, read a portion of the raster image to an array
    #  where each array cell corresponds to a raster pixel
    band = dataset.GetRasterBand(1) # Landsat 8 raster files only have one band per tif, indexing starts at 1
    prePropPix = band.ReadAsArray(minBoxX, maxBoxY, (maxBoxX-minBoxX), (minBoxY-maxBoxY))
    
    #test code start
    fig1 = plt.figure(1)
    plt.imshow(prePropPix, cmap = 'gray', interpolation = 'none')
    plt.xticks([]), plt.yticks([]) # to hide tick values on X and Y axis
    fig1.show()
    #test code end
    
    # Create new affine transform using the bounding box offsets
    detailGeoTrans = [minGeomX - bufferSize, geotransform[1], 0,
                      maxGeomY + bufferSize, 0, geotransform[5]]
    # Set all array cells that are outside of the property to a negative number
    # Create test point object
    testPoint = ogr.Geometry(ogr.wkbPoint)
    testPoint.AddPoint(0,0)
    propPix = np.array(prePropPix)
    for index in np.ndindex(propPix.shape):
        # Determine cell location
        xLoc = (index[1] * detailGeoTrans[1]) + detailGeoTrans[0]
        yLoc = (index[0] * detailGeoTrans[5]) + detailGeoTrans[3]
        # Set point location to (xLoc,yLoc)
        testPoint.SetPoint(0,xLoc,yLoc)
        # if cell is outside of geometry, set to 0
        if not testPoint.Within(workingGeometry):
            propPix[index] = 0
    
    #test code start
    fig2 = plt.figure(2)
    plt.imshow(propPix, cmap = 'gray', interpolation = 'none')
    plt.xticks([]), plt.yticks([]) # to hide tick values on X and Y axis
    fig2.show()
    #test code end    
    
    # Perform cleanup actions
    # Destroy the geometry, band or whatever to free up memory
    workingFeature.Destroy()    
    shapesource.Destroy()
    
    # Return the array
    return propPix
    
def scalePix (propPix, metaFile, colorBand):
    #THIS FUNCTION CURRENTLY ONLY SUPPORTS LANDSAT 8 DATA PRODUCTS    
    #This function takes the array created by isoRasterProp() and scales the
    # pixel values in the array according to the metadata in metaFile to calc
    # the spectral radiance.  The colorBand input defines which color scaling 
    # should be applied.
    #
    # The spectral radiance is computed with the equation:
    #       l_lambda = multFactor * propPix[i,j] + addScale
    # Inputs: ---------------------------------------------------------------------
    # - propPix = the 16-bit integer array produced by isoRasterProp()
    # - metaFile = the Landsat Level 1 product metafile.  The file has naming
    #   convention "*_MTL.txt" and contains datum, projection, reflectance, 
    #   time and date information.  This function scans the metadata file for 
    #   the required info
    # - colorBand = an integer can be entered in this variable.  The only 
    #   acceptable integers are 1-11 and directly correspond to Landsat 8's 
    #   payload bands.  Future functionality may allow for strings as inputs
    #   (i.e. 'red', 'blue', 'nir', etc)
    # Output: ---------------------------------------------------------------------
    # - propPixRad = the input array's values converted to measured spectral
    #   radiance values.  The numpy dtype = float  (64 bit float)
    # Jeff Guido, Jan 2016

    # Create keywords based on colorBand input per 'Landsat8DataUsersHandbook',
    #  pg 69, section 5.1
    multKey = 'RADIANCE_MULT_BAND_' + str(colorBand)   
    addKey = 'RADIANCE_ADD_BAND_' + str(colorBand)   
    
    # Determine values of multFactor and addScale
    file = open(metaFile)
    text = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) # create searchable object
    keyStart = text.find(multKey) # find position of multKey
    text.seek(keyStart) # set mmap pointer to accurate multKey line
    multFactorLine = text.readline()
    # The metafile format for this line is RADIANCE_MULT_BAND_N = X.XXXXE-Y,
    #  we want to isolate 'X.XXXE-Y' and convert it to a float
    multFactorWords = multFactorLine.split(' = ')
    multFactor = float(multFactorWords[1])
    # Perform same operation for addScale
    keyStart = text.find(addKey) # find position of addKey
    text.seek(keyStart) # set mmap pointer to accurate multKey line
    addScaleLine = text.readline()
    addScaleWords = addScaleLine.split(' = ')
    addScale = float(addScaleWords[1])
    text.close() # metaFile is not needed anymore
    
    # Create output array and loop through to generate output
    propPixRad = np.array(propPix, dtype = float)
    for index in np.ndindex(propPixRad.shape):
        if propPix[index] > 0:            
            propPixRad[index] = multFactor*propPix[index] + addScale    

    return propPixRad
    

def getNdvi (redArray, nirArray):
    #This function takes two arrays from isoRasterProp() and computes a 
    # Normalized Difference Vegitation Index (NDVI) for the property rendered 
    # by the arrays.  The inputs are of arrays that are identical in geospatial
    # terms, but were captured in using different parts of the visible and 
    # infrared spectrum.
    #
    # The normalized difference for each pixel is computed with the equation:
    #       NDVI = (nirArray[i,j] - redArray[i,j])/(nirArray[i,j] + redArray[i,j])
    # 
    # After computing the NDVI for each pixel, the pixel values are integrated
    # over the isolated property and a scalar is returned
    # Inputs: ---------------------------------------------------------------------
    # - redArray = the array created by isoRasterProp() that represents the 
    #   property of interest in the red visual band. dtype = float
    # - nirArray = the array created by isoRasterProp() that represents the 
    #   property of interest in the Near Infrared visual band. dtype = float
    # Output: ---------------------------------------------------------------------
    # - ndviVal = the NDVI value for the property of interest.
    # - nvdiArray = the raw array is output so it can be used for visulization
    #   in a GUI.  dtype = float
    # Jeff Guido, Jan 2016

    # Create ndviArray to hold the pixel values   
    ndviArray = np.ndarray(redArray.shape)
    area = 0.0 # property area in pixles
    # Iterate through the red and nir arrays
    for index in np.ndindex(ndviArray.shape):
        # If both arrays have valid values (greater than zero), perform the calc        
        if redArray[index] and nirArray[index] > 0:
            ndviArray[index] = (nirArray[index] - redArray[index])/(nirArray[index] + redArray[index])
            area += 1 # area accumulates the number of pixles in the area of interest
        # If only one array has a valid value in the indexed cell, set a flag
        elif redArray[index] or nirArray[index] > 0:
            ndviArray[index] = -500 # flag set to -500 to indicate an error
            
    # Integrate the values to return
    #ndviVal = np.trapz(np.trapz(ndviArray)) # consider bringing the resolution 30m into this measurement
    # Consider substituting the use of trapz with cumsum and dividing the output by the number of cells   
    ndviVal = (1/(area*30.0))*np.sum(ndviArray)    
    
    return ndviVal, ndviArray