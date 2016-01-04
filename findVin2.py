# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 21:34:16 2015

@author: Jeff
"""

# import needed packages
import math
from osgeo import osr, ogr, gdal
import numpy as np
from matplotlib import pyplot as plt
import cv2

# Enable GDAL/OGR exceptions
gdal.UseExceptions()

# init variables: spatial references, vineyard polygon corner coords
# Thompson Winery coordinate constants
#lat = 34.750945 # Vineyard latitude
#lng = -120.210923 # Vineyard longitude
# Melville Winery coordinate constants
lat = 34.671714
lng = -120.331511
#tifULX = -105.150271116
#tifULY = 39.7278572773
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

def isolateProperty(lat, lng, src, boundary):
#This function takes a latitude, longitude, source raster image and a mask 
# to locate a specific property in a Landsat scene.  Once the property is 
# located the image pixel data is converted to an array and returned for 
# processing in other functions. Currently, this function only utilizes 
# Landsat 8 data and utilizes a spatial reference for WGS84 zone 11N 
# (EPSG: 32611).
# Inputs:
# - lat: Latitude in decimal degrees (dd.dddddd)
# - lng: Longitude in decimal degrees (dd.dddddd)
# - src: Source image.  A GeoTIFF raster image taken by Landsat 8.  Contains 
#   the needed metadata to complete the calcuation 
# - mask: The data that describes the shape of the property being located.  The
#   mask will zero out any data that is outside of the mask after the property
#   is found
# Output:
# - propArray: The property array is an array of size that is dependent on the
#   size of the located property.  Each cell contains a value equivalent to the
#   pixel value the corresponding location in the raster band
# Jeff Guido, Dec 2015

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
    # FUTURE EXTENSIBILITY: make the outputEPSG be set to a value that 
    #  corresponds to the zone the Landsat image was taken from.  The current 
    #  EPSG corresponds only to zone 11. A suggestion would be to read the 
    #  GeoTiff metadata for the zone and set the outputEPSG accordingly
    
    # Create the geometry associated with the property, create a linear ring
    #  that outlines the property in the raster image.  This results in quicker
    #  processing than searching the entire image for the property shape
    propRing = ogr.Geometry(ogr.wkbLinearRing)
#   to be implemented later:
#   it = np.nditer(boundary, flags = ['multi-index'])
#   while not it.finished:        

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
#    print propPoly.ExportToWkt() # print polygon values before transform

#    propPoint = ogr.Geometry(ogr.wkbPoint) # create point object
#    propPoint.AddPoint(lng, lat) # add the point coordinates to object
    # NOTE: the order of lng and lat in the above function are important
    
    propPoly.Transform(coordTransform)  
    
#    print propPoly.ExportToWkt() # print polygon values after transform
    
#    propPoint.Transform(coordTransform) # perform the coordinate transform
#    # DEBUG START: Print for debugging purposes
#    print propPoint.GetX(), propPoint.GetY()
#    # DEBUG END 
    
    # load the source band raster files
    tif = gdal.Open(src, gdal.GA_ReadOnly)
 
    # find the pixel coordinate
    trans = tif.GetGeoTransform()
#    print trans
    
    # Now set the image boundaries (a rectangle that includes the vineyard)
    # First find the max and min boundary values and add a buffer to each side
    maxCoords = np.amax(boundary, axis=0)
    minCoords = np.amin(boundary, axis=0)
#    pixelX = (propPoint.GetX()-trans[0])/trans[1]
#    pixelY = (propPoint.GetY()-trans[3])/trans[5]
#    pixelArray = np.array([int(math.floor(pixelY)),int(math.floor(pixelX))])
#    print pixelArray
#    
#    # create raster channels that only have the vineyard of interest (VOI) (everything 
#    #  else cropped out)
#    preimg = cv2.imread(src,0)
#    zoom = 40
#    img = preimg[pixelArray[0]-zoom:pixelArray[0]+zoom, pixelArray[1]-zoom:pixelArray[1]+zoom]
#    fig = plt.figure(1)
#    plt.imshow(img, cmap = 'gray', interpolation = 'none')
#    plt.xticks([]), plt.yticks([]) # to hide tick values on X and Y axis
#    fig.show()
    return # return nothing at the moment
    
    # Create the vineyard image mask
    # perform edge canny detection on both raster channels to find vineyard boundaries
    # between the two boundaries determined find a way to choose the most appropriate boundary
    # do something to create the mask from the edge detection (fill in the boundaries, etc)
    
    # Use the mask to set all values in each raster channel to zero except for what is 
    #  under the mask
    # convert the rasters to arrays

isolateProperty(lat, lng, src, boundary)

# convert each non-zero pixel value to a reflectance
# perform the RGI calculation between red and green
# sum the remaining array values and divide by the number of non-zero pixels

# return the cumulative RGI value