# -*- coding: utf-8 -*-
# testscript2.py
"""
Created on Sun Jan 03 15:13:21 2016

@author: jeff
"""

from osgeo import ogr
import numpy as np
import agfuse

#boundary = np.array([(34.671431, -120.331162),
#                    (34.671426, -120.327284),
#                    (34.670870, -120.327238),
#                    (34.667834, -120.328600),
#                    (34.665236, -120.329596),
#                    (34.664424, -120.330089),
#                    (34.664159, -120.330776),
#                    (34.664282, -120.331291),
#                    (34.671431, -120.331162)])
#                    
#agfuse.createProperty("Melville", "Vineyard", boundary)
#
#shapefile = ogr.Open("Melville")
#layer = shapefile.GetLayer(0)
#
#for i in range(layer.GetFeatureCount()):
#    feature = layer.GetFeature(i)
#    name = feature.GetField("NAME")
#    use = feature.GetField("USE")
#    geometry = feature.GetGeometryRef()
#    print i, name, use, geometry.ExportToWkt()
    
propArray = agfuse.isoRasterProp('..\LC80420362015194LGN00_B4.TIF', 'Melville\Melville.shp')