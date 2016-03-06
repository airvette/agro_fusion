# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 21:44:48 2016

@author: Jeff
"""

# ANALYSIS SCRIPT

from osgeo import ogr
import numpy as np
import agfuse
from matplotlib import pyplot as plt

nirPropArray = agfuse.isoRasterProp('..\LC80420362015194LGN00_B5.TIF', 'Test\Test.shp')
nirPropArrayRad = agfuse.scalePix(nirPropArray, '..\LC80420362015194LGN00_MTL.txt', 5)
redPropArray = agfuse.isoRasterProp('..\LC80420362015194LGN00_B4.TIF', 'Test\Test.shp')
redPropArrayRad = agfuse.scalePix(redPropArray, '..\LC80420362015194LGN00_MTL.txt', 4)
ndviVal, ndviArray = agfuse.getNdvi(redPropArrayRad, nirPropArrayRad)

#test code start
fig3 = plt.figure(3)
plt.imshow(ndviArray, cmap = 'BrBG', interpolation = 'none')
plt.xticks([]), plt.yticks([]) # to hide tick values on X and Y axis
fig3.show()
#test code end