# Readme for agro_fusion

## Overview
This is a repo that is used to store all the agro fusion code prototyping. I used this project to look into using Geographic Information System tools to perform analytics on vineyards. This repo ingests [Landsat 8](http://landsat.usgs.gov/landsat8.php) products, locates the optical bands of interest, scales the output to represent physical values and performs [NDVI](https://en.wikipedia.org/wiki/Normalized_Difference_Vegetation_Index) calculations.  

### Example output
The below images are the outline and an NDVI example output using the code in this repo. The imagery used for this analysis was taken on 07/13/2015. The coordinates for the north east corner of the property are [[34.627310756333465,-120.27924031019211]](https://www.google.com/maps/place/34%C2%B037'38.3%22N+120%C2%B016'45.3%22W/@34.6273152,-120.281429,788m/data=!3m1!1e3!4m5!3m4!1s0x0:0x0!8m2!3d34.6273108!4d-120.2792403). The property is Ampelos Vineyard near Buellton, CA.

**Vineyard Outline:** The image here is used to indicate the outline of the evaluated vineyard. The image is not used for analytics. 

![Vineyard Outline] (https://github.com/airvette/agro_fusion/blob/master/ampelos_geometry.png)

**NDVI Output:** Output from the boundary above. Brown pixels inside the boundary indicate an NDVI close to -1, white pixels are close to 0 and green pixels indicate values close to +1.

![NDVI Output] (https://github.com/airvette/agro_fusion/blob/master/ampelosTestFig(1611).png)

## Contents 
This repo contains the primary module [agfuse.py](https://github.com/airvette/agro_fusion/blob/master/agfuse.py) that when used with [createSphScript.py](https://github.com/airvette/agro_fusion/blob/master/createSphScript.py) and [analysisScript.py](https://github.com/airvette/agro_fusion/blob/master/analysisScript.py) provide rudimentary insights into NDVI calculations for the identified area. While this repo is used for vineyards, these calculations will deliver NDVI results for any imagery. This repo also contains various test scripts. 

## Workflow
Note that to use this package several Python packages need to be installed onto the host. See the [text file](https://github.com/airvette/agro_fusion/blob/master/Python_Packages_Installed.txt) to see what was installed on my system to run the package. Windows 7 (and later 10) was the operating system.

1. Locate the place on the planet that the analysis should be run. Create the boundary using a tool such as the [gmap pedometer](http://www.mappedometer.com/) and save the boundary points. 

2. Paste the coordinates into the script [createSphScript.py](https://github.com/airvette/agro_fusion/blob/master/createSphScript.py) and run. This creates a shapefile that will be needed to perform the analysis.

3. Download the Landsat data for the date of interest. [EarthExplorer](http://earthexplorer.usgs.gov/) is a good place to start, but some personal exploration will probably be needed to find the desired data. 

4. Run [analysisScript.py](https://github.com/airvette/agro_fusion/blob/master/analysisScript.py). Make sure analysisScript.py points to the files you want to do the NDVI analysis. Note that Band 4 (visual red) and Band 5 (near infrared) are the required bands for NDVI analysis. 
