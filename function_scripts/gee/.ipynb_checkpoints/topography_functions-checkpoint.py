# Begin here
import ee

# Authenticate to the Earth Engine servers
ee.Authenticate(auth_mode="gcloud")

# Initialize the Earth Engine module.
ee.Initialize(project= 'rmrs-wildfire-treatments')

import geemap
import geopandas as gpd
import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import os
import math
from datetime import datetime

def topography_variables(fireID, gee_data_source):
    
    # Load in fire perimeter -- this would be the area to add a buffer
    mtbs_perims = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1')
    
    # Load in DEM 
    nedDEM = ee.Image('USGS/3DEP/10m')
    
    fire_perims = ee.Feature(mtbs_perims.filterMetadata('Event_ID', 'equals', fireID).first())
    fire_perims = fire_perims.buffer(3000)
    demClipped = nedDEM.clip(fire_perims.geometry())

    
    # Calculate topography statisitcs
    elevation = demClipped.select('elevation')
    aspect = ee.Terrain.aspect(elevation)
    northness = aspect.multiply(math.pi/180).cos().rename('northness')
    eastness = aspect.multiply(math.pi/180).sin().rename('eastness')
    slope = ee.Terrain.slope(elevation)
    
    # DEM meta data
    CRS = nedDEM.projection().crs();
    elevationScale = nedDEM.projection().nominalScale();
    
    
    #Calculate heat load index, based on Theobald et al (2015)
    aspectRadians = aspect.multiply(ee.Number(0.01745329252))
    slopeRadians = slope.multiply(ee.Number(0.01745329252))
    foldedAspect = aspectRadians.subtract(ee.Number(4.3196899)).abs().multiply(ee.Number(-1)).add(ee.Number(3.141593)).abs()
    theobaldHLI = (slopeRadians.cos().multiply(ee.Number(1.582)).multiply(ee.Number(0.82886954044))).subtract(foldedAspect.cos().multiply(ee.Number(1.5)).multiply(slopeRadians.sin().multiply(ee.Number(0.55944194062)))).subtract(slopeRadians.sin().multiply(ee.Number(0.262)).multiply(ee.Number(0.55944194062))).add(foldedAspect.sin().multiply(ee.Number(0.607)).multiply(slopeRadians.sin())).add(ee.Number(-1.467));
    theobaldHLI = theobaldHLI.exp().rename('hli')
    
    # Functions are from : https://github.com/zachlevitt/earth-engine/blob/master/topoModules.md#processelevationdata
    # Reduces resolution of image 
    def processElevationData(image, crs, scale): 
        output = image.reduceResolution(
            reducer = ee.Reducer.mean(), 
            maxPixels =  1024).reproject(
            crs = crs,
            scale =  scale)
        return output
    
    # Calculate a neighborhood mean using a pre-defined kernel radius
    def calculateNeighborhoodMean(dem, kernelRadius):
        output = dem.reduceNeighborhood(
            reducer = ee.Reducer.mean(),
            kernel = ee.Kernel.square(kernelRadius,'pixels',False),
            optimization =  'boxcar')
        return output
    
    # Calculate a neighborhood standard deviation using a pre-defined kernel radius    
    def calculateNeighborhoodStdDev (dem, kernelRadius): 
      return dem.reduceNeighborhood(
          reducer =  ee.Reducer.stdDev(),
          kernel = ee.Kernel.square(kernelRadius,'pixels',False),
          optimization = 'window')
    
    # Calculate standardized topographic position index, based on Theobald et al (2015)
    def calculateStandardizedTPI(image, meanImage, stdDevImage):
      return (image.subtract(meanImage)).divide(stdDevImage)
    
    # Calculate mean topographic position index, based on Theobald et al (2015)
    def calculateMeanTPI(image1,image2,image3):
      return (image1.add(image2).add(image3)).divide(ee.Number(3)).rename('meanTPI')
    
    demMean_270m = calculateNeighborhoodMean(demClipped,27)
    demStdDev_270m = calculateNeighborhoodStdDev(demClipped,27)
    stdTPI_270m = calculateStandardizedTPI(demClipped,demMean_270m,demStdDev_270m)
      
    resampledDEM_810m = processElevationData(demClipped,CRS,elevationScale.multiply(3)) # Should I use this?
    demMean_810m = calculateNeighborhoodMean(demClipped,81)
    demStdDev_810m = calculateNeighborhoodStdDev(demClipped,81)
    stdTPI_810m = calculateStandardizedTPI(demClipped,demMean_810m,demStdDev_810m)
          
    resampledDEM_2430m = processElevationData(demClipped,CRS,elevationScale.multiply(9))
    demMean_2430m = calculateNeighborhoodMean(demClipped,243)
    demStdDev_2430m = calculateNeighborhoodStdDev(demClipped,243)
    stdTPI_2430m = calculateStandardizedTPI(demClipped,demMean_2430m,demStdDev_2430m)
            
    meanTPI = calculateMeanTPI(stdTPI_270m,stdTPI_810m,stdTPI_2430m)
    
    
    # Import other topographic variables -- only available with GEE academic accounts
    # Chili and mTPI are only available for no comercial partners
    # chili = ee.Image("CSP/ERGo/1_0/Global/SRTM_CHILI")
    # mtpi = ee.Image("CSP/ERGo/1_0/Global/SRTM_mTPI")
    
    topography = ee.Image([elevation, aspect, northness, eastness, slope, theobaldHLI, meanTPI])#.clip(fire_perims)

    return topography