# Functions to calculate climate indicies and climate normals 
import ee
ee.Authenticate(auth_mode="gcloud")
ee.Initialize(project= 'rmrs-wildfire-treatments')
import geemap
import geopandas as gpd
import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import os
from datetime import datetime


def extract_prism_climate_reducer(fireID, fire_list, reducer):
    # 2. Annual PRISM climate Normals 
    # https://developers.google.com/earth-engine/datasets/catalog/OREGONSTATE_PRISM_Norm91m#description
    # must take mean for annual 30-year Normal
    prismNorm = ee.ImageCollection("OREGONSTATE/PRISM/Norm91m")
    
    # Prism Cliamte Normal Bands
    prism_bands = ['ppt','tmean', 'tmin', 'vpdmax']

    # Load in fire perimeter
    mtbs_perims = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1')
    fire_perims = ee.Feature(mtbs_perims.filterMetadata('Event_ID', 'equals', fireID).first())
   
    # annual mean across bands for all 12 months
    prism_fire = prismNorm.select(prism_bands).reduce(reducer)#.clip(fire_perims)
    # return prism climate normal values 
    return ee.Image(prism_fire)

def climate_data_91(fireID,
                    gee_data_source, 
                    var_crosswalks,  # variable crosswalks sheet
                    reducer = ee.Reducer.mean()): # mean for 30-years
  
    climate_collection = ee.ImageCollection(gee_data_source)
    col_var_df = var_crosswalks[var_crosswalks['file_path'] == gee_data_source]
    bands = col_var_df['Abbreviation'].tolist()
    # Take the 30 - year mean
    climate_mean = climate_collection.filter(ee.Filter.date('1991-01-01', '2020-12-31')).select(bands).mean()#.clip(fire_perims)
    # return prism climate normal values 
    return ee.Image(climate_mean)   
    
def extract_terra_climate_reducer(fireID, reducer):
  
    terraClim = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
    
    # Prism Cliamte Normal Bands
    bands = ['def','swe']

    # Load in fire perimeter
    mtbs_perims = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1')
    fire_perims = ee.Feature(mtbs_perims.filterMetadata('Event_ID', 'equals', fireID).first())
   
    # Take the 30 - year mean
    terraClim = terraClim.filter(ee.Filter.date('1990-01-01', '2020-12-31')).select(bands).reduce(reducer)#.clip(fire_perims)
    # return prism climate normal values 
    return ee.Image(terraClim)   