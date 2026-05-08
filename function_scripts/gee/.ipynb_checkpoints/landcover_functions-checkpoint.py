# FUNCTIONS TO CREATE SPECTRAL INDICES
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

# Ifelse statement to select appropriate NLCD layer
def get_nlcd_layer(fire_year):
    if fire_year < 2021:
        # Import the 2019 NLCD collection
        dataset19 = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')
    else:
        # Import the 2021 NLCD collection
        dataset21 = ee.ImageCollection('USGS/NLCD_RELEASES/2021_REL/NLCD')
    
    if fire_year <= 2005:
        # Filter the collection to the 2019 product
        nlcd2001 = dataset19.filter(ee.Filter.eq('system:index', '2001')).first().select('landcover')
        return nlcd2001
    elif fire_year <= 2010:
        # Filter the collection to the 2019 product
        nlcd2006 = dataset19.filter(ee.Filter.eq('system:index', '2006')).first().select('landcover')
        return nlcd2006
    elif fire_year <= 2015:
        # Filter the collection to the 2019 product
        nlcd2011 = dataset19.filter(ee.Filter.eq('system:index', '2011')).first().select('landcover')
        return nlcd2011
    elif fire_year <= 2018:
        # Filter the collection to the 2019 product
        nlcd2016 = dataset19.filter(ee.Filter.eq('system:index', '2016')).first().select('landcover')
        return nlcd2016
    elif fire_year <= 2020:
        # Filter the collection to the 2019 product
        nlcd2019 = dataset19.filter(ee.Filter.eq('system:index', '2019')).first().select('landcover')
        return nlcd2019
    else:
        # Filter the collection to the 2021 product
        nlcd2021 = dataset21.filter(ee.Filter.eq('system:index', '2021')).first().select('landcover')
        return nlcd2021

# Function to crop NLCD layer to fire perimeter and export
def export_nlcd_for_fire(fireID, fire_pd):
    # select the fire year
    year =  int(str(fire_pd[fire_pd['MTBS_ID'] == fireID]['Year'].iloc[0]))
    # Run function
    nlcd_layer = get_nlcd_layer(year)
    # Clip to fire
    return nlcd_layer

# Function to aquire EPA ecoregions 
def epa_ecoregions(fireID, eco_level):

    # Load EPA ecoregions 
    epa = ee.FeatureCollection("EPA/Ecoregions/2013/L4")
    
    # Select the correct epa ecoregion level
    if eco_level == 3:
        epa_level = epa.select('us_l3name')
        epa_name = 'us_l3name'
        epa_name_id = 'us_l3name_id'
    else: 
        epa_level = epa.select('us_l4name')
        epa_name = 'us_l4name'
        epa_name_id = 'us_l4name_id'
    
    final_nm = 'epa_' + str(eco_level)
        
    # Convert EPA to numeric for rasterization and ee Dictionary
    names = epa.aggregate_array(epa_name).distinct().getInfo()
    name_dict = {name: i+1 for i, name in enumerate(names)}
    name_dict_ee = ee.Dictionary(name_dict)
    
    # add numeric ID
    def add_numeric_id(feature):
        return feature.set(epa_name_id, name_dict_ee.get(feature.get(epa_name)))
    epa_with_id = epa.map(add_numeric_id)
    
    # Rasterize
    epaImg = epa_with_id.reduceToImage(
        properties = [epa_name_id],
        reducer = ee.Reducer.first()).rename([final_nm])
    
    return epaImg 

# Function to gather FRS layer 
def frs(fireID):
    frs_img = ee.Image('projects/rmrs-wildfire-treatments/assets/frs_spatial')
    return(frs_img)

def landfire(code):
    landfire_layer = ee.Image('LANDFIRE/Vegetation/BPS/v1_4_0/CONUS').select(code)
    return(landfire_layer)

def lcms_maxyears(fireID, fire_pd):
     # select the fire year
    fire_year =  int(str(fire_pd[fire_pd['MTBS_ID'] == fireID]['Year'].iloc[0]))
    pre_fire_year = fire_year - 1
    dataset = ee.ImageCollection('USFS/GTAC/LCMS/v2023-9');
    
    lcms = dataset.filterDate('1985', str(pre_fire_year)).filter('study_area == "CONUS"') 
    
    # Function for getting a more general integer year date of change
    def getYrMsk(img):
      yr = ee.Date(img.get('system:time_start')).get('year')
      return ee.Image(yr).int16().updateMask(img.mask())
        
    # Pull apart LCMS fast and slow loss and find the most recent year of each 
    lcms_fast_loss_yr = lcms.select(['Change']).map(lambda img:getYrMsk(img.updateMask(img.eq(3)))).max().rename(['lcms_fast_loss_yr'])
    lcms_slow_loss_yr = lcms.select(['Change']).map(lambda img:getYrMsk(img.updateMask(img.eq(2)))).max().rename(['lcms_slow_loss_yr'])
    lcms_gain_yr = lcms.select(['Change']).map(lambda img:getYrMsk(img.updateMask(img.eq(4)))).max().rename(['lcms_gain_yr'])
    
    return ee.Image([lcms_fast_loss_yr, lcms_slow_loss_yr, lcms_gain_yr])

def mtbs_severity(fireID, fire_pd, gee_data_source):
    
    fire_year =  int(str(fire_pd[fire_pd['MTBS_ID'] == fireID]['Year'].iloc[0]))
    mtbs_severity_str = gee_data_source +'/mtbs_CONUS_' + str(fire_year)
    mtbs_severity = ee.Image(mtbs_severity_str)
    return(mtbs_severity)
