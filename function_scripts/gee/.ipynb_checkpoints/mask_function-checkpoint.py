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
import sys
# !pip install nbimporter
#!pip install math
import math
#!pip install scikit-learn
#!pip install sklearn.preprocessing
import sklearn
from sklearn.preprocessing import StandardScaler

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

def mask_function(image, local_directory, fireID, fire_pd, var_pd, wto_outcome):
    fire_year =  int(str(fire_pd[fire_pd['MTBS_ID'] == fireID]['Year'].iloc[0]))

    # Mask variable dataframe
    var_masks = var_pd[var_pd['Group'] == 'Masks']

    if (var_masks[var_masks['Abbreviation'] == 'forest_mask'][wto_outcome].item()):
        
        # Calculate the pre-fire year
        pre_year = fire_year - 1

        # Load the LCMS land cover dataset select the pre-fire year
        # Define forest classes Trees: 1
        lcms_col = ee.ImageCollection('USFS/GTAC/LCMS/v2024-10')
        lcms = lcms_col.filterDate(str(pre_year), str(fire_year)).filter('study_area == "CONUS"').first().select('Land_Cover')
        forest_mask = lcms.eq(1).And(lcms.lte(5))                              
        
        # # Load the NLCD land cover dataset (e.g., 2016 version). 
        # nlcd = ee.Image('USGS/NLCD_RELEASES/2016_REL/2016').select('landcover') 
        # # Create a forest mask: Define forest classes (41: Deciduous Forest, 42: Evergreen Forest, 43: Mixed Forest).  
        # forest_mask = nlcd.eq(41).Or(nlcd.eq(42)).Or(nlcd.eq(43)) 
        
        # Apply forest mask (include only forest areas). 
        image = image.updateMask(forest_mask)
             
    if (var_masks[var_masks['Abbreviation'] == 'urban_mask'][wto_outcome].item()):
        # Load the pre-fire NLCD land cover dataset
        nlcd = get_nlcd_layer(fire_year)
        
        # Create a forest mask: 21-24 Urban classes and 31 is barren land.  
        developed_classes = [21, 22, 23, 24, 31]
        developed_img = nlcd.remap(developed_classes, [1] * len(developed_classes), 0)
        
        developed_mask = developed_img.eq(0)
        # Apply forest mask (include only forest areas). 
        image = image.updateMask(developed_mask) 
             
    if (var_masks[var_masks['Abbreviation'] == 'water_mask'][wto_outcome].item()):
        # Hansen water mask
        water_mask = ee.Image('UMD/hansen/global_forest_change_2024_v1_12').select(['datamask']).eq(1)
        # Apply water mask (mask out water). 
        image = image.updateMask(water_mask) 
    
    if (var_masks[var_masks['Abbreviation'] == 'landown_mask'][wto_outcome].item()):
        # Load ownership shape to mask out non-federal land -- eventually upload as an asset across WUS
        ownership_fp = local_directory + '/masks/land_ownership/shapes/' + fireID +'_land_own.gpkg'
        gdf = gpd.read_file(ownership_fp)
        ownership_features = geemap.geopandas_to_ee(gdf)
        # Filter to USFS, NPS, and BLM land 
        fed_ownership = ownership_features.filter(ee.Filter.inList('ADMIN_AGENCY_CODE', ee.List(['USFS', 'USBR','FWS' ,'NPS', 'BLM'])))
        image = image.clip(fed_ownership)
             
    return image

# function that exports the masks as tifs
def mask_function_export_images(image, local_directory, fireID, fire_pd, var_pd, wto_outcome):
    fire_year =  int(str(fire_pd[fire_pd['MTBS_ID'] == fireID]['Year'].iloc[0]))

    # Load fire geometry from MTBS for clipping
    mtbs = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1')
    fire_geom = mtbs.filter(ee.Filter.eq('Event_ID', fireID)).geometry()

    # Mask variable dataframe
    var_masks = var_pd[var_pd['Group'] == 'Masks']

    if (var_masks[var_masks['Abbreviation'] == 'forest_mask'][wto_outcome].item()):
        
        # Calculate the pre-fire year
        pre_year = fire_year - 1

        # Load the LCMS land cover dataset select the pre-fire year
        # Define forest classes Trees: 1
        lcms_col = ee.ImageCollection('USFS/GTAC/LCMS/v2024-10')
        lcms = lcms_col.filterDate(str(pre_year), str(fire_year)).filter('study_area == "CONUS"').first().select('Land_Cover')
        forest_mask = lcms.eq(1).clip(fire_geom)
        
        # Export forest mask as TIFF
        task_forest = ee.batch.Export.image.toCloudStorage(
            image=forest_mask,
            description=f'{fireID}_forest_mask',
            bucket='bb-gee-bucket',
            fileNamePrefix=f'mogollon/masks/{fireID}_forest_mask',
            scale=30,
            region=fire_geom,
            fileFormat='GeoTIFF',
            maxPixels=1e9
        )
        task_forest.start()
        
        # Apply forest mask (include only forest areas). 
        image = image.updateMask(forest_mask)
             
    if (var_masks[var_masks['Abbreviation'] == 'urban_mask'][wto_outcome].item()):
        # Load the pre-fire NLCD land cover dataset
        nlcd = get_nlcd_layer(fire_year)
        
        # Create urban mask: 21-24 Urban classes and 31 is barren land.  
        developed_classes = [21, 22, 23, 24, 31]
        developed_img = nlcd.remap(developed_classes, [1] * len(developed_classes), 0)
        
        urban_mask = developed_img.eq(0).clip(fire_geom)
        
        # Export urban mask as TIFF
        task_urban = ee.batch.Export.image.toCloudStorage(
            image=urban_mask,
            description=f'{fireID}_urban_mask',
            bucket='bb-gee-bucket',
            fileNamePrefix=f'mogollon/masks/{fireID}_urban_mask',
            scale=30,
            region=fire_geom,
            fileFormat='GeoTIFF',
            maxPixels=1e9
        )
        task_urban.start()
        
        # Apply urban mask (include only non-urban areas). 
        image = image.updateMask(urban_mask) 
             
    if (var_masks[var_masks['Abbreviation'] == 'water_mask'][wto_outcome].item()):
        # Hansen water mask
        water_mask = ee.Image('UMD/hansen/global_forest_change_2024_v1_12').select(['datamask']).eq(1).clip(fire_geom)
        
        # # Export water mask as TIFF
        # task_water = ee.batch.Export.image.toCloudStorage(
        #     image=water_mask,
        #     description=f'{fireID}_water_mask',
        #     bucket='bb-gee-bucket',
        #     fileNamePrefix=f'klamath/masks/{fireID}_water_mask',
        #     scale=30,
        #     region=fire_geom,
        #     fileFormat='GeoTIFF',
        #     maxPixels=1e9
        # )
        # task_water.start()
        
        # Apply water mask (mask out water). 
        image = image.updateMask(water_mask) 
    
    if (var_masks[var_masks['Abbreviation'] == 'landown_mask'][wto_outcome].item()):
        # Load ownership shape to mask out non-federal land -- eventually upload as an asset across WUS
        ownership_fp = local_directory + '/masks/land_ownership/shapes/' + fireID +'_land_own.gpkg'
        gdf = gpd.read_file(ownership_fp)
        ownership_features = geemap.geopandas_to_ee(gdf)
        # Filter to USFS, NPS, and BLM land 
        fed_ownership = ownership_features.filter(ee.Filter.inList('ADMIN_AGENCY_CODE', ee.List(['USFS', 'USBR','FWS' ,'NPS', 'BLM'])))
        image = image.clip(fed_ownership)
             
    return image