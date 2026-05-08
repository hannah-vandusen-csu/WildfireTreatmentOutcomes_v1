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
# !pip install nbimporter
import nbimporter


# Function to apply scaling factors
def apply_scale_factors(image):
    optical_bands = image.select(['SR_B.']).multiply(0.0000275).add(-0.2)
    return image.addBands(optical_bands, None, True)
# Function to apply old scaling factors
def apply_scale_factors2(image):
    optical_bands = image.select(['SR_B.']).multiply(0.001)#.add(-0.2)
    return image.addBands(optical_bands, None, True)

# Function to compute indices for Landsat 8 and 9
def calculate_indices_ls8_9(image):
    nbr = image.normalizedDifference(['SR_B5', 'SR_B7']).toFloat()
    ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).toFloat()
    ndmi = image.normalizedDifference(['SR_B5', 'SR_B6']).toFloat()

    evi = image.expression(
              '(2.5 * ((B5 - B4) / (B5 + (6 * B4) - (7.5 * B2) + 1)))',
              {'B5': image.select('SR_B5'),
              'B4': image.select('SR_B4'),
              'B2': image.select('SR_B2'),
              }).toFloat()
    mirbi = image.expression(
        '((10 * SR_B6) - (9.8 * SR_B7) + 2)',
        {'SR_B6': image.select('SR_B6'),
         'SR_B7': image.select('SR_B7')}
    ).toFloat()
    qa = image.select(['QA_PIXEL'])
    return nbr.addBands([ndvi, ndmi, evi, mirbi, qa]).select([0, 1, 2, 3, 4, 5], ['nbr', 'ndvi', 'ndmi', 'evi', 'mirbi', 'QA_PIXEL']).copyProperties(image, ['system:time_start'])

# Function to compute indices for Landsat 4, 5, and 7
def calculate_indices_ls4_7(image):
    nbr = image.normalizedDifference(['SR_B4', 'SR_B7']).toFloat()
    ndvi = image.normalizedDifference(['SR_B4', 'SR_B3']).toFloat()
    ndmi = image.normalizedDifference(['SR_B4', 'SR_B5']).toFloat()
    evi = image.expression(
        '2.5 * ((SR_B4 - SR_B3) / (SR_B4 + 6 * SR_B3 - 7.5 * SR_B1 + 1))',
        {'SR_B4': image.select('SR_B4'),
         'SR_B3': image.select('SR_B3'),
         'SR_B1': image.select('SR_B1')}
    ).toFloat()
    mirbi = image.expression(
        '((10 * SR_B5) - (9.8 * SR_B7) + 2)',
        {'SR_B5': image.select('SR_B5'),
         'SR_B7': image.select('SR_B7')}
    ).toFloat()
    qa = image.select(['QA_PIXEL'])
    return nbr.addBands([ndvi, ndmi, evi, mirbi, qa]).select([0, 1, 2, 3, 4, 5], ['nbr', 'ndvi', 'ndmi', 'evi', 'mirbi', 'QA_PIXEL']).copyProperties(image, ['system:time_start'])


# FUNCTION TO MASK CLOUD, WATER, SNOW, ETC.
def lsCfmask(image):
    # Bits 3,4,5,7: cloud,cloud-shadow,snow,water respectively.
    clouds_bit_mask = (1 << 3)
    cloud_shadow_bit_mask = (1 << 4)
    snow_bit_mask = (1 << 5)
    water_bit_mask = (1 << 7)
    # Get the pixel QA band.
    qa = image.select('QA_PIXEL')
    
    # Flags should be set to zero, indicating clear conditions.
    clear = qa.bitwiseAnd(clouds_bit_mask).eq(0).And(qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0)).And(qa.bitwiseAnd(snow_bit_mask).eq(0)).And(qa.bitwiseAnd(water_bit_mask).eq(0))
    return image.updateMask(clear).select([0, 1, 2, 3, 4]).copyProperties(image, ["system:time_start"])


# Main Function for per fire processing
def process_fire(fire_id, fire_pd, local_directory, var_crosswalks, start_day, end_day, mtbs_perim = False ): # true/false for using mtbs or local shapefile
    # 0) FIRE ID ---------------------------------------------------------------------------------
    fire = fire_id
    
    # 1) INITIALIZE FIRE DATAFRAME FOR PROCESSING -------------------------------------------------
    # Add years
    fire_pd['Year'] = pd.to_datetime(fire_pd['Ignition']).dt.year
    
    # Date Manipulation of Dataframe
    # Define fire season for AOI
    # Parks et al 2014 used startday = 91 and endday = 181 for fires in Arizona, New Mexico, and Utah.
    # And startday = 152 and endday = 273 for fires in California, Montana, Washington, and Wyoming.
    
    fire_pd['start_day'] = start_day
    fire_pd['end_day'] = end_day

    # All available bands for export
    col_var_df = var_crosswalks[var_crosswalks['file_path'] == 'parks 2018']
    bandList = col_var_df['Abbreviation'].tolist()
    #bandList = ['dnbr', 'rbr', 'rdnbr', 'dndvi', 'devi', 'dndmi', 'dmirbi', 'post_nbr', 'post_mirbi', 'CBI', 'CBI_bc']

    # 2) RANDOM FOREST SPECIFICATIONS------------------------------------------------------------
    # Bands for the random forest classification
    rf_bands = ['def', 'lat', 'rbr', 'dmirbi', 'dndvi', 'post_mirbi']
    
    # Load training data for Random Forest classification
    cbi = ee.FeatureCollection("users/grumpyunclesean/CBI_predictions/data_for_ee_model")
    
    # Load climatic water deficit variable (def) for random forest classification
    def_img = ee.Image("users/grumpyunclesean/CBI_predictions/def").rename('def').int()
    
    # Create latitude image for random forest classification
    lat = ee.Image.pixelLonLat().select('latitude').rename('lat').round().toInt()
    
    # Parameters for random forest classification
    nrow_training_fold = cbi.size()  # number of training observations
    minLeafPopulation = nrow_training_fold.divide(75).divide(6).round()
    
    # Random forest classifier
    fsev_classifier = ee.Classifier.smileRandomForest(
        numberOfTrees = 500,
        minLeafPopulation = minLeafPopulation,
        seed = 123).train(cbi,'CBI',rf_bands).setOutputMode('REGRESSION')

    #  3) GET LANDSAT COLLECTIONS --------------------------------     
    # Landsat 5, 7, 8 and 9 Surface Reflectance (Level 2) Tier 1 Collection 2 
    
    # Get Landsat collections
    ls_collections = {
        'ls9': ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'),
        'ls8': ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'),
        'ls7': ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'),
        'ls5': ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'),
        'ls4': ee.ImageCollection('LANDSAT/LT04/C02/T1_L2')
    }

    # Create water mask from Hansen's Global Forest Change to use in processing function
    water_mask = ee.Image('UMD/hansen/global_forest_change_2024_v1_12').select(['datamask']).eq(1)
    
    # Apply the functions to the collections
    ls_collections = {key: col.map(apply_scale_factors)
                      .map(calculate_indices_ls8_9 if key in ['ls8', 'ls9'] else calculate_indices_ls4_7)
                      .map(lsCfmask) for key, col in ls_collections.items()}
    
    # Merge Landsat Collections 
    ls_col = ee.ImageCollection(ls_collections['ls9']
                                .merge(ls_collections['ls8'])
                                .merge(ls_collections['ls7'])
                                .merge(ls_collections['ls5'])
                                .merge(ls_collections['ls4']))

    # 4) CALCULATE BURN SEVERITY ON A PER FIRE BASIS -----------------------------------------  
    if mtbs_perim == True:
    # indent lines below  
        mtbs_perims = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1')
        fire_bounds = ee.Feature(mtbs_perims.filterMetadata('Event_ID', 'equals', fire).first())
        fire_bounds = fire_bounds.geometry().buffer(1000)
    else:
        merged_perim_gdf = gpd.read_file(local_directory +'/fire_perimeters/merged/' + fire + '_Extent.shp')
        fire_bounds = geemap.geopandas_to_ee(merged_perim_gdf)
        fire_bounds = fire_bounds.first().geometry().buffer(1000)
    
    fire_year = ee.Date.parse('YYYY', str(fire_pd[fire_pd['MTBS_ID'] == fire]['Year'].iloc[0]))
    start_day = ee.Number.parse(str(fire_pd[fire_pd['MTBS_ID'] == fire]['start_day'].iloc[0]))
    end_day = ee.Number.parse(str(fire_pd[fire_pd['MTBS_ID'] == fire]['end_day'].iloc[0]))


    # Pre-Imagery
    pre_fire_year = fire_year.advance(-1, 'year')
    # Check if imagery is available for 1 year pre fire; if so, get means across pixels; otherwise, output imagery will be masked
    pre_fire_indices = ee.Algorithms.If(
        ls_col.filterBounds(fire_bounds).filterDate(pre_fire_year, fire_year).filter(ee.Filter.dayOfYear(start_day, end_day)).size(),
        ls_col.filterBounds(fire_bounds).filterDate(pre_fire_year, fire_year).filter(ee.Filter.dayOfYear(start_day, end_day)).mean().select([0,1,2,3,4], ['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi']),
        ee.Image.cat(ee.Image(), ee.Image(), ee.Image(), ee.Image(), ee.Image()).rename(['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi'])
    )
    
    # If any pixels within fire have only one 'scene' or less, add additional year backward to fill in behind
    pre_fire_year2 = fire_year.advance(-2, 'year')
    pre_fire_indices2 = ee.Algorithms.If(
        ls_col.filterBounds(fire_bounds).filterDate(pre_fire_year2, fire_year).filter(ee.Filter.dayOfYear(start_day, end_day)).size(),
        ls_col.filterBounds(fire_bounds).filterDate(pre_fire_year2, fire_year).filter(ee.Filter.dayOfYear(start_day, end_day)).mean().select([0,1,2,3,4], ['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi']),
        ee.Image.cat(ee.Image(), ee.Image(), ee.Image(), ee.Image(), ee.Image()).rename(['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi'])
    )
    
    pre_filled = ee.Image(pre_fire_indices).unmask(pre_fire_indices2)
    
    # Post-Imagery
    post_fire_year = fire_year.advance(1, 'year')
    post_fire_year2 = fire_year.advance(2, 'year')
    # Check if imagery is available for 1 year post fire; if so, get means across pixels; otherwise, output imagery will be masked
    post_fire_indices = ee.Algorithms.If(
        ls_col.filterBounds(fire_bounds).filterDate(post_fire_year, post_fire_year2).filter(ee.Filter.dayOfYear(start_day, end_day)).size(),
        ls_col.filterBounds(fire_bounds).filterDate(post_fire_year, post_fire_year2).filter(ee.Filter.dayOfYear(start_day, end_day)).mean().select([0,1,2,3,4], ['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi']),
        ee.Image.cat(ee.Image(), ee.Image(), ee.Image(), ee.Image(), ee.Image()).rename(['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi'])
    )
    
    post_fire_year3 = fire_year.advance(3, 'year')
    # If any pixels within fire have only one 'scene' or less, add additional year forward to fill in behind
    post_fire_indices2 = ee.Algorithms.If(
        ls_col.filterBounds(fire_bounds).filterDate(post_fire_year, post_fire_year3).filter(ee.Filter.dayOfYear(start_day, end_day)).size(),
        ls_col.filterBounds(fire_bounds).filterDate(post_fire_year, post_fire_year3).filter(ee.Filter.dayOfYear(start_day, end_day)).mean().select([0,1,2,3,4], ['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi']),
        ee.Image.cat(ee.Image(), ee.Image(), ee.Image(), ee.Image(), ee.Image()).rename(['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi'])
    )
    
    post_filled = ee.Image(post_fire_indices).unmask(post_fire_indices2)
    fireIndices = pre_filled.addBands(post_filled) #Here
    
    # Calculate dNBR
    burnIndices = fireIndices.expression(
        "(b('pre_nbr') - b('post_nbr')) * 1000"
    ).rename('dnbr').toInt().addBands(fireIndices)
    
    # Calculate RBR
    burnIndices2 = burnIndices.expression(
        "b('dnbr') / (b('pre_nbr') + 1.001)"
    ).rename('rbr').toInt().addBands(burnIndices)
    
    # Calculate RdNBR
    burnIndices3 = burnIndices2.expression(
        "abs(b('pre_nbr')) < 0.001 ? 0.001 : b('pre_nbr')"
    ).abs().sqrt().rename('pre_nbr2').toFloat().addBands(burnIndices2)
   
    burnIndices4 = burnIndices3.expression(
        "b('dnbr') / b('pre_nbr2')"
    ).rename('rdnbr').toInt().addBands(burnIndices3)
    
    # Calculate basal area percent change (Miller 2009a)
    burnIndices45 = burnIndices4.expression(
        "100 * pow(sin((b('rdnbr') - 166.5) / 389), 2)"
    ).rename('per_BA').toInt().addBands(burnIndices4)
    
    # Calculate dNDVI
    burnIndices5 = burnIndices45.expression(
        "(b('pre_ndvi') - b('post_ndvi')) * 1000"
    ).rename('dndvi').toInt().addBands(burnIndices45)
    
    # Calculate dEVI
    burnIndices6 = burnIndices5.expression(
        "(b('pre_evi') - b('post_evi')) * 1000"
    ).rename('devi').toInt().addBands(burnIndices5)
    
    # Calculate dNDMI
    burnIndices7 = burnIndices6.expression(
        "(b('pre_ndmi') - b('post_ndmi')) * 1000"
    ).rename('dndmi').toInt().addBands(burnIndices6)
    
    # Calculate dMIRBI
    burnIndices8 = burnIndices7.expression(
        "(b('pre_mirbi') - b('post_mirbi')) * 1000"
    ).rename('dmirbi').toInt().addBands(burnIndices7)
    
    
    # Multiply post_mirbi band by 1000 to put it on the same scale as CBI plot extractions
    post_mirbi_1000 = burnIndices8.select("post_mirbi").multiply(1000).toInt()
    burnIndices8 = burnIndices8.addBands(post_mirbi_1000, None, True)  # None to copy all bands; True to overwrite original post_mirbi
    
    # Get a list of all metadata properties.
    properties = burnIndices8.propertyNames()
    
    # Add in climatic water deficit variable, i.e. def
    burnIndices9 = burnIndices8.addBands(def_img)
    
    # Add in latitude
    burnIndices10 = burnIndices9.addBands(lat)
    
    # Classify the image with the same bands used to train the Random Forest classifier
    cbi_rf = burnIndices10.select(rf_bands).classify(fsev_classifier).rename('CBI').toFloat().multiply(10**2).floor().divide(10**2)
    
    burnIndices11 = cbi_rf.addBands(burnIndices10)
    
    # Create bias corrected CBI
    def bias_correct(bandName):
        cbi_lo = bandName.expression("((b('CBI') - 1.5) * 1.3)  + 1.5")
        cbi_hi = bandName.expression("((b('CBI') - 1.5) * 1.175) + 1.5")
        cbi_mg = bandName.where(bandName.lte(1.5), cbi_lo).where(bandName.gt(1.5), cbi_hi)
        return cbi_mg.where(cbi_mg.lt(0), 0).where(cbi_mg.gt(3), 3).multiply(10**2).floor().divide(10**2).rename('CBI_bc')
    
    validMask = pre_filled.select('pre_nbr').add(10).add(post_filled.select('post_nbr')).add(10)  # Adding 10 ensures resulting image doesn't have 0 values that would become masked in final output bands
    mask = water_mask.updateMask(validMask)
    burnIndices12 = bias_correct(burnIndices11.select('CBI')).addBands(burnIndices11)
    burnIndices12 = burnIndices12.updateMask(mask).clip(fire_bounds).select(bandList, bandList)

    return burnIndices12.set({'fireID': fire, 'fireYear': fire_year})
    

# PROCESSING

# FUNCTION TO CREATE PREFIRE SPECTRAL IMAGERY FOR EACH FIRE (THIS IS JUST THE FIRST HALF OF THE FIRE SEV FUNCTION
def prefire_spectral_imagery(fireID, fire_pd, local_directory, var_crosswalks, start_day, end_day, mtbs_perim = False ):
    
    # Date Manipulation of Dataframe
    # Define fire season for AOI
    # Parks et al 2014 used startday = 91 and endday = 181 for fires in Arizona, New Mexico, and Utah.
    # And startday = 152 and endday = 273 for fires in California, Montana, Washington, and Wyoming.
    # For Klamath-Trinity, I defined startday = 152 and endday = 273, for other landscapes we can add this as a function input

    fire_pd['start_day'] = start_day
    fire_pd['end_day'] = end_day
    
    # bands
    col_var_df = var_crosswalks[var_crosswalks['file_path'] == 'parks 2018']
    bandList = col_var_df['Abbreviation'].tolist()

    # GET LANDSAT COLLECTIONS --------------------------------     
    # Landsat 5, 7, 8 and 9 Surface Reflectance (Level 2) Tier 1 Collection 2 
    
    # Get Landsat collections
    ls_collections = {
        'ls9': ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'),
        'ls8': ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'),
        'ls7': ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'),
        'ls5': ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'),
        'ls4': ee.ImageCollection('LANDSAT/LT04/C02/T1_L2')
    }

    # Create water mask from Hansen's Global Forest Change to use in processing function
    water_mask = ee.Image('UMD/hansen/global_forest_change_2024_v1_12').select(['datamask']).eq(1)
    
    # Apply the functions to the collections
    ls_collections = {key: col.map(apply_scale_factors)
                      .map(calculate_indices_ls8_9 if key in ['ls8', 'ls9'] else calculate_indices_ls4_7)
                      .map(lsCfmask) for key, col in ls_collections.items()}
    
    # Merge Landsat Collections 
    ls_col = ee.ImageCollection(ls_collections['ls9']
                                .merge(ls_collections['ls8'])
                                .merge(ls_collections['ls7'])
                                .merge(ls_collections['ls5'])
                                .merge(ls_collections['ls4']))

    # Get fire metadata 
    if mtbs_perim == True :
    # indent lines below  
        mtbs_perims = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1')
        fire_bounds = ee.Feature(mtbs_perims.filterMetadata('Event_ID', 'equals', fireID).first())
        fire_bounds = fire_bounds.buffer(1000)
    else:
        merged_perim_gdf = gpd.read_file(local_directory +'/fire_perimeters/merged/' + fireID + '_Extent.gpkg')
        fire_bounds = geemap.geopandas_to_ee(merged_perim_gdf)
        fire_bounds = fire_bounds.first().geometry().buffer(1000)
        
    fire_year = ee.Date.parse('YYYY', str(fire_pd[fire_pd['MTBS_ID'] == fireID]['Year'].iloc[0]))
    pre_fire_year = fire_year.advance(-1, 'year') # Pre-Imagery
    
    # Check if imagery is available for 1 year pre fire; if so, get means across pixels; otherwise, output imagery will be masked
    pre_fire_indices = ee.Algorithms.If(
        ls_col.filterBounds(fire_bounds).filterDate(pre_fire_year, fire_year).filter(ee.Filter.dayOfYear(start_day, end_day)).size(),
        ls_col.filterBounds(fire_bounds).filterDate(pre_fire_year, fire_year).filter(ee.Filter.dayOfYear(start_day, end_day)).mean().select([0,1,2,3,4], ['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi']),
        ee.Image.cat(ee.Image(), ee.Image(), ee.Image(), ee.Image(), ee.Image()).rename(['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi'])
    )
    
    # If any pixels within fire have only one 'scene' or less, add additional year backward to fill in behind
    pre_fire_year2 = fire_year.advance(-2, 'year')
    pre_fire_indices2 = ee.Algorithms.If(
        ls_col.filterBounds(fire_bounds).filterDate(pre_fire_year2, fire_year).filter(ee.Filter.dayOfYear(start_day, end_day)).size(),
        ls_col.filterBounds(fire_bounds).filterDate(pre_fire_year2, fire_year).filter(ee.Filter.dayOfYear(start_day, end_day)).mean()
        .select([0,1,2,3,4], ['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi']),
        ee.Image.cat(ee.Image(), ee.Image(), ee.Image(), ee.Image(), ee.Image())
        .rename(['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi'])
    )
    
    pre_filled = ee.Image(pre_fire_indices).unmask(pre_fire_indices2)
    

    # Mask out water
    validMask = pre_filled.select('pre_ndvi').add(10)
    mask = water_mask.updateMask(validMask)
    pre_filled = pre_filled.updateMask(mask).clip(fire_bounds).select(bandList, bandList)

    return pre_filled


# FUNCTION TO EXTRACT PREFIRE GEDI-BASED CANOPY METRICS
# def prefire_gedi_canopy(fireID, fire_pd, local_directory, var_crosswalks, mtbs_perim = False ):
def prefire_gedi_canopy(fireID, fire_pd, local_directory, var_crosswalks, mtbs_perim = False):
    
    
    fire_year = int(fire_pd[fire_pd['MTBS_ID'] == fireID]['Year'].iloc[0])  # Convert from pandas to integer
    pre_fire_year = fire_year - 1  # Pre-fire imagery year
    
    # GEDI canopy metric asset paths
    base_path = "projects/rmrs-wildfire-treatments/assets/GEDI_veg/"
    Cover_fp = f"{base_path}KT_Canopy_Cover_{pre_fire_year}"
    Height_fp = f"{base_path}KT_Canopy_Height_{pre_fire_year}"
    AGBD_fp = f"{base_path}KT_Canopy_AGBD_{pre_fire_year}"
    PAI_fp = f"{base_path}KT_Canopy_PAI_{pre_fire_year}"

    # Load GEDI images
    Cover = ee.Image(Cover_fp).rename("Cover")
    Height = ee.Image(Height_fp).rename("Height")
    AGBD = ee.Image(AGBD_fp).rename("AGBD")
    PAI = ee.Image(PAI_fp).rename("PAI")

    # Stack images together and clip to fire boundary
    gedi_stack = ee.Image.cat(Cover, Height, AGBD, PAI) #.clip(fire_bounds)
    
    # Load Hansen's Global Forest Change water mask (water = 0, land = 1)
    water_mask = ee.Image('UMD/hansen/global_forest_change_2024_v1_12').select('datamask').eq(1)

    # Apply water mask
    validMask = gedi_stack.select("Cover").add(10)  # Ensure valid values
    mask = water_mask.updateMask(validMask)
    gedi_stack = gedi_stack.updateMask(mask) #.clip(fire_bounds)

    return gedi_stack
