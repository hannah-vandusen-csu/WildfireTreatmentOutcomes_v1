# Begin here
import ee

# Authenticate to the Earth Engine servers
# ee.Authenticate(auth_mode="gcloud")

# Initialize the Earth Engine module.
# ee.Initialize(project= 'rmrs-wildfire-treatments')

import geemap
import geopandas as gpd
import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import os

from datetime import datetime

def new_snow_metrics(fireID, fire_pd):
        
    # Gather fire metadata
    
    wyr =  fire_pd[fire_pd['MTBS_ID'] == fireID]['Year'].iloc[0]
    
    # Load in fire perimeter
    mtbs_perims = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1')
    fire_perims = ee.Feature(mtbs_perims.filterMetadata('Event_ID', 'equals', fireID).first())
    
    # Set the month and day of year for the Water Year to begin
    wyStartDay = '-10-01'
    
    NDSI_threshold = '15'
    
    # Set the lower end of the SCF values for masking
    # This is OPTIONAL. Setting this value to 0 will not mask any values.
    # This value equals about 30 days of snow or less that is being masked. 
    # This is set because without it there are 'snow' designations all across the southern U.S., Mexico, Central America, Africa, etc.
    # My guess is that agricultural water bodies are being mistaken for snow, rice growing, flooded fields, and etc.
    low_value = 0.07
    
    # Set the prior_days parameter for the SDD calculation.
    # This gets fed into the 'min_sdays_prior_to_SDD' image and iable.
    # Setting this value means there are [prior_days] days that snow must be seen/designated in the pixel before the SDD date can be set.
    prior_days = 5
    
    # Set the parameter below for the SDD calculation.
    # This gets fed into the 'minNoSnow' image and iable.
    # Setting this value means there are [minNoSnowdays] days that snow must be absent from the pixel after the SDD date is set.
    minNoSnowdays = 5
    
    # Set the projection for the geotiff export
    projection = 'EPSG:4326'
    
    # These are the min and max corners of the domain for export, in Lat, Long
    # Input the minimum lat, lower left corner
    minLat = 28.0184
    # Input the minimum long, lower left corner
    minLong = -126.5330
    # Input the max lat, upper right corner
    maxLat = 50.1380
    # Input the max Long, upper right corner
    maxLong = -94.6287
    
    # My domain below is the western US and southern Canadian Rockies region
    # my_domain =  {'color': '0b4a8b': ee.Geometry.Polygon(
    #     [[[-126.533, 50.138],
    #       [-126.533, 28.0184],
    #       [-94.6287, 28.0184],
    #       [-94.6287, 50.138]]])}
    
    # Uncomment the lines below to add the my_domain iable to the visualization
    # Map.addLayer(my_domain,'','My Domain')
    # print (my_domain)
    
    # These are the color palettes for the visualization of the metrics below.
    # palette_snow = ['black', '0dffff', '0524ff', 'ffffff']
    # palette_sdd = 'bf812d,dfc27d,f6e8c3,f5f5f5,c7eae5,80cdc1,35978f,01665e'
    
    # Load Function to run for each year
    
    # Take the end of the water year and create a start and end date.
    start = str(wyr - 1) + wyStartDay
    #print(start,'start')
    end = str(wyr) + wyStartDay
    
    # Take the end date of the water year and increase it by 30 days.
    start_adv = end
    end_adv = ee.Date(end).advance(30,'day')
    #print(start_adv, 'start_adv')
    #display(end_adv,'end_adv')
        
    # Advance the days of interest by 30 days into the future
    adv_days = ee.ImageCollection("MODIS/061/MOD10A1").filterDate(start_adv,end_adv).size().getInfo()
    #display(adv_days,'adv_days')
    
    ########## ACCESS the DATA and do some things  ##########
    #############  SATELLITE DATASETS  ######################
    #########################################################
    
    # This makes a list of images where the MOD10A1v6 designation equals cloud(250).
    # See printout to console for more info.
    # This creates a yes/no (0/1) image of cloud locations each day of the year.
    MODIS_cloud_list = ee.ImageCollection("MODIS/061/MOD10A1") \
        .filterDate(start,end_adv) \
        .map(lambda img: 
             img.expression(
                 '(BAND == 250)',{'BAND': img.select('NDSI_Snow_Cover_Class')}
             )) \
        .toList(400)
    #display(MODIS_cloud_list,'MODIS_cloud_list')
    
    # This makes a list of images where the snow value is between 1 and 100 based on the snow 
    # cover band. This creates an yes/no (0/1) image of snow locations each day of the year.
    # See printout to console for more info.
    MODIS_snow_list = ee.ImageCollection("MODIS/061/MOD10A1") \
        .filterDate(start,end_adv) \
        .map(lambda img:img.expression(
            '(BAND>='+ NDSI_threshold +'&&BAND<=100)',{'BAND':img.select('NDSI_Snow_Cover')})) \
        .toList(400)
    #display(MODIS_snow_list,'MODIS_snow_list')
    
    # Access the 30m Landsat datset for water loctions. I don't want these metrics
    # to be calculated over inland water bodies!
    # See printout to console for more info.
    water = ee.ImageCollection("GLCF/GLS_WATER") \
        .map(lambda img: img.expression(
            '(BAND==2)', {'BAND':img.select('water')})).Or().Not()
    #display(water,'This is water from LandSat')
    # Map.addLayer(water,{'palette':palette_snow,'min':0, 'max':1}, 'Inland Water 0=water, 1=not water')
    
    
    # Get number of days in the Water Year
    ndays = MODIS_snow_list.length().getInfo() - adv_days
    #print(ndays,'ndays')
    
    # Get number of days in the Water Year
    ndays = MODIS_snow_list.length().getInfo() - adv_days
    #print(ndays,'ndays')
    
    #% This is a correction for leap years.
    #sddCorrection = (ndays < 365 + (wyr % 4 ? 0 : 1)) ? 365 + (wyr % 4 ? 0 : 1) - ndays : 0
    sddCorrection = (365 + (0 if wyr % 3 else 1)) - ndays if ndays < 365 +(0 if wyr % 4 else 1) else 0
    #print('Processed ' + str(ndays) + ' days in ' + str(wyr))
    
    #####################################/
    # Open a bunch of image iables with zero values, or other defined values
    
    scfImage = 'MODIS_SCF_' + str(wyr)
    #print(scfImage,'scfImage')
    sddImage = 'MODIS_SDD_' + str(wyr)
    #print(sddImage,'sddImage')
    sddImage = ee.Image(0)
    
    future_snow = ee.Image(0)
    Snow = ee.Image(0)
    accSnow = ee.Image(0)
    ephemeral_counter = ee.Image(0)
    noSnowCounter = ee.Image(0)
    minNoSnow = ee.Image(minNoSnowdays)
    min_sdays_prior_to_SDD = ee.Image(prior_days)
    sddDetected = ee.Image(0)
    MODIS_max_ndays_nosnow = ee.Image(0)
    
    #############################
    #######/  Do The STUFFFFFF  ###########/
    #############################
    # Increment backwards from the end of the water year + 30 days
    # REMEMBER: This is complex because it is moving backwards in time. Its really easy to forget. 
    # All the commenting is interpreted through the backwards-through-time lens.
    for current_day  in range(ndays + adv_days - 1, -1,-1):
        
        # During those 30 days after the water year, build some 'memory' of where snow is and is not in the iables
        # Take cloud locations from MODIS and create cloud image of 1=cloud, 0=everything else
        Cloud = ee.Image(MODIS_cloud_list.get(current_day)).unmask()
      
        # During those 30 days after the water year, build some 'memory' of where snow is and is not in the iables
        # Take cloud locations from MODIS and create cloud image of 1=cloud, 0=everything else
        Cloud = ee.Image(MODIS_cloud_list.get(current_day)).unmask()
      
        # If there is a cloud designation TODAY and future_snow iable=1, from tomorrow or snow=1 tomorrow, set to 1. 
        # Essentially, 'remembering' where snow is located in the future
        future_snow = Cloud.And(future_snow.Or(Snow))
      
        # Take snow locations from MODIS TODAY and create snow image of 1=snow, 0=everything else for TODAY
        Snow = ee.Image(MODIS_snow_list.get(current_day)).unmask()
      
        # If TODAY is in the water year of interest do all these things
        # Because of advanced days process above, we have snow 'memory' of the future
        if current_day < ndays:
            
            # If TODAY is snow in the pixel or if there is future snow in the pixel, snow_occurrence = 1 
            snow_occurrence = future_snow.Or(Snow)
            
            # Count the number of days snow has occurred in the future in the pixel and then accumulate it to this iable.
            # This is the last time we use accsnow to create the final SCF_image
            accSnow = accSnow.add(snow_occurrence)
            
            ##################################
            # Everything below in the for loop is to deal with the SDD issues.
            # This accesses the iable defined as prior_days and creates a conditional statement
            # If prior_days is greater than the ephemeral counter, and snow does not occur in the future, reset the value of this iable to ephemeral counter value
            # Purpose:
            min_sdays_prior_to_SDD = min_sdays_prior_to_SDD.where(snow_occurrence.eq(0).And(ephemeral_counter.gt(min_sdays_prior_to_SDD)), ephemeral_counter)
            
            # The number of days after the the TODAY that has snow_occurrence=1 (snow or previously snow)
            ephemeral_counter = (ephemeral_counter.add(snow_occurrence)).multiply(snow_occurrence)
            
            # A running count of the number of days when snow_occurrence was equal to 0 from end of water year to now
            # SDD cannot occur until some number of days without snow (minNoSnow) happen at the end of the year
            noSnowCounter = (noSnowCounter.add(snow_occurrence.eq(0))).multiply(sddDetected.eq(0))
            
            # Day of snow disappearance
            # Number of days of snow prior to CURRENT day must excede a value set by prior_days iable
            # And there must be a minimum number of days without snow after the CUURENT day (minNoSnow)
            sddDetected = (ephemeral_counter.gt(min_sdays_prior_to_SDD)).And(noSnowCounter.gt(minNoSnow))
            
            
            minNoSnow = minNoSnow.where(sddDetected, noSnowCounter)
            
            # Find the day of snow disappearance by adding some values to the CURRENT day
            sdd1 = current_day+sddCorrection
            sdd = min_sdays_prior_to_SDD.add(ee.Image(sdd1))
            sddImage = sddImage.where(sddDetected, sdd)
    
        
    sddImage = sddImage.where(accSnow.eq(ndays-1), ndays + sddCorrection)
    
    # Make the final scf image by dividing the number of snow days by the number of days in the year
    # Also, mask out all the inland water bodies using the Landsat dataset (its more precise than the MODIS land mask, for sure)
    scfImage = accSnow.divide(ndays).updateMask(accSnow).rename('SCF').updateMask(water)
    
    # Make the final scf image
    # Also, mask out all the inland water bodies
    sddImage = sddImage.updateMask(sddImage).rename('SDD').updateMask(water)
    
    ########################################
    ####       Mask the images        ######
    ########################################
    # Mask the low_values in SCF.
    scfimg = scfImage
    SCF = scfimg
    #low_end = ee.Image(low_value)
    #lower_mask = scfimg.gte(low_end)
    #SCF = scfimg.updateMask(lower_mask)
    #SCF = SCF
    #Map.addLayer(SCF,{'palette':palette_snow,'min':low_value, 'max':1},'SCF (post-mask)')
    
    # Mask the SDD image by the low_values image above, just so we have the same coverage
    sddimg = sddImage
    #SDD = sddimg.updateMask(lower_mask)
    # Map.addLayer(SDD,{'palette':palette_sdd,'min':0, 'max':365},'SDD (post-mask)')
    
    SDD = sddimg
    
    return ee.Image([SCF, SDD])
