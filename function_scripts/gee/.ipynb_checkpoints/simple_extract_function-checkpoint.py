
import ee

ee.Authenticate(auth_mode="gcloud")
ee.Initialize(project='rmrs-wildfire-treatments')

import geemap
import geopandas as gpd
import pandas as pd
import numpy as np
from datetime import datetime
from datetime import timedelta
import os
import sys
import math
import sklearn
from sklearn.preprocessing import StandardScaler
import importlib
import re

# import functions
import mask_function
import burn_severity_veg_indices
import landcover_functions
import topography_functions
import snow_functions
import gcs_functions
import daily_weather
import climate_functions


def process_fire_data(
        fireID: str,
        fire_list: list,
        local_directory: str,
        wto_outcome: str,
        variable_list: pd.DataFrame,
        bucket_name: str,
        bucket_fp: str,
        output_directory: str,
        reducer: ee.Reducer,
        filename_suffix: str,
        shape_directory: str = None,
        feature_collection: ee.FeatureCollection = None,
        return_feature: bool = False,
        mask: bool = True,
        timezone: str = 'US/Pacific',
        start_day: int = None,
        end_day: int = None,
):
    """
    Processes fire data for a given fire ID, extracting various environmental variables and exporting the results to Google Cloud Storage.
    It either loads a local shapefile or uses a provided Earth Engine FeatureCollection.

    Args:
        fireID (str): The ID of the fire.
        fire_list (list): A dataframe of all fires with columns: Ignition, MTBS_ID, and Name.
        local_directory (str): The local directory for intermediate files.
        wto_outcome (str): WTO outcome
        variable_list (pd.DataFrame): DataFrame containing variable information from Google Earth Engine.
        bucket_name (str): The name of the Google Cloud Storage bucket.
        bucket_fp (str): The file path within the bucket.
        output_directory (str): The output directory within the bucket for GEE and within the local output directory after download.
        reducer (ee.Reducer): The reducer to use for region extraction.
        filename_suffix (str): The suffix for the sampling shapefile name.
        shape_directory (str): The directory containing the sampling shapefiles.
        feature_collection (ee.FeatureCollection): instead of loading in a load shapefile load feature collection. Only fill in shape_directory or feature_collection.
        return_feature (bool): return an extracted feature for the fire of interest. Default (False)
        mask (bool): Mask pixels before extraction default (True)
        timezone (str): The timezone in TZ identifier format of the fire; defaults to 'US/Arizona'

    """
    # Add year to fire list
    fire_list['Year'] = pd.to_datetime(fire_list['Ignition']).dt.year

    if shape_directory:
        shapefile_path = os.path.join(
            shape_directory, fireID + filename_suffix + ".shp"
        )
        if not os.path.exists(shapefile_path):
            print(f"Skipping {fireID} fire: Shapefile not found at {shapefile_path}")
            return
        print(f"Processing {fireID} fire using local shapefile.")

        shape_gdf = gpd.read_file(shapefile_path)
        shape_feat = geemap.geopandas_to_ee(shape_gdf)  # shapes_dob is now a FeatureCollection

    elif feature_collection:
        print(f"Processing {fireID} fire using provided FeatureCollection.")
        shape_feat = feature_collection  # Use the FeatureCollection directly
    else:
        print(f"Skipping {fireID} fire")
        raise ValueError(
            "Either shape_directory and filename_suffix or a feature_collection must be provided."
        )

    sys.setrecursionlimit(16385)

    # Filter Variable dataframe
    var_outcome_df = variable_list[variable_list[wto_outcome] == True]
    vars_only_df = var_outcome_df[var_outcome_df['Group'] != 'Masks']
    var_gee_df = vars_only_df[vars_only_df['extraction_source'] == 'gee']

    # Extract DOB dependent variables
    var_dob_df = var_gee_df[var_gee_df["dob_dependent"] == "yes"]

    # Only calculate point dob if needed
    if len(var_dob_df) > 0:
        dob_datasources = var_dob_df["file_path"].unique().tolist()
        # Calculate Day of Burn and Variables that are dependent on dob like fire weather
        shapes_dob = daily_weather.extract_dob(
            fireID, shape_feat, fire_list, local_directory, pixel_area=False
        )

        # Extract DOB dependent variables
        var_dob_df = var_gee_df[var_gee_df["dob_dependent"] == "yes"]
        dob_datasources = var_dob_df["file_path"].unique().tolist()

        for data in dob_datasources:
            shapes_dob = daily_weather.daily_weather(
                sample_dob_fc=shapes_dob,
                var_crosswalks=var_dob_df,
                timezone=timezone
            )
    else:
        shapes_dob = shape_feat

    # Process variables not dependent on DOB, grouped by CRS
    var_gee_df2 = var_gee_df[var_gee_df["dob_dependent"] == "no"]
    unique_crs = var_gee_df2["crs"].unique().tolist()

    for crs in unique_crs:
        var_proj = var_gee_df2[var_gee_df2["crs"] == crs]
        proj_img_ls = []  # empty list to store data

        # Helper function to extract and append images
        def extract_images(group_name, data_processor):
            if len(var_proj[var_proj["Group"] == group_name]) > 0:
                datasources = var_proj[var_proj["Group"] == group_name][
                    "file_path"
                ].unique().tolist()
                for data in datasources:
                    img = data_processor(data, var_proj)
                    proj_img_ls.append(img)

        # Burn Severity
        def process_burn_severity(data, var_proj):
            if data == "parks 2018":
                return burn_severity_veg_indices.process_fire(
                    fireID, fire_list, local_directory, var_proj[var_proj["Group"] == "Burn Severity"], mtbs_perim=False,
                    start_day = start_day, 
                    end_day = end_day
                )
            elif data == "USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1":
                return landcover_functions.mtbs_severity(fireID, fire_list, data)
            return None  # Handle unexpected data values

        extract_images("Burn Severity", process_burn_severity)

        # Climate Variables
        def process_climate(data, var_proj):
            if data == "crumley 2020":
                return snow_functions.new_snow_metrics(fireID, fire_list)
            else:
                return climate_functions.climate_data_91(fireID, data, var_proj)

        extract_images("Climate", process_climate)

        # Topography variables
        def process_topography(data, var_proj):
            return topography_functions.topography_variables(fireID, data)

        extract_images("Topography", process_topography)

        # Vegetation Variables
        def process_vegetation(data, var_proj):
            if data == "USGS/NLCD_RELEASES":
                return landcover_functions.export_nlcd_for_fire(fireID, fire_list)
            elif data == "parks 2018":
                return burn_severity_veg_indices.prefire_spectral_imagery(
                    fireID, fire_list, local_directory, var_proj[var_proj["Group"] == "Vegetation"]
                )
            elif "EPA/Ecoregions" in data:
                epa_level = re.search(r"\d+$", data).group(0)
                return landcover_functions.epa_ecoregions(fireID, int(epa_level))
            elif "LANDFIRE/Vegetation/BPS" in data:
                return landcover_functions.landfire("BPS")
            elif "vegetation/GEDI_Canopy_Products" in data:
                return burn_severity_veg_indices.prefire_gedi_canopy(
                    fireID, fire_list, local_directory, var_proj[var_proj["Group"] == "Vegetation"]
                )
            return None

        extract_images("Vegetation", process_vegetation)

        # Condense into one image
        all_img = ee.Image.cat(proj_img_ls)


        # Utilize masking function if mask is true
        if mask:
            all_img = mask_function.mask_function(all_img, local_directory, fireID, fire_list, variable_list,
                                                  wto_outcome)

        if all_img:
            # Extract at the minimum scale of this projection
            var_scale = var_gee_df2[var_gee_df2["crs"] == crs]["scale"].min()

            # # --- NEW: apply bilinear interpolation + 90 m focal mean ---
            # all_img_bilinear = (
            #     all_img
            #     .resample('bilinear')  # interpolate pixel values
            #     .focal_mean(radius=90, units='meters')  # average over 90 m neighborhood
            # )

        shapes_dob = shapes_dob.map(
            lambda f: f.setMulti(
                all_img.reduceRegion(
                    reducer=reducer,
                    geometry=f.geometry(),
                    scale=var_scale,
                    maxPixels=1e13
                )
            ).copyProperties(f, f.propertyNames())
        )

    export_task = ee.batch.Export.table.toCloudStorage(
        collection=shapes_dob,
        description="exporting " + fireID + " table",
        bucket=bucket_name,
        fileNamePrefix=bucket_fp
                       + "/"
                       + output_directory
                       + "/"
                       + fireID
                       + filename_suffix
                       + "_df",
        fileFormat="CSV",
        priority=500
    )
    export_task.start()

    if (return_feature):
        return (shapes_dob)