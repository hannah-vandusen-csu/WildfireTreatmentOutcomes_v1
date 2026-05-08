# updated to use UNIX timestamp format burn rasters

import ee
import geemap
import pandas as pd
from datetime import datetime, timedelta


# Create date object from year, month, and day
def add_date(feature):
    date = ee.Date.fromYMD(feature.get('year'), feature.get('month'), feature.get('day'))
    return feature.set('ee_date', date)


def extract_dob(
    fireID,
    shape_features,      # ee.FeatureCollection
    fire_list_pd,        # pandas DataFrame with ignition year, MTBS_ID, etc.
    local_directory,     # not directly used (kept for compatibility)
    pixel_area = False,    # if True, computes pixel area burned by date
    timezone ='US/Pacific'  # <-- NEW: specify timezone for Unix conversion; default is US/Arizona
):
    """
    Extracts date of burn (DOB) from GEE raster.

    Args:
        fireID (str): Fire identifier.
        shape_features (ee.FeatureCollection): Polygons or points of interest.
        fire_list_pd (pd.DataFrame): Metadata table with fire info.
        local_directory (str): Placeholder (kept for compatibility).
        pixel_area (bool): If True, returns burned area by DOB.
        timezone (str): Timezone for Unix time conversion (default = 'US/Arizona').
    """

    # Get fire year from metadata
    year = str(fire_list_pd[fire_list_pd['MTBS_ID'] == fireID]['Year'].iloc[0])

    # Load Date of Burn (DOB) raster asset for the given fire
    dob_raw_fp = f"projects/rmrs-wildfire-treatments/assets/klamath/dob_parks/{fireID}_1_dob"
    dob_raw = ee.Image(dob_raw_fp)

    # The DOB raster stores Unix timestamps (seconds since 1970-01-01 UTC)
    dob_unix = dob_raw.select('b1').toInt64()

    # Extract the minimum (earliest) burn timestamp under each feature
    shape_dob = dob_unix.reduceRegions(
        reducer=ee.Reducer.min(),
        collection=shape_features,
        scale=30
    )

    # Convert DOB to a Feature Collection and a Date Object the date strings a 6-digit integer representing YYMMDD for years >= 2010 or a 5 digit integer (YMMDD) for 2000 <= years < 2010 
    if int(year) >= 2010:
        shape_ymd = shape_dob.map(lambda feature:
                                  feature.set('day', ee.Number.parse(ee.String(feature.get('min')).slice(4, 6)))
                                  .set('year', ee.Number.parse(ee.String(feature.get('min')).slice(0, 2)).add(2000))
                                  .set('month', ee.Number.parse(ee.String(feature.get('min')).slice(2, 4))))
    else:
        shape_ymd = shape_dob.map(lambda feature:
                                feature.set('day', ee.Number.parse(ee.String(feature.get('min')).slice(3, 5)))
                                .set('year', ee.Number.parse(ee.String(feature.get('min')).slice(0, 1)).add(2000))
                                .set('month', ee.Number.parse(ee.String(feature.get('min')).slice(1, 3))))
    
    dob_date = shape_ymd.map(add_date)

    # Default behavior: return FeatureCollection with ee_date
    if not pixel_area:
        return dob_date

    # Optional: compute pixel area by burn date
    else:
        areaImage = ee.Image.pixelArea().addBands(dob_unix)
        areas = areaImage.reduceRegion(
            reducer=ee.Reducer.sum().group(groupField=1, groupName='class'),
            geometry=shape_features.geometry(),
            scale=30,
            maxPixels=1e10
        )
        dob_areas = ee.List(areas.get('groups'))

        def map_function(item):
            area_dict = ee.Dictionary(item)
            class_number = ee.Number(area_dict.get('class')).format()
            area = ee.Number(area_dict.get('sum'))
            return ee.List([class_number, area])

        dob_area_lists = dob_areas.map(map_function)
        dob_dict = ee.Dictionary(dob_area_lists.flatten())
        return pd.DataFrame.from_dict(dob_dict.getInfo().items(), orient='columns')



# Helper: Extract weather for a given feature/date)

def extract_weather_data(weather_collection, advance_day, feature, band, scale, reducer, sub_daily_reducer, sub_daily_reducer_string):
    date = ee.Date(feature.get('ee_date'))

    if sub_daily_reducer == '':
        img = weather_collection.filterDate(date, date.advance(advance_day, 'day')).first().select(band)
    else:
        img = (weather_collection
               .filterDate(date, date.advance(advance_day, 'day'))
               .select(band)
               .reduce(sub_daily_reducer))
        band = band + '_' + sub_daily_reducer_string

    vector_value = img.reduceRegion(
        reducer=reducer,
        geometry=feature.geometry(),
        scale=scale,
        bestEffort=True
    ).get(band)

    return feature.set(band, ee.Number(vector_value))


# Main Function: Extract weather variables (daily or subdaily)
def daily_weather(
    sample_dob_fc,           # FeatureCollection with ee_date (from extract_dob)
    var_crosswalks,          # Variable list DataFrame with file_path, Abbreviation, etc.
    reducer=ee.Reducer.mean(),
    sub_daily_reducer='',
    sub_daily_reducer_string='',
    timezone='US/Pacific'    # Optional timezone argument
):
    """
    Extracts weather variables for features with known burn date/time.
    Supports multiple datasets (e.g., GRIDMET and RTMA) in a single run.

    Args:
        sample_dob_fc (ee.FeatureCollection): Input features with ee_date.
        var_crosswalks (pd.DataFrame): Variable list with at least
            ['file_path', 'Abbreviation', 'scale'] columns.
        reducer (ee.Reducer): Spatial reducer (default = mean).
        sub_daily_reducer (ee.Reducer): Optional temporal reducer (for multi-image composites).
        sub_daily_reducer_string (str): Suffix to append when sub-daily reducer applied.
        timezone (str): Timezone for RTMA window extraction (default = 'US/Arizona').

    Returns:
        ee.FeatureCollection: Input features with new weather variable attributes.
    """

    # Get all unique datasets from the variable crosswalk
    datasets = var_crosswalks['file_path'].unique().tolist()

    # Iterate through each dataset (can include both GRIDMET and RTMA)
    for gee_data_source in datasets:

        # Subset the variable list for this dataset
        col_var_df = var_crosswalks[var_crosswalks['file_path'] == gee_data_source]
        bands = col_var_df['Abbreviation'].tolist()
        col_scale = col_var_df['scale'].unique()[0]


        # Set advance_day for non-RTMA datasets
        if 'GRIDMET/DROUGHT' in gee_data_source:
            advance_day = 5
        else:
            advance_day = 1

        # Load the GEE ImageCollection
        weather_collection = ee.ImageCollection(gee_data_source)


        # GRIDMET 
        for band in bands:
            sample_dob_fc = sample_dob_fc.map(
                lambda f, b=band: extract_weather_data(
                    weather_collection=weather_collection,
                    advance_day=advance_day,
                    feature=f,
                    band=b,
                    scale=col_scale,
                    reducer=reducer,
                    sub_daily_reducer=sub_daily_reducer,
                    sub_daily_reducer_string=sub_daily_reducer_string
                )
            )

    return sample_dob_fc