# please just look at the other notebooks that does all of this separately, I just created this file in case the
# assignment auto grades

from pystac.extensions.eo import EOExtension as eo
import pystac_client
import planetary_computer
import glob
import rioxarray as rxr
from rioxarray.merge import merge_arrays
import re
import datetime
import pandas as pd
from shapely.geometry import box
import odc.stac
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import rasterio as rio
from os.path import basename, exists, expanduser, join
import shutil
import shapely
import rioxarray as rxr
import xarray as xr
import time
import gzip
from urllib.request import urlretrieve
import traceback

def cop30_for_aso(aso_raster_fn):
    aso_raster = rxr.open_rasterio(aso_raster_fn).squeeze()
    aso_raster = aso_raster.where(aso_raster>=0, drop=True)
    aso_raster = aso_raster.interpolate_na(dim='x')
    bounds_latlon = box(*aso_raster.rio.transform_bounds("EPSG:4326"))
    
    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1")

    search = catalog.search(
        collections=["cop-dem-glo-30"],
        intersects=bounds_latlon)

    # Check how many items were returned
    items = search.item_collection()
    print(f"Returned {len(items)} Items")
    
    data = []
    for item in items:
        dem_path = planetary_computer.sign(item.assets['data']).href
        data.append(rxr.open_rasterio(dem_path))
    cop30_da = merge_arrays(data)  
    
    # clip to ASO extent
    cop30_da = cop30_da.squeeze().rio.reproject_match(aso_raster, resampling=rio.enums.Resampling.bilinear)
    
    cop30_da.to_netcdf(rf"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\cop30\cop30_for_{aso_raster_fn.split("/")[-1][:-4]}.nc")
    
    #return cop30_da, aso_raster
def cop30_for_aso_all(dir_path):
    error_paths = []
    raster_paths = glob.glob(f'{dir_path}/*/ASO_50M_SD*.tif')[:1]
    for i, path in enumerate(raster_paths):
        print(f'----\nworking on {path.split("/")[-1]}, {i+1}/{len(raster_paths)}\n----')
        
        try:
            cop30_for_aso(path)
        except Exception as e:
            print(e)
            print('encountered error downloading, skipping for now')
            error_paths.append(path)

dir_path = r"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\ASO_50m_SD_cleaned"
cop30_for_aso_all(dir_path)

def url_download(url, out_fp, overwrite = False):
    # check if file already exists
    if not exists(out_fp) or overwrite == True:
            urlretrieve(url, out_fp)
    # if already exists. skip download.
    else:
        print('file already exists, skipping')
def download_fcf(out_fp):
    # this is the url from Lievens et al. 2021 paper
    fcf_url = 'https://zenodo.org/record/3939050/files/PROBAV_LC100_global_v3.0.1_2019-nrt_Tree-CoverFraction-layer_EPSG-4326.tif'
    # download just forest cover fraction to out file
    url_download(fcf_url, out_fp)
fcf_path =r"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\fcf\global\fcf_global.tif"
download_fcf(fcf_path)
def fcf_for_aso(aso_raster_fn, fcf_path):
    aso_raster = rxr.open_rasterio(aso_raster_fn).squeeze()
    aso_raster = aso_raster.where(aso_raster>=0, drop=True)
    aso_raster = aso_raster.interpolate_na(dim='x')
    
    # open as dataArray and return
    fcf = rxr.open_rasterio(fcf_path)

    # clip first to avoid super long reproject processes
    fcf = fcf.rio.clip_box(*aso_raster.rio.transform_bounds("EPSG:4326"))
    # reproject FCF to match dataset
    fcf = fcf.rio.reproject_match(aso_raster, resampling=rio.enums.Resampling.bilinear)
    # remove band dimension as it only has one band
    fcf = fcf.squeeze('band')
    # if max is greater than 1 set to 0-1
    if fcf.max() >= 1:
        fcf = fcf / 100
    
    assert fcf.max() <= 1, "Forest cover fraction must be bounded 0-1"
    assert fcf.min() >= 0, "Forest cover fraction must be bounded 0-1"
    
    fcf.to_netcdf(rf"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\fcf\fcf_for_{aso_raster_fn.split("\\")[-1][:-4]}.nc")

    #return fcf, aso_raster
def fcf_for_aso_all(dir_path, fcf_path):
    raster_paths = glob.glob(rf'{dir_path}\*\ASO_50M_SD*.tif')
    for i, path in enumerate(raster_paths):
        print(f'----\nworking on {path.split("/")[-1]}, {i+1}/{len(raster_paths)}\n----')

        fcf_for_aso(path, fcf_path)
dir_path = r"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\ASO_50m_SD_cleaned"
fcf_for_aso_all(dir_path, fcf_path)

def rtc_for_aso_snowon_mean(aso_raster_fn):
    time = pd.to_datetime(re.search(r"(\d{4}\d{2}\d{2})", aso_raster_fn).group())
    week_before = (time - datetime.timedelta(weeks=2)).strftime('%Y-%m-%d')
    week_after = (time + datetime.timedelta(weeks=2)).strftime('%Y-%m-%d')
    time_of_interest = f'{week_before}/{week_after}'

    aso_raster = rxr.open_rasterio(aso_raster_fn).squeeze()
    aso_raster = aso_raster.where(aso_raster>=0, drop=True)
    bounds_latlon = box(*aso_raster.rio.transform_bounds("EPSG:4326"))

    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace)

    search = catalog.search(
        collections=["sentinel-1-rtc"],
        intersects=bounds_latlon,
        datetime=time_of_interest)

    # Check how many items were returned
    items = search.item_collection()

    rtc_stac = odc.stac.load(items,chunks={"x": 2048, "y": 2048},resolution=50, groupby='sat:absolute_orbit')
    print(f"Returned {len(rtc_stac.time)} acquisitions")
    rtc_stac_clipped = rtc_stac.rio.clip_box(*bounds_latlon.bounds,crs="EPSG:4326")

    rel_orbits = [scene.properties['sat:relative_orbit'] for scene in items.items]
    ac_times = [scene.properties['datetime'] for scene in items.items]
    ac_times = [np.datetime64(item) for item in ac_times]

    # clip to ASO extent
    rtc_stac_clipped = rtc_stac_clipped.rio.reproject_match(aso_raster, resampling=rio.enums.Resampling.bilinear)

    # limit to morning acquisitions
    rtc_ds = rtc_stac_clipped.where(rtc_stac_clipped.time.dt.hour > 11, drop=True)
    if 'vv' not in list(rtc_ds.keys()) or 'vh' not in list(rtc_ds.keys()):
        print('missing polarization')
        return None

    if len(rtc_ds.time) == 0:
        print('no morning acquisitions')
        return None

    # calculate percent vh coverage of each acquisition
    perc_cover = (rtc_ds.vh.where(aso_raster >= 0) > 0).sum(dim=['x', 'y'])/(rtc_ds.vh.where(aso_raster >= 0) >= -1000000000).sum(dim=['x', 'y'])

    # if multiple with full coverage, grab nearest in time with full coverage
    if perc_cover.values.tolist().count(1) > 1:
        print('total snow-on coverage available')
        rtc_ds = rtc_ds.where(perc_cover == 1, drop=True).sortby('time')
        rtc_ds = rtc_ds.sel(time=time, method='nearest')
        # should probably redo with largest number of scenes with full coverage

    # exit if no rasters have good vh coverage
    elif perc_cover.max() < 0.01:
        print('max vh coverage is < 1%--recommend skipping ASO raster')
        return None

    # otherwise, grab max coverage 
    else:
        if perc_cover.max() == 1:
            print('total snow-on coverage available')
        else: 
            print(f'{perc_cover.max().item()} snow-on coverage')
        rtc_ds = rtc_ds.sel(time=perc_cover.idxmax())

    # get relative orbit of scene
    rel_orbit = rel_orbits[ac_times.index(rtc_ds.time)]

    orbit_dict = {}
    for i, orbit in enumerate(rel_orbits):
        if orbit not in orbit_dict.keys():
            orbit_dict[orbit] = [ac_times[i]]
        else:
            orbit_dict[orbit].append(ac_times[i])

    rtc_ds = rtc_stac_clipped.where(rtc_stac_clipped.time.isin(orbit_dict[rel_orbit]), drop=True)
    
    # mask negative areas
    rtc_ds = rtc_ds.where(rtc_ds.vh > 0)
    rtc_ds = rtc_ds.where(rtc_ds.vv > 0)

    # take mean of all acquisitions
    print(f'taking mean of {rtc_ds.time.size} snow-on rasters')
    rtc_ds = rtc_ds.mean(dim='time', skipna=True)

    rtc_ds = rtc_ds.compute()

    rtc_ds.to_netcdf(rf"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\S1_rtc_mean\S1_snow-on_orbit{rel_orbit}_for_{aso_raster_fn.split("\\")[-1][:-4]}.nc")
    
    return rel_orbit
def rtc_for_aso_snowoff_mean(aso_raster_fn, rel_orbit):
    year = pd.to_datetime(re.search(r"(\d{4}\d{2}\d{2})", aso_raster_fn).group()).year
    time = pd.to_datetime(f'{year-1}0910')
    week_before = (time - datetime.timedelta(weeks=2)).strftime('%Y-%m-%d')
    week_after = (time + datetime.timedelta(weeks=2)).strftime('%Y-%m-%d')
    time_of_interest = f'{week_before}/{week_after}'

    aso_raster = rxr.open_rasterio(aso_raster_fn).squeeze()
    aso_raster = aso_raster.where(aso_raster>=0, drop=True)
    bounds_latlon = box(*aso_raster.rio.transform_bounds("EPSG:4326"))

    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace)

    search = catalog.search(
        collections=["sentinel-1-rtc"],
        intersects=bounds_latlon,
        datetime=time_of_interest)

    # Check how many items were returned
    items = search.item_collection()

    rel_orbits = [scene.properties['sat:relative_orbit'] for scene in items.items]
    ac_times = [scene.properties['datetime'] for scene in items.items]
    ac_times = [np.datetime64(item) for item in ac_times]

    rtc_stac = odc.stac.load(items,chunks={"x": 2048, "y": 2048},resolution=50, groupby='sat:absolute_orbit')
    print(f"Returned {len(rtc_stac.time)} acquisitions")
    rtc_stac_clipped = rtc_stac.rio.clip_box(*bounds_latlon.bounds,crs="EPSG:4326")

    orbit_dict = {}
    for i, orbit in enumerate(rel_orbits):
        if orbit not in orbit_dict.keys():
            orbit_dict[orbit] = [ac_times[i]]
        else:
            orbit_dict[orbit].append(ac_times[i])

    if rel_orbit not in orbit_dict.keys():
        print('no acquisitons from same orbit, skipping')
        return

    rtc_stac_clipped = rtc_stac_clipped.where(rtc_stac_clipped.time.isin(orbit_dict[rel_orbit]), drop=True)

    # clip to ASO extent
    rtc_ds = rtc_stac_clipped.rio.reproject_match(aso_raster, resampling=rio.enums.Resampling.bilinear)

    if 'vv' not in list(rtc_ds.keys()) or 'vh' not in list(rtc_ds.keys()):
        print('missing polarization, skipping')
        return
    
    if len(rtc_ds.time) == 0:
        print('no morning acquisitions')
        return 

    # calculate percent vh coverage of each acquisition
    perc_cover = (rtc_ds.vh.where(aso_raster >= 0) > 0).sum(dim=['x', 'y'])/(rtc_ds.vh.where(aso_raster >= 0) >= -1000000000).sum(dim=['x', 'y'])
    
    # if multiple with full coverage, grab nearest in time with full coverage
    if perc_cover.values.tolist().count(1) > 1:
        print('total snow-on coverage available')

    # exit if no rasters have good vh coverage
    elif perc_cover.max() < 0.01:
        print('max vh coverage is < 1%--recommend skipping ASO raster')
        return
    
    else:
        if perc_cover.max() == 1:
            print('total snow-on coverage available')
        else: 
            print(f'{perc_cover.max().item()} snow-on coverage')
            
    # mask negative areas
    rtc_ds = rtc_ds.where(rtc_ds.vh > 0)
    rtc_ds = rtc_ds.where(rtc_ds.vv > 0)
    
    # take mean of all acquisitions
    print(f'taking mean of {rtc_ds.time.size} snow-off rasters')
    rtc_ds = rtc_ds.mean(dim='time', skipna=True)

    rtc_ds = rtc_ds.compute()
    
    rtc_ds.to_netcdf(rf"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\S1_rtc_mean\S1_snow-off_orbit{rel_orbit}_for_{aso_raster_fn.split("\\")[-1][:-4]}.nc")
def rtc_for_aso_all(dir_path):
    raster_paths = glob.glob(rf'{dir_path}\*\ASO_50M_SD*.tif')
    for i, path in enumerate(raster_paths):
        print(f'----\nworking on {path.split("/")[-1]}, {i+1}/{len(raster_paths)}\n----')
        
        try:
            relative_orbit = rtc_for_aso_snowon_mean(path)
            if relative_orbit == None:
                continue
            rtc_for_aso_snowoff_mean(path, relative_orbit)
        except Exception as exc:
            print(traceback.format_exc())
            print(exc)
            print('encountered error, skipping')
dir_path = r"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\ASO_50m_SD_cleaned"
test = rtc_for_aso_all(dir_path)

def rtc_for_aso_snowon(aso_raster_fn):
    time = pd.to_datetime(re.search("(\d{4}\d{2}\d{2})", aso_raster_fn).group())
    week_before = (time - datetime.timedelta(weeks=2)).strftime('%Y-%m-%d')
    week_after = (time + datetime.timedelta(weeks=2)).strftime('%Y-%m-%d')
    time_of_interest = f'{week_before}/{week_after}'
    
    aso_raster = rxr.open_rasterio(aso_raster_fn).squeeze()
    aso_raster = aso_raster.where(aso_raster>=0, drop=True)
    bounds_latlon = box(*aso_raster.rio.transform_bounds("EPSG:4326"))
    
    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace)

    search = catalog.search(
        collections=["sentinel-1-rtc"],
        intersects=bounds_latlon,
        datetime=time_of_interest)
    
    # Check how many items were returned
    items = search.item_collection()
    
    rtc_stac = odc.stac.load(items,chunks={"x": 2048, "y": 2048},resolution=50, groupby='sat:absolute_orbit')
    print(f"Returned {len(rtc_stac.time)} acquisitions")
    rtc_stac_clipped = rtc_stac.rio.clip_box(*bounds_latlon.bounds,crs="EPSG:4326")
    
    rel_orbits = [scene.properties['sat:relative_orbit'] for scene in items.items]
    ac_times = [scene.properties['datetime'] for scene in items.items]
    ac_times = [np.datetime64(item) for item in ac_times]
    
    # clip to ASO extent
    rtc_stac_clipped = rtc_stac_clipped.rio.reproject_match(aso_raster, resampling=rio.enums.Resampling.bilinear)

    # limit to morning acquisitions
    rtc_ds = rtc_stac_clipped.where(rtc_stac_clipped.time.dt.hour > 11, drop=True)
    if 'vv' not in list(rtc_ds.keys()) or 'vh' not in list(rtc_ds.keys()):
        print('missing polarization')
        return None
    
    if len(rtc_ds.time) == 0:
        print('no morning acquisitions')
        return None
        
    # calculate percent vh coverage of each acquisition
    perc_cover = (rtc_ds.vh.where(aso_raster >= 0) > 0).sum(dim=['x', 'y'])/(rtc_ds.vh.where(aso_raster >= 0) >= -1000000000).sum(dim=['x', 'y'])
    
    # if multiple with full coverage, grab nearest in time with full coverage
    if perc_cover.values.tolist().count(1) > 1:
        print('total snow-on coverage available')
        rtc_ds = rtc_ds.where(perc_cover == 1, drop=True).sortby('time')
        rtc_ds = rtc_ds.sel(time=time, method='nearest')

    # exit if no rasters have good vh coverage
    elif perc_cover.max() < 0.01:
        print('max vh coverage is < 1%--recommend skipping ASO raster')
        return None

    # otherwise, grab max coverage 
    else:
        if perc_cover.max() == 1:
            print('total snow-on coverage available')
        else: 
            print(f'{perc_cover.max().item()} snow-on coverage')
        rtc_ds = rtc_ds.sel(time=perc_cover.idxmax())

    # get relative orbit of scene
    rel_orbit = rel_orbits[ac_times.index(rtc_ds.time)]
    
    # mask negative areas
    rtc_ds = rtc_ds.where(rtc_ds.vh > 0)
    rtc_ds = rtc_ds.where(rtc_ds.vv > 0)
    
    rtc_ds.to_netcdf(rf"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\S1_rtc\S1_snow-on_{rtc_ds.time.dt.strftime("%Y%m%d").item()}_for_{aso_raster_fn.split("\\")[-1][:-4]}.nc")
    
    return rel_orbit
def rtc_for_aso_snowoff(aso_raster_fn, rel_orbit):
    year = pd.to_datetime(re.search("(\d{4}\d{2}\d{2})", aso_raster_fn).group()).year
    time = pd.to_datetime(f'{year-1}0910')
    week_before = (time - datetime.timedelta(weeks=2)).strftime('%Y-%m-%d')
    week_after = (time + datetime.timedelta(weeks=2)).strftime('%Y-%m-%d')
    time_of_interest = f'{week_before}/{week_after}'

    aso_raster = rxr.open_rasterio(aso_raster_fn).squeeze()
    aso_raster = aso_raster.where(aso_raster>=0, drop=True)
    bounds_latlon = box(*aso_raster.rio.transform_bounds("EPSG:4326"))

    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace)

    search = catalog.search(
        collections=["sentinel-1-rtc"],
        intersects=bounds_latlon,
        datetime=time_of_interest)

    # Check how many items were returned
    items = search.item_collection()

    rel_orbits = [scene.properties['sat:relative_orbit'] for scene in items.items]
    ac_times = [scene.properties['datetime'] for scene in items.items]
    ac_times = [np.datetime64(item) for item in ac_times]

    rtc_stac = odc.stac.load(items,chunks={"x": 2048, "y": 2048},resolution=50, groupby='sat:absolute_orbit')
    print(f"Returned {len(rtc_stac.time)} acquisitions")
    rtc_stac_clipped = rtc_stac.rio.clip_box(*bounds_latlon.bounds,crs="EPSG:4326")

    orbit_dict = {}
    for i, orbit in enumerate(rel_orbits):
        if orbit not in orbit_dict.keys():
            orbit_dict[orbit] = [ac_times[i]]
        else:
            orbit_dict[orbit].append(ac_times[i])

    if rel_orbit not in orbit_dict.keys():
        print('no acquisitons from same orbit, skipping')
        return

    rtc_stac_clipped = rtc_stac_clipped.where(rtc_stac_clipped.time.isin(orbit_dict[rel_orbit]), drop=True)

    # clip to ASO extent
    rtc_stac_clipped = rtc_stac_clipped.rio.reproject_match(aso_raster, resampling=rio.enums.Resampling.bilinear)

    # limit to morning acquisitions
    rtc_ds = rtc_stac_clipped.where(rtc_stac_clipped.time.dt.hour > 11, drop=True)
    if 'vv' not in list(rtc_ds.keys()) or 'vh' not in list(rtc_ds.keys()):
            print('missing polarization, skipping')
            return
    
    if len(rtc_ds.time) == 0:
        print('no morning acquisitions')
        return None

    # calculate percent vh coverage of each acquisition
    perc_cover = (rtc_ds.vh.where(aso_raster >= 0) > 0).sum(dim=['x', 'y'])/(rtc_ds.vh.where(aso_raster >= 0) >= -1000000000).sum(dim=['x', 'y'])

    # if multiple with full coverage, grab nearest in time with full coverage
    if perc_cover.values.tolist().count(1) > 1:
        print('total snow-off coverage available')
        rtc_ds = rtc_ds.where(perc_cover == 1, drop=True).sortby('time')
        rtc_ds = rtc_ds.sel(time=time, method='nearest')

    # exit if no rasters have good vh coverage
    elif perc_cover.max() < 0.01:
        print('max vh coverage is < 1%--recommend skipping ASO raster')
        return

    # otherwise, grab max coverage 
    else:
        if perc_cover.max() == 1:
            print('total snow-off coverage available')
        else: 
            print(f'{perc_cover.max().item()} snow-off coverage')
        rtc_ds = rtc_ds.sel(time=perc_cover.idxmax())
    
    # mask negative areas
    rtc_ds = rtc_ds.where(rtc_ds.vh > 0)
    rtc_ds = rtc_ds.where(rtc_ds.vv > 0)
    
    rtc_ds.to_netcdf(rf"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\S1_rtc\S1_snow-off_{rtc_ds.time.dt.strftime("%Y%m%d").item()}_for_{aso_raster_fn.split("\\")[-1][:-4]}.nc")

def rtc_for_aso_all(dir_path):
    raster_paths = glob.glob(rf'{dir_path}\utm11n\ASO_50M_SD*.tif')
    for i, path in enumerate(raster_paths):
        print(f'----\nworking on {path.split("/")[-1]}, {i+1}/{len(raster_paths)}\n----')
        
        try:
            relative_orbit = rtc_for_aso_snowon(path)
            if relative_orbit == None:
                continue
            rtc_for_aso_snowoff(path, relative_orbit)
        except Exception as exc:
            print(traceback.format_exc())
            print(exc)
            print('encountered error, skipping')
dir_path = r"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\ASO_50m_SD_cleaned"
test = rtc_for_aso_all(dir_path)

def sentinel2_for_aso(aso_raster_fn):
    
    time = pd.to_datetime(re.search(r"(\d{4}\d{2}\d{2})", aso_raster_fn).group())
    week_before = (time - datetime.timedelta(weeks=0.2)).strftime('%Y-%m-%d')
    week_after = (time + datetime.timedelta(weeks=0.2)).strftime('%Y-%m-%d')
    time_of_interest = f'{week_before}/{week_after}'
    
    aso_raster = rxr.open_rasterio(aso_raster_fn).squeeze()
    bounds_latlon = box(*aso_raster.rio.transform_bounds("EPSG:4326"))
    
    catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace)

    search = catalog.search(
        collections=["sentinel-2-l2a"],
        intersects=bounds_latlon,
        datetime=time_of_interest)

    # Check how many items were returned
    items = search.item_collection()
    print(f"Returned {len(items)} Items")
    
    sentinel2_stac = odc.stac.load(items,chunks={"x": 2048, "y": 2048},resolution=50, groupby='solar_day')
    sentinel2_stac_clipped = sentinel2_stac.rio.clip_box(*bounds_latlon.bounds,crs="EPSG:4326")
    scl = sentinel2_stac_clipped['SCL'].rio.reproject_match(aso_raster, resampling=rio.enums.Resampling.bilinear).where(aso_raster>=0)
    classes = [ #SCL classes here: https://custom-scripts.sentinel-hub.com/custom-scripts/sentinel-2/scene-classification/
    #0,   #No Data (Missing data)	#000000	
    #1,   #Saturated or defective pixel	#ff0000	
    #2,   #Topographic casted shadows (called "Dark features/Shadows" for data before 2022-01-25)	#2f2f2f	
    #3,   #Cloud shadows	#643200	
    4,   #Vegetation	#00a000	
    5,   #Not-vegetated	#ffe65a	
    #6,   #Water	#0000ff	
    #7,   #Unclassified	#808080	
    #8,   #Cloud medium probability	#c0c0c0	
    #9,   #Cloud high probability	#ffffff	
    #10,   #Thin cirrus	#64c8ff	
    11    #Snow or ice      
    ]
    
    idx_least_clouds = scl.where(scl.isin(classes)).sum(dim=['x','y']).idxmax()
    sentinel2_best_lowcloud = sentinel2_stac_clipped.sel(time=idx_least_clouds)
    sentinel2_best_lowcloud.to_netcdf(rf'C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\sentinel_2\{pd.to_datetime(idx_least_clouds.values).strftime("%Y%m%d")}_for_{aso_raster_fn.split("\\")[-1][:-4]}.nc')    
    #return sentinel2_best_lowcloud
aso_raster_fns = glob.glob(r"C:\Users\JackE\uw\courses\aut24\ml_geo\final_data\ASO_50m_SD_cleaned\utm12n\ASO*.tif")
for i, aso_raster_fn in enumerate(aso_raster_fns):
    error_list = []
    print(f'----\nworking on {aso_raster_fn.split("/")[-1]}, {i+1}/{len(aso_raster_fns)}\n----')
    try: 
        sentinel2_for_aso(aso_raster_fn)
    except:
        print('error, skipping')
        error_list.append(aso_raster_fn)