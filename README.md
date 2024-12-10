# ESS 569 A Au 24: Machine Learning In Geosciences Final Project
It should be noted that the basis of this project is forked from a [repo](https://github.com/geo-smart/deep-snow) developed by Quinn Brencher (gbrench@uw.edu) and Eric Gagliano (egagli@uw.edu) that documents some of their incredible contributions to snow research. 
The data used in this analysis is private and was used for an exploration of geospatial data processing pipelines and CNN architectures for our final project in ESS 569. Real progress in this scientific problem can be seen through Brencher et. al's work.

### Group members
Jack Hayes, Ayush Gupta, Sry Wei

*old repo with contributions can be found at https://github.com/Jack-Hayes/mlgeo-2024-deep-snow

### Motivation
Understanding snow depth is crucial in hydrological risk assessment, water resource management, climate change modeling, and more. Remote sensing technologies such as light detection and ranging (LiDAR), synthetic aperture radar (SAR), and optical imagery allow for measurements of snow depth, land cover, and topography across a spatial scale unachievable by traditional manual measurements and models. Current machine learning models that use remote sensing data to measure snow depth are making great strides, but struggle in terms of accuracy at a large spatial scale. The incorporation of spatially-sparse, highly-precise snow depth stations into these models to improve this accuracy is a challenge many snow scientists are facing today. We hope to develop methodology that efficiently encodes point and raster data into machine learning architectures, using Quinn and Eric's "deep-snow" data and models (CNNs).

### Installation
```
$ conda install mamba -n base -c conda-forge
```
Clone the repo and set up the environment
```
$ git clone https://github.com/UW-MLGEO/MLGEO2024_SnowDepth.git
$ cd ./MLGEO2024_SnowDepth
$ mamba env create -f environment.yml
$ conda activate deep-snow
```
Install the package locally
```
$ pip install -e .
```

### Data
*The below is copied from Quinn and Eric's repo
- Sentinel-1 RTC backscatter data (snow on and snow off)
- Sentinel-2 imagery (snow on)
- Fractional forest cover
- COP30 digital elevation model
- Airborne Snow Observatory (ASO) lidar snow depth maps
- SNOTEL snowpack monitoring stations

Snow-on Sentinel-1 and 2 data were collected nearby in time to corresponding ASO acquistions. All products were reprojected to the appropriate UTM zone and resampled to a matching 50 m grid. Products were divided up spatially into training, testing, and validation tiles and subset to produce a machine-learning ready dataset. Our training dataset includes ~37,000 image stacks, each of which includes all of the above listed inputs.  

### MLGeo-AI-Ready-DS
It should be noted that we unfortunately cannot share our cleaned ASO LiDAR snow depth data due to ASO's data sharing policy. The data downloading notebooks do rely on this data so this is not entirely replicable (though accessing the other data sources through the API endpoints with any geometry is still possible, you just won't have the 'ground truth' data to test against). Furthermore, we have too much data to host on our repo so the notebooks below (particularly [notebooks/visualizations/ai-ready_viz.ipynb](notebooks/visualizations/ai-ready_viz.ipynb) will showcase the data).
The AI-ready dataset is split into tiles corresponding to ASO LiDAR flight extents. Each tile (32x32km or 1024km2) will have multiple subdivided rasters of size 128x128 pixels where each pixel is 50m (! 6.4x6.4km). Each of these pixels will contain data listed from the various remote sensing sources below. Tiles do not all have the same amount of rasters within due to missing data. In total, there are 451 testing rasters, 3366 training rasters, and 589 validation rasters. There are 16 testing tiles, 125 training tiles, and 16 validation tiles. The data covers the Sierra Nevada, Rocky Mountain, and Olympic mountain ranges from 2016-2023.
- data
  - [Airborne Snow Observatory (ASO) lidar snow depth maps](https://nsidc.org/data/aso_3m_sd/versions/1)      
    - ASO lidar snow depth (target dataset)   
    - gaps in ASO data 

  - [Sentinel-1 radar backscatter](https://sentinel.esa.int/web/sentinel/missions/sentinel-1)
    - snow on Sentinel-1 VV polarization backscatter (snowon_vv)
    - snow on Sentinel-1 VH polarization backscatter (snowon_vh)
    - snow off Sentinel-1 VV polarization backscatter (snowoff_vv)
    - snow off Sentinel-1 VH polarization backscatter (snowoff_vh)
    - snow on Sentinel-1 VV polarization backscatter, 4-week mean (snowon_vv_mean)
    - snow on Sentinel-1 VH polarization backscatter, 4-week mean (snowon_vh_mean)
    - snow off Sentinel-1 VV polarization backscatter, 4-week mean (snowoff_vv_mean)
    - snow off Sentinel-1 VH polarization backscatter, 4-week mean (snowoff_vh_mean)
    - snow on Sentinel-1 cross ratio (VH - VV) (snowon_cr)
    - snow off Sentinel-1 cross ratio (VH - VV) (snowoff_cr)
    - change in cross ratio, snow on vs snow off (delta_cr)
    - gaps in Sentinel-1 data (rtc_gap_map)
    - gaps in Sentinel-1 mean data (rtc_mean_gap_map)

  - [Sentinel-2 multispectral imagery](https://sentinel.esa.int/web/sentinel/missions/sentinel-2)
    - aerosol optical thickness (snow on aerosol_optical_thickness)
    - coastal aerosol band (snow on coastal_aerosol)
    - blue band (snow on blue)
    - green band (snow on green)
    - red band (snow on red)
    - red edge 1 band (snow on red_edge1)
    - red edge 2 band (snow on red_edge2)
    - red edge 3 band (snow on red_edge3)
    - near infrared (NIR) band (snow on nir)
    - water vapor band (snow on water_vapor)
    - shortwave infrared 1 band (snow on swir1)
    - shortwave infrared 2 band (snow on swir2)
    - scene classification map (snow on scene_class_map)
    - water vapor product (snow on water_vapor_product)
    - Normalized Difference Vegetation Index (NDVI)
    - Normalized Difference Snow Index (NDSI)
    - Normalized Difference Water Index (NDWI)
    - gaps in Sentinel-2 data (s2_gap_map)

  - [PROBA-V Global Land Cover dataset](https://doi.org/10.5281/zenodo.3939050)
    - fractional forest cover (fcf)

  - [COP30 Digital Elevation Model](https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model)
    - elevation
    - slope
    - aspect
    - curvature
    - topographic position index (tpi)
    - terrain ruggedness index (tri)

  - Geospatial Information
    - latitude
    - longitude

  - Derived Variables
    - day of water year (dowy)

- Data Download and Raw Data Organization
  -  The assignment-required script is at [scripts/clean_data.py](scripts/clean_data.py) but please just look through the [notebooks/dataset_prep](notebooks/dataset_prep) notebooks as they're much better organized and easier to understand
- Basic Data Cleaning and Manipulation
  - [notebooks/Data_Cleaning.ipynb](notebooks/Data_Cleaning.ipynb) 
- Organizing Data into AI-Ready Format
  - [notebooks/Prepare_AI_Ready_Data.ipynb](notebooks/Prepare_AI_Ready_Data.ipynb)
  - note that the subdirectories in [data/ai_ready/](data/ai_ready/) are empty due to the ASO data sharing policies described above. the data are stored in .nc files so all of the training, testing, and validation files contain this sensitive data
  - See [notebooks/visualizations/ai-ready_viz.ipynb](notebooks/visualizations/ai-ready_viz.ipynb) for data visualizations which will serve as our [data/ai_ready/](data/ai_ready/) replacement
  - Note the the tile boundary polygons were generated in [QGIS](https://www.qgis.org/) using the "Create grid" tool.
  - The SNOTEL snow depth station data we want is easily queryable as seen at [notebooks/dataset_prep/snotel_exploration.ipynb](notebooks/dataset_prep/snotel_exploration.ipynb) (note that the interactive maps aren't viewable on the repo so see [this notebook](https://notebooksharing.space/view/0a79b23e82c6e156fddbcb0ae7e7727016099edce27b87b94757317c59b9d910#displayOptions=) instead). This data is not technically AI-ready, but a main goal of this project is to how best make this data AI-ready and how that ties in with encoding the point data most effectively.
- Exploratory Data Analysis (EDA)
  - [notebooks/EDA.ipynb](notebooks/EDA.ipynb)
- Dimensionality Discussion and Reduction
  - [notebooks/Dimensionality_Reduction.ipynb](notebooks/Dimensionality_Reduction.ipynb)

### Additional resources or background reading
- [spicy-snow background](https://github.com/SnowEx/spicy-snow/blob/main/contrib/brencher/tutorial/01background.ipynb)
- [spicy-snow paper](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-1018/egusphere-2024-1018.pdf)
- [Lievens et al. (2022) paper](https://tc.copernicus.org/articles/16/159/2022/) 
- [SAR basics](https://asf.alaska.edu/information/sar-information/what-is-sar/)
- [More SAR basics](https://www.earthdata.nasa.gov/learn/backgrounders/what-is-sar)
- [Sentinel-1 SAR](https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-1-sar)
- [More on ASF HyP3 RTC](https://hyp3-docs.asf.alaska.edu/guides/rtc_product_guide/)
- [SAR theory from 2022 UNAVCO InSAR class (more advanced)](https://nbviewer.org/github/parosen/Geo-SInC/blob/main/UNAVCO2022/0.8_SAR_Theory_Phenomenology/SAR.ipynb)
