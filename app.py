import streamlit as st
import pandas as pd
import osmnx as ox
import ee
from geopy.geocoders import Nominatim
import geemap.foliumap as geemap
import geemap.colormaps as cm
import folium
import pathlib
import json

service_account_info = st.secrets["ee"]

credentials = ee.ServiceAccountCredentials(
    email=service_account_info["client_email"],
    key_data=json.dumps(dict(service_account_info))  # Fix: convert AttrDict to JSON
)
try:
    ee.Initialize(credentials)
except Exception as e:
    st.error(f"Failed to initialize Earth Engine: {e}")
    st.stop()

def CalculateNDVI(image, sensor):
    bands = {'89': ['SR_B4', 'SR_B5'], '457': ['SR_B3', 'SR_B4'], 'S2': ['B8', 'B4']}
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    NDVI = image.normalizedDifference(bands[sensor]).rename('NDVI')
    return image.addBands(NDVI)

def CalculateNDBI(image, sensor):
    bands = {'89': ['SR_B6', 'SR_B5'], '457': ['SR_B5', 'SR_B4'], 'S2': ['B11', 'B8']}
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    NDBI = image.normalizedDifference(bands[sensor]).rename('NDBI')
    return image.addBands(NDBI)

def CalculateNDWI(image, sensor):
    bands = {'89': ['SR_B3', 'SR_B5'], '457': ['SR_B2', 'SR_B4'], 'S2': ['B3', 'B8']}
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    NDWI = image.normalizedDifference(bands[sensor]).rename('NDWI')
    return image.addBands(NDWI)

def CalculateNDSI(image, sensor):
    bands = {'89': ['SR_B3', 'SR_B6'], '457': ['SR_B2', 'SR_B5'], 'S2': ['B3', 'B11']}
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    NDSI = image.normalizedDifference(bands[sensor]).rename('NDSI')
    return image.addBands(NDSI)

def CalculateSAVI(image, sensor):
    bands = {
        '89': {'NIR': 'SR_B5', 'RED': 'SR_B4'},
        '457': {'NIR': 'SR_B4', 'RED': 'SR_B3'},
        'S2': {'NIR': 'B8', 'RED': 'B4'},
    }
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    SAVI = image.expression(
        '(NIR - RED) / (NIR + RED + 0.5) * 1.5',
        {
            'NIR': image.select(bands[sensor]['NIR']),
            'RED': image.select(bands[sensor]['RED']),
        }
    ).rename('SAVI')
    return image.addBands(SAVI)

def CalculateMSAVI(image, sensor):
    bands = {
        '89': {'NIR': 'SR_B5', 'RED': 'SR_B4'},
        '457': {'NIR': 'SR_B4', 'RED': 'SR_B3'},
        'S2': {'NIR': 'B8', 'RED': 'B4'},
    }
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    MSAVI = image.expression(
        '(2 * NIR + 1 - sqrt((2 * NIR + 1) ** 2 - 8 * (NIR - RED))) / 2',
        {
            'NIR': image.select(bands[sensor]['NIR']),
            'RED': image.select(bands[sensor]['RED']),
        }
    ).rename('MSAVI')
    return image.addBands(MSAVI)

def CalculateGNDVI(image, sensor):
    bands = {
        '89': {'NIR': 'SR_B5', 'GREEN': 'SR_B3'},
        '457': {'NIR': 'SR_B4', 'GREEN': 'SR_B2'},
        'S2': {'NIR': 'B8', 'GREEN': 'B3'},
    }
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    GNDVI = image.normalizedDifference([bands[sensor]['NIR'], bands[sensor]['GREEN']]).rename('GNDVI')
    return image.addBands(GNDVI)

def CalculateGCI(image, sensor):
    bands = {
        '89': {'NIR': 'SR_B5', 'GREEN': 'SR_B3'},
        '457': {'NIR': 'SR_B4', 'GREEN': 'SR_B2'},
        'S2': {'NIR': 'B8', 'GREEN': 'B3'},
    }
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    GCI = image.expression(
        '(NIR / GREEN) - 1',
        {
            'NIR': image.select(bands[sensor]['NIR']),
            'GREEN': image.select(bands[sensor]['GREEN']),
        }
    ).rename('GCI')
    return image.addBands(GCI)

def CalculateNBR(image, sensor):
    bands = {'89': ['SR_B5', 'SR_B7'], '457': ['SR_B4', 'SR_B7'], 'S2': ['B8', 'B12']}
    if sensor not in bands:
        raise ValueError(f"Unknown sensor type: {sensor}")
    NBR = image.normalizedDifference(bands[sensor]).rename('NBR')
    return image.addBands(NBR)

def CalculateEVI(image, sensor):
    bands = {
        '89': {'NIR': 'SR_B5', 'RED': 'SR_B4', 'BLUE': 'SR_B2'},
        '457': {'NIR': 'SR_B4', 'RED': 'SR_B3', 'BLUE': 'SR_B1'},
        'S2': {'NIR': 'B8', 'RED': 'B4', 'BLUE': 'B2'},
    }
    EVI = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
        {
            'NIR': image.select(bands[sensor]['NIR']),
            'RED': image.select(bands[sensor]['RED']),
            'BLUE': image.select(bands[sensor]['BLUE']),
        }
    ).rename('EVI')
    return image.addBands(EVI)

st.title("🛰️ Satellite / Index Viewer 🛰️")
start_date = st.text_input("Start Date", "2017-01-01")
end_date = st.text_input("End Date", "2025-01-01")

# Use OSMnx for location suggestions
def get_location_suggestions(query, max_results=5):
    if not query:
        return []
    try:
        gdf = ox.geocoder.geocode_to_gdf(query, which_result=None)
        if gdf.empty:
            return []
        return gdf['display_name'].tolist()[:max_results]
    except Exception as e:
        st.warning(f"Location suggestion error: {e}")
        return []

location_query = st.text_input("Location", "", placeholder="Enter a place name or coordinates")
location_suggestions = get_location_suggestions(location_query, 5) if location_query else []
selected_location = st.selectbox("Suggestions", location_suggestions) if location_suggestions else location_query

# Select options
index = st.selectbox("Index", ['RGB', 'NDVI', 'NDBI', 'NBR', 'EVI', 'NDSI', 'NDWI', 'SAVI', 'MSAVI', 'GNDVI', 'GCI'])
stretch_type = st.selectbox("Stretch", ["Default (-1, 1)", "Percentile (like Earth Engine)"])
stretch_min = st.number_input("Stretch Min Percentile", min_value=0, max_value=100, value=2, step=1)
stretch_max = st.number_input("Stretch Max Percentile", min_value=0, max_value=100, value=98, step=1)
cloud_percent = st.slider("Cloud %", 0, 100, 10)

satellite_options = {
    'Landsat 4': 'LANDSAT/LT04/C02/T1_L2',
    'Landsat 5': 'LANDSAT/LT05/C02/T1_L2',
    'Landsat 7': 'LANDSAT/LE07/C02/T1_L2',
    'Landsat 8': 'LANDSAT/LC08/C02/T1_L2',
    'Landsat 9': 'LANDSAT/LC09/C02/T1_L2',
    'Sentinel 2': 'COPERNICUS/S2_SR_HARMONIZED',
}
satellite_name = st.selectbox("Image", list(satellite_options.keys()))
satellite = satellite_options[satellite_name]

# Get ROI with OSMnx

def get_roi_geom(location_name):
    try:
        gdf = ox.geocoder.geocode_to_gdf(location_name)
        if gdf.empty:
            st.error("Location not found.")
            return None
        centroid = gdf.geometry.centroid.iloc[0]
        lon, lat = centroid.x, centroid.y
        return ee.Geometry.Point([lon, lat]).buffer(5000).bounds()
    except Exception as e:
        st.error(f"Error locating geometry: {e}")
        return None

# Visualization palettes/index functions (as in your original)
index_function = {'NDVI': CalculateNDVI, 'NDBI': CalculateNDBI, 'NBR': CalculateNBR, 'EVI': CalculateEVI, 'NDSI': CalculateNDSI, 'NDWI': CalculateNDWI, 'SAVI': CalculateSAVI, 'MSAVI': CalculateMSAVI, 'GNDVI': CalculateGNDVI, 'GCI': CalculateGCI}
index_palette = {
    'NDVI': ['#440154', '#3b528b', '#21908d', '#5dc962', '#fde725'],  # Viridis
    'NDBI': ['#f7fcf5', '#c7e9c0', '#74c476', '#238b45', '#00441b'],  # Green shades for built-up
    'NBR':  ['#d73027', '#f46d43', '#fdae61', '#fee08b', '#d9ef8b', '#a6d96a', '#1a9850'],  # Fire severity scale
    'EVI':  ['#ffffe5', '#f7fcb9', '#d9f0a3', '#addd8e', '#78c679', '#31a354', '#006837'],  # Vegetation
    'NDSI': ['#f7fcfd', '#e0ecf4', '#bfd3e6', '#9ebcda', '#8c96c6', '#8c6bb1', '#88419d', '#6e016b'],  # Snow index (Blue-Purple)
    'NDWI': ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#3182bd', '#08519c'],  # Water
    'SAVI': ['#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837'],  # Similar to NDVI but with less saturation
    'MSAVI': ['#f7fcf5', '#e5f5e0', '#c7e9c0', '#a1d99b', '#74c476', '#41ab5d', '#238b45'],  # Moderate SAVI
    'GNDVI': ['#ffffe5', '#f7fcb9', '#d9f0a3', '#addd8e', '#78c679', '#31a354', '#006837'],  # Green NDVI
    'GCI':   ['#edf8fb', '#b2e2e2', '#66c2a4', '#2ca25f', '#006d2c'],  # Chlorophyll concentration
}


# Display map

def show_map(roi_geom, start_date, end_date, satellite, index, cloud_percent):
    dataset = ee.ImageCollection(satellite).filterBounds(roi_geom).filterDate(start_date, end_date)
    try:
        if dataset.size().getInfo() == 0:
            st.error("No images found for the selected location and time range.")
            return

        sensor_type = 'S2' if 'S2' in satellite else ('89' if '89' in satellite else '457')
        dataset = dataset.filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE' if 'S2' in satellite else 'CLOUD_COVER', cloud_percent))

        rgb_vis_params = {'min': 0, 'max': 3000, 'bands': ['B4', 'B3', 'B2']} if 'S2' in satellite else {'min': 0, 'max': 3000, 'bands': ['SR_B4', 'SR_B3', 'SR_B2']}
        rgb_bands = ['B4', 'B3', 'B2'] if 'S2' in satellite else ['SR_B4', 'SR_B3', 'SR_B2']

        # --- Fix: Always select only the bands needed for the chosen index ---
        if index == 'RGB':
            dataset = dataset.select(rgb_bands)
            img = dataset.mosaic()
            img = img.select(rgb_bands)
            vis_params = rgb_vis_params
            index_name = 'RGB'
        else:
            # For indices, select only the bands needed for the calculation
            # Get all bands needed for the index
            if index in ['NDVI', 'NDBI', 'NBR', 'NDWI', 'NDSI']:
                # These use a pair of bands
                bands_needed = {
                    'NDVI': ['B8', 'B4'] if sensor_type == 'S2' else (['SR_B4', 'SR_B5'] if sensor_type == '89' else ['SR_B3', 'SR_B4']),
                    'NDBI': ['B11', 'B8'] if sensor_type == 'S2' else (['SR_B6', 'SR_B5'] if sensor_type == '89' else ['SR_B5', 'SR_B4']),
                    'NBR': ['B8', 'B12'] if sensor_type == 'S2' else (['SR_B5', 'SR_B7'] if sensor_type == '89' else ['SR_B4', 'SR_B7']),
                    'NDWI': ['B3', 'B8'] if sensor_type == 'S2' else (['SR_B3', 'SR_B5'] if sensor_type == '89' else ['SR_B2', 'SR_B4']),
                    'NDSI': ['B3', 'B11'] if sensor_type == 'S2' else (['SR_B3', 'SR_B6'] if sensor_type == '89' else ['SR_B2', 'SR_B5']),
                }[index]
                dataset = dataset.select(bands_needed)
            elif index in ['EVI', 'SAVI', 'MSAVI', 'GNDVI', 'GCI']:
                # These use more or different bands, so select all possible bands needed for all indices
                if sensor_type == 'S2':
                    bands_needed = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
                elif sensor_type == '89':
                    bands_needed = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
                else:
                    bands_needed = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']
                dataset = dataset.select(bands_needed)
            img = dataset.mosaic()
            img = index_function[index](img, sensor_type).select(index)
            vis_params = {'palette': index_palette[index]}
            stats_geom = roi_geom.buffer(2000)
            if stretch_type == "Percentile (like Earth Engine)":
                percentiles = img.reduceRegion(
                    reducer=ee.Reducer.percentile([stretch_min, stretch_max]),
                    geometry=stats_geom,
                    scale=30,
                    maxPixels=1e13
                )
                vis_params['min'] = percentiles.get(f"{index}_p{int(stretch_min)}").getInfo()
                vis_params['max'] = percentiles.get(f"{index}_p{int(stretch_max)}").getInfo()
            else:
                stats = img.reduceRegion(
                    reducer=ee.Reducer.minMax(),
                    geometry=stats_geom,
                    scale=30,
                    maxPixels=1e13
                )
                vis_params['min'] = stats.get(index + '_min').getInfo()
                vis_params['max'] = stats.get(index + '_max').getInfo()
            index_name = index

        m = geemap.Map()
        m.add_basemap('HYBRID')
        m.centerObject(roi_geom, 9)
        m.addLayer(img, vis_params, index_name)
        # --- Add HTML colorbar gradient for indices ---
        if index != 'RGB':
            # Build CSS gradient from palette
            palette = index_palette[index]
            gradient = f"linear-gradient(to top, {', '.join(palette)})"
            colorbar_html = f"""
            <div class="colorbar-panel" id="colorbar-panel" style="position: absolute; bottom: 35px; right: 20px; z-index: 99999; background-color: rgba(89, 99, 88, 0.80); padding: 10px; border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.2); font-size: 12px; text-align: center;">
                <b>{index_name} Value</b><br>
                <div style="display: flex; flex-direction: row; align-items: center; gap: 16px;">
                    <div class="colorbar-gradient" style="width: 20px; height: 120px; background: {gradient}; border-radius: 10px; margin-bottom: 5px;"></div>
                    <div class="colorbar-labels" style="display: flex; flex-direction: column; align-items: center; height: 120px; justify-content: space-between; margin-left: 8px;">
                        <span style="font-size: 12px;">{vis_params['max']:.2f}</span>
                        <span style="font-size: 12px;">{vis_params['min']:.2f}</span>
                    </div>
                </div>
            </div>
            """
            m.get_root().html.add_child(folium.Element(colorbar_html))
        st.components.v1.html(m.to_html(), height=600, width=800)
    except Exception as e:
        st.error(f"Map error: {e}")

# Button trigger
if st.button("Update Map"):
    if not selected_location:
        st.warning("Please enter a location.")
    else:
        roi_geom = get_roi_geom(selected_location)
        if roi_geom:
            show_map(roi_geom, start_date, end_date, satellite, index, cloud_percent)
