import streamlit as st
import pandas as pd
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

def get_location_suggestions(query, max_results=5):
    if not query:
        return []
    geolocator = Nominatim(user_agent="MyGeocoder")
    try:
        locations = geolocator.geocode(query, exactly_one=False, limit=max_results, addressdetails=True)
        if not locations:
            return []
        suggestions = []
        for loc in locations:
            display_name = loc.address
            suggestions.append(display_name)
        return suggestions
    except Exception:
        return []

# Location search with suggestions
location_query = st.text_input("Location", "", placeholder="Enter a place name or coordinates")
location_suggestions = get_location_suggestions(location_query, 5) if location_query else []
if location_suggestions:
    selected_location = st.selectbox("Suggestions", location_suggestions)
else:
    selected_location = location_query

index = st.selectbox("Index", ['RGB', 'NDVI', 'NDBI', 'NBR', 'EVI'])
# Add stretch option
stretch_type = st.selectbox("Stretch", ["Default (-1, 1)", "Percentile (like Earth Engine)"])
stretch_min = st.number_input("Stretch Min Percentile", min_value=0, max_value=100, value=2, step=1)
stretch_max = st.number_input("Stretch Max Percentile", min_value=0, max_value=100, value=98, step=1)
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
cloud_percent = st.slider("Cloud %", 0, 100, 10)

index_function = {'NDVI': CalculateNDVI, 'NDBI': CalculateNDBI, 'NBR': CalculateNBR, 'EVI': CalculateEVI}
index_palette = {
    'NDVI': ['#654321', '#D2B48C', '#FFFF00', '#7FFF00', '#006400'],
    'NDBI': ['#FFFFFF', '#D3D3D3', '#A9A9A9', '#696969', '#000000'],
    'NBR': ['#006400', '#7CFC00', '#FFBF00', '#FF0000', '#8B0000'],
    'EVI': ['#ADD8E6', '#00FFFF', '#008080', '#00FF00', '#006400'],
}

def get_roi_geom(location):
    geolocator = Nominatim(user_agent="MyGeocoder")
    loc = geolocator.geocode(location)
    if not loc:
        st.error("Invalid location. Please enter a valid place.")
        return None
    return ee.Geometry.Point([loc.longitude, loc.latitude])

def show_map(roi_geom, start_date, end_date, satellite, index, cloud_percent):
    dataset = ee.ImageCollection(satellite).filterBounds(roi_geom).filterDate(str(start_date), str(end_date))
    try:
        if dataset.size().getInfo() > 0:
            sensor_type = 'S2' if 'S2' in satellite else ('89' if '89' in satellite else '457')
            if 'S2' in satellite:
                dataset = dataset.filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cloud_percent))
                rgb_vis_params = {'min': 0, 'max': 3000, 'bands': ['B4', 'B3', 'B2']}
            else:
                dataset = dataset.filter(ee.Filter.lte('CLOUD_COVER', cloud_percent))
                rgb_vis_params = {'min': 0, 'max': 3000, 'bands': ['SR_B4', 'SR_B3', 'SR_B2']}

            img = dataset.first()
            if img is None or img.getInfo() is None:
                st.error("No images found for the selected location and time range.")
                return

            if index == 'RGB':
                img = img.select(rgb_vis_params['bands'])
                vis_params = rgb_vis_params
                index_name = 'RGB'
            else:
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
                    min_val = percentiles.get(f"{index}_p{int(stretch_min)}").getInfo()
                    max_val = percentiles.get(f"{index}_p{int(stretch_max)}").getInfo()
                else:
                    stats = img.reduceRegion(
                        reducer=ee.Reducer.minMax(),
                        geometry=stats_geom,
                        scale=30,
                        maxPixels=1e13
                    )
                    min_val = stats.get(index + '_min').getInfo()
                    max_val = stats.get(index + '_max').getInfo()

                vis_params['min'] = min_val
                vis_params['max'] = max_val
                index_name = index

            # Display map
            m = geemap.Map()
            m.add_basemap('HYBRID')
            m.centerObject(roi_geom, 9)
            m.addLayer(img, vis_params, index_name)
            st.components.v1.html(m.to_html(), height=600, width=800)

            # Download GeoTIFF if index
            if index != 'RGB':
                try:
                    clipped = img.clip(roi_geom)
                    download_params = {
                        'scale': 30 if 'LANDSAT' in satellite else 10,
                        'region': roi_geom.buffer(1000).bounds().getInfo()['coordinates'],
                        'fileFormat': 'GeoTIFF',
                        'formatOptions': {'cloudOptimized': True}
                    }
                    url = clipped.getDownloadURL(download_params)

                    if st.button(f"\U0001F4E5 Generate {index} GeoTIFF"):
                        st.success("Download link is ready:")
                        st.markdown(f"[Click here to download {index} as GeoTIFF]({url})", unsafe_allow_html=True)

                except Exception as e:
                    st.error(f"Failed to generate GeoTIFF download: {e}")

        else:
            st.error("No images found for the selected location and time range.")

    except Exception as e:
        st.error(f"Error processing image: {e}")

if st.button("Update Map"):
    if not selected_location:
        st.warning("Please enter a location.")
    else:
        roi_geom = get_roi_geom(selected_location)
        if roi_geom:
            show_map(roi_geom, start_date, end_date, satellite, index, cloud_percent)