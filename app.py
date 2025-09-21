"""
Satellite Index Viewer - A Streamlit application for viewing satellite imagery and indices
Enhanced with RGB band selection and percentile-based min/max controls
"""

import streamlit as st
import geopandas as gpd
import osmnx as ox
import ee
from geopy.geocoders import Nominatim
import geemap.foliumap as geemap
import folium
import json
import pandas as pd

# ============================================================================
# PAGE CONFIGURATION - MUST BE FIRST
# ============================================================================

st.set_page_config(
    page_title="🛰️ Satellite Index Visualizer 🛰️",
    page_icon="🌎",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# CONFIGURATION AND INITIALIZATION
# ============================================================================


def initialize_earth_engine():
    """Initialize Google Earth Engine with service account credentials"""
    try:
        service_account_info = st.secrets["ee"]
        credentials = ee.ServiceAccountCredentials(
            email=service_account_info["client_email"],
            key_data=json.dumps(dict(service_account_info))
        )
        ee.Initialize(credentials)
    except Exception as e:
        st.error(f"Failed to initialize Earth Engine: {e}")
        st.stop()


# Satellite datasets configuration
SATELLITE_OPTIONS = {
    'Landsat 4': 'LANDSAT/LT04/C02/T1_L2',
    'Landsat 5': 'LANDSAT/LT05/C02/T1_L2',
    'Landsat 7': 'LANDSAT/LE07/C02/T1_L2',
    'Landsat 8': 'LANDSAT/LC08/C02/T1_L2',
    'Landsat 9': 'LANDSAT/LC09/C02/T1_L2',
    'Sentinel 2': 'COPERNICUS/S2_SR_HARMONIZED',
}

# Band information for each satellite
SATELLITE_BANDS = {
    'LANDSAT/LT04/C02/T1_L2': {
        'SR_B1': 'Blue',
        'SR_B2': 'Green',
        'SR_B3': 'Red',
        'SR_B4': 'NIR',
        'SR_B5': 'SWIR1',
        'SR_B7': 'SWIR2'
    },
    'LANDSAT/LT05/C02/T1_L2': {
        'SR_B1': 'Blue',
        'SR_B2': 'Green',
        'SR_B3': 'Red',
        'SR_B4': 'NIR',
        'SR_B5': 'SWIR1',
        'SR_B7': 'SWIR2'
    },
    'LANDSAT/LE07/C02/T1_L2': {
        'SR_B1': 'Blue',
        'SR_B2': 'Green',
        'SR_B3': 'Red',
        'SR_B4': 'NIR',
        'SR_B5': 'SWIR1',
        'SR_B7': 'SWIR2'
    },
    'LANDSAT/LC08/C02/T1_L2': {
        'SR_B2': 'Blue',
        'SR_B3': 'Green',
        'SR_B4': 'Red',
        'SR_B5': 'NIR',
        'SR_B6': 'SWIR1',
        'SR_B7': 'SWIR2'
    },
    'LANDSAT/LC09/C02/T1_L2': {
        'SR_B2': 'Blue',
        'SR_B3': 'Green',
        'SR_B4': 'Red',
        'SR_B5': 'NIR',
        'SR_B6': 'SWIR1',
        'SR_B7': 'SWIR2'
    },
    'COPERNICUS/S2_SR_HARMONIZED': {
        'B2': 'Blue',
        'B3': 'Green',
        'B4': 'Red',
        'B8': 'NIR',
        'B11': 'SWIR1',
        'B12': 'SWIR2'
    }
}

# Index visualization palettes
INDEX_PALETTES = {
    'NDVI': ['#d73027', '#fdae61', '#fee08b', '#d9ef8b', '#91cf60', '#1a9850'],
    'NDBI': ['#f7f7f7', '#cccccc', '#969696', '#525252', '#252525'],
    'NBR':  ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee08b', '#d9ef8b', '#66bd63', '#1a9850', '#006837'],
    'EVI':  ['#ffffe5', '#f7fcb9', '#d9f0a3', '#addd8e', '#78c679', '#31a354', '#006837'],
    'NDSI': ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#3182bd', '#08519c'],
    'NDWI': ['#f7fcfd', '#e0ecf4', '#9ecae1', '#3182bd', '#08519c'],
    'SAVI': ['#f7f4f9', '#d0d1e6', '#a6bddb', '#74a9cf', '#2b8cbe', '#045a8d'],
    'MSAVI': ['#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837'],
    'GNDVI': ['#edf8e9', '#c7e9c0', '#a1d99b', '#74c476', '#31a354', '#006d2c'],
    'GCI':   ['#ffffe5', '#f7fcb9', '#d9f0a3', '#78c679', '#238443', '#004529']
}

# Default value ranges for indices
INDEX_DEFAULTS = {
    'NDVI': (-1, 1),
    'NDBI': (-1, 1),
    'NBR': (-1, 1),
    'NDWI': (-1, 1),
    'NDSI': (-1, 1),
    'EVI': (-1, 2),
    'SAVI': (-1, 1),
    'MSAVI': (-1, 1),
    'GNDVI': (-1, 1),
    'GCI': (0, 10)
}

# ============================================================================
# INDEX CALCULATION FUNCTIONS
# ============================================================================


def get_sensor_bands():
    """Return band mappings for different sensors"""
    return {
        'NIR_RED': {
            '89': ['SR_B5', 'SR_B4'],
            '457': ['SR_B4', 'SR_B3'],
            'S2': ['B8', 'B4']
        },
        'SWIR_NIR': {
            '89': ['SR_B6', 'SR_B5'],
            '457': ['SR_B5', 'SR_B4'],
            'S2': ['B11', 'B8']
        },
        'GREEN_NIR': {
            '89': ['SR_B3', 'SR_B5'],
            '457': ['SR_B2', 'SR_B4'],
            'S2': ['B3', 'B8']
        },
        'GREEN_SWIR': {
            '89': ['SR_B3', 'SR_B6'],
            '457': ['SR_B2', 'SR_B5'],
            'S2': ['B3', 'B11']
        },
        'NIR_SWIR2': {
            '89': ['SR_B5', 'SR_B7'],
            '457': ['SR_B4', 'SR_B7'],
            'S2': ['B8', 'B12']
        },
        'RGB': {
            '89': {'NIR': 'SR_B5', 'RED': 'SR_B4', 'BLUE': 'SR_B2'},
            '457': {'NIR': 'SR_B4', 'RED': 'SR_B3', 'BLUE': 'SR_B1'},
            'S2': {'NIR': 'B8', 'RED': 'B4', 'BLUE': 'B2'},
        },
        'RGB_WITH_GREEN': {
            '89': {'NIR': 'SR_B5', 'RED': 'SR_B4', 'GREEN': 'SR_B3'},
            '457': {'NIR': 'SR_B4', 'RED': 'SR_B3', 'GREEN': 'SR_B2'},
            'S2': {'NIR': 'B8', 'RED': 'B4', 'GREEN': 'B3'},
        }
    }


def calculate_normalized_difference(image, bands, index_name):
    """Calculate normalized difference index"""
    return image.normalizedDifference(bands).rename(index_name)


def calculate_ndvi(image, sensor):
    """Calculate Normalized Difference Vegetation Index"""
    bands = get_sensor_bands()['NIR_RED'][sensor]
    ndvi = calculate_normalized_difference(image, bands, 'NDVI')
    return image.addBands(ndvi)


def calculate_ndbi(image, sensor):
    """Calculate Normalized Difference Built-up Index"""
    bands = get_sensor_bands()['SWIR_NIR'][sensor]
    ndbi = calculate_normalized_difference(image, bands, 'NDBI')
    return image.addBands(ndbi)


def calculate_ndwi(image, sensor):
    """Calculate Normalized Difference Water Index"""
    bands = get_sensor_bands()['GREEN_NIR'][sensor]
    ndwi = calculate_normalized_difference(image, bands, 'NDWI')
    return image.addBands(ndwi)


def calculate_ndsi(image, sensor):
    """Calculate Normalized Difference Snow Index"""
    bands = get_sensor_bands()['GREEN_SWIR'][sensor]
    ndsi = calculate_normalized_difference(image, bands, 'NDSI')
    return image.addBands(ndsi)


def calculate_nbr(image, sensor):
    """Calculate Normalized Burn Ratio"""
    bands = get_sensor_bands()['NIR_SWIR2'][sensor]
    nbr = calculate_normalized_difference(image, bands, 'NBR')
    return image.addBands(nbr)


def calculate_gndvi(image, sensor):
    """Calculate Green Normalized Difference Vegetation Index"""
    bands_map = get_sensor_bands()['RGB_WITH_GREEN'][sensor]
    gndvi = image.normalizedDifference(
        [bands_map['NIR'], bands_map['GREEN']]).rename('GNDVI')
    return image.addBands(gndvi)


def calculate_savi(image, sensor):
    """Calculate Soil Adjusted Vegetation Index"""
    bands_map = get_sensor_bands()['RGB_WITH_GREEN'][sensor]
    savi = image.expression(
        '(NIR - RED) / (NIR + RED + 0.5) * 1.5',
        {
            'NIR': image.select(bands_map['NIR']),
            'RED': image.select(bands_map['RED']),
        }
    ).rename('SAVI')
    return image.addBands(savi)


def calculate_msavi(image, sensor):
    """Calculate Modified Soil Adjusted Vegetation Index"""
    bands_map = get_sensor_bands()['RGB_WITH_GREEN'][sensor]
    msavi = image.expression(
        '(2 * NIR + 1 - sqrt((2 * NIR + 1) ** 2 - 8 * (NIR - RED))) / 2',
        {
            'NIR': image.select(bands_map['NIR']),
            'RED': image.select(bands_map['RED']),
        }
    ).rename('MSAVI')
    return image.addBands(msavi)


def calculate_gci(image, sensor):
    """Calculate Green Chlorophyll Index"""
    bands_map = get_sensor_bands()['RGB_WITH_GREEN'][sensor]
    gci = image.expression(
        '(NIR / GREEN) - 1',
        {
            'NIR': image.select(bands_map['NIR']),
            'GREEN': image.select(bands_map['GREEN']),
        }
    ).rename('GCI')
    return image.addBands(gci)


def calculate_evi(image, sensor):
    """Calculate Enhanced Vegetation Index"""
    bands_map = get_sensor_bands()['RGB'][sensor]
    evi = image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
        {
            'NIR': image.select(bands_map['NIR']),
            'RED': image.select(bands_map['RED']),
            'BLUE': image.select(bands_map['BLUE']),
        }
    ).rename('EVI')
    return image.addBands(evi)


# Index function mapping
INDEX_FUNCTIONS = {
    'NDVI': calculate_ndvi,
    'NDBI': calculate_ndbi,
    'NBR': calculate_nbr,
    'EVI': calculate_evi,
    'NDSI': calculate_ndsi,
    'NDWI': calculate_ndwi,
    'SAVI': calculate_savi,
    'MSAVI': calculate_msavi,
    'GNDVI': calculate_gndvi,
    'GCI': calculate_gci
}

# ============================================================================
# LOCATION AND GEOMETRY UTILITIES
# ============================================================================


def get_location_suggestions(query, limit=5):
    """Get location suggestions using Nominatim geocoder"""
    if not query:
        return []

    try:
        geolocator = Nominatim(user_agent="sat-index-viewer")
        hits = geolocator.geocode(
            query, exactly_one=False, limit=limit, addressdetails=True
        )
        if not hits:
            return []
        return [{"name": h.address, "lat": h.latitude, "lon": h.longitude} for h in hits]
    except Exception as e:
        st.warning(f"Suggestion lookup error: {e}")
        return []


def get_boundary_geometry(name):
    """Get polygon/multipolygon boundary as a GeoDataFrame"""
    try:
        # Special handling for Algerian wilayas
        if "wilaya" in name.lower():
            gdf = ox.geocode_to_gdf({
                "state": name.replace("wilaya", "").strip(),
                "country": "Algeria"
            })
        else:
            gdf = ox.geocode_to_gdf(name, which_result=None)

        if gdf is None or gdf.empty:
            return None

        gdf = gdf.to_crs(epsg=4326)

        # Prefer polygon/multipolygon
        polys = gdf[gdf.geometry.geom_type.isin(["Polygon", "MultiPolygon"])]
        if not polys.empty:
            return polys.iloc[[0]][["geometry"]]

        # If only point -> create buffer
        pt_rows = gdf[gdf.geometry.geom_type == "Point"]
        if not pt_rows.empty:
            pt = pt_rows.geometry.iloc[0]
            rough = gpd.GeoDataFrame(
                geometry=[pt.buffer(0.02).envelope], crs="EPSG:4326"
            )
            return rough

        return gdf.iloc[[0]][["geometry"]]

    except Exception as e:
        st.error(f"Error locating geometry: {e}")
        return None

# ============================================================================
# RGB BAND SELECTION AND PERCENTILE FUNCTIONS
# ============================================================================


def get_band_statistics(image, geometry, bands, percentiles, scale=30):
    """Calculate percentile statistics for selected bands"""
    try:
        stats = image.select(bands).reduceRegion(
            reducer=ee.Reducer.percentile(percentiles),
            geometry=geometry,
            scale=scale,
            maxPixels=1e13,
            bestEffort=True,
            tileScale=2
        ).getInfo()

        return stats
    except Exception as e:
        st.warning(f"Error calculating band statistics: {e}")
        return {}


def create_rgb_band_selector(satellite, available_bands):
    """Create RGB band selection interface"""
    st.sidebar.subheader("🎨 RGB Band Selection")

    band_options = list(available_bands.keys())
    band_labels = [
        f"{band} ({available_bands[band]})" for band in band_options]

    # Default RGB combinations
    default_rgb = get_default_rgb_bands(satellite)

    # Find indices for default bands
    try:
        red_idx = band_options.index(default_rgb['red'])
        green_idx = band_options.index(default_rgb['green'])
        blue_idx = band_options.index(default_rgb['blue'])
    except (ValueError, KeyError):
        red_idx, green_idx, blue_idx = 0, 1, 2

    col1, col2, col3 = st.sidebar.columns(3)
    with col1:
        red_band = st.selectbox("Red", band_options, index=red_idx,
                                format_func=lambda x: f"{x}\n({available_bands[x]})")
    with col2:
        green_band = st.selectbox("Green", band_options, index=green_idx,
                                  format_func=lambda x: f"{x}\n({available_bands[x]})")
    with col3:
        blue_band = st.selectbox("Blue", band_options, index=blue_idx,
                                 format_func=lambda x: f"{x}\n({available_bands[x]})")

    return [red_band, green_band, blue_band]


def get_default_rgb_bands(satellite):
    """Get default RGB band combination for satellite"""
    if 'S2' in satellite:
        return {'red': 'B4', 'green': 'B3', 'blue': 'B2'}
    elif any(x in satellite for x in ['LC08', 'LC09']):
        return {'red': 'SR_B4', 'green': 'SR_B3', 'blue': 'SR_B2'}
    else:  # Landsat 4,5,7
        return {'red': 'SR_B3', 'green': 'SR_B2', 'blue': 'SR_B1'}


def create_percentile_controls():
    """Create percentile-based min/max controls"""
    st.sidebar.subheader("📊 Percentile Stretch")

    stretch_method = st.sidebar.selectbox(
        "Stretch Method",
        ["Percentile", "Manual", "Standard Deviation"]
    )

    if stretch_method == "Percentile":
        col1, col2 = st.sidebar.columns(2)
        with col1:
            min_percentile = st.number_input("Min %", 0.0, 50.0, 2.0, 0.1)
        with col2:
            max_percentile = st.number_input("Max %", 50.0, 100.0, 98.0, 0.1)
        return stretch_method, min_percentile, max_percentile, None, None

    elif stretch_method == "Manual":
        col1, col2 = st.sidebar.columns(2)
        with col1:
            min_val = st.number_input("Min Value", 0, 10000, 0, 50)
        with col2:
            max_val = st.number_input("Max Value", 0, 10000, 3000, 50)
        return stretch_method, None, None, min_val, max_val

    else:  # Standard Deviation
        std_devs = st.sidebar.slider("Std Deviations", 1.0, 4.0, 2.0, 0.1)
        return stretch_method, None, None, std_devs, None


# ============================================================================
# IMAGE PROCESSING AND VISUALIZATION
# ============================================================================


def get_required_bands(index, sensor_type):
    """Get the bands needed for a specific index calculation"""
    band_requirements = {
        'NDVI': get_sensor_bands()['NIR_RED'][sensor_type],
        'NDBI': get_sensor_bands()['SWIR_NIR'][sensor_type],
        'NBR': get_sensor_bands()['NIR_SWIR2'][sensor_type],
        'NDWI': get_sensor_bands()['GREEN_NIR'][sensor_type],
        'NDSI': get_sensor_bands()['GREEN_SWIR'][sensor_type],
    }

    # For complex indices, return all possible bands
    complex_indices = ['EVI', 'SAVI', 'MSAVI', 'GNDVI', 'GCI']
    if index in complex_indices:
        if sensor_type == 'S2':
            return ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
        elif sensor_type == '89':
            return ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
        else:
            return ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']

    return band_requirements.get(index, [])


def calculate_image_statistics(image, geometry, scale, stretch_type, stretch_min, stretch_max, index):
    """Calculate min/max statistics for image visualization"""
    try:
        # --- Step 1: Try absolute min/max first ---
        stats = image.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=geometry,
            scale=scale,
            maxPixels=1e13,
            bestEffort=True,
            tileScale=2
        ).getInfo()

        min_val = stats.get(f"{index}_min") or stats.get("min")
        max_val = stats.get(f"{index}_max") or stats.get("max")

        # Fallback: look for any 'min'/'max' keys
        if min_val is None or max_val is None:
            min_keys = [k for k in stats.keys() if "min" in k.lower()]
            max_keys = [k for k in stats.keys() if "max" in k.lower()]
            if min_keys and max_keys:
                min_val = stats[min_keys[0]]
                max_val = stats[max_keys[0]]

        # --- Step 2: If min/max invalid or too narrow, try percentiles ---
        if (
            min_val is None or max_val is None
            or (max_val - min_val) < 0.2  # suspiciously narrow range
            or stretch_type == "Percentile (like Earth Engine)"
        ):
            stats = image.reduceRegion(
                reducer=ee.Reducer.percentile([stretch_min, stretch_max]),
                geometry=geometry,
                scale=scale,
                maxPixels=1e13,
                bestEffort=True,
                tileScale=2
            ).getInfo()

            min_key = f"{index}_p{int(stretch_min)}"
            max_key = f"{index}_p{int(stretch_max)}"
            min_val = stats.get(min_key) or stats.get(f"p{int(stretch_min)}")
            max_val = stats.get(max_key) or stats.get(f"p{int(stretch_max)}")

        # --- Step 3: Final fallback to defaults ---
        if min_val is None or max_val is None:
            min_val, max_val = INDEX_DEFAULTS.get(index, (-1, 1))

        st.info(f"Retrieved {index} min={min_val:.2f}, max={max_val:.2f}")

        return float(min_val), float(max_val)

    except Exception as e:
        print(f"Statistics calculation failed: {e}")
        return INDEX_DEFAULTS.get(index, (-1, 1))


def calculate_rgb_stretch_values(image, geometry, rgb_bands, stretch_method,
                                 min_percentile, max_percentile,
                                 min_val=None, max_val=None, scale=30):
    """Calculate stretch values for RGB visualization using percentiles"""
    try:
        if stretch_method == "Manual":
            return min_val, max_val

        # Use percentile method for RGB bands
        stats = get_band_statistics(image, geometry, rgb_bands,
                                    [min_percentile, max_percentile], scale)

        # Extract percentile values for each band and take average
        min_vals = []
        max_vals = []

        for band in rgb_bands:
            min_key = f"{band}_p{min_percentile:.1f}".replace('.', '_')
            max_key = f"{band}_p{max_percentile:.1f}".replace('.', '_')

            # Try different key formats
            possible_min_keys = [f"{band}_p{min_percentile:.1f}".replace('.', '_'),
                                 f"{band}_p{int(min_percentile)}"]
            possible_max_keys = [f"{band}_p{max_percentile:.1f}".replace('.', '_'),
                                 f"{band}_p{int(max_percentile)}"]

            min_found = False
            max_found = False

            for key in possible_min_keys:
                if key in stats:
                    min_vals.append(stats[key])
                    min_found = True
                    break

            for key in possible_max_keys:
                if key in stats:
                    max_vals.append(stats[key])
                    max_found = True
                    break

            if not min_found or not max_found:
                st.warning(f"Could not find percentile stats for band {band}")

        if min_vals and max_vals:
            return int(sum(min_vals) / len(min_vals)), int(sum(max_vals) / len(max_vals))
        else:
            return 0, 3000

    except Exception as e:
        st.warning(f"Error calculating RGB stretch values: {e}")
        return 0, 3000


def calculate_index_stretch_values(image, geometry, index_name, stretch_method,
                                   min_percentile, max_percentile,
                                   min_val=None, max_val=None, scale=30):
    """Calculate stretch values for index visualization using percentiles"""
    try:
        if stretch_method == "Manual":
            return min_val, max_val

        # Use percentile method for indices
        stats = image.reduceRegion(
            reducer=ee.Reducer.percentile([min_percentile, max_percentile]),
            geometry=geometry,
            scale=scale,
            maxPixels=1e13,
            bestEffort=True,
            tileScale=2
        ).getInfo()

        # Try different key formats for the index
        possible_min_keys = [
            f"{index_name}_p{min_percentile:.1f}".replace('.', '_'),
            f"{index_name}_p{int(min_percentile)}",
            f"p{min_percentile:.1f}".replace('.', '_'),
            f"p{int(min_percentile)}"
        ]

        possible_max_keys = [
            f"{index_name}_p{max_percentile:.1f}".replace('.', '_'),
            f"{index_name}_p{int(max_percentile)}",
            f"p{max_percentile:.1f}".replace('.', '_'),
            f"p{int(max_percentile)}"
        ]

        min_val_found = None
        max_val_found = None

        for key in possible_min_keys:
            if key in stats:
                min_val_found = stats[key]
                break

        for key in possible_max_keys:
            if key in stats:
                max_val_found = stats[key]
                break

        if min_val_found is not None and max_val_found is not None:
            return float(min_val_found), float(max_val_found)
        else:
            # Fallback to default ranges
            return INDEX_DEFAULTS.get(index_name, (-1, 1))

    except Exception as e:
        st.warning(f"Error calculating index stretch values: {e}")
        return INDEX_DEFAULTS.get(index_name, (-1, 1))


def create_colorbar_html(index_name, min_val, max_val, palette):
    """Create HTML for the colorbar visualization"""
    gradient = f"linear-gradient(to top, {', '.join(palette)})"

    def format_val(val):
        try:
            return f"{float(val):.2f}"
        except:
            return "N/A"

    min_str = format_val(min_val)
    max_str = format_val(max_val)

    return f"""
    <div class="colorbar-panel" style="position: absolute; bottom: 35px; right: 20px; z-index: 99999; background-color: rgb(255, 255, 255); padding: 10px; border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.2); font-size: 12px; text-align: center;">
        <b>{index_name} Value</b><br>
        <div style="display: flex; flex-direction: row; align-items: center; gap: 16px;">
            <div class="colorbar-gradient" style="width: 20px; height: 120px; background: {gradient}; border-radius: 10px; margin-bottom: 5px;"></div>
            <div class="colorbar-labels" style="display: flex; flex-direction: column; align-items: center; height: 120px; justify-content: space-between; margin-left: 8px;">
                <span style="font-size: 12px;">{max_str}</span>
                <span style="font-size: 12px;">{min_str}</span>
            </div>
        </div>
    </div>
    """


def show_map(roi_geom, start_date, end_date, satellite, index, cloud_percent,
             stretch_method, min_percentile, max_percentile, min_val=None, max_val=None,
             rgb_bands=None):
    """Main function to process and display satellite imagery with unified percentile controls"""
    try:
        # Get dataset and filter
        dataset = ee.ImageCollection(
            satellite).filterDate(start_date, end_date)

        if dataset.size().getInfo() == 0:
            st.error("No images found for the selected location and time range.")
            return

        # Determine sensor type
        sensor_type = 'S2' if 'S2' in satellite else (
            '89' if '89' in satellite else '457')
        scale = 10 if sensor_type == 'S2' else 30

        # Filter by cloud cover
        cloud_property = 'CLOUDY_PIXEL_PERCENTAGE' if 'S2' in satellite else 'CLOUD_COVER'
        dataset = dataset.filter(ee.Filter.lte(cloud_property, cloud_percent))

        # Process RGB or Index
        if index == 'RGB':
            # Use custom RGB bands if provided, otherwise use defaults
            if rgb_bands:
                selected_bands = rgb_bands
            else:
                default_rgb = get_default_rgb_bands(satellite)
                selected_bands = [default_rgb['red'],
                                  default_rgb['green'], default_rgb['blue']]

            img = dataset.select(selected_bands).median().clip(roi_geom)

            # Calculate stretch values using percentiles
            min_stretch, max_stretch = calculate_rgb_stretch_values(
                img, roi_geom, selected_bands, stretch_method,
                min_percentile, max_percentile, min_val, max_val, scale
            )

            vis_params = {
                'min': min_stretch,
                'max': max_stretch,
                'bands': selected_bands
            }

            st.info(
                f"RGB stretch ({stretch_method}): {min_stretch} - {max_stretch}")
            index_name = f'RGB ({"-".join(selected_bands)})'

        else:
            # Select required bands and process index
            required_bands = get_required_bands(index, sensor_type)
            dataset = dataset.select(required_bands)
            img = dataset.median().clip(roi_geom)

            # Calculate index
            index_img = INDEX_FUNCTIONS[index](img, sensor_type)

            # Select the index band
            try:
                index_img = index_img.select(index)
            except:
                index_img = index_img.select([0]).rename(index)

            # Calculate stretch values using percentiles
            min_stretch, max_stretch = calculate_index_stretch_values(
                index_img, roi_geom, index, stretch_method,
                min_percentile, max_percentile, min_val, max_val, scale
            )

            vis_params = {
                'palette': INDEX_PALETTES[index],
                'min': min_stretch,
                'max': max_stretch
            }

            st.info(
                f"{index} stretch ({stretch_method}): {min_stretch:.3f} - {max_stretch:.3f}")
            img = index_img
            index_name = index

        # Create and display map
        m = geemap.Map()
        m.add_basemap('HYBRID')
        m.centerObject(roi_geom, 12)
        m.addLayer(img, vis_params, index_name)

        # Add colorbar for indices (not for RGB)
        if index != 'RGB':
            colorbar_html = create_colorbar_html(
                index_name, vis_params['min'], vis_params['max'], INDEX_PALETTES[index]
            )
            m.get_root().html.add_child(folium.Element(colorbar_html))

        # Display in Streamlit
        st.components.v1.html(m.to_html(), height=600, width=1400)

        return img, vis_params

    except Exception as e:
        st.error(f"Map error: {e}")
        return None, None

# ============================================================================
# STREAMLIT USER INTERFACE
# ============================================================================

def create_sidebar_controls():
    """Create sidebar controls with unified percentile approach for both RGB and indices"""
    st.sidebar.header("🛰️ Satellite Parameters")

    # Date inputs
    start_date = st.sidebar.date_input(
        "Start Date", value=pd.to_datetime("2020-01-01")).strftime("%Y-%m-%d")
    end_date = st.sidebar.date_input(
        "End Date", value=pd.to_datetime("2024-01-01")).strftime("%Y-%m-%d")

    # Satellite selection
    satellite_name = st.sidebar.selectbox(
        "Satellite", list(SATELLITE_OPTIONS.keys()))
    satellite = SATELLITE_OPTIONS[satellite_name]

    # Index selection
    index = st.sidebar.selectbox(
        "Index/Band", ['RGB'] + list(INDEX_FUNCTIONS.keys()))

    # Cloud cover
    cloud_percent = st.sidebar.slider("Max Cloud Cover (%)", 0, 100, 10)

    # RGB band selector (always show for RGB mode)
    rgb_bands = None
    if index == 'RGB':
        available_bands = SATELLITE_BANDS[satellite]
        rgb_bands = create_rgb_band_selector(satellite, available_bands)

    # Unified percentile controls for both RGB and indices
    stretch_method, min_percentile, max_percentile, min_val, max_val = create_percentile_controls()

    return (start_date, end_date, satellite, index, cloud_percent,
            rgb_bands, stretch_method, min_percentile, max_percentile, min_val, max_val)


def create_location_input():
    """Create location input and suggestion interface"""
    st.header("📍 Location Selection")

    col1, col2 = st.columns([3, 1])
    with col1:
        location_query = st.text_input(
            "Search for a location:", placeholder="Search for your location")

    selected_geometry = None
    selected_name = None

    if location_query:
        suggestions = get_location_suggestions(location_query)
        if suggestions:
            suggestion_names = [s['name'] for s in suggestions]
            selected_idx = st.selectbox("Select from suggestions:", range(len(suggestion_names)),
                                        format_func=lambda i: suggestion_names[i])

            if selected_idx is not None:
                selected_name = suggestions[selected_idx]['name']
                gdf = get_boundary_geometry(selected_name)
                if gdf is not None:
                    selected_geometry = geemap.geopandas_to_ee(gdf)
                    st.success(f"✅ Location found: {selected_name}")
                else:
                    st.warning("⚠️ Could not find boundary for this location")
        else:
            st.info("No suggestions found. Try a different search term.")

    return selected_geometry, selected_name


def display_band_info(satellite, rgb_bands=None):
    """Display information about selected bands"""
    if rgb_bands:
        st.subheader("📡 Selected Bands Information")
        available_bands = SATELLITE_BANDS[satellite]

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric(
                "Red Band", f"{rgb_bands[0]} ({available_bands[rgb_bands[0]]})")
        with col2:
            st.metric("Green Band",
                      f"{rgb_bands[1]} ({available_bands[rgb_bands[1]]})")
        with col3:
            st.metric(
                "Blue Band", f"{rgb_bands[2]} ({available_bands[rgb_bands[2]]})")


def main():
    """Main application function"""
    # Initialize Earth Engine
    initialize_earth_engine()

    # Title and description
    st.title("🛰️ Enhanced Satellite Index Viewer")
    st.markdown("""
    Visualize satellite imagery with customizable RGB band combinations and unified 
    percentile-based stretching for both RGB composites and spectral indices.
    """)

    # Create sidebar controls
    (start_date, end_date, satellite, index, cloud_percent,
     rgb_bands, stretch_method, min_percentile, max_percentile, min_val, max_val) = create_sidebar_controls()

    # Create location input
    selected_geometry, selected_name = create_location_input()

    # Display band information for RGB
    if index == 'RGB' and rgb_bands:
        display_band_info(satellite, rgb_bands)

    # Display stretch information
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Visualization Method", stretch_method)
    with col2:
        if stretch_method == "Percentile":
            st.metric("Percentile Range",
                      f"{min_percentile}% - {max_percentile}%")
        else:
            st.metric("Manual Range", f"{min_val} - {max_val}")
    with col3:
        st.metric("Selected Index/Bands", index)

    # Display map button and processing
    if st.button("🗺️ Generate Map", type="primary"):
        if not selected_name or not selected_geometry:
            st.error("⚠️ Please enter a location and select a valid suggestion.")
        else:
            with st.spinner(f"Processing {index} for {selected_name}..."):
                show_map(
                    selected_geometry, start_date, end_date, satellite,
                    index, cloud_percent, stretch_method, min_percentile, max_percentile,
                    min_val, max_val, rgb_bands
                )

    # Information sections
    col1, col2 = st.columns(2)

    with col1:
        with st.expander("ℹ️ About Indices"):
            st.markdown("""
            - **NDVI**: Normalized Difference Vegetation Index - vegetation health
            - **NDBI**: Normalized Difference Built-up Index - urban areas
            - **NDWI**: Normalized Difference Water Index - water bodies
            - **NDSI**: Normalized Difference Snow Index - snow cover
            - **NBR**: Normalized Burn Ratio - fire damage assessment
            - **EVI**: Enhanced Vegetation Index - improved vegetation analysis
            - **SAVI**: Soil Adjusted Vegetation Index - vegetation with soil correction
            - **MSAVI**: Modified Soil Adjusted Vegetation Index
            - **GNDVI**: Green Normalized Difference Vegetation Index
            - **GCI**: Green Chlorophyll Index - chlorophyll content
            """)

    with col2:
        with st.expander("🎨 RGB Band Combinations"):
            st.markdown("""
            **Popular Band Combinations:**
            - **Natural Color**: Red-Green-Blue (true color)
            - **False Color**: NIR-Red-Green (vegetation in red)
            - **SWIR Composite**: SWIR2-NIR-Red (geology/moisture)
            - **Urban**: SWIR2-SWIR1-Red (urban features)
            
            **Percentile Stretching:**
            - **2%-98%**: Default Earth Engine style (removes extreme outliers)
            - **1%-99%**: More aggressive stretch (higher contrast)
            - **5%-95%**: Conservative stretch (preserves more detail in extremes)
            - **Manual Override**: Set exact min/max values when needed
            """)

    # Advanced features section
    with st.expander("🔧 Advanced Features & Tips"):
        st.markdown("""
        **Unified Percentile Controls:**
        - Same percentile settings work for both RGB composites and spectral indices
        - Percentile method automatically calculates optimal display ranges
        - Manual override available when you need precise control
        
        **RGB Band Selection:**
        - Choose any three bands for Red, Green, and Blue channels
        - Create false color composites for different analysis purposes
        - Band labels show spectral regions (Blue, Green, Red, NIR, SWIR1, SWIR2)
        
        **Best Practices:**
        - Use 2%-98% percentiles for most visualizations (Earth Engine default)
        - Try NIR-Red-Green for vegetation analysis
        - Use SWIR bands for geological and moisture analysis
        - Lower percentiles (1%-99%) for higher contrast in urban areas
        """)


if __name__ == "__main__":
    main()
