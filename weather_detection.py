import warnings
warnings.filterwarnings('ignore')
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import boto3
import pytz
from datetime import datetime, timedelta
from botocore import UNSIGNED
from botocore.config import Config
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
import shutil
import zipfile
import json
from geopy.geocoders import Nominatim
from pyproj import Proj
import argparse
from pathlib import Path
import eumdac
from scipy.interpolate import griddata
from matplotlib.patches import Rectangle
from matplotlib.path import Path as MplPath
from shapely.geometry import Polygon as ShapelyPolygon, box as shapely_box
from scipy.spatial import cKDTree

# ==========================================
# CONFIGURATION
# ==========================================
class SatelliteConfig:
    SATELLITES = {
        'GOES-18': {
            'bucket': 'noaa-goes18', 'product': 'ABI-L2-ACMF', 'variable': 'ACM',
            'coverage': {'lon_min': -180, 'lon_max': -100, 'lat_min': -60, 'lat_max': 60},
            'path_format': '{product}/{year}/{day_of_year:03d}/{hour:02d}/',
            'engine': 'h5netcdf', 'projection_var': 'goes_imager_projection',
            'data_source': 'aws', 'standard': 'NOAA', 'sat_longitude': -137.2,
            'dqf_var': 'DQF'
        },
        'GOES-19': {
            'bucket': 'noaa-goes19', 'product': 'ABI-L2-ACMF', 'variable': 'ACM',
            'coverage': {'lon_min': -110, 'lon_max': -25, 'lat_min': -60, 'lat_max': 60},
            'path_format': '{product}/{year}/{day_of_year:03d}/{hour:02d}/',
            'engine': 'h5netcdf', 'projection_var': 'goes_imager_projection',
            'data_source': 'aws', 'standard': 'NOAA', 'sat_longitude': -75.0,
            'dqf_var': 'DQF'
        },
        'MTG-0': {
            'product': 'FCI-2-CLM', 'variable': 'cloud_state',
            'coverage': {'lon_min': -45, 'lon_max': 30, 'lat_min': -60, 'lat_max': 60},
            'engine': 'h5netcdf', 'collection_id': 'EO:EUM:DAT:0678',
            'projection_var': 'mtg_geos_projection', 'sat_longitude': 0.0,
            'data_source': 'eumdac', 'file_ext': '.nc', 'standard': 'EUMETSAT'
        },
        'MSG-2': {
            'product': 'MSGCLMK', 'variable': 'p260537',
            'coverage': {'lon_min': 20, 'lon_max': 95, 'lat_min': -60, 'lat_max': 60},
            'engine': 'cfgrib', 'collection_id': 'EO:EUM:DAT:MSG:CLM-IODC',
            'data_source': 'eumdac', 'file_ext': '.grb', 'standard': 'EUMETSAT',
            'sat_longitude': 45.5
        },
        'Himawari-8': {
            'bucket': 'noaa-himawari8', 'product': 'AHI-L2-FLDK-Clouds',
            'file_pattern': 'CMSK', 'variable': 'CloudMask',
            'coverage': {'lon_min': 85, 'lon_max': 180, 'lat_min': -60, 'lat_max': 60},
            'path_format': '{product}/{year}/{month:02d}/{day:02d}/',
            'engine': 'h5netcdf', 'lat_var': 'Latitude', 'lon_var': 'Longitude',
            'data_source': 'aws', 'standard': 'NOAA', 'sat_longitude': 140.7,
            'dqf_var': 'DQF'
        },
        'Himawari-9': {
            'bucket': 'noaa-himawari9', 'product': 'AHI-L2-FLDK-Clouds',
            'file_pattern': 'CMSK', 'variable': 'CloudMask',
            'coverage': {'lon_min': 85, 'lon_max': 180, 'lat_min': -60, 'lat_max': 60},
            'path_format': '{product}/{year}/{month:02d}/{day:02d}/',
            'engine': 'h5netcdf', 'lat_var': 'Latitude', 'lon_var': 'Longitude',
            'data_source': 'aws', 'standard': 'NOAA', 'sat_longitude': 140.7,
            'dqf_var': 'DQF'
        },
        'GOES-16': {
            'bucket': 'noaa-goes16', 'product': 'ABI-L2-ACMF', 'variable': 'ACM',
            'coverage': {'lon_min': -110, 'lon_max': -25, 'lat_min': -60, 'lat_max': 60},
            'path_format': '{product}/{year}/{day_of_year:03d}/{hour:02d}/',
            'engine': 'h5netcdf', 'projection_var': 'goes_imager_projection',
            'data_source': 'aws', 'standard': 'NOAA', 'sat_longitude': -75.2,
            'dqf_var': 'DQF'
        },
        'MSG-3': {
            # Meteosat-10, SEVIRI prime service at 0°E
            # Same GRIB structure as MSG-2 (3712x3712, paramId=260537)
            # Only difference: longitudeOfSubSatellitePoint = 0 (not 45.5)
            'product': 'MSGCLMK', 'variable': 'p260537',
            'coverage': {'lon_min': -75, 'lon_max': 75, 'lat_min': -75, 'lat_max': 75},
            'engine': 'cfgrib', 'collection_id': 'EO:EUM:DAT:MSG:CLM',
            'sat_longitude': 0.0,
            'data_source': 'eumdac', 'file_ext': '.grb', 'standard': 'EUMETSAT'
        }
    }

# GRIB Code Table 4.217
MSG2_VALUE_MAP = {0: 'Clear over water', 1: 'Clear over land',
                  2: 'Cloud',            3: 'No data'}
MSG2_VALID = {0, 1, 2} 


# VZA reliability threshold (degrees).
# Above this angle the pixel footprint is too distorted and cloud parallax
# displacement becomes too large for reliable cloud fraction estimates.
# EUMETSAT/NOAA operational guidance: 70° for GEO instruments.
# (MODIS uses 60° but that is a LEO sensor at 700 km altitude — GEO at 36000 km
# uses the same pixel-size scaling law but operationally flags at 70°.)
VZA_UNRELIABLE_THRESHOLD = 70.0

# ==========================================
# UTILITIES
# ==========================================
def get_region_details(region):
    clat = region['centroid_lat']
    clon = region['centroid_lon']
    try:
        loc = Nominatim(user_agent="geo_forensic_tool").reverse(
            f"{clat}, {clon}", timeout=10)
        return loc.address if loc else f"Lat:{clat:.4f}, Lon:{clon:.4f}"
    except:
        return f"Lat:{clat:.4f}, Lon:{clon:.4f}"

def get_satellite_style(satellite):
    std = SatelliteConfig.SATELLITES[satellite]['standard']
    if std == 'NOAA':
        colors  = ['#FFFF00', '#ADD8E6', '#CBC3E3', '#FFFFFF', '#000000']
        labels  = ['Clear (0)', 'Prob Clear (1)', 'Prob Cloudy (2)',
                   'Cloudy (3)', 'No Data']
        mapping = {'nodata_val': 4}
    else:
        colors  = ['#0055FF', '#00AA00', '#FFFFFF', '#444444']
        labels  = ['Clear Water (0)', 'Clear Land (1)',
                   'Cloud (2)',       'No Data (3)']
        mapping = {'nodata_val': 3}
    return ListedColormap(colors), labels, mapping


# ==========================================
# SATELLITE PRIORITY / VIEWING GEOMETRY
# ==========================================

def compute_viewing_zenith_angle(lat_c, lon_c, sat_lon_0):
    """
    Returns the satellite viewing zenith angle (VZA) in degrees at the
    given ground point.

    VZA is measured at the target: angle between the local vertical and
    the line of sight to the satellite. Lower VZA = satellite closer to
    overhead = better geometry, smaller pixel, less parallax.

    Uses the Earth-center / satellite / target triangle:
      beta = arccos(cos(lat) * cos(lon - sat_lon_0))
      VZA  = arctan(r_s * sin(beta) / (r_s * cos(beta) - R_E))
    where r_s = 42164 km (GEO orbital radius), R_E = 6378 km.
    """
    R_E = 6378.0
    r_s = 42164.0
    lat_r  = np.deg2rad(lat_c)
    dlon_r = np.deg2rad(lon_c - sat_lon_0)
    beta   = np.arccos(np.clip(np.cos(lat_r) * np.cos(dlon_r), -1.0, 1.0))
    vza    = np.arctan2(r_s * np.sin(beta), r_s * np.cos(beta) - R_E)
    return float(np.rad2deg(vza))


def assign_priority_order(results):
    """
    Assigns priority_order (1, 2, or 3) to each result in-place and sorts
    the list so index 0 is Priority 1.

    Ranking rules:
      1. Lower nodata_pct wins first: a satellite with 0% missing data
         always ranks above one with any missing data, regardless of VZA.
      2. Among satellites with equal nodata, lower VZA wins: the satellite
         closer to overhead has less pixel distortion and parallax.
      3. Priority is capped at 3; any extra satellites also get 3.
    """
    results.sort(key=lambda r: (r.get('nodata_pct', 0.0), r.get('vza_deg', 90.0)))
    for rank, r in enumerate(results, start=1):
        r['priority_order'] = min(rank, 3)

# --- Exact constants from GRIB keys ---
_NI      = 3712
_NJ      = 3712
_DX      = 3622          # GRIB key dx (angular step units, see below)
_DY      = 3610          # GRIB key dy
_XP      = 1856000       # GRIB key Xp  (sub-satellite column * 1000)
_YP      = 1856000       # GRIB key Yp  (sub-satellite row    * 1000)
_NR      = 6610700       # GRIB key Nr  (H/R_eq * 10^6)
_R_EQ_KM = 6378.14       # GRIB key earthMajorAxis (km)
_R_PO_KM = 6356.755      # GRIB key earthMinorAxis (km)
_SAT_LON = 45.5          # GRIB key longitudeOfSubSatellitePoint (deg)

# Derived
_H_KM    = (_NR / 1e6) * _R_EQ_KM   # satellite distance from Earth centre (km)
_MSG2_CFAC = 781648399
_DX_RAD  = (2**16) / _MSG2_CFAC   # ≈ 8.384e-5 rad/pixel (column step)
_DY_RAD  = (2**16) / _MSG2_CFAC   # ≈ 8.384e-5 rad/pixel (row step, same CFAC)


def _msg2_latlon_grid(sat_lon=_SAT_LON):
    """
    Build (lats2d, lons2d) shaped (_NJ, _NI) from exact GRIB key values.
    Works for any SEVIRI prime/IODC product; only sat_lon differs.

    Scan order from GRIB keys (scanningMode=192):
      iScansNegatively=1 - column index increases West (pixel 0 = easternmost)
      jScansPositively=1 - row    index increases North (pixel 0 = southernmost)

    WMO GDT 3.90 angular offset formula:
      x =  (Xp/1000 - c) * dx_rad   [positive = West of sub-satellite point]
      y =  (l - Yp/1000) * dy_rad   [positive = North of sub-satellite point]
    """
    Xp_px = _XP / 1000.0   # 1856.0
    Yp_px = _YP / 1000.0   # 1856.0

    c = np.arange(_NI, dtype=np.float64)   # column index 0…3711
    l = np.arange(_NJ, dtype=np.float64)   # row    index 0…3711

    # iScansNegatively=1: col 0 is East edge; col increases going West
    # EUMETSAT formula: Se = -Sd*sin(X_ang), so X_ang must be NEGATIVE for East.
    # X_ang = (c - Xp) * DX_RAD: at c=0 (East) -> negative ; at c=3711 (West) -> positive 
    # jScansPositively=1: row 0 is South edge; Y_ang = (l - Yp)*DY -> negative at row 0 (South) 
    C2d, L2d = np.meshgrid(np.arange(_NI, dtype=np.float64),
                           np.arange(_NJ, dtype=np.float64))
    X_ang = (C2d - Xp_px) * _DX_RAD
    Y_ang = (L2d - Yp_px) * _DY_RAD

    # EUMETSAT geostationary inverse projection (EUM/MSG/SPE/057 eq. 5.2)
    R_eq = _R_EQ_KM
    R_po = _R_PO_KM
    H    = _H_KM

    a = (np.sin(X_ang)**2 +
         np.cos(X_ang)**2 * (np.cos(Y_ang)**2 +
                              (R_eq/R_po)**2 * np.sin(Y_ang)**2))
    b = -2.0 * H * np.cos(X_ang) * np.cos(Y_ang)
    c_coeff = H**2 - R_eq**2

    disc = b**2 - 4.0*a*c_coeff
    off_earth = disc < 0
    disc = np.where(off_earth, 0.0, disc)

    Sd = (-b - np.sqrt(disc)) / (2.0 * a)
    Sn =  H    - Sd * np.cos(X_ang) * np.cos(Y_ang)
    Se = -Sd * np.sin(X_ang)
    Sz =  Sd * np.cos(X_ang) * np.sin(Y_ang)

    lat_rad = np.arctan((R_eq/R_po)**2 * Sz / np.sqrt(Sn**2 + Se**2))
    lon_rad = np.arctan2(Se, Sn) + np.deg2rad(sat_lon)

    lats2d = np.rad2deg(lat_rad)
    lons2d = np.rad2deg(lon_rad)
    lons2d = (lons2d + 180) % 360 - 180

    lats2d[off_earth] = np.nan
    lons2d[off_earth] = np.nan

    valid_pct = (~off_earth).sum() / (_NI * _NJ) * 100
    print(f"  [SEVIRI GRID lon_0={sat_lon}°] {_NJ}x{_NI} pixels decoded, "
          f"{valid_pct:.1f}% on-Earth")
    print(f"  [SEVIRI GRID] Lat: {np.nanmin(lats2d):.2f}–{np.nanmax(lats2d):.2f}  "
          f"Lon: {np.nanmin(lons2d):.2f}–{np.nanmax(lons2d):.2f}")
    return lats2d, lons2d   # shape (_NJ, _NI)


def get_msg2_scatter(ds, sat_lon=_SAT_LON, sat_label='SEVIRI'):
    var_name = 'p260537'
    print(f"  [{sat_label}] Reconstructing coordinates from GRIB projection keys (lon_0={sat_lon}°)...")
    lats2d, lons2d = _msg2_latlon_grid(sat_lon)

    # cfgrib loads as flat 1-D; reshape to 2-D (NJ rows × NI cols)
    raw_vals = ds[var_name].values.reshape(_NJ, _NI).astype(np.float64)

    # Flatten all arrays
    lats_1d = lats2d.flatten()
    lons_1d = lons2d.flatten()
    vals_1d = raw_vals.flatten()

    # Valid = on-Earth AND actual cloud-mask value (0–3), not GRIB fill
    valid = (
        ~np.isnan(lats_1d) &
        (vals_1d >= 0) & (vals_1d <= 3)
    )

    v_lats = lats_1d[valid]
    v_lons = lons_1d[valid]
    v_vals = vals_1d[valid]

    print(f"  [{sat_label}] Valid scatter points: {len(v_vals):,}")
    u, c = np.unique(v_vals, return_counts=True)
    print(f"  [{sat_label}] Full-disk distribution:")
    for val, cnt in zip(u, c):
        print(f"    {int(val)} ({MSG2_VALUE_MAP[int(val)]}): "
              f"{cnt:,}  ({cnt/len(v_vals)*100:.1f}%)")
    return v_lats, v_lons, v_vals


def extract_msg2_bbox(v_lats, v_lons, v_vals, region, sat_label='SEVIRI'):
    b_weights_full = _msg_scatter_pixel_weights(v_lats, v_lons, region)
    in_region = b_weights_full > 0
    b_lats    = v_lats[in_region]
    b_lons    = v_lons[in_region]
    b_vals    = v_vals[in_region]
    b_weights = b_weights_full[in_region]
    eff       = b_weights.sum()
    n_full    = int((b_weights == 1.0).sum())
    n_part    = int((b_weights < 1.0).sum())
    print(f"  [{sat_label} BBOX] Pixels in region: {len(b_vals):,}  "
          f"(full: {n_full:,}, partial: {n_part:,}, effective: {eff:.1f})")
    if len(b_vals) > 0:
        u, c = np.unique(b_vals, return_counts=True)
        vd = dict(zip(u.astype(int), c))
        for val, cnt in zip(u, c):
            print(f"    {int(val)} ({MSG2_VALUE_MAP[int(val)]}): "
                  f"{cnt:,}  ({cnt/len(b_vals)*100:.1f}%)")
        if vd.get(0, 0) > 0 and vd.get(1, 0) > 0:
            water_pct = vd.get(0, 0) / len(b_vals) * 100
            if water_pct > 5:
                print(f"  [NOTE] {water_pct:.0f}% Clear Water (0) alongside Clear Land (1) —")
                print(f"         SEVIRI uses a static land-sea mask; some land pixels near")
                print(f"         coastlines or in desert regions may be misclassified.")
                print(f"         Both Value 0 and Value 1 mean 'clear' — combined clear = "
                      f"{(vd.get(0,0)+vd.get(1,0))/len(b_vals)*100:.1f}%")
    return b_lats, b_lons, b_vals, b_weights


def build_msg2_display_grid(v_lats, v_lons, v_vals, dlon=0.15, dlat=0.15):
    valid = np.isin(v_vals, list(MSG2_VALID))
    lon_grid = np.arange(-20, 110, dlon)
    lat_grid = np.arange(-60,  60, dlat)
    LON, LAT = np.meshgrid(lon_grid, lat_grid)
    print(f"  [MSG-2 DISPLAY] Gridding {valid.sum():,} points "
          f"-> {len(lat_grid)}x{len(lon_grid)}...")
    full_grid = griddata(
        (v_lons[valid], v_lats[valid]),
        v_vals[valid],
        (LON, LAT), method='nearest'
    )
    return lon_grid, lat_grid, full_grid


# ==========================================
# REGION HELPERS (polygon input + fractional pixel weights)
# ==========================================

def _build_region(corners):
    """
    Build a unified region dict from 4 (lat, lon) corner tuples.
    Input order does not matter — corners are automatically sorted via convex hull
    so that any permutation of the 4 points produces the correct polygon.
    Works correctly for convex quadrilaterals (any tilted rectangle / parallelogram).
    Returns dict with:
      polygon      => shapely Polygon (lon, lat coords = x, y)
      corners      => corners in convex-hull order [(lat,lon), ...]
      bbox         => axis-aligned envelope {north, south, east, west}
      centroid_lat => polygon centroid latitude
      centroid_lon => polygon centroid longitude
    """
    from shapely.geometry import MultiPoint
    # Convex hull of the 4 points always gives the correct outer polygon
    # regardless of the order the user provides them in.
    hull = MultiPoint([(lon, lat) for lat, lon in corners]).convex_hull
    if hull.geom_type != 'Polygon':
        raise ValueError("The 4 corners do not form a valid quadrilateral "
                         "(are some points collinear or duplicate?)")
    poly = hull
    centroid = poly.centroid
    lats = [c[0] for c in corners]
    lons = [c[1] for c in corners]
    # Rebuild corners list in the convex-hull winding order for consistent plotting
    hull_coords = list(poly.exterior.coords)[:-1]   # drop repeated closing point
    ordered_corners = [(lat, lon) for lon, lat in hull_coords]
    return {
        'polygon':      poly,
        'corners':      ordered_corners,
        'bbox':         {'north': max(lats), 'south': min(lats),
                         'east':  max(lons),  'west':  min(lons)},
        'centroid_lat': centroid.y,
        'centroid_lon': centroid.x,
    }


def _compute_pixel_half_extents(lats, lons):
    """
    Estimate per-pixel half-extents in lat and lon from local grid spacing.
    Uses np.gradient (central differences interior, one-sided at edges).
    Clamps to physically plausible range; NaN/inf replaced with array median.
    Returns dlat_half, dlon_half arrays (same shape as lats/lons).
    """
    dlat_row, dlat_col = np.gradient(lats)
    dlon_row, dlon_col = np.gradient(lons)
    # Conservative axis-aligned bounding box of the (possibly tilted) pixel footprint
    dlat_half = 0.5 * (np.abs(dlat_row) + np.abs(dlat_col))
    dlon_half = 0.5 * (np.abs(dlon_row) + np.abs(dlon_col))
    for arr in (dlat_half, dlon_half):
        finite = arr[np.isfinite(arr)]
        med = float(np.nanmedian(finite)) if len(finite) > 0 else 0.05
        arr[~np.isfinite(arr)] = med
    dlat_half = np.clip(dlat_half, 0.001, 2.0)
    dlon_half = np.clip(dlon_half, 0.001, 2.0)
    return dlat_half, dlon_half


def _compute_pixel_weights(lats, lons, region):
    """
    Compute fractional overlap weight [0,1] for every pixel with the region polygon.
    Strategy:
      1. Prescreen with the axis-aligned envelope (fast numpy).
      2. Vectorised point-in-polygon (matplotlib Path) to classify interior vs edge.
      3. Eroded polygon: pixels whose centres fall inside -> guaranteed fully covered
         -> weight = 1.0 (avoids per-pixel shapely call for the bulk of pixels).
      4. Remaining candidates (edge pixels) -> exact shapely box-intersection.
    """
    orig_shape = lats.shape
    lats_f = lats.ravel()
    lons_f = lons.ravel()

    weights = np.zeros(len(lats_f), dtype=float)
    poly = region['polygon']
    bbox = region['bbox']

    dlat_h, dlon_h = _compute_pixel_half_extents(lats, lons)
    dlat_h_f = dlat_h.ravel()
    dlon_h_f = dlon_h.ravel()

    med_half = max(float(np.nanmedian(dlat_h_f)), float(np.nanmedian(dlon_h_f)))

    # Step 1: envelope prescreen
    in_env = (
        np.isfinite(lats_f) & np.isfinite(lons_f) &
        (lats_f >= bbox['south'] - med_half) & (lats_f <= bbox['north'] + med_half) &
        (lons_f >= bbox['west']  - med_half) & (lons_f <= bbox['east']  + med_half)
    )
    cand_idx = np.where(in_env)[0]
    if len(cand_idx) == 0:
        return weights.reshape(orig_shape)

    c_lons = lons_f[cand_idx]
    c_lats = lats_f[cand_idx]
    centers = np.column_stack([c_lons, c_lats])  # (x=lon, y=lat)

    # Step 2: vectorised PIP
    poly_coords = np.array(list(poly.exterior.coords))
    mpl_path = MplPath(poly_coords[:, :2])
    inside_poly = mpl_path.contains_points(centers)

    # Step 3: eroded polygon for bulk-interior classification
    eroded = poly.buffer(-med_half)
    if not eroded.is_empty and eroded.geom_type == 'Polygon':
        ec = np.array(list(eroded.exterior.coords))
        fully_interior = MplPath(ec[:, :2]).contains_points(centers)
    elif not eroded.is_empty and eroded.geom_type == 'MultiPolygon':
        fully_interior = np.zeros(len(cand_idx), dtype=bool)
        for geom in eroded.geoms:
            ec = np.array(list(geom.exterior.coords))
            fully_interior |= MplPath(ec[:, :2]).contains_points(centers)
    else:
        fully_interior = np.zeros(len(cand_idx), dtype=bool)

    # Interior: weight = 1.0
    interior_mask = inside_poly & fully_interior
    weights[cand_idx[interior_mask]] = 1.0

    # Step 4: edge pixels – exact shapely intersection
    for k in np.where(~interior_mask)[0]:
        gi = cand_idx[k]
        if not (np.isfinite(lats_f[gi]) and np.isfinite(lons_f[gi])):
            continue
        dh_lat = dlat_h_f[gi]
        dh_lon = dlon_h_f[gi]
        pix = shapely_box(lons_f[gi] - dh_lon, lats_f[gi] - dh_lat,
                          lons_f[gi] + dh_lon, lats_f[gi] + dh_lat)
        pix_area = pix.area
        if pix_area > 0:
            try:
                weights[gi] = poly.intersection(pix).area / pix_area
            except Exception:
                weights[gi] = 1.0 if inside_poly[k] else 0.0

    return weights.reshape(orig_shape)


def _msg_scatter_pixel_weights(v_lats, v_lons, region):
    """
    Fractional overlap weights for MSG scatter data (1-D decoded pixel centres).
    Pixel half-size estimated via nearest-neighbour spacing (captures limb growth).
    """
    poly = region['polygon']
    bbox = region['bbox']
    n = len(v_lats)
    weights = np.zeros(n, dtype=float)

    MARGIN = 0.2  # degrees – larger than any SEVIRI pixel footprint
    in_env = (
        np.isfinite(v_lats) & np.isfinite(v_lons) &
        (v_lats >= bbox['south'] - MARGIN) & (v_lats <= bbox['north'] + MARGIN) &
        (v_lons >= bbox['west']  - MARGIN) & (v_lons <= bbox['east']  + MARGIN)
    )
    cand_idx = np.where(in_env)[0]
    if len(cand_idx) == 0:
        return weights

    c_lats = v_lats[cand_idx]
    c_lons = v_lons[cand_idx]

    # Local pixel half-size from nearest-neighbour distance in (lat,lon) space
    pts = np.column_stack([c_lats, c_lons])
    nn_dists, _ = cKDTree(pts).query(pts, k=2)
    nn_half = nn_dists[:, 1] * 0.5  # half of spacing to nearest neighbour

    # Vectorised PIP
    poly_coords = np.array(list(poly.exterior.coords))
    mpl_path = MplPath(poly_coords[:, :2])
    centers = np.column_stack([c_lons, c_lats])
    inside_poly = mpl_path.contains_points(centers)

    # Eroded polygon for interior classification
    med_half = float(np.median(nn_half))
    eroded = poly.buffer(-med_half)
    if not eroded.is_empty and eroded.geom_type == 'Polygon':
        ec = np.array(list(eroded.exterior.coords))
        fully_interior = MplPath(ec[:, :2]).contains_points(centers)
    elif not eroded.is_empty and eroded.geom_type == 'MultiPolygon':
        fully_interior = np.zeros(len(cand_idx), dtype=bool)
        for geom in eroded.geoms:
            ec = np.array(list(geom.exterior.coords))
            fully_interior |= MplPath(ec[:, :2]).contains_points(centers)
    else:
        fully_interior = np.zeros(len(cand_idx), dtype=bool)

    # Interior: weight = 1.0
    interior_mask = inside_poly & fully_interior
    weights[cand_idx[interior_mask]] = 1.0

    # Edge: exact shapely intersection
    for k in np.where(~interior_mask)[0]:
        gi = cand_idx[k]
        dh = nn_half[k]
        pix = shapely_box(c_lons[k] - dh, c_lats[k] - dh,
                          c_lons[k] + dh, c_lats[k] + dh)
        pix_area = pix.area
        if pix_area > 0:
            try:
                weights[gi] = poly.intersection(pix).area / pix_area
            except Exception:
                weights[gi] = 1.0 if inside_poly[k] else 0.0

    return weights


# ==========================================
# GENERAL EXTRACTION (non-MSG2)
# ==========================================
def extract_bounding_box_data(satellite, ds, region):
    config = SatelliteConfig.SATELLITES[satellite]
    bbox   = region['bbox']
    print(f"\n--- EXTRACTING DATA FOR {satellite} ---")
    lat_var = config.get('lat_var', 'latitude')
    lon_var = config.get('lon_var', 'longitude')
    if lat_var not in ds and 'Latitude' in ds: lat_var = 'Latitude'
    if lon_var not in ds and 'Longitude' in ds: lon_var = 'Longitude'
    try:
        if lat_var in ds and lon_var in ds:
            lats = ds[lat_var].values
            lons = ds[lon_var].values
            offset = config.get('lon_offset', 0.0)
            if abs(offset) > 0.1: lons = lons + offset
            if lats.ndim == 1 and lons.ndim == 1:
                if not (len(lats) == len(lons) and len(lats) > 10000):
                    lons, lats = np.meshgrid(lons, lats)
            lons = (lons + 180) % 360 - 180
        elif 'x' in ds and 'y' in ds:
            sat_lon = config.get('sat_longitude', 0.0)
            sat_h   = 35786023.0   # metres above Earth surface
            r_eq    = 6378137.0
            r_pol   = 6356752.0
            sweep   = 'x'          # GOES: sweep_angle_axis=x; FCI: sweep_angle_axis=y

            pv_name = config.get('projection_var', 'goes_imager_projection')
            if pv_name in ds:
                pvar = ds[pv_name]
                if hasattr(pvar, 'perspective_point_height'):
                    sat_h = float(pvar.perspective_point_height)
                if hasattr(pvar, 'longitude_of_projection_origin'):
                    sat_lon = float(pvar.longitude_of_projection_origin)
                if hasattr(pvar, 'semi_major_axis'):
                    r_eq = float(pvar.semi_major_axis)
                if hasattr(pvar, 'semi_minor_axis'):
                    r_pol = float(pvar.semi_minor_axis)
                if hasattr(pvar, 'sweep_angle_axis'):
                    sweep = str(pvar.sweep_angle_axis)

            # Scale detection: radians (|max| < 1.0) - multiply by sat_h to get metres
            x_raw = ds['x'].values
            y_raw = ds['y'].values
            xx = x_raw * sat_h if np.nanmax(np.abs(x_raw)) < 1.0 else x_raw
            yy = y_raw * sat_h if np.nanmax(np.abs(y_raw)) < 1.0 else y_raw

            # sweep_angle_axis=y means the projection sweeps in y (FCI convention)
            proj_str = (f"+proj=geos +lon_0={sat_lon} +h={sat_h} "
                        f"+a={r_eq} +b={r_pol} +sweep={sweep} +units=m +no_defs")
            p = Proj(proj_str)
            X, Y = np.meshgrid(xx, yy)
            lons, lats = p(X, Y, inverse=True)
            off = (np.abs(lons) > 360) | (np.abs(lats) > 90)
            lons[off] = np.nan
            lats[off] = np.nan
            lons = (lons + 180) % 360 - 180
            lons[off] = np.nan

        else:
            print("  [ERROR] Cannot determine coordinates."); return None, None, None

        weights = _compute_pixel_weights(lats, lons, region)
        n_sel  = int(np.sum(weights > 0))
        n_full = int(np.sum(weights == 1.0))
        n_part = n_sel - n_full
        eff    = float(weights.sum())
        print(f"  Pixels in region: {n_sel:,}  "
              f"(full: {n_full:,}, partial: {n_part:,}, effective: {eff:.1f})")
        return weights, lats, lons
    except Exception as e:
        print(f"  [CRITICAL ERROR] {e}"); return None, None, None


# ==========================================
# ANALYSIS
# ==========================================
def analyze_centroid(region, lats_in, lons_in, data_in, satellite):
    clat = region['centroid_lat']
    clon = region['centroid_lon']
    print(f"\n=== CENTROID ANALYSIS FOR {satellite} ===")
    flat_lats = np.asarray(lats_in).flatten()
    flat_lons = np.asarray(lons_in).flatten()
    flat_vals = np.asarray(data_in).flatten()
    valid = ~np.isnan(flat_vals)
    if not np.any(valid): print("  No valid pixels."); return
    flat_lats, flat_lons, flat_vals = (flat_lats[valid], flat_lons[valid],
                                       flat_vals[valid])
    dists = np.sqrt((flat_lats - clat)**2 + (flat_lons - clon)**2)
    i     = np.argmin(dists)
    is_eu = SatelliteConfig.SATELLITES[satellite]['standard'] == 'EUMETSAT'
    lbl   = (MSG2_VALUE_MAP.get(int(flat_vals[i]), '?') if is_eu
             else {0:'Clear',1:'Prob Clear',2:'Prob Cloudy',3:'Cloudy'
                   }.get(int(flat_vals[i]), '?'))
    print(f"  Target Centroid: {clat:.4f} N, {clon:.4f} E")
    print(f"  Closest Pixel:   {flat_lats[i]:.4f} N, {flat_lons[i]:.4f} E")
    print(f"  Distance:        {dists[i]:.4f} deg")
    print(f"  Pixel Value:     {int(flat_vals[i])} ({lbl})")


def calculate_comprehensive_statistics(data_in, satellite, weights_in=None):
    """
    weights_in: per-pixel fractional overlap weight [0,1].  None -> all pixels weight 1.
    Counts are replaced by weight sums so edge pixels contribute proportionally.
    total_pixels is the sum of weights (effective pixel count, may be float).
    """
    config = SatelliteConfig.SATELLITES[satellite]
    std    = config['standard']
    flat   = np.asarray(data_in).flatten()
    w      = (np.asarray(weights_in).flatten()
              if weights_in is not None else np.ones(len(flat), dtype=float))

    # Capture total weight before any filtering (for nodata_pct denominator)
    total_w_raw = float(np.nansum(w))

    # Drop NaN values (xarray-masked fill values land here)
    valid  = ~np.isnan(flat)
    nan_w  = float(w[~valid].sum())   # weight of NaN-masked fill pixels
    flat   = flat[valid]
    w      = w[valid]

    print(f"\n=== COMPREHENSIVE STATISTICS FOR {satellite} ({std}) ===")
    print(f"  Valid pixels (all values): {len(flat):,}  (effective: {w.sum():.1f})")
    if len(flat) == 0:
        print("  [WARNING] No valid data.")
        return {"cloud_pct": 0, "clear_pct": 0, "total_pixels": 0.0, "nodata_pct": 100}

    # Weighted value distribution
    u_all = np.unique(flat)
    print("  Value distribution (weighted):")
    for val in u_all:
        mask_v = flat == val
        wsum   = w[mask_v].sum()
        lbl    = (MSG2_VALUE_MAP.get(int(val), '?') if std == 'EUMETSAT'
                  else {0:'Clear',1:'Prob Clear',2:'Prob Cloudy',3:'Cloudy'
                        }.get(int(val), '?'))
        print(f"    {int(val)} ({lbl}): {mask_v.sum():,} px  "
              f"(eff {wsum:.1f}, {wsum/w.sum()*100:.1f}%)")

    def wsum(val):
        return float(w[flat == val].sum())

    if std == 'EUMETSAT':
        v0_w = wsum(0); v1_w = wsum(1); v2_w = wsum(2); v3_w = wsum(3)
        total_w = v0_w + v1_w + v2_w + v3_w
        pct = lambda n: (n / total_w * 100) if total_w > 0 else 0

        nodata_pct      = pct(v3_w)
        clear_water_pct = pct(v0_w)
        clear_land_pct  = pct(v1_w)
        cloud_pct       = pct(v2_w)
        clear_pct       = pct(v0_w + v1_w)

        print(f"\n  ABSOLUTE COVERAGE (% of effective pixels in region):")
        if nodata_pct > 50:
            print(f"  [WARNING] >50% No-data — result may be unreliable "
                  f"(possible night-time / edge-of-disk / processing gap)")
        print(f"  Clear over water (0): {v0_w:.1f}  ({clear_water_pct:.2f}%)")
        print(f"  Clear over land  (1): {v1_w:.1f}  ({clear_land_pct:.2f}%)")
        print(f"  Cloud            (2): {v2_w:.1f}  ({cloud_pct:.2f}%)")
        print(f"  No-data          (3): {v3_w:.1f}  ({nodata_pct:.2f}%)")
        print(f"  Effective total:       {total_w:.1f}")
        print(f"  >> CLEAR: {clear_pct:.2f}%   CLOUD: {cloud_pct:.2f}%   NO-DATA: {nodata_pct:.2f}%")
        total_all = total_w
    else:
        valid_mask = flat <= 3
        fill_w = float(w[~valid_mask].sum())   # weight of fill pixels (>3, e.g. 128/255)
        fv = flat[valid_mask]; wv = w[valid_mask]
        def wv_sum(val): return float(wv[fv == val].sum())
        v0_w,v1_w,v2_w,v3_w = wv_sum(0),wv_sum(1),wv_sum(2),wv_sum(3)
        total_w = v0_w + v1_w + v2_w + v3_w
        pct = lambda n: (n / total_w * 100) if total_w > 0 else 0
        # nodata_pct: fraction of ALL pixels (including NaN and fill) with no valid value
        nodata_pct = ((nan_w + fill_w) / total_w_raw * 100) if total_w_raw > 0 else 0
        print(f"  Fill/no-data (>3):   {fill_w:.1f} + NaN: {nan_w:.1f}  ({nodata_pct:.2f}%)")
        print(f"  Clear (0):           {v0_w:.1f}  ({pct(v0_w):.2f}%)")
        print(f"  Probably Clear (1):  {v1_w:.1f}  ({pct(v1_w):.2f}%)")
        print(f"  Probably Cloudy (2): {v2_w:.1f}  ({pct(v2_w):.2f}%)")
        print(f"  Cloudy (3):          {v3_w:.1f}  ({pct(v3_w):.2f}%)")
        print(f"  Effective valid:     {total_w:.1f}")
        if nodata_pct > 50:
            print(f"  [WARNING] >50% fill/no-data — result unreliable "
                  f"(edge-of-disk / DQF failure / processing gap)")
        # Probability-weighted clear/cloud (same coefficients as before)
        wc_clear = 1.0*v0_w + 0.67*v1_w + 0.34*v2_w + 0.0*v3_w
        wc_cloud = 0.0*v0_w + 0.33*v1_w + 0.66*v2_w + 1.0*v3_w
        clear_pct = (wc_clear / total_w * 100) if total_w > 0 else 0
        cloud_pct = (wc_cloud / total_w * 100) if total_w > 0 else 0
        print(f"  >> WEIGHTED CLEAR: {clear_pct:.2f}%   CLOUD: {cloud_pct:.2f}%   NO-DATA: {nodata_pct:.2f}%")
        total_all  = total_w

    return {"cloud_pct": cloud_pct, "clear_pct": clear_pct,
            "total_pixels": total_all, "nodata_pct": nodata_pct}


# ==========================================
# VISUALIZATION
# ==========================================
def plot_msg2_dual_view(sat, region, address, utc_dt, s_dir,
                        vals2d_full,
                        b_lats, b_lons, b_vals,
                        sat_lon=_SAT_LON):
    """
    Two-panel SEVIRI GRIB plot (used for MSG-2 and MSG-3).
    vals2d_full: raw (3712,3712) value array — used for native GEOS full-disk.
    sat_lon: sub-satellite longitude (45.5 for MSG-2 IODC, 0.0 for MSG-3 prime).
    """
    cmap, labels, _ = get_satellite_style(sat)
    plt.style.use('dark_background')

    bbox = region['bbox']
    valid_disp = np.isin(b_vals, list(MSG2_VALID))
    if valid_disp.sum() > 0:
        gx = np.linspace(bbox['west'],  bbox['east'],  800)
        gy = np.linspace(bbox['south'], bbox['north'], 800)
        GX, GY = np.meshgrid(gx, gy)
        zoom_grid = griddata(
            (b_lons[valid_disp], b_lats[valid_disp]),
            b_vals[valid_disp],
            (GX, GY), method='nearest')
    else:
        zoom_grid = None

    # Native geostationary projection centred on sub-satellite longitude
    geo_proj = ccrs.Geostationary(central_longitude=sat_lon,
                                   satellite_height=35786000)
    fig = plt.figure(figsize=(20, 8))
    ax1 = fig.add_subplot(1, 2, 1, projection=geo_proj)
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())

    # Full-disk: render raw 2D array subsampled to ~1400×1400
    step_fd = max(1, _NI // 1400)
    raw_ss  = vals2d_full[::step_fd, ::step_fd].astype(float)
    half_ext = 5568000.0
    ax1.imshow(np.fliplr(raw_ss),
               origin='lower',
               extent=[-half_ext, half_ext, -half_ext, half_ext],
               transform=geo_proj,
               cmap=cmap, vmin=0, vmax=3, interpolation='nearest')
    # region outline (polygon corners)
    _corners_c = region['corners'] + [region['corners'][0]]
    lons_b = [c[1] for c in _corners_c]
    lats_b = [c[0] for c in _corners_c]
    ax1.plot(lons_b, lats_b, color='red', linewidth=3,
             transform=ccrs.PlateCarree(), zorder=10)
    ax1.coastlines(color='white', linewidth=0.8)
    ax1.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')
    ax1.set_title(f"{sat} Full Disk View\n{utc_dt.strftime('%Y-%m-%d %H:%M')} UTC")
    ax1.set_global()

    if zoom_grid is not None:
        img2 = ax2.imshow(
            zoom_grid,
            extent=[bbox['west'], bbox['east'], bbox['south'], bbox['north']],
            origin='lower', transform=ccrs.PlateCarree(),
            cmap=cmap, vmin=0, vmax=3, interpolation='nearest')
    else:
        img2 = ax2.imshow(
            np.full((2,2), 3.0),
            extent=[bbox['west'], bbox['east'], bbox['south'], bbox['north']],
            origin='lower', transform=ccrs.PlateCarree(),
            cmap=cmap, vmin=0, vmax=3)
        ax2.text(0.5, 0.5, 'No Data', transform=ax2.transAxes,
                 ha='center', va='center', color='white', fontsize=14)

    ax2.set_extent([bbox['west'], bbox['east'], bbox['south'], bbox['north']],
                   crs=ccrs.PlateCarree())
    # draw polygon outline on zoom panel too
    ax2.plot(lons_b, lats_b, color='red', linewidth=2,
             transform=ccrs.PlateCarree(), zorder=10)
    ax2.coastlines(color='white', linewidth=2)
    ax2.add_feature(cfeature.BORDERS, linestyle=':', linewidth=1)
    ax2.set_title(f"Target Region Analysis\n{address}")
    ax2.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)

    cbar = plt.colorbar(img2, ax=[ax1, ax2], shrink=0.6, pad=0.02, aspect=30)
    cbar.set_ticks([0.375, 1.125, 1.875, 2.625])
    cbar.set_ticklabels(labels)

    save_path = s_dir / f"{sat}_{utc_dt.strftime('%Y%m%d_%H%M')}_dual_analysis.png"
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [Output] Saved: {save_path}")


def _detect_lat_direction(lats_array):
    flat = np.asarray(lats_array).flatten()
    flat = flat[~np.isnan(flat)]
    if len(flat) < 2:
        return 'NS'  # default
    return 'NS' if flat[0] > flat[-1] else 'SN'


def plot_split_view(sat, ds, mask, lats, lons, data_vals,
                    region, address, utc_dt, s_dir):
    """
    Two-panel plot:
      Left  : full-disk context using PlateCarree
      Right : scatter re-gridded bbox zoom with correct
    """
    config = SatelliteConfig.SATELLITES[sat]
    cmap, labels, mapping = get_satellite_style(sat)
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(22, 10))

    sat_lon_fd = config.get('sat_longitude', 0.0)
    if 'GOES-18'   in sat: sat_lon_fd = -137.2
    if 'GOES-16'   in sat: sat_lon_fd =  -75.2
    if 'GOES-19'   in sat: sat_lon_fd =  -75.0
    if 'Himawari'  in sat: sat_lon_fd =  140.7

    geo_proj   = ccrs.Geostationary(central_longitude=sat_lon_fd,
                                     satellite_height=35786000)
    ax_global  = fig.add_subplot(1, 2, 1, projection=geo_proj)

    try:
        raw2d = ds[config['variable']].values
        if raw2d.ndim == 1:
            side = int(np.sqrt(len(raw2d)))
            raw2d = raw2d.reshape(side, side)

        # Subsample for rendering speed (every Nth row/col)
        N_side = raw2d.shape[0]
        step   = max(1, N_side // 1400)   # target ~1400×1400 render
        raw_ss = raw2d[::step, ::step].astype(float)

        if 'x' in ds and 'y' in ds:
            x_vals = ds['x'].values
            y_vals = ds['y'].values
            fd_origin   = 'lower' if y_vals[0] < y_vals[-1] else 'upper'
            need_fliplr = x_vals[0] > x_vals[-1]   # True when x decreases (col0=East)
            # Compute half_ext from actual coordinate range.
            # GOES stores x/y in radians (|max|<1); multiply by satellite height to get metres.
            # Read sat_h from projection variable if available, else use the same default
            # already passed to ccrs.Geostationary above (35786000 m).
            _sat_h_fd = 35786000.0
            _pv_name = config.get('projection_var', 'goes_imager_projection')
            if _pv_name in ds and hasattr(ds[_pv_name], 'perspective_point_height'):
                _sat_h_fd = float(ds[_pv_name].perspective_point_height)
            scale = _sat_h_fd if np.nanmax(np.abs(x_vals)) < 1.0 else 1.0
            half_ext = max(np.nanmax(np.abs(x_vals)),
                           np.nanmax(np.abs(y_vals))) * scale
        else:
            fd_origin   = 'upper'
            need_fliplr = False
            half_ext    = 5568000.0   # fallback: SEVIRI full-disk extent

        render_ss = np.fliplr(raw_ss) if need_fliplr else raw_ss
        ax_global.imshow(render_ss,
                         origin=fd_origin,
                         extent=[-half_ext, half_ext, -half_ext, half_ext],
                         transform=geo_proj,
                         cmap=cmap, vmin=0, vmax=len(labels)-1,
                         interpolation='nearest')
    except Exception as e:
        print(f"  Full-disk render failed: {e}")

    ax_global.set_global()
    bbox = region['bbox']
    ax_global.coastlines(color='cyan', linewidth=0.8)
    ax_global.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.4, edgecolor='gray')
    _corners_c = region['corners'] + [region['corners'][0]]
    lons_b = [c[1] for c in _corners_c]
    lats_b = [c[0] for c in _corners_c]
    ax_global.plot(lons_b, lats_b, color='magenta', linewidth=1,
                   transform=ccrs.PlateCarree(), zorder=10)
    ax_global.set_title(f"{sat} Full Disk View\n{utc_dt.strftime('%Y-%m-%d %H:%M')} UTC",
                        fontsize=14, pad=10)

    # ── region zoom ─────────────────────────────────
    lat_dir = _detect_lat_direction(lats[mask] if hasattr(lats, '__len__') else lats)

    ax_zoom = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
    gx = np.linspace(bbox['west'],  bbox['east'],  1200)
    gy = np.linspace(bbox['south'], bbox['north'], 1200)
    GX, GY = np.meshgrid(gx, gy)
    pts   = np.column_stack((lons[mask], lats[mask]))
    vals  = np.nan_to_num(data_vals.astype(float), nan=mapping['nodata_val'])
    gdata = griddata(pts, vals, (GX, GY), method='nearest')

    img = ax_zoom.imshow(gdata,
                         extent=[bbox['west'], bbox['east'],
                                 bbox['south'], bbox['north']],
                         transform=ccrs.PlateCarree(), cmap=cmap,
                         origin='lower', vmin=0, vmax=len(labels)-1,
                         interpolation='nearest')
    ax_zoom.plot(lons_b, lats_b, color='magenta', linewidth=2,
                 transform=ccrs.PlateCarree(), zorder=10)
    ax_zoom.set_extent([bbox['west'],  bbox['east'],
                        bbox['south'], bbox['north']],
                       crs=ccrs.PlateCarree())
    ax_zoom.add_feature(cfeature.COASTLINE, edgecolor='white', linewidth=2)
    ax_zoom.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='white',
                        linewidth=1.5)
    ax_zoom.add_feature(cfeature.STATES, linestyle='--', edgecolor='gray',
                        alpha=0.5)
    gl = ax_zoom.gridlines(draw_labels=True, linewidth=0.5,
                           color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    ax_zoom.set_title(f"Target Region Analysis\n{address}", fontsize=12, pad=10)

    cbar = plt.colorbar(img, ax=ax_zoom, ticks=range(len(labels)),
                        fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels(labels)
    plt.suptitle(
        f"Global Cloud Forensic Analysis: {sat} | "
        f"{utc_dt.strftime('%Y-%m-%d %H:%M')} UTC",
        fontsize=20, y=0.98, weight='bold')
    save_path = s_dir / f"{sat}_{utc_dt.strftime('%Y%m%d_%H%M')}_forensic.png"
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [Output] Saved: {save_path}")


# ==========================================
# DOWNLOADERS
# ==========================================
def download_aws_data(satellite, utc_dt, download_dir):
    config = SatelliteConfig.SATELLITES[satellite]
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    if 'GOES' in satellite:
        year, doy = utc_dt.year, int(utc_dt.strftime('%j'))
        prefix = config['path_format'].format(
            product=config['product'], year=year,
            day_of_year=doy, hour=utc_dt.hour)
    else:
        prefix = config['path_format'].format(
            product=config['product'], year=utc_dt.year,
            month=utc_dt.month, day=utc_dt.day)
    try:
        print(f"  Searching S3: {config['bucket']} | {prefix}")
        candidates = []
        for page in s3.get_paginator('list_objects_v2').paginate(
                Bucket=config['bucket'], Prefix=prefix):
            if 'Contents' not in page: continue
            for obj in page['Contents']:
                if config.get('file_pattern', '') in obj['Key']:
                    candidates.append(obj['Key'])
        if not candidates: print("  No files found."); return None
        best_key, min_diff = None, float('inf')
        for key in candidates:
            try:
                parts = os.path.basename(key).split('_')
                ts = next((p for p in parts if p.startswith('s') and
                           len(p) >= 12 and p[1].isdigit()), None)
                if ts:
                    fmt = "%Y%j%H%M" if 'GOES' in satellite else "%Y%m%d%H%M"
                    end = 12 if 'GOES' in satellite else 13
                    f_dt = pytz.utc.localize(datetime.strptime(ts[1:end], fmt))
                    diff = abs((f_dt - utc_dt).total_seconds())
                    if diff < min_diff: min_diff, best_key = diff, key
            except: continue
        if not best_key: best_key = candidates[-1]
        local = os.path.join(download_dir, os.path.basename(best_key))
        print(f"  Best: {best_key}  (dt={min_diff/60:.1f} min)")
        if os.path.exists(local): print("  Cached."); return local
        s3.download_file(config['bucket'], best_key, local)
        return local
    except Exception as e:
        print(f"  S3 Error: {e}"); return None


def download_eumdac_data(satellite, utc_dt, download_dir):
    config = SatelliteConfig.SATELLITES[satellite]
    try:
        cred = Path.home() / '.eumdac' / 'credentials'
        if not cred.exists():
            print("  [ERROR] EUMDAC credentials not found."); return None
        creds    = cred.read_text().strip().split(',')
        token    = eumdac.AccessToken((creds[0], creds[1]))
        coll     = eumdac.DataStore(token).get_collection(config['collection_id'])
        if utc_dt.tzinfo is None: utc_dt = pytz.utc.localize(utc_dt)
        products = coll.search(dtstart=utc_dt-timedelta(minutes=45),
                               dtend=utc_dt+timedelta(minutes=45))
        if not products: print("  No products found."); return None
        best, min_diff = None, float('inf')
        for p in products:
            s, e = p.sensing_start, p.sensing_end
            if s.tzinfo is None: s = pytz.utc.localize(s)
            if e.tzinfo is None: e = pytz.utc.localize(e)
            diff = abs((s+(e-s)/2 - utc_dt).total_seconds())
            if diff < min_diff: min_diff, best = diff, p
        if not best: best = products.first()
        print(f"  Selected: {best}  (dt={min_diff/60:.1f} min)")
        out_zip = Path(download_dir) / f"{str(best)}.zip"
        if not out_zip.exists():
            print("  Downloading...")
            with best.open() as fsrc, open(out_zip, 'wb') as fdst:
                shutil.copyfileobj(fsrc, fdst)
        else:
            print("  Zip cached.")
        target = None
        with zipfile.ZipFile(out_zip, 'r') as z:
            for f in z.namelist():
                if f.endswith(config['file_ext']):
                    z.extract(f, download_dir)
                    target = os.path.join(download_dir, f)
                    print(f"  Extracted: {f}"); break
        return target
    except Exception as e:
        print(f"  EUMDAC Error: {e}")
        import traceback; traceback.print_exc()
        return None


# ==========================================
# MAIN
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Multi-Satellite Cloud Analysis")
    # Axis-aligned bbox (original interface)
    parser.add_argument('--north',    type=float, default=None)
    parser.add_argument('--south',    type=float, default=None)
    parser.add_argument('--west',     type=float, default=None)
    parser.add_argument('--east',     type=float, default=None)
    # Tilted polygon (4 corners, overrides north/south/east/west)
    parser.add_argument('--corners',  type=str,   default=None,
        help="4 corner points: 'lat1,lon1 lat2,lon2 lat3,lon3 lat4,lon4' (any winding)")
    parser.add_argument('--date',     type=str,   required=True, help="YYYY-MM-DD")
    parser.add_argument('--time_utc', type=str,   required=True, help="HH:MM")
    parser.add_argument('--plots',    action='store_true')
    args = parser.parse_args()

    if args.corners:
        pairs = args.corners.strip().split()
        if len(pairs) != 4:
            parser.error("--corners requires exactly 4 'lat,lon' pairs")
        corners = [tuple(float(x) for x in p.split(',')) for p in pairs]
    else:
        if any(v is None for v in [args.north, args.south, args.east, args.west]):
            parser.error("Provide either --corners or all of --north/--south/--east/--west")
        if args.north < args.south: args.north, args.south = args.south, args.north
        if args.east  < args.west:  args.east,  args.west  = args.west,  args.east
        corners = [
            (args.north, args.west),
            (args.north, args.east),
            (args.south, args.east),
            (args.south, args.west),
        ]

    region = _build_region(corners)
    bbox   = region['bbox']   # axis-aligned envelope (used for coverage checks + display)
    utc_dt = pytz.utc.localize(
        datetime.strptime(f"{args.date} {args.time_utc}", "%Y-%m-%d %H:%M"))

    import sys
    d_dir_early, s_dir_early = Path("Downloads"), Path("Clouds_coverage_c")
    d_dir_early.mkdir(exist_ok=True)
    s_dir_early.mkdir(exist_ok=True)

    class _Tee:
        def __init__(self, *files): self.files = files
        def write(self, obj):
            for f in self.files: f.write(obj); f.flush()
        def flush(self):
            for f in self.files: f.flush()

    log_path = s_dir_early / f"run_{utc_dt.strftime('%Y%m%d_%H%M')}.log"
    _log_fh  = open(log_path, 'w', buffering=1)
    sys.stdout = _Tee(sys.__stdout__, _log_fh)
    print(f"[LOG] Writing to {log_path}")

    print(f"\n{'='*60}")
    print("GLOBAL MULTI-SATELLITE CLOUD FORENSIC ANALYSIS")
    print(f"{'='*60}")
    print(f"Target : {get_region_details(region)}")
    print(f"Time   : {utc_dt.strftime('%Y-%m-%d %H:%M UTC')}")
    print(f"Region : corners {['({:.3f},{:.3f})'.format(lat,lon) for lat,lon in corners]}")
    print(f"BBox   : {bbox['south']:.3f}S–{bbox['north']:.3f}N  "
          f"{bbox['west']:.3f}W–{bbox['east']:.3f}E")

    if abs(bbox['north']) > 81 or abs(bbox['south']) > 81:
        print("[WARNING] Target near/above Pole — geostationary coverage limited.")

    d_dir, s_dir = d_dir_early, s_dir_early   # already created for logging

    relevant_sats = []
    print(f"\n{'='*60}\nSATELLITE COVERAGE:\n{'='*60}")
    for name, cfg in SatelliteConfig.SATELLITES.items():
        cov     = cfg['coverage']
        overlap = not (bbox['east'] < cov['lon_min'] or bbox['west'] > cov['lon_max'])
        print(f"{name:12} | {cov['lon_min']:6.1f}–{cov['lon_max']:6.1f}E "
              f"| {'YES' if overlap else 'NO'}")
        if overlap: relevant_sats.append(name)

    if not relevant_sats:
        print("[ERROR] No satellites cover this region."); return

    print(f"\nProcessing: {', '.join(relevant_sats)}")
    results = []

    for idx, sat in enumerate(relevant_sats, 1):
        print(f"\n{'='*60}\nSATELLITE {idx}/{len(relevant_sats)}: {sat}\n{'='*60}")
        config = SatelliteConfig.SATELLITES[sat]

        print("[Step 1] Downloading...")
        path = (download_aws_data(sat, utc_dt, str(d_dir))
                if config['data_source'] == 'aws'
                else download_eumdac_data(sat, utc_dt, str(d_dir)))
        if not path:
            print(f"  [SKIP] No data for {sat}"); continue

        try:
            print(f"[Step 2] Opening: {os.path.basename(path)}")
            ds = xr.open_dataset(path, engine=config['engine'])

            # ----------------------------------------------------------
            # MSG-2 / MSG-3: GRIB SEVIRI — coordinate decode + scatter stats
            # (identical structure; only sub-satellite longitude differs)
            # ----------------------------------------------------------
            if sat in ('MSG-2', 'MSG-3'):
                _sat_lon = SatelliteConfig.SATELLITES[sat].get('sat_longitude', _SAT_LON)
                print(f"[Step 3] Decoding {sat} SEVIRI coordinates from GRIB keys "
                      f"(lon_0={_sat_lon}°)...")
                v_lats, v_lons, v_vals = get_msg2_scatter(ds, sat_lon=_sat_lon,
                                                         sat_label=sat)

                print("[Step 4] Extracting region (ground-truth scatter)...")
                b_lats, b_lons, b_vals, b_weights = extract_msg2_bbox(
                    v_lats, v_lons, v_vals, region, sat_label=sat)

                if len(b_vals) == 0:
                    print("  [SKIP] Zero pixels in region"); ds.close(); continue

                print("[Step 5] Centroid analysis...")
                analyze_centroid(region, b_lats, b_lons, b_vals, sat)

                print("[Step 6] Statistics (raw scatter, no interpolation)...")
                stats = calculate_comprehensive_statistics(b_vals, sat,
                                                           weights_in=b_weights)

                _clat = region['centroid_lat']
                _clon = region['centroid_lon']
                _vza  = compute_viewing_zenith_angle(_clat, _clon, _sat_lon)
                print(f"  [VZA] {sat} viewing zenith angle to bbox centroid: {_vza:.1f}°"
                      + (f"  *** UNRELIABLE: VZA>{VZA_UNRELIABLE_THRESHOLD:.0f}° ***"
                         if _vza > VZA_UNRELIABLE_THRESHOLD else ""))

                results.append({
                    'satellite':      sat,
                    'cloud_pct':      stats['cloud_pct'],
                    'clear_pct':      stats['clear_pct'],
                    'total_pixels':   stats['total_pixels'],
                    'nodata_pct':     stats.get('nodata_pct', 0),
                    'dqf_masked_pct': 0.0,   # MSG GRIB: no DQF variable; value=3 encodes no-data
                    'vza_deg':        _vza,
                })

                if args.plots:
                    print("[Step 7] Plotting full-disk + region...")
                    vals2d_raw = ds['p260537'].values.reshape(_NJ, _NI).astype(float)
                    plot_msg2_dual_view(
                        sat, region, get_region_details(region), utc_dt, s_dir,
                        vals2d_raw,
                        b_lats, b_lons, b_vals,
                        sat_lon=_sat_lon)

            # ----------------------------------------------------------
            # All other satellites
            # ----------------------------------------------------------
            else:
                print("[Step 3] Extracting region...")
                weights, lats, lons = extract_bounding_box_data(sat, ds, region)
                if weights is None or not np.any(weights > 0):
                    print(f"  [SKIP] No pixels in region"); ds.close(); continue

                sel       = weights > 0
                data_vals = ds[config['variable']].values[sel].astype(float)
                pix_w     = weights[sel]

                # DQF masking: set bad-quality pixels to NaN so they count as no-data
                _dqf_masked_pct = 0.0
                dqf_var = config.get('dqf_var')
                if dqf_var:
                    _dqf_candidates = [dqf_var, 'DQF', 'QualityFlag', 'quality_flag']
                    _dqf_found = next((v for v in _dqf_candidates if v in ds), None)
                    if _dqf_found:
                        dqf_vals = ds[_dqf_found].values[sel]
                        dqf_bad  = dqf_vals >= 2   # 0=Good, 1=Conditionally usable, ≥2=invalid
                        n_dqf    = int(dqf_bad.sum())
                        n_total  = len(data_vals)
                        _dqf_masked_pct = n_dqf / n_total * 100 if n_total > 0 else 0.0
                        print(f"  [DQF] Using '{_dqf_found}': "
                              f"{n_dqf:,}/{n_total:,} pixels masked (DQF>=2, "
                              f"{_dqf_masked_pct:.1f}% of region bbox — treated as no-data)")
                        data_vals[dqf_bad] = np.nan
                    else:
                        print(f"  [DQF] Variable '{dqf_var}' not found in dataset — skipping DQF filter")

                print("[Step 4] Centroid analysis...")
                analyze_centroid(region, lats[sel], lons[sel], data_vals, sat)

                print("[Step 5] Statistics...")
                stats = calculate_comprehensive_statistics(data_vals, sat,
                                                           weights_in=pix_w)

                _sat_lon_vza = config.get('sat_longitude', 0.0)
                _clat = region['centroid_lat']
                _clon = region['centroid_lon']
                _vza  = compute_viewing_zenith_angle(_clat, _clon, _sat_lon_vza)
                print(f"  [VZA] {sat} viewing zenith angle to region centroid: {_vza:.1f}°"
              + (f"  *** UNRELIABLE: VZA>{VZA_UNRELIABLE_THRESHOLD:.0f}° — pixel distortion "
                 f"too high for reliable cloud fraction ***" if _vza > VZA_UNRELIABLE_THRESHOLD else ""))

                results.append({
                    'satellite':      sat,
                    'cloud_pct':      stats['cloud_pct'],
                    'clear_pct':      stats['clear_pct'],
                    'total_pixels':   stats['total_pixels'],
                    'nodata_pct':     stats.get('nodata_pct', 0),
                    'dqf_masked_pct': _dqf_masked_pct,
                    'vza_deg':        _vza,
                })

                if args.plots:
                    print("[Step 6] Plotting...")
                    plot_split_view(sat, ds, sel, lats, lons, data_vals,
                                    region, get_region_details(region), utc_dt, s_dir)

            ds.close()

        except Exception as e:
            print(f"  [ERROR] {e}")
            import traceback; traceback.print_exc()

    # ---- Summary -------------------------------------------------------
    print(f"\n{'='*80}\nMULTI-SATELLITE SUMMARY\n{'='*80}")
    print(f"Target : {get_region_details(region)}")
    print(f"Time   : {utc_dt.strftime('%Y-%m-%d %H:%M UTC')}")

    if results:
        # Assign priority order based on VZA + no-data penalty
        assign_priority_order(results)   # sorts results in-place; adds priority_order

        # Reliable = no-data ≤ 30%
        reliable = [r for r in results if r.get('nodata_pct', 0) <= 30]
        if not reliable:
            reliable = results
            print("  [NOTE] All results have >30% No-data; using all anyway.")

        print(f"\n{'Pri':<4} {'Satellite':<12} {'Cloud%':<8} {'Clear%':<8} "
              f"{'Pixels':<10} {'NoData%':<9} {'VZA°'}")
        print('-'*60)
        for r in results:
            nd  = r.get('nodata_pct', 0)
            vza = r.get('vza_deg', 90.0)
            pri = r.get('priority_order', '-')
            print(f"#{pri:<3} {r['satellite']:<12} {r['cloud_pct']:<8.2f} "
                  f"{r['clear_pct']:<8.2f} {r['total_pixels']:<10.1f} "
                  f"{nd:<9.1f} {vza:.1f}")
        print('-'*60)

        p1 = results[0]
        print(f"\n>>> PRIORITY ASSESSMENT <<<")
        print(f"Primary source : {p1['satellite']}  "
              f"VZA={p1['vza_deg']:.1f}°  nodata={p1['nodata_pct']:.1f}%")
        print(f"Cloud coverage : {p1['cloud_pct']:.2f}%")
        print(f"Clear coverage : {p1['clear_pct']:.2f}%")
        print(f"Effective pixels: {p1['total_pixels']:.1f}")

        # Export to GeoJSON
        export_to_geojson(region, results, utc_dt, s_dir)
    else:
        print("[ERROR] No successful analyses.")

    print(f"\n{'='*80}\nDone. Outputs in '{s_dir}'.\n{'='*80}")


def export_to_geojson(region, results, utc_dt, s_dir):
    """
    Export analysis results as GeoJSON FeatureCollection.
    Includes the analysis polygon (from region corners) and centroid point.
    """
    clat = region['centroid_lat']
    clon = region['centroid_lon']

    serializable_results = []
    for r in results:
        nd_pct  = float(r.get('nodata_pct', 0))
        dqf_pct = round(float(r.get('dqf_masked_pct', 0.0)), 2)
        vza     = float(r.get('vza_deg', 0.0))
        vza_ok  = vza <= VZA_UNRELIABLE_THRESHOLD
        status  = "Reliable" if (nd_pct <= 30.0 and vza_ok) else (
                  f"Unreliable (VZA={vza:.1f}°>{VZA_UNRELIABLE_THRESHOLD:.0f}°)" if not vza_ok
                  else "Unreliable (>30% No-Data)")
        if not vza_ok and nd_pct > 30.0:
            status = f"Unreliable (VZA={vza:.1f}°>{VZA_UNRELIABLE_THRESHOLD:.0f}° AND >30% No-Data)"
        if dqf_pct > 0:
            status += f"; {dqf_pct:.1f}% excluded by DQF>=2"
        serializable_results.append({
            "satellite":          str(r['satellite']),
            "priority_order":     int(r.get('priority_order', 3)),
            "cloud_pct":          float(r['cloud_pct']),
            "clear_pct":          float(r['clear_pct']),
            "effective_pixels":   round(float(r['total_pixels']), 1),
            "nodata_pct":         nd_pct,
            "dqf_masked_pct":     dqf_pct,
            "dqf_note":           (f"{dqf_pct:.2f}% of bbox pixels had DQF>=2 "
                                   f"(out-of-range/no-value) and were excluded from analysis "
                                   f"— counted as no-data in priority selection"
                                   if dqf_pct > 0 else "All pixels DQF-good (0 or 1)"),
            "vza_deg":            round(vza, 1),
            "vza_reliable":       vza_ok,
            "vza_threshold_deg":  VZA_UNRELIABLE_THRESHOLD,
            "status":             status,
        })

    # Region Polygon Feature (uses actual corner points, not axis-aligned bbox)
    corners_closed = region['corners'] + [region['corners'][0]]
    bbox_feature = {
        "type": "Feature",
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[float(lon), float(lat)] for lat, lon in corners_closed]]
        },
        "properties": {
            "feature_type": "Analysis Region",
            "analysis_time_utc": utc_dt.isoformat(),
            "satellites": serializable_results
        }
    }

    # Centroid Point Feature
    centroid_feature = {
        "type": "Feature",
        "geometry": {
            "type": "Point",
            "coordinates": [float(clon), float(clat)]
        },
        "properties": {
            "feature_type": "Region Centroid",
            "analysis_time_utc": utc_dt.isoformat()
        }
    }

    # Write FeatureCollection
    output_path = s_dir / f"analysis_{utc_dt.strftime('%Y%m%d_%H%M')}.geojson"
    fc = {
        "type": "FeatureCollection",
        "features": [bbox_feature, centroid_feature]
    }
    
    with open(output_path, 'w') as f:
        json.dump(fc, f, indent=2)
    
    print(f"  [Output] GeoJSON saved: {output_path}")


if __name__ == "__main__":
    main()