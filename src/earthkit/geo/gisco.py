# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import os
import warnings
from pathlib import Path
from zipfile import ZipFile

import requests

GISCO_RESOLUTIONS = ["01M", "03M", "10M", "20M", "60M"]

GISCO_GEOMETRY_TYPES = {
    "polygons": "RG",
    "lines": "BN",
    "points": "LB",
}

GISCO_URL_TEMPLATE = "https://gisco-services.ec.europa.eu" "/distribution/v2/{category}/shp"

AVAILABLE_YEARS = {
    "nuts": [2024, 2021, 2016, 2013, 2010, 2006, 2003],
    "countries": [2024, 2020, 2016, 2013, 2010, 2006, 2001],
}


_DEFAULT_DATA_DIR = Path.home() / ".cache" / "earthkit-geo"
DATA_DIR = Path(os.environ.get("EARTHKIT_GEO_DATA_DIR", _DEFAULT_DATA_DIR))


def _get_gisco_cache_dir():
    """Return (and create if missing) the base directory where GISCO shapefiles will be stored."""
    cache_dir = DATA_DIR / "gisco-data"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def _load_with_earthkit_data(path):
    """Load shapefile with earthkit.data if available, warn otherwise."""
    try:
        import earthkit.data
    except ImportError:
        warnings.warn("The 'earthkit.data' package is not installed; returning the file path instead.")
        return path
    else:
        return earthkit.data.from_source("file", path)


class DownloadWarning(Warning):
    """Warning for download operations."""


def _download_and_cache_gisco(name, category, geometry_type, resolution, suffix="", year=2024):
    """
    Download (if needed) and unzip the GISCO shapefile for a given
    country code `name` and resolution. Returns the path to the .shp file.
    """
    # Convert geometry_type if needed
    if geometry_type not in GISCO_GEOMETRY_TYPES.values():
        try:
            geometry_type = GISCO_GEOMETRY_TYPES[geometry_type]
        except KeyError:
            raise ValueError(
                f"Invalid geometry type '{geometry_type}'. "
                f"Valid types are: {list(GISCO_GEOMETRY_TYPES.keys())}"
            )

    # Validate year
    if year not in AVAILABLE_YEARS.get(category, []):
        raise ValueError(
            f"Year {year} is not available for category '{category}'. "
            f"Available years: {AVAILABLE_YEARS.get(category, [])}"
        )

    # Only include resolution if not LB
    resolution_str = resolution if geometry_type != "LB" and resolution else ""
    filename_parts = [name, geometry_type]
    if resolution_str:
        filename_parts.append(resolution_str)
    filename_parts.append(str(year))
    filename_parts.append("4326")  # CRS
    if suffix:
        filename_parts.append(suffix.lstrip("_"))
    filename = "_".join(filename_parts) + ".shp.zip"

    # Set up paths
    cache_dir = _get_gisco_cache_dir()
    zip_path = cache_dir / filename
    extract_dir = cache_dir / "_".join([name, resolution_str]) if resolution_str else cache_dir / name
    shp_path = extract_dir / filename.replace(".zip", "")

    # Download and extract if needed
    if not shp_path.exists():
        from filelock import FileLock

        lock = str(zip_path) + ".lock"
        with FileLock(lock):
            # Check again, another thread/process may have created the file
            if not shp_path.exists():
                if not zip_path.exists():
                    url = f"{GISCO_URL_TEMPLATE.format(category=category)}/{filename}"
                    warnings.warn(f"Downloading: {url}", DownloadWarning)
                    resp = requests.get(url)
                    resp.raise_for_status()
                    zip_path.write_bytes(resp.content)
                extract_dir.mkdir(parents=True, exist_ok=True)
                with ZipFile(zip_path, "r") as zf:
                    zf.extractall(extract_dir)

    return str(shp_path)


def nuts_regions(level, resolution="10M", year=2024, geometry_type="polygons"):
    """
    Retrieve NUTS regions of the given level.

    Please see https://ec.europa.eu/eurostat/web/nuts for more information.

    Parameters
    ----------
    level : int
        The NUTS level (0, 1, 2, or 3).
    resolution : str, optional
        The desired resolution of the shapefile. Must be one of the following:
        '01M', '03M', '10M', '20M', or '60M'. Default is '10M'.
    geometry_type : str, optional
        The type of geometry to retrieve. Must be one of 'polygons', 'lines', or 'points'.
        Default is 'polygons'.
    year : int, optional
        The year of the NUTS data to retrieve. Default is 2024 (latest available).
    """
    suffix = f"_LEVL_{level}"
    path = _download_and_cache_gisco(
        name="NUTS",
        category="nuts",
        geometry_type=geometry_type,
        resolution=resolution,
        suffix=suffix,
        year=year,
    )
    return _load_with_earthkit_data(path)


def countries(resolution="10M", year=2024, geometry_type="polygons"):
    """
    Retrieve country boundaries.

    Parameters
    ----------
    resolution : str, optional
        The desired resolution of the shapefile. Must be one of the following:
        '01M', '03M', '10M', '20M', or '60M'. Default is '10M'.
    geometry_type : str, optional
        The type of geometry to retrieve. Must be one of 'polygons', 'lines', or 'points'.
        Default is 'polygons'.
    year : int, optional
        The year of the country data to retrieve. Default is 2024 (latest available).
    """
    path = _download_and_cache_gisco(
        name="CNTR",
        category="countries",
        geometry_type=geometry_type,
        resolution=resolution,
        year=year,
    )
    return _load_with_earthkit_data(path)
