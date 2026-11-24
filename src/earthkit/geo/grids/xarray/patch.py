# (C) Copyright 2024 Anemoi contributors.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import logging
from typing import Any

import xarray as xr

LOG = logging.getLogger(__name__)


def patch_attributes(ds: xr.Dataset, attributes: dict[str, dict[str, Any]]) -> xr.Dataset:
    """Patch the attributes of the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to patch.
    attributes : Dict[str, Dict[str, Any]]
        The attributes to patch.

    Returns
    -------
    Any
        The patched dataset.
    """
    for name, value in attributes.items():
        variable = ds[name]
        variable.attrs.update(value)

    return ds


def patch_coordinates(ds: xr.Dataset, coordinates: list[str]) -> xr.Dataset:
    """Patch the coordinates of the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to patch.
    coordinates : List[str]
        The coordinates to patch.

    Returns
    -------
    Any
        The patched dataset.
    """
    for name in coordinates:
        ds = ds.assign_coords({name: ds[name]})

    return ds


def patch_rename(ds: xr.Dataset, renames: dict[str, str]) -> xr.Dataset:
    """Rename variables in the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to patch.
    renames : dict[str, str]
        Mapping from old variable names to new variable names.

    Returns
    -------
    Any
        The patched dataset.
    """
    return ds.rename(renames)


def patch_sort_coordinates(ds: xr.Dataset, sort_coordinates: list[str]) -> xr.Dataset:
    """Sort the coordinates of the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to patch.
    sort_coordinates : List[str]
        The coordinates to sort.

    Returns
    -------
    Any
        The patched dataset.
    """
    for name in sort_coordinates:
        ds = ds.sortby(name)
    return ds


def patch_subset_dataset(ds: xr.Dataset, selection: dict[str, Any]) -> xr.Dataset:
    """Select a subset of the dataset using xarray's sel method.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to patch.
    selection : dict[str, Any]
        Dictionary mapping dimension names to selection criteria.
        Keys must be existing dimension names in the dataset.
        Values can be any type accepted by xarray's sel method, including:
        - Single values (int, float, str, datetime)
        - Lists or arrays of values
        - Slices (using slice() objects)
        - Boolean arrays

    Returns
    -------
    xr.Dataset
        The patched dataset containing only the selected subset.

    Examples
    --------
    >>> # Select specific time and pressure level
    >>> patch_subset_dataset(ds, {"time": "2020-01-01", "pressure": 500})

    >>> # Select a range using slice
    >>> patch_subset_dataset(ds, {"lat": slice(-90, 90), "lon": slice(0, 180)})
    """
    ds = ds.sel(selection)

    return ds


def patch_latlon(ds: xr.Dataset, patch: dict[str, dict[str, Any]]) -> Any:
    """Patch the latitude and longitude coordinates of the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to patch.
    patch : dict[str, dict[str, Any]]
        The patch to apply to the latitude and longitude coordinates.

    Returns
    -------
    Any
        The patched dataset.
    """
    # Implement the patching logic for latitude and longitude here
    if (
        "latitude" not in ds.coords
        and "latitude" in ds.data_vars
        and "longitude" not in ds.coords
        and "longitude" in ds.data_vars
    ):
        ds = patch_coordinates(ds, ["latitude", "longitude"])
    return ds


PATCHES = {
    "attributes": patch_attributes,
    "coordinates": patch_coordinates,
    "latlon": patch_latlon,
    "rename": patch_rename,
    "sort_coordinates": patch_sort_coordinates,
    "subset_dataset": patch_subset_dataset,
}


def patch_dataset(ds: xr.Dataset, patch: dict[str, dict[str, Any]]) -> Any:
    """Patch the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to patch.
    patch : dict[str, dict[str, Any]]
        The patch to apply.

    Returns
    -------
    Any
        The patched dataset.
    """
    ORDER = [
        "coordinates",
        "attributes",
        "rename",
        "sort_coordinates",
        "subset_dataset",
        "latlon",
    ]
    for what, values in sorted(patch.items(), key=lambda x: ORDER.index(x[0])):
        if what not in PATCHES:
            raise ValueError(f"Unknown patch type {what!r}")

        ds = PATCHES[what](ds, values)

    return ds
