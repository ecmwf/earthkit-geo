# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


_NE_DATASET_NAME_KEYS = {
    "admin_0_map_units": "NAME_LONG",
    "admin_0_countries": "NAME_LONG",
    "admin_1_states_provinces": "name_en",
}

_NE_RESOLUTIONS = {
    "10m": 10e6,
    "50m": 50e6,
    "110m": 110e6,
}


def multipolygon_to_coordinates(multipolygon):
    """
    Convert a shapely MultiPolygon object to a list of coordinates.

    Parameters
    ----------
    multipolygon : shapely.geometry.MultiPolygon
        A MultiPolygon object to convert to a list of coordinates.

    Returns
    -------
    list
        A list of lists of coordinates, where each sublist represents a polygon
        in the combined geometry of the specified countries.
    """
    coordinates = []
    for polygon in multipolygon.geoms:
        exterior_coords = list(polygon.exterior.coords)
        coordinates.append([[x, y] for y, x in exterior_coords])

    return coordinates


def _closest_resolution(resolution):
    """
    Get the closest available resolution from Natural Earth to the specified
    resolution in metres.

    Parameters
    ----------
    resolution : float
        The desired resolution, in metres, to find the closest available
        resolution to.
    """
    return min(
        _NE_RESOLUTIONS,
        key=lambda res: (abs(_NE_RESOLUTIONS[res] - resolution), -_NE_RESOLUTIONS[res]),
    )


def country_polygons(country_names, resolution=110e6):
    """
    Get the combined geometry of one or more countries by name from Natural
    Earth's shapefiles.

    Parameters
    ----------
    country_names : str or list of str
        The name(s) of the country or countries to get the geometry for.
    resolution : float, optional
        The desired resolution, in metres, of the shapefile to use. Will return
        the closest available resolution from Natural Earth. Default is 110e6.

    Returns
    -------
    list
        A list of lists of coordinates, where each sublist represents a polygon
        in the combined geometry of the specified countries.
    """
    import cartopy.io.shapereader as shpreader
    from shapely.geometry import MultiPolygon
    from shapely.ops import unary_union

    ne_resolution = _closest_resolution(resolution)

    if isinstance(country_names, str):
        country_names = [country_names]

    country_names = [name.lower() for name in country_names]
    geometries = []

    for source, attribute in _NE_DATASET_NAME_KEYS.items():
        shpfilename = shpreader.natural_earth(
            resolution=ne_resolution,
            category="cultural",
            name=source,
        )
        reader = shpreader.Reader(shpfilename)
        for record in reader.records():
            name = record.attributes.get(attribute) or ""
            name = name.replace("\x00", "").lower()
            if name in country_names:
                geometries.append(record.geometry)
                if len(geometries) == len(country_names):
                    break
        if len(geometries) == len(country_names):
            break

    # Check for any missing countries
    missing_countries = [
        name
        for name in country_names
        if name not in [record.attributes.get(attribute).lower() for record in reader.records()]
    ]
    if missing_countries:
        raise ValueError(
            f"No countries or states named {missing_countries} found in Natural " "Earth's shapefiles"
        )

    combined_geometry = unary_union(geometries)

    if combined_geometry.geom_type == "Polygon":
        combined_geometry = MultiPolygon([combined_geometry])

    return multipolygon_to_coordinates(combined_geometry)
