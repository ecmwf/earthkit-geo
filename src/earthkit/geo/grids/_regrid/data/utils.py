# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import logging

LOG = logging.getLogger(__name__)


def create_grid_object(grid):
    """Convert a grid or grid spec into a ``Grid`` object.

    Parameters
    ----------
    grid : eckit.geo.Grid, earthkit.geo.grids.Grid, dict, str or None
        The grid or grid spec to convert. Left unchanged if already a
        ``Grid`` instance.

    Returns
    -------
    eckit.geo.Grid, earthkit.geo.grids.Grid or None
        ``grid`` unchanged if it is already a ``Grid`` instance, None if
        ``grid`` is None, otherwise a new ``Grid`` built from the
        (normalised) grid spec.

    Raises
    ------
    Exception
        Whatever ``eckit.geo.Grid`` raises if ``grid`` cannot be parsed
        (logged before being re-raised).
    """
    if grid is None:
        return None

    from eckit.geo import Grid

    from earthkit.geo.grids import Grid as EarthKitGrid

    if isinstance(grid, (Grid, EarthKitGrid)):
        return grid
    else:
        try:
            return Grid(normalise_grid_spec(grid))
        except Exception as e:
            LOG.exception(f"Cannot create Grid object from grid_spec: {e}")
            raise


def get_grid_spec(grid):
    """Return the grid spec dict for a grid or grid spec.

    Parameters
    ----------
    grid : eckit.geo.Grid, earthkit.geo.grids.Grid, dict, str or None
        The grid or grid spec to inspect.

    Returns
    -------
    dict or None
        None if ``grid`` is None, ``grid.spec`` if ``grid`` is a ``Grid``
        instance, otherwise the normalised grid spec (see
        :func:`normalise_grid_spec`).
    """
    if grid is None:
        return None

    from eckit.geo import Grid

    from earthkit.geo.grids import Grid as EarthKitGrid

    if isinstance(grid, (Grid, EarthKitGrid)):
        return grid.spec
    else:
        return normalise_grid_spec(grid)


def normalise_grid_spec(grid_spec):
    """Return a normalised grid spec.

    Parameters
    ----------
    grid_spec : Any
        The grid spec to normalise. Left untouched if not a dict.

    Returns
    -------
    Any
        ``grid_spec`` unchanged if it is not a dict, otherwise a shallow
        copy with its ``"reference"`` values (if any) converted to float.

    Raises
    ------
    ValueError
        If a "reference" value is a str that cannot be converted to float.

    Notes
    -----
        The "reference" values are converted to float when they are given as str.

    """
    if not isinstance(grid_spec, dict):
        return grid_spec

    reference = grid_spec.get("reference", None)
    if isinstance(reference, (list, tuple)):
        grid_spec = grid_spec.copy()
        r = []
        for i, v in enumerate(reference):
            if isinstance(v, str):
                try:
                    v = float(v)
                except ValueError:
                    raise ValueError(
                        f"Invalid value={v!r} at index={i} in grid_spec['reference']={reference!r}. "
                        "Cannot be converted to float."
                    )
            r.append(v)
        grid_spec["reference"] = r
    return grid_spec


def create_unstructured_grid_object(lats, lons):
    import numpy as np
    from eckit.geo import Grid

    lats = np.asarray(lats).flatten()
    lons = np.asarray(lons).flatten()
    return Grid({"latitudes": lats, "longitudes": lons})


def _is_proj_source(projection) -> bool:
    """Whether a grid spec projection is described by a PROJ source string."""
    if not isinstance(projection, dict) or projection.get("type") != "proj":
        return False
    return isinstance(projection.get("source"), str)


def scatter_projected_grid(grid):
    """Convert a PROJ-source projected grid into the unstructured lat/lon grid of its points.

    Parameters
    ----------
    grid : eckit.geo.Grid, earthkit.geo.grids.Grid, dict, str or None
        The grid to scatter, as a ``Grid`` object, a plain grid spec dict or a
        JSON string. A ``regular_xy`` grid whose projection is a PROJ source
        string, such as ``{"type": "proj", "source": "+proj=laea ..."}`` (the
        way the xarray loader spells a projection detected from a CF
        ``proj4_params`` attribute), is converted into the unstructured lat/lon
        grid of its points. A plain dict spec of that kind, with or without its
        top-level ``"type"``, is turned into a ``Grid`` first.

    Returns
    -------
    eckit.geo.Grid, earthkit.geo.grids.Grid, dict, str or None
        ``grid.to_unstructured_ll()`` for a PROJ-source projected grid,
        otherwise ``grid`` itself, unchanged and the same object: any other
        grid, dict spec or JSON string is passed through as it is, and None
        stays None.

    Raises
    ------
    Exception
        Whatever ``eckit.geo.Grid`` raises if a dict PROJ spec cannot be
        parsed into a grid (see :func:`create_grid_object`).

    Notes
    -----
    The conversion keeps the point order of the projected grid, i.e. its x axis
    scanned fastest, which is the order data laid out on the grid is in, so the
    values still match their points. Backends that cannot parse a ``+proj=...``
    source string as an input grid can use this to take the points as a grid
    they do read.
    """
    from eckit.geo import Grid

    original = grid
    if isinstance(grid, dict) and _is_proj_source(grid.get("projection")):
        # a plain spec of the PROJ kind: eckit decides which grid type it describes
        grid = create_grid_object(grid)

    if isinstance(grid, Grid) and grid.type == "regular_xy" and _is_proj_source(grid.spec.get("projection")):
        return grid.to_unstructured_ll()
    # anything not scattered is handed back as it came in, even a dict spec that was parsed above
    return original
