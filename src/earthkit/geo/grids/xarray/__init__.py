# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""xarray geography: coordinate classification, grid inference and variable discovery.

Public package holding the machinery that turns the coordinates and CF
metadata of an ``xarray.Dataset`` into geographical grids. It classifies each
coordinate (see :mod:`~.coordinates` and the guessers in :mod:`~.guesser`),
builds the grid a variable lives on (see :mod:`~.grid` and the grid specs in
:mod:`~.grid_spec`), and discovers the geographical variables of a dataset
with their grids (see :func:`~.loader.variables`).

It is used by the regrid xarray data handler (see
:mod:`earthkit.geo.grids._regrid.data.xarray.handler`) and is meant to be
imported by other earthkit packages. Unlike the rest of ``earthkit.geo``,
importing it does import ``xarray``; it never imports ``earthkit.data``.
"""

from .coordinates import (
    Coordinate,
    LatitudeCoordinate,
    LongitudeCoordinate,
    PointCoordinate,
    ScalarCoordinate,
    UnsupportedCoordinate,
    XCoordinate,
    YCoordinate,
    datetime64_to_datetime,
    extract_single_value,
    is_scalar,
    timedelta64_to_timedelta,
)
from .grid import (
    Grid,
    LatLonGrid,
    MeshedGrid,
    MeshedXYGrid,
    MeshProjectionGrid,
    ProjectionGrid,
    RawXYGrid,
    ReducedGaussianGrid,
    UnstructuredGrid,
    UnstructuredProjectionGrid,
    UnstructuredXYGrid,
    XarrayGrid,
    XYGrid,
)
from .guesser import CoordinateAttributes, CoordinateGuesser, DefaultCoordinateGuesser, RulesCoordinateGuesser
from .loader import variables

__all__ = [
    "Coordinate",
    "CoordinateAttributes",
    "CoordinateGuesser",
    "DefaultCoordinateGuesser",
    "Grid",
    "LatLonGrid",
    "LatitudeCoordinate",
    "LongitudeCoordinate",
    "MeshProjectionGrid",
    "MeshedGrid",
    "MeshedXYGrid",
    "PointCoordinate",
    "ProjectionGrid",
    "RawXYGrid",
    "ReducedGaussianGrid",
    "RulesCoordinateGuesser",
    "ScalarCoordinate",
    "UnstructuredGrid",
    "UnstructuredProjectionGrid",
    "UnstructuredXYGrid",
    "UnsupportedCoordinate",
    "XCoordinate",
    "XarrayGrid",
    "XYGrid",
    "YCoordinate",
    "datetime64_to_datetime",
    "extract_single_value",
    "is_scalar",
    "timedelta64_to_timedelta",
    "variables",
]
