# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

try:
    # NOTE: the `version.py` file must not be present in the git repository
    #   as it is generated by setuptools at install time
    from ._version import __version__
except ImportError:  # pragma: no cover
    # Local copy or not installed with setuptools
    __version__ = "999"


from earthkit.geo.distance import (
    GeoKDTree,
    haversine_distance,
    nearest_point_haversine,
    nearest_point_kdtree,
)
from earthkit.geo.figure import IFS_SPHERE, UNIT_SPHERE

__all__ = [
    "GeoKDTree",
    "haversine_distance",
    IFS_SPHERE,
    "nearest_point_haversine",
    "nearest_point_kdtree",
    UNIT_SPHERE,
    "__version__",
]
