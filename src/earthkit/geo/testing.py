# (C) Copyright 2020 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import logging
import os

from earthkit.geo.wind import _normalise_longitude

LOG = logging.getLogger(__name__)

_ROOT_DIR = top = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
)
if not os.path.exists(os.path.join(_ROOT_DIR, "tests", "data")):
    _ROOT_DIR = "./"


def earthkit_test_data_file(*args):
    return os.path.join(_ROOT_DIR, "tests", "data", *args)


def normalise_lon(lon):
    if isinstance(lon, (int, float)):
        return _normalise_longitude(lon, -180)
    else:
        import numpy as np

        lon = np.asarray(lon)
        return np.array([_normalise_longitude(x, -180) for x in lon])
