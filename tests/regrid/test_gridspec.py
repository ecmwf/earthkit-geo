# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


from earthkit.geo.grids._regrid.gridspec import normalise_grid_spec


def test_gridspec_normaliser_1():
    gs = {"area": [71.9834273, -25, 25, 49.9834833], "grid": [0.0166667, 0.0166667], "reference": ["5e-05", 0.0166167]}

    r = normalise_grid_spec(gs)
    assert r == {
        "area": [71.9834273, -25, 25, 49.9834833],
        "grid": [0.0166667, 0.0166667],
        "reference": [5e-05, 0.0166167],
    }


def test_gridspec_normaliser_2():
    gs = {"area": [71.9834273, -25, 25, 49.9834833], "grid": [0.0166667, 0.0166667], "reference": [5e-05, 0.0166167]}

    r = normalise_grid_spec(gs)
    assert r == {
        "area": [71.9834273, -25, 25, 49.9834833],
        "grid": [0.0166667, 0.0166667],
        "reference": [5e-05, 0.0166167],
    }
