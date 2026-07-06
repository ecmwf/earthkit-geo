# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


from earthkit.geo.utils.testing import earthkit_test_data_path


def tests_eccodes_grid_spec():
    import os

    os.environ["ECCODES_ECKIT_GEO"] = "1"
    import eccodes

    fname = earthkit_test_data_path("o32.grib2")

    with open(fname, "rb") as f:
        while True:
            handle = eccodes.codes_new_from_file(f, eccodes.CODES_PRODUCT_GRIB)
            if handle is not None:
                v = eccodes.codes_get(handle, "gridSpec")
                print(f"gridSpec={v}", flush=True)
                assert isinstance(v, str)
                break
