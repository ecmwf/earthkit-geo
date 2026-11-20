# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""earthkit-data non-GRIB (array) fieldlist regridding tests (``earthkit.geo.regrid``) against MIR references.

The non-GRIB half of the old ``test_regrid_fieldlist.py``; its GRIB half is
``test_regrid_fieldlist_grib.py``. What the two files share -- the fieldlist loader, the
10x10 reference check and the ``..._to_ll`` case table and test body -- is
``regrid_fieldlist_common.py``, while the shared pieces they have with the rest of the
folder, the kwarg alias parametrisation, the ``mir`` / ``earthkit-data`` skip marks and
the loading of the legacy npz references, still come from ``regrid_common``.

The use cases are the ones of the old file, unchanged: the same two inputs -- a regular
lat-lon and an octahedral Gaussian one -- regridded to a 10x10 global lat-lon grid as a
fieldlist and as a single field (the ``field_type="array"`` cases of the parametrised
``..._to_ll`` test), and a field made with ``set()``, which is decoupled from the GRIB
message it came from. The inputs are still GRIB files read from the test-data URL with
``get_test_data_path``, converted to array fields by
``regrid_fieldlist_common.create_test_fieldlist`` with ``to_fieldlist()``, so what these
tests exercise is the regridding of fields that hold no GRIB message.

The reference comparisons deliberately keep the strict ``compare_global_ll_results``
helper of ``earthkit.geo.utils.testing`` (through ``regrid_fieldlist_common``) instead of
the registry compare of ``regrid_common.get_interpolation(...).compare`` (used by the grid
matrix tests): these checks keep the tolerance and the pole-row handling they always had.
"""

import numpy as np
import pytest
from regrid_common import KWARG_ALIASES, legacy_values, requires_ekd, requires_mir
from regrid_fieldlist_common import TO_LL_CASES, create_test_fieldlist, run_to_ll

from earthkit.geo import regrid
from earthkit.geo.grids import Grid


@requires_mir
@requires_ekd
@pytest.mark.parametrize("_kwarg,interpolation", KWARG_ALIASES)
@pytest.mark.parametrize("case", TO_LL_CASES)
def test_regrid_fieldlist_array_to_ll(case, _kwarg, interpolation):
    run_to_ll(case, _kwarg, interpolation, field_type="array")


@requires_mir
@requires_ekd
@pytest.mark.download
@pytest.mark.tmp_cache
def test_regrid_fieldlist_array_single_field_non_grib():
    interpolation = "linear"
    out_grid = {"grid": [10, 10]}

    fl = create_test_fieldlist("5x5.grib", field_type="array")
    field = fl[0].set({"parameter.variable": "msl"}, labels={"my_label": "my_value"})

    # the field is now decoupled from the original grib message
    assert field.message() is None

    v_ref = legacy_values(f"out_5x5_10x10_{interpolation}.npz")

    r = regrid(field, out_grid=out_grid, interpolation=interpolation)

    assert r.geography.shape() == (19, 36)
    assert np.allclose(r.values, v_ref)
    assert r.get("parameter.variable") == "msl"
    assert r.get("labels.my_label") == "my_value"

    grid_ref = Grid({"grid": [10, 10]}).spec

    assert r.geography.grid_spec() == grid_ref
