# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Shared pieces of the two earthkit-data fieldlist test modules of this folder.

This module is not a test file. ``test_regrid_fieldlist_grib.py`` and
``test_regrid_fieldlist_array.py`` split the use cases of the old
``test_regrid_fieldlist.py`` between them -- the GRIB fields and fieldlists in one file,
their non-GRIB (array) counterparts in the other -- and what both need is here: the
fieldlist loader :func:`create_test_fieldlist`, the 10x10 reference check
:func:`check_regrid_to_ll`, and the case table :data:`TO_LL_CASES` with the test body
:func:`run_to_ll` each file runs as one parametrised test with its own field type. The
GRIB-only parts (the round-trip and points metadata constants, the ecCodes grid-spec
probe) stayed in ``test_regrid_fieldlist_grib.py``.

The folder has no ``__init__.py``: pytest puts it on ``sys.path``, and the two test files
import this module with ``from regrid_fieldlist_common import ...``, like
``regrid_common``. Everything the fieldlist tests share with the rest of the folder -- the
``mir`` / ``earthkit-data`` skip marks, the ``interpolation=`` kwarg spellings and the
loading of the legacy npz references -- is in ``regrid_common``, not duplicated here.

The module has no side effects at import and imports without earthkit-data installed:
``from_source`` is guarded by ``NO_EKD``, and the tests of both files skip through
``regrid_common.requires_ekd`` before any function here is called.
"""

import numpy as np
import pytest
from regrid_common import legacy_values

from earthkit.geo import regrid
from earthkit.geo.utils.testing import NO_EKD, compare_global_ll_results, get_test_data_path

if not NO_EKD:
    from earthkit.data import from_source

# Metadata keys compared between the input and the regridded result. The regular lat-lon
# input case also compares the gridType of the input message; for the Gaussian input the
# gridType of the result is compared on its own ("regular_ll", regridding a reduced
# Gaussian grid to a regular one).
METADATA_KEYS = ["param", "level", "date", "time"]
METADATA_KEYS_GRID_TYPE = [*METADATA_KEYS, "gridType"]

# Direction increments of the 10x10 global lat-lon output grid.
_INCREMENT_REF = {"iDirectionIncrementInDegrees": 10.0, "jDirectionIncrementInDegrees": 10.0}

# The cases of the ``..._to_ll`` test of both files: the input file, the stem of its
# reference file, the metadata keys to compare, whether a single field (instead of the
# whole fieldlist) is regridded, and whether the gridType of the result is also compared
# with "regular_ll". Each file runs them for its own field type.
TO_LL_CASES = [
    # a regular lat-lon fieldlist, the gridType of the input takes part in the metadata comparison
    pytest.param(("5x5.grib", "out_5x5_10x10", METADATA_KEYS_GRID_TYPE, False, False), id="reg_ll"),
    # an octahedral Gaussian fieldlist, its result compared against the regular lat-lon gridType
    pytest.param(("O32.grib", "out_O32_10x10", METADATA_KEYS, False, True), id="gg"),
    # the same Gaussian input, regridding a single field instead of the whole fieldlist
    pytest.param(("O32.grib", "out_O32_10x10", METADATA_KEYS, True, True), id="single_field"),
]


def create_test_fieldlist(filename, subfolder="global_0_360", field_type="grib"):
    """The earthkit-data fieldlist of *filename*, as a GRIB or an array fieldlist.

    The input is read from the test-data URL with ``get_test_data_path``; *field_type*
    "grib" keeps the GRIB fieldlist, "array" converts it with ``to_fieldlist()``.
    """
    ds = from_source("url", get_test_data_path(filename, subfolder=subfolder)).to_fieldlist()
    if field_type == "array":
        return ds.to_fieldlist()
    elif field_type == "grib":
        return ds
    else:
        raise ValueError(f"Unknown field type: {field_type}")


def check_regrid_to_ll(result, v_ref, interpolation, metadata_keys, metadata_ref, is_field, check_grid_type):
    """Check a result regridded to the 10x10 global lat-lon grid against its reference.

    Parameters
    ----------
    result : field or fieldlist
        The regridded field, or a fieldlist holding one field.
    v_ref : numpy.ndarray
        Reference values of the regridded field.
    interpolation : str
        Name of the interpolation used, for the reference comparison.
    metadata_keys : list of str
        Metadata keys compared with ``metadata_ref``.
    metadata_ref : object
        The same metadata, read from the input before regridding.
    is_field : bool
        Whether ``result`` is a single field (instead of a fieldlist).
    check_grid_type : bool
        Whether the gridType of the result is also compared with "regular_ll".
    """
    if is_field:
        fields = [result]
    else:
        assert len(result) == 1
        fields = list(result)

    for f in fields:
        assert f.shape == (19, 36)
        compare_global_ll_results(f.to_numpy(), v_ref, interpolation, rtol=1e-4)
        for k, v in _INCREMENT_REF.items():
            assert np.isclose(f.metadata(k), v), k

    assert result.metadata(list(metadata_keys)) == metadata_ref

    if check_grid_type:
        assert result.metadata("gridType") == ("regular_ll" if is_field else ["regular_ll"])


def run_to_ll(case, kwarg, interpolation, field_type):
    """Run one :data:`TO_LL_CASES` case of one field type and check it.

    Parameters
    ----------
    case : tuple
        One entry of :data:`TO_LL_CASES`: ``(input_file, ref_stem, metadata_keys,
        is_field, check_grid_type)``.
    kwarg : dict
        The ``interpolation=`` spelling under test, as an ``earthkit.geo.regrid`` kwargs
        dict (empty for the default); an entry of ``regrid_common.KWARG_ALIASES``.
    interpolation : str
        The canonical name that ``kwarg`` must select, for the reference comparison.
    field_type : str
        "grib" or "array": the kind of fieldlist :func:`create_test_fieldlist` builds.
    """
    input_file, ref_stem, metadata_keys, is_field, check_grid_type = case

    ds = create_test_fieldlist(input_file, field_type=field_type)
    source = ds[0] if is_field else ds

    v_ref = legacy_values(f"{ref_stem}_{interpolation}.npz")
    metadata_ref = source.metadata(list(metadata_keys))

    r = regrid(source, out_grid={"grid": [10, 10]}, **kwarg)

    check_regrid_to_ll(r, v_ref, interpolation, metadata_keys, metadata_ref, is_field, check_grid_type)
