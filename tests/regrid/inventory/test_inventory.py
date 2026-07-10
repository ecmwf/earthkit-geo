# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import pytest

from earthkit.geo.grids._regrid.backends.db import SYS_DB
from earthkit.geo.grids._regrid.gridspec import _GridSpec


@pytest.mark.parametrize(
    "gs_in, gs_out, docs_in, docs_out",
    [
        (
            {"grid": "O32"},
            {"grid": [0.1, 0.1]},
            {"grid": "O32", "_type": {"reduced_gg", "reduced-gg"}, "octahedral": True},
            {"grid": [0.1, 0.1], "_type": {"regular_ll", "regular-ll"}},
        ),
        (
            {"grid": "H8"},
            {"grid": [0.5, 0.5]},
            {"grid": "H8", "_type": {"healpix", "HEALPix"}, "order": "ring"},
            {"grid": [0.5, 0.5], "_type": {"regular_ll", "regular-ll"}},
        ),
        (
            {"grid": "H8", "order": "nested"},
            {"grid": [0.5, 0.5]},
            {"grid": "H8", "_type": {"healpix", "HEALPix"}, "order": "nested"},
            {"grid": [0.5, 0.5], "_type": {"regular_ll", "regular-ll"}},
        ),
    ],
)
def test_inventory_build(gs_in, gs_out, docs_in, docs_out):
    """Test objects and methods used to build the matrix inventory documentation for regridding."""
    r = SYS_DB.find_entry(gs_in, gs_out, "linear")

    assert r

    r_in = _GridSpec.from_dict(r["input"]).inventory_docs_spec
    r_out = _GridSpec.from_dict(r["output"]).inventory_docs_spec

    assert r_in
    for k, v in docs_in.items():
        assert k in r_in
        if isinstance(v, set):
            assert r_in[k] in v
        else:
            assert r_in[k] == v

    for k, v in docs_out.items():
        assert k in r_out
        if isinstance(v, set):
            assert r_out[k] in v
        else:
            assert r_out[k] == v
