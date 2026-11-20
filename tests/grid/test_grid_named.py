#!/usr/bin/env python3

# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""``Grid`` object tests for the named grids of the eckit catalogues.

Every grid is built from a gridspec dict which names it (``Grid({"grid":
"O32"})``, ``Grid({"grid": "ORCA2_T"})``, ``Grid({"grid": "0.5/0.5"})``) and
checked against expected values worked out independently of the code: from the
grid formulas in ``docs/concepts/grid/types/*.rst``, from the arithmetic of the
catalogue entries themselves (a ``[start, end, step]`` range, ``n`` cells per
90 degrees), or from the data a catalogue entry declares (the shape, uid and
bounding box of an ICON, ORCA or FESOM grid). Where the two disagree, that is a
finding: the case is marked BROKEN below, never "fixed" by copying what the
code returns. The plain string form, ``Grid("O32")``, is kept in one small
test of its own, :func:`test_string_form_is_the_same_grid`.

The names which are executed come from the written-out tables in this file: the
allowlist of what the probe covered
(``.ecmwf-multi-agents/notes/grid-probe.md``). Nothing is discovered from the
catalogue for execution. The one test which reads the catalogue files,
:func:`test_every_catalogue_name_is_reviewed`, only checks that every catalogue
name is either tested here or listed in :data:`NOT_TESTED` with a reason, so a
new name fails until it has been reviewed.

Two rules keep the run safe and quick:

* ``grid.order`` is never read for the types which abort the interpreter: ORCA,
  ICON and FESOM here, ``sh`` and ``unstructured_ll`` in the other files of this
  folder. The crash is an uncaught C++ ``eckit::NotImplemented``: SIGABRT, which
  Python cannot catch, so those cases are ``xfail(run=False)`` and are never
  executed (finding F1).
* Coordinates are only materialised (``to_latlons()``) below
  :data:`MAX_POINTS` points. ``CAMA-01min`` has 233,280,000 points, which is
  gigabytes of Python floats. Above the gate only metadata is checked, which the
  catalogues serve offline; that is a resource limit, not a failure, so those
  names are simply not in the points list.

Everything which does not work is
``pytest.param(..., marks=pytest.mark.xfail(run=False, reason="BROKEN: ..."))``:
it shows up as an expected failure, nothing is skipped in silence. The module is
marked ``long_test``: a default run skips it, ``-E long`` runs it.
"""

import math

import pytest
import yaml

from earthkit.geo.grids import Grid
from earthkit.geo.sphinxext import generate_catalogue_rst as catalogue

pytestmark = pytest.mark.long_test

# to_latlons() and any other coordinate-materialising call is only used below
# this many points.
MAX_POINTS = 1_000_000

# The tolerance for a value of a gridspec: eckit writes the numbers of a spec
# with about 15 significant digits, so a derived value is never bit-identical.
REL = 1e-9
ABS = 1e-12

# bounding_box() of a grid which is not regional is the global placeholder
# (north, west, south, east), whatever the points reach (finding F5).
GLOBAL = (90.0, 0.0, -90.0, 360.0)

# How far a projected grid's points may stick out of its bounding_box(). The box
# is made of the corners of the projected grid while its largest latitude is on
# an edge: swisslv03 exceeds its own box by 4e-7 degrees (a few centimetres).
AREA_TOL = 1e-4

# Sentinel for a grid whose ``order`` aborts the interpreter: never read it.
CRASH = "not implemented"

# The abort messages, one per type (finding F1), as BROKEN reasons.
CRASH_ORCA = "BROKEN: F1 grid.order aborts the interpreter (SIGABRT): uncaught eckit::NotImplemented: Not implemented"
CRASH_ICON = (
    "BROKEN: F1 grid.order aborts the interpreter (SIGABRT): uncaught eckit::NotImplemented: "
    "Not implemented: Grid: order() is not implemented for type 'ICON'"
)
CRASH_FESOM = (
    "BROKEN: F1 grid.order aborts the interpreter (SIGABRT): uncaught eckit::NotImplemented: "
    "Not implemented: Grid: order() is not implemented for type 'FESOM'"
)


# ---------------------------------------------------------------------------
# the expected values, from the grid formulas and the catalogue data
# ---------------------------------------------------------------------------


def octahedral_points(n):
    """Points of the octahedral reduced Gaussian grid ``O<n>``: 4n² + 36n.

    docs/concepts/grid/types/reduced_gg.rst: the latitude i counted from the pole
    has 4n + 2i - 2 points, which adds up to 4n² + 36n over the 2n latitudes.
    """
    return 4 * n * n + 36 * n


def linear_points(n):
    """Points of the classical ("linear") reduced Gaussian grid ``N<n>``.

    These grids have no closed formula: the points per latitude come from a
    table, and only some ``n`` are in it (see
    :func:`test_classical_grid_names_start_at_n16`). The values are the standard
    counts of the ECMWF linear grids, the N16/N32/N64 grids used by the shallow
    water test case community.
    """
    return {16: 1688, 32: 6114, 64: 23112, 1280: 8505906}[n]


def regular_gg_shape(n):
    """Shape of the regular Gaussian grid ``F<n>``: 2n latitudes, 4n longitudes."""
    return (2 * n, 4 * n)


def healpix_points(nside):
    """Pixels of a HEALPix grid: 12 * nside²."""
    return 12 * nside * nside


def regular_ll_shape(dlon, dlat):
    """Shape of the global ``<dlon>/<dlat>`` grid: poles included, 0 E once."""
    return (round(180 / dlat) + 1, round(360 / dlon))


def range_points(triple):
    """Number of values of a catalogue ``[start, end, step]`` range."""
    start, end, step = triple
    return round((end - start) / step) + 1


def xy_order(x, y):
    """Scan order of a regular_xy grid given by ``[start, end, step]`` ranges.

    ``i`` follows x and ``j`` follows y, and ``+`` means increasing. The default
    order, ``i+j-``, is left out of a canonical spec.
    """
    return ("i+" if x[2] > 0 else "i-") + ("j+" if y[2] > 0 else "j-")


def hydro_cells(n):
    """Cell size in degrees of an arakawa_c grid with ``n`` cells per 90 degrees."""
    return 90 / n


def global_hydro_spec(n):
    """Canonical spec of the global arakawa_c grids of hydro.yaml (``west: -180``).

    The points are the centres of a 2n by 4n field of cells: 2n latitude rows
    from half a cell inside the north pole, 4n columns, so the area of the grid
    keeps the poles and wraps once past 180 in longitude.
    """
    d = hydro_cells(n)
    return {"area": [90, -180 + d / 2, -90, 180 + d / 2], "grid": [d, d], "reference": [d / 2, d / 2]}


def area_hydro_spec(n, north, west, south, east):
    """Canonical spec of a regional arakawa_c grid given its area.

    The area is the outer edge of the border cells, so the grid of cell centres
    covers it shrunk by half a cell on each side, with (north - south) / d rows
    and (east - west) / d columns.
    """
    d = hydro_cells(n)
    return {
        "area": [north - d / 2, west + d / 2, south + d / 2, east - d / 2],
        "grid": [d, d],
        "reference": [d / 2, d / 2],
    }


def hydro_shape(n, north, west, south, east):
    """Shape of the regional arakawa_c grid of :func:`area_hydro_spec`."""
    d = hydro_cells(n)
    return (round((north - south) / d), round((east - west) / d))


class Named:
    """One grid that can be named in a gridspec, and what it should look like."""

    def __init__(
        self,
        name,
        grid_type,
        shape,
        spec,
        order,
        bbox=GLOBAL,
        area_bbox=False,
        crash="",
        broken_points="",
    ):
        self.name = name
        self.grid_type = grid_type
        self.shape = tuple(shape)
        # The canonical gridspec the name resolves to. Compared with _same(),
        # numbers with a tolerance.
        self.spec = spec
        # The scan order, or CRASH when reading grid.order aborts the process.
        self.order = order
        # The expected bounding_box(), or None for a regional grid whose area is
        # not known independently: then only its invariants are checked.
        self.bbox = bbox
        # True when bounding_box() is the area of the grid, not the globe.
        self.area_bbox = area_bbox
        # The BROKEN reason for the order() abort of this grid's type.
        self.crash = crash
        # The BROKEN reason of a to_latlons() failure, empty when it works.
        self.broken_points = broken_points

    @property
    def size(self):
        """The expected number of points: the product of the expected shape."""
        return math.prod(self.shape)

    @property
    def points(self):
        """Whether to_latlons() may be called, and is expected to work."""
        return self.size < MAX_POINTS and not self.broken_points

    def __str__(self):
        return self.name


def orca_grids():
    """The ORCA grids of the catalogue: every base name, every arrangement.

    The shapes are the ``shape`` fields of ORCA.yaml; the five arrangements of a
    base share one. A catalogue name always carries the arrangement suffix:
    ``ORCA2`` on its own is not a name.
    """
    bases = {
        "ORCA2": (149, 182),
        "ORCA1": (292, 362),
        "eORCA1": (332, 362),
        "ORCA025": (1021, 1442),
        "eORCA025": (1207, 1442),
        "ORCA12": (3059, 4322),
        "eORCA12": (3606, 4322),
    }
    grids = []
    for base, shape in bases.items():
        for arrangement in "FTUVW":
            name = f"{base}_{arrangement}"
            grids.append(
                Named(
                    name=name,
                    grid_type="ORCA",
                    shape=shape,
                    spec={"grid": name},
                    order=CRASH,
                    crash=CRASH_ORCA,
                )
            )
    return grids


def icon_grids():
    """The ICON grids of the catalogue: name, points, uid and, for MeteoSwiss, area.

    The points and uids are the ``shape`` and ``uid`` fields of ICON-CH.yaml,
    ICON-MPIM.yaml and ICON-DWD.yaml. An ICON canonical spec is the name with the
    arrangement suffix ``_C``, the catalogue uid, and the catalogue bounding_box
    as ``area`` where one is declared. The 5th field is the name the canonical
    spec uses when it is not the name which was asked for: the ``-V1`` names of
    ICON-CH are aliases, and resolve to the grid of the plain name.
    """
    ch = (50.6, -0.9, 41.9, 17.8)
    entries = [
        # MeteoSwiss: the one family with a bounding_box in the catalogue.
        ("ICON-CH1", 1147980, "17643da2574959b644d254a3cd6e2bc0", ch, None),
        ("ICON-CH1-V1", 1147980, "17643da2574959b644d254a3cd6e2bc0", ch, "ICON-CH1"),
        ("ICON-CH2", 283876, "bbbd5a09855499243c7a4aa4c8762920", ch, None),
        ("ICON-CH2-V1", 283876, "bbbd5a09855499243c7a4aa4c8762920", ch, "ICON-CH2"),
        # Max Planck Institute for Meteorology: all 7 names.
        ("icon-grid-0016-r02b09-o", 14886338, "375cb0cc637e11e89d6f8f41a9b9ff4b", None, None),
        ("icon-grid-0020-r02b05-o", 59261, "e96d69e2ab5811e89e299566f2744a3c", None, None),
        ("icon-grid-0024-r02b07-o", 933398, "28f3a5e2de9c11e88555cbf786e364f5", None, None),
        ("icon-grid-0036-r02b04-o", 15105, "5bd948e8ac1a11eaa6b1d317264fdca9", None, None),
        ("icon-grid-0038-r02b11-o", 237316813, "9b323b62670411eb85afdd9149101e1a", None, None),
        ("icon-grid-0040-r02b10-o", 59359799, "99c03a3c657811eb8e42f565c1ad0089", None, None),
        ("icon-grid-0046-r02b06-o", 235403, "f4ed57f6b2ea11e9ae92c52a3fa37d96", None, None),
        # DWD: a written-out sample of 14 of the 56 names, covering every
        # icon_type class, every suffix family (-g -r -n -n02 -l -lr -ln02), the
        # base refinements r02, r03 and r19, the smallest grid of the catalogue
        # and grids on both sides of the resource gate. The other 42 names are in
        # NOT_TESTED.
        ("icon-grid-0001-r02b05-r", 81920, "eee9e7881dd111b2aedbc760324270de", None, None),
        ("icon-grid-0002-r02b06-g", 327680, "fc9f45081dd111b2aedbc760324270de", None, None),
        ("icon-grid-0004-r02b07-g", 1310720, "8bb13be81dd211b2bdd32f117dd15b92", None, None),
        ("icon-grid-0009-r02b03-r", 5120, "65dadcba1dcf11b2880c0f1645f3d1dc", None, None),
        ("icon-grid-0010-r02b04-g", 20480, "66a341d21dcf11b2880c0f1645f3d1dc", None, None),
        ("icon-grid-0020-r03b06-g", 737280, "365801d21dd611b2b98a6b5c1380b5e3", None, None),
        ("icon-grid-0022-r03b07-g", 2949120, "eb8c0fa91dd611b2837c051661b31a74", None, None),
        ("icon-grid-0027-r03b08-n02", 659156, "ec13b8bcb82d11e4b13f4d55411d42e6", None, None),
        ("icon-grid-0031-r02b06-n02", 19744, "a668eacafe2e11e49a47e3e5371a87a7", None, None),
        ("icon-grid-0042-r02b10-l", 253520, "d378d127082af71232e62a2093c92ba0", None, None),
        ("icon-grid-0043-r19b06-lr", 134541, "ac1190a0a846e31092a72a2bf62329c0", None, None),
        ("icon-grid-0045-r19b08-ln02", 716980, "b258b3391435078ecc92d51b8a3f2ac0", None, None),
        ("icon-grid-0054-r02b04-g", 20480, "87397f993e6f913d00fc361367e725a0", None, None),
        ("icon-grid-0055-r02b05-n", 2656, "e234e01a8556e9a84bcb42361d2f24e0", None, None),
    ]

    grids = []
    for name, points, uid, box, alias in entries:
        spec = {"grid": f"{alias or name}_C", "uid": uid}
        if box is not None:
            spec = {"area": list(box), **spec}
        grids.append(
            Named(
                name=name,
                grid_type="ICON",
                shape=(points,),
                spec=spec,
                order=CRASH,
                bbox=box if box is not None else GLOBAL,
                area_bbox=box is not None,
                crash=CRASH_ICON,
            )
        )
    return grids


def fesom_grids():
    """The FESOM grids of FESOM.yaml: the cell and node arrangements of 4 meshes.

    The shapes are the ``shape`` fields of the catalogue.
    """
    meshes = {
        "pi": (5839, 3140),
        "CORE2": (244659, 126858),
        "DART": (6262485, 3160340),
        "NG5": (14741520, 7402886),
    }
    grids = []
    for mesh, (cell, node) in meshes.items():
        for arrangement, points in (("C", cell), ("N", node)):
            name = f"{mesh}_{arrangement}"
            grids.append(
                Named(
                    name=name,
                    grid_type="FESOM",
                    shape=(points,),
                    spec={"grid": name},
                    order=CRASH,
                    crash=CRASH_FESOM,
                )
            )
    return grids


# The x, y and false ranges of grid.yaml, as the catalogue file gives them.
LAEA_EFAS_X = (0, 4995000, 5000)
LAEA_EFAS_Y = (0, -4745000, -5000)
LAEA_EFAS_FALSE = (1818500.021, -2287499.955)
SMUFF_X = (0, 3798000, 2000)
SMUFF_Y = (0, -4398000, -2000)
SMUFF_FALSE = (1949000, -2101000)
SWISSLV95_X = (2439000, 2867000, 1000)
SWISSLV95_Y = (1040500, 1333500, 1000)
SWISSLV03_X = (255500, 964500, 1000)
SWISSLV03_Y = (-159500, 479500, 1000)


def grid_yaml_grids():
    """The four regional projected grids of grid.yaml.

    Each entry gives x and y as ``[start, end, step]`` and a projection: the
    shape is (values of y, values of x), and the canonical spec is the catalogue
    entry, with an ``order`` key only where it is not the default ``i+j-``. The
    bounding box follows from the projection, which is not spelled out here: no
    value is known independently, so its invariants are checked, and in the
    points test that it holds the points. Note that the ``false`` key of the
    catalogue is a boolean to YAML, and eckit keeps it as ``false`` in the spec,
    not ``false_xy`` (finding F12).
    """
    grids = []
    entries = [
        (
            "LAEA-EFAS-5km",
            LAEA_EFAS_X,
            LAEA_EFAS_Y,
            dict(false=LAEA_EFAS_FALSE),
            dict(type="laea", lat_0=52, lon_0=10, figure="grs80", reverse=True),
        ),
        (
            "SMUFF-OPERA-2km",
            SMUFF_X,
            SMUFF_Y,
            dict(false=SMUFF_FALSE),
            dict(type="laea", lat_0=55, lon_0=10, figure="wgs84", reverse=True),
        ),
        ("swisslv95", SWISSLV95_X, SWISSLV95_Y, dict(), "swisslv95"),
        ("swisslv03", SWISSLV03_X, SWISSLV03_Y, dict(), "swisslv03"),
    ]
    for name, x, y, extra, projection in entries:
        spec = {"x": list(x), "y": list(y), "projection": projection}
        spec.update(extra)
        order = xy_order(x, y)
        if order != "i+j-":
            spec["order"] = order
        grids.append(
            Named(
                name=name,
                grid_type="regular_xy",
                shape=(range_points(y), range_points(x)),
                spec=spec,
                order=order,
                bbox=None,
                area_bbox=True,
            )
        )
    return grids


# The catalogue area of EFAS-5, and the point counts it gives over cells of
# 1/60 degree: 49.5° of latitude is 2970 rows, 75.5° of longitude is 4530 columns.
EFAS_AREA = (72.25, -25.25, 22.75, 50.25)

# The catalogue area of GLOFAS-4.
GLOFAS_AREA = (89.975, -179.975, -59.975, 179.975)


def hydro_grids():
    """The seven hydrological grids of hydro.yaml, all of type arakawa_c.

    arakawa_c is resolved into a cell-centred grid, and the canonical spec says
    ``regular_ll`` (finding F12). The five CAMA grids are global: the catalogue
    gives ``west: -180`` and no area. EFAS-5 covers its catalogue area.
    """
    grids = []
    cama = {"CAMA-15min": 360, "CAMA-06min": 900, "CAMA-05min": 1080, "CAMA-03min": 1800, "CAMA-01min": 5400}
    for name, n in cama.items():
        spec = global_hydro_spec(n)
        grids.append(
            Named(
                name=name,
                grid_type="regular_ll",
                shape=(2 * n, 4 * n),
                spec=spec,
                order="i+j-",
                bbox=tuple(spec["area"]),
            )
        )

    spec = area_hydro_spec(5400, *EFAS_AREA)
    grids.append(
        Named(
            name="EFAS-5",
            grid_type="regular_ll",
            shape=hydro_shape(5400, *EFAS_AREA),
            spec=spec,
            order="i+j-",
            bbox=tuple(spec["area"]),
        )
    )

    # GLOFAS-4 (n=1800, so cells of 0.05 degrees) reaches the north edge within
    # half a cell and its east edge 0.025 degrees short of 180. The builder
    # treats the latitude range as the extent of the points, which are then
    # 89.975 down to -59.975, 3000 rows of them, and wraps the longitudes once
    # past 180, 7200 columns. That snapping is not documented anywhere, so these
    # values are not derived from a formula: they are written out as the area of
    # the catalogue entry implies them, and the report notes the convention.
    grids.append(
        Named(
            name="GLOFAS-4",
            grid_type="regular_ll",
            shape=(3000, 7200),
            spec={
                "area": [90, -179.975, -59.975, 180.025],
                "grid": [0.05, 0.05],
                "reference": [0.025, 0.025],
            },
            order="i+j-",
            bbox=(90.0, -179.975, -59.975, 180.025),
        )
    )
    return grids


def pattern_grids():
    """The names built from a pattern: O<N>, N<N>, F<N>, H<N>, <dlon>/<dlat>."""
    grids = []

    for n in (8, 16, 32, 64, 128, 1280):
        grids.append(
            Named(
                name=f"O{n}",
                grid_type="reduced_gg",
                shape=(octahedral_points(n),),
                spec={"grid": f"O{n}"},
                order="i+j-",
            )
        )

    for n in (16, 32, 64, 1280):
        grids.append(
            Named(
                name=f"N{n}",
                grid_type="reduced_gg",
                shape=(linear_points(n),),
                spec={"grid": f"N{n}"},
                order="i+j-",
            )
        )

    for n in (16, 32, 48, 1600):
        grids.append(
            Named(
                name=f"F{n}",
                grid_type="regular_gg",
                shape=regular_gg_shape(n),
                spec={"grid": f"F{n}"},
                order="i+j-",
            )
        )

    for n in (4, 8, 12, 16, 31, 1024):
        grids.append(
            Named(
                name=f"H{n}",
                grid_type="HEALPix",
                shape=(healpix_points(n),),
                spec={"grid": f"H{n}"},
                order="ring",
            )
        )

    # HR<n> spells out the ring order, and canonicalises to the plain H<n>.
    grids.append(
        Named(
            name="HR16",
            grid_type="HEALPix",
            shape=(healpix_points(16),),
            spec={"grid": "H16"},
            order="ring",
        )
    )
    # HN<n> is the nested order, which the canonical spec keeps as an order key.
    grids.append(
        Named(
            name="HN16",
            grid_type="HEALPix",
            shape=(healpix_points(16),),
            spec={"grid": "H16", "order": "nested"},
            order="nested",
        )
    )
    # HN10: a nested grid whose Nside is not a power of 2 builds and reports its
    # metadata, but its points cannot be made (finding F4).
    grids.append(
        Named(
            name="HN10",
            grid_type="HEALPix",
            shape=(healpix_points(10),),
            spec={"grid": "H10", "order": "nested"},
            order="nested",
            broken_points=(
                "BROKEN: F4 to_latlons() raises RuntimeError: Assertion failed: "
                "Nside must be a power of 2, for a grid which builds and reports 1200 pixels"
            ),
        )
    )

    for dlon, dlat in ((1, 1), (0.5, 0.5), (1.5, 0.5), (2, 2), (0.25, 0.25)):
        grids.append(
            Named(
                name=f"{dlon}/{dlat}",
                grid_type="regular_ll",
                shape=regular_ll_shape(dlon, dlat),
                spec={"grid": [dlon, dlat]},
                order="i+j-",
            )
        )

    return grids


ORCA_GRIDS = orca_grids()
ICON_GRIDS = icon_grids()
FESOM_GRIDS = fesom_grids()
GRID_YAML_GRIDS = grid_yaml_grids()
HYDRO_GRIDS = hydro_grids()
PATTERN_GRIDS = pattern_grids()

NAMED = tuple(ORCA_GRIDS + ICON_GRIDS + FESOM_GRIDS + GRID_YAML_GRIDS + HYDRO_GRIDS + PATTERN_GRIDS)

# The names of the catalogues which this file builds: the allowlist of
# :func:`test_every_catalogue_name_is_reviewed`. The pattern names are not in it,
# they are not in a catalogue file.
CATALOGUE_TESTED = {grid.name for grid in ORCA_GRIDS + ICON_GRIDS + FESOM_GRIDS + GRID_YAML_GRIDS + HYDRO_GRIDS}

# The catalogue names which are not tested, each with its reason. A name which
# is neither tested nor listed here fails
# :func:`test_every_catalogue_name_is_reviewed`.
NOT_TESTED = {
    name: "not probed: ICON-DWD is covered by the 14-name sample of ICON_GRIDS"
    for name in (
        "icon-grid-0003-r02b06-r",
        "icon-grid-0005-r03b06-r",
        "icon-grid-0006-r03b07-g",
        "icon-grid-0007-r02b04-r",
        "icon-grid-0008-r02b05-g",
        "icon-grid-0011-r02b03-r",
        "icon-grid-0012-r02b04-g",
        "icon-grid-0013-r02b04-r",
        "icon-grid-0014-r02b05-g",
        "icon-grid-0015-r02b05-r",
        "icon-grid-0016-r02b06-g",
        "icon-grid-0017-r02b06-r",
        "icon-grid-0018-r02b07-g",
        "icon-grid-0019-r03b05-r",
        "icon-grid-0021-r03b06-r",
        "icon-grid-0023-r02b05-r",
        "icon-grid-0024-r02b06-g",
        "icon-grid-0025-r03b06-r",
        "icon-grid-0026-r03b07-g",
        "icon-grid-0028-r02b07-n02",
        "icon-grid-0029-r02b04-r",
        "icon-grid-0030-r02b05-g",
        "icon-grid-0032-r03b04-r",
        "icon-grid-0033-r03b05-g",
        "icon-grid-0034-r03b06-n02",
        "icon-grid-0035-r03b05-r",
        "icon-grid-0036-r03b06-g",
        "icon-grid-0037-r03b07-n02",
        "icon-grid-0038-r02b06-r",
        "icon-grid-0039-r02b07-g",
        "icon-grid-0040-r02b08-n02",
        "icon-grid-0041-r02b09-lr",
        "icon-grid-0044-r19b07-l",
        "icon-grid-0046-r19b06-lr",
        "icon-grid-0047-r19b07-l",
        "icon-grid-0048-r03b07-lr",
        "icon-grid-0049-r03b08-l",
        "icon-grid-0050-r02b07-n02",
        "icon-grid-0051-r03b07-n02",
        "icon-grid-0052-r03b07-lr",
        "icon-grid-0053-r03b08-l",
        "icon-grid-0056-r02b05-n",
    )
}

# Invalid names, each with why it is invalid. What each of them raises is the
# finding (F6): ``SpecError: cannot build grid without 'type'``, as if the
# gridspec had no type, instead of an error naming the unknown grid name. The
# rejection itself is right — none of these is a grid — the message is not, so
# it is recorded as a BROKEN case which is never executed, never asserted as
# correct behaviour.
INVALID_NAMES = {
    # Not catalogue names: an ORCA name always carries the arrangement suffix.
    "ORCA2": "no ORCA entry of this name; the suffix _F/_T/_U/_V/_W is part of it",
    "eORCA1": "no ORCA entry of this name; the suffix _F/_T/_U/_V/_W is part of it",
    "ORCA025": "no ORCA entry of this name; the suffix _F/_T/_U/_V/_W is part of it",
    # The name quoted in the plan; ICON-CH.yaml has ICON-CH1, ICON-CH2 and their
    # -V1 aliases, not this.
    "ICON-CH1-EPS": "not in any catalogue; ICON-CH.yaml has ICON-CH1, ICON-CH2 and the -V1 aliases",
    # Cited in docs/concepts/grid/types/projected.rst but not installed with eckit.
    "meps-bris-nordic": "documented in projected.rst but absent from the installed catalogue",
    "wibble-42": "not a grid of any catalogue or name pattern",
    # Matches the HEALPix pattern H<N>, but Nside 0 has no pixels.
    "H0": "an H<N> name with Nside 0, which is not a grid",
}

NAME_ERROR_PARAMS = [
    pytest.param(
        name,
        id=name,
        marks=pytest.mark.xfail(
            run=False,
            reason=f'BROKEN: F6 misleading error: {name!r} ({why}) is rejected with "cannot '
            "build grid without 'type',\" as if the gridspec had no type; the expected "
            "error is a clear \"unknown grid name\" which names the value of 'grid'",
        ),
    )
    for name, why in INVALID_NAMES.items()
]

# Representative names, one per family, for which the plain string still builds
# the same grid as the gridspec dict (the tests above all use the dict).
STRING_FORM_NAMES = (
    "O32",
    "N32",
    "F32",
    "H8",
    "HR16",
    "HN16",
    "0.5/0.5",
    "ORCA2_T",
    "ICON-CH1-V1",
    "icon-grid-0009-r02b03-r",
    "pi_C",
    "swisslv95",
    "LAEA-EFAS-5km",
    "EFAS-5",
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _same(actual, expected, path="spec"):
    """Compare a value of a gridspec with the expected one, numbers with a tolerance."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict), f"{path}: expected a mapping, got {actual!r}"
        assert sorted(map(str, actual)) == sorted(map(str, expected)), (
            f"{path}: keys {sorted(map(str, actual))} != {sorted(map(str, expected))}"
        )
        for key in expected:
            _same(actual[key], expected[key], f"{path}[{key!r}]")
    elif isinstance(expected, (list, tuple)):
        assert isinstance(actual, (list, tuple)), f"{path}: expected a sequence, got {actual!r}"
        assert len(actual) == len(expected), f"{path}: {len(actual)} values, expected {len(expected)}"
        for i, (got, want) in enumerate(zip(actual, expected)):
            _same(got, want, f"{path}[{i}]")
    elif isinstance(expected, bool):
        assert actual is expected, f"{path}: {actual!r} != {expected!r}"
    elif isinstance(expected, (int, float)):
        assert actual == pytest.approx(expected, rel=REL, abs=ABS), f"{path}: {actual!r} != {expected!r}"
    else:
        assert actual == expected, f"{path}: {actual!r} != {expected!r}"


def _assert_area(bbox, name):
    """bounding_box() of a regional grid, for which no value is known independently."""
    north, west, south, east = bbox
    assert 90.0 >= north > south >= -90.0, f"{name}: latitudes of the area {bbox}"
    assert -180.0 <= west < east <= 360.0, f"{name}: longitudes of the area {bbox}"
    assert tuple(bbox) != GLOBAL, f"{name}: bounding_box() is the global placeholder, not the area"


def _params(grids):
    return [pytest.param(grid, id=grid.name) for grid in grids]


METADATA_PARAMS = _params(NAMED)

ORDER_BROKEN_PARAMS = [
    pytest.param(grid.name, id=grid.name, marks=pytest.mark.xfail(run=False, reason=grid.crash))
    for grid in NAMED
    if grid.order == CRASH
]

POINTS_PARAMS = _params([grid for grid in NAMED if grid.points])

POINTS_BROKEN_PARAMS = [
    pytest.param(grid, id=grid.name, marks=pytest.mark.xfail(run=False, reason=grid.broken_points))
    for grid in NAMED
    if grid.broken_points
]


# ---------------------------------------------------------------------------
# the tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("grid", METADATA_PARAMS)
def test_named_grid_metadata(grid):
    """spec, type, shape, size, bounding_box, order and the spec round trip."""
    named = Grid({"grid": grid.name})

    _same(named.spec, grid.spec)

    assert named.type == grid.grid_type, f"{grid.name}: type"
    assert tuple(named.shape) == grid.shape, f"{grid.name}: shape"
    assert named.size() == grid.size, f"{grid.name}: size"
    assert named.size() == math.prod(named.shape), f"{grid.name}: size is not the product of shape"

    if grid.bbox is not None:
        assert tuple(named.bounding_box()) == pytest.approx(grid.bbox), f"{grid.name}: bounding_box"
    else:
        assert grid.area_bbox, f"{grid.name}: no expected bounding_box, and not a regional grid"
        _assert_area(named.bounding_box(), grid.name)

    # The sentinel is never read: for these types it aborts the interpreter
    # (finding F1), which test_named_grid_order_is_not_implemented records.
    if grid.order != CRASH:
        assert named.order == grid.order, f"{grid.name}: order"

    assert Grid(named.spec) == named, f"{grid.name}: the canonical spec does not rebuild the grid"


@pytest.mark.parametrize("name", ORDER_BROKEN_PARAMS)
def test_named_grid_order_is_not_implemented(name):
    """BROKEN: grid.order aborts the process for every ORCA, ICON and FESOM grid.

    eckit throws a C++ ``eckit::NotImplemented`` which nothing catches, so the
    interpreter dies with SIGABRT (exit code 134) and no Python test run can
    survive it. Every other member of these grids works, which
    :func:`test_named_grid_metadata` checks. Verified again for each type on
    2026-11-18.
    """
    assert Grid({"grid": name}).order


@pytest.mark.parametrize("grid", POINTS_PARAMS)
def test_named_grid_points(grid):
    """to_latlons() gives one point per grid point, on the globe, inside the area."""
    named = Grid({"grid": grid.name})
    assert named.size() < MAX_POINTS, f"{grid.name}: above the resource gate of this file"

    latlons = named.to_latlons()
    assert len(latlons) == 2
    lat, lon = latlons
    assert len(lat) == named.size(), f"{grid.name}: {len(lat)} latitudes for {named.size()} points"
    assert len(lon) == named.size(), f"{grid.name}: {len(lon)} longitudes for {named.size()} points"

    assert -90.0 <= min(lat), f"{grid.name}: latitude {min(lat)}"
    assert max(lat) <= 90.0, f"{grid.name}: latitude {max(lat)}"
    assert -360.0 <= min(lon), f"{grid.name}: longitude {min(lon)}"
    assert max(lon) <= 360.0, f"{grid.name}: longitude {max(lon)}"

    if grid.area_bbox:
        # The area of the grid must hold its points. A catalogue grid whose
        # bounding box is the global placeholder is not checked: there the box is
        # the area of the grid, not the extent of its points (finding F5).
        north, west, south, east = named.bounding_box()
        assert west - AREA_TOL <= min(lon), f"{grid.name}: longitude {min(lon)} west of {west}"
        assert max(lon) <= east + AREA_TOL, f"{grid.name}: longitude {max(lon)} east of {east}"
        assert south - AREA_TOL <= min(lat), f"{grid.name}: latitude {min(lat)} south of {south}"
        assert max(lat) <= north + AREA_TOL, f"{grid.name}: latitude {max(lat)} north of {north}"


@pytest.mark.parametrize("grid", POINTS_BROKEN_PARAMS)
def test_named_grid_points_are_broken(grid):
    """BROKEN: HN10 builds, says it has 1200 pixels, and has no coordinates."""
    assert len(Grid({"grid": grid.name}).to_latlons()[0]) == grid.size


def test_broken_params_are_marked():
    """The BROKEN cases are xfail with run=False and a reason, nothing else.

    run=False is what keeps them out of the process: an xfail which is executed
    would abort the test run instead of reporting a failure, and the misleading
    error of F6 would be asserted as if it were correct.
    """
    assert ORDER_BROKEN_PARAMS and POINTS_BROKEN_PARAMS and NAME_ERROR_PARAMS

    for param in ORDER_BROKEN_PARAMS + POINTS_BROKEN_PARAMS + NAME_ERROR_PARAMS:
        assert len(param.marks) == 1, param.id
        mark = param.marks[0]
        assert mark.name == "xfail", param.id
        assert mark.kwargs.get("run") is False, param.id
        assert str(mark.kwargs.get("reason")).startswith("BROKEN: "), param.id


def test_classical_grid_sizes_are_between_octahedral_and_regular():
    """For the same n, N<n> has more points than O<n> and fewer than F<n>.

    The bounds of the tabulated linear grids: the octahedral grid is the sparsest
    and the regular Gaussian grid the densest way of placing points on n
    latitudes. No grid is built.
    """
    for n in (16, 32, 64, 1280):
        assert octahedral_points(n) < linear_points(n) < math.prod(regular_gg_shape(n)), f"N{n}"


@pytest.mark.parametrize("name", INVALID_NAMES)
def test_invalid_name_does_not_build(name):
    """An invalid name does not build, whatever its error message says.

    Only that the build raises is asserted. The message these names get is the
    finding F6, recorded as a non-running BROKEN case by
    :func:`test_unknown_name_error_is_clear_broken`: asserting today's wrong
    message as the expected one would assert a bug as correct behaviour.
    """
    with pytest.raises(RuntimeError):
        Grid({"grid": name})


@pytest.mark.parametrize("name", NAME_ERROR_PARAMS)
def test_unknown_name_error_is_clear_broken(name):
    """BROKEN (F6): rejecting an invalid name should name the problem.

    Today ``Grid({"grid": name})`` raises ``SpecError: [Grid: cannot build grid
    without 'type']`` for every one of these names, as if the gridspec had no
    type at all — re-verified per name in a subprocess on 2026-11-18; the dict
    form and the plain string form give the same message. The spec does have a
    resolved type; it is the *name* which is not a grid, and the expected
    behaviour is a clear "unknown grid name" error which says so. The body
    asserts what should happen; it is never executed (run=False), so the finding
    is recorded without the wrong message ever passing as correct. That the
    build raises at all is checked, message-free, by
    :func:`test_invalid_name_does_not_build`.
    """
    with pytest.raises(RuntimeError) as raised:
        Grid({"grid": name})
    assert "unknown grid" in str(raised.value).lower(), f"{name}: {raised.value}"


def test_string_form_is_the_same_grid():
    """The one small test of the plain string form: Grid(name) still works.

    Every other test in this file builds the gridspec dict. For a representative
    name of each family, the string and the dict form give the same grid.
    """
    for name in STRING_FORM_NAMES:
        assert Grid(name) == Grid({"grid": name}), f"{name}: string form differs"


def test_classical_grid_names_start_at_n16():
    """N4, N8, N12 and N21 are not in the table of classical reduced grids.

    ``N<n>`` is only a name for the values of n of the classical (linear) point
    count table; the others raise ``Bad value: reduced_classical_pl: unknown
    N=<n>``. This message is a clear rejection, unlike the finding F6 above, so
    it is asserted as-is. Worth recording, because an arbitrary ``pl`` list is
    canonicalised into an invented ``N<len(pl)/2>`` name which then cannot be
    built again (finding F2, see test_grid_types.py).
    """
    for n in (4, 8, 12, 21):
        with pytest.raises(RuntimeError) as raised:
            Grid({"grid": f"N{n}"})
        assert f"reduced_classical_pl: unknown N={n}" in str(raised.value), f"N{n}: {raised.value}"


# ---------------------------------------------------------------------------
# the catalogue allowlist
# ---------------------------------------------------------------------------


def _key(value):
    """The name of a YAML key: grid.yaml uses the key ``false``, a boolean to YAML."""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _entry(item):
    """Return the (name, fields) of one item of a ``grid_names`` list.

    The catalogues use two layouts. In grid.yaml and hydro.yaml an item is
    ``{NAME: {fields}}``. In the ORCA, ICON and FESOM files the YAML anchor sits
    on the entry name, so an item is one flat dict in which the name is the key
    with a null value next to the fields. An alias, as ICON-CH has, is a name
    with nothing else, because safe_load resolves the anchor to that null.
    """
    if not isinstance(item, dict):
        return _key(item), {}

    if len(item) == 1:
        name, fields = next(iter(item.items()))
        return _key(name), fields if isinstance(fields, dict) else {}

    names = [key for key, value in item.items() if value is None]
    if names:
        first = names[0]
        return _key(first), {key: value for key, value in item.items() if key is not first}

    return "?", dict(item)


def _catalogues():
    """{file name: {entry name: fields}} for the grid entries of every catalogue."""
    folder = catalogue.catalogue_dir()
    assert folder is not None, "the catalogue folder of the installed eckit was not found"

    catalogues = {}
    for path in sorted(folder.glob("*.yaml")):
        data = yaml.safe_load(path.read_text())
        # Only grid_names: area.yaml and projection.yaml list other things, and
        # the grid_uids lists of the grid catalogues hold uids, which are not
        # names.
        items = data.get("grid_names") if isinstance(data, dict) else None
        if isinstance(items, list):
            catalogues[path.name] = dict(_entry(item) for item in items)

    assert catalogues, f"no catalogue file with grid_names in {folder}"
    return catalogues


def test_every_catalogue_name_is_reviewed():
    """Every name of every catalogue grid is tested here or listed with a reason.

    The catalogue is read, no grid is built: a new name in a catalogue file fails
    until it is either added to the tables above or to :data:`NOT_TESTED`, so
    nothing goes untested in silence. A name which leaves the catalogue fails
    too, to keep the allowlist honest.
    """
    catalogues = _catalogues()
    listed = CATALOGUE_TESTED | set(NOT_TESTED)
    in_catalogue = {name for entries in catalogues.values() for name in entries}

    unreviewed = sorted(name for name in in_catalogue if name not in listed)
    assert not unreviewed, f"catalogue names which are neither tested nor in NOT_TESTED: {unreviewed}"

    stale = sorted(name for name in listed if name not in in_catalogue)
    assert not stale, f"names tested or listed here which are no longer in the catalogue: {stale}"


def test_expected_data_matches_the_catalogue():
    """The shapes, uids and bounding boxes the tables expect are the declared ones.

    This compares this file with the catalogue data, not with the code: it
    catches a table which was written for another version of the catalogues.
    """
    catalogues = _catalogues()
    entries = {name: fields for catalogued in catalogues.values() for name, fields in catalogued.items()}

    # Catalogue names which declare a shape, and are therefore checkable: every
    # ORCA, ICON and FESOM entry except the ICON-CH aliases, which declare
    # nothing and are checked through the entries they alias.
    checkable = {name for name, fields in entries.items() if "shape" in fields and name in CATALOGUE_TESTED}

    checked = set()
    for grid in NAMED:
        # grid.yaml and hydro.yaml declare no shape: those grids are checked
        # against the formulas above, not against their entry.
        if grid.name not in checkable:
            continue

        declared = entries[grid.name]
        assert tuple(declared["shape"]) == grid.shape, f"{grid.name}: shape in the catalogue"
        checked.add(grid.name)

        uid = grid.spec.get("uid")
        if uid is not None:
            assert declared.get("uid") == uid, f"{grid.name}: uid in the catalogue"

        box = declared.get("bounding_box")
        if box is not None:
            assert tuple(box) == grid.bbox, f"{grid.name}: bounding_box in the catalogue"

    assert checked == checkable, f"catalogue shapes not checked: {sorted(checkable - checked)}"
    # A floor, so that the flat layout of the grid catalogues cannot quietly
    # parse into entries without fields: every ORCA and FESOM name declares one.
    assert {grid.name for grid in ORCA_GRIDS + FESOM_GRIDS} <= checkable, "catalogue shapes were not parsed"
