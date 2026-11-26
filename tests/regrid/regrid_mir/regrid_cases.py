# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""The grid table of the regrid matrix tests: sources, targets and per-interpolation selections.

This module is data, not a test: it is imported by ``test_regrid_numpy_grids.py``, by
``generate_regrid_data.py`` and by nothing in ``tests/grid/``. It imports only the standard
library at module level; ``earthkit.geo`` is imported inside :func:`check_table` alone, so the
table can be read (and the matrix walked) without building a single grid.

The table holds EVERY grid used by the grid tests -- ``tests/grid/test_grid_types.py``
(``ALL_CASES`` and the specs of its BROKEN lists which do build), ``test_grid_named.py``
(all of ``NAMED``, not a sample), ``test_grid_sh.py`` and ``test_grid.py`` -- copied from
those files as gridspec dicts, so that a change in a grid test never silently changes a
regrid test. Grids which cannot be regrid sources stay in the table with ``covered=False``
and the reason, see :data:`NOT_COVERED`. Equivalent spellings of one grid (same points, same
order, same uid) are listed once; the merged spellings are named in the comment of the entry
they were merged into.

Contents
--------
* :class:`Source` / :data:`SOURCES`: every grid of the grid tests as a regrid source.
* :class:`Target` / :data:`TARGETS`: the four grids the sources are regridded to.
* :data:`SELECTION`: canonical interpolation name -> ``(source ids, target ids)``.
* :data:`NOT_COVERED`: source id -> why that source is not in the matrix.
* :data:`BROKEN` / :data:`CLIPPED`: per-combination findings, filled in from the probe.
* :func:`matrix`: the ``(source, target)`` pairs of one interpolation, in a stable order.
* :func:`check_table`: builds the grids and reports where the table is wrong.

Adding a grid
-------------
Copy its gridspec from the grid test that uses it, give it a short id of letters, digits,
``-``, ``_`` and ``.`` only (ids become file names: ``input/<source>.npz``), add one
``_source(id, spec, shape, ...)`` line to the type section it belongs to, and run
``python -c "import regrid_cases as c; print(c.check_table())"`` from this folder: it must
return ``[]``. Coverage needs no attention: a source of :data:`MAX_POINTS` points or more is
flagged ``too large`` by ``_source`` itself, so raising ``MAX_POINTS`` is the only change
needed to put a large grid into the matrix. A grid which cannot be a source at all (no
points, no coordinates) gets an explicit ``reason=``. Add the id to a selection list to have
the cheaper interpolations cover it, and list one spelling only: two specs which give the same
grid, same order and same uid, are one entry. Ids of grids which differ only in their scan
order spell ``+`` and ``-`` as ``p`` and ``m``: ``ll-order-ipp`` is the ``i+j+`` grid.

Adding a target
---------------
Append ``Target(id, spec, shape)`` to :data:`TARGETS`. ``nearest-neighbour`` picks new targets
up automatically (its selection is every covered source times every target); the selections
of the other interpolations name their targets and need to be edited.

Adding an interpolation
-----------------------
Register it in ``INTERPOLATIONS`` in ``regrid_common.py`` under its canonical name, then add
one key to :data:`SELECTION` with the source ids and target ids it should cover -- every
covered source and every target for an interpolation which must run on all of them, a hand-
picked subset otherwise. Nothing else changes: the tests are parametrised from
:func:`matrix`, and :data:`BROKEN` and :data:`CLIPPED` are keyed by that canonical name.

Never read ``Grid.order`` and never call ``to_latlons()`` in this module or from it: for
``sh``, ``unstructured_ll``, ``ORCA``, ``ICON`` and ``FESOM`` grids reading the order aborts
the interpreter (an uncaught C++ ``eckit::NotImplemented``, SIGABRT -- finding F1 of
``.ecmwf-multi-agents/notes/grid-probe.md``), and above :data:`MAX_POINTS` points the
coordinates are gigabytes. :func:`check_table` only builds grids and reads shape and size.
"""

import math
from collections.abc import Iterator
from dataclasses import dataclass

# to_latlons() and any other coordinate-materialising call is only possible below this many
# points; it is the same resource gate the grid tests use (SIZE_GATE / MAX_POINTS there).
MAX_POINTS = 1_000_000

# The reasons a source is not in the matrix, apart from being too large.
NO_POINTS = "sh has no grid points"
NO_COORDINATES = "to_latlons raises: Nside must be a power of 2"
NO_LATLONS = "no lat/lon: Cartesian grid"

GLOBAL = "global"
REGIONAL = "regional"


@dataclass(frozen=True)
class Source:
    """One grid of the grid tests, as a regrid source.

    Parameters
    ----------
    id: str
        Short filesystem-safe name, unique among the sources and used in file names.
    spec: dict
        The gridspec, as the grid test builds it.
    shape: tuple of int
        The expected ``Grid(spec).shape``.
    size: int
        The expected ``Grid(spec).size()``: the product of :attr:`shape`.
    region: str
        ``"global"`` or ``"regional"``, whether the grid declares an area.
    covered: bool
        False when the grid cannot be a source of the matrix; :attr:`reason` says why.
    reason: str
        Why the source is not covered, empty when it is covered.
    """

    id: str
    spec: dict
    shape: tuple
    size: int
    region: str
    covered: bool = True
    reason: str = ""


@dataclass(frozen=True)
class Target:
    """One grid the sources are regridded to.

    Parameters
    ----------
    id: str
        Short filesystem-safe name, unique among the targets and used in file names.
    spec: dict
        The gridspec of the target grid.
    shape: tuple of int
        The expected ``Grid(spec).shape`` of the regridded result.
    """

    id: str
    spec: dict
    shape: tuple


def _source(id, spec, shape, region=GLOBAL, reason=""):
    """One row of :data:`SOURCES`, with the size derived from the shape and the coverage from it.

    A source of :data:`MAX_POINTS` points or more is flagged ``too large`` here, so that the
    constant is the only knob; a grid which cannot be a source for another reason passes an
    explicit ``reason``.
    """
    shape = tuple(shape)
    size = math.prod(shape)
    if not reason and size >= MAX_POINTS:
        reason = f"too large: {size} points"
    return Source(id=id, spec=spec, shape=shape, size=size, region=region, covered=not reason, reason=reason)


# ---------------------------------------------------------------------------
# the sources: every grid of tests/grid, copied from the test that uses it
# ---------------------------------------------------------------------------

# test_grid_types.py REGULAR_LL_CASES, plus the <dlon>/<dlat> pattern names of
# test_grid_named.py. The pattern names "1/1", "0.5/0.5", "1.5/0.5" name the grids of
# ll-1x1, ll-0.5x0.5 and ll-1.5x0.5 (same uid), "2/2" and "0.25/0.25" are grids of their own.
REGULAR_LL_SOURCES = (
    # ll-1x1 is also test_grid.py:19 (dict(grid=[1, 1])) and the pattern name "1/1"
    _source("ll-1x1", {"grid": [1, 1]}, (181, 360)),
    _source("ll-0.5x0.5", {"grid": [0.5, 0.5]}, (361, 720)),
    _source("ll-1.5x0.5", {"grid": [1.5, 0.5]}, (361, 240)),
    _source("ll-2x2", {"grid": [2, 2]}, (91, 180)),  # pattern name "2/2"
    _source("ll-0.25x0.25", {"grid": [0.25, 0.25]}, (721, 1440)),  # pattern name "0.25/0.25"
    _source("ll-area-0.1", {"grid": [0.1, 0.1], "area": [72, -25, 25, 45]}, (471, 701), region=REGIONAL),
    # the string area spelling: the same grid as the target ll-1x1-area, kept as a source
    _source("ll-area-string", {"grid": [1, 1], "area": "72/-25/25/45"}, (48, 71), region=REGIONAL),
    # the ac-T-n90 spec {"type": "arakawa_c", "n": 90} is this grid (same uid), listed once
    _source("ll-reference", {"grid": [1, 1], "reference": [0.5, 0.5]}, (180, 360)),
    # scan orders: "+" and "-" of the order name are "p" and "m" in the ids (file system safe)
    _source("ll-order-ipp", {"grid": [1, 1], "order": "i+j+"}, (181, 360)),
    _source("ll-order-mmm", {"grid": [1, 1], "order": "i-j-"}, (181, 360)),
    _source("ll-order-mpp", {"grid": [1, 1], "order": "i-j+"}, (181, 360)),
    _source("ll-rotation", {"grid": [1, 1], "rotation": [40, -30]}, (181, 360)),
    _source("ll-rotated_ll-alias", {"type": "rotated_ll", "grid": [1, 1], "rotation": [30, 90]}, (181, 360)),
    _source(
        "ll-rotated-area",
        {
            "grid": [0.02, 0.02],
            "area": [3.36, -6.82, -4.42, 4.8],
            "order": "i+j+",
            "projection": {"type": "rotation", "south_pole": [10, -43]},
        },
        (390, 582),
        region=REGIONAL,
    ),
    _source("ll-figure-sun", {"grid": [2, 2], "figure": "sun"}, (91, 180)),
    # XY_NO_TYPE_CASES of test_grid_types.py (finding F10): these regular_xy specs silently
    # build a regular_ll. {"grid": [2, 2], "bounding_box_xy": ...} is then the grid of ll-2x2
    # (same uid, checked), so it is not listed again; the laea one is a grid of its own.
    _source(
        "xy-no-type-laea",
        {
            "grid": [10, 10],
            "bounding_box_xy": [-100, -100, 100, 100],
            "false_lonlat": [0, 52],
            "projection": {"type": "laea", "lat_0": 52, "lon_0": 0},
        },
        (19, 36),
    ),
)

# test_grid_types.py ARAKAWA_C_CASES. ac-T-n90 is the grid of ll-reference above.
ARAKAWA_C_SOURCES = (
    _source("ac-U-n90", {"type": "arakawa_c", "n": 90, "arrangement": "U"}, (180, 360)),
    _source("ac-V-n90", {"type": "arakawa_c", "n": 90, "arrangement": "V"}, (181, 360)),
    _source("ac-grid_factor-n90", {"type": "arakawa_c", "n": 90, "grid_factor": [2, 1.5]}, (120, 180)),
    _source("ac-area-n600", {"type": "arakawa_c", "n": 600, "area": [72, -25, 25, 45]}, (313, 467), region=REGIONAL),
    _source("ac-order-ipp-n90", {"type": "arakawa_c", "n": 90, "order": "i+j+"}, (180, 360)),
    _source("ac-rotation-n90", {"type": "arakawa_c", "n": 90, "rotation": [40, -30]}, (180, 360)),
)

# test_grid_types.py ARAKAWA_C_UM_CASES, plus the grids of test_grid.py. All the equivalent
# spellings of the N96 T grid there (kwargs, yaml strings, explicit grid and reference or
# area) are ac-um-T-N96; only the 8-decimal grid_factor of test_grid.py:67 is another grid.
ARAKAWA_C_UM_SOURCES = (
    _source("ac-um-T-N96", {"type": "arakawa_c_um", "N": 96}, (144, 192)),
    _source("ac-um-U-N96", {"type": "arakawa_c_um", "N": 96, "arrangement": "U"}, (144, 192)),
    _source("ac-um-V-N96", {"type": "arakawa_c_um", "N": 96, "arrangement": "V"}, (145, 192)),
    _source("ac-um-T-N48", {"type": "arakawa_c_um", "N": 48}, (72, 96)),
    _source(
        "ac-n96-gf-2x1.33333333",
        {"type": "arakawa_c", "N": 96, "grid_factor": [2, 1.33333333], "order": "i+j+"},
        (144, 192),
    ),
)

# test_grid_types.py REGULAR_GG_CASES and the F<N> pattern names: {"type": "regular_gg",
# "N": 16} and {"type": "regular_gg", "N": 32} are the grids named F16 and F32 (same uid).
REGULAR_GG_SOURCES = (
    _source("F16", {"grid": "F16"}, (32, 64)),
    _source("F32", {"grid": "F32"}, (64, 128)),
    _source("F48", {"grid": "F48"}, (96, 192)),
    _source("F1600", {"grid": "F1600"}, (3200, 6400)),
    _source("gg-rotated-N16", {"type": "rotated_gg", "N": 16, "rotation": [40, -30]}, (32, 64)),
    _source("gg-area-N32", {"type": "regular_gg", "N": 32, "area": [80, -20, 20, 40]}, (22, 22), region=REGIONAL),
    _source("gg-order-ipp-N16", {"type": "regular_gg", "N": 16, "order": "i+j+"}, (32, 64)),
)

# test_grid_types.py REDUCED_GG_CASES and the O<N>/N<N> pattern names. {"type": "reduced_gg",
# "N": 32} is O32 (the N shortcut means octahedral, finding F12) and
# {"type": "reduced_rotated_gg", "N": 32, "rotation": [40, -30]} is O32-rotated (same uids).
# The order variant of the pl grids, {"pl": [8, 8, 8, 8], "order": "i+j+"}, is rd-pl-4lat:
# order is silently ignored on reduced_gg (finding F3), the uid is the same.
REDUCED_GG_SOURCES = (
    _source("O8", {"grid": "O8"}, (544,)),
    _source("O16", {"grid": "O16"}, (1600,)),
    _source("O32", {"grid": "O32"}, (5248,)),
    _source("O64", {"grid": "O64"}, (18688,)),
    _source("O128", {"grid": "O128"}, (70144,)),
    _source("O1280", {"grid": "O1280"}, (6599680,)),
    _source("N16", {"grid": "N16"}, (1688,)),
    _source("N32", {"grid": "N32"}, (6114,)),
    _source("N64", {"grid": "N64"}, (23112,)),
    _source("N1280", {"grid": "N1280"}, (8505906,)),
    _source("O32-rotated", {"grid": "O32", "rotation": [40, -30]}, (5248,)),
    _source("O32-area", {"grid": "O32", "area": [80, -20, 20, 40]}, (274,), region=REGIONAL),
    _source("rd-pl-8lat", {"type": "reduced_gg", "pl": [4, 8, 12, 16, 16, 12, 8, 4]}, (80,)),
    _source("rd-pl-4lat", {"type": "reduced_gg", "pl": [8, 8, 8, 8]}, (32,)),
)

# test_grid_types.py HEALPIX_CASES and the H<N>/HR<N>/HN<N> pattern names: {"type": "healpix",
# "Nside": 4} is H4, {"type": "HEALPix", "Nside": 12} is H12, HR16 (hp-HR16 and the pattern
# name) is the ring order of H16, and HN16 (hp-HN16 and the pattern name) is its nested order,
# listed as H16-nested. H8 is the same grid as the target H8-ring.
HEALPIX_SOURCES = (
    _source("H4", {"grid": "H4"}, (192,)),
    _source("H8", {"grid": "H8"}, (768,)),
    _source("H8-nested", {"grid": "H8", "order": "nested"}, (768,)),
    _source("H12", {"grid": "H12"}, (1728,)),
    _source("H16", {"grid": "H16"}, (3072,)),
    _source("H16-nested", {"grid": "HN16"}, (3072,)),
    _source("H31", {"grid": "H31"}, (11532,)),
    _source("H1024", {"grid": "H1024"}, (12582912,)),
    # nested with a Nside which is not a power of 2: metadata works, points do not (finding F4)
    _source("HN10", {"grid": "HN10"}, (1200,), reason=NO_COORDINATES),
)

# test_grid_types.py REGULAR_XY_CASES. {"x": [0, 10, 1], "y": [0, 10, 1], "order": "i+j+"} is
# xy-plain (same uid) and the plain {"x": [0, 100, 10], "y": [0, 100, 10]} of the false_xy
# BROKEN test is xy-false_xy: a dropped false_xy does not change the grid (a T3 finding).
# Only grids WITHOUT a projection are out: their points are Cartesian, not (lon, lat).
REGULAR_XY_SOURCES = (
    _source("xy-plain", {"x": [0, 10, 1], "y": [0, 10, 1]}, (11, 11), reason=NO_LATLONS),
    _source("xy-descending-x", {"x": [10, 0, -1], "y": [0, 10, 1]}, (11, 11), reason=NO_LATLONS),
    _source(
        "xy-false_xy", {"x": [0, 100, 10], "y": [0, 100, 10], "false_xy": [1000, 2000]}, (11, 11), reason=NO_LATLONS
    ),
    _source(
        "xy-laea-grid-bbox-reverse",
        {
            "type": "regular_xy",
            "grid": [5000, 5000],
            "bounding_box_xy": [-500000, -400000, 500000, 400000],
            "projection": {"type": "laea", "lat_0": 52, "lon_0": 10, "reverse": True},
        },
        (161, 201),
        region=REGIONAL,
    ),
    # the coarse named-projection grid of the types test; the catalogue grid.yaml name
    # swisslv03 below is the same area at 1 km: a different, finer grid
    _source(
        "xy-swisslv03",
        {"projection": "swisslv03", "x": [260000, 960000, 50000], "y": [-150000, 450000, 50000]},
        (13, 15),
        region=REGIONAL,
    ),
)

# test_grid_types.py UNSTRUCTURED_CASES. The {"type": "unstructured", ...} alias and all three
# specs of the order-BROKEN list are these two grids (same uids, checked); reading .order on
# this type aborts the interpreter (finding F1), which is why the grids are here but their
# order never is.
UNSTRUCTURED_SOURCES = (
    _source(
        "uns-4pt",
        {"latitudes": [0.0, 10.5, -20.25, 45.0], "longitudes": [0.0, 90.0, 180.0, 270.5]},
        (4,),
    ),
    _source(
        "uns-named",
        {
            "type": "unstructured_ll",
            "latitudes": [0.0, 1.0, 2.0, 3.0, 4.0],
            "longitudes": [10.0, 11.0, 12.0, 13.0, 14.0],
            "name": "EARTHKIT-GEO-T3-UNS",
        },
        (5,),
    ),
)

# test_grid_sh.py TRUNCATIONS: an sh "grid" is a spectral truncation with coefficients, not
# coordinates, so none of these is a source (finding F8). Shape and size are the flat
# coefficient count (T + 1)(T + 2) of docs/concepts/grid/types/spectral.rst.
SH_SOURCES = tuple(
    _source(f"T{t}", {"grid": f"T{t}"}, ((t + 1) * (t + 2),), reason=NO_POINTS)
    for t in (10, 21, 42, 63, 85, 159, 319, 639)
)

# test_grid_named.py orca_grids(): the 7 base grids, all 5 arrangements, shapes from the
# catalogue. The *_W grids carry the coordinates of their *_T grid (finding F9) but have
# their own uid, so they are listed as the separate grids of the grid tests.
ORCA_SOURCES = tuple(
    _source(f"{base}_{arrangement}", {"grid": f"{base}_{arrangement}"}, shape)
    for base, shape in (
        ("ORCA2", (149, 182)),
        ("ORCA1", (292, 362)),
        ("eORCA1", (332, 362)),
        ("ORCA025", (1021, 1442)),
        ("eORCA025", (1207, 1442)),
        ("ORCA12", (3059, 4322)),
        ("eORCA12", (3606, 4322)),
    )
    for arrangement in "FTUVW"
)

# test_grid_named.py icon_grids(): name and point count, as the catalogue declares them. The
# aliases ICON-CH1-V1 and ICON-CH2-V1 are ICON-CH1 and ICON-CH2 (same uid), listed once. The
# two 20480-point grids 0010-r02b04-g and 0054-r02b04-g are different grids (different uid).
ICON_ENTRIES = (
    ("ICON-CH1", 1147980, REGIONAL),
    ("ICON-CH2", 283876, REGIONAL),
    ("icon-grid-0016-r02b09-o", 14886338, GLOBAL),
    ("icon-grid-0020-r02b05-o", 59261, GLOBAL),
    ("icon-grid-0024-r02b07-o", 933398, GLOBAL),
    ("icon-grid-0036-r02b04-o", 15105, GLOBAL),
    ("icon-grid-0038-r02b11-o", 237316813, GLOBAL),
    ("icon-grid-0040-r02b10-o", 59359799, GLOBAL),
    ("icon-grid-0046-r02b06-o", 235403, GLOBAL),
    ("icon-grid-0001-r02b05-r", 81920, GLOBAL),
    ("icon-grid-0002-r02b06-g", 327680, GLOBAL),
    ("icon-grid-0004-r02b07-g", 1310720, GLOBAL),
    ("icon-grid-0009-r02b03-r", 5120, GLOBAL),
    ("icon-grid-0010-r02b04-g", 20480, GLOBAL),
    ("icon-grid-0020-r03b06-g", 737280, GLOBAL),
    ("icon-grid-0022-r03b07-g", 2949120, GLOBAL),
    ("icon-grid-0027-r03b08-n02", 659156, GLOBAL),
    ("icon-grid-0031-r02b06-n02", 19744, GLOBAL),
    ("icon-grid-0042-r02b10-l", 253520, GLOBAL),
    ("icon-grid-0043-r19b06-lr", 134541, GLOBAL),
    ("icon-grid-0045-r19b08-ln02", 716980, GLOBAL),
    ("icon-grid-0054-r02b04-g", 20480, GLOBAL),
    ("icon-grid-0055-r02b05-n", 2656, GLOBAL),
)

ICON_SOURCES = tuple(_source(name, {"grid": name}, (points,), region=region) for name, points, region in ICON_ENTRIES)

# test_grid_named.py fesom_grids(): the cell (C) and node (N) arrangements of 4 meshes, with
# the points each arrangement has in the catalogue.
FESOM_SOURCES = tuple(
    _source(name, {"grid": name}, (points,))
    for name, points in (
        ("pi_C", 5839),
        ("pi_N", 3140),
        ("CORE2_C", 244659),
        ("CORE2_N", 126858),
        ("DART_C", 6262485),
        ("DART_N", 3160340),
        ("NG5_C", 14741520),
        ("NG5_N", 7402886),
    )
)

# test_grid_named.py grid_yaml_grids(): the four regional projected grids of grid.yaml, as
# catalogue names (the grid tests build Grid({"grid": name})).
GRID_YAML_SOURCES = (
    _source("LAEA-EFAS-5km", {"grid": "LAEA-EFAS-5km"}, (950, 1000), region=REGIONAL),
    _source("SMUFF-OPERA-2km", {"grid": "SMUFF-OPERA-2km"}, (2200, 1900), region=REGIONAL),
    _source("swisslv95", {"grid": "swisslv95"}, (294, 429), region=REGIONAL),
    _source("swisslv03", {"grid": "swisslv03"}, (640, 710), region=REGIONAL),
)

# test_grid_named.py hydro_grids(): the seven hydrological grids of hydro.yaml. All are large;
# EFAS-5 and GLOFAS-4 cover a region, the CAMA grids are global.
HYDRO_SOURCES = (
    _source("CAMA-15min", {"grid": "CAMA-15min"}, (720, 1440)),
    _source("CAMA-06min", {"grid": "CAMA-06min"}, (1800, 3600)),
    _source("CAMA-05min", {"grid": "CAMA-05min"}, (2160, 4320)),
    _source("CAMA-03min", {"grid": "CAMA-03min"}, (3600, 7200)),
    _source("CAMA-01min", {"grid": "CAMA-01min"}, (10800, 21600)),
    _source("EFAS-5", {"grid": "EFAS-5"}, (2970, 4530), region=REGIONAL),
    _source("GLOFAS-4", {"grid": "GLOFAS-4"}, (3000, 7200)),
)

SOURCES = (
    REGULAR_LL_SOURCES
    + ARAKAWA_C_SOURCES
    + ARAKAWA_C_UM_SOURCES
    + REGULAR_GG_SOURCES
    + REDUCED_GG_SOURCES
    + HEALPIX_SOURCES
    + REGULAR_XY_SOURCES
    + UNSTRUCTURED_SOURCES
    + SH_SOURCES
    + ORCA_SOURCES
    + ICON_SOURCES
    + FESOM_SOURCES
    + GRID_YAML_SOURCES
    + HYDRO_SOURCES
)

# ---------------------------------------------------------------------------
# the targets: four grids every source is regridded to
# ---------------------------------------------------------------------------

TARGETS = (
    Target("ll-10x10", {"grid": [10, 10]}, (19, 36)),
    Target("O32", {"grid": "O32"}, (5248,)),
    Target("H8-ring", {"grid": "H8", "order": "ring"}, (768,)),
    Target("ll-1x1-area", {"grid": [1, 1], "area": [72, -25, 25, 45]}, (48, 71)),
)

# ---------------------------------------------------------------------------
# what runs against what
# ---------------------------------------------------------------------------

# The small grids the interpolations which are not nearest-neighbour are checked on: a regular
# global and a regular regional ll, a regular Gaussian, an octahedral and a classical reduced
# Gaussian, HEALPix in both orders, an ORCA grid and a FESOM grid.
SMALL_SOURCES = (
    "ll-1x1",
    "ll-area-string",
    "F32",
    "O32",
    "N32",
    "H8",
    "H8-nested",
    "ORCA2_T",
    "pi_C",
)

# Two targets for the cheaper selections: a global regular grid and a reduced Gaussian one.
SMALL_TARGETS = ("ll-10x10", "O32")

# Canonical interpolation name -> (source ids, target ids). Nearest-neighbour covers every
# covered source against every target; the others a hand-picked subset. The names must match
# the keys of INTERPOLATIONS in regrid_common.py.
SELECTION = {
    "nearest-neighbour": (tuple(s.id for s in SOURCES if s.covered), tuple(t.id for t in TARGETS)),
    "linear": (SMALL_SOURCES, SMALL_TARGETS),
    "grid-box-average": (SMALL_SOURCES, SMALL_TARGETS),
}

# Source id -> why it is not a source of the matrix, derived from the table.
NOT_COVERED = {s.id: s.reason for s in SOURCES if not s.covered}

# ---------------------------------------------------------------------------
# the findings per combination, filled in from the probe of the matrix
# ---------------------------------------------------------------------------

# (interpolation, source id, target id) -> reason, for a combination whose result cannot be
# produced at all: MIR raises, aborts the process or hangs. The tests turn these into
# pytest.param(..., marks=pytest.mark.xfail(run=False, reason=...)), so they are listed and
# never silently skipped. Every reason starts with "BROKEN: " and says what happens.
# The 7 reasons are the verbatim MIR errors of the probe (_dev/probe/broken.json).
BROKEN = {
    # Finding F1 of .ecmwf-multi-agents/notes/regrid-matrix-probe.md: the three Swiss grids fail
    # only against the global ll-10x10 target, whose 10 degree lattice misses their source box
    # entirely; they succeed against O32 (clipped to 1-4 points), H8-ring (NaN-filled, see the
    # note under CLIPPED) and ll-1x1-area (clipped).
    ("nearest-neighbour", "swisslv03", "ll-10x10"): (
        "BROKEN: RuntimeError: Serious bug: KNearestNeighbours: failed to interpolate"
    ),
    ("nearest-neighbour", "swisslv95", "ll-10x10"): (
        "BROKEN: RuntimeError: Serious bug: KNearestNeighbours: failed to interpolate"
    ),
    ("nearest-neighbour", "xy-swisslv03", "ll-10x10"): (
        "BROKEN: RuntimeError: Serious bug: KNearestNeighbours: failed to interpolate"
    ),
    # Finding F2: grid-box-average needs grid boxes the unstructured meshes do not have. ORCA2_T
    # and pi_C are the only unstructured sources of its selection and both fail on both global
    # targets; the same sources are fine with nearest-neighbour and linear. Expect this for every
    # ORCA / FESOM / ICON source if the grid-box-average selection ever grows.
    ("grid-box-average", "ORCA2_T", "O32"): (
        'BROKEN: RuntimeError: Serious bug: Representation::gridBoxes() not implemented for {"grid":"ORCA2_T"}'
    ),
    ("grid-box-average", "ORCA2_T", "ll-10x10"): (
        'BROKEN: RuntimeError: Serious bug: Representation::gridBoxes() not implemented for {"grid":"ORCA2_T"}'
    ),
    ("grid-box-average", "pi_C", "O32"): (
        'BROKEN: RuntimeError: Serious bug: Representation::gridBoxes() not implemented for {"grid":"pi_C"}'
    ),
    ("grid-box-average", "pi_C", "ll-10x10"): (
        'BROKEN: RuntimeError: Serious bug: Representation::gridBoxes() not implemented for {"grid":"pi_C"}'
    ),
}

# (interpolation, source id, target id) -> the gridspec the regridded result describes, for a
# combination where MIR clips the output to the area of the source instead of returning the
# target grid (an ICON-CH source onto a global target gives the small European grid). The
# tests assert that gridspec instead of the target, and it is written by hand from the source
# area, never copied blindly from the run. A combination which is neither the target nor a
# sensible clip is BROKEN, not CLIPPED.
#
# All 35 entries come from the probe (_dev/probe/clipped.json and _dev/probe/atypical.json) and
# all 34 clips are produced by the 12 covered regional sources; H8-ring is never clipped at all.
# The last entry is not a clip in the "sub-area" sense: it keeps the target shape but carries the
# source's rotation into the returned spec (finding F3); it lives here because its values are
# valid and the odd spec should stay visible as a regression test, not as an xfail.
CLIPPED = {
    # 32 hand-verified clips (probe section 4): for each of these, keeping the target points whose
    # (lat, lon) lie inside the source Grid.bounding_box() reproduces the returned gridspec
    # exactly, in point count and bounding box (_dev/probe/check_clipped.py). The three F4
    # entries below count among the 32 -- their point counts are verified, only the gridspec
    # spelling is degenerate. The ll-1x1-area clips are always the intersection of the target
    # area [72, -25, 25, 45] with the source box, so a Swiss source gives 6x10 of the 48x71 grid.
    ("grid-box-average", "ll-area-string", "O32"): {
        "grid": "O32",
        "area": [71.1577520115873, -25, 26.510769325211, 45],
    },
    ("grid-box-average", "ll-area-string", "ll-10x10"): {"grid": [10, 10], "area": [70, -20, 30, 40]},
    ("linear", "ll-area-string", "O32"): {"grid": "O32", "area": [71.1577520115873, -25, 26.510769325211, 45]},
    ("linear", "ll-area-string", "ll-10x10"): {"grid": [10, 10], "area": [70, -20, 30, 40]},
    ("nearest-neighbour", "ICON-CH2", "O32"): {
        "grid": "O32",
        "area": [48.8352409662506, 0, 43.2541946653509, 17.1428571428571],
    },
    ("nearest-neighbour", "ICON-CH2", "ll-1x1-area"): {"grid": [1, 1], "area": [50, 0, 42, 17]},
    ("nearest-neighbour", "LAEA-EFAS-5km", "O32"): {
        "grid": "O32",
        "area": [71.1577520115873, -35, 26.510769325211, 73.6363636363636],
    },
    ("nearest-neighbour", "LAEA-EFAS-5km", "ll-10x10"): {"grid": [10, 10], "area": [70, -30, 30, 70]},
    ("nearest-neighbour", "O32-area", "ll-10x10"): {"grid": [10, 10], "area": [70, -20, 30, 40]},
    ("nearest-neighbour", "O32-area", "ll-1x1-area"): {"grid": [1, 1], "area": [72, -20, 25, 40]},
    ("nearest-neighbour", "ac-area-n600", "O32"): {
        "grid": "O32",
        "area": [71.1577520115873, -24.5454545454545, 26.510769325211, 43.3333333333333],
    },
    ("nearest-neighbour", "ac-area-n600", "ll-10x10"): {"grid": [10, 10], "area": [70, -20, 30, 40]},
    ("nearest-neighbour", "ac-area-n600", "ll-1x1-area"): {"grid": [1, 1], "area": [71, -24, 26, 44]},
    ("nearest-neighbour", "gg-area-N32", "ll-10x10"): {"grid": [10, 10], "area": [70, -10, 30, 30]},
    ("nearest-neighbour", "gg-area-N32", "ll-1x1-area"): {"grid": [1, 1], "area": [72, -19, 25, 39]},
    ("nearest-neighbour", "ll-area-0.1", "O32"): {"grid": "O32", "area": [71.1577520115873, -25, 26.510769325211, 45]},
    ("nearest-neighbour", "ll-area-0.1", "ll-10x10"): {"grid": [10, 10], "area": [70, -20, 30, 40]},
    ("nearest-neighbour", "ll-area-string", "O32"): {
        "grid": "O32",
        "area": [71.1577520115873, -25, 26.510769325211, 45],
    },
    ("nearest-neighbour", "ll-area-string", "ll-10x10"): {"grid": [10, 10], "area": [70, -20, 30, 40]},
    ("nearest-neighbour", "ll-rotated-area", "O32"): {
        "grid": "O32",
        "area": [48.8352409662506, 0, 43.2541946653509, 17.1428571428571],
    },
    ("nearest-neighbour", "ll-rotated-area", "ll-1x1-area"): {"grid": [1, 1], "area": [50, 0, 43, 17]},
    ("nearest-neighbour", "swisslv03", "O32"): {
        "grid": "O32",
        "area": [48.8352409662506, 4.5, 46.0447266311017, 9.47368421052632],
    },
    ("nearest-neighbour", "swisslv03", "ll-1x1-area"): {"grid": [1, 1], "area": [49, 3, 44, 12]},
    ("nearest-neighbour", "swisslv95", "O32"): {"grid": "O32", "area": [46.0447266311017, 9, 46.0447266311017, 9]},
    ("nearest-neighbour", "swisslv95", "ll-1x1-area"): {"grid": [1, 1], "area": [48, 6, 46, 11]},
    ("nearest-neighbour", "xy-laea-grid-bbox-reverse", "O32"): {
        "grid": "O32",
        "area": [54.4161995260862, 4.73684210526316, 48.8352409662506, 15.8823529411765],
    },
    ("nearest-neighbour", "xy-laea-grid-bbox-reverse", "ll-1x1-area"): {"grid": [1, 1], "area": [55, 3, 49, 17]},
    ("nearest-neighbour", "xy-swisslv03", "O32"): {
        "grid": "O32",
        "area": [48.8352409662506, 4.5, 46.0447266311017, 9.47368421052632],
    },
    ("nearest-neighbour", "xy-swisslv03", "ll-1x1-area"): {"grid": [1, 1], "area": [49, 3, 44, 12]},
    # Finding F4 (3 of the 32 above): a clip too small for the target resolution produces a
    # degenerate gridspec -- a spacing of 0 is not a legal regular-ll resolution (probable MIR
    # bug). The point counts are exactly the target points inside the source box (1x2 = 2 points
    # for ICON-CH2 and ll-rotated-area, 1x1 = 1 point for xy-laea-grid-bbox-reverse) and the
    # values are valid, so the odd spelling stays visible as a regression test.
    ("nearest-neighbour", "ICON-CH2", "ll-10x10"): {"grid": [10, 0], "area": [50, 0, 50, 10]},
    ("nearest-neighbour", "ll-rotated-area", "ll-10x10"): {"grid": [10, 0], "area": [50, 0, 50, 10]},
    ("nearest-neighbour", "xy-laea-grid-bbox-reverse", "ll-10x10"): {"grid": [0, 0], "area": [50, 10, 50, 10]},
    # Finding F6 (2 of the 34, the exceptions of the hand check): for a reduced-Gaussian target
    # the clip is not exactly target ∩ source box -- O32-area and gg-area-N32 onto O32 return
    # 274 and 267 points while only 250 and 243 O32 points lie inside the source boxes; the
    # returned area equals the source box and the extra points sit on the rows the cut keeps
    # whole. Self-consistent and finite, so still a CLIPPED entry.
    ("nearest-neighbour", "O32-area", "O32"): {"grid": "O32", "area": [79.5256065726594, -20, 20.9295742544895, 40]},
    ("nearest-neighbour", "gg-area-N32", "O32"): {
        "grid": "O32",
        "area": [79.5256065726594, -19.5652173913043, 20.9295742544895, 39.375],
    },
    # Finding F3 (the 3 atypical pairs of _dev/probe/atypical.json, not clips of an area): the
    # result has the target shape but the returned gridspec carries the source's rotation into
    # the unrotated target, so its uid differs from the target's (probable MIR bug). Every value
    # check passes, so the pair stays visible here instead of as an xfail. It is the only one
    # left: since the backend passes a target's area inside the grid spec (2026-11-26), the two
    # rotated ll-1x1-area pairs of the probe -- ll-rotation and ll-rotated_ll-alias -- return
    # the target spec itself and are ordinary matrix pairs with no entry here; their references
    # were regenerated for the new (correct) result.
    ("nearest-neighbour", "O32-rotated", "O32"): {
        "grid": "O32",
        "projection": {"south_pole": [-30, 40], "type": "rotation"},
    },
}

# Not clipped, but worth a note before CLIPPED is mistaken for the whole story of the regional
# sources: 12 nearest-neighbour pairs onto H8-ring return the WHOLE target grid with 698-767 of
# the 768 values NaN (finding F5 of the probe) -- MIR's NaN-fill outside the source area, since
# a HEALPix target is never clipped. They are ordinary matrix pairs, not CLIPPED entries; the
# only consequence is that the nearest reference comparison in regrid_common must require the
# NaN positions of result and reference to be identical and compare the finite values only.
# The 12 sources: ICON-CH2, LAEA-EFAS-5km, O32-area, ac-area-n600, gg-area-N32, ll-area-0.1,
# ll-area-string, ll-rotated-area, swisslv03, swisslv95, xy-laea-grid-bbox-reverse, xy-swisslv03.

_SOURCES_BY_ID = {s.id: s for s in SOURCES}
_TARGETS_BY_ID = {t.id: t for t in TARGETS}


def matrix(interpolation: str) -> Iterator[tuple[Source, Target]]:
    """Yield the ``(source, target)`` pairs of one interpolation, sources in table order.

    Pure table walking: no grid is built and no coordinate is computed here, so the matrix can
    be parametrised without touching eckit. The pairs are the covered sources and the targets
    of the selection of ``interpolation``, in the order they are listed there.
    """
    try:
        source_ids, target_ids = SELECTION[interpolation]
    except KeyError:
        raise KeyError(
            f"unknown interpolation {interpolation!r}: SELECTION has {sorted(SELECTION)}, "
            "add it there (and to INTERPOLATIONS in regrid_common.py)"
        ) from None

    for source_id in source_ids:
        source = _SOURCES_BY_ID[source_id]
        for target_id in target_ids:
            yield source, _TARGETS_BY_ID[target_id]


def check_table() -> list[str]:
    """Check the table against the grids: return a list of problems, empty when it is right.

    Run by hand or from a test, never at import: this is the only function here which touches
    ``earthkit.geo``. Builds ``Grid(spec)`` for EVERY source, covered or not, and for every
    target, and compares ``shape`` and ``size()`` with the table; a grid which does not build is
    a problem. Only metadata is read, so even the grids above :data:`MAX_POINTS` are cheap here.
    The rest of the checks are pure table walking: unique ids,
    file-system-safe ids, a reason for every uncovered source and none for a covered one, and a
    selection of ids which exist and are covered. Never call ``to_latlons()`` or read
    ``order`` from here.
    """
    from earthkit.geo.grids import Grid

    problems = []

    def check(id, spec, shape, size):
        """Compare a built grid with the table; the build itself may raise."""
        grid = Grid(spec)
        if tuple(grid.shape) != tuple(shape):
            problems.append(f"{id}: shape {tuple(grid.shape)}, the table says {tuple(shape)}")
        if grid.size() != size:
            problems.append(f"{id}: size() {grid.size()}, the table says {size}")

    for source in SOURCES:
        if source.size != math.prod(source.shape):
            problems.append(f"{source.id}: size {source.size} is not the product of shape {source.shape}")

        if source.covered:
            try:
                check(source.id, source.spec, source.shape, source.size)
            except Exception as error:  # every grid of the table must build: report, do not crash
                problems.append(f"{source.id}: does not build: {error}")

    for target in TARGETS:
        try:
            check(target.id, target.spec, target.shape, math.prod(target.shape))
        except Exception as error:  # every target must build: report, do not crash
            problems.append(f"{target.id}: does not build: {error}")

    for source in SOURCES:
        # An uncovered grid must still build, and still have the shape and size the table claims.
        # Not being a source only means its coordinates cannot be used (no points, no lat/lon,
        # too many points to materialise): the grid itself exists, and building it reads the
        # metadata the catalogues serve offline. A grid that does not build is a problem
        # whatever reason the table gives for not covering it -- that is how a typo in a spec or
        # a grid a catalogue dropped stops being invisible.
        if source.covered:
            continue
        try:
            check(source.id, source.spec, source.shape, source.size)
        except Exception as error:  # every grid of the table must build: report, do not crash
            problems.append(f"{source.id}: does not build: {error}")

    problems.extend(_check_table_itself())
    return problems


def _check_table_itself() -> list[str]:
    """The checks which need no grid: ids, reasons and the consistency of the selections."""
    problems = []
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.")

    for kind, entries in (("source", SOURCES), ("target", TARGETS)):
        ids = [entry.id for entry in entries]
        for id in ids:
            if set(id) - allowed:
                problems.append(f"{kind} {id!r}: not a file system safe id ({sorted(set(id) - allowed)})")
        for id in sorted({id for id in ids if ids.count(id) > 1}):
            problems.append(f"{id}: {kind} id used twice")

    for source in SOURCES:
        if source.covered and source.reason:
            problems.append(f"{source.id}: covered, but it has a reason {source.reason!r}")
        if not source.covered and not source.reason:
            problems.append(f"{source.id}: not covered without a reason")
        if source.region not in (GLOBAL, REGIONAL):
            problems.append(f"{source.id}: unknown region {source.region!r}")

    if NOT_COVERED != {s.id: s.reason for s in SOURCES if not s.covered}:
        problems.append("NOT_COVERED does not match the uncovered entries of SOURCES")

    for interpolation, (source_ids, target_ids) in SELECTION.items():
        for id in source_ids:
            source = _SOURCES_BY_ID.get(id)
            if source is None:
                problems.append(f"{interpolation}: selection names unknown source {id!r}")
            elif not source.covered:
                problems.append(f"{interpolation}: selection names uncovered source {id!r}")
        for id in target_ids:
            if id not in _TARGETS_BY_ID:
                problems.append(f"{interpolation}: selection names unknown target {id!r}")
        if len(set(source_ids)) != len(source_ids) or len(set(target_ids)) != len(target_ids):
            problems.append(f"{interpolation}: selection lists an id twice")

    nearest_sources, nearest_targets = SELECTION["nearest-neighbour"]
    if tuple(nearest_sources) != tuple(s.id for s in SOURCES if s.covered):
        problems.append("nearest-neighbour does not cover every covered source")
    if tuple(nearest_targets) != tuple(t.id for t in TARGETS):
        problems.append("nearest-neighbour does not cover every target")

    for key in BROKEN | CLIPPED:
        interpolation, source_id, target_id = key
        if interpolation not in SELECTION:
            problems.append(f"{key}: BROKEN/CLIPPED for unknown interpolation {interpolation!r}")
            continue
        source_ids, target_ids = SELECTION[interpolation]
        if source_id not in source_ids or target_id not in target_ids:
            problems.append(f"{key}: BROKEN/CLIPPED combination is not in the selection")
    for key, reason in BROKEN.items():
        if not str(reason).startswith("BROKEN: "):
            problems.append(f"{key}: BROKEN reason does not start with 'BROKEN: '")

    return problems
