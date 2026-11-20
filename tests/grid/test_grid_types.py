# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Grid-object tests for all non-named grid types built from gridspec dicts.

Every expected value is derived independently from the formulas and rules in
``docs/concepts/grid/`` (regular_ll: shape = (180/dlat + 1, 360/dlon) globally,
area rows/cols = span/step + 1; F<N>: 2N Gaussian latitudes x 4N longitudes, the
latitudes being the roots of the Legendre polynomial of degree 2N and the clipped
area following from them -- see :func:`gaussian_latitudes` and :func:`gg_area_nodes`;
O<N>: 4N^2 + 36N; H<N>: 12N^2; arakawa_c: cells 90/n, T/U/V reference staggering;
arakawa_c_um: arakawa_c with grid_factor [2, 4/3] and order i+j+; canonical-spec
normalisation rules), never copied from whatever the code returns. Where a value is
not independently derivable (projected bounding boxes, uids and other content hashes,
the tabulated ``pl`` of the classical ``N<N>`` grids, how a projection is normalised
in a canonical spec) that check is skipped for the case (``None`` in the case field,
with a comment saying why); no value observed from the code is written into an
expectation.

Cases the T2 probe (``.ecmwf-multi-agents/notes/grid-probe.md``) found broken are
``xfail(run=False)`` parameters whose reason starts with ``BROKEN: `` -- they are
listed but never executed, see test names below and the report. Each of them carries
the complete data its test needs, so that lifting the mark checks the intended
behaviour instead of failing on a missing piece of case data. ``.order`` is never
read for grids where it aborts the interpreter (F1: unstructured_ll; also
sh/ORCA/ICON/FESOM, covered in their own files). ``to_latlons()`` is only called on
grids below 1,000,000 points (resource gate); for those grids only the length is
compared with ``size()``, except where the coordinates are derived above. All tests
are marked ``long_test``.
"""

import math

import numpy as np
import pytest

from earthkit.geo.grids import Grid

pytestmark = pytest.mark.long_test

# Reference bbox of a grid with no regional area (and the placeholder that
# area-less grids such as Gaussian/HEALPix/reduced/unstructured report, F5).
GLOBAL_BBOX = (90.0, 0.0, -90.0, 360.0)

# Point-count resource gate: to_latlons() materialises the coordinates.
SIZE_GATE = 1_000_000

_MISSING = object()


def case(
    id,
    spec,
    type,
    order,
    shape,
    size,
    bbox=GLOBAL_BBOX,
    canonical=_MISSING,
    canon="full",
    rt=True,
    latlons="yes",
    broken=None,
):
    """A parametrised grid case.

    ``canonical`` is the independently derived canonical ``g.spec`` (defaults to
    ``spec``); ``canon`` selects the spec check: ``full`` (dict equality),
    ``pl`` (only the ``pl`` list is preserved; the fabricated name key is the
    BROKEN round-trip case below), ``uid`` (unstructured: keys and uid only, the
    uid is a content hash that cannot be written by hand). ``None`` in
    ``canonical``/``order``/``shape``/``size``/``bbox`` skips that check (value
    not independently derivable -- see the per-case comments). ``latlons`` is
    ``yes``, ``no`` (grid has no geographic meaning) or ``big`` (resource gate).
    ``broken`` is the ``BROKEN: ...`` reason of a case that does not work: the
    parameter becomes ``xfail(run=False)``, listed but never executed, with all the
    data the check needs so that lifting the mark tests the intended behaviour.
    """
    if canonical is _MISSING:
        canonical = spec
    marks = [] if broken is None else [pytest.mark.xfail(run=False, reason=broken)]
    return pytest.param(
        dict(
            spec=spec,
            canonical=canonical,
            canon=canon,
            type=type,
            order=order,
            shape=shape,
            size=size,
            bbox=bbox,
            rt=rt,
            latlons=latlons,
        ),
        id=id,
        marks=marks,
    )


def check_eq(actual, expected):
    """Compare nested specs/bboxes; floats with a tolerance, everything else exactly."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict)
        assert set(actual) == set(expected)
        for key, value in expected.items():
            check_eq(actual[key], value)
    elif isinstance(expected, (list, tuple)):
        assert len(actual) == len(expected)
        for a, e in zip(actual, expected):
            check_eq(a, e)
    elif isinstance(expected, float):
        assert actual == pytest.approx(expected)
    else:
        assert actual == expected


def gaussian_latitudes(n):
    """The 2N Gaussian latitudes of an ``F<N>``/``O<N>`` grid, north to south, in degrees.

    Docs (``types/gaussian.rst``): the latitudes are the Gaussian quadrature points,
    the roots of the Legendre polynomial of degree 2N, where N is the number of
    latitudes between a pole and the equator. ``leggauss`` returns those roots in
    ascending order as values of sin(lat), so the latitude is ``asin(root)``.
    """
    nodes, _weights = np.polynomial.legendre.leggauss(2 * n)
    return sorted((math.degrees(math.asin(float(x))) for x in nodes), reverse=True)


def gaussian_longitudes(n):
    """The 4N longitudes of the global ``F<N>`` grid, in degrees east of Greenwich.

    Docs: an ``F<N>`` grid has 4N points on every latitude, so the increment is
    360/(4N) = 90/N; as for ``regular_ll`` the eastmost point repeats the westmost
    one, so there is no 360 (``types/regular_ll.rst``).
    """
    dlon = 90.0 / n
    return [i * dlon for i in range(4 * n)]


def wrap_lon(lon):
    """A longitude normalised onto ``[-180, 180)`` (docs: ``shared/area.rst``)."""
    return (float(lon) + 180.0) % 360.0 - 180.0


def gg_area_nodes(n, area):
    """The points of the global ``F<N>`` grid inside ``area: [north, west, south, east]``.

    Docs (``shared/area.rst``): "for regular latitude-longitude and Gaussian grids,
    the area selects the points of the global grid inside it", with the west edge
    normalised into ``[-180, 180)``. Returns ``(latitudes, longitudes)`` north to
    south and west to east; the clipped grid's own area is the box spanned by them,
    which is what its canonical spec and ``bounding_box()`` must report.
    """
    north, west, south, east = (float(v) for v in area)
    tol = 1e-9  # nodes are exact; this only absorbs rounding at the area edges
    lats = [x for x in gaussian_latitudes(n) if south - tol <= x <= north + tol]
    lons = sorted(wrap_lon(x) for x in gaussian_longitudes(n) if west - tol <= wrap_lon(x) <= east + tol)
    assert lats and lons, f"F{n}: the area {area} leaves no grid points"
    return lats, lons


REGULAR_LL_CASES = [
    # global: rows = 180/dlat + 1, cols = 360/dlon (the eastmost longitude repeats the westmost)
    case("ll-1x1", {"grid": [1, 1]}, "regular_ll", "i+j-", (181, 360), 181 * 360),
    case("ll-0.5x0.5", {"grid": [0.5, 0.5]}, "regular_ll", "i+j-", (361, 720), 361 * 720),
    # grid is [dlon, dlat]: 360/1.5 = 240 columns, 180/0.5 + 1 = 361 rows
    case("ll-1.5x0.5", {"grid": [1.5, 0.5]}, "regular_ll", "i+j-", (361, 240), 361 * 240),
    # area [N, W, S, E]; both spans inclusive: (72-25)/0.1 + 1 = 471, (45+25)/0.1 + 1 = 701
    case(
        "ll-area-0.1",
        {"grid": [0.1, 0.1], "area": [72, -25, 25, 45]},
        "regular_ll",
        "i+j-",
        (471, 701),
        471 * 701,
        bbox=(72.0, -25.0, 25.0, 45.0),
    ),
    # string areas normalise to lists (docs: "n/w/s/e"); (72-25)/1 + 1 = 48, (45+25)/1 + 1 = 71
    case(
        "ll-area-string",
        {"grid": [1, 1], "area": "72/-25/25/45"},
        "regular_ll",
        "i+j-",
        (48, 71),
        48 * 71,
        canonical={"area": [72, -25, 25, 45], "grid": [1, 1]},
        bbox=(72.0, -25.0, 25.0, 45.0),
    ),
    # reference [0.5, 0.5]: latitudes 89.5..-89.5 -> 180 rows, longitudes 0.5..359.5 -> 360 cols
    case(
        "ll-reference",
        {"grid": [1, 1], "reference": [0.5, 0.5]},
        "regular_ll",
        "i+j-",
        (180, 360),
        180 * 360,
    ),
    case("ll-order-i+j+", {"grid": [1, 1], "order": "i+j+"}, "regular_ll", "i+j+", (181, 360), 181 * 360),
    case("ll-order-i-j-", {"grid": [1, 1], "order": "i-j-"}, "regular_ll", "i-j-", (181, 360), 181 * 360),
    case("ll-order-i-j+", {"grid": [1, 1], "order": "i-j+"}, "regular_ll", "i-j+", (181, 360), 181 * 360),
    # rotation [lat, lon] of the south pole normalises to projection/south_pole [lon, lat] (docs)
    case(
        "ll-rotation",
        {"grid": [1, 1], "rotation": [40, -30]},
        "regular_ll",
        "i+j-",
        (181, 360),
        181 * 360,
        canonical={"grid": [1, 1], "projection": {"type": "rotation", "south_pole": [-30, 40]}},
    ),
    case(
        "ll-rotated_ll-alias",
        {"type": "rotated_ll", "grid": [1, 1], "rotation": [30, 90]},
        "regular_ll",
        "i+j-",
        (181, 360),
        181 * 360,
        canonical={"grid": [1, 1], "projection": {"type": "rotation", "south_pole": [90, 30]}},
    ),
    # docs example: rotated fine-area grid; the geographic bbox of a rotated area
    # is not derivable by hand, so bbox is not checked here.
    case(
        "ll-rotated-area",
        {
            "grid": [0.02, 0.02],
            "area": [3.36, -6.82, -4.42, 4.8],
            "order": "i+j+",
            "projection": {"type": "rotation", "south_pole": [10, -43]},
        },
        "regular_ll",
        "i+j+",
        # in the rotated frame: (3.36+4.42)/0.02 + 1 = 390, (4.8+6.82)/0.02 + 1 = 582
        (390, 582),
        390 * 582,
        bbox=None,
    ),
    case(
        "ll-figure-sun",
        {"grid": [2, 2], "figure": "sun"},
        "regular_ll",
        "i+j-",
        (91, 180),
        91 * 180,
        canonical={"figure": "sun", "grid": [2, 2]},
    ),
]

ARAKAWA_C_CASES = [
    # n = cells per 90 deg -> d = 90/n = 1; T: half-cell reference in x and y
    case(
        "ac-T-n90",
        {"type": "arakawa_c", "n": 90},
        "regular_ll",  # canonicalises as the regular_ll it builds (docs)
        "i+j-",
        (180, 360),
        180 * 360,
        canonical={"grid": [1, 1], "reference": [0.5, 0.5]},
    ),
    # U: nodes on cell edges in x (offset 0), centres in y
    case(
        "ac-U-n90",
        {"type": "arakawa_c", "n": 90, "arrangement": "U"},
        "regular_ll",
        "i+j-",
        (180, 360),
        180 * 360,
        canonical={"grid": [1, 1], "reference": [0, 0.5]},
    ),
    # V: centres in x, nodes in y -> rows 180/1 + 1 = 181 including both poles
    case(
        "ac-V-n90",
        {"type": "arakawa_c", "n": 90, "arrangement": "V"},
        "regular_ll",
        "i+j-",
        (181, 360),
        181 * 360,
        canonical={"grid": [1, 1], "reference": [0.5, 0]},
    ),
    # grid_factor [2, 1.5] -> d = [2, 1.5]; rows = (180 - 1.5)/1.5 + 1 = 120, cols = 360/2 = 180
    case(
        "ac-grid_factor-n90",
        {"type": "arakawa_c", "n": 90, "grid_factor": [2, 1.5]},
        "regular_ll",
        "i+j-",
        (120, 180),
        120 * 180,
        canonical={"grid": [2, 1.5], "reference": [1, 0.75]},
    ),
    # n=600 -> d = 0.15; area clipped to the grid points inside it:
    # lat points are 0.075 mod 0.15 from 90 -> 72 -> 71.925, 25 -> 25.125;
    # lon points -25 -> -24.975, 45 -> 44.925; rows (71.925-25.125)/0.15 + 1 = 313,
    # cols (44.925+24.975)/0.15 + 1 = 467
    case(
        "ac-area-n600",
        {"type": "arakawa_c", "n": 600, "area": [72, -25, 25, 45]},
        "regular_ll",
        "i+j-",
        (313, 467),
        313 * 467,
        canonical={"area": [71.925, -24.975, 25.125, 44.925], "grid": [0.15, 0.15], "reference": [0.075, 0.075]},
        bbox=(71.925, -24.975, 25.125, 44.925),
    ),
    case(
        "ac-order-i+j+-n90",
        {"type": "arakawa_c", "n": 90, "order": "i+j+"},
        "regular_ll",
        "i+j+",
        (180, 360),
        180 * 360,
        canonical={"grid": [1, 1], "order": "i+j+", "reference": [0.5, 0.5]},
    ),
    # rotated arakawa_c: the canonical spec gains an area adjustment (89.617...) in the
    # rotated frame that is not hand-derivable, so spec and bbox are not checked here;
    # the underlying grid is still the T-arranged 1-degree grid: 180 x 360
    case(
        "ac-rotation-n90",
        {"type": "arakawa_c", "n": 90, "rotation": [40, -30]},
        "regular_ll",
        "i+j-",
        (180, 360),
        180 * 360,
        canonical=None,
        bbox=None,
    ),
]

ARAKAWA_C_UM_CASES = [
    # arakawa_c with grid_factor [2, 4/3] and order i+j+ (docs):
    # n=96 -> d = 90/96 = 0.9375; [dx, dy] = [1.875, 1.25]
    # rows = 180/1.25 = 144 (cell centres), cols = 360/1.875 = 192
    case(
        "ac-um-T-N96",
        {"type": "arakawa_c_um", "N": 96},
        "regular_ll",
        "i+j+",
        (144, 192),
        144 * 192,
        canonical={"grid": [1.875, 1.25], "order": "i+j+", "reference": [0.9375, 0.625]},
    ),
    case(
        "ac-um-U-N96",
        {"type": "arakawa_c_um", "N": 96, "arrangement": "U"},
        "regular_ll",
        "i+j+",
        (144, 192),
        144 * 192,
        canonical={"grid": [1.875, 1.25], "order": "i+j+", "reference": [0, 0.625]},
    ),
    # V: nodes in y -> 180/1.25 + 1 = 145 rows
    case(
        "ac-um-V-N96",
        {"type": "arakawa_c_um", "N": 96, "arrangement": "V"},
        "regular_ll",
        "i+j+",
        (145, 192),
        145 * 192,
        canonical={"grid": [1.875, 1.25], "order": "i+j+", "reference": [0.9375, 0]},
    ),
    # n=48 -> d = 1.875; [dx, dy] = [3.75, 2.5]; rows 72, cols 96
    case(
        "ac-um-T-N48",
        {"type": "arakawa_c_um", "N": 48},
        "regular_ll",
        "i+j+",
        (72, 96),
        72 * 96,
        canonical={"grid": [3.75, 2.5], "order": "i+j+", "reference": [1.875, 1.25]},
    ),
]

# The F32 grid clipped to the area [80, -20, 20, 40], derived by gg_area_nodes: the
# rows it keeps are the roots of P_64 between 20 and 80 degrees, the columns the
# multiples of 360/(4*32) = 2.8125 whose longitude wrapped into [-180, 180) lies
# inside [-20, 40], and the grid's own area is bounded by the outermost of them
GG_N = 32
GG_AREA = [80, -20, 20, 40]
GG_LATS, GG_LONS = gg_area_nodes(GG_N, GG_AREA)
GG_ROWS, GG_COLS = len(GG_LATS), len(GG_LONS)
GG_SNAPPED = (GG_LATS[0], GG_LONS[0], GG_LATS[-1], GG_LONS[-1])

REGULAR_GG_CASES = [
    # F<N>: 2N latitudes x 4N longitudes
    case(
        "gg-N16",
        {"type": "regular_gg", "N": 16},
        "regular_gg",
        "i+j-",
        (32, 64),
        32 * 64,
        canonical={"grid": "F16"},
    ),
    case(
        "gg-N32",
        {"type": "regular_gg", "N": 32},
        "regular_gg",
        "i+j-",
        (64, 128),
        64 * 128,
        canonical={"grid": "F32"},
    ),
    case(
        "gg-F48",
        {"grid": "F48"},
        "regular_gg",
        "i+j-",
        (96, 192),
        96 * 192,
    ),
    case(
        "gg-rotated-N16",
        {"type": "rotated_gg", "N": 16, "rotation": [40, -30]},
        "regular_gg",
        "i+j-",
        (32, 64),
        32 * 64,
        canonical={"grid": "F16", "projection": {"type": "rotation", "south_pole": [-30, 40]}},
    ),
    # area clip on F32: rows, columns and the snapped area edges come from the Gaussian
    # nodes and the longitude grid (GG_* above), not from anything the code returned
    case(
        "gg-area-N32",
        {"type": "regular_gg", "N": GG_N, "area": list(GG_AREA)},
        "regular_gg",
        "i+j-",
        (GG_ROWS, GG_COLS),
        GG_ROWS * GG_COLS,
        canonical={"area": list(GG_SNAPPED), "grid": f"F{GG_N}"},
        bbox=GG_SNAPPED,
    ),
    case(
        "gg-order-i+j+-N16",
        {"type": "regular_gg", "N": 16, "order": "i+j+"},
        "regular_gg",
        "i+j+",  # regular_gg honours order (unlike reduced_gg, see the BROKEN case)
        (32, 64),
        32 * 64,
        canonical={"grid": "F16", "order": "i+j+"},
    ),
    case(
        "gg-F1600",
        {"grid": "F1600"},
        "regular_gg",
        "i+j-",
        (3200, 6400),
        3200 * 6400,
        latlons="big",
    ),
]


def o_size(n):
    """Octahedral reduced Gaussian grid: 4N^2 + 36N points (docs)."""
    return 4 * n * n + 36 * n


REDUCED_GG_CASES = [
    case("rd-O8", {"grid": "O8"}, "reduced_gg", "i+j-", (o_size(8),), o_size(8)),
    case("rd-O16", {"grid": "O16"}, "reduced_gg", "i+j-", (o_size(16),), o_size(16)),
    case("rd-O32", {"grid": "O32"}, "reduced_gg", "i+j-", (o_size(32),), o_size(32)),
    case("rd-O64", {"grid": "O64"}, "reduced_gg", "i+j-", (o_size(64),), o_size(64)),
    case("rd-O128", {"grid": "O128"}, "reduced_gg", "i+j-", (o_size(128),), o_size(128)),
    case("rd-O1280", {"grid": "O1280"}, "reduced_gg", "i+j-", (o_size(1280),), o_size(1280), latlons="big"),
    # classical N<N> pl is tabulated (docs), no formula: sizes/shapes are not checked
    # independently, only consistency and the point count vs to_latlons()
    case("rd-N16", {"grid": "N16"}, "reduced_gg", "i+j-", None, None),
    case("rd-N32", {"grid": "N32"}, "reduced_gg", "i+j-", None, None),
    case("rd-N64", {"grid": "N64"}, "reduced_gg", "i+j-", None, None),
    # the ``N: <N>`` shortcut builds the OCTAHEDRAL grid, not the classical one (T2/F12;
    # the docs shortcut table is ambiguous about this)
    case(
        "rd-typeN32->O32",
        {"type": "reduced_gg", "N": 32},
        "reduced_gg",
        "i+j-",
        (o_size(32),),
        o_size(32),
        canonical={"grid": "O32"},
    ),
    case(
        "rd-O32-rotated",
        {"grid": "O32", "rotation": [40, -30]},
        "reduced_gg",
        "i+j-",
        (o_size(32),),
        o_size(32),
        canonical={"grid": "O32", "projection": {"type": "rotation", "south_pole": [-30, 40]}},
    ),
    case(
        "rd-reduced_rotated_gg-alias",
        {"type": "reduced_rotated_gg", "N": 32, "rotation": [40, -30]},
        "reduced_gg",
        "i+j-",
        (o_size(32),),
        o_size(32),
        canonical={"grid": "O32", "projection": {"type": "rotation", "south_pole": [-30, 40]}},
    ),
    # reduced_gg area: every latitude of a reduced grid has its own longitude increment
    # and its own placement inside [W, E] (docs: "pl[j] points on latitude j"), and the
    # docs give no rule for aligning them with the area, so the number of points that
    # survive the clip cannot be derived here; shape, size, spec and bbox are skipped
    case(
        "rd-O32-area",
        {"grid": "O32", "area": [80, -20, 20, 40]},
        "reduced_gg",
        "i+j-",
        None,
        None,
        canonical=None,
        bbox=None,
    ),
    # arbitrary pl: point count = sum(pl) (docs: pl[j] points on latitude j); the
    # canonical spec wrongly adds a fabricated classical name (see BROKEN cases)
    case(
        "rd-pl-8lat",
        {"type": "reduced_gg", "pl": [4, 8, 12, 16, 16, 12, 8, 4]},
        "reduced_gg",
        "i+j-",
        (sum([4, 8, 12, 16, 16, 12, 8, 4]),),
        sum([4, 8, 12, 16, 16, 12, 8, 4]),
        canon="pl",
        rt=False,
    ),
    case(
        "rd-pl-4lat",
        {"type": "reduced_gg", "pl": [8, 8, 8, 8]},
        "reduced_gg",
        "i+j-",
        (32,),
        32,
        canon="pl",
        rt=False,
    ),
]

HEALPIX_CASES = [
    # H<Nside>: 12 * Nside^2 pixel centres
    case(
        "hp-Nside4",
        {"type": "healpix", "Nside": 4},
        "HEALPix",  # type strings are inconsistently cased (F13); match the docs section name
        "ring",
        (12 * 4 * 4,),
        12 * 4 * 4,
        canonical={"grid": "H4"},
    ),
    case("hp-H8", {"grid": "H8"}, "HEALPix", "ring", (12 * 8 * 8,), 12 * 8 * 8),
    case(
        "hp-H8-nested",
        {"grid": "H8", "order": "nested"},
        "HEALPix",
        "nested",
        (12 * 8 * 8,),
        12 * 8 * 8,
    ),
    # HR<N>/HN<N> names normalise to {grid: H<N>[, order: nested]} (docs)
    case(
        "hp-HR16",
        {"grid": "HR16"},
        "HEALPix",
        "ring",
        (12 * 16 * 16,),
        12 * 16 * 16,
        canonical={"grid": "H16"},
    ),
    case(
        "hp-HN16",
        {"grid": "HN16"},
        "HEALPix",
        "nested",
        (12 * 16 * 16,),
        12 * 16 * 16,
        canonical={"grid": "H16", "order": "nested"},
    ),
    case(
        "hp-typeHEALPix-Nside12",
        {"type": "HEALPix", "Nside": 12},
        "HEALPix",
        "ring",
        (12 * 12 * 12,),
        12 * 12 * 12,
        canonical={"grid": "H12"},
    ),
    case("hp-H31", {"grid": "H31"}, "HEALPix", "ring", (12 * 31 * 31,), 12 * 31 * 31),
    case("hp-H1024", {"grid": "H1024"}, "HEALPix", "ring", (12 * 1024 * 1024,), 12 * 1024 * 1024, latlons="big"),
]

REGULAR_XY_CASES = [
    # x/y ranges [start, stop, step] with stop included -> 11 points each
    case(
        "xy-plain",
        {"x": [0, 10, 1], "y": [0, 10, 1]},
        "regular_xy",
        "i+j+",  # both ranges ascend (docs: the directions determine the order)
        (11, 11),
        121,
        canonical={"order": "i+j+", "x": [0, 10, 1], "y": [0, 10, 1]},
        latlons="no",  # plain Cartesian, to_latlons() raises (see the negative test)
    ),
    case(
        "xy-descending-x",
        {"x": [10, 0, -1], "y": [0, 10, 1]},
        "regular_xy",
        "i-j+",
        (11, 11),
        121,
        canonical={"order": "i-j+", "x": [10, 0, -1], "y": [0, 10, 1]},
        latlons="no",
    ),
    case(
        "xy-order-i+j+",
        {"x": [0, 10, 1], "y": [0, 10, 1], "order": "i+j+"},
        "regular_xy",
        "i+j+",
        (11, 11),
        121,
        canonical={"order": "i+j+", "x": [0, 10, 1], "y": [0, 10, 1]},
        latlons="no",
    ),
    # false_xy without a projection: docs say the canonical form keeps a non-default
    # false_xy, but the key is dropped and the uid ignores it -> spec not asserted here,
    # see the BROKEN case in test_false_xy_dropped_from_spec_broken
    case(
        "xy-false_xy",
        {"x": [0, 100, 10], "y": [0, 100, 10], "false_xy": [1000, 2000]},
        "regular_xy",
        "i+j+",
        (11, 11),
        121,
        canonical=None,
        latlons="no",
    ),
    # grid+bbox form with a reverse projection: 201 x 161 points (step included at both ends);
    # y is normalised north-to-south with order i+j-; the projected bbox is not hand-derivable
    case(
        "xy-laea-grid-bbox-reverse",
        {
            "type": "regular_xy",
            "grid": [5000, 5000],
            "bounding_box_xy": [-500000, -400000, 500000, 400000],
            "projection": {"type": "laea", "lat_0": 52, "lon_0": 10, "reverse": True},
        },
        "regular_xy",
        "i+j-",
        (161, 201),
        161 * 201,
        canonical={
            "projection": {"lat_0": 52, "lon_0": 10, "reverse": True, "type": "laea"},
            "x": [-500000, 500000, 5000],
            "y": [400000, -400000, -5000],
        },
        bbox=None,
    ),
    # named projection: cols (960000-260000)/50000 + 1 = 15, rows (450000+150000)/50000 + 1 = 13
    case(
        "xy-swisslv03",
        {"projection": "swisslv03", "x": [260000, 960000, 50000], "y": [-150000, 450000, 50000]},
        "regular_xy",
        "i+j+",
        (13, 15),
        13 * 15,
        canonical={
            "order": "i+j+",
            "projection": "swisslv03",
            "x": [260000, 960000, 50000],
            "y": [-150000, 450000, 50000],
        },
        bbox=None,
    ),
]

UNSTRUCTURED_CASES = [
    # canonical form {type: unstructured_ll, uid: ...} (docs); the uid is a content
    # hash, so only its presence/consistency is checked (canon="uid")
    case(
        "uns-latitudes-longitudes",
        {"latitudes": [0.0, 10.5, -20.25, 45.0], "longitudes": [0.0, 90.0, 180.0, 270.5]},
        "unstructured_ll",
        None,  # reading .order aborts the process (F1), see test_unstructured_order_aborts_broken
        (4,),
        4,
        canon="uid",
    ),
    case(
        "uns-unstructured-alias",
        {"type": "unstructured", "latitudes": [0.0, 10.5, -20.25, 45.0], "longitudes": [0.0, 90.0, 180.0, 270.5]},
        "unstructured_ll",
        None,
        (4,),
        4,
        canon="uid",
    ),
    # a named unstructured grid normalises to {grid: <name>} (docs) and resolves in-process
    case(
        "uns-named",
        {
            "type": "unstructured_ll",
            "latitudes": [0.0, 1.0, 2.0, 3.0, 4.0],
            "longitudes": [10.0, 11.0, 12.0, 13.0, 14.0],
            "name": "EARTHKIT-GEO-T3-UNS",
        },
        "unstructured_ll",
        None,
        (5,),
        5,
        canonical={"grid": "EARTHKIT-GEO-T3-UNS"},
    ),
]

ALL_CASES = (
    REGULAR_LL_CASES
    + ARAKAWA_C_CASES
    + ARAKAWA_C_UM_CASES
    + REGULAR_GG_CASES
    + REDUCED_GG_CASES
    + HEALPIX_CASES
    + REGULAR_XY_CASES
    + UNSTRUCTURED_CASES
)

# Specs the T2 probe's "Broken" list found do not build at all, but should. They stay in
# the members parametrisation as xfail(run=False) -- listed, never executed -- and carry
# the complete expected data, so that lifting the mark re-checks the grid members instead
# of failing on a missing piece of case data. The two remaining BROKEN specs of that list
# have expectations the members test cannot express; they have their own tests below.
BROKEN_BUILD_CASES = [
    # N4 is the name the library itself invents for an arbitrary 8-latitude pl, and then
    # refuses to read; no pl (hence no shape/size, and no canonical spec) is derivable
    # for a grid that is neither octahedral nor one of the tabulated classical ones
    case(
        "BROKEN-construct:rd-N4-name",
        {"grid": "N4"},
        "reduced_gg",
        "i+j-",  # the default order of a reduced Gaussian grid (docs)
        None,
        None,
        canonical=None,
        broken="BROKEN: (F2/F6) the build fails with RuntimeError 'Bad value: "
        "reduced_classical_pl: unknown N=4'. N4 is not a classical Gaussian grid (docs tabulate "
        "N=16, 24, ...), but it is the name the library INVENTS in the canonical spec of an "
        "arbitrary 8-latitude pl grid (see the BROKEN round trips below): it rejects a spec it "
        "produced itself, and a name it generates must be a name it can read",
    ),
    # the documented grid + bounding_box_xy form of regular_xy: rows = (400000+400000)/5000
    # + 1 = 161, cols = (500000+500000)/5000 + 1 = 201 (docs: shared/area.rst gives
    # bounding_box_xy as [min_x, min_y, max_x, max_y]); how this form normalises the y
    # direction, the ranges and the projection is not derivable, and neither is the
    # geographic extent of a projected grid, so those checks are skipped
    case(
        "BROKEN-construct:xy-laea-no-reverse",
        {
            "type": "regular_xy",
            "grid": [5000, 5000],
            "bounding_box_xy": [-500000, -400000, 500000, 400000],
            "projection": {"type": "laea", "lat_0": 52, "lon_0": 10},
        },
        "regular_xy",
        None,
        (161, 201),
        161 * 201,
        canonical=None,
        bbox=None,
        broken="BROKEN: (F11, T2 Broken list) the build fails with RuntimeError 'SpecError: "
        "[RegularXY: 'projection' from (x, y) expected (e.g. with reverse: true)]'. "
        "docs/concepts/grid/types/projected.rst lists ``grid`` + ``bounding_box_xy`` as a form of "
        "regular_xy, so it must build the 161 x 201 point grid; if rejecting a built-in "
        "projection without 'reverse' is by design, the docs table has to say so",
    ),
]


def build(c):
    return Grid(c["spec"])


@pytest.mark.parametrize("c", ALL_CASES + BROKEN_BUILD_CASES)
def test_grid_members(c):
    grid = build(c)

    # canonical spec
    if c["canon"] == "full":
        if c["canonical"] is not None:
            check_eq(grid.spec, c["canonical"])
    elif c["canon"] == "pl":
        # only the pl list is expected to survive (docs: "Any other pl is accepted,
        # and kept in the canonical spec"); the fabricated name key is the BROKEN case
        assert grid.spec["pl"] == c["spec"]["pl"]
    elif c["canon"] == "uid":
        assert set(grid.spec) == {"type", "uid"}
        assert grid.spec["type"] == "unstructured_ll"
        assert grid.spec["uid"] == grid.uid

    assert grid.type == c["type"]

    # order: None skips the check -- either reading it aborts the process (F1, see the
    # unstructured cases) or the value is not derivable for that case (see its comment)
    if c["order"] is not None:
        assert grid.order == c["order"]

    if c["shape"] is not None:
        assert tuple(grid.shape) == c["shape"]
        if c["size"] is not None:
            assert c["size"] == math.prod(c["shape"])
    assert grid.size() == math.prod(grid.shape)

    if c["bbox"] is None:
        # not independently derivable (rotated/projected extents): shape of the value only
        bbox = grid.bounding_box()
        assert len(bbox) == 4
        assert bbox[0] >= bbox[2]
    else:
        check_eq(tuple(grid.bounding_box()), c["bbox"])

    # uid format and stability across a rebuild from the canonical spec
    assert isinstance(grid.uid, str)
    assert len(grid.uid) == 32

    if c["rt"]:
        rebuilt = Grid(grid.spec)
        assert rebuilt == grid
        assert rebuilt.uid == grid.uid


def test_case_data_is_complete():
    """Every parametrised case carries all the fields its test reads.

    A BROKEN case left in a parametrisation with only a spec would fail on a missing
    key as soon as the library could build it, instead of checking the behaviour it
    documents, so the case data is asserted to be complete here.
    """
    members = {"spec", "canonical", "canon", "type", "order", "shape", "size", "bbox", "rt", "latlons"}
    for param in ALL_CASES + BROKEN_BUILD_CASES:
        assert set(param.values[0]) == members, param.id
    for param in XY_NO_TYPE_CASES:
        assert {"spec", "type", "shape", "size"} <= set(param.values[0]), param.id
    for param in PL_ROUNDTRIP_CASES:
        assert "pl" in param.values[0], param.id
    for param in BROKEN_NAME_ERROR_CASES:
        assert {"spec", "hint"} == set(param.values[0]), param.id
    for param in BROKEN_UID_CASES + UNSTRUCTURED_ORDER_CASES:
        assert isinstance(param.values[0], dict), param.id


# --------------------------------------------------------------------------
# BROKEN cases (T2 probe): marked xfail(run=False), listed but never executed.
# --------------------------------------------------------------------------

# The two BROKEN specs whose intended behaviour test_grid_members cannot express keep
# their ids so that the report and the list above point at the same cases. One must
# build from a value only the library knows (a grid's own uid, taken at run time rather
# than copied from the probe); the other must not build at all and only has to give the
# right reason for refusing.
BROKEN_UID_CASES = [
    pytest.param(
        {"grid": [1, 1]},
        id="BROKEN-construct:ll-own-uid",
        marks=pytest.mark.xfail(
            run=False,
            reason="BROKEN: (F7) building a grid from its OWN uid fails with RuntimeError "
            "'GridUnknownError: [Grid: unknown grid uid ...]' -- only catalogue and unstructured grids are "
            "registered, so a regular grid's uid cannot be used to rebuild it, although docs list uid as a "
            "spec key ('a registered grid')",
        ),
    ),
]


@pytest.mark.parametrize("spec", BROKEN_UID_CASES)
def test_broken_build_from_own_uid(spec):
    # the uid is read from the grid itself, so nothing here comes from the probe output
    ref = Grid(spec)
    grid = Grid({"uid": ref.uid})
    assert grid == ref
    assert grid.uid == ref.uid
    assert grid.type == "regular_ll"
    assert tuple(grid.shape) == (181, 360)  # 180/1 + 1 rows by 360/1 columns (docs)
    assert Grid(grid.spec) == grid  # how that grid is then canonicalised is not derived here


BROKEN_NAME_ERROR_CASES = [
    pytest.param(
        dict(spec={"grid": "H0"}, hint=r"Nside|unknown|invalid|zero"),
        id="BROKEN-construct:hp-H0-name",
        marks=pytest.mark.xfail(
            run=False,
            reason="BROKEN: (F6) Grid({grid: H0}) fails with 'SpecError: [Grid: cannot build grid without "
            "'type']' -- the spec is a name, and a name never carries a 'type' key, so the message blames a "
            "missing key that is not needed; the real problem is Nside=0, since H<N> is 12*N*N points, so "
            "H0 is no grid at all",
        ),
    ),
]


@pytest.mark.parametrize("c", BROKEN_NAME_ERROR_CASES)
def test_broken_unknown_name_message(c):
    # a rejected name must say what is wrong with it (and not blame a missing 'type' key)
    with pytest.raises(RuntimeError, match=c["hint"]) as excinfo:
        Grid(c["spec"])
    assert "'type'" not in str(excinfo.value)


PL_ROUNDTRIP_CASES = [
    pytest.param(
        {"type": "reduced_gg", "pl": [4, 8, 12, 16, 16, 12, 8, 4]},
        id="BROKEN:rd-pl-8lat-roundtrip",
        marks=pytest.mark.xfail(
            run=False,
            reason="BROKEN: (F2) Grid(g.spec) raises RuntimeError 'Bad value: reduced_classical_pl: unknown N=4' "
            "-- the canonical spec is {'grid': 'N4', 'pl': [...]} but docs say a non-octahedral, "
            "non-classical grid is canonicalised with pl only, no fabricated classical name",
        ),
    ),
    pytest.param(
        {"type": "reduced_gg", "pl": [8, 8, 8, 8]},
        id="BROKEN:rd-pl-4lat-roundtrip",
        marks=pytest.mark.xfail(
            run=False,
            reason="BROKEN: (F2) Grid(g.spec) raises RuntimeError 'Bad value: reduced_classical_pl: unknown N=2' "
            "-- the canonical spec invents the non-classical name N2 (same cause as the 8-latitude pl)",
        ),
    ),
    pytest.param(
        {"type": "reduced_gg", "pl": [8, 8, 8, 8], "order": "i+j+"},
        id="BROKEN:rd-pl-4lat-order-roundtrip",
        marks=pytest.mark.xfail(
            run=False,
            reason="BROKEN: (F2+F3) the canonical spec drops order (it is silently ignored on reduced_gg) "
            "and adds the non-classical name N2, so Grid(g.spec) raises 'unknown N=2' (verified: the uid "
            "is identical to the same grid without order)",
        ),
    ),
]


@pytest.mark.parametrize("spec", PL_ROUNDTRIP_CASES)
def test_reduced_gg_pl_roundtrip_broken(spec):
    grid = Grid(spec)
    # docs: canonical form keeps pl for grids that are neither octahedral nor classical
    check_eq(grid.spec, {"pl": spec["pl"]})
    assert Grid(grid.spec) == grid


@pytest.mark.xfail(
    run=False,
    reason="BROKEN: (F3) order is silently ignored on reduced_gg -- Grid({grid: O32, order: i+j+}) reads "
    "order 'i+j-', its spec drops the key and its uid and points equal the default (verified in T2 and by "
    "hand); regular_gg honours order",
)
def test_reduced_gg_order_ignored_broken():
    grid = Grid({"grid": "O32", "order": "i+j+"})
    default = Grid({"grid": "O32"})
    # docs list four orders for reduced_gg; regular_gg honours order (tested above)
    assert grid.order == "i+j+"
    assert grid.uid != default.uid
    assert grid.spec == {"grid": "O32", "order": "i+j+"}


UNSTRUCTURED_ORDER_CASES = [
    pytest.param(
        {"latitudes": [0.0, 10.5, -20.25, 45.0], "longitudes": [0.0, 90.0, 180.0, 270.5]},
        id="BROKEN:uns-latitudes-longitudes-order",
    ),
    pytest.param(
        {"type": "unstructured", "latitudes": [0.0, 10.5, -20.25, 45.0], "longitudes": [0.0, 90.0, 180.0, 270.5]},
        id="BROKEN:uns-alias-order",
    ),
    pytest.param(
        {
            "type": "unstructured_ll",
            "latitudes": [0.0, 1.0, 2.0, 3.0, 4.0],
            "longitudes": [10.0, 11.0, 12.0, 13.0, 14.0],
            "name": "EARTHKIT-GEO-T3-UNS-ORDER",
        },
        id="BROKEN:uns-named-order",
    ),
]


@pytest.mark.xfail(
    run=False,
    reason="BROKEN: (F1) reading .order on unstructured_ll grids ABORTS the whole Python process "
    '(SIGABRT, rc -6; uncaught C++ eckit::NotImplemented "Grid: order() is not implemented for type '
    "'unstructured_ll'\") -- re-verified in a subprocess (rc=134); never read .order for this type "
    "in process. Every other member works",
)
@pytest.mark.parametrize("spec", UNSTRUCTURED_ORDER_CASES)
def test_unstructured_order_aborts_broken(spec):
    grid = Grid(spec)
    # reading .order must raise a Python exception, not abort the interpreter
    with pytest.raises(RuntimeError):
        grid.order


XY_NO_TYPE_CASES = [
    pytest.param(
        {
            "spec": {"grid": [2, 2], "bounding_box_xy": [0, 0, 10, 10]},
            "type": "regular_xy",
            "shape": (6, 6),  # dx=dy=2 over 0..10 inclusive -> 6 x 6 points
            "size": 36,
        },
        id="BROKEN:xy-grid-bbox-without-type",
    ),
    pytest.param(
        {
            "spec": {
                "grid": [10, 10],
                "bounding_box_xy": [-100, -100, 100, 100],
                "false_lonlat": [0, 52],
                "projection": {"type": "laea", "lat_0": 52, "lon_0": 0},
            },
            "type": "regular_xy",
            "shape": (21, 21),  # dx=dy=10 over -100..100 inclusive -> 21 x 21 points
            "size": 441,
        },
        id="BROKEN:xy-false_lonlat-without-type",
    ),
]


@pytest.mark.xfail(
    run=False,
    reason="BROKEN: (F10) {grid, bounding_box_xy} without type: regular_xy silently builds a plain "
    "regular_ll and drops bounding_box_xy: {grid: [2,2], ...} gives 91 x 180 ll instead of 6 x 6 xy; "
    "with false_lonlat/laea it gives a 19 x 36 regular_ll whose to_latlons() is nonsense (all lat "
    "~52.0008, metre steps in degrees). The explicit-type form works (see xy-laea-grid-bbox-reverse)",
)
@pytest.mark.parametrize("c", XY_NO_TYPE_CASES)
def test_xy_without_regular_xy_type_broken(c):
    grid = Grid(c["spec"])
    # docs: "Required (or grid and bounding_box_xy)" -- this form must build regular_xy
    assert grid.type == c["type"]
    assert tuple(grid.shape) == c["shape"]
    assert grid.size() == c["size"]


@pytest.mark.xfail(
    run=False,
    reason="BROKEN: (F4) HN10 (nested, Nside not a power of 2) builds, reports spec {grid: H10, order: nested}, "
    "1200 points and round-trips, but to_latlons() raises RuntimeError 'Assertion failed: Nside must be a "
    "power of 2' -- a deferred failure; the build should reject it (docs: nested requires Nside a power of 2)",
)
def test_hn10_to_latlons_broken():
    # nested order requires Nside to be a power of 2 (docs) -- this must not build
    grid = Grid({"grid": "HN10"})
    lat, lon = grid.to_latlons()
    assert len(lat) == grid.size() == 12 * 10 * 10


@pytest.mark.xfail(
    run=False,
    reason="BROKEN: (new T3 finding, not in the T2 list) false_xy without a projection is silently ignored -- "
    "the canonical spec drops the key (returns {order, x, y}) and the uid is IDENTICAL to the same grid "
    "without false_xy (verified in a subprocess), so the two grids cannot be distinguished by spec or uid; "
    "docs say the canonical form keeps a non-default false_xy",
)
def test_false_xy_dropped_from_spec_broken():
    plain = Grid({"x": [0, 100, 10], "y": [0, 100, 10]})
    grid = Grid({"x": [0, 100, 10], "y": [0, 100, 10], "false_xy": [1000, 2000]})
    # docs: canonical form is {x, y} "with false_xy, projection and order when not the defaults"
    assert grid.spec == {"false_xy": [1000, 2000], "order": "i+j+", "x": [0, 100, 10], "y": [0, 100, 10]}
    assert grid.uid != plain.uid


# --------------------------------------------------------------------------
# to_latlons(): point count vs size(), only for grids below the 1M gate, plus the
# Gaussian coordinates the clipped gg-area-N32 expectations are derived from.
# --------------------------------------------------------------------------


def latlons_cases():
    return [p for p in ALL_CASES if p.values[0]["latlons"] == "yes"]


@pytest.mark.parametrize("c", latlons_cases(), ids=[p.id for p in latlons_cases()])
def test_grid_to_latlons(c):
    grid = build(c)
    assert grid.size() < SIZE_GATE  # the resource gate
    lat, lon = grid.to_latlons()
    assert len(lat) == grid.size()
    assert len(lon) == grid.size()


@pytest.mark.parametrize("n", [8, 16, 32])
def test_regular_gg_points_are_gaussian_nodes(n):
    # the global F<N> grid really is the roots of P_2N by i*90/N degrees, which is what
    # gg_area_nodes() derives the clipped expectations from (docs: types/gaussian.rst)
    lat, lon = Grid({"type": "regular_gg", "N": n}).to_latlons()
    check_eq(sorted({float(v) for v in lat}, reverse=True), gaussian_latitudes(n))
    check_eq(sorted({float(v) for v in lon}), gaussian_longitudes(n))


def test_regular_gg_area_nodes():
    # the clipped grid holds exactly the rows and columns gg_area_nodes() keeps, so the
    # shape, area and bbox asserted for gg-area-N32 are the same derivation, not a copy
    lat, lon = Grid({"type": "regular_gg", "N": GG_N, "area": list(GG_AREA)}).to_latlons()
    check_eq(sorted({float(v) for v in lat}, reverse=True), GG_LATS)
    check_eq(sorted({wrap_lon(float(v)) for v in lon}), GG_LONS)


def test_unstructured_to_latlons_exact():
    # an unstructured grid returns exactly the points it was built from
    lats = [0.0, 10.5, -20.25, 45.0]
    lons = [0.0, 90.0, 180.0, 270.5]
    grid = Grid({"latitudes": lats, "longitudes": lons})
    lat, lon = grid.to_latlons()
    check_eq(list(lat), lats)
    check_eq(list(lon), lons)


# --------------------------------------------------------------------------
# Negative controls (expected errors that DO work -- probed by T2).
# --------------------------------------------------------------------------


def test_bad_order_raises():
    with pytest.raises(RuntimeError, match=r"OrderError: \[Scan invalid order: 'nonsense'\]"):
        Grid({"grid": [1, 1], "order": "nonsense"})


def test_latitudes_without_longitudes_raises():
    with pytest.raises(
        RuntimeError,
        match=r"SpecError: \[Grid: both 'latitudes' and 'longitudes' are required\]",
    ):
        Grid({"latitudes": [0.0, 1.0]})


@pytest.mark.parametrize("spec", [{"x": [0, 10, 1], "y": [0, 10, 1]}, {"x": [10, 0, -1], "y": [0, 10, 1]}])
def test_cartesian_xy_to_latlons_raises(spec):
    # documented: without a projection the points are plain Cartesian, no geographic meaning
    with pytest.raises(RuntimeError, match=r"GridError: \[Grid::to_latlons: grid points are not \(lon, lat\)\]"):
        Grid(spec).to_latlons()
