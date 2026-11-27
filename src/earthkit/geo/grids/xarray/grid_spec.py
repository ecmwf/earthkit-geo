# (C) Copyright 2024 Anemoi contributors.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


"""Functions building gridspecs from distinct coordinate values.

See :ref:`gridspec <grid-spec>` for the gridspec format.
"""

import warnings
from fractions import Fraction
from functools import lru_cache
from typing import Any, Optional

import numpy as np


def _offset(value: float, increment: float, atol: float = 1e-9) -> float:
    """Get ``value`` modulo ``increment``, snapping values near 0 or ``increment`` to 0."""
    offset = value % increment
    if np.isclose(offset, 0, atol=atol) or np.isclose(offset, increment, atol=atol):
        return 0.0
    return offset


def _is_single_precision(values: np.ndarray) -> bool:
    """Whether ``values`` is stored with float32 or coarser precision."""
    return np.issubdtype(values.dtype, np.floating) and np.finfo(values.dtype).eps >= np.finfo(np.float32).eps


def _round6(value: float) -> float:
    """Re-read ``value`` at 6 significant digits, the reliable precision of a float32."""
    return float(f"{value:.6g}")


def _round12(value: float) -> float:
    """Re-read ``value`` at 12 decimals, the precision a gridspec carries its increments at."""
    return float(f"{value:.12f}")


# how far a point of a float32 axis may sit from the lattice its cleaned increment and ends
# describe, in ulps of the coordinate magnitude. Axes built by arange, by linspace or in
# float32 stay within 1.4 ulps of it, over grids of up to 3600 points; an axis that is only
# regular on average, its points wandering from side to side, is off by orders of magnitude
# more, and is not a grid
_LATTICE_ULPS = 4.0

# the cleaned increment and ends are only used when the area they describe spans exactly the
# points of the axis: eckit counts the increments that fit in the area, so an area a fraction
# of an increment too short loses a line of points
_SPAN_ATOL = 1e-9

# an increment with more than 6 significant digits, such as the 1/12 degree of a 5-minute
# ocean grid, is beyond what rounding a float32 axis recovers: it is recovered by snapping
# the mean step to the fraction k/n it sits on, denominators up to this bound
_SNAP_MAX_DENOM = 240

# the mean step of a float32 axis is biased by the quantisation of its two ends, spread over
# all its steps; a fraction further than this relative distance from the mean step was not
# the increment the axis was built on, and is not offered to the axis check
_SNAP_RTOL = 1e-6


def _max_drift(values: np.ndarray, first: float, step: float) -> float:
    """How far the points in ``values`` are from the regular axis starting at ``first`` with increment ``step``."""
    lattice = first + step * np.arange(len(values))
    return float(np.max(np.abs(np.asarray(values, dtype=np.float64) - lattice)))


def _axis_atol(values: np.ndarray) -> float:
    """How far the points of a single-precision axis may sit from the lattice they belong to.

    A coordinate stored in float32 is quantised to an ulp of its magnitude, about 1e-5 degrees
    at the size of the Earth and much coarser than the absolute tolerance ``np.allclose`` uses.
    For magnitudes below one degree the ulp of a degree is used instead, a floor that no real
    coordinate reaches.
    """
    return _LATTICE_ULPS * np.finfo(values.dtype).eps * max(float(np.max(np.abs(values))), 1.0)


def _reproduces_axis(values: np.ndarray, step: float, low: float, high: float, span_atol: float = _SPAN_ATOL) -> bool:
    """Whether the increment ``step`` and the ends ``low``/``high`` describe ``values``.

    The ends must span exactly ``len(values) - 1`` increments, so that the grid built from
    them has as many points as the axis, and the axis rebuilt from the end the increment
    starts at must match the points. ``span_atol`` is how far the span may sit from those
    increments; a snapped axis is re-read at 12 decimals, which loosens it, see
    :func:`_snapped_axis`.
    """
    if abs((high - low) - (len(values) - 1) * abs(step)) > span_atol:
        return False
    return _max_drift(values, low if step > 0 else high, step) <= _axis_atol(values)


def _within_bounds(values: np.ndarray, bounds: tuple) -> bool:
    """Whether a single-precision coordinate stays within ``bounds``, allowing for its precision.

    Integer and float64 coordinates are not checked: their values are taken as they are.
    """
    if not _is_single_precision(values):
        return True
    low, high = bounds
    atol = _axis_atol(values)
    return float(np.min(values)) >= low - atol and float(np.max(values)) <= high + atol


def _snapped_axis(values: np.ndarray, step64: float) -> Optional[tuple]:
    """Fit a regular axis whose increment is a simple fraction of a degree.

    Returns ``(step, low, high)``, or ``None`` when no fraction reproduces the axis.

    An increment with more than 6 significant digits, such as the 1/12 degree of a
    5-minute ocean grid, is lost by the 6-digit cleaning of :func:`_regular_axis`: the
    cleaning rounds it to 0.083333, whose multiples are so far from the true lattice that
    they miss the far end of the axis by hundreds of ulps. It is recovered here instead by
    snapping the mean step to the closest fraction k/n with a small denominator, taking the
    near end of the axis (its first point, or its last when descending) snapped onto that
    fraction's lattice as the low end and the point exactly ``len(values) - 1`` increments
    beyond it as the high end. Step and ends are re-read at the 12 decimals a gridspec
    carries its increments at, which is what eckit round-trips them through. As for the
    cleaned axis, the candidate is only kept when :func:`_reproduces_axis` finds it
    describes the axis: an axis that merely averages out near a simple fraction while its
    points wander off its lattice is rejected. The span check is loosened by an ulp of the
    12-decimal rounding per increment, the price of re-reading step and ends.
    """
    if not np.isfinite(step64) or step64 == 0.0:
        return None

    magnitude = abs(step64)
    candidates = set()
    for denom in range(2, _SNAP_MAX_DENOM + 1):
        # k/n is only worth reducing when k is the multiple of 1/n the mean step sits on
        fraction = Fraction(round(magnitude * denom), denom)
        if fraction.numerator and abs(float(fraction) - magnitude) <= magnitude * _SNAP_RTOL:
            candidates.add(fraction)

    start = float(values[0]) if step64 > 0 else float(values[-1])
    for fraction in sorted(candidates, key=lambda f: abs(float(f) - magnitude)):
        unit = float(fraction)
        lattice = round(start / unit) * unit
        step = _round12(unit) if step64 > 0 else -_round12(unit)
        low = _round12(lattice)
        high = _round12(lattice + (len(values) - 1) * unit)
        span_atol = _SPAN_ATOL + len(values) * 1e-12
        if _reproduces_axis(values, step, low, high, span_atol=span_atol):
            return step, low, high

    return None


def _regular_axis(values: np.ndarray) -> Optional[tuple]:
    """Fit a regular axis to the coordinate ``values``.

    Returns ``(step, low, high)``: the signed increment and the two ends of the axis, or
    ``None`` when ``values`` does not describe a regular axis.

    A coordinate stored in float32 (e.g. read from a netCDF file) is quantised to an ulp of its
    magnitude, about 1e-5 degrees at the size of the Earth, so its steps are uniform only to
    within a few of those ulps and it does not meet the bar a float64 axis is held to: np.allclose
    with its default tolerances. Its increment and its ends carry about 7 significant decimal
    digits, so they are re-read at 6: this drops the last digit of quantisation noise and recovers
    the nominal values the file was meant to hold (0.2, 89.8, ...). Without it, an increment such
    as ``float32(0.2)`` = 0.20000000298... leaves the ends off the increment lattice by a few ulps,
    which makes the gridspec round-trip through eckit into a different, shifted grid. Those cleaned
    values are only used when they still describe the axis, which :func:`_reproduces_axis` checks:
    the same number of points, and every point of the axis within a few ulps of the lattice they
    build. An axis that is only regular on average, one whose points wander from side to side
    along it, fails that check, and so does one whose ends really are off the increment lattice
    (10.03125, say).

    A nominal value with more than 6 significant digits, the increment of a simple fraction such
    as the 1/12 degree of a 5-minute grid, survives neither that rounding nor the strict bar;
    :func:`_snapped_axis` recovers it from the lattice of the closest small fraction instead, and
    is held to the same check.

    Both of those keep the raw increment and ends, but only when the steps are uniform to
    np.allclose's default tolerances: that is the bar this fitting was held to before single
    precision was taken into account, so an axis that used to be detected still is, with the very
    same gridspec. An axis that neither cleans up nor snaps onto a simple fraction nor meets it
    is not a regular axis.
    """
    steps = np.diff(values)
    step = float(np.mean(steps))
    low = float(np.min(values))
    high = float(np.max(values))

    # the bar this fitting was held to before single precision: np.allclose's defaults
    strict = bool(np.allclose(steps, step))

    if _is_single_precision(values):
        # the increment counts the points of the grid, so it is averaged in float64: a float32
        # step carries too few digits for that average to be reliable over a whole axis
        step64 = float(np.mean(np.diff(values.astype(np.float64))))
        cleaned = (_round6(step64), _round6(low), _round6(high))
        if _reproduces_axis(values, *cleaned):
            return cleaned
        snapped = _snapped_axis(values, step64)
        if snapped is not None:
            return snapped

    if not strict:
        return None

    return step, low, high


def regular_ll_grid_spec(lats: Any, lons: Any) -> Optional[dict]:
    """Build the gridspec for a regular latitude-longitude grid.

    Parameters
    ----------
    lats : Any
        The distinct latitudes.
    lons : Any
        The distinct longitudes.

    Returns
    -------
    dict or None
        The gridspec as a dict, or ``None`` if ``lats``/``lons`` do not
        describe a regular latitude-longitude grid.
    """
    lats = np.asarray(lats)
    lons = np.asarray(lons)

    if len(lats) < 2 or len(lons) < 2:
        return None

    # the increment of each axis is the average of all its steps, which smooths out the
    # per-step floating-point noise a single step (e.g. the first) can carry; see
    # _regular_axis for how the float32 quantisation noise is dealt with
    lat_fit = _regular_axis(lats)
    lon_fit = _regular_axis(lons)

    if lat_fit is None or lon_fit is None:
        return None

    # a float32 latitude axis reaching past the poles is not one: the gridspec that fits it
    # would be clamped by eckit into a grid with fewer points than the axis has. Float64
    # axes are taken as they are, as they always have been
    if not _within_bounds(lats, (-90.0, 90.0)):
        return None

    dlat_step, south, north = lat_fit
    dlon_step, west, east = lon_fit
    dlat = abs(dlat_step)
    dlon = abs(dlon_step)

    spec = {
        "grid": [dlon, dlat],
        "area": [north, west, south, east],
    }

    if dlat_step > 0:
        spec["order"] = "i+j+"

    # the offset of the corners from the default (0, 0) reference is only snapped to 0
    # when it really is 0: snapping a small but genuine residual, a few ulps of the
    # coordinate, would move the lattice off the corner, which is enough for eckit to
    # drop a line of points. The cleaned corners of a float32 axis are on the lattice to
    # within float64 rounding, so they still end up with no "reference"
    lon_offset = _offset(west, dlon)
    lat_offset = _offset(north, dlat)
    if lon_offset or lat_offset:
        spec["reference"] = [lon_offset, lat_offset]

    return spec


def rotation_projection(mapping: Any) -> Optional[dict]:
    """Build the eckit-geo rotation projection of a CF rotated latitude-longitude grid.

    Parameters
    ----------
    mapping : Any
        The CF grid mapping attributes of the variable, e.g. the attributes
        of its ``grid_mapping`` variable.

    Returns
    -------
    dict or None
        The ``projection`` entry for a gridspec, or ``None`` if ``mapping``
        is not a usable ``rotated_latitude_longitude`` description.

    Notes
    -----
    eckit describes the rotation by its south pole. In this eckit the
    ``south_pole`` pair is ``[longitude, latitude]`` (verified numerically
    against the points of a rotated grid), the antipode of the CF north
    pole: ``grid_north_pole_longitude + 180``, ``-grid_north_pole_latitude``.
    A non-zero ``north_pole_grid_longitude`` (rotating the grid about its
    own pole, an angle eckit's rotation projection cannot carry) is
    rejected.
    """
    if not isinstance(mapping, dict):
        return None

    if mapping.get("grid_mapping_name") != "rotated_latitude_longitude":
        return None

    try:
        pole_latitude = float(mapping["grid_north_pole_latitude"])
        pole_longitude = float(mapping["grid_north_pole_longitude"])
        # absent means 0.0, the CF default
        angle = float(mapping.get("north_pole_grid_longitude", 0.0) or 0.0)
    except (KeyError, TypeError, ValueError):
        return None

    if not np.isfinite(angle) or not np.isclose(angle, 0.0):
        return None

    south_pole_longitude = (pole_longitude + 180.0) % 360.0
    return {"type": "rotation", "south_pole": [south_pole_longitude, -pole_latitude]}


def _proj_string(value: Any) -> Optional[str]:
    """The PROJ definition in ``value``, a CF attribute holding one, or ``None`` if it holds none."""
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, str) and value.strip():
        return value.strip()
    return None


def _proj_string_from_wkt(value: Any) -> Optional[str]:
    """The PROJ definition of a CF ``crs_wkt`` attribute, read with pyproj, or ``None``.

    Anything pyproj makes of it that is not a usable definition (an unreadable WKT, a CRS
    whose PROJ translation PROJ itself refuses, ...) yields ``None``, as detection has
    nothing to check the grid's points against then.

    pyproj warns that a PROJ string cannot carry all a WKT says (its datum, mainly). That
    loss is what the point check of
    :meth:`~.grid.UnstructuredGrid._projected_xy_grid_spec` is there to catch, so the
    general warning is silenced here; the grid keeps its unstructured treatment if the
    translated definition does not reproduce the dataset's own coordinates.
    """
    wkt = _proj_string(value)
    if wkt is None:
        return None

    try:
        from pyproj import CRS

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return _proj_string(CRS.from_wkt(wkt).to_proj4())
    except Exception:
        return None


def proj_projection(mapping: Any) -> Optional[dict]:
    """Build the eckit-geo PROJ projection of a CF grid mapping.

    Parameters
    ----------
    mapping : Any
        The CF grid mapping attributes of the variable, e.g. the attributes
        of its ``grid_mapping`` variable.

    Returns
    -------
    dict or None
        The ``projection`` entry for a gridspec, or ``None`` if ``mapping``
        carries no usable PROJ description of the projection.

    Notes
    -----
    The definition is read from the CF ``proj4_params`` attribute and, when
    that is absent, from ``crs_wkt`` translated with pyproj. Nothing is
    assembled from the individual projection attributes (a centre latitude,
    an ellipsoid, ...): a PROJ string that merely looks right would describe a
    grid no point of the file confirms, and no gridspec at all, which is what
    detection falls back to, is the safer answer. A
    ``rotated_latitude_longitude`` mapping is left to :func:`rotation_projection`.
    """
    if not isinstance(mapping, dict):
        return None

    if mapping.get("grid_mapping_name") == "rotated_latitude_longitude":
        return None

    source = _proj_string(mapping.get("proj4_params"))
    if source is None:
        source = _proj_string_from_wkt(mapping.get("crs_wkt"))
    if source is None:
        return None

    return {"type": "proj", "source": source}


def _xy_range(step: float, low: float, high: float) -> list:
    """The ``[start, stop, step]`` range of a fitted axis, in the direction it runs."""
    return [low, high, step] if step > 0 else [high, low, step]


def regular_xy_grid_spec(x: Any, y: Any, projection: Any) -> Optional[dict]:
    """Build the gridspec for a regular grid in the plane of a projection.

    Parameters
    ----------
    x : Any
        The distinct x coordinates of the grid, in the units of the projection plane.
    y : Any
        The distinct y coordinates of the grid, in the same units.
    projection : Any
        The ``projection`` entry the grid carries (see :func:`proj_projection`).

    Returns
    -------
    dict or None
        The gridspec as a dict, or ``None`` if ``x``/``y`` do not describe two
        regular axes.

    Notes
    -----
    The axes are fitted by :func:`_regular_axis`, the same fitting a regular
    lat/lon grid uses, so a single precision axis is cleaned the same way. Its
    tolerance is a few ulps of the coordinate, which for coordinates in metres
    is a couple of metres rather than the hundred-thousandths of a degree it
    means for latitudes and longitudes.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    if len(x) < 2 or len(y) < 2:
        return None

    x_fit = _regular_axis(x)
    y_fit = _regular_axis(y)

    if x_fit is None or y_fit is None:
        return None

    return {
        "type": "regular_xy",
        "x": _xy_range(*x_fit),
        "y": _xy_range(*y_fit),
        "projection": projection,
    }


def _unique_in_order(values: np.ndarray) -> np.ndarray:
    """Get the unique values in ``values``, in their order of first appearance."""
    _, index = np.unique(values, return_index=True)
    return values[np.sort(index)]


def regular_ll_grid_spec_from_points(lats: Any, lons: Any) -> Optional[dict]:
    """Detect a regular latitude-longitude grid from per-point coordinates.

    Unlike :func:`regular_ll_grid_spec`, ``lats`` and ``lons`` here hold one
    value per grid point (as for an unstructured grid), not the distinct
    latitudes/longitudes. They may be 1D or 2D (e.g. ``(values,)`` or
    ``(y, x)``), as long as they have the same shape.

    Parameters
    ----------
    lats : Any
        The latitude of every grid point.
    lons : Any
        The longitude of every grid point.

    Returns
    -------
    dict or None
        The gridspec as a dict, or ``None`` if ``lats``/``lons`` do not
        describe a regular latitude-longitude grid.
    """
    lats = np.asarray(lats).reshape(-1)
    lons = np.asarray(lons).reshape(-1)

    if lats.shape != lons.shape or lats.size < 4:
        return None

    lat_distinct = _unique_in_order(lats)
    lon_distinct = _unique_in_order(lons)

    if len(lat_distinct) < 2 or len(lon_distinct) < 2:
        return None

    if len(lat_distinct) * len(lon_distinct) != lats.size:
        return None

    # every (lat, lon) combination must appear exactly once for the points to
    # form a complete regular latitude-longitude grid
    num_distinct_points = len(set(zip(lats.tolist(), lons.tolist())))
    if num_distinct_points != lats.size:
        return None

    return regular_ll_grid_spec(lat_distinct, lon_distinct)


# For each supported Nside, the latitudes/longitudes of the first and last 10
# points of the "ring"- and "nested"-ordered grids, as (first_lat, first_lon,
# last_lat, last_lon). Used as a cheap fingerprint to identify a HEALPix grid
# without needing to build the full eckit.geo.Grid and its ~12 * Nside**2
# points, which gets expensive for the larger resolutions.
HEALPIX_REFERENCE_POINTS = {
    2: {
        "ring": (
            [
                66.4435356909,
                66.4435356909,
                66.4435356909,
                66.4435356909,
                41.8103148958,
                41.8103148958,
                41.8103148958,
                41.8103148958,
                41.8103148958,
                41.8103148958,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -41.8103148958,
                -41.8103148958,
                -41.8103148958,
                -41.8103148958,
                -41.8103148958,
                -41.8103148958,
                -66.4435356909,
                -66.4435356909,
                -66.4435356909,
                -66.4435356909,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                19.4712206345,
                41.8103148958,
                41.8103148958,
                66.4435356909,
                19.4712206345,
                41.8103148958,
                41.8103148958,
                66.4435356909,
                19.4712206345,
                41.8103148958,
            ],  # noqa: E501
            [45.0, 67.5, 22.5, 45.0, 135.0, 157.5, 112.5, 135.0, 225.0, 247.5],
            [
                -41.8103148958,
                -19.4712206345,
                -66.4435356909,
                -41.8103148958,
                -41.8103148958,
                -19.4712206345,
                -66.4435356909,
                -41.8103148958,
                -41.8103148958,
                -19.4712206345,
            ],  # noqa: E501
            [112.5, 135.0, 225.0, 247.5, 202.5, 225.0, 315.0, 337.5, 292.5, 315.0],
        ),
    },
    4: {
        "ring": (
            [
                78.2841476051,
                78.2841476051,
                78.2841476051,
                78.2841476051,
                66.4435356909,
                66.4435356909,
                66.4435356909,
                66.4435356909,
                66.4435356909,
                66.4435356909,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -66.4435356909,
                -66.4435356909,
                -66.4435356909,
                -66.4435356909,
                -66.4435356909,
                -66.4435356909,
                -78.2841476051,
                -78.2841476051,
                -78.2841476051,
                -78.2841476051,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                9.5940682269,
                19.4712206345,
                19.4712206345,
                30.0,
                30.0,
                41.8103148958,
                41.8103148958,
                54.3409123039,
                30.0,
                41.8103148958,
            ],  # noqa: E501
            [45.0, 56.25, 33.75, 45.0, 67.5, 78.75, 56.25, 75.0, 22.5, 33.75],
            [
                -41.8103148958,
                -30.0,
                -54.3409123039,
                -41.8103148958,
                -41.8103148958,
                -30.0,
                -30.0,
                -19.4712206345,
                -19.4712206345,
                -9.5940682269,
            ],  # noqa: E501
            [326.25, 337.5, 285.0, 303.75, 281.25, 292.5, 315.0, 326.25, 303.75, 315.0],
        ),
    },
    8: {
        "ring": (
            [
                84.1497329363,
                84.1497329363,
                84.1497329363,
                84.1497329363,
                78.2841476051,
                78.2841476051,
                78.2841476051,
                78.2841476051,
                78.2841476051,
                78.2841476051,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -78.2841476051,
                -78.2841476051,
                -78.2841476051,
                -78.2841476051,
                -78.2841476051,
                -78.2841476051,
                -84.1497329363,
                -84.1497329363,
                -84.1497329363,
                -84.1497329363,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                4.7801918472,
                9.5940682269,
                9.5940682269,
                14.4775121859,
                14.4775121859,
                19.4712206345,
                19.4712206345,
                24.6243183522,
                14.4775121859,
                19.4712206345,
            ],  # noqa: E501
            [45.0, 50.625, 39.375, 45.0, 56.25, 61.875, 50.625, 56.25, 33.75, 39.375],
            [
                -19.4712206345,
                -14.4775121859,
                -24.6243183522,
                -19.4712206345,
                -19.4712206345,
                -14.4775121859,
                -14.4775121859,
                -9.5940682269,
                -9.5940682269,
                -4.7801918472,
            ],  # noqa: E501
            [320.625, 326.25, 303.75, 309.375, 298.125, 303.75, 315.0, 320.625, 309.375, 315.0],
        ),
    },
    16: {
        "ring": (
            [
                87.0758196429,
                87.0758196429,
                87.0758196429,
                87.0758196429,
                84.1497329363,
                84.1497329363,
                84.1497329363,
                84.1497329363,
                84.1497329363,
                84.1497329363,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -84.1497329363,
                -84.1497329363,
                -84.1497329363,
                -84.1497329363,
                -84.1497329363,
                -84.1497329363,
                -87.0758196429,
                -87.0758196429,
                -87.0758196429,
                -87.0758196429,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                2.3880154633,
                4.7801918472,
                4.7801918472,
                7.1807557815,
                7.1807557815,
                9.5940682269,
                9.5940682269,
                12.0246991806,
                7.1807557815,
                9.5940682269,
            ],  # noqa: E501
            [45.0, 47.8125, 42.1875, 45.0, 50.625, 53.4375, 47.8125, 50.625, 39.375, 42.1875],
            [
                -9.5940682269,
                -7.1807557815,
                -12.0246991806,
                -9.5940682269,
                -9.5940682269,
                -7.1807557815,
                -7.1807557815,
                -4.7801918472,
                -4.7801918472,
                -2.3880154633,
            ],  # noqa: E501
            [317.8125, 320.625, 309.375, 312.1875, 306.5625, 309.375, 315.0, 317.8125, 312.1875, 315.0],
        ),
    },
    32: {
        "ring": (
            [
                88.5380288374,
                88.5380288374,
                88.5380288374,
                88.5380288374,
                87.0758196429,
                87.0758196429,
                87.0758196429,
                87.0758196429,
                87.0758196429,
                87.0758196429,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -87.0758196429,
                -87.0758196429,
                -87.0758196429,
                -87.0758196429,
                -87.0758196429,
                -87.0758196429,
                -88.5380288374,
                -88.5380288374,
                -88.5380288374,
                -88.5380288374,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                1.1937484371,
                2.3880154633,
                2.3880154633,
                3.5833216985,
                3.5833216985,
                4.7801918472,
                4.7801918472,
                5.9791567963,
                3.5833216985,
                4.7801918472,
            ],  # noqa: E501
            [45.0, 46.40625, 43.59375, 45.0, 47.8125, 49.21875, 46.40625, 47.8125, 42.1875, 43.59375],
            [
                -4.7801918472,
                -3.5833216985,
                -5.9791567963,
                -4.7801918472,
                -4.7801918472,
                -3.5833216985,
                -3.5833216985,
                -2.3880154633,
                -2.3880154633,
                -1.1937484371,
            ],  # noqa: E501
            [316.40625, 317.8125, 312.1875, 313.59375, 310.78125, 312.1875, 315.0, 316.40625, 313.59375, 315.0],
        ),
    },
    64: {
        "ring": (
            [
                89.2690292916,
                89.2690292916,
                89.2690292916,
                89.2690292916,
                88.5380288374,
                88.5380288374,
                88.5380288374,
                88.5380288374,
                88.5380288374,
                88.5380288374,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -88.5380288374,
                -88.5380288374,
                -88.5380288374,
                -88.5380288374,
                -88.5380288374,
                -88.5380288374,
                -89.2690292916,
                -89.2690292916,
                -89.2690292916,
                -89.2690292916,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                0.5968418305,
                1.1937484371,
                1.1937484371,
                1.7907846593,
                1.7907846593,
                2.3880154633,
                2.3880154633,
                2.9855060057,
                1.7907846593,
                2.3880154633,
            ],  # noqa: E501
            [45.0, 45.703125, 44.296875, 45.0, 46.40625, 47.109375, 45.703125, 46.40625, 43.59375, 44.296875],
            [
                -2.3880154633,
                -1.7907846593,
                -2.9855060057,
                -2.3880154633,
                -2.3880154633,
                -1.7907846593,
                -1.7907846593,
                -1.1937484371,
                -1.1937484371,
                -0.5968418305,
            ],  # noqa: E501
            [315.703125, 316.40625, 313.59375, 314.296875, 312.890625, 313.59375, 315.0, 315.703125, 314.296875, 315.0],  # noqa: E501
        ),
    },
    128: {
        "ring": (
            [
                89.6345165048,
                89.6345165048,
                89.6345165048,
                89.6345165048,
                89.2690292916,
                89.2690292916,
                89.2690292916,
                89.2690292916,
                89.2690292916,
                89.2690292916,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -89.2690292916,
                -89.2690292916,
                -89.2690292916,
                -89.2690292916,
                -89.2690292916,
                -89.2690292916,
                -89.6345165048,
                -89.6345165048,
                -89.6345165048,
                -89.6345165048,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                0.2984168675,
                0.5968418305,
                0.5968418305,
                0.8952829866,
                0.8952829866,
                1.1937484371,
                1.1937484371,
                1.4922462896,
                0.8952829866,
                1.1937484371,
            ],  # noqa: E501
            [45.0, 45.3515625, 44.6484375, 45.0, 45.703125, 46.0546875, 45.3515625, 45.703125, 44.296875, 44.6484375],
            [
                -1.1937484371,
                -0.8952829866,
                -1.4922462896,
                -1.1937484371,
                -1.1937484371,
                -0.8952829866,
                -0.8952829866,
                -0.5968418305,
                -0.5968418305,
                -0.2984168675,
            ],  # noqa: E501
            [
                315.3515625,
                315.703125,
                314.296875,
                314.6484375,
                313.9453125,
                314.296875,
                315.0,
                315.3515625,
                314.6484375,
                315.0,
            ],  # noqa: E501
        ),
    },
    256: {
        "ring": (
            [
                89.8172584848,
                89.8172584848,
                89.8172584848,
                89.8172584848,
                89.6345165048,
                89.6345165048,
                89.6345165048,
                89.6345165048,
                89.6345165048,
                89.6345165048,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -89.6345165048,
                -89.6345165048,
                -89.6345165048,
                -89.6345165048,
                -89.6345165048,
                -89.6345165048,
                -89.8172584848,
                -89.8172584848,
                -89.8172584848,
                -89.8172584848,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                0.1492079278,
                0.2984168675,
                0.2984168675,
                0.447627831,
                0.447627831,
                0.5968418305,
                0.5968418305,
                0.7460598782,
                0.447627831,
                0.5968418305,
            ],  # noqa: E501
            [
                45.0,
                45.17578125,
                44.82421875,
                45.0,
                45.3515625,
                45.52734375,
                45.17578125,
                45.3515625,
                44.6484375,
                44.82421875,
            ],  # noqa: E501
            [
                -0.5968418305,
                -0.447627831,
                -0.7460598782,
                -0.5968418305,
                -0.5968418305,
                -0.447627831,
                -0.447627831,
                -0.2984168675,
                -0.2984168675,
                -0.1492079278,
            ],  # noqa: E501
            [
                315.17578125,
                315.3515625,
                314.6484375,
                314.82421875,
                314.47265625,
                314.6484375,
                315.0,
                315.17578125,
                314.82421875,
                315.0,
            ],  # noqa: E501
        ),
    },
    512: {
        "ring": (
            [
                89.9086292714,
                89.9086292714,
                89.9086292714,
                89.9086292714,
                89.8172584848,
                89.8172584848,
                89.8172584848,
                89.8172584848,
                89.8172584848,
                89.8172584848,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -89.8172584848,
                -89.8172584848,
                -89.8172584848,
                -89.8172584848,
                -89.8172584848,
                -89.8172584848,
                -89.9086292714,
                -89.9086292714,
                -89.9086292714,
                -89.9086292714,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                0.0746039007,
                0.1492079278,
                0.1492079278,
                0.2238122079,
                0.2238122079,
                0.2984168675,
                0.2984168675,
                0.373022033,
                0.2238122079,
                0.2984168675,
            ],  # noqa: E501
            [
                45.0,
                45.087890625,
                44.912109375,
                45.0,
                45.17578125,
                45.263671875,
                45.087890625,
                45.17578125,
                44.82421875,
                44.912109375,
            ],  # noqa: E501
            [
                -0.2984168675,
                -0.2238122079,
                -0.373022033,
                -0.2984168675,
                -0.2984168675,
                -0.2238122079,
                -0.2238122079,
                -0.1492079278,
                -0.1492079278,
                -0.0746039007,
            ],  # noqa: E501
            [
                315.087890625,
                315.17578125,
                314.82421875,
                314.912109375,
                314.736328125,
                314.82421875,
                315.0,
                315.087890625,
                314.912109375,
                315.0,
            ],  # noqa: E501
        ),
    },
    1024: {
        "ring": (
            [
                89.9543146393,
                89.9543146393,
                89.9543146393,
                89.9543146393,
                89.9086292714,
                89.9086292714,
                89.9086292714,
                89.9086292714,
                89.9086292714,
                89.9086292714,
            ],  # noqa: E501
            [45.0, 135.0, 225.0, 315.0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5],
            [
                -89.9086292714,
                -89.9086292714,
                -89.9086292714,
                -89.9086292714,
                -89.9086292714,
                -89.9086292714,
                -89.9543146393,
                -89.9543146393,
                -89.9543146393,
                -89.9543146393,
            ],  # noqa: E501
            [112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.0, 135.0, 225.0, 315.0],
        ),
        "nested": (
            [
                0.0373019424,
                0.0746039007,
                0.0746039007,
                0.1119058905,
                0.1119058905,
                0.1492079278,
                0.1492079278,
                0.1865100283,
                0.1119058905,
                0.1492079278,
            ],  # noqa: E501
            [
                45.0,
                45.0439453125,
                44.9560546875,
                45.0,
                45.087890625,
                45.1318359375,
                45.0439453125,
                45.087890625,
                44.912109375,
                44.9560546875,
            ],  # noqa: E501
            [
                -0.1492079278,
                -0.1119058905,
                -0.1865100283,
                -0.1492079278,
                -0.1492079278,
                -0.1119058905,
                -0.1119058905,
                -0.0746039007,
                -0.0746039007,
                -0.0373019424,
            ],  # noqa: E501
            [
                315.0439453125,
                315.087890625,
                314.912109375,
                314.9560546875,
                314.8681640625,
                314.912109375,
                315.0,
                315.0439453125,
                314.9560546875,
                315.0,
            ],  # noqa: E501
        ),
    },
}

# the supported HEALPix resolutions, each with 12 * Nside**2 points
HEALPIX_NSIDES = sorted(HEALPIX_REFERENCE_POINTS)


def healpix_grid_spec_from_points(lats: Any, lons: Any) -> Optional[dict]:
    """Detect a HEALPix grid from per-point coordinates.

    ``lats`` and ``lons`` hold one value per grid point, with no further
    information about the grid (e.g. ``Nside`` or point ordering). A HEALPix
    grid is always stored as a 1D array of points, so ``lats``/``lons`` must
    be 1D. The points' count is used to guess the ``Nside``, and for each
    candidate ``order`` (``"ring"`` or ``"nested"``) the first and last 10
    points are compared against :data:`HEALPIX_REFERENCE_POINTS`. This avoids
    building the full ``eckit.geo.Grid``, which for the larger resolutions
    means generating and comparing millions of points.

    Parameters
    ----------
    lats : Any
        The latitude of every grid point, as a 1D array.
    lons : Any
        The longitude of every grid point, as a 1D array.

    Returns
    -------
    dict or None
        The gridspec as a dict, or ``None`` if ``lats``/``lons`` do not
        describe a HEALPix grid.
    """
    lats = np.asarray(lats)
    lons = np.asarray(lons)

    if lats.ndim != 1 or lons.ndim != 1 or lats.shape != lons.shape or lats.size < 10:
        return None

    nside = round((lats.size / 12) ** 0.5)
    if nside not in HEALPIX_REFERENCE_POINTS or 12 * nside * nside != lats.size:
        return None

    first_lat, first_lon = lats[:10], lons[:10]
    last_lat, last_lon = lats[-10:], lons[-10:]

    for order, (ref_first_lat, ref_first_lon, ref_last_lat, ref_last_lon) in HEALPIX_REFERENCE_POINTS[nside].items():
        if (
            np.allclose(first_lat, ref_first_lat)
            and np.allclose(first_lon, ref_first_lon)
            and np.allclose(last_lat, ref_last_lat)
            and np.allclose(last_lon, ref_last_lon)
        ):
            spec = {"grid": f"H{nside}"}
            if order != "ring":
                spec["order"] = order
            return spec

    return None


# For each supported reduced Gaussian grid ("O" for octahedral, "N" for classical,
# followed by the Gaussian number N), keyed by its point count: the grid name and the
# latitudes/longitudes of its first and last 10 points. Used as a cheap fingerprint to
# identify a reduced Gaussian grid without needing to build the full eckit.geo.Grid.
REDUCED_GG_REFERENCE_POINTS = {
    1600: (
        "O16",
        [
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    1688: (
        "N16",
        [
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
            85.7605871204,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
            -85.7605871204,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    3168: (
        "O24",
        [
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    3528: (
        "N24",
        [
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
            87.1590945559,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
            -87.1590945559,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    5248: (
        "O32",
        [
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    6114: (
        "N32",
        [
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
            87.8637988392,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
            -87.8637988392,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    10944: (
        "O48",
        [
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    13280: (
        "N48",
        [
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
            88.572168514,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
            -88.572168514,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    18688: (
        "O64",
        [
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    23112: (
        "N64",
        [
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
            88.9277353523,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
            -88.9277353523,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    28480: (
        "O80",
        [
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    35718: (
        "N80",
        [
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
            89.1415194265,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
            -89.1415194265,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    40320: (
        "O96",
        [
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    50662: (
        "N96",
        [
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
            89.2842275325,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
            -89.2842275325,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    70144: (
        "O128",
        [
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    88838: (
        "N128",
        [
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
            89.4628215686,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
            -89.4628215686,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    108160: (
        "O160",
        [
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    138346: (
        "N160",
        [
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
            89.5700895506,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
            -89.5700895506,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    167200: (
        "O200",
        [
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    213988: (
        "N200",
        [
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
            89.6559642469,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
            -89.6559642469,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    271360: (
        "O256",
        [
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    348528: (
        "N256",
        [
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
            89.7311486184,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
            -89.7311486184,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    421120: (
        "O320",
        [
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    542080: (
        "N320",
        [
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
            89.7848769072,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
            -89.7848769072,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    654400: (
        "O400",
        [
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    843490: (
        "N400",
        [
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
            89.8278746459,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
            -89.8278746459,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    1067008: (
        "O512",
        [
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    1373624: (
        "N512",
        [
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
            89.8655086877,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
            -89.8655086877,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    1661440: (
        "O640",
        [
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    2140702: (
        "N640",
        [
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
            89.8923964456,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
            -89.8923964456,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    2588800: (
        "O800",
        [
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    3336946: (
        "N800",
        [
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
            89.9139104326,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
            -89.9139104326,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    4231168: (
        "O1024",
        [
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    5447118: (
        "N1024",
        [
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
            89.9327379285,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
            -89.9327379285,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    6599680: (
        "O1280",
        [
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    8505906: (
        "N1280",
        [
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
            89.9461877157,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
            -89.9461877157,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    10297600: (
        "O1600",
        [
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    13260400: (
        "N1600",
        [
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
            89.9569484911,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
            -89.9569484911,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
    16072000: (
        "O2000",
        [
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
        ],  # noqa: E501
        [0.0, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0],
        [
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
        ],  # noqa: E501
        [180.0, 198.0, 216.0, 234.0, 252.0, 270.0, 288.0, 306.0, 324.0, 342.0],
    ),
    20696844: (
        "N2000",
        [
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
            89.9655577166,
        ],  # noqa: E501
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0],
        [
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
            -89.9655577166,
        ],  # noqa: E501
        [160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0],
    ),
}


def reduced_gg_grid_spec_from_points(lats: Any, lons: Any) -> Optional[dict]:
    """Detect a reduced Gaussian grid from per-point coordinates.

    ``lats`` and ``lons`` hold one value per grid point, with no further
    information about the grid. Unlike a regular lat/lon or HEALPix grid, a
    reduced Gaussian grid cannot be stored as a 2D array (each latitude row
    has a different number of longitude points), so ``lats``/``lons`` must
    be 1D. The points' count uniquely identifies a candidate grid (octahedral
    "O" or classical "N", at one of the supported Gaussian numbers), looked
    up in :data:`REDUCED_GG_REFERENCE_POINTS`, and the first/last 10 points
    are compared against that candidate's reference points to confirm the
    match.

    Parameters
    ----------
    lats : Any
        The latitude of every grid point, as a 1D array.
    lons : Any
        The longitude of every grid point, as a 1D array.

    Returns
    -------
    dict or None
        The gridspec as a dict, or ``None`` if ``lats``/``lons`` do not
        describe a reduced Gaussian grid.
    """
    lats = np.asarray(lats)
    lons = np.asarray(lons)

    if lats.ndim != 1 or lons.ndim != 1 or lats.shape != lons.shape or lats.size < 10:
        return None

    reference = REDUCED_GG_REFERENCE_POINTS.get(lats.size)
    if reference is None:
        return None

    grid_name, ref_first_lat, ref_first_lon, ref_last_lat, ref_last_lon = reference

    if (
        np.allclose(lats[:10], ref_first_lat)
        and np.allclose(lons[:10], ref_first_lon)
        and np.allclose(lats[-10:], ref_last_lat)
        and np.allclose(lons[-10:], ref_last_lon)
    ):
        return {"grid": grid_name}

    return None


def _scalar_number(value: Any) -> Optional[float]:
    """``value``, a single finite number, as a Python float, or ``None`` if it is not one.

    A netCDF attribute is read as a numpy scalar, which :func:`float` handles. A vector
    attribute, a string or a missing value is not the single number a grid parameter
    should be, and yields ``None`` rather than an exception.
    """
    try:
        number = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError):
        return None

    if number.ndim != 0 or not bool(np.isfinite(number)):
        return None

    return float(number)


def _whole_number(value: Any) -> Optional[int]:
    """``value``, a single finite number, as a Python int when it is a whole one, or ``None``."""
    number = _scalar_number(value)
    if number is None or number != float(int(number)):
        return None
    return int(number)


def _integer_row(values: Any) -> Optional[np.ndarray]:
    """``values`` as a 1D array of whole numbers in int64, or ``None`` if it is not one.

    Only numpy walks the values: the array a point index is held in has as many
    entries as the grid has points, over 6.5 million for an O1280, which no
    Python loop over its values should be asked to step through.
    """
    try:
        array = np.asarray(values)
    except Exception:
        return None

    if array.ndim != 1 or array.size == 0 or np.issubdtype(array.dtype, np.bool_):
        return None

    if np.issubdtype(array.dtype, np.integer):
        low, high = int(array.min()), int(array.max())
    elif np.issubdtype(array.dtype, np.floating):
        # a coordinate which packing or a float dtype turned into floats is still
        # usable as long as every value is a whole number
        if not bool(np.all(np.isfinite(array))) or not bool(np.all(array == np.floor(array))):
            return None
        low, high = int(array.min()), int(array.max())
    else:
        return None

    if low < np.iinfo(np.int64).min or high > np.iinfo(np.int64).max:
        return None

    return array.astype(np.int64, copy=False)


# an octahedral reduced Gaussian grid has 20 points on the latitude line closest to each
# pole, and 4 more on every line towards the equator, its southern half mirroring the
# northern one, so the table of points per latitude says which grid it is
_OCTAHEDRAL_FIRST_ROW = 20
_OCTAHEDRAL_ROW_STEP = 4


@lru_cache(maxsize=None)
def _classical_pl(rows: int) -> Optional[tuple]:
    """The points per latitude of the classical reduced Gaussian grid ``N<rows>``.

    A classical grid does not follow a rule a table can be checked against; the
    one eckit builds ``N<rows>`` from is read off the grid itself, whose points
    come west to east and north to south, each latitude line holding its points
    together. eckit only knows the classical resolutions of the ECMWF tables, so
    a ``rows`` that is not one of them raises: the result is then ``None``, as
    for a grid that cannot be built. Building the grid means generating all its
    points, which for the larger resolutions takes about a second, so the table
    is cached per resolution, and only asked for once a candidate table has the
    point count of the grid it claims.
    """
    try:
        from eckit.geo import Grid as EckitGrid

        lat, _ = EckitGrid({"grid": f"N{rows}"}).to_latlons()
    except Exception:
        return None

    lat = np.asarray(lat, dtype=np.float64)
    ends = np.flatnonzero(np.diff(lat) != 0) + 1
    return tuple(int(size) for size in np.diff(np.concatenate(([0], ends, [lat.size]))))


def _cf_reduced_gaussian_table(
    attrs: Any,
    pl: Any,
    accumulated: bool = False,
) -> Optional[tuple]:
    """Check a CF reduced Gaussian grid mapping and its global points-per-latitude table.

    These are the checks that a mapping and a table can answer on their own, before
    any of the file's points are looked at: the mapping is a ``reduced_gaussian`` one
    naming an "octahedral" or "normal" grid whose first meridian is Greenwich, and
    the table is a 1D even number of mirrored lines of whole counts, following the
    octahedral rule when it claims to. The classical table of an ``N<n>`` grid
    follows no rule to check against, so it is only the mirroring that is settled
    here; the comparison against the table eckit builds ``N<n>`` from (a lookup that
    means generating that grid) is left to the caller, once the cheaper point counts
    and index checks have narrowed the candidates.

    Returns
    -------
    tuple or None
        ``(subtype, counts)``, the normalised grid subtype and the points per
        latitude (``pl`` un-accumulated when ``accumulated``), or ``None`` if the
        metadata does not describe the global table of a reduced Gaussian grid.
    """
    if not isinstance(attrs, dict):
        return None

    if attrs.get("grid_mapping_name") != "reduced_gaussian":
        return None

    subtype = attrs.get("grid_subtype")
    if isinstance(subtype, bytes):
        subtype = subtype.decode("utf-8", errors="replace")
    if not isinstance(subtype, str):
        return None
    subtype = subtype.strip().lower()
    if subtype not in ("octahedral", "normal"):
        return None

    # an attribute the caller holds as None is an attribute the dataset does not carry,
    # whose CF default is the Greenwich meridian eckit builds
    meridian = attrs.get("longitude_of_first_meridian", 0.0)
    meridian = _scalar_number(0.0 if meridian is None else meridian)
    if meridian is None or not np.isclose(meridian, 0.0, rtol=0.0, atol=1e-9):
        return None

    counts = _integer_row(pl)
    if counts is None or counts.size % 2:
        return None

    if accumulated:
        # the accumulated points per latitude is the cumulative sum of the points per
        # latitude, along the latitude dimension, whichever end of it it starts at
        counts = np.diff(np.concatenate(([0], counts)))

    rows = counts.size // 2
    if subtype == "octahedral":
        north = _OCTAHEDRAL_FIRST_ROW + _OCTAHEDRAL_ROW_STEP * np.arange(rows, dtype=np.int64)
        if not np.array_equal(counts, np.concatenate((north, north[::-1]))):
            return None
    elif not bool(np.array_equal(counts, counts[::-1])):
        # a classical table mirrors itself like an octahedral one, but its rows follow
        # no rule: which grid it belongs to, if any, is only known against the tables
        # eckit has, which the caller looks up once the cheap checks are exhausted
        return None

    return subtype, counts


def cf_reduced_gaussian_grid_spec(
    attrs: Any,
    pl: Any,
    n_points: Optional[int] = None,
    index_values: Any = None,
    accumulated: bool = False,
) -> Optional[dict]:
    """Build the gridspec of a reduced Gaussian grid from its CF grid mapping.

    The CF conventions describe a reduced Gaussian grid (Appendix F, "Reduced
    Gaussian") with a ``reduced_gaussian`` grid mapping whose parameters give the
    number of longitude points on every latitude line, and a
    ``reduced_gaussian_index`` coordinate numbering the points west to east and
    north to south, from 0 to the last point of the grid. A grid so described is
    named by its subtype: an "octahedral" one is the eckit grid ``O<n>``, a
    "normal", classical one the eckit grid ``N<n>``, ``n`` the number of latitude
    lines of one hemisphere, which is all the gridspec needs to carry::

        {"grid": "O<n>"}
        {"grid": "N<n>"}

    An ``O<n>`` grid has a rule for its rows, 20 on the line closest to each pole
    and 4 more on every line towards the equator, so any table following it is
    such a grid. A classical grid has no rule to follow: it is only recognised
    when its table is exactly the one eckit builds ``N<n>`` from, a table which
    merely looks Gaussian belongs to no grid eckit can name, and is not detected.
    A classical grid whose points are in the file can still be found from them,
    as :func:`reduced_gg_grid_spec_from_points` does.

    Parameters
    ----------
    attrs : Any
        The CF grid mapping attributes of the variable, e.g. the attributes of its
        ``grid_mapping`` variable. Only a ``reduced_gaussian`` mapping with an
        "octahedral" or a "normal" ``grid_subtype`` is a candidate for a gridspec.
    pl : Any
        The number of longitude points on each latitude line, along the dimension the
        ``latitude_dimension`` attribute names, as a 1D array. Pass the
        ``accumulated_points_per_latitude`` array the same way, with ``accumulated``
        set.
    n_points : Optional[int]
        The size of the dimension the points are indexed by, i.e. the number of points
        of the grid. Defaults to the size of ``index_values``.
    index_values : Any
        The ``reduced_gaussian_index`` of every point, as a 1D array.
    accumulated : bool
        Whether ``pl`` holds the accumulated points per latitude, the cumulative sum
        of the points per latitude, rather than those points per latitude themselves.

    Returns
    -------
    dict or None
        The gridspec as a dict, or ``None`` if the metadata does not describe a global
        reduced Gaussian grid, octahedral or classical.

    Notes
    -----
    A gridspec is only returned when every parameter is one the grid named by it
    carries: the table of points per latitude has an even number of lines and its
    southern half mirrors the northern one, following the rule of an ``O<n>`` grid
    or being exactly the table of the ``N<n>`` grid of the same number of lines,
    and summing up to ``n_points``; ``longitude_of_first_meridian`` is absent or 0,
    as eckit starts both grids at the Greenwich meridian; and the index is exactly
    0 to ``n_points`` - 1, which the index of a regional or land-only subset is
    not, as it keeps the numbers the points have on the global grid. Anything else
    yields ``None``, never an exception for data which turns out to be unusable,
    and never a gridspec which does not describe the points the file holds.

    The direction the table of points per latitude runs in is not checked, and does
    not matter: a Gaussian table reads the same from either end, and it is the
    index, not the table, which says in which order the points come.
    """
    table = _cf_reduced_gaussian_table(attrs, pl, accumulated)
    if table is None:
        return None
    subtype, counts = table
    rows = counts.size // 2

    index = _integer_row(index_values)
    if index is None:
        return None

    if n_points is None:
        n_points = int(index.size)
    else:
        n_points = _whole_number(n_points)
        if n_points is None:
            return None

    # the point count is checked before the index is built, so that a count which is not
    # one cannot send numpy looking for an array of points no file has; it is also checked
    # before any eckit grid is built, a table which does not even hold the points of the
    # grid it claims is no grid, and should not cost a grid build
    if n_points < 1 or index.size != n_points or int(counts.sum(dtype=np.int64)) != n_points:
        return None

    if not np.array_equal(index, np.arange(n_points, dtype=np.int64)):
        return None

    if subtype == "octahedral":
        return {"grid": f"O{rows}"}

    # the classical table of N<rows> says whether this is that grid: a table which is
    # well shaped and well counted but belongs to no resolution eckit knows describes
    # no grid a name can carry, and is left undetected
    reference = _classical_pl(rows)
    if reference is None or not np.array_equal(counts, reference):
        return None

    return {"grid": f"N{rows}"}
