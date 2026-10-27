# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import numbers

import numpy as np

from .coord import latlon_to_xyz, xyz_to_latlon
from .figure import IFS_SPHERE

METHODS = ("great_circle", "rhumb", "linear")

# Tolerance (in radians) below which two points are considered to be coincident
# (great circle) or at the same latitude (rhumb).
_EPSILON = 1e-12


def _wrap_lon(lon, minimum):
    """Wrap longitudes (degrees) into the range [minimum, minimum + 360)."""
    upper = minimum + 360.0
    wrapped = minimum + np.mod(lon - minimum, 360.0)
    # np.mod can return exactly the modulus for tiny negative values
    return np.where(wrapped >= upper, minimum, wrapped)


def _mercator(phi):
    return np.log(np.tan(np.pi / 4 + phi / 2))


def _great_circle_prepare(lat, lon):
    """Return per-segment (angle, start vector, tangent vector) for the great circle path.

    The tangent is the unit vector at the start pointing towards the end, so that the
    points of the segment are ``cos(t * angle) * start + sin(t * angle) * tangent``.
    This is numerically stable for any angle, unlike the usual slerp formula which
    divides by ``sin(angle)``.
    """
    x, y, z = latlon_to_xyz(lat, lon)
    vec = np.stack([x, y, z], axis=-1)
    a = vec[:-1]
    b = vec[1:]
    cross_vec = np.cross(a, b)
    cross = np.linalg.norm(cross_vec, axis=-1)
    dot = np.sum(a * b, axis=-1)

    coincident = cross < _EPSILON
    if np.any(coincident & (dot < 0)):
        i = int(np.flatnonzero(coincident & (dot < 0))[0])
        raise ValueError(
            f"polyline: segment {i} of the great circle polyline has antipodal ends "
            f"({lat[i]}, {lon[i]}) and ({lat[i + 1]}, {lon[i + 1]}); the path is not unique. "
            "Add an intermediate waypoint"
        )

    angle = np.where(coincident, 0.0, np.arctan2(cross, dot))
    normal = cross_vec / np.where(coincident, 1.0, cross)[:, None]
    tangent = np.where(coincident[:, None], 0.0, np.cross(normal, a))
    return angle, a, tangent


def _rhumb_prepare(lat, lon):
    """Return per-segment lengths and the quantities needed to interpolate the rhumb line."""
    phi = np.deg2rad(lat)
    dphi = np.diff(phi)
    dlam = np.deg2rad(np.diff(lon))
    psi = _mercator(phi)
    dpsi = np.diff(psi)

    flat = np.abs(dphi) < _EPSILON
    safe_dpsi = np.where(flat, 1.0, dpsi)
    q = np.where(flat, np.cos(phi[:-1]), dphi / safe_dpsi)
    length = np.sqrt(dphi**2 + (q * dlam) ** 2)
    return length, flat, psi, dpsi


def _distribute(lengths, extra):
    """Share ``extra`` points among segments in proportion to ``lengths``.

    Largest-remainder method; ties go to the earlier segment. Zero-length
    segments get no points unless all segments have zero length, in which case
    the points are spread evenly by count.
    """
    nseg = len(lengths)
    total = float(np.sum(lengths))
    if total > 0:
        quota = extra * lengths / total
        eligible = lengths > 0
    else:
        quota = np.full(nseg, extra / nseg)
        eligible = np.ones(nseg, dtype=bool)

    # Plain largest remainder: round every quota down and give out the leftovers.
    # Float noise needs no special care: a quota of 2.9999999999999996 floors to 2
    # with a remainder of almost 1, so it gets one of the leftovers and ends at 3.
    counts = np.floor(quota).astype(int)
    left = extra - int(counts.sum())
    if left < 0:  # pragma: no cover
        raise RuntimeError(f"polyline: internal error in the sharing of the points (left={left})")

    candidates = np.flatnonzero(eligible)
    if left >= len(candidates):  # pragma: no cover
        # cannot happen without floating point drift: give whole rounds first
        counts[candidates] += left // len(candidates)
        left %= len(candidates)

    if left > 0:
        # Largest remainder first. Remainders are snapped to a grid of size ``tol``
        # (machine scale) so that ties up to rounding noise go to the earlier segment,
        # while real differences still decide the order.
        tol = 1e-12 * max(extra, 1)
        grid = np.round((quota[candidates] - counts[candidates]) / tol)
        # stable sort: the last key of lexsort is the primary one
        order = np.lexsort((candidates, -grid))
        counts[candidates[order[:left]]] += 1
    return counts


def _validate_waypoints(name, what, lats, lons):
    """Check the common input conditions of the polyline functions.

    Returns the latitudes and longitudes as float arrays and the number of vertices.
    ``name`` prefixes the error messages and ``what`` names the vertices in them.
    """
    lat = np.asarray(lats, dtype=float)
    lon = np.asarray(lons, dtype=float)

    if lat.ndim != 1 or lon.ndim != 1:
        raise ValueError(f"{name}: lats and lons must be 1D! lats.ndim={lat.ndim}, lons.ndim={lon.ndim}")

    if lat.shape != lon.shape:
        raise ValueError(f"{name}: lats and lons must have the same length! {lat.shape[0]} != {lon.shape[0]}")

    k = lat.shape[0]
    if k < 2:
        raise ValueError(f"{name}: at least 2 {what} are required! Got {k}")

    if not (np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))):
        raise ValueError(f"{name}: lats and lons must be finite")

    if np.any(np.abs(lat) > 90):
        raise ValueError(f"{name}: latitudes must be in the range [-90, 90]")

    return lat, lon, k


def polyline(lats, lons, n=100, method="great_circle"):
    """Sample points along a polyline defined by a sequence of waypoints.

    Parameters
    ----------
    lats: array-like
        Latitudes of the waypoints (degrees). Must be 1D and contain at least two
        values in the range [-90, 90].
    lons: array-like
        Longitudes of the waypoints (degrees). Must have the same shape as ``lats``.
    n: int, optional
        Total number of points in the result, including the waypoints (default: 100).
        Must be at least the number of waypoints.
    method: str, optional
        Path between consecutive waypoints. One of ``"great_circle"`` (default),
        ``"rhumb"`` or ``"linear"``.

    Returns
    -------
    ndarray
        Latitudes (degrees) of the ``n`` points.
    ndarray
        Longitudes (degrees) of the ``n`` points.

    Raises
    ------
    ValueError
        If the inputs are not 1D, have different lengths or have fewer than two
        waypoints; if a latitude is outside [-90, 90] or a value is not finite; if ``n``
        is not an int or is smaller than the number of waypoints; if ``method`` is
        unknown; if a great circle segment has antipodal ends; or if a waypoint is at a
        pole when ``method`` is ``"rhumb"``.

    Notes
    -----
    The result starts with the first waypoint and ends with the last one. Every waypoint
    is included exactly as given (it is copied, never recomputed or normalised), in the
    order: waypoint 0, the points inside segment 0, waypoint 1, and so on.

    The model is a sphere, so the result does not depend on its radius.

    The methods are:

    - ``"great_circle"``: each segment follows the shortest great circle arc between
      its ends (spherical linear interpolation). The path is not unique if the ends are
      antipodal, which raises an error. Longitude differences are not used to choose
      the direction, so a global great circle needs at least three waypoints, each less
      than 180° apart from the next one (see below).
    - ``"rhumb"``: each segment is a rhumb line (loxodrome), a path of constant bearing.
      The longitude difference is taken as given, with no wrapping and no
      shortest-way choice, e.g. going from 170° to -170° runs westwards through 0° (340°
      of longitude); use 190° as the end to cross the dateline eastwards. Waypoints at
      the poles are not allowed.
    - ``"linear"``: latitude and longitude values are interpolated linearly as given,
      with no wrapping.

    Sharing the points: the ``n`` minus the number of waypoints extra points are shared
    among the segments in proportion to their length measured along the chosen path
    (central angle for ``"great_circle"``, rhumb line length for ``"rhumb"``
    and the Euclidean distance in degrees for ``"linear"``). Shares are rounded with the
    largest remainder method, so the total is exactly ``n``; ties go to the earlier
    segment. A zero-length segment (repeated waypoint) gets no extra points. If all
    segments have zero length, the extra points are spread evenly by count. Within a
    segment the extra points are equally spaced along the path.

    Longitude range: the interior points (the ones that are not waypoints) of the
    ``"great_circle"`` and ``"rhumb"`` methods are returned in [0, 360) if all the input
    longitudes are non-negative and in [-180, 180) otherwise. The ``"linear"`` method
    does not change longitudes. Interpolation always uses the longitudes as given.

    Global lines:

    - ``"rhumb"`` and ``"linear"``: a full parallel needs the longitudes to span 360°,
      e.g. ``lons=[-180, 180]`` or ``lons=[0, 360]``, with constant latitude.
    - ``"great_circle"``: use three or more waypoints, each less than 180° apart from the
      next one, e.g. ``lats=[0, 0, 0, 0]`` and ``lons=[-180, -60, 60, 180]``.

    The first and the last point of a full circle are at the same place.

    Examples
    --------
    Points along the great circle between Reading and Bologna. The first and the last
    points are the waypoints.

    >>> from earthkit.geo.polyline import polyline
    >>> lats, lons = polyline([51.45, 44.49], [-0.97, 11.34], n=5)
    >>> lats.shape
    (5,)
    >>> float(lats[0]), float(lons[0])
    (51.45, -0.97)
    >>> float(lats[-1]), float(lons[-1])
    (44.49, 11.34)

    Sample the full equator every 90° (the last point is the same place as the first
    one) with a rhumb line:

    >>> lats, lons = polyline([0, 0], [-180, 180], n=5, method="rhumb")
    >>> lats.tolist()
    [0.0, 0.0, 0.0, 0.0, 0.0]
    >>> lons.tolist()
    [-180.0, -90.0, 0.0, 90.0, 180.0]

    The same with great circle segments needs more than two waypoints:

    >>> lats, lons = polyline([0, 0, 0, 0], [-180, -60, 60, 180], n=7)
    >>> [round(float(x), 6) for x in lats]
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    >>> [round(float(x), 6) for x in lons]
    [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0]

    """
    lat, lon, k = _validate_waypoints("polyline", "waypoints", lats, lons)

    if isinstance(n, bool) or not isinstance(n, numbers.Integral):
        raise ValueError(f"polyline: n must be an int! Got {type(n).__name__}")

    n = int(n)
    if n < k:
        raise ValueError(f"polyline: n must be at least the number of waypoints! n={n} < {k}")

    if method not in METHODS:
        raise ValueError(f"polyline: unknown method={method!r}! Must be one of {METHODS}")

    if method == "rhumb" and np.any(np.abs(lat) == 90):
        raise ValueError("polyline: waypoints at the poles are not allowed with method='rhumb'")

    # segment lengths along the chosen path
    if method == "great_circle":
        lengths, gc_a, gc_v = _great_circle_prepare(lat, lon)
    elif method == "rhumb":
        lengths, flat, psi, dpsi = _rhumb_prepare(lat, lon)
    else:
        lengths = np.hypot(np.diff(lat), np.diff(lon))

    counts = _distribute(lengths, n - k)

    out_lat = np.empty(n, dtype=float)
    out_lon = np.empty(n, dtype=float)
    interior = np.zeros(n, dtype=bool)

    pos = 0
    for i in range(k - 1):
        out_lat[pos] = lat[i]
        out_lon[pos] = lon[i]
        e = int(counts[i])
        if e > 0:
            t = np.arange(1, e + 1) / (e + 1)
            sl = slice(pos + 1, pos + 1 + e)

            if method == "great_circle":
                # rotation of the start vector towards the end vector in their plane;
                # for coincident ends the angle is 0 and all the points are at the start
                angle = t * lengths[i]
                p = np.cos(angle)[:, None] * gc_a[i] + np.sin(angle)[:, None] * gc_v[i]
                lat_i, lon_i = xyz_to_latlon(p[:, 0], p[:, 1], p[:, 2])
            elif method == "rhumb":
                lat_i = lat[i] + t * (lat[i + 1] - lat[i])
                dlon = lon[i + 1] - lon[i]
                if flat[i]:
                    frac = t
                else:
                    frac = (_mercator(np.deg2rad(lat_i)) - psi[i]) / dpsi[i]
                lon_i = lon[i] + dlon * frac
            else:
                lat_i = lat[i] + t * (lat[i + 1] - lat[i])
                lon_i = lon[i] + t * (lon[i + 1] - lon[i])

            out_lat[sl] = lat_i
            out_lon[sl] = lon_i
            interior[sl] = True
        pos += 1 + e

    out_lat[pos] = lat[-1]
    out_lon[pos] = lon[-1]

    # normalise the longitude of the interior points only, as the last step
    if method != "linear" and np.any(interior):
        minimum = 0.0 if np.all(lon >= 0) else -180.0
        out_lon[interior] = _wrap_lon(out_lon[interior], minimum)

    return out_lat, out_lon


def _local_basis(lat, lon):
    """Return the unit vectors pointing east and north at the given points (``(..., 3)`` arrays)."""
    phi = np.deg2rad(lat)
    lam = np.deg2rad(lon)
    sin_phi, cos_phi = np.sin(phi), np.cos(phi)
    sin_lam, cos_lam = np.sin(lam), np.cos(lam)
    east = np.stack([-sin_lam, cos_lam, np.zeros_like(lam)], axis=-1)
    north = np.stack([-sin_phi * cos_lam, -sin_phi * sin_lam, cos_phi], axis=-1)
    return east, north


def _great_circle_tangents(lat, lon):
    """Return the per-segment (east, north) tangents at the start and at the end of the segment.

    Result: ``(valid, start, end)``, with ``start`` and ``end`` of shape ``(k - 1, 2)``.
    Invalid segments are the zero-length ones.
    """
    x, y, z = latlon_to_xyz(lat, lon)
    vec = np.stack([x, y, z], axis=-1)
    a = vec[:-1]
    b = vec[1:]
    cross_vec = np.cross(a, b)
    cross = np.linalg.norm(cross_vec, axis=-1)
    dot = np.sum(a * b, axis=-1)

    coincident = cross < _EPSILON
    antipodal = coincident & (dot < 0)
    if np.any(antipodal):
        i = int(np.flatnonzero(antipodal)[0])
        raise ValueError(
            f"polyline_normal: segment {i} of the great circle polyline has antipodal ends "
            f"({lat[i]}, {lon[i]}) and ({lat[i + 1]}, {lon[i + 1]}); the path is not unique. "
            "Add an intermediate waypoint"
        )

    normal = cross_vec / np.where(coincident, 1.0, cross)[:, None]
    east, north = _local_basis(lat, lon)

    ends = []
    for t, sl in ((np.cross(normal, a), slice(None, -1)), (np.cross(normal, b), slice(1, None))):
        ends.append(np.stack([np.sum(t * east[sl], axis=-1), np.sum(t * north[sl], axis=-1)], axis=-1))
    return ~coincident, ends[0], ends[1]


def _mercator_pole_safe(lat):
    """Return the isometric latitude ``psi = ln(tan(pi/4 + phi/2))`` of latitudes (degrees, not at a pole).

    It is computed from the colatitude ``90 - |lat|``, which is exact in degrees, so it stays
    accurate (in the relative sense) next to a pole, where ``pi/4 + phi/2`` or ``sin(phi)``
    would round to a value that is too close to the pole.
    """
    colat = np.deg2rad(90.0 - np.abs(lat))
    return -np.sign(lat) * np.log(np.tan(colat / 2))


def _mercator_difference(lat1, lat2):
    """Return ``psi(phi2) - psi(phi1)``, accurate for close latitudes and next to the poles (degrees, not at a pole).

    There are two cases, chosen by ``x = tanh(dpsi)``:

    - ``|x| < 0.5`` (``|dpsi| < 0.55``): ``psi = atanh(sin(phi))``, so
      ``x = (s2 - s1) / (1 - s1 * s2)``, where ``s2 - s1 = 2 cos(m) sin(d / 2)`` and
      ``1 - s1 * s2 = sin(d / 2)**2 + cos(m)**2``, with ``d = phi2 - phi1`` and
      ``m = (phi1 + phi2) / 2``. No subtraction of close values is involved, so the sign of
      the result is the sign of ``d`` and it is not zero if ``d`` is not zero (``d`` is
      computed in degrees, where the difference of close values is exact). ``x`` is
      far from the saturation of ``atanh``, which would happen for a large ``|dpsi|``.
    - Otherwise the direct difference of the pole safe ``psi`` values is used. The
      absolute error of each value is about ``1e-16 * |psi|`` (``|psi| < 40``), so the
      cancellation costs less than 1e-13 of relative accuracy for ``|dpsi| >= 0.55``.
      The threshold 0.5 on ``x`` is where the two errors are both small.

    Segments with both ends at least 45° from the equator, in the same hemisphere, use the
    colatitudes instead (see the code), to stay accurate when both ends are very close to a pole.
    """
    half = np.sin(np.deg2rad(lat2 - lat1) / 2)
    cos_mid = np.cos(np.deg2rad((lat1 + lat2) / 2))
    x = 2 * cos_mid * half / (half**2 + cos_mid**2)

    small = np.abs(x) < 0.5
    stable = np.arctanh(np.where(small, x, 0.0))
    direct = _mercator_pole_safe(lat2) - _mercator_pole_safe(lat1)
    general = np.where(small, stable, direct)

    # Both ends at least 45° from the equator, in the same hemisphere: use the colatitudes
    # c = 90 - |lat|. They are exact in floating point there (Sterbenz), while the formula
    # above loses the distance to the pole through ``cos(m)``, ``sin(phi)`` or ``phi``.
    polar = ((lat1 >= 45) & (lat2 >= 45)) | ((lat1 <= -45) & (lat2 <= -45))
    c1 = np.deg2rad(np.where(polar, 90.0 - np.abs(lat1), 1.0))
    c2 = np.deg2rad(np.where(polar, 90.0 - np.abs(lat2), 2.0))
    # in the north psi = -ln(tan(c / 2)), so dpsi = ln(tan(c1 / 2) / tan(c2 / 2)); the
    # argument of log1p, ``ratio - 1`` or ``1 / ratio - 1`` for the opposite sign, is
    # never negative, which keeps it accurate for tiny and for large steps
    dc = np.deg2rad(np.where(polar, np.abs(lat2) - np.abs(lat1), 1.0))  # c1 - c2, exact in degrees
    towards_pole = dc >= 0
    big, small_c = np.where(towards_pole, c1, c2), np.where(towards_pole, c2, c1)
    arg = np.sin(np.abs(dc) / 2) / (np.cos(big / 2) * np.sin(small_c / 2))
    north = np.where(towards_pole, 1.0, -1.0) * np.log1p(arg)
    return np.where(polar, np.where(lat1 >= 0, 1.0, -1.0) * north, general)


def _rhumb_tangents(lat, lon):
    """Return ``(valid, start, end)`` as in ``_great_circle_tangents`` for rhumb line segments."""
    dphi = np.diff(lat)
    dlam = np.diff(lon)
    dpsi = _mercator_difference(lat[:-1], lat[1:])
    beta = np.arctan2(np.deg2rad(dlam), dpsi)
    tangent = np.stack([np.sin(beta), np.cos(beta)], axis=-1)
    return (dphi != 0) | (dlam != 0), tangent, tangent


def _linear_tangents(lat, lon):
    """Return ``(valid, start, end)`` as in ``_great_circle_tangents`` for linear segments."""
    dphi = np.diff(lat)
    dlam = np.diff(lon)
    cos_phi = np.cos(np.deg2rad(lat))
    valid = (dphi != 0) | (dlam != 0)

    ends = []
    for c in (cos_phi[:-1], cos_phi[1:]):
        te = dlam * c
        norm = np.hypot(te, dphi)
        norm = np.where(valid, norm, 1.0)
        ends.append(np.stack([te / norm, dphi / norm], axis=-1))
    return valid, ends[0], ends[1]


def polyline_normal(lats, lons, method="great_circle"):
    """Compute the unit normal vectors at the vertices of a polyline.

    Parameters
    ----------
    lats: array-like
        Latitudes of the vertices (degrees). Must be 1D and contain at least two
        values in the range [-90, 90]. No vertex can be at a pole.
    lons: array-like
        Longitudes of the vertices (degrees). Must have the same shape as ``lats``.
    method: str, optional
        Path between consecutive vertices. One of ``"great_circle"`` (default),
        ``"rhumb"`` or ``"linear"``. It is the same as in :py:func:`polyline`.

    Returns
    -------
    ndarray
        East component of the unit normal at each vertex.
    ndarray
        North component of the unit normal at each vertex.

    Raises
    ------
    ValueError
        If the inputs are not 1D, have different lengths or have fewer than two
        vertices; if a latitude is outside [-90, 90] or a value is not finite; if
        ``method`` is unknown; if a vertex is at a pole; if a great circle segment has
        antipodal ends; if all the segments have zero length; or if the line reverses
        its direction at a vertex (U-turn).

    Notes
    -----
    The normal is a unit vector in the local horizontal plane, given by its ``(east, north)``
    components (like wind ``u`` and ``v``), perpendicular to the line. It points to the
    **left** of the direction of travel, which goes from the first vertex to the last one.
    If the unit tangent is ``(te, tn)`` the normal is ``(-tn, te)``, so a line heading
    east has the normal ``(0, 1)`` (north).

    At the first and the last vertex the tangent is the one of the only segment there.
    At an interior vertex it is the bisector ``normalize(t_in + t_out)``, where ``t_in`` is
    the direction of travel arriving at the vertex and ``t_out`` the one leaving it. If the
    line reverses at the vertex (``|t_in + t_out|`` is nearly zero) the bisector is
    undefined and an error is raised.

    The methods are the same as in :py:func:`polyline`, and the normals follow the path
    that it samples:

    - ``"great_circle"``: the normal at both ends of a segment lies in the plane
      tangent to the sphere and is perpendicular to the plane of the great circle arc.
      Antipodal ends are an error.
    - ``"rhumb"``: the bearing is constant along the segment, ``atan2(dlon, dpsi)``, where
      ``dpsi`` is the difference of the isometric (Mercator) latitudes. The longitude
      difference is taken as given, with no wrapping, so ``lons=[0, 360]`` is a full
      parallel heading east.
    - ``"linear"``: latitude and longitude change linearly in degrees (longitude difference
      as given). The tangent at a vertex with latitude ``lat`` is the normalised
      ``(dlon * cos(lat), dlat)``, evaluated at each end of the segment.

    Zero-length segments (repeated vertices) are skipped: a vertex uses the nearest segment
    with a non-zero length before it for ``t_in`` and the nearest one after it for ``t_out``,
    so repeated vertices get the same normal. For ``"great_circle"`` two vertices are
    the same if they are at the same place, so ``lons=[0, 360]`` is a zero-length segment,
    while for ``"rhumb"`` and ``"linear"`` they are the same only if both latitudes
    and longitudes are equal. An error is raised if all the segments have zero length.

    Longitudes are used as given and the output is not normalised in any way.

    Examples
    --------
    The normal of a line heading east along the equator points north:

    >>> from earthkit.geo.polyline import polyline_normal
    >>> east, north = polyline_normal([0, 0], [0, 10])
    >>> [round(float(x), 6) + 0.0 for x in east]
    [0.0, 0.0]
    >>> [round(float(x), 6) + 0.0 for x in north]
    [1.0, 1.0]

    The normals of a line heading north along a meridian point west:

    >>> east, north = polyline_normal([0, 10, 20], [0, 0, 0], method="rhumb")
    >>> [round(float(x), 6) + 0.0 for x in east]
    [-1.0, -1.0, -1.0]
    >>> [round(float(x), 6) + 0.0 for x in north]
    [0.0, 0.0, 0.0]

    """
    lat, lon, k = _validate_waypoints("polyline_normal", "vertices", lats, lons)

    if method not in METHODS:
        raise ValueError(f"polyline_normal: unknown method={method!r}! Must be one of {METHODS}")

    pole = np.abs(lat) == 90
    if np.any(pole):
        i = int(np.flatnonzero(pole)[0])
        raise ValueError(f"polyline_normal: vertex {i} is at a pole (lat={lat[i]}), where east and north are undefined")

    if method == "great_circle":
        valid, start, end = _great_circle_tangents(lat, lon)
    elif method == "rhumb":
        valid, start, end = _rhumb_tangents(lat, lon)
    else:
        valid, start, end = _linear_tangents(lat, lon)

    if not np.any(valid):
        raise ValueError("polyline_normal: all the segments have zero length, the direction is undefined")

    # nearest non-degenerate segment before (i < j) and after (i >= j) each vertex j
    nseg = k - 1
    index = np.arange(nseg)
    before = np.maximum.accumulate(np.where(valid, index, -1))
    after = np.minimum.accumulate(np.where(valid, index, nseg)[::-1])[::-1]

    seg_in = np.concatenate([[-1], before])  # for the vertices 0..k-1
    seg_out = np.concatenate([after, [nseg]])
    has_in = seg_in >= 0
    has_out = seg_out < nseg

    t_in = np.where(has_in[:, None], end[np.where(has_in, seg_in, 0)], 0.0)
    t_out = np.where(has_out[:, None], start[np.where(has_out, seg_out, 0)], 0.0)
    total = t_in + t_out

    norm = np.hypot(total[:, 0], total[:, 1])
    uturn = has_in & has_out & (norm < _EPSILON)
    if np.any(uturn):
        i = int(np.flatnonzero(uturn)[0])
        raise ValueError(f"polyline_normal: the line reverses its direction at vertex {i} (U-turn)")

    tangent = total / norm[:, None]
    return -tangent[:, 1], tangent[:, 0]


def _real_scalar(name, value):
    """Return ``value`` as a finite float. It must be a real-valued numeric scalar (not a bool)."""
    arr = np.asarray(value)
    if arr.ndim != 0 or arr.dtype.kind not in "iuf":
        raise ValueError(
            f"circle: {name} must be a real-valued numeric scalar! Got {type(value).__name__}"
            + (f" with ndim={arr.ndim}" if arr.ndim else "")
        )
    out = float(arr)
    if not np.isfinite(out):
        raise ValueError(f"circle: {name} must be finite! Got {out}")
    return out


def circle(lat, lon, radius, n=100, figure=IFS_SPHERE):
    """Compute points on a circle of a given radius around a location.

    Parameters
    ----------
    lat: float
        Latitude of the centre (degrees). Must be in the range [-90, 90].
    lon: float
        Longitude of the centre (degrees).
    radius: float
        Radius of the circle (km), measured along the surface of the Earth. Must be
        positive and smaller than the distance to the antipode of the centre (about
        20015.8 km on the default sphere).
    n: int, optional
        Total number of points in the result, including the repeated last point
        (default: 100). Must be at least 4.
    figure: :class:`geo.figure.Figure`, optional
        Figure of the spheroid (default: :obj:`geo.figure.IFS_SPHERE`). Only its radius
        (in metres) is used.

    Returns
    -------
    ndarray
        Latitudes (degrees) of the ``n`` points.
    ndarray
        Longitudes (degrees) of the ``n`` points.

    Raises
    ------
    ValueError
        If ``lat``, ``lon`` or ``radius`` is not a real-valued numeric scalar (a bool, a
        complex number, a string and an array with dimensions are not accepted) or is not
        finite; if ``lat`` is outside [-90, 90]; if ``radius`` is not positive; if the
        radius reaches the antipode of the centre; if ``n`` is not an int; or if ``n`` is
        smaller than 4.

    Notes
    -----
    The radius is the distance along the surface of the sphere (the great circle
    distance), not the straight line distance through the Earth. The angular radius is
    ``radius * 1000 / figure.radius``.

    The circle is closed: the last point is an exact copy of the first one. The ``n - 1``
    distinct points are equally spaced in bearing, ``360 * j / (n - 1)`` degrees for
    ``j = 0 .. n - 2``, as seen from the centre. The first point is due north of the
    centre and the points go clockwise seen from above: north, east, south, west.

    Longitude range: the longitudes of all the points are returned in [0, 360) if the
    longitude of the centre is non-negative and in [-180, 180) otherwise, as in
    :py:func:`polyline`. A circle that crosses the dateline therefore jumps in longitude.

    The circle may contain a pole. The centre can be at a pole, where the bearings are
    undefined: they are taken as the limit of approaching the pole along the meridian
    ``lon``. The first point (due north) is then on the meridian ``lon + 180`` at the
    north pole, as north continues beyond the pole, and on the meridian ``lon`` at the
    south pole. All the points have the same latitude, ``90 - delta`` (north pole) or
    ``-90 + delta`` (south pole), where ``delta`` is the angular radius in degrees.

    The model is a sphere. The points are computed with vectors, which is accurate for
    any radius and next to the poles.

    Because the points travel clockwise, :py:func:`polyline_normal` on the result gives
    normals pointing outwards (away from the centre), wherever it is defined. It is not
    defined if a vertex is exactly at a pole (e.g. a circle around (0, 0) with a radius
    of a quarter of the circumference starts at the north pole) or if the segments are
    degenerate (smaller than the coincidence threshold of :py:func:`polyline_normal`).

    Examples
    --------
    Points at 1000 km from the point (0, 0), going north, east, south and west and back
    to north. The longitude of the centre is non-negative, so the longitudes are in
    [0, 360) and the western point is at 351.0071°:

    >>> from earthkit.geo.polyline import circle
    >>> lats, lons = circle(0, 0, 1000, n=5)
    >>> [round(float(x), 4) + 0.0 for x in lats]
    [8.9929, 0.0, -8.9929, 0.0, 8.9929]
    >>> [round(float(x), 4) + 0.0 for x in lons]
    [0.0, 8.9929, 0.0, 351.0071, 0.0]

    The first and the last point are the same, and the result has ``n`` points:

    >>> lats, lons = circle(51.5, -0.1, 500)
    >>> lats.shape
    (100,)
    >>> bool(lats[0] == lats[-1] and lons[0] == lons[-1])
    True

    """
    lat = _real_scalar("lat", lat)
    lon = _real_scalar("lon", lon)
    radius = _real_scalar("radius", radius)

    if abs(lat) > 90:
        raise ValueError(f"circle: lat must be in the range [-90, 90]! Got {lat}")

    if radius <= 0:
        raise ValueError(f"circle: radius must be positive! Got {radius}")

    delta = radius * 1000.0 / figure.radius
    if delta >= np.pi:
        raise ValueError(
            f"circle: radius={radius} km reaches the antipode of the centre "
            f"(at most {np.pi * figure.radius / 1000.0} km is allowed)"
        )

    if isinstance(n, bool) or not isinstance(n, numbers.Integral):
        raise ValueError(f"circle: n must be an int! Got {type(n).__name__}")

    n = int(n)
    if n < 4:
        raise ValueError(f"circle: n must be at least 4! Got {n}")

    centre = np.stack(latlon_to_xyz(lat, lon))
    east, north = (v[0] for v in _local_basis(np.array([lat]), np.array([lon])))

    theta = np.deg2rad(360.0 * np.arange(n - 1) / (n - 1))
    tangent = np.cos(theta)[:, None] * north + np.sin(theta)[:, None] * east
    p = np.cos(delta) * centre + np.sin(delta) * tangent

    lats = np.rad2deg(np.arctan2(p[:, 2], np.hypot(p[:, 0], p[:, 1])))
    lons = _wrap_lon(np.rad2deg(np.arctan2(p[:, 1], p[:, 0])), 0.0 if lon >= 0 else -180.0)

    # close the circle with an exact copy of the first point
    return np.append(lats, lats[0]), np.append(lons, lons[0])
