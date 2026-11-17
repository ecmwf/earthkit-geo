.. _grid-projected:

Projected grids: x, y and the false origin
==========================================

A projected grid in the plane of a map projection, described by coordinate ranges
and a false origin that places it on the Earth.

Coordinates in the projection plane
-----------------------------------

A projected grid (``regular_xy``) is a regular grid in the plane of a projection.
Its coordinates are given as ranges, ``x: [start, stop, step]`` and
``y: [start, stop, step]``, in the units of the projection plane (usually metres);
``stop`` is included, and a negative step runs backwards. Alternatively,
``grid: [dx, dy]`` with ``bounding_box_xy: [min_x, min_y, max_x, max_y]`` and
``order`` give the same grid. The directions of ``x`` and ``y`` determine the order,
which can also make ``j`` the fastest-varying index.

The false origin
----------------

The **false origin** relates the grid's coordinates to those of the projection. With
``false_xy: [x0, y0]`` (false easting and northing), grid point ``(x, y)`` is at
``(x − x0, y − y0)`` in the projection plane. With ``false_lonlat: [lon, lat]``, the
grid's origin, ``(x, y) = (0, 0)``, is placed at that geographic location, and the
projection computes the offsets.

The projection
--------------

The projection must map ``(x, y)`` to ``(longitude, latitude)``: built-in map
projections are used with ``reverse: true``, and ``proj`` projections with a
projected ``source``. Without a projection, the points are plain Cartesian
coordinates, with no geographic meaning.

.. code-block:: python

    Grid(
        x=[0, 2120000, 2500],
        y=[0, 2420000, 2500],
        false_lonlat=[1.60491, 51.633816],
        projection={"type": "lcc", "lat_0": 63.3, "lat_1": 63.3, "lon_0": 15, "figure": {"R": 6371000}, "reverse": True},
    )  # the MEPS grid, also available as Grid("meps-bris-nordic")

.. seealso::

    :ref:`grid-ref-grids` for the spec keys of ``regular_xy``, and
    :ref:`grid-projection` for projections.
