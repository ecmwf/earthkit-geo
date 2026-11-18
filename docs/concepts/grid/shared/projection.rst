.. _grid-projection:

Projection
==========

A *projection* maps points from a source to a target coordinate type, forwards (``fwd``) and backwards (``inv``). The coordinate types are geographic ``(longitude, latitude)``, planar ``(x, y)`` and geocentric ``(x, y, z)``. Map projections such as ``lcc`` or ``laea`` map ``(longitude, latitude)`` to ``(x, y)``, the ``rotation`` maps ``(longitude, latitude)`` to ``(longitude, latitude)`` on a rotated sphere, ``ll-to-xyz`` maps them to geocentric coordinates, and PROJ-based projections map between whatever coordinate reference systems they are given.

**The default projection** is the equidistant cylindrical (plate carrée) projection with ``x`` the longitude and ``y`` the latitude in degrees, which is the identity for all practical purposes. It is omitted from specs, and a projection that does not change anything (such as a rotation with the south pole at the South Pole and no rotation angle) becomes the default projection.

**How grids use projections.** A grid defines its points in its own frame, and its projection maps them to geographic coordinates. There are two cases:

- *Latitude–longitude grids* (regular, Gaussian, HEALPix, ...) can be rotated: their points are defined on a rotated sphere, and the rotation maps them to geographic coordinates.

  .. code-block:: python

      Grid(grid="F48", rotation=[30, 30]).spec
      # {'grid': 'F48', 'projection': {'south_pole': [30, 30], 'type': 'rotation'}}

- *Projected grids* (``regular_xy``) define their points in a plane, in metres, and their projection maps ``(x, y)`` to ``(longitude, latitude)``. A map projection is the other way round, so it is used reversed, with ``reverse: true``:

  .. code-block:: python

      Grid(
          x=[0, 4995000, 5000],
          y=[0, -4745000, -5000],
          false_xy=[1818500.021, -2287499.955],
          projection={"type": "laea", "lat_0": 52, "lon_0": 10, "figure": "grs80", "reverse": True},
      )  # the EFAS 5 km grid, also available as Grid("LAEA-EFAS-5km")

**``reverse: true``** swaps the forward and backward mappings, and the source and target coordinate types, of any projection, and is kept in its canonical spec. The types ``rotation``, ``stretch`` and ``proj`` interpret it themselves: a reversed rotation, for example, is the inverse rotation, and is written canonically as such.

**False easting and northing** are properties of the grid, not of the projection: a projection is pure mathematics, while the offset of a grid's origin describes the grid. A projection spec therefore rejects PROJ's ``x_0`` and ``y_0``, and a projected grid gives its offsets as ``false_xy: [x, y]``, or as ``false_lonlat: [lon, lat]``, the geographic location of the grid's origin, from which the projection computes the offsets (:ref:`grid-projected`).

**Named projections** are registered in ``projection.yaml`` and can be used wherever a projection spec is expected: ``projection: swisslv95``.

**PROJ strings.** Every projection type with an equivalent in `PROJ <https://proj.org>`_ has a ``proj_str``, computed by eckit without requiring PROJ itself (:ref:`grid-ref-projections` lists the correspondence):

.. code-block:: python

    Projection(type="laea", lat_0=52, lon_0=10, figure="grs80").proj_str
    # '+proj=laea +lat_0=52 +lon_0=10 +ellps=GRS80'

    Projection(type="rotation", south_pole=[10, -40], rotation_angle=25).proj_str
    # '+proj=ob_tran +lon_0=10 +o_lat_p=40 +o_lon_p=25 +o_proj=longlat +R=6371229'

The tests check that each built-in projection gives the same coordinates as PROJ given its ``proj_str``, and that a projected grid gives the same points with either projection.

**Using a projection directly.** ``fwd`` and ``inv`` take coordinates in order or by name, as scalars (one point) or as 1-d arrays (many points). Points outside the domain of the projection result in NaN:

.. code-block:: python

    lcc = Projection(type="lcc", lat_1=63.3, lat_0=63.3, lon_0=15)

    lcc.source_point_coordinates        # ('longitude', 'latitude')
    lcc.target_point_coordinates        # ('x', 'y')

    x, y = lcc.fwd(longitude=[10, 20], latitude=[60, 65])
    lon, lat = lcc.inv(x, y)

.. _grid-projection-proj:

PROJ-based projections
----------------------

The projection type ``proj`` interfaces with the `PROJ <https://proj.org>`_ library. Rather than reimplementing PROJ's projection descriptions, it adopts them: a ``proj`` projection is given by two coordinate reference systems, and PROJ computes the transformation between them.

.. code-block:: yaml

    projection:
      type: proj
      source: "EPSG:2056"       # Swiss CH1903+ / LV95, (x, y) in metres
      target: "EPSG:4326"       # WGS 84, (lon, lat) in degrees; the default

- **Source and target.** ``source`` and ``target`` accept any coordinate reference system definition that PROJ accepts, such as ``AUTHORITY:CODE``, a PROJ string, WKT or PROJJSON (see `proj_create_crs_to_crs <https://proj.org/en/stable/development/reference/functions.html#c.proj_create_crs_to_crs>`_ and `PROJ strings <https://proj.org/en/stable/usage/projections.html>`_). Both default to ``EPSG:4326``, and defaults are omitted from canonical specs. The shorthand ``{proj: "<definition>"}`` gives the target.
- **Axis order** is normalised to longitude first and easting first (`proj_normalize_for_visualization <https://proj.org/en/stable/development/reference/functions.html#c.proj_normalize_for_visualization>`_).
- **Coordinate types** follow from the coordinate systems: ellipsoidal and spherical systems give ``(longitude, latitude)``, two-dimensional Cartesian systems ``(x, y)``, and three-dimensional Cartesian systems ``(x, y, z)``. For example, ``EPSG:4326`` to ``EPSG:3857`` maps ``(longitude, latitude)`` to ``(x, y)``, and ``EPSG:4326`` to ``EPSG:4978`` to geocentric ``(x, y, z)``.
- **Longitudes** are returned in ``[lon_minimum, lon_minimum + 360)``, by default ``[0, 360)``; set ``lon_minimum: -180`` for ``[-180, 180)``.
- **Points** that PROJ fails to transform result in NaN.

**The figure comes from PROJ.** A ``proj`` projection takes its figure from the ellipsoid of the target coordinate reference system, as PROJ defines it. It is a working ``Figure`` with a spec, semi-axes and a PROJ string, which can be used anywhere else a figure is expected:

.. code-block:: python

    etrs89 = Projection(type="proj", target="EPSG:4258")    # ETRS89, on the GRS 1980 ellipsoid
    etrs89.figure.spherical                                 # False
    etrs89.figure.a                                         # 6378137.0
    etrs89.figure.spec                                      # the figure spec, reusable elsewhere

**Projected grids on any CRS.** A ``regular_xy`` grid can be defined in the plane of any projected coordinate reference system. This is how the Swiss national grids are defined in the catalogue:

.. code-block:: yaml

    # projection.yaml
    projection_names:
      - swisslv95:
          type: proj
          source: "EPSG:2056"

    # grid.yaml
    grid_names:
      - swisslv95:
          type: regular_xy
          x: [2439000, 2867000, 1000]
          y: [1040500, 1333500, 1000]
          projection: swisslv95

**Built-in or PROJ?** The built-in projections (:ref:`grid-ref-projections`) need no external library and are available everywhere. The ``proj`` type opens up every coordinate reference system PROJ knows. The two agree: a grid with a built-in projection gives the same points as the same grid with ``projection: {type: proj, source: <its proj_str>}``.

**Availability.** The ``proj`` type requires PROJ and its database, ``proj.db`` (`resource files <https://proj.org/en/stable/resource_files.html>`_); ``eckit.geo.projdb_is_available()`` tells whether it is usable.

.. seealso::

    :ref:`grid-projected` describes projected grids, and :ref:`grid-ref-projections` lists the built-in projection types and their PROJ equivalents.
