.. _grid-ref-projections:

Projection specs
================

The keys of every projection type: the keys common to all projections first, then the projection types and their PROJ equivalents, and the named projections. :ref:`grid-projection` explains how projections work.

Keys common to projections
--------------------------

Every projection understands these keys, whatever its type:

.. list-table::
   :widths: 12 30 18 40
   :header-rows: 1

   * - Key
     - Value
     - Default
     - Meaning
   * - ``type``
     - Projection type
     - ``eqc``
     - The builder
   * - ``figure``
     - Figure spec, or the figure keys ``R``, ``a``, ``b``, ...
     - ``earth``
     - :ref:`grid-figure`
   * - ``reverse``
     - Boolean
     - ``false``
     - Swap ``fwd`` and ``inv`` (:ref:`grid-projection`)

A projection can also be given by name: ``projection: swisslv95``. PROJ's false easting and northing (``x_0``, ``y_0``) are rejected; projected grids use ``false_xy`` instead.

Projection types
----------------

.. list-table::
   :widths: 18 22 42 18
   :header-rows: 1

   * - Type (aliases)
     - Maps
     - Keys (defaults)
     - PROJ equivalent
   * - ``eqc`` (``plate-carree``)
     - (lon, lat) → (x, y) in degrees; the default projection
     - ``lat_ts`` (0), ``lat_0`` (0)
     - `eqc <https://proj.org/en/stable/operations/projections/eqc.html>`_, with ``to_meter`` scaling
   * - ``rotation``
     - (lon, lat) → (lon, lat) on the rotated sphere
     - ``south_pole: [lon, lat]``, or ``rotation: [lat, lon]``, or ``south_pole_lon`` and ``south_pole_lat`` (default ``[0, -90]``); ``rotation_angle`` (0)
     - `ob_tran <https://proj.org/en/stable/operations/projections/ob_tran.html>`_
   * - ``lcc`` (``lambert``, ``lambert_lam``)
     - (lon, lat) → (x, y)
     - ``lat_1``, ``lat_2`` (``lat_1``), ``lon_0`` (0), ``lat_0`` (0)
     - `lcc <https://proj.org/en/stable/operations/projections/lcc.html>`_
   * - ``aea`` (``albers``)
     - (lon, lat) → (x, y)
     - ``lon_0``, ``lat_0``, ``lat_1``, ``lat_2`` (all required)
     - `aea <https://proj.org/en/stable/operations/projections/aea.html>`_
   * - ``laea`` (``lambert_azimuthal_equal_area``)
     - (lon, lat) → (x, y)
     - ``lat_0`` (0), ``lon_0`` (0)
     - `laea <https://proj.org/en/stable/operations/projections/laea.html>`_
   * - ``stere`` (``polar_stereographic``)
     - (lon, lat) → (x, y)
     - ``lat_0`` (90 or −90, required), ``lat_ts`` (``lat_0``), ``lon_0`` (0)
     - `stere <https://proj.org/en/stable/operations/projections/stere.html>`_, polar aspect
   * - ``merc`` (``mercator``)
     - (lon, lat) → (x, y)
     - ``lon_0`` (0), ``lat_ts`` (0)
     - `merc <https://proj.org/en/stable/operations/projections/merc.html>`_
   * - ``stretch`` (``schmidt``)
     - (lon, lat) → (lon, lat), Schmidt transform: :math:`\tan(\mathrm{colat}'/2) = c\,\tan(\mathrm{colat}/2)`
     - ``factor`` (c > 0, required)
     - None
   * - ``ll-to-xyz``
     - (lon, lat) → geocentric (x, y, z) on the figure
     - Figure only
     - `cart <https://proj.org/en/stable/operations/conversions/cart.html>`_
   * - ``proj``
     - Between two coordinate reference systems, computed by PROJ
     - ``source`` (``EPSG:4326``), ``target`` (``EPSG:4326``), ``lon_minimum`` (0); shorthand ``proj: <definition>``
     - Itself (:ref:`grid-projection-proj`)

The map projections follow Snyder, *Map Projections: A Working Manual* (USGS, 1987), on the sphere and on the ellipsoid.

Named projections
-----------------

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Name
     - Spec
   * - ``swisslv95``
     - ``{type: proj, source: "EPSG:2056"}`` (CH1903+ / LV95)
   * - ``swisslv03``
     - ``{type: proj, source: "EPSG:21781"}`` (CH1903 / LV03)

.. seealso::

    :ref:`grid-projection` for how grids use projections, and :ref:`grid-projection-proj` for the ``proj`` type.
