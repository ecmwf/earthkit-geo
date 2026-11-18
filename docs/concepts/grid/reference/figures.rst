.. _grid-ref-figures:

Figure specs
============

The names of the predefined figures, with their shapes and dimensions, and the keys that define a figure by its dimensions. :ref:`grid-figure` explains how figures are used.

The named figures are:

.. list-table::
   :widths: 22 33 45
   :header-rows: 1

   * - Name
     - Shape
     - Dimensions [m]
   * - ``earth`` (default)
     - Sphere, the IFS Earth
     - R = 6 371 229
   * - ``grib1``
     - Sphere, the GRIB edition 1 Earth
     - R = 6 367 470
   * - ``wgs84_sphere``
     - Sphere
     - R = 6 371 200
   * - ``wgs84``
     - WGS 84 ellipsoid
     - a = 6 378 137, b = 6 356 752.314245
   * - ``grs80``
     - GRS 1980 ellipsoid
     - a = 6 378 137, b = 6 356 752.314140
   * - ``iau1965``
     - IAU 1965 ellipsoid
     - a = 6 378 160, b = 6 356 775
   * - ``sun``
     - Sphere
     - R = 695 990 000
   * - ``unit-sphere``
     - Sphere
     - R = 1

A figure can also be given by its dimensions:

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Key
     - Value
   * - ``figure``
     - A name, or a nested figure spec
   * - ``R`` (or ``r``, ``radius``)
     - Radius of a sphere
   * - ``a`` and ``b`` (or ``semi_major_axis`` and ``semi_minor_axis``)
     - Semi-axes of an oblate spheroid, ``a ≥ b``

Canonical form: ``{figure: <name>}`` when a named figure matches, otherwise ``{figure: {r: R}}`` or ``{figure: {a: a, b: b}}``. The default figures, ``earth`` and ``grib1``, are omitted from canonical grid and projection specs. See also `PROJ ellipsoids <https://proj.org/en/stable/usage/ellipsoids.html>`_.
