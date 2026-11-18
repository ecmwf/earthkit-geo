.. _grid-figure:

Figure
======

The *figure* is the shape and size of the reference surface: a sphere or an oblate spheroid (an ellipsoid of revolution). Geographic coordinates are geodetic latitude and longitude on the figure, and the figure determines everything that depends on the size and shape of the Earth: geocentric coordinates (``to_xyz()``), the ellipsoidal formulas of map projections, surface areas, and distances.

A figure is given by name, or by its dimensions:

.. code-block:: python

    Figure(figure="wgs84")                                  # by name
    Figure(R=6371229)                                       # a sphere; also r or radius
    Figure(a=6378137, b=6356752.314245)                     # an oblate spheroid; also semi_major_axis, semi_minor_axis
    Figure(figure={"a": 6378137, "b": 6356752.314245})      # nested

Dimensions are in metres. A spheroid with ``a == b`` is a sphere, and prolate spheroids (``b > a``) are rejected. The named figures are listed in :ref:`grid-ref-figures`; the default is ``earth``, the sphere of radius 6 371 229 m used by the IFS.

The figure belongs to the projection: a projection spec carries its figure, and ``grid.figure`` is ``grid.projection.figure``. When the projection is the default, the figure can be given at the top level of the grid spec, as a shorthand:

.. code-block:: python

    Grid(grid=[30, 30], figure="wgs84").spec
    # {'figure': 'wgs84', 'grid': [30, 30]}

    Grid(grid=[30, 30], figure="earth").spec
    # {'grid': [30, 30]}                     the default figure is omitted

In a canonical spec a figure is written by name when it matches a named figure, and otherwise as ``{r: ...}`` or ``{a: ..., b: ...}``. Figures are equal when their semi-axes are equal.

A ``Figure`` object offers its dimensions (``R``, ``a``, ``b``), ``spherical``, ``eccentricity`` and ``flattening``, the surface area of the whole figure or of a bounding box (``area()``, ``area(bbox)``), its spec, and its PROJ representation (``proj_str``, e.g. ``+ellps=WGS84`` or ``+R=6371229``). PROJ-based projections take their figure from PROJ, and the result is a figure like any other (see :ref:`PROJ-based projections <grid-projection-proj>`).

.. seealso::

    :ref:`grid-ref-figures` lists the named figures, and :ref:`grid-projection` explains how a figure is attached to a projection.
