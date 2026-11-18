.. _grid-shared:

Concepts shared by all grids
============================

A grid is made of points in a given order (determined by its type), on a surface (the *figure*), possibly mapped through a *projection*, and possibly restricted to an *area*. Each of these concepts has its own spec, nested in the grid's:

.. code-block:: text

    grid spec
    ├── type                     the grid type, often implied
    ├── type-specific keys       grid, pl, Nside, x, y, …
    ├── order                    the numbering of the points
    ├── area                     a bounding box
    ├── projection               a projection spec
    │   ├── type                 the projection type
    │   ├── type-specific keys   lat_0, lat_1, south_pole, source, …
    │   ├── reverse
    │   └── figure               a figure spec
    └── figure                   shorthand, when the projection is the default

The grid type is described in :ref:`grid-spec-type`, its type-specific keys in :ref:`grid-types` and :ref:`grid-ref-grids`, and the projection's type-specific keys in :ref:`grid-ref-projections`. The shared concepts each have their own page:

.. toctree::
   :maxdepth: 1

   figure
   projection
   area
   order
