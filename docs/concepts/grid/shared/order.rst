.. _grid-order:

Order
=====

The points of a grid are numbered, and the *order* says how. Iteration, ``to_latlons()``, ``to_points()`` and the values of any field on the grid all follow this order, so it is what ties a data array to the points it belongs to.

**Regular grids** (``regular_ll``, ``regular_gg``, ``regular_xy``) support every scanning mode. An order such as ``i+j-`` reads: ``i`` (along x, or longitude) varies fastest, in increasing direction (``+``), and ``j`` (along y, or latitude) varies slowest, in decreasing direction (``-``). The default, ``i+j-``, numbers the rows from north to south, each from west to east, which is ECMWF's and GRIB's default scanning mode. Writing ``j`` first makes columns the fastest-varying index, and ``+-`` or ``-+`` makes consecutive rows (or columns) alternate in direction, starting in the direction of the first sign, like a boustrophedon:

.. list-table:: Scanning modes for regular grids
   :header-rows: 1
   :stub-columns: 1

   * -
     - Rows (``i`` fastest)
     - Columns (``j`` fastest)
   * - One direction
     - ``i+j-`` (default), ``i+j+``, ``i-j-``, ``i-j+``
     - ``j+i+``, ``j+i-``, ``j-i+``, ``j-i-``
   * - Alternating
     - ``i+-j-``, ``i+-j+``, ``i-+j-``, ``i-+j+``
     - ``j+-i+``, ``j+-i-``, ``j-+i+``, ``j-+i-``

**Reduced Gaussian grids** (``reduced_gg``) are numbered by rows, in one direction: ``i+j-`` (default), ``i+j+``, ``i-j-`` or ``i-j+``.

**HEALPix grids** are numbered in ``ring`` order (the default) or ``nested`` order (:ref:`grid-healpix`).

**Grids defined by coordinates** are numbered in the order of their coordinates.

Non-default orders appear in the canonical spec, and ``grid.reorder(order)`` returns the renumbering to another order: point ``k`` of the grid is point ``ren[k]`` in that order.

.. code-block:: python

    grid = Grid(grid="H2")
    ren = grid.reorder("nested")

    Grid(grid=[1, 1], order="i+j+").spec
    # {'grid': [1, 1], 'order': 'i+j+'}

.. seealso::

    :ref:`grid-ref-grids` lists the order accepted by each grid type, and :ref:`grid-spec-canonical` explains how non-default orders appear in canonical specs.
