.. _grid-regular-ll:

Regular latitude–longitude grids
================================

A regular lattice of latitudes and longitudes, given by its increments, plus the
specialised Arakawa C-grid types built on the same geometry.

Increments and reference point
------------------------------

A regular latitude–longitude grid is given by its increments, ``grid: [dlon, dlat]``
in degrees. Its points lie on the lines through a **reference** point, by default
``(0, 0)``: the grid includes the reference, and every point a whole number of
increments away from it. A reference offset by half an increment gives a
cell-centred grid:

.. code-block:: python

    Grid(grid=[1, 1], reference=[0.5, 0.5])     # points at ..., -0.5, 0.5, 1.5, ...

Given an area without a reference, the grid passes through the area's south-west
corner (see :ref:`grid-spec-canonical`).

Arakawa C-grids
---------------

Two specialised types describe **Arakawa C-grids** on regular latitude–longitude
geometries, where different variables are staggered by half a cell:

- ``arakawa_c`` is given by ``n``, the number of cells per 90° (like the Gaussian
  ``N``; see :ref:`grid-gaussian`), so that the increments are 90°/n, optionally
  scaled by ``grid_factor: [fx, fy]``. The ``arrangement`` places the points at the
  cell centres (``T``, the default), or at the midpoints of the cell sides (``U`` on
  the west and east sides, ``V`` on the south and north sides).

- ``arakawa_c_um`` describes the grids of the Met Office Unified Model, with
  ``grid_factor: [2, 4/3]`` and order ``i+j+``. For example, the UM N96 grid has
  increments of 1.875° by 1.25°.

The hydrological grids of the catalogue are Arakawa C-grids, such as ``GLOFAS-4``
(3 arcminutes, near-global). Being written in terms of ``regular_ll``, these are
ordinary regular grids once built (see :ref:`grid-spec-canonical`).

.. _grid-staggering:

Staggering: arrangement
-----------------------

Models often place different variables at different positions of the same mesh
(Arakawa staggering). The **arrangement** selects which positions a grid's points
are at, and ``grid.arrangement`` returns it:

.. list-table::
   :header-rows: 1

   * - Grids
     - Arrangements
   * - ``arakawa_c``, ``arakawa_c_um``
     - ``T`` (cell centres, default), ``U``, ``V``
   * - ORCA
     - ``T``, ``U``, ``V``, ``F``, ``W``

Arrangement letters are case-insensitive. For ORCA grids the arrangement is part of
the name, as in ``eORCA025_T``.

.. seealso::

    :ref:`grid-ref-grids` for the spec keys of ``regular_ll``, ``arakawa_c`` and
    ``arakawa_c_um``, and :ref:`grid-orca` for the ORCA arrangements.
