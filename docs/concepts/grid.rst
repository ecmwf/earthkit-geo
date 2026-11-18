.. _grid-object:

The Grid object
===============

This section describes the ``Grid`` object of earthkit-geo (``earthkit.geo.grids.Grid``), the
specification, or *spec*, that fully describes a grid, and the concepts a grid is built from: the
figure of the Earth, the projection, the area and the order of the points. It complements the
:ref:`MIR gridspecs <mir-gridspec>` and :ref:`precomputed gridspecs <precomputed-gridspec>` pages,
which list the gridspecs the regridding backends accept.

``earthkit.geo.grids.Grid`` is currently identical to ``eckit.geo.Grid``, the Python interface to
eckit's geometry library, which implements everything described here. The figures, projections and
areas of :ref:`grid-shared` are available from the same module, ``eckit.geo``. The string and dict
grid specs used in :ref:`regridding <mir-backend>` are all converted to a ``Grid`` object
internally.

.. warning::

    ``Grid`` is an **experimental feature** and its interface is still evolving: it may change in
    future releases. This section describes what is implemented today.

What a Grid is
--------------

A grid is a set of points on a surface, in a well-defined order. earthkit-geo represents a grid
with a ``Grid`` object, which is:

- **Completely described by its spec.** From the spec alone, any library, process or machine
  reconstructs the same points, in the same order.
- **Functional.** A ``Grid`` is not a passive record of metadata. It computes the coordinates of its
  points (geographic, geocentric or in its native projected plane), its bounding box, cropped and
  reordered versions of itself, and nearest-neighbour searches; it also hands out its projection and
  its figure as working objects.
- **Immutable.** A grid never changes once built; a cropped or reordered grid is a new object.
  Expensive derived data, such as coordinates or search trees, are computed on first use, so
  creating a grid is cheap, even a very large one.

A first example
~~~~~~~~~~~~~~~

A ``Grid`` is created from its spec, which can be given as keyword arguments, a dict or a string,
or just from a grid name. See :ref:`grid-creating` for all the ways to create a grid.

.. code-block:: python

    from earthkit.geo.grids import Grid

    grid = Grid(grid="O320")

    grid.type           # 'reduced_gg'
    grid.shape          # (421120,)
    grid.spec           # {'grid': 'O320'}
    grid.uid            # 32 hexadecimal digits, the same in every process

    lat, lon = grid.to_latlons()

The spec ``{"grid": "O320"}`` names the octahedral reduced Gaussian grid with 320 latitudes between
pole and equator. The name is enough: the grid knows its type, its 640 latitudes, the number of
points on each, the order in which they are numbered, and its identifier.

What a Grid offers
~~~~~~~~~~~~~~~~~~

.. list-table::
    :header-rows: 1

    * - Property or method
      - Returns
    * - ``spec``, ``spec_str``
      - The canonical spec, as a dict or as compact JSON (:ref:`grid-spec-canonical`)
    * - ``uid``
      - A 32-hexadecimal-digit identifier (:ref:`grid-identity`)
    * - ``type``, ``name``, ``arrangement``
      - The type of the grid, and for catalogue grids their name and staggering arrangement
        (:ref:`grid-staggering`)
    * - ``shape``, ``size()``, ``len(grid)``
      - The shape of the grid, e.g. ``(181, 360)`` for a regular grid or ``(421120,)`` for a reduced
        one, and its number of points
    * - ``order``
      - The order of the points (:ref:`grid-order`)
    * - ``point_coordinates``
      - The names of the native coordinates, ``("longitude", "latitude")`` or ``("x", "y")``
    * - ``to_points()``
      - The native coordinates, as a dict of coordinate name to NumPy array
    * - ``to_latlons()``
      - Geographic coordinates on the figure, as ``(latitudes, longitudes)``
    * - ``to_xyz()``
      - Geocentric Cartesian coordinates on the figure, as ``(x, y, z)``
    * - ``lat()``, ``lon()``, ``x()``, ``y()``
      - Distinct coordinate values, for grids built from ranges (e.g. regular grids)
    * - ``bounding_box()``
      - ``(north, west, south, east)``
    * - ``figure``, ``projection``
      - Working ``Figure`` and ``Projection`` objects (:ref:`grid-figure` and :ref:`grid-projection`)
    * - ``crop(area)``, ``cropped(area)``
      - Indices of the points in an area, or the cropped grid itself (:ref:`grid-area`)
    * - ``reorder(order)``
      - The renumbering to another order (:ref:`grid-order`)
    * - ``to_unstructured_ll(name="")``, ``cropped_unstructured_ll(area)``
      - The same points, or those in an area, as an unstructured grid (:ref:`grid-as-unstructured`)
    * - ``search(**spec)``
      - A search on the points: nearest neighbour, k nearest neighbours, or within a radius
    * - ``catalog``
      - What the catalogue knows about the grid, if anything (:ref:`grid-catalogues`)

.. _grid-identity:

Identity: spec, uid and equality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A grid has two identities, both derived from its points.

The **canonical spec** (``grid.spec``) describes the grid in its simplest form. Two grids are equal
when their canonical specs are equal, however each was created:

.. code-block:: python

    Grid(grid=[1, 1]) == Grid("1/1")                                       # True
    Grid(grid=[5, 5], rotation=[-40, -20]) == Grid(
        grid=[5, 5], projection={"type": "rotation", "south_pole": [-20, -40]}
    )                                                                      # True

The **uid** (``grid.uid``) is a 128-bit identifier written as 32 hexadecimal digits. For grids
generated from parameters it is a digest of the canonical spec. For grids defined by coordinate
data, such as unstructured and ORCA grids, it is a digest of the coordinates themselves, so it
identifies the points regardless of how they were obtained. Every grid registers its uid when it is
computed, so ``Grid(uid=...)`` rebuilds the grid in the same process, and catalogue grids, such as
ORCA grids, are registered under their published uids.

Where to go next
----------------

To learn how a grid is described and created, start with :ref:`grid-spec` and :ref:`grid-creating`.
The concepts that grids are built from are explained in :ref:`grid-shared` and :ref:`grid-types`.
When you need the details of a particular spec, look it up in :ref:`grid-reference`.

.. toctree::
    :maxdepth: 2

    grid/spec
    grid/creating
    grid/shared/index
    grid/types/index
    grid/reference/index
