.. _regrid-xarray:

Regridding xarray Dataset and DataArray
=========================================

:func:`regrid` also accepts an :class:`xarray.DataArray` or :class:`xarray.Dataset` as its ``data``
argument. A ``DataArray`` input yields a ``DataArray`` output; a ``Dataset`` is regridded variable by
variable and a ``Dataset`` is returned.

Determining the input grid
----------------------------

Unless ``in_grid`` is given explicitly, the input grid of each variable is determined as follows:

- if the dataset carries an earthkit-data ``"earthkit.grid_spec"`` attribute, that grid spec is used
  as-is ; otherwise
- the grid spec is inferred from the coordinates and CF metadata of the variable, when they describe a
  grid earthkit-geo can build and whose points they reproduce: a regular or rotated latitude-longitude
  grid, a regular ``x``/``y`` grid in the plane of a projection (:ref:`grid-projected`), a HEALPix or a
  reduced Gaussian grid ; otherwise
- as a fallback, the latitude/longitude values of every point are extracted and used to build a generic :ref:`unstructured <gridspec-unstructured>` grid regardless of the original coordinate layout. See :ref:`regridding-unstructured` for more details.

The inference needs the metadata that states the grid to be part of the object handed to :func:`regrid`:
a single ``DataArray`` cut out of a CF dataset no longer carries the grid mapping variable of the file,
so its projected or rotated grid cannot be inferred from it. The fallback needs its points to be
available too: either latitude/longitude coordinates, or a usable projection to derive them from. A
``DataArray`` of a file which holds nothing but the projected ``x``/``y`` axes has neither once its grid
mapping is lost, and cannot be regridded: keep the grid mapping variable as a scalar coordinate of the
``DataArray``, or pass the ``Dataset``.

A projected grid held by nothing but the 1D ``x``/``y`` axes of its plane, with no latitude/longitude
coordinates, is detected when its grid mapping variable carries ``crs_wkt``, or ``proj4_params``
together with the individual CF parameters: pyproj does not read ``proj4_params``, so a mapping
carrying only that one takes the fallback.

A reduced Gaussian grid held by nothing but a ``reduced_gaussian_index`` coordinate, with no
latitude/longitude coordinates of its own, is detected from its CF ``reduced_gaussian`` grid mapping
variable (``grid_mapping_name = "reduced_gaussian"``) and the table of points per latitude it names
(``points_per_latitude`` or ``accumulated_points_per_latitude``). An octahedral grid (``grid_subtype =
"octahedral"``) gives ``{"grid": "O<N>"}``. A classical one (``grid_subtype = "normal"``) gives ``{"grid":
"N<N>"}``, but only when its points per latitude are those of the classical grid of that resolution.
Mapping and table always describe the *global* grid, whichever part of it the data is on: it is the
``reduced_gaussian_index``, the number each point has on the global grid, west to east within a latitude
line and north to south from line to line, which says which points the data holds. The direction of the
latitude table does not matter, but the index has to be a strictly increasing sequence of numbers of the
global grid, and ``longitude_of_first_meridian`` has to be 0. What the index gives depends on the points
it numbers:

- every point of the global grid, numbered ``0, 1, 2, ...`` with no gaps: that grid, in the point order of
  the CF convention;
- the points of a rectangle of longitudes on uninterrupted latitude lines: a sub-area of it, as
  ``{"grid": "O<N>", "area": [north, west, south, east]}`` (:ref:`grid-area`). The area is an equivalent
  covering the same points, not necessarily the one the file was written with, and the values are
  reordered when it crosses the Greenwich meridian: eckit runs every latitude line of an area from its
  west edge, where the CF file runs it from the meridian. The MIR backend clips the result to the
  bounding box of the area of its input, so a sub-area input gives a regional output;
- any other subset of points, the land points only for instance: the :ref:`unstructured
  <gridspec-unstructured>` grid of their latitudes and longitudes, computed from the table of the global
  grid (see :ref:`regridding-unstructured`).

The mapping takes precedence over the per-point latitude/longitude coordinates CF also allows: a variable
holding both is detected from the mapping, and its values number the points the index names, in the order
of the file. Both the mapping and the table have to reach :func:`regrid` for that: a ``DataArray`` cut out
of a CF dataset holds the index but neither of them, and regridding it raises an error unless it also
holds the latitude and longitude of its points. Pass the ``Dataset``, or the grid yourself in ``in_grid``
with the point dimension in ``in_dims``: the values are then regridded as they stand, with no reordering.

An explicitly provided ``in_grid`` always takes precedence over what can be inferred from the dataset.

Determining the input dimensions
----------------------------------

Similarly, unless ``in_dims`` is given explicitly, the geographical dimensions of each variable are
determined alongside its input grid, from whichever coordinates (grid spec attribute, latitude/longitude,
or x/y) were used to establish that grid.

An explicitly provided ``in_dims`` always takes precedence over what can be inferred from the dataset.

Regridding each variable
--------------------------

Each variable is regridded independently with ``xarray.apply_ufunc`` (dask-aware where
applicable), which delegates the actual point-to-point interpolation to the selected backend. The
output geography (dimensions and coordinate arrays) is rebuilt from the resulting output grid and
attached to the result; the earthkit-data grid spec attribute is updated too, where present.

Output grids
-------------------------------------

The output grid is always given by the user via ``out_grid`` (unlike the input grid, it cannot be
inferred from the dataset). By default, the dimension names attached to the regridded result depend
on the shape of that output grid, as described below. This default can be overridden by giving
``out_dims`` explicitly, in which case those dimension names are used instead, regardless of the
output grid's shape.

**Meshed grid.** When the output grid is two-dimensional (a regular mesh), the result gets
separate ``latitude`` and ``longitude`` dimensions, each with its own 1D coordinate array.

If distinct 1D latitude/longitude coordinate arrays cannot be recovered for the output grid, this
falls back to generic ``y``/``x`` dimensions instead, carrying the full 2D latitude/longitude
coordinate arrays.

**1D grid.** When the output grid is one-dimensional (e.g. an unstructured grid), the result is
instead given a single ``values`` dimension, carrying 1D latitude/longitude coordinates indexed by
that same dimension.

**Reduced Gaussian grid.** A one-dimensional output grid which is a plain octahedral ``O<N>`` or
classical ``N<N>`` reduced Gaussian grid, global or with an ``area``, is written with the CF
``reduced_gaussian`` grid mapping rather than as ``values``: its points are numbered by a
``reduced_gaussian_index`` dimension coordinate holding their index on the global grid, and the values are
gathered into that increasing order. A ``Dataset`` carries the whole encoding, located by the global grid
alone, with no per-point latitude/longitude. The result of regridding a field onto ``{"grid": "O16",
"area": [60, 350, 30, 20]}``, as xarray prints it with the byte counts and the attributes left out:

.. code-block:: text

    Dimensions:                 (reduced_gaussian_index: 26, lat: 32)
    Coordinates:
      * reduced_gaussian_index  (reduced_gaussian_index) int32 140 141 142 ... 438 439
      * lat                     (lat) float64 85.76 80.27 74.74 ... -80.27 -85.76
    Data variables:
        t                       (reduced_gaussian_index) float64 200.5 200.5 ... 201.0
        pl                      (lat) int32 20 24 28 32 36 ... 32 28 24 20
        reduced_gaussian        int64 0

``pl`` counts the points of every latitude line of the global grid, ``lat`` gives their latitudes, and
the scalar ``reduced_gaussian`` variable is the grid mapping (``grid_mapping_name = "reduced_gaussian"``,
``grid_subtype`` ``"octahedral"`` or ``"normal"``, ``latitude_dimension = "lat"`` and
``points_per_latitude = "pl"``) which every regridded variable refers to with its ``grid_mapping``
attribute. Such a result written with :meth:`~xarray.Dataset.to_netcdf` and read back is therefore
detected again as the grid it was regridded onto.

A ``DataArray`` can hold neither the table nor the mapping, so its encoding is only partial: the index
coordinate, the latitude and longitude of every point, and no ``grid_mapping`` attribute — enough to
regrid it again, which is not the case of an index-only ``DataArray`` cut out of a CF file. Any other
one-dimensional target keeps the ``values`` form above, including a reduced Gaussian grid with a
rotation, a projection or a table of its own. In every case the ``grid_mapping`` attribute the input
carried is not passed on to the result. Since the encoding writes ``lat``, ``pl`` and
``reduced_gaussian`` into the result, regridding onto such a target raises a ``ValueError`` when the data
already takes one of those names.


Examples
--------

- :ref:`/tutorials/regrid/mir_regrid_xarray.ipynb`
- :ref:`/tutorials/regrid/mir_regrid_xarray_unstructured.ipynb`

See also
--------

- :ref:`mir-regrid-high`
- :ref:`precomputed-regrid-high`
- :ref:`grid-spec`
