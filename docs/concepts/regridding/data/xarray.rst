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

**1D grid.** When the output grid is one-dimensional (e.g. an unstructured or reduced grid), the
result is instead given a single ``values`` dimension, carrying 1D latitude/longitude coordinates
indexed by that same dimension.


Examples
--------

- :ref:`/tutorials/regrid/mir_regrid_xarray.ipynb`
- :ref:`/tutorials/regrid/mir_regrid_xarray_unstructured.ipynb`

See also
--------

- :ref:`mir-regrid-high`
- :ref:`precomputed-regrid-high`
- :ref:`grid-spec`
