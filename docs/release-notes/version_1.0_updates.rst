..  _release-notes-1.0:

Version 1.0 Updates
/////////////////////


Version 1.0.0
===============

This is the first stable release of ``earthkit-geo``.

Fixes
------

- Fixed import examples in docstrings for :func:`~earthkit.geo.distance.haversine_distance`, :func:`~earthkit.geo.distance.nearest_point_haversine` and :func:`~earthkit.geo.distance.nearest_point_kdtree` to use ``from earthkit.geo.distance import ...`` (:pr:`76`).
- Updated regridding how-to notebooks to use the ``out_grid`` kwarg (:pr:`77`).


Regridding interface
---------------------

Added the ``in_grid`` and ``out_grid`` kwargs to the high-level :func:`regrid` methods :ref:`regrid() <mir-regrid-high>` and :ref:`regrid() <precomputed-regrid-high>` (:pr:`73`). These kwargs are used to specify the input and output grid specifications for regridding. The previously used ``grid`` kwarg is now deprecated and replaced with ``out_grid``. However, it can still be used for backward compatibility but it will raise a deprecation warning. The ``in_grid`` kwarg can be used to specify the input grid specification if it cannot be inferred from the input data, e.g. for Xarray input without proper metadata.

An example of how to use the new interface:

.. code-block:: python

    import earthkit.data as ekd
    import earthkit.geo as ekg

    # get fieldlist from a sample GRIB file
    ds = ekd.from_source("sample", "O32_t2.grib2").to_fieldlist()

    # the target is a regular latitude-longitude grid
    out_grid = {"grid": [5, 5]}

    ds_res = ekg.regrid(ds, out_grid=out_grid)


Breaking changes for Xarray regridding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the input grid cannot be automatically inferred from the input data, it can now be specified via the new ``in_grid`` kwarg (:pr:`73`). Previously, it required adding the grid specification to the Xarray data as a custom attribute ``earthkit_grid_spec``, which was not straightforward for users. This attribute is no longer supported.

See the notebook :ref:`/tutorials/mir_regrid_xarray.ipynb` for an example.

.. code-block:: python

    import earthkit.geo as ekg
    import xarray as xr

    # create a sample Xarray DataArray without grid metadata
    da = xr.DataArray(...)

    # specify the input grid specification via the new in_grid kwarg
    in_grid = {"grid": [0.25, 0.25]}  # regular latitude-longitude grid

    # specify the output grid specification via the new out_grid kwarg
    out_grid = {"grid": "O320"}  # octahedral reduced Gaussian grid

    da_res = ekg.regrid(da, in_grid=in_grid, out_grid=out_grid)


Importing
----------

Changed the way methods/objects can be imported from ``earthkit.geo`` (:pr:`68`):

- The high-level :func:`regrid` method can still be imported from the top level, but other methods have to be imported from their submodules.
- The ``grids`` submodule was added to contain all grid-related methods.
- As a new feature, the ``eckit.geo.Grid`` object is available as ``earthkit.geo.grids.Grid``. Please note this is an experimental feature.

.. code-block:: python
    :caption: The new way to import methods/objects from earthkit.geo

    # high-level regrid method
    from earthkit.geo import regrid

    # array level regrid methods
    from earthkit.geo.grids.array import regrid

    # new feature: eckit.geo.Grid object
    from earthkit.geo.grids import Grid

    # other examples
    from earthkit.geo.distance import haversine_distance
    from earthkit.geo.figure import Sphere

Breaking changes
+++++++++++++++++

The array-level :func:`regrid` method is now available as :func:`earthkit.geo.grids.array.regrid`. The old way of importing it from the top level is no longer supported.

.. code-block:: python
   :caption: How to import the array-level regrid method from earthkit.geo

    # new way
    from earthkit.geo.grids.array import regrid

    # old way - no longer supported
    from earthkit.geo.regrid.array import regrid

The following methods are no longer available at the top level and should be imported from their respective submodules:

- methods in the ``distance`` module (e.g. :func:`haversine_distance`) should be imported from :mod:`earthkit.geo.distance`:

   - GeoKDTree
   - haversine_distance
   - nearest_point_haversine
   - nearest_point_kdtree

- methods/objects in the ``figure`` module (e.g. :class:`Sphere`) should be imported from :mod:`earthkit.geo.figure`:

   - IFS_SPHERE
   - UNIT_SPHERE
   - Sphere


Regridding
-----------

Added support for :ref:`regridding` with :ref:`MIR <mir-backend>` or :ref:`precomputed weights <precomputed-backend>` backends (:pr:`rc1`).

- Regridding GRIB data requires earthkit-data version 1.0.0rc4 or later.
- Implemented regridding for non-GRIB fields/fieldlists (:pr:`62`).
- Ensured that field metadata is preserved during regridding (:pr:`62`).

See the following guides for more details:

Regridding with the MIR backend:

- :ref:`/how-tos/mir/mir_numpy_array.ipynb`
- :ref:`/how-tos/mir/mir_healpix_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_octahedral_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_interpolation_types.ipynb`

Regridding with the precomputed weights backend:

- :ref:`/how-tos/precomputed/precomp_numpy_array.ipynb`
- :ref:`/how-tos/precomputed/precomp_healpix_fieldlist.ipynb`
- :ref:`/how-tos/precomputed/precomp_octahedral_fieldlist.ipynb`


Other fixes
-----------

- Fixed issue when URL paths were not properly encoded in the precomputed weights backend, which prevented downloading index files and weights on Windows (:pr:`49`).
- Ensured that no eckit.geo related warnings appear when an index file is loaded for the "precomputed" regridding backend (:pr:`60`).
- Added missing dependency (:pr:`55`).
