..  _release-notes-1.0.0rc:

Version 1.0.0 Release Candidate Updates
/////////////////////////////////////////


Version 1.0.0rc7
==================

- Use ``earthkit_grid_spec`` instead of ``ek_grid_spec`` as the attribute name for the grid specification in Xarray input for :func:`regrid` (:pr:`69`).


Version 1.0.0rc6
==================

Importing
------------------

Changed the way methods/objects can be imported from `earthkit.geo` (:pr:`68`):

- The high-level :func:`regrid` method can still be imported from the top level, but other methods have to be imported from their submodules.
- The `grids` submodule was added to contain all the grid related methods.
- As a new feature the `eckit.geo.Grid` object is available as `earthkit.geo.grids.Grid`. Please note this is an experimental feature.

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

Breaking changes 1
++++++++++++++++++++

The array-level :func:`regrid` method is now available as :func:`earthkit.geo.grids.array.regrid`. The old way of importing it from the top level is no longer supported.

.. code-block:: python
   :caption: How to import the array-level regrid method from earthkit.geo

    # new way
    from earthkit.geo.grids.array import regrid

    # old way - no longer supported
    from earthkit.geo.regrid.array import regrid

Breaking changes 2
++++++++++++++++++++

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


Version 1.0.0rc5
==================

- Regridding GRIB data requires earthkit-data version 1.0.0rc4 or later
- Implemented regridding for non-GRIB fields/fieldlists (:pr:`62`)
- Ensured that field metadata is preserved during regridding (:pr:`62`)


Version 1.0.0rc4
==================

- Ensured that no eckit.geo related warnings appear when index file is loaded for the "precomputed" regridding backend (:pr:`60`)


Version 1.0.0rc3
==================

- Added missing dependency (:pr:`55`)


Version 1.0.0rc2
==================

- Fixed issue when URL paths were not properly encoded in the precomputed weights backend, which prevented downloading index files and weights on Windows (:pr:`49`)


Version 1.0.0rc1
==================

Regridding
++++++++++++++++

Added support for :ref:`regridding` with :ref:`MIR <mir-backend>` or :ref:`precomputed weights <precomputed-backend>` backends.

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
