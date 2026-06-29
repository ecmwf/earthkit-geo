..  _release-notes-1.0:

Version 1.0 Updates
/////////////////////

Version 1.0.0
===============

Regridding
-----------

Added support for :ref:`regridding` with :ref:`MIR <mir-backend>` or :ref:`precomputed weights <precomputed-backend>` backends. The input can be a Numpy array, an earthkit-data :py:class:`~earthkit.data.core.field.Field` or :py:class:`~earthkit.data.core.fieldlist.FieldList` object, or an Xarray DataArray or Dataset.


See the following guides for more details:

Regridding with the MIR backend:

- :ref:`/how-tos/mir/mir_numpy_array.ipynb`
- :ref:`/how-tos/mir/mir_healpix_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_octahedral_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_interpolation_types.ipynb`
- :ref:`/tutorials/mir_regrid_xarray.ipynb`

Regridding with the precomputed weights backend:

- :ref:`/how-tos/precomputed/precomp_numpy_array.ipynb`
- :ref:`/how-tos/precomputed/precomp_healpix_fieldlist.ipynb`
- :ref:`/how-tos/precomputed/precomp_octahedral_fieldlist.ipynb`


Importing
----------

Changed the way methods/objects can be imported from ``earthkit.geo``. The following methods are no longer available at the top level and should be imported from their respective submodules:

- methods in the ``distance`` module (e.g. :func:`haversine_distance`) should be imported from :mod:`earthkit.geo.distance`:

   - GeoKDTree
   - haversine_distance
   - nearest_point_haversine
   - nearest_point_kdtree

- methods/objects in the ``figure`` module (e.g. :class:`Sphere`) should be imported from :mod:`earthkit.geo.figure`:

   - IFS_SPHERE
   - UNIT_SPHERE
   - Sphere
