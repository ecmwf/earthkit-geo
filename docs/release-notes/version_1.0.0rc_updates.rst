..  _release-notes-1.0.0rc:

Version 1.0.0 Release Candidate Updates
/////////////////////////////////////////

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
