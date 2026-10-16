.. _release-notes-1.2:

Version 1.2 Updates
/////////////////////////


Version 1.2.0
===============

New features
++++++++++++++++

- Added support for :ref:`regridding <regridding>` to/from :ref:`unstructured <gridspec-unstructured>` grids (:pr:`113`, :pr:`114`, :pr:`118`, :pr:`121`, :pr:`124`, :pr:`125`).
- Automatically detects the geography related dimensions and extracts the corresponding latitude/longitude values for a wider range of Xarray data for :ref:`regridding <regridding>`.
- Added fallback mechanism for :ref:`fieldlist <regrid-fieldlist>` and :ref:`Xarray <regrid-xarray>` when the input grid spec cannot be determined. In this case the latitude/longitude values of every point are extracted and used to build a generic :ref:`unstructured <gridspec-unstructured>` grid regardless of the original coordinate layout. See :ref:`regridding-unstructured` for more details and also the following notebook example:

    - :ref:`/tutorials/mir_regrid_xarray_unstructured.ipynb`

Fixes
++++++++++++++++

- Removed unwanted error-level logging emitted while guessing the grid mapping
  for Xarray :ref:`regridding <regridding>`; these messages are now logged at
  debug level (:pr:`118`, :pr:`122`).

Other changes
++++++++++++++++

- Refactored GRIB regridding and reorganised the Xarray regridding code into
  a dedicated subpackage (:pr:`112`, :pr:`117`).
- Added new concept pages for each input type for regridding:

   - :ref:`regrid-xarray`
   - :ref:`regrid-fieldlist`
   - :ref:`regrid-grib-message`
   - :ref:`regrid-array`

- Various documentation improvements for :ref:`regridding <regridding>`
  (:pr:`93`).
