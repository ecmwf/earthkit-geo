.. _release-notes-1.1:

Version 1.1 Updates
/////////////////////////

Version 1.1.1
===============

Fixes
++++++++++++++++

- Fixed an issue in Xarray :ref:`regridding <regridding>` where dataset and
  variable attributes were not propagated to the output, causing failures with
  Xarray 2025.10.1 (:pr:`100`).



Version 1.1.0
===============

New features
++++++++++++++++

- Added support for non-trailing geographical dimensions in Xarray
  :ref:`regridding <regridding>` (:pr:`89`). Geographical dimensions are now
  detected by name rather than position. Input arrays are made C-contiguous before being passed to MIR,
  since dimension reordering can produce non-contiguous views that MIR rejects;
  this is a no-op for arrays that are already contiguous.
