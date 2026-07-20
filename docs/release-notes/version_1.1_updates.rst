.. _release-notes-1.1:

Version 1.1 Updates
/////////////////////////

Version 1.1.0
===============

New features
++++++++++++++++

- Added support for non-trailing geographical dimensions in Xarray
  :ref:`regridding <regridding>` (:pr:`89`). Geographical dimensions are now
  detected by name rather than position. Input arrays are made C-contiguous before being passed to MIR,
  since dimension reordering can produce non-contiguous views that MIR rejects;
  this is a no-op for arrays that are already contiguous.
