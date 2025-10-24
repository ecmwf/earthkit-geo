Version 0.4 Updates
/////////////////////////

Version 0.4.0
===============

New features
++++++++++++++++

- Added the :py:mod:`earthkit.geo.gisco` module for retrieving, caching, and loading shapefiles from the Eurostat GISCO API (:pr:`22`)


Installation
++++++++++++++++++++++++

- Made `cartopy` an optional dependency as "[cartography]". It is required for the :py:meth:`earthkit.geo.cartography.country_polygons` method. See :ref:`install` for details on how to install ``earthkit-geo`` with optional dependencies.
