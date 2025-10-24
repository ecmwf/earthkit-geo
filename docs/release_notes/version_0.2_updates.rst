Version 0.2 Updates
/////////////////////////

Version 0.2.0
===============

New features
++++++++++++++++

- added methods :py:meth:`earthkit.geo.rotate.rotate` and :py:meth:`earthkit.geo.rotate.unrotate` to perform spehrical rotation. See the notebook example: :ref:`/examples/rotate.ipynb`
- added methods :py:meth:`earthkit.geo.rotate.rotate_vector` and :py:meth:`earthkit.geo.rotate.unrotate_vector` to perform local rotation of vectors.
- added methods :py:meth:`earthkit.geo.coord.latlon_to_xyz` and :py:meth:`earthkit.geo.coord.xyz_to_latlon` to convert between [ECEF]_ and geodetic coordinates.
- renamed :py:attr:`earthkit.geo.constants.NORTH` to :py:attr:`earthkit.geo.constants.NORTH_POLE_LAT` and added :py:attr:`earthkit.geo.constants.SOUTH_POLE_LAT`
