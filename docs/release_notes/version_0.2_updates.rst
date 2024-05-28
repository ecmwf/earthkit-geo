Version 0.2 Updates
/////////////////////////

Version 0.2.0
===============

New features
++++++++++++++++

- added methods :py:meth:`geo.rotate.rotate` and :py:meth:`geo.rotate.unrotate` to perform spehrical rotation. See the notebook example: :ref:`/examples/rotate.ipynb`
- added methods :py:meth:`geo.rotate.rotate_vector` and :py:meth:`geo.rotate.unrotate_vector` to perform local rotation of vectors.
- added methods :py:meth:`geo.coord.latlon_to_xyz` and :py:meth:`geo.coord.xyz_to_latlon` to convert between [ECEF]_ and geodetic coordinates.
- renamed :py:attr:`geo.constants.NORTH` to :py:attr:`geo.constants.NORTH_POLE_LAT` and added :py:attr:`geo.constants.SOUTH_POLE_LAT`