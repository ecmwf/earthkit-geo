.. _grids:

Grid Types
==========

This page provides an overview of all grid types supported by earthkit-geo.

The :class:`~earthkit.geo.grids.Grid` class provides factory methods that automatically return the appropriate grid subclass:

.. code-block:: python

    from earthkit.geo import Grid

    # Automatically returns HEALPix instance
    grid = Grid.from_string("H32")

    # Automatically returns Equirectangular instance
    grid = Grid.from_dict({"grid": [1, 1]})

.. autoclass:: earthkit.geo.grids.Equirectangular

----

.. autoclass:: earthkit.geo.grids.HEALPix

----

.. autoclass:: earthkit.geo.grids.ReducedGaussian

----

.. autoclass:: earthkit.geo.grids.Octahedral

----

.. autoclass:: earthkit.geo.grids.ICON_CH1_EPS

----

.. autoclass:: earthkit.geo.grids.ICON_CH2_EPS
