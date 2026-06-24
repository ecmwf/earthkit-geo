.. _grid-object:

The Grid object
===========================

earthkit-geo uses the ``earthkit.geo.grids.Grid`` object to represent a grid. This is currently identical to the ``eckit.geo.Grid`` object from the ``eckit`` library. The str and dict grid specs used in `regridding <mir-backend>`_ are all converted to a ``Grid`` object internally.

Currently, ``Grid`` is an **experimental feature**  and its interface may change in future releases. For that reason, it is not documented here. If you are interested in using it, you can instantiate a ``Grid`` object from a grid spec and explore its interface.

.. code-block:: python

    from earthkit.geo.grids import Grid

    g = Grid({"grid": "O320"})
    g = Grid("O320")
