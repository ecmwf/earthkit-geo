API Reference Guide
===================

Regidding with MIR
-------------------------

.. list-table::
   :header-rows: 1

   * - Function
     - Description
   * - :ref:`regrid() <mir-regrid-high>`
     - Regrid high-level data (FieldList, Xarray etc.) with the MIR backend.
   * - :ref:`regrid() <mir-regrid-array>`
     - Regrid array data with the MIR backend.


Regidding with precomputed weights
--------------------------------------

.. list-table::
   :header-rows: 1

   * - Function
     - Description
   * - :ref:`regrid() <precomputed-regrid-high>`
     - Regrid high-level data (FieldList, Xarray etc.) with precomputed weights.
   * - :ref:`regrid() <precomputed-regrid-array>`
     - Regrid array data with precomputed weights.


The full API reference guide
-------------------------------

Traverse through the API reference guide to find out more about the available functions
and classes in `earthkit-geo`.

- :doc:`autoapi/earthkit/geo/cartography/index`
- :doc:`autoapi/earthkit/geo/constants/index`
- :doc:`autoapi/earthkit/geo/coord/index`
- :doc:`autoapi/earthkit/geo/distance/index`
- :doc:`autoapi/earthkit/geo/figure/index`
- :doc:`autoapi/earthkit/geo/gisco/index`
- :doc:`autoapi/earthkit/geo/grids/index`
- :doc:`autoapi/earthkit/geo/grids/xarray/index`
- :doc:`autoapi/earthkit/geo/polyline/index`
- :doc:`autoapi/earthkit/geo/rotate/index`

.. The toctree below puts the generated API pages into the navigation. It starts at
   the ``earthkit.geo`` package index, which lists the modules linked above. The
   ``_regrid`` subpackage is an implicit namespace package, so autoapi generates no
   index page for it and its packages are listed here instead.

.. toctree::
   :hidden:
   :maxdepth: 1

   autoapi/earthkit/geo/index
   autoapi/earthkit/geo/grids/_regrid/array/index
   autoapi/earthkit/geo/grids/_regrid/backends/index
   autoapi/earthkit/geo/grids/_regrid/data/index
   autoapi/earthkit/geo/grids/_regrid/input/index
   autoapi/earthkit/geo/grids/_regrid/regrid/index
