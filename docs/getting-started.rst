
Installation and Getting Started
================================

Installing from PyPi
--------------------

Install **earthkit-geo** with python3 (>= 3.10) and ``pip`` as follows:

.. code-block:: bash

    python3 -m pip install earthkit-geo


For more details on the installation please see :ref:`install`.


Import and use
--------------

For regridding, a different interface is available depending on the input data type.

High-level interface
//////////////////////

Use it for data containing geographical information, e.g. earthkit-data :xref:`fieldlist` objects, Xarray DataArrays or Datasets.

.. code-block:: python

    import earthkit.data as ekd
    from earthkit.geo.regrid import regrid

    # get fieldlist from a sample GRIB file
    ds = ekd.from_source("sample", "O32_t2.grib2")

    # the target is a regular latitude-longitude grid
    grid = {"grid": [5, 5]}

    ds_res = regrid(ds, grid=grid)


Array-level interface
////////////////////////

Use it for raw data arrays, e.g. Numpy ndarrays.


.. code-block:: python

    from earthkit.geo.regrid.array import regrid
    import numpy as np

    vals = np.random.rand(320, 640)
    in_grid = {"grid": [0.25, 0.25]}  # regular latitude-longitude grid
    out_grid = {"grid": "O320"}  # octahedral reduced Gaussian grid

    res_vals, res_grid = regrid(vals, in_grid=in_grid, out_grid=out_grid)
