Frequently asked questions
==========================

How can I interpolate data from O2560 to O1280 grid?
-----------------------------------------------------------

It is supported by :py:func:`regrid` when used with the :ref:`MIR backend <mir-backend>`. Here is an example of how to use it for GRIB data:

.. code-block:: python

    from earthkit.geo.regrid import regrid

    # assume fl is a GRIB fieldlist/field on the O2560 grid

    # the target gridspec
    grid = {"grid": "O1280"}

    # regrid the fieldlist/field to the target grid, the default backend is MIR. The output
    # is stored in memory.
    fl_res = regrid(fl, grid=grid)


How can I use conservative regridding?
--------------------------------------

It is supported by :py:func:`regrid` when used with the :ref:`MIR backend <mir-backend>` and the ``interpolation=grid-box-average`` option. Here is an example of how to use it for GRIB data:

.. code-block:: python

    import earthkit.data as ekd
    from earthkit.geo.regrid import regrid

    # get fieldlist from a sample GRIB file
    ds = ekd.from_source("sample", "O32_t2.grib2")

    # the target gridspec (regular latitude-longitude grid)
    grid = {"grid": [5, 5]}

    # regrid the fieldlist/field to the target grid, the default backend is MIR. The output
    # is stored in memory.
    ds_res = regrid(ds, grid=grid, interpolation="grid-box-average")


How can I specify the cache directory for MIR?
--------------------------------------------------

The cache directory for MIR is managed by MIR itself and not by ``earthkit-geo``. It can be specified by setting the ``MIR_CACHE_PATH`` environment variable to the desired path.
