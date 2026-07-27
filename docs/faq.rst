Frequently Asked Questions
==========================

How can I interpolate data from O2560 to O1280 grid?
-----------------------------------------------------------

It is supported by :py:func:`regrid` when used with the :ref:`MIR backend <mir-backend>`. Here is an example of how to use it for GRIB data:

.. code-block:: python

    import earthkit.geo as ekg

    # assume fl is a GRIB fieldlist/field on the O2560 grid

    # the target gridspec
    out_grid = {"grid": "O1280"}

    # regrid the fieldlist/field to the target grid, the default backend is MIR. The output
    # is stored in memory.
    fl_res = ekg.regrid(fl, out_grid=out_grid)


How can I use conservative regridding?
--------------------------------------

It is supported by :py:func:`regrid` when used with the :ref:`MIR backend <mir-backend>` and the ``interpolation=grid-box-average`` option. Here is an example of how to use it for GRIB data:

.. code-block:: python

    import earthkit.data as ekd
    import earthkit.geo as ekg

    # get fieldlist from a sample GRIB file
    ds = ekd.from_source("sample", "O32_t2.grib2")

    # the target gridspec (regular latitude-longitude grid)
    out_grid = {"grid": [5, 5]}

    # regrid the fieldlist/field to the target grid, the default backend is MIR. The output
    # is stored in memory.
    ds_res = ekg.regrid(ds, out_grid=out_grid, interpolation="grid-box-average")


How can I specify the cache directory for MIR?
--------------------------------------------------

The cache directory for MIR is managed by MIR itself and not by ``earthkit-geo``. It can be specified by setting the ``MIR_CACHE_PATH`` environment variable to the desired path.
