.. _mir-regrid-high:

regrid (high-level) with MIR
===============================================

*New in version 1.0.0.*

.. py:function:: regrid(data, grid=None, *, interpolation='linear',  backend="mir")
    :noindex:

    Regrid the high-level ``data`` object (with geography information) using **MIR** (Meteorological Interpolation and Regridding).

    The ``backend`` parameter is set to "mir" by default so it is not necessary to specify it explicitly.

    :param data: the following input data types are supported:

        - an earthkit-data GRIB :py:class:`~earthkit.data.core.fieldlist.FieldList` (requires :xref:`earthkit-data` >= 1.0.0).
        - an earthkit-data GRIB :py:class:`~earthkit.data.core.field.Field` (requires :xref:`earthkit-data` >= 1.0.0).
        - a GRIB message as a bytes or :class:`io.BytesIO` object
        - an :class:`xarray.DataArray` or :class:`xarray.Dataset`. Please note this is an **experimental** feature and
          only works if the input :ref:`gridspec <gridspec>` can be automatically inferred from the data. If the
          Xarray was created with ``earthkit-data`` it most probably contains the necessary metadata to be automatically
          inferred. Otherwise, the Xarray dataset/dataarray must have the ``ek_grid_spec`` attribute containing
          the :ref:`gridspec <gridspec>` as a dict or a JSON string.

    :type data: :py:class:`~earthkit.data.core.fieldlist.FieldList`, :py:class:`~earthkit.data.core.field.Field`, bytes, :class:`io.BytesIO`, :class:`xarray.DataArray`, :class:`xarray.Dataset`
    :param grid: the :ref:`gridspec <gridspec>` describing the target grid that ``data`` will be interpolated onto.
    :type grid: dict, str, :class:`eckit.geo.Grid`
    :param interpolation: the interpolation method. Please note not all the interpolation methods support all possible grid types. The possible values are as follows:

        - "linear": Finite Element based interpolation with linear base functions with supporting triangular mesh
        - "grid-box-average": input/output grid box intersections interpolation preserving input value integrals (conservative interpolation).
        - "nearest-neighbour": choose a nearest neighbouring input point to define output point value

    :type interpolation: str


    The interpolation only works if both the input and output grid are supported. For the list of supported grids, please refer to the :ref:`gridspec <grid_spec>` documentation.

    :return: Return the regridded data with the same type as ``data`` but with the grid changed to the output grid.
    :rtype:  The same type of data as ``data``.
Examples
--------

- :ref:`/how-tos/mir/mir_healpix_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_octahedral_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_interpolation_types.ipynb`
