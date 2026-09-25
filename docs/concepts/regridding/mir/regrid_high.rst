.. _mir-regrid-high:

regrid (high-level) with MIR
===============================================

.. py:function:: regrid(data, in_grid=None, out_grid=None, in_dims=None, out_dims=None, interpolation='linear',  backend="mir")
    :noindex:

    Regrid the high-level ``data`` object (with geography information) using **MIR** (Meteorological Interpolation and Regridding).

    The ``backend`` parameter is set to "mir" by default so it is not necessary to specify it explicitly.

    :param data: The data to be regridded. The following input data types are supported:

        - earthkit-data :py:class:`~earthkit.data.core.fieldlist.FieldList` (see details :ref:`here <regrid-fieldlist>`)
        - earthkit-data :py:class:`~earthkit.data.core.field.Field` (see details :ref:`here <regrid-fieldlist>`)
        - :class:`xarray.DataArray` or :class:`xarray.Dataset` (see details :ref:`here <regrid-xarray>`)
        - GRIB message as a bytes or :class:`io.BytesIO` object (see details :ref:`here <regrid-grib-message>`)

    :type data: :py:class:`~earthkit.data.core.fieldlist.FieldList`, :py:class:`~earthkit.data.core.field.Field`, :class:`xarray.DataArray`, :class:`xarray.Dataset`, bytes, :class:`io.BytesIO`
    :param in_grid: The :ref:`gridspec <gridspec>` describing the input grid. When None (the default), the input grid is inferred from the input data if possible. If grid information cannot be inferred, but the latitudes and longitudes are available the input grid is treated as an unstructured lat/lon grid. For the list of supported grids, please refer to the :ref:`gridspec <gridspec>` documentation. Ignored when ``data`` is a GRIB message.
    :type in_grid: dict, str, :py:class:`Grid`
    :param out_grid: The :ref:`gridspec <gridspec>` describing the target grid that ``data`` will be interpolated onto. For the list of supported grids, please refer to the :ref:`gridspec <gridspec>` documentation.
    :type out_grid: dict, str, :py:class:`Grid`
    :param in_dims: The names of the geographical dimensions in the Xarray input data. It is only needed when the dimension names cannot be automatically inferred. When it is provided, it takes precedence over the metadata of the input data.
    :type in_dims: list, tuple, None
    :param out_dims: The names of the geographical dimensions in the Xarray output data. It is only needed when the dimension names cannot be automatically inferred. When it is provided, it takes precedence over the metadata of the input data.
    :type out_dims: list, tuple, None
    :param interpolation: The interpolation method. Please note not all the interpolation methods support all possible grid types. The possible values are as follows:

        - "linear": Finite Element based interpolation with linear base functions with supporting triangular mesh
        - "grid-box-average": input/output grid box intersections interpolation preserving input value integrals (conservative interpolation).
        - "nearest-neighbour": choose a nearest neighbouring input point to define output point value

    :type interpolation: str


    :return: The regridded data with the same type as ``data`` but with the grid changed to the output grid.
    :rtype:  The same type of data as ``data``.


Notes
-----

The interpolation only works if both the input and output grid are supported. For the list of supported grids, please refer to the :ref:`gridspec <gridspec>` documentation.


Examples
--------

- :ref:`/how-tos/mir/mir_healpix_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_octahedral_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_interpolation_types.ipynb`
- :ref:`/tutorials/mir_regrid_xarray.ipynb`
