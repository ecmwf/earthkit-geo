.. _mir-regrid-high:

regrid (high-level) with MIR
===============================================

.. py:function:: regrid(data, in_grid=None, out_grid=None, in_dims=None, out_dims=None, interpolation='linear',  backend="mir")
    :noindex:

    Regrid the high-level ``data`` object (with geography information) using **MIR** (Meteorological Interpolation and Regridding).

    The ``backend`` parameter is set to "mir" by default so it is not necessary to specify it explicitly.

    :param data: The data to be regridded. The following input data types are supported:

        - earthkit-data :py:class:`~earthkit.data.core.fieldlist.FieldList`
        - earthkit-data :py:class:`~earthkit.data.core.field.Field`
        - :class:`xarray.DataArray` or :class:`xarray.Dataset`
        - GRIB message as a bytes or :class:`io.BytesIO` object

    :type data: :py:class:`~earthkit.data.core.fieldlist.FieldList`, :py:class:`~earthkit.data.core.field.Field`, :class:`xarray.DataArray`, :class:`xarray.Dataset`, bytes, :class:`io.BytesIO`
    :param in_grid: the :ref:`gridspec <gridspec>` describing the input grid. It is only needed when the input grid
        cannot be automatically inferred from the input data. This can be the case for
        Xarray, where at present the grid information can only be accessed via the "earthkit.grid_spec" attribute,
        and if it is missing or does not contain the necessary information ``in_grid`` needs to be provided. When ``in_grid`` is provided, it takes precedence over the metadata of the
        input data. For the list of supported grids, please refer to the :ref:`gridspec <gridspec>` documentation.
    :type in_grid: dict, str, :py:class:`Grid`
    :param out_grid: the :ref:`gridspec <gridspec>` describing the target grid that ``data`` will be interpolated onto. For the list of supported grids, please refer to the :ref:`gridspec <gridspec>` documentation.
    :type out_grid: dict, str, :py:class:`Grid`
    :param in_dims: the names of the geographical dimensions in the Xarray input data. It is only needed when the dimension names cannot be automatically inferred. When it is provided, it takes precedence over the metadata of the input data.
    :type in_dims: list, tuple, None
    :param out_dims: the names of the geographical dimensions in the Xarray output data. It is only needed when the dimension names cannot be automatically inferred. When it is provided, it takes precedence over the metadata of the input data.
    :type out_dims: list, tuple, None
    :param interpolation: the interpolation method. Please note not all the interpolation methods support all possible grid types. The possible values are as follows:

        - "linear": Finite Element based interpolation with linear base functions with supporting triangular mesh
        - "grid-box-average": input/output grid box intersections interpolation preserving input value integrals (conservative interpolation).
        - "nearest-neighbour": choose a nearest neighbouring input point to define output point value

    :type interpolation: str


    :return: Return the regridded data with the same type as ``data`` but with the grid changed to the output grid.
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
