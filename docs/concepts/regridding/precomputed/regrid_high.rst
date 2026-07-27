.. _precomputed-regrid-high:

regrid (high-level) with precomputed weights
=============================================================

.. py:function:: regrid(data, in_grid=None, out_grid=None, in_dims=None, out_dims=None, interpolation="linear", backend="precomputed", inventory="ecmwf")
    :noindex:

    Regrid the high-level ``data`` object (with geography information) using precomputed weights.

    :param data:  The data to be regridded. The following input data types are supported:

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
    :type in_grid: dict, str, :class:`eckit.geo.Grid`

    :param out_grid: the :ref:`gridspec <gridspec>` describing the target grid that ``data`` will be interpolated onto. For the list of supported grids, please refer to the :ref:`gridspec <gridspec>` documentation.
    :type out_grid: dict, str, :class:`eckit.geo.Grid`
    :param in_dims: the names of the geographical dimensions in the Xarray input data. It is only needed when the dimension names cannot be automatically inferred. When it is provided, it takes precedence over the metadata of the input data.
    :type in_dims: list, tuple, None
    :param out_dims: the names of the geographical dimensions in the Xarray output data. It is only needed when the dimension names cannot be automatically inferred. When it is provided, it takes precedence over the metadata of the input data.
    :type out_dims: list, tuple, None
    :param interpolation: the interpolation method. Possible values are ``linear`` and ``nearest-neighbour``. For ``nearest-neighbour`` the following aliases are also supported: ``nn``, ``nearest-neighbor``.
    :type interpolation: str
    :param inventory: the path to the inventory of the precomputed weights. The interpolation only works when the weights are available for the given input grid (automatically determined from the data), target ``grid`` and ``interpolation`` combination. At present, two inventory types are available:

       - If ``inventory`` is "ecmwf" on None, the remote inventory managed by ECMWF is used. In this case the weights are automatically downloaded and stored in a local cache (at ``"~/.cache/earthkit-geo"``) and when it is needed again the cached version is used. See the :ref:`inventory <precomputed_inventory>` for the list of supported grid to grid combinations with this backend.
       - If ``inventory`` is a local path, a local inventory is used. Please note this in experimental feature only used for development purposes.
    :type inventory: str
    :return: Return the regridded data with the same type as ``data`` but with the grid changed to the output grid.
    :rtype:  The same type of data as ``data``.

    :raises ValueError: if the precomputed weights are not available


    The regridding is performed by multiplying the ``data`` vector with the interpolation weights, which forms a sparse matrix (sparse matrix-vector multiplication).



Notes
-------

The interpolation only works when the weights are available for the given input grid (automatically determined from the data), target ``grid`` and ``interpolation`` combination. See the :ref:`inventory <precomputed_inventory>` for the available combinations with the "ecmwf" inventory.

The regridding is performed by multiplying the ``data`` vector with the interpolation weights, which forms a sparse matrix (sparse matrix-vector multiplication).


Examples
--------

- :ref:`/how-tos/precomputed/precomp_healpix_fieldlist.ipynb`
- :ref:`/how-tos/precomputed/precomp_octahedral_fieldlist.ipynb`
- :ref:`/how-tos/precomputed/memory_cache.ipynb`
