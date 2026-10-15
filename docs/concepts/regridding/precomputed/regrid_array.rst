.. _precomputed-regrid-array:

regrid (array-level) with precomputed weights
=============================================

.. py:function:: regrid(date, in_grid=None, out_grid=None, interpolation='linear', backend="precomputed", inventory="ecmwf")
    :noindex:

    Regrid array ``data`` using precomputed weights.

    :param data: The data to be regridded, represented as a NumPy array, defining a single field on the ``in_grid``.
    :type data: ndarray
    :param in_grid: The :ref:`gridspec <gridspec>` describing the grid that ``data`` are defined on.
    :type in_grid: dict
    :param out_grid: The :ref:`gridspec <gridspec>` describing the target grid that ``data`` will be interpolated onto
    :type out_grid: dict
    :param interpolation: The interpolation method. Possible values are ``linear`` and ``nearest-neighbour``. For ``nearest-neighbour`` the following aliases are also supported: ``nn``, ``nearest-neighbor``.
    :type interpolation: str
    :param inventory: the path to the inventory of the precomputed weights. The interpolation only works when the weights are available for the given ``in_grid``, ``out_grid`` and ``interpolation`` combination. At present, two inventory types are available:

       - If ``inventory`` is "ecmwf" or None, the remote inventory managed by ECMWF is used. In this case the weights are automatically downloaded and stored in a local cache (at ``"~/.cache/earthkit-geo"``) and when it is needed again the cached version is used. See the :ref:`inventory <precomputed_inventory>` for the list of supported grid to grid combinations with this backend.
       - If ``inventory`` is a local path, a local inventory is used. Please note this an experimental feature only used for development purposes.
    :type inventory: str
    :return: Return a tuple with the interpolated values and the :ref:`gridspec <gridspec>` of the output grid.
    :rtype: tuple of ndarray and dict
    :raises ValueError: if the precomputed weights are not available


Notes
----------

The interpolation only works when the weights are available for the given input grid (automatically determined from the data), target ``grid`` and ``interpolation`` combination. See the :ref:`inventory <precomputed_inventory>` for the available combinations with the "ecmwf" inventory.

The regridding is performed by multiplying the ``values`` vector with the interpolation weights, which forms a sparse matrix (sparse matrix) -vector multiplication).


Examples
--------

- :ref:`/how-tos/precomputed/precomp_numpy_array.ipynb`
- :ref:`/how-tos/precomputed/memory_cache.ipynb`
