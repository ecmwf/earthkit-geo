.. _mir-regrid-array:

regrid (array-level) with MIR
===============================================

*New in version 1.0.0.*

.. py:function:: regrid(data, in_grid=None, out_grid=None, *, interpolation='linear', backend="mir")
    :noindex:

    Regrid array ``data`` using **MIR** (Meteorological Interpolation and Regridding).

    The ``backend`` parameter is set to "mir" by default so it is not necessary to specify it explicitly.

    :param data: array representing a single field defined on the ``in_grid``.
    :type data: ndarray
    :param in_grid: the :ref:`gridspec <gridspec>` describing the grid that ``data`` are defined on. Ignored when ``data`` is not an ndarray.
    :type in_grid: dict, str,  :obj:`Grid <earthkit.geo.grid.Grid>`
    :param out_grid: the :ref:`gridspec <gridspec>` describing the target grid that ``data`` will be interpolated onto
    :type out_grid: dict, str, :obj:`Grid <earthkit.geo.grid.Grid>`
    :param interpolation: the interpolation method. There is a high degree of customisation available to parametrise the available interpolation methods. Please note ot all the interpolation methods support all possible grid types. The possible values are  as follows:

        - "linear": Finite Element based interpolation with linear base functions with supporting triangular mesh
        - "grid-box-average": input/output grid box (see [model_grid_box]_) intersections interpolation preserving input value integrals (conservative interpolation).
        - "nearest-neighbour": choose a nearest neighbouring input point to define output point value

    :type interpolation: str

    :return: Return a tuple with the interpolated values and the :ref:`gridspec <gridspec>` of the output grid.
    :rtype: tuple of ndarray and dict

Examples
--------

- :ref:`/how-tos/mir/mir_numpy_array.ipynb`
