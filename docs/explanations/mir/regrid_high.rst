.. _mir-regrid-high:

regrid (high-level) with MIR
===============================================

*New in version 1.0.0.*

.. py:function:: regrid(data, grid=None, *, interpolation='linear',  backend="mir")
    :noindex:

    Regrid the high-level ``data`` object (with geography information) using **MIR** (Meteorological Interpolation and Regridding).

    The ``backend`` parameter is set to "mir" by default so it is not necessary to specify it explicitly.

    :param data: the following input data types are supported:

        - an earthkit-data GRIB :xref:`fieldlist` (requires :xref:`earthkit-data` >= 1.0.0).
        - an earthkit-data GRIB :xref:`field` (requires :xref:`earthkit-data` >= 1.0.0).
        - a GRIB message as a bytes or :class:`io.BytesIO` object
        - an :class:`xarray.DataArray` or :class:`xarray.Dataset`

    :type data: :xref:`fieldlist`, :xref:`field`, bytes, or :class:`io.BytesIO`
    :param grid: the :ref:`grid_spec <grid_spec>` describing the target grid that ``data`` will be interpolated onto
    :type grid: dict, str, :class:`Grid`
    :param interpolation: the interpolation method. There is a high degree of customisation available to parametrise the available interpolation methods. Please note ot all the interpolation methods support all possible grid types. The possible values are as follows:

        - "linear": Finite Element based interpolation with linear base functions with supporting triangular mesh
        - "grid-box-average": input/output grid box (see [model_grid_box]_) intersections interpolation preserving input value integrals (conservative interpolation).
        - "nearest-neighbour": choose a nearest neighbouring input point to define output point value

    :type interpolation: str


Examples
--------

- :ref:`/how-tos/mir/mir_numpy_array.ipynb`
- :ref:`/how-tos/mir/mir_healpix_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_octahedral_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_interpolation_types.ipynb`
