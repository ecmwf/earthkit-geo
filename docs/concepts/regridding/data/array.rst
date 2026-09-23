.. _regrid-array:

Regridding NumPy arrays
=========================

For plain NumPy array data — a single field of values with no attached geography metadata — use the
dedicated array-level entry point, :func:`earthkit.geo.grids.array.regrid`, rather than the generic,
high-level :func:`earthkit.geo.grids.regrid`.

.. code-block::

    from earthkit.geo.grids import array

    out_values, out_grid = array.regrid(values, in_grid={"grid": [1, 1]}, out_grid={"grid": "O320"})

Since a raw array carries no geography information of its own, both ``in_grid`` and ``out_grid`` must
always be given explicitly as :ref:`gridspecs <gridspec>`. The function returns a tuple of the
interpolated values and the output grid spec; the latter may differ from the requested ``out_grid``
since the backend can normalise or otherwise adjust it during regridding.

.. note::

    Passing a raw NumPy array to the generic :func:`earthkit.geo.grids.regrid` function raises a
    ``ValueError`` pointing to this array-level entry point instead, since that function's grid
    inference and dispatch mechanism only applies to data types with attached geography metadata
    (see :ref:`regrid-fieldlist`, :ref:`regrid-grib` and :ref:`regrid-xarray`).

Examples
--------

Using the "mir" backend:

- :ref:`/how-tos/mir/mir_numpy_array.ipynb`

Using the "precomputed" backend:

- :ref:`/how-tos/precomputed/precomp_numpy_array.ipynb`

See also
--------

- :ref:`mir-regrid-array`
- :ref:`precomputed-regrid-array`
- :ref:`gridspec`
