.. _regrid-fieldlist:

Regridding FieldList and Field data
====================================

:func:`regrid` accepts an earthkit-data :py:class:`~earthkit.data.core.fieldlist.FieldList` or
:py:class:`~earthkit.data.core.field.Field` as its ``data`` argument. A ``FieldList`` is regridded
field by field; a single ``Field`` is regridded by internally wrapping it in a one-field ``FieldList``.

Determining the input grid
---------------------------

When ``in_grid`` is not specified, the input grid is inferred independently for each field:

1. First, the field's ``geography.grid_spec()`` metadata is used, when available.
2. If no grid spec can be produced this way, the field's latitude/longitude arrays
   (``geography.latlons()``) are used instead to build an unstructured lat/lon grid.

Each field is regridded differently depending on whether it is backed by an actual GRIB message.

GRIB fields
-----------

If a field has an associated GRIB message, GRIB-specific regridding is used so that, as far as
possible, the original GRIB metadata is preserved in the result (only the geometry-related keys
change). See :ref:`regrid-grib` for further details on how a GRIB message itself is regridded.

Whenever possible, the GRIB message is regridded directly by the backend, without going through an
intermediate NumPy array. This fast path is only used when:

- the **input** grid is structured. When the input grid is unstructured (for example, a set of
  points without a regular grid description), regridding falls back to extracting the field's values
  as a flat NumPy array, interpolating them with the backend, and re-encoding the result into a new
  GRIB message; and
- the backend implements GRIB message regridding. At present, only the MIR backend
  (:ref:`mir-backend`) provides this; with other backends, every GRIB field is regridded through the
  array-based fallback described above.

The **output** grid can be unstructured on either path, but the two paths handle it differently:

- On the GRIB-specific fast path, the backend re-encodes the interpolated values into a new GRIB
  message as usual (for MIR this changes the GRIB edition from 1 to 2), and the output grid spec is
  then additionally set on the resulting field.
- On the array-based fallback, the interpolated values currently cannot be re-encoded into a
  new GRIB message. In that case, instead of producing a new GRIB-backed field, the values and the
  output grid spec are set directly on the (copied) original field, which keeps the original GRIB
  message and metadata as-is.

.. note::

    In GRIB, only edition 2 can encode an arbitrary set of latitudes/longitudes as a grid, using the
    ecCodes GRIB key ``gridType="unstructured"``. Even then, the actual latitude/longitude values
    themselves cannot be encoded into the GRIB message. For this reason, earthkit-data represents
    such data as a field composed of the GRIB message together with a separate grid spec describing
    the custom grid, rather than embedding the coordinates in the message itself.

    This limitation does not apply to standard, named unstructured grids (for example ORCA grids):
    those are fully described by the GRIB message alone. It only affects an arbitrary, custom set of
    latitudes/longitudes.

Non-GRIB fields
----------------

If a field has no associated GRIB message, generic array-based regridding is used unconditionally:
the field's values are extracted as a flat NumPy array, interpolated by the backend, and a new field
is created from the interpolated values and the output grid spec, with the rest of the original
field's metadata (bar the grid) preserved. Since there is no GRIB message to preserve or re-encode,
this works the same way regardless of whether the input or output grid is structured or unstructured.

Examples
--------

- :ref:`/how-tos/mir/mir_healpix_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_octahedral_fieldlist.ipynb`
- :ref:`/how-tos/mir/mir_interpolation_types.ipynb`
- :ref:`/how-tos/precomputed/precomp_healpix_fieldlist.ipynb`
- :ref:`/how-tos/precomputed/precomp_octahedral_fieldlist.ipynb`

See also
--------

- :ref:`mir-regrid-high`
- :ref:`precomputed-regrid-high`
- :ref:`gridspec`
