.. _regrid-grib-message:

Regridding raw GRIB messages
=============================

:func:`regrid` also accepts a raw, already-encoded GRIB message as its ``data`` argument, given as a
:py:class:`bytes` object or an :class:`io.BytesIO`. This is a lower-level entry point than passing an
earthkit-data :py:class:`~earthkit.data.core.fieldlist.FieldList`/:py:class:`~earthkit.data.core.field.Field`
(see :ref:`regrid-fieldlist`): there is no field object involved here, only the encoded message itself.

.. code-block::

    from earthkit.geo import regrid

    # assuming "data.grib" contains a single GRIB message
    with open("data.grib", "rb") as f:
        out = regrid(f.read(), out_grid={"grid": [1, 1]})

``out`` is returned as a new, re-encoded GRIB message with the same type (``bytes`` or ``BytesIO``) as the
input.

Handling
--------

If the selected backend implements GRIB-specific message regridding (currently only the MIR backend,
:ref:`mir-backend`), the message is passed to it directly, and ``in_grid`` is not used since the backend
reads the input grid from the GRIB metadata itself. Otherwise (for example, with the precomputed-weights
backend, :ref:`precomputed-backend`), the message is transparently wrapped into an earthkit-data
``Field`` and regridded exactly as described in :ref:`regrid-fieldlist` — including its handling of
unstructured input/output grids — before being re-encoded back into a new GRIB message.

.. note::

    If the backend has no GRIB-specific regridding support and earthkit-data is not installed, regridding
    a raw GRIB message raises a ``ValueError``.

See also
--------

- :ref:`regrid-fieldlist`
- :ref:`mir-regrid-high`
- :ref:`gridspec`
