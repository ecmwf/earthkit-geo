.. _regridding-unstructured:

Unstructured grid fallback
============================

Regridding is driven by the :ref:`gridspec <gridspec>`: the input and output grids are each described
as a gridspec, and the selected backend interpolates between them accordingly.

Not all grids can currently be described by a gridspec. The gridspec, as an object/framework, is not
yet fully formulated, so for some grids there is no gridspec representation available at all. A common
example is data defined on a projection (for example, a Lambert Azimuthal Equal Area or a rotated pole
grid), described by 1D ``x``/``y`` coordinates together with 2D ``latitude``/``longitude`` coordinates.

For high-level data (such as GRIB fields or xarray objects), when the gridspec of the input data
cannot be specified or inferred, regridding automatically falls back to treating the input grid as
:ref:`unstructured <gridspec-unstructured>`. In earthkit-geo, "unstructured" means an arbitrary set of
latitude/longitude points, with no assumption about their layout, spacing or ordering — only the
per-point coordinates are used. Concretely, this means extracting the flat latitude/longitude arrays
for every point in the field and building a gridspec of the form::

    {"grid": "unstructured", "latitudes": [LAT1, LAT2, ...], "longitudes": [LON1, LON2, ...]}

This gridspec is then used as ``in_grid`` instead of a spec describing the original grid's actual
structure (which may not exist).

.. note::

    Interpolation methods other than nearest-neighbour (e.g. ``"linear"``) can be very slow on
    unstructured input grids, since no regular structure can be exploited by the backend to speed up
    neighbour search. Nearest-neighbour interpolation (``interpolation="nn"``) is recommended whenever
    the input ends up being treated as unstructured.

GRIB data
----------

For GRIB data, the gridspec can be inferred automatically from the data's own metadata. When this is
not possible, the reason is that a gridspec is not yet defined for that particular grid type — not that
anything is wrong with the data itself.

When using :func:`regrid` with earthkit-data :py:class:`~earthkit.data.core.fieldlist.FieldList`/
:py:class:`~earthkit.data.core.field.Field` objects, or with raw GRIB messages, the input grid is
determined as follows (see :ref:`regrid-fieldlist` and :ref:`regrid-grib-message` for full details):

1. the field's ``geography.grid_spec()`` metadata is used, when available;
2. otherwise, as a fallback, the field's latitude/longitude arrays (``geography.latlons()``) are used
   to build an unstructured gridspec.

.. code-block:: python

    import earthkit.data as ekd
    from earthkit.geo import regrid

    ds = ekd.from_source("file", "data.grib")
    field = ds[0]

    # None when the field's grid has no gridspec representation yet
    field.metadata().geography.grid_spec()

    # regridding still works: it falls back to an unstructured input grid automatically
    out = regrid(field, out_grid={"grid": [1, 1]}, interpolation="nn")

An explicitly provided ``in_grid`` always takes precedence over what can be inferred from the field,
so it can be used to bypass the fallback if the grid's gridspec is in fact known::

    out = regrid(field, in_grid={"grid": "O320"}, out_grid={"grid": [1, 1]})

Xarray data
------------

For xarray data, the gridspec can only be inferred automatically if the dataset carries the
earthkit-data ``"earthkit.grid_spec"`` attribute (set, for example, when the data was produced from
GRIB via earthkit-data's xarray engine).

- If this attribute is not present, but the gridspec is nonetheless known by other means, it can be
  passed explicitly via the ``in_grid`` keyword argument to :func:`regrid`.
- Otherwise — if the gridspec is not known, or not yet available for this kind of grid — regridding
  automatically falls back to treating the input as an unstructured grid, built from the dataset's
  latitude/longitude coordinate values, regardless of their original layout (1D, 2D, curvilinear, ...).

See :ref:`regrid-xarray` for the full grid- and dimension-inference rules.

Example: a field on a projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following loads a river discharge field defined on a Lambert Azimuthal Equal Area projection
(1D ``x``/``y`` plus 2D ``latitude``/``longitude`` coordinates). This kind of grid has no gridspec
representation yet, so the earthkit-data grid spec attribute is unavailable:

.. code-block:: python

    import earthkit.data as ekd
    from earthkit.geo import regrid

    ds = ekd.from_source("sample", "efas.nc").to_xarray()

    ds.earthkit.grid_spec  # None: no gridspec available for this grid

    # the target grid: a regular 0.1x0.1 degree latitude-longitude grid over a fixed area
    out_grid = {"grid": [0.1, 0.1], "area": [75, -30, 30, 40]}

    # since no gridspec can be inferred for the input, earthkit-geo falls back to treating
    # it as unstructured, extracting the per-point latitude/longitude values automatically
    r = regrid(ds.dis06, out_grid=out_grid, interpolation="nn")

Because the input had no gridspec attribute to begin with, the output does not get one either::

    r.earthkit.grid_spec  # still None

If the input grid's gridspec happens to be known even though it is not attached to the dataset, it can
be supplied explicitly to avoid the fallback::

    r = regrid(ds.dis06, in_grid={"grid": "O320"}, out_grid=out_grid)

See the full worked example in :ref:`/tutorials/mir_regrid_xarray_unstructured.ipynb`.

See also
--------

- :ref:`gridspec`
- :ref:`regrid-fieldlist`
- :ref:`regrid-grib-message`
- :ref:`regrid-xarray`
- :ref:`/tutorials/mir_regrid_xarray_unstructured.ipynb`
