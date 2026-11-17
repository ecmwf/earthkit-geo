.. _precomputed-gridspec:

Precomputed gridspecs
=====================

A gridspec is the spec of a ``Grid`` object (:ref:`grid-object`): a dict of
keys that describes a grid completely (:ref:`grid-spec`). With the
precomputed backend, the source and target grids of a regridding,
``in_grid`` and ``out_grid``, are each given as a gridspec.

Regridding only works for the source-target gridspec pairs listed in the
:ref:`weights inventory <precomputed_inventory>`, which is the authoritative
list of what is available. In addition:

- only global grids are supported
- for regular latitude-longitude grids, only the default scanning mode used
  at ECMWF is supported (values start at North-West and follow in
  consecutive rows from West to East)
- for ORCA grids, only the ``linear`` and ``nearest-neighbour``
  interpolation methods are supported

Supported gridspec forms
------------------------

The gridspecs in the :ref:`inventory <precomputed_inventory>` use the
following forms. The "Role" column says whether the form appears as a
source grid, as a target grid, or as both.

.. list-table:: Gridspec forms in the precomputed inventory
   :widths: 32 30 14 24
   :header-rows: 1

   * - Form
     - Example
     - Role
     - Grid documentation
   * - ``O<N>`` reduced Gaussian, octahedral
     - ``{"grid": "O320"}``
     - source, target
     - :ref:`grid-gaussian`
   * - ``N<N>`` reduced Gaussian, non-octahedral
     - ``{"grid": "N320"}``
     - source, target
     - :ref:`grid-gaussian`
   * - ``[dlon, dlat]`` regular latitude-longitude
     - ``{"grid": [0.25, 0.25]}``
     - source, target
     - :ref:`grid-regular-ll`
   * - ``H<N>`` HEALPix, ring order (the default)
     - ``{"grid": "H1024"}``
     - source
     - :ref:`grid-healpix`
   * - ``H<N>`` HEALPix, nested order
     - ``{"grid": "H1024", "order": "nested"}``
     - source
     - :ref:`grid-healpix`
   * - ORCA grid
     - ``{"grid": "eORCA025_T"}``
     - source
     - :ref:`grid-orca`

.. seealso::

    :ref:`grid-ref-grids` for the spec keys of these grids, and the
    :ref:`inventory pages <precomputed_inventory>` for the exact
    source-target combinations with available weights.
