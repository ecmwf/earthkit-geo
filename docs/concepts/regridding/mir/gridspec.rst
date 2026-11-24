.. _mir-gridspec:

MIR gridspecs
=============

A gridspec is the spec of a ``Grid`` object (:ref:`grid-object`): a dict of keys that describes a
grid completely (:ref:`grid-spec`). With the MIR backend, the input and output grids of a
regridding, ``in_grid`` and ``out_grid``, can be given as a gridspec: as a dict, as a string or
as a ``Grid`` object (``eckit.geo.Grid``), in :ref:`regrid <mir-regrid-high>` and
:ref:`regrid_array <mir-regrid-array>`. How specs are written — names, patterns and shortcuts — is
described in :ref:`grid-spec-names`.

Supported specs
---------------

The MIR backend accepts the grids described in the Grid documentation: the grid types of
:ref:`grid-types` and every spec key listed in :ref:`grid-reference`. The one exception is
spherical harmonics (:ref:`grid-sh`: the ``sh`` type, named ``T<T>``), which describes spectral
coefficients rather than points, and is not supported.

Besides the grid type, a spec can contain components that refine the grid. They are explained in
:ref:`grid-shared`:

- **area**: the extent of a regional grid, given as ``[north, west, south, east]`` in degrees
  (:ref:`grid-area`).
- **figure**: the shape and size of the reference surface, a sphere or an ellipsoid, such as
  ``wgs84`` (:ref:`grid-figure`).
- **projection**: how the points of the grid are mapped to longitude and latitude, for example a
  rotation or a map projection (:ref:`grid-projection`).
- **order**: how the points of the grid are numbered, such as the scanning mode of a regular grid
  (:ref:`grid-order`).

Examples by grid family
-----------------------

Each section below gives the accepted ``grid`` patterns and a few examples for a family of grids.
How these grids work is described in :ref:`grid-types`, and their spec keys in
:ref:`grid-ref-grids`.

Gaussian grids
~~~~~~~~~~~~~~

Three patterns for ``grid``, case-insensitive::

    1. OXXX (global octahedral reduced Gaussian grid)
    2. NXXX (global non-octahedral reduced Gaussian grid)
    3. FXXX (global regular Gaussian grid)

where *XXX* is the Gaussian number, the number of parallels between the pole and the equator. The
regular grid, ``FXXX``, has a constant number of points per parallel, while the reduced grids,
``NXXX`` and ``OXXX``, vary it to approach a constant local resolution. See :ref:`grid-gaussian`
for the details, including the pre-defined ``NXXX`` grids and the octahedral point counts, and
:ref:`grid-ref-grids` for the spec keys.

Examples:

.. code-block::

    {"grid": "O320"}
    {"grid": "N320"}
    {"grid": "F320"}


Regular latitude-longitude grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``grid`` format is::

    [DLON, DLAT]

where *DLON* and *DLAT* are the increments in degrees in longitudes and latitudes, respectively.
This grid is global, includes the origin (0°, 0°), and follows ECMWF's default scan mode. See
:ref:`grid-regular-ll` for the details, including the ``reference`` point, :ref:`grid-order` for
the scanning modes, and :ref:`grid-ref-grids` for the spec keys.

Examples:

.. code-block::

    {"grid": [1, 1]}
    {"grid": [0.1, 0.1]}


Arakawa C-grids
~~~~~~~~~~~~~~~

Two related ``type``-based gridspecs describe Arakawa C-grids:

- ``arakawa_c_um`` for Met Office Unified Model grids, which supports a range of UM *N* values,
  including::

      N48, N96, N144, N216, N320, N400, N512, N640, N768, N1280

- ``arakawa_c`` for general equirectangular projection-based grids, a compact way to define a
  C-grid on a regular latitude-longitude geometry

For both types, ``arrangement`` can be ``T``, ``U`` or ``V`` (case-insensitive); the default is
``T`` (:ref:`grid-staggering`). The two types map to progressively more foundational descriptions
(``arakawa_c_um`` -> ``arakawa_c`` with explicit UM-specific grid stretching and order,
``arakawa_c`` -> ``regular_ll`` with explicit ``grid`` and ``reference``), as shown by these
equivalent constructions (for UM N96, arrangement ``T``):

.. code-block::

    {"type": "arakawa_c_um", "N": 96}
    {"type": "arakawa_c", "N": 96, "grid_factor": [2, 1.3333333333], "order": "i+j+"}
    {"grid": [1.875, 1.25], "reference": [0.9375, 0.625], "order": "i+j+"}

See :ref:`grid-regular-ll` for the details and :ref:`grid-ref-grids` for the spec keys.


HEALPix grids
~~~~~~~~~~~~~

The ``grid`` is case-insensitive, in the formats::

    1. HXXX
    2. HRXXX
    3. HNXXX

where *XXX* is the HEALPix *Nside*. The second letter gives the ring (``R``, the default) or
nested (``N``) order of the points; the ``order`` key does the same. See :ref:`grid-healpix` for
the details and :ref:`grid-ref-grids` for the spec keys.

Examples (the last two are in nested order):

.. code-block::

    {"grid": "H512"}
    {"grid": "HR512"}
    {"grid": "H512", "order": "ring"}
    {"grid": "HN512"}
    {"grid": "H512", "order": "nested"}


ORCA grids
~~~~~~~~~~

These grids are the global ocean grids associated with the NEMO model, taken from the catalogue.
The ``grid`` format is case-insensitive, in the formats::

    1. ORCAXXX_[FTUVW]
    2. eORCAYYY_[FTUVW]

where the leading ``e`` of the second pattern stands for *extended*, indicating coverage closer to
the South pole. The suffix is required: *F*, *T*, *U*, *V* or *W*, the staggered position of the
points on NEMO's Arakawa C-grid (:ref:`grid-staggering`). The horizontal resolution numbers *XXX*
and *YYY* are written ``2``, ``1``, ``025`` and ``12``, approximately corresponding to 2°, 1°,
0.25° and 1/12°. See :ref:`grid-orca` for the details and :ref:`grid-ref-grids` for the spec keys.

Examples:

.. code-block::

    {"grid": "eORCA025_T"}
    {"grid": "ORCA1_F"}


.. _gridspec-unstructured:

Unstructured latitude-longitude grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These grids consist of an arbitrary set of latitudes and longitudes. They are often used to
represent observational data or model outputs on non-uniform grids. The ``grid`` in this case is
optional and can be omitted from the spec, which takes the following forms::

    {"grid": "unstructured", "latitudes": [LAT1, LAT2, LAT3, ...], "longitudes": [LON1, LON2, LON3, ...]}
    {"latitudes": [LAT1, LAT2, LAT3, ...], "longitudes": [LON1, LON2, LON3, ...]}

See :ref:`grid-unstructured` and :ref:`grid-as-unstructured` for the details, and
:ref:`regridding-unstructured` for how regridding falls back to unstructured grids when no other
specification is available.

Examples:

.. code-block::

    {"latitudes": [34.05, 36.12, 40.71], "longitudes": [-118.24, -115.17, -74.01]}


Rotated and projected grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A gridspec is not limited to a grid name: the keys of :ref:`grid-reference` can be combined. A
rotated regular latitude-longitude grid:

.. code-block::

    {"grid": [5, 5], "rotation": [-40, -20]}

where ``rotation`` is a shortcut for a rotation projection, and gives the latitude and then the
longitude of its south pole (:ref:`grid-projection`). Like ``area``, ``rotation`` is passed to MIR
as a job option for now, and works for GRIB input only (an array output cannot describe a rotated
grid, so array input raises ``ValueError``). A projected grid, named from the catalogue:

.. code-block::

    {"grid": "LAEA-EFAS-5km"}

See :ref:`grid-projected` for grids in a projection plane and :ref:`grid-ref-grids` for their spec
keys.


Related Grid concepts
---------------------

- :ref:`grid-spec` — what a spec is: its keys, values and rules
- :ref:`grid-spec-names` — the names, patterns and shortcuts a spec can use
- :ref:`grid-types` — the concepts specific to each family of grids
- :ref:`grid-reference` — the spec keys of every grid type
- :ref:`grid-area` — the ``area`` key, and cropping grids
- :ref:`grid-figure` — the shape and size of the reference surface
- :ref:`grid-projection` — the projection mapping a grid's points to longitude and latitude
- :ref:`grid-order` — how the points of a grid are numbered
