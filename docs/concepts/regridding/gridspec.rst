.. _gridspec:

Gridspecs
=========

A gridspec is earthkit's way of describing a spatial grid. It takes the form of a dict.

.. warning::

    The gridspec format is not finalised yet and may change in future releases. Area specification is currently supported only as a provisional feature (subject to change or removal) and may emit deprecation warnings.

.. note::

    Not all gridspecs are supported for all backends. Please read :ref:`here <precomputed_inventory>` for information regarding support for the precomputed backend.

Summary of support
------------------

.. list-table:: Currently supported grids in earthkit-geo
    :widths: 10 30 10 10
    :header-rows: 1
    :stub-columns: 1

    * - Gridspec
      - Description
      - MIR Support
      - Precomputed Support
    * - OXXX
      - Global octahedral reduced Gaussian grid
      - ✔
      - ✔
    * - NXXX
      - Global non-octahedral reduced Gaussian grid
      - ✔
      - ✔
    * - FXXX
      - Regular Gaussian grid
      - ✔
      - ✘
    * - [DLON, DLAT]
      - Regular latitude-longitude grid
      - ✔
      - ✔
    * - arakawa_c
      - General equirectangular projection-based grids
      - ✔
      - ✘
    * - arakawa_c_um
      - Met Office Unified Model grids
      - ✔
      - ✘
    * - HXXX
      - HEALPix grid
      - ✔
      - ✔
    * - HRXXX
      - Ring HEALPix grid
      - ✔
      - ✘
    * - HNXXX
      - Nested HEALPix grid
      - ✔
      - ✘
    * - ORCAXXX_[FTUVW]
      - ORCA grid
      - ✔
      - ✘
    * - eORCAYYY_[FTUVW]
      - Extended ORCA grid
      - ✔
      - ✔

Read below for further information on any of these grid types.

Gaussian grids
--------------

Gaussian grids are a family of global grids defined by the number of constant-latitude lines (parallels) between the pole and equator.
Each parallel's latitude is at Gaussian quadrature points so integrals can be calculated over the sphere.
With increasing resolution these parallels approach but do not match the poles and equator.
Longitudinally there are two descriptions: regular with constant number of points per parallel, and reduced with varying number of points per parallel, increasingly higher near the equator so as to approach a constant local horizontal resolution over the globe.

There are three patterns for ``grid``, case-insensitive::

    1. OXXX (global octahedral reduced Gaussian grid)
    2. NXXX (global non-octahedral reduced Gaussian grid)
    3. FXXX (regular Gaussian grid)

where *XXX* is known as the Gaussian number, and represents the number of latitude lines between the pole and equator.
For details about these grids, see `here <https://confluence.ecmwf.int/display/FCST/Gaussian+grids>`_ and `here <https://confluence.ecmwf.int/spaces/FCST/pages/47300374/Introducing+the+octahedral+reduced+Gaussian+grid>`_.
The regular Gaussian grid, ``[Ff]XXX``, has a constant number of points per parallel, while the reduced Gaussian grid, ``[NnOo]XXX``, has a varying number of points per parallel.

There are a few pre-defined reduced Gaussian grids (eg. ``N160``, ``N320``, ``N640``) which have a pre-defined number of points per parallel.
The octahedral Gaussian grid (``[Oo]XXX``) has a specific number of points per parallel starting from 20 near the poles and increasing 4 points per parallel towards the equator, so it is defined for any Gaussian number.

Examples:

.. code-block::

    {"grid": "O320"}
    {"grid": "o320"}
    {"grid": "N320"}
    {"grid": "n320"}
    {"grid": "F320"}
    {"grid": "f320"}


Regular latitude-longitude grids
--------------------------------

The ``grid`` format is::

    [DLON, DLAT]

where *DLON* and *DLAT* are the increments in degrees in longitudes and latitudes, respectively. This grid is global and includes the origin (0°, 0°), and the ECMWF's default scan mode: first point at North-West with longitude the fastest index from West to East.

Example:

.. code-block::

    {"grid": [1, 1]}
    {"grid": [0.1, 0.1]}


Arakawa C-grids
---------------

Two related ``type``-based gridspecs describing Arakawa C-grids are supported:

- ``arakawa_c_um`` for Met Office Unified Model grids
- ``arakawa_c`` for general equirectangular projection-based grids

The ``arakawa_c_um`` type is specific to UM and supports a range of UM *N* values, including::

    N48, N96, N144, N216, N320, N400, N512, N640, N768, N1280

For both types, ``arrangement`` can be ``T``, ``U`` or ``V`` (case-insensitive). The default arrangement is ``T``.

The ``arakawa_c`` type provides a compact way to define an Arakawa C-grid on a regular latitude-longitude geometry. This is useful when you want explicit control of the grid construction while keeping the C-grid arrangement semantics.

The two types map to progressively more foundational descriptions:

- ``arakawa_c_um`` -> ``arakawa_c`` (with explicit UM-specific grid stretching and order)
- ``arakawa_c`` -> ``regular_ll`` (with explicit ``grid`` and ``reference``)

Example equivalent constructions (for UM N96, arrangement ``T``):

.. code-block::

    {"type": "arakawa_c_um", "N": 96}
    {"type": "arakawa_c", "N": 96, "grid_factor": [2, 1.3333333333], "order": "i+j+"}
    {"grid": [1.875, 1.25], "reference": [0.9375, 0.625], "order": "i+j+"}

You can then explore related grids by setting ``arrangement``.

References:

- `Met Office Unified Model <https://www.metoffice.gov.uk/research/approach/modelling-systems/unified-model>`_
- `Arakawa grids overview <https://en.wikipedia.org/wiki/Arakawa_grids>`_



HEALPix grids
-------------

The ``grid`` is case-insensitive, in the format::

    1. HXXX
    2. HRXXX
    3. HNXXX

The HEALPix grid is a global, hierarchical equal area isolatitude pixelisation of the sphere, where *XXX* is the HEALPix *Nside* representing the number of points per pixel side.
There are 12 base pixels, so the total number of points is :math:`12 N_\mathrm{side}^2`.
The points can be in ring (default) or nested order, indicated by the second letter or ``order`` (default ``ring``); if the second letter is ``[Nn]`` or ``order`` is ``nested``, the points order is nested. For details about this grid, see `here  <https://en.wikipedia.org/wiki/HEALPix>`_.

Examples (the last two are in nested order):

.. code-block::

    {"grid": "Hr512"}
    {"grid": "h512"}
    {"grid": "H512", "order": "ring"}
    {"grid": "Hn512"}
    {"grid": "h512", "order": "nested"}


ORCA grids
----------

These grids are designed for global ocean coverage, and are associated with the NEMO model. The ``grid`` format is case-insensitive, in the formats::

    1. ORCAXXX_[FTUVW]
    2. eORCAYYY_[FTUVW]

where the first letter stands for *extended* indicating coverage closer to the South pole.
The suffix is required: *F*, *T*, *U*, *V* or *W*. These letters define the point location respective to the supporting mesh elements, respectively in cell-centred, vertex and ``u``/``v`` edges point arrangements.
Use *F* for the cell-centred arrangement.

Horizontal resolution numbers *XXX* and *YYY* are in increasing order *2* (*XXX* only), *1*, *025*, *12*, approximately corresponding to resolution in degrees 2°, 1°, 0.25° and 1/12°.

Examples:

.. code-block::

    {"grid": "eORCA025_T"}
    {"grid": "ORCA1_F"}
