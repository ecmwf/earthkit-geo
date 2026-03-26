.. _gridspec:

Gridspecs (for MIR)
====================================

A gridspec describes spatial grids in the form of a dict.

.. warning::

    The gridspec format is not finalised yet and may change in future releases. Area specification is currently supported only as a provisional feature (subject to change or removal) and may emit deprecation warnings.

The gridspecs supported by the ``in_grid`` and ``out_grid`` options in :ref:`regrid() <mir-regrid>` with the (default) MIR backend are summarised below:


Gaussian grids
--------------

Gaussian grids are a family of global grids defined by the number of constant-latitude lines (parallels) between the pole and equator. Each parallel's latitude is at Gaussian quadrature points so integrals can be calculated over the sphere. With increasing resolution these parallels approach but do not match the poles and equator. Longitudinally there are two descriptions: regular with constant number of points per parallel, and reduced with varying number of points per parallel, increasingly higher near the equator so as to approach a constant local horizontal resolution over the globe.

There are three patterns for ``grid``, case-insensitive::

    OXXX
    NXXX
    FXXX

where *XXX* is known as the Gaussian number, and represents the number of latitude lines between the pole and equator. For details about these grids, see `here <https://confluence.ecmwf.int/display/FCST/Gaussian+grids>`_. The regular Gaussian grid, ``[Ff]XXX``, has a constant number of points per parallel, while the reduced Gaussian grid, ``[NnOo]XXX``, has a varying number of points per parallel.

There are a few pre-defined reduced Gaussian grids (eg. ``N160``, ``N320``, ``N640``) which have a pre-defined number of points per parallel. The octahedral Gaussian grid (``[Oo]XXX``) has a specific number of points per parallel starting from 20 near the poles and increasing 4 points per parallel towards the equator, so it is defined for any Gaussian number.

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



HEALPix grids
-------------

The ``grid`` is case-insensitive, in the format::

    HXXX
    HRXXX
    HNXXX

The HEALPix grid is a global, hierarchical equal area isolatitude pixelisation of the sphere, where *XXX* is the HEALPix *Nside* representing the number of points per pixel side. There are 12 base pixels, so the total number of points is :math:`12 N_\mathrm{side}^2`. The points can be in ring (default) or nested order, indicated by the second letter or ``order`` (default ``ring``); if the second letter is ``[Nn]`` or ``order`` is ``nested``, the points order is nested. For details about this grid, see `here  <https://en.wikipedia.org/wiki/HEALPix>`_.

Examples (the last two are in nested order):

.. code-block::

    {"grid": "Hr512"}
    {"grid": "h512"}
    {"grid": "H512", "order": "ring"}
    {"grid": "Hn512"}
    {"grid": "h512", "order": "nested"}


FESOM grids
----------

These grids are designed for global ocean coverage, and are associated with the FESOM model. The ``grid`` format is case-insensitive, in the formats::

     CORE2
     CORE2_[CN]
     DART
     DART_[CN]
     NG5
     NG5_[CN]

where *CORE2*, *DART* and *NG5* are distinct grids in increasing horizontal resolution, up to approximately 5 km. *C* and *N* define the point location respective to the supporting mesh elements, respectively in cell-centred and vertex point arrangements. The default arrangement is cell-centred or *C*.

Examples:

.. code-block::

    {"grid": "core2_n"}
    {"grid": "NG5"}


ICON grids
----------

ICON is a model used for both atmospheric and ocean modelling. ICON grids are identified by name (case-insensitive). Three
sets of pre-defined grids are available:

- **MeteoSwiss** (ICON-CH) — regional grids covering the Alpine region
- **Deutscher Wetterdienst** (ICON-DWD) — global and regional grids for operational NWP
- **Max Planck Institute for Meteorology** (ICON-MPIM) — global ocean-only grids for the ICON-O ocean model

MeteoSwiss (ICON-CH)
~~~~~~~~~~~~~~~~~~~~

Regional operational grids covering the Alpine region. The version suffix ``-V1`` is optional, version 1 the default. The grid names are case-insensitive:

- ``ICON-CH1``, ``ICON-CH1-V1``, 1° resolution
- ``ICON-CH2``, ``ICON-CH2-V1``, 2° resolution

Deutscher Wetterdienst (ICON-DWD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DWD grids follow the naming convention::

    icon-grid-NNNN-rXXbYY-Z

where *NNNN* is a sequential grid number, *rXXbYY* is the refinement level (root division *XX*,
bisection level *YY*), and *Z* indicates the grid type (``g`` global, ``n`` regional nest,
``l`` limited area, ``r`` radiation auxiliary). Many grids are rotated 36° around the z-axis.

*Global grids*:

- ``icon-grid-0002-r02b06-g``, 40 km
- ``icon-grid-0004-r02b07-g``, 20 km
- ``icon-grid-0006-r03b07-g``, 13 km
- ``icon-grid-0008-r02b05-g``, 80 km
- ``icon-grid-0010-r02b04-g``, 160 km
- ``icon-grid-0012-r02b04-g``, 160 km
- ``icon-grid-0014-r02b05-g``, 80 km
- ``icon-grid-0016-r02b06-g``, 40 km
- ``icon-grid-0018-r02b07-g``, 20 km
- ``icon-grid-0020-r03b06-g``, 26 km
- ``icon-grid-0022-r03b07-g``, 13 km
- ``icon-grid-0024-r02b06-g``, 40 km
- ``icon-grid-0026-r03b07-g``, 13 km
- ``icon-grid-0030-r02b05-g``, 80 km
- ``icon-grid-0033-r03b05-g``, 53 km
- ``icon-grid-0036-r03b06-g``, 26.5 km
- ``icon-grid-0039-r02b07-g``, 20 km
- ``icon-grid-0054-r02b04-g``, 160 km (idealized test grid)

*Regional grids* (Europe):

- ``icon-grid-0027-r03b08-n02``, 6.5 km
- ``icon-grid-0028-r02b07-n02``, 20 km
- ``icon-grid-0031-r02b06-n02``, 40 km
- ``icon-grid-0034-r03b06-n02``, 26.5 km
- ``icon-grid-0037-r03b07-n02``, 13 km
- ``icon-grid-0040-r02b08-n02``, 10 km
- ``icon-grid-0048-r03b07-lr``, 13 km
- ``icon-grid-0049-r03b08-l``, 6.5 km
- ``icon-grid-0050-r02b07-n02``, 20 km (Europe, North Atlantic and North Africa)
- ``icon-grid-0051-r03b07-n02``, 13 km (Europe, North Atlantic and North Africa)
- ``icon-grid-0052-r03b07-lr``, 13 km
- ``icon-grid-0053-r03b08-l``, 6.5 km

*Regional grids* (Germany):

- ``icon-grid-0041-r02b09-lr``, 5 km (COSMO-DE)
- ``icon-grid-0042-r02b10-l``, 2.5 km (COSMO-DE)
- ``icon-grid-0043-r19b06-lr``, 4 km (COSMO-D2)
- ``icon-grid-0044-r19b07-l``, 2 km (COSMO-D2)
- ``icon-grid-0045-r19b08-ln02``, 1 km (COSMO-D2)
- ``icon-grid-0046-r19b06-lr``, 4 km (enhanced COSMO-D2)
- ``icon-grid-0047-r19b07-l``, 2 km (enhanced COSMO-D2)

Max Planck Institute for Meteorology (ICON-MPIM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Global ocean-only grids for the ICON-O ocean model, rotated 37° around the z-axis and
symmetrical about the equator. Each grid covers ocean areas only:

- ``icon-grid-0036-r02b04-o``, 160 km
- ``icon-grid-0020-r02b05-o``, 80 km
- ``icon-grid-0046-r02b06-o``, 40 km
- ``icon-grid-0024-r02b07-o``, 20 km
- ``icon-grid-0016-r02b09-o``, 5 km
- ``icon-grid-0040-r02b10-o``, 2.5 km
- ``icon-grid-0038-r02b11-o``, 1.2 km

Examples:

.. code-block::

   {"grid": "icon-ch2"}
   {"grid": "ICON-CH1-V1"}
   {"grid": "icon-grid-0022-r03b07-g"}
   {"grid": "icon-grid-0044-r19b07-l"}
   {"grid": "icon-grid-0024-r02b07-o"}


ORCA grids
----------

These grids are designed for global ocean coverage, and are associated with the NEMO model. The ``grid`` format is case-insensitive, in the formats::

     ORCAXXX(|_[FTUVW])
     eORCAYYY(_|[FTUVW])

where the first letter stands for *extended* indicating coverage closer to the South pole. *T*, *U*, *V* or *W* define the point location respective to the supporting mesh elements, respectively in cell-centred, vertex and ``u``/``v`` edges point arrangements. The default arrangement is cell-centred or *F*.

Horizontal resolution numbers *XXX* and *YYY* are in increasing order *2* (*XXX* only), *1*, *025*, *12*, approximately corresponding to resolution in degrees 2°, 1°, 0.25° and 1/12°.

Examples:

.. code-block::

    {"grid": "eORCA025_T"}
    {"grid": "ORCA1"}
