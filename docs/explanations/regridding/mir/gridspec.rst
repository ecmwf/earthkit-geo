.. _gridspec:

Gridspecs (for MIR)
====================================

A gridspec describes spatial grids in the form of a dict.

.. warning::

    The gridspec format is not finalised yet and may change in future releases. Area specification is currently supported only as a provisional feature (subject to change or removal) and may emit deprecation warnings.

The gridspecs supported when :ref:`regridding <mir-backend>` with the MIR backend are summarised below:


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
-------------

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

    ICON-GRID-NNNN-RXXBYY-Z

where *NNNN* is a sequential grid number, *rXXbYY* is the refinement level (root division *XX*,
bisection level *YY*), and *Z* indicates the grid type (``g`` global, ``n`` regional nest,
``l`` limited area, ``r`` radiation auxiliary). Many grids are rotated 36° around the z-axis.

*Global grids*:

- ``ICON-GRID-0002-R02B06-G``, 40 km
- ``ICON-GRID-0004-R02B07-G``, 20 km
- ``ICON-GRID-0006-R03B07-G``, 13 km
- ``ICON-GRID-0008-R02B05-G``, 80 km
- ``ICON-GRID-0010-R02B04-G``, 160 km
- ``ICON-GRID-0012-R02B04-G``, 160 km
- ``ICON-GRID-0014-R02B05-G``, 80 km
- ``ICON-GRID-0016-R02B06-G``, 40 km
- ``ICON-GRID-0018-R02B07-G``, 20 km
- ``ICON-GRID-0020-R03B06-G``, 26 km
- ``ICON-GRID-0022-R03B07-G``, 13 km
- ``ICON-GRID-0024-R02B06-G``, 40 km
- ``ICON-GRID-0026-R03B07-G``, 13 km
- ``ICON-GRID-0030-R02B05-G``, 80 km
- ``ICON-GRID-0033-R03B05-G``, 53 km
- ``ICON-GRID-0036-R03B06-G``, 26.5 km
- ``ICON-GRID-0039-R02B07-G``, 20 km
- ``ICON-GRID-0054-R02B04-G``, 160 km (idealized test grid)

*Regional grids* (Europe):

- ``ICON-GRID-0027-R03B08-N02``, 6.5 km
- ``ICON-GRID-0028-R02B07-N02``, 20 km
- ``ICON-GRID-0031-R02B06-N02``, 40 km
- ``ICON-GRID-0034-R03B06-N02``, 26.5 km
- ``ICON-GRID-0037-R03B07-N02``, 13 km
- ``ICON-GRID-0040-R02B08-N02``, 10 km
- ``ICON-GRID-0048-R03B07-LR``, 13 km
- ``ICON-GRID-0049-R03B08-L``, 6.5 km
- ``ICON-GRID-0050-R02B07-N02``, 20 km (Europe, North Atlantic and North Africa)
- ``ICON-GRID-0051-R03B07-N02``, 13 km (Europe, North Atlantic and North Africa)
- ``ICON-GRID-0052-R03B07-LR``, 13 km
- ``ICON-GRID-0053-R03B08-L``, 6.5 km

*Regional grids* (Germany):

- ``ICON-GRID-0041-R02B09-LR``, 5 km (COSMO-DE)
- ``ICON-GRID-0042-R02B10-L``, 2.5 km (COSMO-DE)
- ``ICON-GRID-0043-R19B06-LR``, 4 km (COSMO-D2)
- ``ICON-GRID-0044-R19B07-L``, 2 km (COSMO-D2)
- ``ICON-GRID-0045-R19B08-LN02``, 1 km (COSMO-D2)
- ``ICON-GRID-0046-R19B06-LR``, 4 km (enhanced COSMO-D2)
- ``ICON-GRID-0047-R19B07-L``, 2 km (enhanced COSMO-D2)

Max Planck Institute for Meteorology (ICON-MPIM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Global ocean-only grids for the ICON-O ocean model, rotated 37° around the z-axis and
symmetrical about the equator. Each grid covers ocean areas only:

- ``ICON-GRID-0036-R02B04-O``, 160 km
- ``ICON-GRID-0020-R02B05-O``, 80 km
- ``ICON-GRID-0046-R02B06-O``, 40 km
- ``ICON-GRID-0024-R02B07-O``, 20 km
- ``ICON-GRID-0016-R02B09-O``, 5 km
- ``ICON-GRID-0040-R02B10-O``, 2.5 km
- ``ICON-GRID-0038-R02B11-O``, 1.2 km

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

     ORCAXXX
     ORCAXXX_[FTUVW]
     eORCAYYY
     eORCAYYY_[FTUVW]

where the first letter stands for *extended* indicating coverage closer to the South pole. *T*, *U*, *V* or *W* define the point location respective to the supporting mesh elements, respectively in cell-centred, vertex and ``u``/``v`` edges point arrangements. The default arrangement is cell-centred or *F*.

Horizontal resolution numbers *XXX* and *YYY* are in increasing order *2* (*XXX* only), *1*, *025*, *12*, approximately corresponding to resolution in degrees 2°, 1°, 0.25° and 1/12°.

Examples:

.. code-block::

    {"grid": "eORCA025_T"}
    {"grid": "ORCA1"}
