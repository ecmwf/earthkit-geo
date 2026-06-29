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


ORCA grids
----------

These grids are designed for global ocean coverage, and are associated with the NEMO model. The ``grid`` format is case-insensitive, in the formats::

    ORCAXXX_[FTUVW]
    eORCAYYY_[FTUVW]

where the first letter stands for *extended* indicating coverage closer to the South pole. The suffix is required: *F*, *T*, *U*, *V* or *W*. These letters define the point location respective to the supporting mesh elements, respectively in cell-centred, vertex and ``u``/``v`` edges point arrangements. Use *F* for the cell-centred arrangement.

Horizontal resolution numbers *XXX* and *YYY* are in increasing order *2* (*XXX* only), *1*, *025*, *12*, approximately corresponding to resolution in degrees 2°, 1°, 0.25° and 1/12°.

Examples:

.. code-block::

    {"grid": "eORCA025_T"}
    {"grid": "ORCA1_F"}
