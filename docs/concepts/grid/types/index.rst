.. _grid-types:

Grid types
==========

Some concepts belong to grids in general: see :ref:`grid-shared`. Others only make
sense for a particular family of grids, and are more specialised concepts of
``Grid``. The pages below explain them; :ref:`grid-ref-grids` lists the keys.

.. list-table:: Grid types and the page that describes them
   :widths: 30 35 35
   :header-rows: 1

   * - Page
     - Grid types
     - Concepts
   * - :ref:`grid-gaussian`
     - ``regular_gg``, ``reduced_gg``
     - Gaussian ``N``, the latitude point counts ``pl``
   * - :ref:`grid-regular-ll`
     - ``regular_ll``, ``arakawa_c``, ``arakawa_c_um``
     - increments, reference point, staggering and ``arrangement``
   * - :ref:`grid-projected`
     - ``regular_xy``
     - ``x`` and ``y`` ranges, the false origin
   * - :ref:`grid-healpix`
     - ``HEALPix``
     - ``Nside``, ring and nested ordering
   * - :ref:`grid-unstructured`
     - ``unstructured_ll``, ORCA
     - grids from coordinate arrays, any grid as an unstructured grid
   * - :ref:`grid-sh`
     - ``sh``
     - spherical harmonics truncation

.. toctree::
   :hidden:
   :maxdepth: 2

   gaussian
   regular_ll
   projected
   healpix
   coordinates
   spectral
