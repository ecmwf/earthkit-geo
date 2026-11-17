.. _grid-healpix:

HEALPix: Nside and ordering
===========================

HEALPix (Hierarchical Equal Area isoLatitude Pixelisation) divides the sphere into
12 base pixels, each divided into :math:`N_\mathrm{side} \times N_\mathrm{side}`
pixels of equal area. The grid points are the :math:`12 N_\mathrm{side}^2` pixel
centres, on :math:`4 N_\mathrm{side} - 1` rings of constant latitude. Equal areas
and isolatitude rings make HEALPix well suited to global analysis, spectral
transforms, and machine learning.

Ordering
--------

The points can be numbered in two ways:

- **``ring``** (the default): ring by ring from north to south, each from west to
  east;

- **``nested``**: hierarchically, following the subdivision of the base pixels, so
  that nearby pixels have nearby indices. This requires Nside to be a power of 2.

.. code-block:: python

    Grid("H512")                         # ring order
    Grid("HN512")                        # nested order
    Grid(grid="H512", order="nested")    # the same

The HEALPix grid is ECMWF's own implementation, independent of PROJ.

.. seealso::

    :ref:`grid-ref-grids` for the spec keys of the HEALPix grid, and :ref:`grid-order`
    for the orders of grids in general.
