.. _grid-sh:

Spherical harmonics: truncation
===============================

Spectral fields are described by a spherical harmonics "grid" (``sh``, named
``T<T>``), given by its triangular **truncation**. It has no points: its size is the
number of real coefficients, :math:`(T + 1)(T + 2)`. Including it in the same
framework lets spectral and grid-point fields be described, compared and exchanged
in the same way.

.. seealso::

    :ref:`grid-ref-grids` for the spec keys of ``sh``, and :ref:`grid-gaussian` for
    the Gaussian grids of spectral models.
