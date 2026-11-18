.. _grid-gaussian:

Gaussian grids: N and pl
========================

Gaussian grids place their latitudes at the Gaussian quadrature points, the roots
of the Legendre polynomial of degree 2N, so that integrals over the sphere, and in
particular the spherical harmonic transforms of spectral models such as the IFS,
are computed accurately from the grid values. **N** is the number of latitudes
between a pole and the equator; the latitudes approach, but never reach, the poles
and the equator.

- A **regular Gaussian grid** (``F<N>``, type ``regular_gg``) has 4N points on
  every latitude.

- A **reduced Gaussian grid** (type ``reduced_gg``) has ``pl[j]`` points on
  latitude ``j``, fewer towards the poles, to approach a uniform resolution over
  the sphere. ``pl`` gives the numbers of points on the 2N latitudes, from north
  to south.

  - The **octahedral** grid (``O<N>``) has 20 points on the latitude nearest each
    pole and 4 more on each latitude towards the equator. It is defined for any N
    and has :math:`4N^2 + 36N` points.

  - The **classical** grids (``N<N>``) use the tabulated ``pl`` of ECMWF's
    operational history, known for N = 16, 24, 32, 48, 64, 80, 96, 128, 160, 200,
    256, 320, 400, 512, 576, 640, 800, 1024, 1280, 1600, 2000, 4000 and 8000.

  - Any other ``pl`` is accepted, and kept in the canonical spec.

Gaussian grids can be regional (``area``), rotated (``rotation`` or
``projection``), and reordered.

.. seealso::

    :ref:`grid-ref-grids` for the spec keys of the Gaussian grid types.
