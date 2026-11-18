.. _grid-ref-grids:

Grid specs
==========

The keys of every grid type: the keys common to all grids first, then one section per type with its aliases, names, shortcuts, type-specific keys and canonical form. :ref:`grid-types` explains how these grids work.

Keys common to grids
--------------------

Every grid understands these keys, whatever its type:

.. list-table::
   :widths: 12 28 22 38
   :header-rows: 1

   * - Key
     - Value
     - Default
     - Meaning
   * - ``type``
     - Grid type
     - Implied by names and shortcuts
     - The builder (:ref:`grid-spec-type`)
   * - ``grid``
     - A name, a uid, or a type-specific value
     - —
     - The grid (:ref:`grid-spec-names`)
   * - ``uid``
     - 32 hexadecimal digits
     - —
     - A registered grid (:ref:`grid-object`)
   * - ``area``
     - ``[north, west, south, east]``, or ``"n/w/s/e"``
     - Global
     - The extent of the grid (:ref:`grid-area`)
   * - ``order``
     - Order name
     - Type-specific
     - The numbering of the points (:ref:`grid-order`)
   * - ``projection``
     - Projection spec, or projection name
     - Default projection
     - :ref:`grid-projection`
   * - ``rotation``
     - ``[lat, lon]`` of the south pole
     - No rotation
     - Shortcut for a rotation (:ref:`grid-spec-names`)
   * - ``figure``
     - Figure spec
     - ``earth``
     - Shorthand for the default projection's figure (:ref:`grid-figure`)

Which keys a type uses is listed with each type below.

``regular_ll`` — regular latitude–longitude grid
------------------------------------------------

See :ref:`grid-regular-ll`.

Alias ``rotated_ll``. Names ``<dlon>/<dlat>``; shortcut ``grid: [dlon, dlat]``.

.. list-table::
   :widths: 25 40 35
   :header-rows: 1

   * - Key
     - Value
     - Default
   * - ``grid``
     - ``[dlon, dlat]``, positive
     - Required
   * - ``reference``
     - ``[lon, lat]``, a point of the grid
     - ``[0, 0]``, or the area's south-west corner
   * - ``area``, ``order``, ``projection``, ``rotation``, ``figure``
     -
     -

.. code-block:: python

    Grid(grid=[0.25, 0.25])
    Grid(grid="0.1/0.1", area=[72, -25, 25, 45])
    Grid(grid=[0.02, 0.02], area=[3.36, -6.82, -4.42, 4.8], order="i+j+", projection={"type": "rotation", "south_pole": [10, -43]})

Canonical form: ``{grid: [dlon, dlat]}``, with ``area``, ``reference``, ``order`` and ``projection`` when not the defaults. The area of a rotated grid is in the grid's own (rotated) frame.

``arakawa_c`` — Arakawa C-grid on a regular latitude–longitude grid
-------------------------------------------------------------------

See :ref:`grid-regular-ll`, and :ref:`grid-staggering` for the arrangements.

.. list-table::
   :widths: 25 40 35
   :header-rows: 1

   * - Key
     - Value
     - Default
   * - ``n``
     - Cells per 90°
     - Required
   * - ``arrangement``
     - ``T``, ``U`` or ``V``
     - ``T``
   * - ``grid_factor``
     - ``[fx, fy]``, positive scaling of the increments
     - ``[1, 1]``
   * - ``area``, ``order``, ``projection``, ``rotation``, ``figure``
     - As ``regular_ll``
     -

Canonical form: as the ``regular_ll`` grid it builds (:ref:`grid-regular-ll`). Catalogue example: ``GLOFAS-4``.

``arakawa_c_um`` — Met Office Unified Model grid
------------------------------------------------

See :ref:`grid-regular-ll`.

As ``arakawa_c``, with ``grid_factor: [2, 4/3]`` and ``order: i+j+``. Example: ``{type: arakawa_c_um, N: 96, arrangement: U}``.

``regular_gg`` — regular Gaussian grid
--------------------------------------

See :ref:`grid-gaussian`.

Alias ``rotated_gg``. Names ``F<N>``.

.. list-table::
   :widths: 25 40 35
   :header-rows: 1

   * - Key
     - Value
     - Default
   * - ``N``
     - Latitudes between pole and equator
     - Required (or the name)
   * - ``area``, ``order``, ``projection``, ``rotation``, ``figure``
     -
     -

Canonical form: ``{grid: F<N>}``.

``reduced_gg`` — reduced Gaussian grid
--------------------------------------

See :ref:`grid-gaussian`.

Alias ``reduced_rotated_gg``. Names ``O<N>`` (octahedral) and ``N<N>`` (classical); shortcuts ``pl: [...]`` and ``N: <N>``.

.. list-table::
   :widths: 25 40 35
   :header-rows: 1

   * - Key
     - Value
     - Default
   * - ``pl``
     - Points per latitude, 2N values from north to south
     - From the name
   * - ``area``, ``projection``, ``rotation``, ``figure``
     -
     -
   * - ``order``
     - ``i+j-``, ``i+j+``, ``i-j-`` or ``i-j+``
     - ``i+j-``

Canonical form: ``{grid: O<N>}`` or ``{grid: N<N>}``, with ``pl`` when neither octahedral nor classical.

``HEALPix`` — HEALPix grid
--------------------------

See :ref:`grid-healpix`.

Alias ``healpix``. Names ``H<Nside>``, ``HR<Nside>``, ``HN<Nside>`` (and ``H<Nside>R``, ``H<Nside>N``).

.. list-table::
   :widths: 25 40 35
   :header-rows: 1

   * - Key
     - Value
     - Default
   * - ``Nside``
     - Pixels per side of a base pixel
     - Required (or the name)
   * - ``order``
     - ``ring`` or ``nested`` (Nside a power of 2)
     - ``ring``
   * - ``projection``, ``rotation``, ``figure``
     -
     -

Canonical form: ``{grid: H<Nside>}``, with ``order: nested`` when nested.

``regular_xy`` — regular grid in a projection plane
---------------------------------------------------

See :ref:`grid-projected`.

Implied by ``x`` and ``y``, or by a ``type`` that names a projection.

.. list-table::
   :widths: 25 40 35
   :header-rows: 1

   * - Key
     - Value
     - Default
   * - ``x``, ``y``
     - ``[start, stop, step]``, stop included
     - Required (or ``grid`` and ``bounding_box_xy``)
   * - ``grid``
     - ``[dx, dy]``, positive
     - —
   * - ``bounding_box_xy``
     - ``[min_x, min_y, max_x, max_y]``
     - —
   * - ``false_xy``
     - ``[x0, y0]``, false easting and northing
     - ``[0, 0]``
   * - ``false_lonlat``
     - ``[lon, lat]`` of the grid origin
     - —
   * - ``projection``
     - From ``(x, y)`` to ``(lon, lat)``
     - None (plain Cartesian points)
   * - ``order``
     - Consistent with the directions of ``x`` and ``y``
     - From the directions

Canonical form: ``{x: [...], y: [...]}``, with ``false_xy``, ``projection`` and ``order`` when not the defaults. Catalogue examples: ``LAEA-EFAS-5km``, ``swisslv95``.

``unstructured_ll`` — unstructured grid
---------------------------------------

See :ref:`grid-unstructured`.

Alias ``unstructured``. Shortcut ``latitudes`` and ``longitudes``. Also the result of ``to_unstructured_ll()`` and ``cropped_unstructured_ll()`` (:ref:`grid-as-unstructured`).

.. list-table::
   :widths: 25 40 35
   :header-rows: 1

   * - Key
     - Value
     - Default
   * - ``latitudes``, ``longitudes``
     - Arrays of the same size
     - Required (or ``uid``)
   * - ``uid``
     - A registered unstructured grid
     - —
   * - ``name``
     - Registers the grid under this name
     - —

Canonical form: ``{type: unstructured_ll, uid: ...}``, or ``{grid: <name>}`` when named.

``ORCA`` — NEMO ocean grids
---------------------------

See :ref:`grid-orca`.

Names ``<name>_<arrangement>``, from ``ORCA2`` to ``eORCA12``, in arrangements ``F``, ``T``, ``U``, ``V`` and ``W``, e.g. ``eORCA025_T``.

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Key
     - Value
   * - ``name``, ``arrangement``, ``uid``
     - From the catalogue
   * - ``shape``, ``url``, ``intgrid``
     - Catalogue information

Canonical form: ``{grid: <name>_<arrangement>}``, e.g. ``{grid: eORCA1_T}``.

``sh`` — spherical harmonics
----------------------------

See :ref:`grid-sh`.

Names ``T<T>``.

.. list-table::
   :widths: 25 40 35
   :header-rows: 1

   * - Key
     - Value
     - Default
   * - ``truncation``
     - Triangular truncation, positive
     - Required (or the name)

Canonical form: ``{grid: T<T>}``.
