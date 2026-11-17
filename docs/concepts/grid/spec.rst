.. _grid-spec:

The spec
========

A grid is created from a *spec*: a mapping from keys to values that describes it completely. This
page describes the spec of a ``Grid`` — the keys it uses, how names, patterns and shortcuts are
resolved, and the canonical form a grid writes back. :ref:`grid-reference` lists the keys of every
type.

A hierarchical description
--------------------------

A spec is a mapping from keys to values. A value is a scalar (string, integer, real or boolean), a
list of scalars, or a nested mapping, which is itself a spec:

.. code-block:: yaml

    grid: [0.1, 0.1]
    area: [72, -25, 25, 45]
    projection:
      type: rotation
      south_pole: [10, -40]
      figure: wgs84

A few rules apply throughout:

- **Keys are case-insensitive.** They are converted to lower case on input, and canonical specs use
  lower case. Values are case-sensitive, with two exceptions: grid names are matched
  case-insensitively (``O320`` and ``o320`` name the same grid), and so are some enumerated values,
  such as the staggering arrangement.
- **Keys belong to concepts and types.** Some keys apply to every object of a concept, whatever its
  type, and others only to a particular type (:ref:`grid-spec-type`). Nested mappings carry the specs
  of components (``projection: {...}``, ``figure: {...}``), each interpreted by its own concept.
- **Keys a builder does not use are ignored**, and do not appear in the canonical spec. There are
  deliberate exceptions where a key would be misleading, such as a false easting inside a projection
  (:ref:`grid-projection`).

.. _grid-spec-type:

The ``type`` key
----------------

Exactly one key is treated specially: ``type``. This holds for every spec, not only for grid specs:
every object built from a spec, whether a grid, a projection or an area, is built by the builder
that its ``type`` names. That builder receives the whole spec and interprets the other keys
according to its own rules.

The other keys of a spec are therefore of two kinds:

- **Keys common to a concept**, understood by every type of that concept: ``area``, ``projection``,
  ``figure`` and ``order`` for grids, or ``figure`` and ``reverse`` for projections.
- **Type-specific keys**, understood only by the type that defines them: ``pl`` for reduced Gaussian
  grids, ``Nside`` for HEALPix, and ``x`` and ``y`` for projected grids; ``lat_1`` and ``lat_2`` for
  Lambert conformal conic projections, ``south_pole`` for rotations, and ``source`` and ``target``
  for PROJ-based projections; ``latitudes`` and ``longitudes`` for polygons, and ``name_field`` for
  shapefile regions.

The same key can even mean different things to different types: ``grid: [1, 1]`` gives the increments
in degrees of a ``regular_ll`` grid, while ``grid: [2500, 2500]`` gives the increments in metres of a
``regular_xy`` grid. Nested specs follow the same rule at their own level, independently of the spec
around them. In the EFAS grid below, ``x``, ``y`` and ``false_xy`` are type-specific keys of the grid
(a ``regular_xy`` grid, implied by ``x`` and ``y``), and ``projection`` is a key common to grids. Its
value is a projection spec, in which ``type`` selects the projection builder, ``lat_0`` and ``lon_0``
are type-specific keys of ``laea``, and ``figure`` and ``reverse`` are keys common to projections:

.. code-block:: yaml

    x: [0, 4995000, 5000]
    y: [0, -4745000, -5000]
    false_xy: [1818500.021, -2287499.955]
    projection:
      type: laea
      lat_0: 52
      lon_0: 10
      figure: grs80
      reverse: true

Each concept has its own builders and its own set of types: grid types (``regular_ll``, ``HEALPix``,
...), projection types (``lcc``, ``rotation``, ``proj``, ...) and area types (``bounding_box``,
``polygon``, ...). Type names are matched exactly. For each concept, :ref:`grid-reference` lists the
common keys, and then the types with their aliases and type-specific keys. Figures, the simplest
component, need no ``type``: they are selected by name or by their dimensions.

Users rarely write ``type`` themselves, because it is implied when absent. For grids, names and
shortcuts imply it (:ref:`grid-spec-names`), and canonical specs omit it when it is implied:
``{grid: O320}`` is a ``reduced_gg`` grid, ``{grid: [1, 1]}`` is a ``regular_ll`` grid, and a spec
with ``x`` and ``y`` is a ``regular_xy`` grid. A projection spec without a ``type`` is the default
projection, ``eqc``, and an area spec without a ``type`` is a bounding box.

.. _grid-spec-names:

Names, patterns and shortcuts
-----------------------------

The value of the ``grid`` key can be a name. A name is either registered literally, in a catalogue
(:ref:`grid-catalogues`), or matches a registered pattern. Names are case-insensitive.

.. list-table::
    :header-rows: 1

    * - Name
      - Grid
      - Examples
    * - ``O<N>``
      - Octahedral reduced Gaussian grid (any N)
      - ``O320``, ``O1280``
    * - ``N<N>``
      - Classical reduced Gaussian grid
      - ``N320``, ``N640``
    * - ``F<N>``
      - Regular Gaussian grid
      - ``F320``
    * - ``<dlon>/<dlat>``
      - Regular latitude–longitude grid
      - ``1/1``, ``0.25/0.25``, ``1e-1/1e-1``
    * - ``H<Nside>``, ``HR<Nside>``, ``HN<Nside>`` (also ``H<Nside>R``, ``H<Nside>N``)
      - HEALPix grid, in ring (``R``, the default) or nested (``N``) order
      - ``H512``, ``HN512``
    * - ``T<T>``
      - Spherical harmonics, triangular truncation T
      - ``T1279``
    * - ``<name>_<arrangement>``
      - ORCA grids from the catalogue
      - ``eORCA025_T``, ``ORCA1_F``
    * - Catalogue names
      - Regional projected and hydrological grids
      - ``LAEA-EFAS-5km``, ``GLOFAS-4``
    * - A uid
      - The grid registered under that uid
      - ``d5bde4f52ff3a9bea5629cd9ac514410`` (``ORCA2_T``)

Some keys are also recognised as **shortcuts**, which imply other keys:

.. list-table::
    :header-rows: 1

    * - Shortcut
      - Implies
    * - ``N: <N>``
      - ``grid: O<N>``, the octahedral reduced Gaussian grid (unless ``type`` says otherwise)
    * - ``pl: [...]``
      - ``type: reduced_gg``
    * - ``grid: [dlon, dlat]``
      - ``type: regular_ll``
    * - ``latitudes: [...]`` and ``longitudes: [...]``
      - ``type: unstructured_ll``
    * - ``rotation: [lat, lon]``
      - ``projection: {type: rotation, ...}``, with the south pole at latitude ``lat`` and longitude
        ``lon``
    * - ``x: [...]`` and ``y: [...]``
      - ``type: regular_xy``
    * - ``type: <projection name>``
      - ``type: regular_xy`` with that projection (e.g. ``type: swisslv95``)
    * - ``uid: <uid>``
      - The registered grid with that uid

Note the order of the coordinates in ``rotation``: the latitude of the south pole comes first, as in
the rotation keyword of MARS requests. The explicit form,
``projection: {type: rotation, south_pole: [lon, lat]}``, gives the longitude first, as every other
point in a spec does.

.. _grid-spec-resolution:

How a spec is resolved
----------------------

Building a grid from a spec takes two steps. First the spec is *resolved*: names, uids and shortcuts
are expanded into the keys a builder needs. Then the ``type`` key selects the builder, which
constructs and validates the grid. Resolution proceeds as follows:

1. The catalogues are loaded, once per process (:ref:`grid-catalogues`).
2. If the spec gives a ``uid``, or a ``grid`` that is a uid, the spec registered under that uid is
   looked up. Since a uid identifies a grid exactly, the registered spec takes precedence over the
   keys given with it.
3. The shortcuts of :ref:`grid-spec-names` add the keys they imply.
4. If ``grid`` is a name, the spec registered for that name is added.
5. A spec with ``x`` and ``y`` but no ``type`` becomes a ``regular_xy`` grid.
6. A ``type`` that names a projection rather than a grid becomes a ``regular_xy`` grid with that
   projection.

Resolution never discards a key the user gave. It layers the specs contributed by names and shortcuts
behind the user's spec, so explicitly given keys take precedence over the keys a name or shortcut
contributes (uids being the exception, as noted). For example, in ``{grid: O320, area: [60, -10, 30,
40]}`` the name contributes ``{type: reduced_gg, N: 320, pl: [...]}`` and the user's area is kept.

.. _grid-spec-canonical:

The canonical spec
------------------

Once built, a grid writes its own spec, the *canonical spec*. It is minimal (defaults are omitted),
uses names where possible, normalises values, and sorts its keys. Two properties make it the form to
store and exchange:

- **Round trip.** Building a grid from its canonical spec gives an equal grid with the same points in
  the same order: ``Grid(grid.spec) == grid``. The test suites check this across grid types and
  orders.
- **Identity.** Equality, and the uid of grids generated from parameters, are defined by the
  canonical spec.

The table shows how the canonical spec normalises what users write (examples from the eckit tests):

.. list-table::
    :header-rows: 1

    * - Spec as given
      - Canonical spec
    * - ``{grid: 10/10}``
      - ``{"grid":[10,10]}``
    * - ``{pl: [20, 24, 24, 20]}``
      - ``{"grid":"O2"}``
    * - ``{grid: hr2}``
      - ``{"grid":"H2"}``
    * - ``{grid: h2n}``
      - ``{"grid":"H2","order":"nested"}``
    * - ``{pl: [20, 24, 28, 32, 32, 28, 24, 99]}``
      - ``{"grid":"N4","pl":[20,24,28,32,32,28,24,99]}``
    * - ``{grid: O96, area: [89.2842275325138, 0, -89.2842275325138, 359.1]}``
      - ``{"grid":"O96"}``
    * - ``{grid: [0.05, 0.05], area: 89.975/-179.975/-89.975/179.975}``
      - ``{"area":[90,-179.975,-90,180.025],"grid":[0.05,0.05],"reference":[0.025,0.025]}``
    * - ``{grid: [5, 5], rotation: [-40, -20]}``
      - ``{"grid":[5,5],"projection":{"south_pole":[-20,-40],"type":"rotation"}}``
    * - ``{grid: [5, 5], rotation: [-90, 0]}``
      - ``{"grid":[5,5]}``
    * - ``{grid: [30, 30], figure: earth}``
      - ``{"grid":[30,30]}``
    * - ``{grid: [30, 30], figure: wgs84}``
      - ``{"figure":"wgs84","grid":[30,30]}``
    * - ``{type: arakawa_c_um, N: 96}``
      - ``{"grid":[1.875,1.25],"order":"i+j+","reference":[0.9375,0.625]}``

A few of these deserve a comment. A ``pl`` that matches a known grid becomes its name, and an unknown
``pl`` is kept. An area that contains every point of a global grid is the global area, and
disappears. In the 0.05° example the area selects the points, and the grid passes through the area's
south-west corner. The canonical spec records the resulting offset of the grid as a *reference* point
near the origin, ``[0.025, 0.025]``, and the area the points actually span, here from pole to pole
and once around the globe. A rotation that does not rotate, and the default figure, are omitted.
Finally, a specialised type such as ``arakawa_c_um`` is written in the more fundamental terms of the
grid it builds.

.. _grid-catalogues:

Catalogues and named grids
--------------------------

Named grids, projections and areas are registered from YAML catalogues, distributed with eckit (in
``share/eckit/geo``). They include:

.. list-table::
    :header-rows: 1

    * - File
      - Contents
    * - ``grid.yaml``
      - Regional projected grids, such as the EFAS grid (``LAEA-EFAS-5km``) and the Swiss national
        grid (``swisslv95``)
    * - ``hydro.yaml``
      - Hydrological grids, such as GloFAS (``GLOFAS-4``)
    * - ``ORCA.yaml``
      - NEMO ORCA ocean grids, such as ``eORCA025_T``
    * - ``projection.yaml``
      - Named projections, such as ``swisslv95`` (EPSG:2056)
    * - ``area.yaml``, ``shapefile.yaml``
      - Areas, and libraries of region polygons (Eurostat GISCO country boundaries)

A grid catalogue lists entries under ``grid_names`` (by name) and ``grid_uids`` (by uid); YAML
anchors register one entry under several names and uids. An entry is an ordinary spec, which may
carry information beyond what identifies the grid:

.. code-block:: yaml

    grid_names:
      - ORCA2_T: &ORCA2_T
          type: ORCA
          name: ORCA2
          intgrid: O42
          arrangement: T
          shape: [149, 182]
          uid: "d5bde4f52ff3a9bea5629cd9ac514410"
          url: https://sites.ecmwf.int/repository/atlas/grids/orca/v1/ORCA2_T.atlas

    grid_uids:
      - "d5bde4f52ff3a9bea5629cd9ac514410": *ORCA2_T

The ``url`` of an entry points to its coordinates, which are retrieved when they are first needed.

The spec and the catalogue answer different questions. ``grid.spec`` says what identifies the grid;
``grid.catalog`` returns the catalogue entry that describes it (looked up by uid, then by name), with
whatever else is known about it, such as its shape or the URL of its coordinates. For a grid that is
not in a catalogue, ``grid.catalog`` is empty.

.. seealso::

    :ref:`grid-creating` for the forms a spec can take, and :ref:`grid-reference` for the specs of
    every grid type.
