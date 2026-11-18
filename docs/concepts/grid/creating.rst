.. _grid-creating:

Creating a Grid
===============

A ``Grid`` is created from a spec: :ref:`grid-spec` describes what a spec can contain, this page how
it is passed in, stored and read back.

Keyword arguments, a dict or a string
-------------------------------------

A grid can be created from keyword arguments, from a dict, or from a string in YAML (and so also
JSON, which is valid YAML). The four lines below create the same grid:

.. code-block:: python

    from earthkit.geo.grids import Grid

    a = Grid(grid="O320", area=[60, -10, 30, 40])                # keyword arguments
    b = Grid({"grid": "O320", "area": [60, -10, 30, 40]})        # a dict
    c = Grid("{grid: O320, area: [60, -10, 30, 40]}")            # a YAML string
    d = Grid('{"grid": "O320", "area": [60, -10, 30, 40]}')      # a JSON string

    assert a == b == c == d

All forms end up as the same structured spec and are interpreted by exactly the same code: keyword
arguments and dicts are serialised to a YAML flow string, which eckit parses as it would any other
string. Choose whichever is convenient: keyword arguments interactively, dicts in programs, strings in
configuration files, on the command line, or in storage.

A string that is not a mapping is taken as the value of ``grid``, so a grid can be created from its
name alone:

.. code-block:: python

    Grid("O320")          # same as Grid(grid="O320")
    Grid("eORCA025_T")
    Grid("H512")
    Grid("1/1")

A few rules apply:

- Pass either a spec or keyword arguments, not both.
- Values can be Python scalars, lists, tuples or dicts, NumPy arrays and scalars (converted to plain
  lists and numbers), or bytes (decoded as UTF-8).
- Unstructured grids given only by ``latitudes``, ``longitudes`` and optionally ``name`` are built
  directly from the arrays, without going through a string.

Round trips and persistence
---------------------------

The canonical spec is the natural way to store a grid, and a grid can always be rebuilt from it:

.. code-block:: python

    spec = grid.spec                    # a dict, e.g. stored as JSON alongside data
    assert Grid(spec) == grid
    assert Grid(grid.spec_str) == grid  # the same, from the JSON string

    same = Grid(uid=grid.uid)           # by uid, in the same process

Components on their own
-----------------------

Figures, projections and areas have their own specs, and their own classes in ``eckit.geo``, which
accept keyword arguments, dicts and mapping strings in the same way:

.. code-block:: python

    from eckit.geo import Area, Figure, Projection

    wgs84 = Figure(figure="wgs84")
    sphere = Figure(R=6371229)
    lcc = Projection(type="lcc", lat_1=63.3, lat_0=63.3, lon_0=15)
    lv95 = Projection("{type: proj, source: 'EPSG:2056'}")
    europe = Area(area=[72, -25, 25, 45])

The shorthand of a bare name applies to grids only; a figure is named with ``Figure(figure="wgs84")``.
The ``figure`` and ``projection`` properties of a grid return the grid's own objects, which keep the
grid alive as long as they are in use.

When a spec is not understood
-----------------------------

A spec that cannot be interpreted raises a ``RuntimeError``. Error messages state what was expected as
a spec fragment, for example ``RegularLL: {grid: [dlon, dlat]} expected`` or
``Rotation: {south_pole: [lon, lat]} expected``.

.. seealso::

    :ref:`grid-spec` for the keys a spec can hold and how they are resolved, and
    :ref:`grid-spec-canonical` for the spec a grid writes back.
