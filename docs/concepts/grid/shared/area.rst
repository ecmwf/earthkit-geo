.. _grid-area:

Area
====

An *area* is a region of the surface. Areas serve two purposes: they give the extent of regional grids, and they crop grids.

**The extent of a grid.** A grid's own area is always a bounding box, ``area: [north, west, south, east]`` in degrees, also accepted as a string ``area: north/west/south/east`` or as separate ``north``, ``west``, ``south`` and ``east`` keys. The default area is global, ``[90, 0, -90, 360]``, and is omitted from specs. For regular latitude–longitude and Gaussian grids, the area selects the points of the global grid inside it, as in the examples of :ref:`grid-spec-canonical`. Longitudes are normalised (west in ``[-180, 180)``, east in ``[west, west + 360]``), and a bounding box spanning 360° of longitude is periodic.

**Cropping.** Any area can crop any grid:

.. code-block:: python

    from eckit.geo import Area

    grid = Grid(grid=[1, 1])
    europe = Area(area=[72, -25, 25, 45])

    indices = grid.crop(europe)            # indices of the points inside, in the grid's order
    cropped = grid.cropped(europe)         # the cropped grid

The cropped grid is of the same type when possible, as for unrotated regular and reduced Gaussian grids cropped by a bounding box; otherwise it is an unstructured grid of the points inside (:ref:`grid-as-unstructured`).

**Area types.** Besides bounding boxes there are rectangles in a projection plane, polygons, and regions from shapefiles (:ref:`grid-ref-areas`):

.. code-block:: python

    Area(type="bounding_box_xy", bounding_box_xy=[0, 0, 1000000, 1000000])
    Area(type="polygon", latitudes=[0.5, 0.5, 10.5, 10.5], longitudes=[0.5, 10.5, 10.5, 0.5])
    Area(type="shapefile", library="cntr-rg-60m-2024", name="Portugal")

A polygon is a list of rings, each an unstructured grid of at least three points. Its canonical spec refers to the rings by their uids, ``{type: polygon, uid: ...}``, so a polygon is as compact to exchange as a grid. Shapefile regions are read from a file, a URL or a library registered in ``shapefile.yaml``, and selected by name or by index.

.. seealso::

    :ref:`grid-ref-areas` gives the full area specification reference, and :ref:`grid-as-unstructured` explains what cropping returns.
