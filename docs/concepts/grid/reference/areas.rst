.. _grid-ref-areas:

Area specs
==========

The types of area, with their keys and canonical forms. :ref:`grid-area` explains how areas are used.

.. list-table::
   :widths: 15 55 30
   :header-rows: 1

   * - Type
     - Keys
     - Canonical form
   * - ``bounding_box`` (default)
     - ``area: [north, west, south, east]``, or ``area: "n/w/s/e"``, or ``north``, ``west``, ``south``, ``east`` (defaults global; ``west`` alone implies ``east = west + 360``, and vice versa)
     - ``{area: [n, w, s, e]}``
   * - ``bounding_box_xy``
     - ``bounding_box_xy: [min_x, min_y, max_x, max_y]``, in a projection plane
     - ``{type: bounding_box_xy, bounding_box_xy: [...]}``
   * - ``polygon``
     - ``latitudes`` and ``longitudes`` (one ring), or ``rings: [...]`` (each a dict of coordinates or a ``Grid``), or ``uid`` (one or a list, of unstructured grids of 3 points or more)
     - ``{type: polygon, uid: ...}``
   * - ``shapefile``
     - ``file`` (or ``shp``, with optional ``dbf``), or ``url``, or ``library``; ``name`` with ``name_field``, or ``index``
     - As given

Shapefile regions are an optional feature (``"shapefile" in eckit.geo.features()``). The libraries in ``shapefile.yaml`` are Eurostat GISCO country polygons (``cntr-rg-<scale>-<year>``), at scales from 1:1 million to 1:60 million.
