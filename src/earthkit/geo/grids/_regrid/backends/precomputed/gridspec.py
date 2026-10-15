# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Grid spec normalisation and comparison for the precomputed matrix inventory.

The matrix inventory (see :mod:`.backends.db`) needs to compare grid specs
for equality (to look up a precomputed interpolation matrix for a given
input/output grid pair) and to serialise them back to a canonical spec/string
form. :class:`_GridWrapper` and its subclasses wrap a grid spec (dict, JSON
string, or ``eckit.geo.Grid``) to provide this comparison/serialisation
uniformly, dispatching most grids to ``eckit.geo.Grid`` (:class:`_EckitGridWrapper`)
while handling a few grid families that ``eckit.geo.Grid`` cannot parse, or
parses inconsistently with the matrix inventory, as special cases
(:class:`_OrcaGridWrapper`, :class:`_CustomGridWrapper`).
"""

import json
import logging
import re

LOG = logging.getLogger(__name__)

HEALPIX_PATTERN = re.compile(r"[Hh]\d+")


class ShapeDoesNotMatchError(Exception):
    """Raised when a grid spec declares a ``"shape"`` that the actual grid does not have."""

    pass


# Temporary code to support gridspecs for the precomputed matrix inventory.


# Wrapper around eckit.geo.Grid so that matrix inventory can be used.
class _GridWrapper(dict):
    """Base class for grid spec wrappers used by the matrix inventory.

    A ``_GridWrapper`` is itself a dict holding the grid's canonical spec, plus
    ``type``/``grid``/``shape`` properties and an equality comparison
    suitable for matching grids in the matrix inventory index (see
    :class:`.backends.db.MatrixIndex`). Use :meth:`from_dict` or
    :meth:`from_any` rather than instantiating a subclass directly, as they
    pick the appropriate subclass for the given grid.
    """

    @staticmethod
    def from_dict(d):
        """Wrap a grid spec dict in the appropriate :class:`_GridWrapper` subclass.

        Parameters
        ----------
        d : dict
            The grid spec to wrap.

        Returns
        -------
        _GridWrapper
            :class:`_OrcaGridWrapper` for ORCA grids, :class:`_CustomGridWrapper`
            for a few other grid families ``eckit.geo.Grid`` cannot parse
            reliably, and :class:`_EckitGridWrapper` otherwise.
        """
        # in some environments eckit.geoGrid generate a lot of erros when
        # initializing with some specific gridspecs. We use dedicated classes
        # for these gridspecs to avoid the issue.
        if isinstance(d, dict):
            grid = d.get("grid", "")
            if isinstance(grid, str):
                grid = grid.upper()
                if "ORCA" in grid:
                    return _OrcaGridWrapper(d)
                elif "ICON" in grid or "CORE2" in grid or "NG5" in grid or "DART" in grid:
                    return _CustomGridWrapper(d)
        return _EckitGridWrapper(d)

    @staticmethod
    def from_any(d):
        """Wrap a grid spec, whatever form it is given in.

        Parameters
        ----------
        d : _GridWrapper, dict, str or Any
            The grid spec to wrap: an existing ``_GridWrapper`` (returned
            as-is), a dict, a JSON-encoded string, or anything
            ``eckit.geo.Grid`` itself can parse (e.g. a grid name).

        Returns
        -------
        _GridWrapper
            See :meth:`from_dict`.
        """
        if isinstance(d, _GridWrapper):
            return d
        elif isinstance(d, dict):
            return _GridWrapper.from_dict(d)
        elif isinstance(d, str):
            try:
                d = json.loads(d)
                return _GridWrapper.from_dict(d)
            except Exception:
                pass
        return _EckitGridWrapper(d)

    @property
    def grid_object(self):
        """eckit.geo.Grid or None: The underlying ``Grid`` object, when available."""
        return None

    @property
    def spec(self):
        """Dict or None: The grid spec as a dictionary, when available."""
        return None

    @property
    def type(self):
        """Str or None: The eckit-geo grid type, e.g. ``"regular-ll"``, ``"healpix"``."""
        return None

    @property
    def grid(self):
        """str, list or None: The eckit-geo grid, when available."""
        return None

    @property
    def shape(self):
        """Tuple[int, ...] or None: The grid's shape, when available."""
        return None

    @property
    def spec_str(self):
        """Str or None: The grid spec as a JSON-encoded string, when available."""
        return None

    def __eq__(self, value):
        return False


class _EckitGridWrapper(_GridWrapper):
    """Grid spec wrapper backed by ``eckit.geo.Grid``.

    Used for every grid ``eckit.geo.Grid`` can parse (i.e. everything except
    the special-cased grids handled by :class:`_OrcaGridWrapper` and
    :class:`_CustomGridWrapper`).
    """

    def __init__(self, gs_in):
        """Build the wrapper by parsing ``gs_in`` with ``eckit.geo.Grid``.

        Parameters
        ----------
        gs_in : eckit.geo.Grid, dict or Any
            The grid spec to parse: an existing ``Grid`` instance (used
            as-is), a dict (patched via :meth:`_patch` and checked against
            an optional ``"shape"`` entry), or anything else
            ``eckit.geo.Grid`` accepts directly (e.g. a grid name).

        Raises
        ------
        ValueError
            If ``eckit.geo.Grid`` cannot parse ``gs_in``.
        ShapeDoesNotMatchError
            If ``gs_in`` is a dict declaring a ``"shape"`` that does not
            match the shape of the parsed grid.
        """
        from eckit.geo import Grid

        expected_shape = None
        if isinstance(gs_in, Grid):
            try:
                self._grid_object = gs_in
                spec = self._grid_object.spec
            except Exception as e:
                raise ValueError(f"Cannot parse gridspec with eckit.geo: {e}")
        else:
            if isinstance(gs_in, dict):
                gs = gs_in.copy()
                self._patch(gs)
                expected_shape = gs.pop("shape", None)
                if gs.pop("global", None) in ("1", 1, True):
                    gs.pop("area", None)
            else:
                gs = gs_in

            try:
                self._grid_object = Grid(gs)
                spec = self._grid_object.spec
            except Exception as e:
                raise ValueError(f"Cannot parse gridspec with eckit.geo: {e}")

            if expected_shape is not None:
                shape = expected_shape
                if not isinstance(shape, (list, tuple)):
                    shape = [shape]
                if isinstance(shape, list):
                    shape = tuple(shape)
                if self._grid_object.shape != shape:
                    raise ShapeDoesNotMatchError(
                        f"Grid shape {self._grid_object.shape} does not match expected shape in {gs_in}"
                    )

        super().__init__(spec)

    def _patch(self, d):
        """Adjust a grid spec dict in place to work around known issues.

        Parameters
        ----------
        d : dict
            The grid spec dict to patch, modified in place.
        """
        grid = d.get("grid", "")
        if isinstance(grid, str):
            # the offical key is "order" in HEALPix, we still support "ordering" for compatibility
            if HEALPIX_PATTERN.match(grid):
                if "ordering" in d:
                    d["order"] = d.pop("ordering")
            # in the matrix inventory this orca grid has a 1D shape, but eckit.geo.Grid
            # generates a 2D shape for it, so we patch it here
            elif grid == "eORCA025_T":
                d["shape"] = (1442, 1207)

    @property
    def grid_object(self):
        """eckit.geo.Grid: The underlying ``Grid`` object."""
        return self._grid_object

    @property
    def spec(self):
        """dict: The grid's canonical spec, as returned by ``eckit.geo.Grid.spec``."""
        return self._grid_object.spec

    @property
    def type(self):
        """str: The eckit-geo grid type, e.g. ``"regular-ll"``, ``"healpix"``."""
        return self._grid_object.type

    @property
    def grid(self):
        """str, list, or None: The grid name or value, as returned by ``eckit.geo.Grid.grid``."""
        return self._grid_object.grid

    @property
    def shape(self):
        """Tuple[int, ...]: The grid's shape."""
        return self._grid_object.shape

    @property
    def spec_str(self):
        """str: The grid's canonical spec, as a string."""
        return self._grid_object.spec_str

    def __eq__(self, o):
        """Compare two grids by their ``eckit.geo.Grid`` uid.

        Parameters
        ----------
        o : _GridWrapper
            The grid wrapper to compare against.

        Returns
        -------
        bool
            True if both grids have a ``grid`` and their uids match, False
            otherwise.
        """
        if self._grid_object is not None and o.grid_object is not None:
            return self._grid_object.uid == o.grid_object.uid
        return False

    @property
    def inventory_docs_spec(self):
        """dict: The grid spec used for matching entries in the matrix inventory.

        Same as :attr:`spec` but with default values filled in for a couple
        of grid types (``order`` for HEALPix, ``octahedral`` for reduced
        Gaussian grids) and an added ``"_type"`` key.
        """
        # return a dict spec that can be used for inventory matching
        spec = dict(self.spec)
        if self.type in {"healpix", "HEALPix"} and "order" not in spec:
            spec["order"] = "ring"
        if self.type in {"reduced-gg", "reduced_gg"} and "octahedral" not in spec:
            spec["octahedral"] = True if self["grid"][0].lower() == "o" else False
        spec["_type"] = self.type
        return spec


class _OrcaGridWrapper(_GridWrapper):
    """Grid spec wrapper for ORCA grids in the precomputed matrix inventory.

    ORCA grids are handled separately from :class:`_EckitGridWrapper` because,
    in the matrix inventory, they are described with a 1D shape while
    ``eckit.geo.Grid`` reports a 2D shape for them.
    """

    def __init__(self, d):
        """Build the wrapper from a grid spec dict.

        Parameters
        ----------
        d : dict
            The grid spec to wrap. Only its ``"grid"`` and ``"shape"`` keys
            are used.
        """
        super().__init__(d)
        d = dict(d)
        self._patch(d)
        self._spec = {"grid": d.get("grid", None)}

        self._type = d.get("grid", None)
        self._shape = d.get("shape", None)

    @property
    def grid_object(self):
        """None: ORCA grids have no ``eckit.geo.Grid`` object."""
        return None

    @property
    def type(self):
        """Str or None: The grid name (e.g. ``"eORCA025_T"``), used as the type."""
        return self._type

    @property
    def grid(self):
        """Str or None: The grid name (e.g. ``"eORCA025_T"``)."""
        return self._type

    @property
    def shape(self):
        """Tuple[int, ...] or None: The grid's shape, as given in the input spec."""
        return self._shape

    @property
    def spec(self):
        """dict: The grid's canonical spec: ``{"grid": <name>}``."""
        return dict(self._spec)

    @property
    def spec_str(self):
        """str: The grid's canonical spec, JSON-encoded."""
        return json.dumps(self.spec)

    def _patch(self, d):
        """Adjust a grid spec dict in place for known ORCA grid shape mismatches.

        Parameters
        ----------
        d : dict
            The grid spec dict to patch, modified in place.
        """
        grid = d.get("grid", "")
        if isinstance(grid, str):
            # in the matrix inventory this orca grid has a 1D shape, but eckit.geo.Grid
            # generates a 2D shape for it, so we patch it here
            if grid == "eORCA025_T":
                d["shape"] = (1442, 1207)

    def __eq__(self, o):
        """Compare two ORCA grids by their grid name.

        Parameters
        ----------
        o : _GridWrapper
            The grid wrapper to compare against.

        Returns
        -------
        bool
            True if both grids have a non-None ``type`` (grid name) and they
            match, False otherwise.
        """
        if self._type is not None and o.type is not None:
            return self._type == o.type

        return False

    @property
    def inventory_docs_spec(self):
        """dict: :attr:`spec` with an added ``"_type": "ORCA"`` key."""
        # return a dict spec that can be used for inventory matching
        spec = dict(self.spec)
        spec["_type"] = "ORCA"
        return spec


class _CustomGridWrapper(_GridWrapper):
    """Grid spec wrapper for a few specific grids ``eckit.geo.Grid`` cannot parse.

    Used for grid names containing ``"ICON"``, ``"CORE2"``, ``"NG5"`` or
    ``"DART"`` (see :meth:`_GridWrapper.from_dict`). The spec is stored as-is,
    with no interpretation of the grid beyond wrapping it.
    """

    def __init__(self, d):
        """Build the wrapper from a grid spec dict.

        Parameters
        ----------
        d : dict
            The grid spec to wrap, stored as-is.
        """
        super().__init__(d)
        self._spec = dict(d)
        super().__init__(self._spec)

    @property
    def grid_object(self):
        """None: Custom grids have no ``eckit.geo.Grid`` object."""
        return None

    @property
    def type(self):
        """None: Custom grids have no known eckit-geo grid type."""
        return None

    @property
    def grid(self):
        """None: Custom grids have no ``eckit.geo.Grid`` grid name."""
        return None

    @property
    def shape(self):
        """None: Custom grids have no known shape."""
        return None

    @property
    def spec(self):
        """dict: The grid's spec, as originally given."""
        return dict(self._spec)

    @property
    def spec_str(self):
        """str: The grid's spec, JSON-encoded."""
        return json.dumps(self._spec)

    def __eq__(self, o):
        """Compare two custom grids by their spec dict.

        Parameters
        ----------
        o : Any
            The value to compare against.

        Returns
        -------
        bool
            True if ``o`` is a :class:`_CustomGridWrapper` with the same
            spec, False otherwise.
        """
        if isinstance(o, _CustomGridWrapper):
            return self._spec == o.spec
        return False

    @property
    def inventory_docs_spec(self):
        """dict: :attr:`spec` with an added ``"_type": None`` key."""
        spec = dict(self.spec)
        spec["_type"] = None
        return spec
