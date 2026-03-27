# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import logging
import re
from abc import abstractmethod

LOG = logging.getLogger(__name__)

HEALPIX_PATTERN = re.compile(r"[Hh]\d+")


class ShapeDoesNotMatchError(Exception):
    pass


# Wrapper around eckit.geo.Grid so that matrix inventory can be used.
# TODO: add logic here to support custom gridspecs that cannot be parsed by eckit.geo.
class GridSpec(dict):
    """Base class for grid specifications."""

    @staticmethod
    def from_dict(d):
        return EckitGridSpec(d)

    @staticmethod
    def from_any(d):
        if isinstance(d, GridSpec):
            return d
        elif isinstance(d, dict):
            return GridSpec.from_dict(d)

        return EckitGridSpec(d)

    @property
    @abstractmethod
    def grid(self):
        pass

    @property
    @abstractmethod
    def shape(self):
        pass

    def __eq__(self, o):
        if self.grid is not None and o.grid is not None:
            return self.grid.uid == o.grid.uid

        return False


class EckitGridSpec(GridSpec):
    """GridSpec that is parsed by eckit.geo.Grid."""

    def __init__(self, gs_in):
        from eckit.geo import Grid

        expected_shape = None
        if isinstance(gs_in, Grid):
            try:
                self._grid = gs_in
                spec = gs_in.spec
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
                self._grid = Grid(gs)
                spec = self._grid.spec
            except Exception as e:
                raise ValueError(f"Cannot parse gridspec with eckit.geo: {e}")

            if expected_shape is not None:
                shape = expected_shape
                if not isinstance(shape, (list, tuple)):
                    shape = [shape]
                if isinstance(shape, list):
                    shape = tuple(shape)
                if self._grid.shape != shape:
                    raise ShapeDoesNotMatchError(
                        f"Grid shape {self._grid.shape} does not match expected shape in {gs_in}"
                    )

        super().__init__(spec)

    def _patch(self, d):
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
    def grid(self):
        return self._grid

    @property
    def shape(self):
        return self._grid.shape

    @property
    def spec(self):
        return self._grid.spec


class CustomGridSpec(GridSpec):
    """Custom grid spec that cannot be parsed by eckit.geo.

    This is a placeholder for now, and can be extended in the future to suppor
    custom gridspecs.
    """

    @property
    def grid(self):
        return None

    @property
    def shape(self):
        return None

    @property
    def spec(self):
        return dict(self)
