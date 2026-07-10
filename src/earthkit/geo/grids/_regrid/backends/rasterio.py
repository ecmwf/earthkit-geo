# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


import numpy as np

from . import Backend


def to_rasterio_kwargs(gs, nodata, prefix=""):
    return {f"{prefix}crs": gs.crs, f"{prefix}transform": gs.transform, f"{prefix}nodata": nodata}


def _resolution(shape, bounds, outer=True):
    # Based on: rioxarray.rioxarray.XRasterBase.resolution
    ny, nx = shape
    left, bottom, right, top = bounds
    dx = (right - left) / (nx - (0 if outer else 1))
    dy = (bottom - top) / (ny - (0 if outer else 1))
    return dx, dy


class CRSAffineGridSpec:
    """Affine transformation mapping into a coordinate reference system."""

    def __init__(self, crs, transform, shape):
        self.crs = crs
        self.transform = transform
        self.shape = shape  # array order (y, x)

    @classmethod
    def from_shape(cls, crs, bounds, shape):
        """Grid from a CRS, bounding box and shape.

        Parameters
        ----------
        crs : rasterio.crs.CRS
            Coordinate reference system.
        bounds : tuple[number, number, number, number]
            Bounding coordinates of the grid in CRS coordinates. Order is
            left, bottom, right, top. References grid boxes, grid points are
            generated at the center of the boxes.
        shape : tuple[int, int]
            Number of grid points in y and x directions (array order).
        """
        from affine import Affine

        left, _, _, top = bounds  # reference corner
        dx, dy = _resolution(shape, bounds, outer=True)
        affine = Affine.translation(left, top) * Affine.scale(dx, dy)
        return cls(crs, affine, shape)

    @classmethod
    def from_gridspec(cls, spec):
        from rasterio.crs import CRS

        if "area" in spec:
            north, west, south, east = spec["area"]
            bounds = (
                west,
                south,
                east,
                north,
            )
        else:
            bounds = (-180, -90, 180, 90)

        crs = CRS.from_epsg(4326)

        from earthkit.geo.grids import Grid

        shape = Grid(spec).shape

        return cls.from_shape(crs, bounds, shape)


class RasterioBackend(Backend):
    name = "rasterio"

    def regrid(self, data, in_grid, out_grid, interpolation="linear"):

        from rasterio.enums import Resampling
        from rasterio.warp import reproject

        if "area" not in out_grid and "area" in in_grid:
            out_grid["area"] = in_grid["area"]

        in_grid_obj = CRSAffineGridSpec.from_gridspec(in_grid)
        out_grid_obj = CRSAffineGridSpec.from_gridspec(out_grid)

        reproject_kwargs = {
            **to_rasterio_kwargs(in_grid_obj, np.nan, prefix="src_"),
            **to_rasterio_kwargs(out_grid_obj, np.nan, prefix="dst_"),
        }

        interpolation_key_map = {"linear": "bilinear", "nearest-neighbour": "nearest", "grid-box-average": "average"}

        if interpolation is not None:
            reproject_kwargs["resampling"] = Resampling[interpolation_key_map[interpolation]]

        # Flatten leading dimensions and restore after reprojection
        # (rasterio only handles 2- and 3-dimensional arrays)
        extra_shape = data.shape[:-2]
        flat_size = np.prod(extra_shape, dtype=int)
        data = data.reshape((flat_size, *in_grid_obj.shape[-2:]))
        regridded = np.empty((flat_size, *out_grid_obj.shape), dtype=data.dtype)
        reproject(source=data, destination=regridded, **reproject_kwargs)
        return regridded.reshape((*extra_shape, *out_grid_obj.shape)), out_grid


backend = RasterioBackend
