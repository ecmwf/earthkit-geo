# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""MIR-backed regrid backend.

Provides :class:`MirBackend`, which delegates point-to-point interpolation
(of arrays and GRIB messages) to the MIR C++ library via its Python bindings.
"""

from warnings import warn

from . import Backend


class MirBackend(Backend):
    """Regrid backend delegating interpolation to the MIR library."""

    name = "mir"

    @staticmethod
    def normalise_area(area):
        """Normalise an ``area`` value into the string format MIR expects.

        Parameters
        ----------
        area : str, list or tuple
            The area, either already a string, or a 4-element
            ``[north, west, south, east]`` sequence.

        Returns
        -------
        str
            ``area`` unchanged if it was already a string, otherwise its
            4 elements joined with ``"/"``.

        Raises
        ------
        ValueError
            If ``area`` is neither a string nor a 4-element list/tuple.
        """
        if isinstance(area, str):
            return area
        if isinstance(area, (list, tuple)):
            if len(area) == 4:
                return "/".join(map(str, area))

        raise ValueError(f"Invalid area format: {area}")

    @staticmethod
    def adjust_options(grid, kwargs):
        """Move a temporary ``"area"`` grid-spec key into the MIR job options.

        This is a temporary workaround for representing ``area`` in a grid
        spec: MIR expects it as a job option, not a grid-spec key.

        Parameters
        ----------
        grid : dict
            The (output) grid spec, possibly containing an ``"area"`` key.
        kwargs : dict
            The extra MIR job options to update.

        Returns
        -------
        Tuple[dict, dict]
            ``(grid, kwargs)`` unchanged if ``grid`` has no ``"area"`` key,
            otherwise copies of both with ``"area"`` removed from ``grid``
            and added (normalised via :meth:`normalise_area`) to ``kwargs``.
        """
        # TODO: remove this once we have a better way to handle area in gridspec
        if "area" in grid:
            warn(
                "The area key is a temporary workaround for area in gridspec",
                DeprecationWarning,
                stacklevel=2,
            )
            grid = grid.copy()
            kwargs = kwargs.copy()
            area = grid.pop("area")
            kwargs["area"] = MirBackend.normalise_area(area)
        return grid, kwargs

    @staticmethod
    def get_grid_spec(grid):
        """Return the grid spec dict for ``grid``.

        Parameters
        ----------
        grid : eckit.geo.Grid or dict
            The grid, either as a ``Grid`` object or an already-plain spec.

        Returns
        -------
        dict
            ``grid.spec`` if ``grid`` is a ``Grid`` instance, otherwise
            ``grid`` unchanged.
        """
        from eckit.geo import Grid

        if isinstance(grid, Grid):
            return grid.spec
        return grid

    def regrid(
        self,
        data,
        in_grid,
        out_grid,
        interpolation="linear",
    ):
        """Interpolate an array from ``in_grid`` onto ``out_grid`` via MIR.

        Parameters
        ----------
        data : numpy.ndarray
            The values to interpolate, defined on ``in_grid``.
        in_grid : eckit.geo.Grid or dict
            The input grid spec.
        out_grid : eckit.geo.Grid or dict
            The output grid spec.
        interpolation : str, default="linear"
            The interpolation method (e.g. ``"linear"``, ``"grid-box-average"``,
            ``"nearest-neighbour"``).

        Returns
        -------
        Tuple[numpy.ndarray, dict]
            The interpolated values and the resulting output grid spec, as
            reported by MIR.
        """
        import mir
        import numpy as np

        kwargs = {
            "interpolation": interpolation,
        }

        in_grid = self.get_grid_spec(in_grid)
        out_grid = self.get_grid_spec(out_grid)

        out_grid, kwargs = self.adjust_options(out_grid, {})

        # mir.ArrayInput requires a C-contiguous array but apply_ufunc's dim
        # reordering can produce non-contiguous views. No-op if already contiguous.
        data = np.ascontiguousarray(data)

        input = mir.ArrayInput(data, in_grid)
        out = mir.ArrayOutput()

        job = mir.Job()
        job.set("grid", out_grid)
        job.set("interpolation", interpolation)  # NOTE: needs generalisation
        for k, v in kwargs.items():
            job.set(k, v)

        job.execute(input, out)

        return out.values(), out.spec

    # TODO: remove this once the gridspec can be written into the GRIB message
    def regrid_grib(
        self,
        message,
        out_grid,
        interpolation="linear",
    ):
        """Interpolate a GRIB message directly onto ``out_grid`` via MIR.

        Parameters
        ----------
        message : eccodes GRIB message
            The GRIB message to interpolate.
        out_grid : eckit.geo.Grid or dict
            The output grid spec.
        interpolation : str, default="linear"
            The interpolation method.

        Returns
        -------
        bytes
            The regridded output, encoded as a new GRIB message.
        """
        from io import BytesIO

        import mir

        out_grid = self.get_grid_spec(out_grid)

        kwargs = {
            "interpolation": interpolation,
        }

        # no 'automatic' necessary
        kremove = [k for k, v in kwargs.items() if v == "automatic"]
        for k in kremove:
            del kwargs[k]

        out_grid, kwargs = self.adjust_options(out_grid, kwargs)

        in_data = mir.GribMemoryInput(message)
        out = BytesIO()

        job = mir.Job()
        job.set("grid", out_grid)
        for k, v in kwargs.items():
            job.set(k, v)

        job.execute(in_data, out)

        return out.getvalue()


backend = MirBackend
