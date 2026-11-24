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
    def normalise_rotation(rotation):
        """Normalise a ``rotation`` value into the string format MIR expects.

        Parameters
        ----------
        rotation : str, list or tuple
            The south pole of the rotation, either already a string, or a
            2-element ``[south_pole_lat, south_pole_lon]`` sequence.

        Returns
        -------
        str
            ``rotation`` unchanged if it was already a string, otherwise its
            2 elements joined with ``"/"``.

        Raises
        ------
        ValueError
            If ``rotation`` is neither a string nor a 2-element list/tuple.
        """
        if isinstance(rotation, str):
            return rotation
        if isinstance(rotation, (list, tuple)):
            if len(rotation) == 2:
                return "/".join(map(str, rotation))

        raise ValueError(f"Invalid rotation format: {rotation}")

    @staticmethod
    def adjust_options(grid, kwargs):
        """Move temporary ``"area"`` and ``"rotation"`` grid-spec keys into the MIR job options.

        This is a temporary workaround for representing ``area`` and ``rotation``
        in a grid spec: MIR expects them as job options, not grid-spec keys.

        Parameters
        ----------
        grid : dict
            The (output) grid spec, possibly containing an ``"area"`` key.
        kwargs : dict
            The extra MIR job options to update.

        Returns
        -------
        Tuple[dict, dict]
            ``(grid, kwargs)`` unchanged if ``grid`` has neither key, otherwise
            copies of both with the keys removed from ``grid`` and added
            (normalised via :meth:`normalise_area` / :meth:`normalise_rotation`)
            to ``kwargs``.
        """
        # TODO: remove this once we have a better way to handle area and rotation in gridspec
        normalise = {"area": MirBackend.normalise_area, "rotation": MirBackend.normalise_rotation}
        keys = [k for k in normalise if k in grid]
        if keys:
            warn(
                f"The {'/'.join(keys)} key is a temporary workaround for {'/'.join(keys)} in gridspec",
                DeprecationWarning,
                stacklevel=2,
            )
            grid = grid.copy()
            kwargs = kwargs.copy()
            for k in keys:
                kwargs[k] = normalise[k](grid.pop(k))
        return grid, kwargs

    @staticmethod
    def get_grid_spec(grid):
        """Return the grid spec dict for ``grid``.

        Parameters
        ----------
        grid : eckit.geo.Grid|dict|str
            The grid, either as a ``Grid`` object or an already-plain spec.

        Returns
        -------
        dict | str
            ``grid.spec`` if ``grid`` is a ``Grid`` instance, otherwise
            ``grid`` unchanged. There is a special treatment for unstructured
            grids (see Notes). If ``grid`` is a dict or JSON string, it is
            returned as-is.


        Notes
        -----
        MIR (>= 1.30.2) accepts the grid specification either as a ``Grid`` object, as
        a plain dictionary or as a JSON string. In former versions the ``Grid`` object
        was not supported and the spec was extracted from it and passed to MIR as a
        plain dictionary. We keep this behaviour for backward compatibility, but will
        likely change it in future versions.

        Unstructured grids, at present, require a special treatment. For these grids, the
        most performant option is to pass the "uid" wrapped in a plain dictionary.
        This may be subject to change in future versions of MIR.
        """
        from eckit.geo import Grid

        if isinstance(grid, Grid):
            if grid.type in ("unstructured_ll", "unstructured"):
                grid = {"uid": grid.uid}
            else:
                return grid.spec
        return grid

    def regrid(
        self,
        data,
        in_grid,
        out_grid,
        interpolation="linear",
        rotation=None,
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
        rotation : str, list or tuple, optional
            The south pole of the output rotation, passed on by the field list
            handler. Not supported for arrays: raises ``ValueError``.

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

        if rotation is not None or (isinstance(out_grid, dict) and "rotation" in out_grid):
            raise ValueError("rotation in the output gridspec is only supported for GRIB input, not for arrays")

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
        rotation=None,
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
        rotation : str, list or tuple, optional
            The south pole of the output rotation, passed on by the field list
            handler. Added to the MIR job options like ``area``.

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
        if rotation is not None:
            kwargs["rotation"] = self.normalise_rotation(rotation)

        in_data = mir.GribMemoryInput(message)
        out = BytesIO()

        job = mir.Job()
        job.set("grid", out_grid)
        for k, v in kwargs.items():
            job.set(k, v)

        job.execute(in_data, out)

        return out.getvalue()


backend = MirBackend
