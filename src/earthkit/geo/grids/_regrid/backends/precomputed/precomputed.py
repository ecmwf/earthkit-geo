# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


from .. import Backend


class MatrixBackend(Backend):
    """Regrid backend interpolating with precomputed sparse matrices.

    Looks up a precomputed interpolation matrix for a given input/output
    grid pair and interpolation method in a :class:`~.db.MatrixDb` matrix
    inventory (the built-in system inventory by default, or a local
    directory/URL), and applies it as a sparse matrix-vector product.
    """

    name = "precomputed"
    system_inventory_id = "ecmwf"

    def __init__(self, inventory=None):
        """Initialise the backend with a matrix inventory.

        Parameters
        ----------
        inventory : str, optional
            Selects the matrix inventory to use: None or
            :attr:`system_inventory_id` for the built-in system inventory
            (downloaded/cached remotely), an ``http(s)://`` URL for a
            remote inventory, or a local directory path.
        """
        self.path_or_url = inventory
        self.db = self.get_db(inventory)

    def regrid(self, data, in_grid, out_grid, interpolation):
        """Interpolate ``data`` from ``in_grid`` onto ``out_grid`` using a precomputed matrix.

        Parameters
        ----------
        data : numpy.ndarray
            The values to interpolate, defined on ``in_grid``.
        in_grid : Any
            The input grid spec, in a form accepted by
            :class:`~.gridspec._GridWrapper.from_any`.
        out_grid : Any
            The output grid spec, in the same form as ``in_grid``.
        interpolation : str
            The interpolation method (e.g. ``"linear"``, ``"nearest-neighbour"``).

        Returns
        -------
        Tuple[numpy.ndarray, dict]
            The interpolated values, reshaped to the output grid's shape,
            and the output grid spec.

        Raises
        ------
        ValueError
            If no precomputed matrix is found for the given grid pair and
            interpolation method.
        """
        from .gridspec import _GridWrapper

        z, shape = self.db.find(in_grid, out_grid, interpolation)

        if z is None:
            raise ValueError(f"No precomputed weights found! {in_grid=} {out_grid=} {interpolation=}")

        # This should check for 1D (GG) and 2D (LL) matrices
        data = data.reshape(-1, 1)

        data = z @ data
        data = data.reshape(shape)

        _out_grid = _GridWrapper.from_any(out_grid)

        return data, _out_grid.spec

    def prepare_grid_object(self, grid_spec):
        """Wrap ``grid_spec`` for use in matrix inventory lookups.

        Parameters
        ----------
        grid_spec : Any
            The grid spec to wrap, in a form accepted by
            :class:`~.gridspec._GridWrapper.from_any`.

        Returns
        -------
        _GridWrapper or None
            The wrapped grid spec, or None if ``grid_spec`` is None.
        """
        from .gridspec import _GridWrapper

        if grid_spec is not None:
            return _GridWrapper.from_any(grid_spec)
        return None

    def get_db(self, path_or_url):
        """Get the :class:`~.db.MatrixDb` matrix inventory to use.

        Parameters
        ----------
        path_or_url : str or None
            None or :attr:`system_inventory_id` for the built-in system
            inventory, an ``http(s)://`` URL for a remote inventory, or a
            local directory path.

        Returns
        -------
        MatrixDb
            The corresponding matrix inventory.

        Raises
        ------
        ValueError
            If ``path_or_url`` is falsy and not None (e.g. an empty string).
        """
        if path_or_url is None or path_or_url == self.system_inventory_id:
            from .db import SYS_DB

            return SYS_DB
        elif path_or_url.startswith("http://") or path_or_url.startswith("https://"):
            from .db import MatrixDb

            return MatrixDb.from_url(path_or_url)
        elif path_or_url:
            from .db import MatrixDb

            return MatrixDb.from_path(path_or_url)
        else:
            raise ValueError(f"Invalid path_or_url={path_or_url} for backend={self.name}")
