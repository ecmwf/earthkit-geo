# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""NumPy array backend for the regrid data-handler dispatch mechanism.

Exposes :data:`handler`, the :class:`NumpyDataHandler` class used by
:func:`earthkit.geo.grids._regrid.array.regrid` to regrid a plain NumPy
array.
"""

import logging

from .handler import DataHandler

LOG = logging.getLogger(__name__)


class NumpyDataHandler(DataHandler):
    """Regrid data handler for a plain NumPy ``ndarray``."""

    @staticmethod
    def match(data):
        """Return True if ``data`` is a NumPy ``ndarray``.

        Parameters
        ----------
        data : Any
            The data to check.

        Returns
        -------
        bool
            True if ``data`` is a NumPy ``ndarray``.

        Notes
        -----
        This handler is not registered in
        :data:`earthkit.geo.grids._regrid.data.DATA_HANDLERS` (only
        ``fieldlist``, ``grib`` and ``xarray`` are), so ``match`` is not
        exercised by :func:`~earthkit.geo.grids._regrid.data.get_data_handler`;
        :func:`earthkit.geo.grids._regrid.array.regrid` uses
        :class:`NumpyDataHandler` directly instead.
        """
        import earthkit.geo.grids._regrid.data.numpy as np

        return isinstance(data, np.ndarray)

    def regrid(self, data, in_grid=None, out_grid=None, **kwargs):
        """Regrid a NumPy array by delegating directly to the backend.

        Parameters
        ----------
        data : numpy.ndarray
            The values to regrid, defined on ``in_grid``.
        in_grid : Any, optional
            The input grid or grid spec.
        out_grid : Any, optional
            The output grid or grid spec.
        **kwargs : dict
            Must include ``backend`` (name or object of the regrid backend
            to use). Any other keyword arguments are forwarded to the
            backend's ``regrid`` method.

        Returns
        -------
        Tuple[numpy.ndarray, Any]
            The regridded values and the (possibly backend-adjusted) output
            grid spec.
        """
        backend = self.backend_from_kwargs(kwargs)
        return backend.regrid(data, in_grid, out_grid, **kwargs)


handler = NumpyDataHandler
