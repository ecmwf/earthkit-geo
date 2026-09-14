# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""xarray backend for the regrid data-handler dispatch mechanism.

Exposes :data:`handler`, the :class:`~.handler.XarrayDataHandler` class used
to regrid ``xarray.Dataset``/``xarray.DataArray`` objects.
"""

from .handler import XarrayDataHandler

handler = XarrayDataHandler
