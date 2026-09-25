# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


"""Precomputed-matrix regrid backend.

Exposes :data:`backend`, the :class:`~.precomputed.MatrixBackend` class used
to interpolate values with precomputed sparse interpolation matrices from
the regrid matrix inventory (see :mod:`.db`).
"""

from .precomputed import MatrixBackend

backend = MatrixBackend
