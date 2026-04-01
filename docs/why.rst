Why earthkit-geo?
==========================


**earthkit-geo** is a Python package providing geospatial computations such as regridding and rotation. It is part of the :xref:`earthkit` ecosystem.


Regridding
--------------

Regridding is available for inputs of Numpy arrays, earthkit-data GRIB :py:class:`~earthkit.data.core.field.Field` or :py:class:`~earthkit.data.core.fieldlist.FieldList` objects, and as an experimental feature for Xarray. By default it is using the :ref:`mir <mir-backend>` based on ECMWF's **MIR (Meteorological Interpolation and Regridding)** library.
