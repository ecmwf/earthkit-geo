# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import logging

from .handler import DataHandler

LOG = logging.getLogger(__name__)


class FieldListDataHandler(DataHandler):
    @staticmethod
    def match(data):
        from earthkit.geo.utils import is_module_loaded

        if not is_module_loaded("earthkit.data"):
            return False

        try:
            import earthkit.data

            return isinstance(data, earthkit.data.FieldList) and hasattr(data, "get")
        except Exception:
            return False

    @staticmethod
    def input_gridspec(field, index):
        try:
            in_grid = field.geography.grid_spec()
        except Exception as e:
            LOG.exception(f"Cannot get input grid_spec from metadata for field[{index}]: {e}")
            raise
        if in_grid is None:
            raise ValueError(f"Cannot get input grid_spec from metadata for field[{index}]")

        return in_grid

    def regrid(self, data, **kwargs):
        backend = self.backend_from_kwargs(kwargs)

        grid = kwargs.pop("grid", None)
        if grid is None:
            raise ValueError("Missing 'grid' argument")

        if hasattr(backend, "regrid_grib"):
            # TODO: remove this when ecCodes supports setting the gridSpec on a GRIB handle
            return self._regrid_grib(data, backend, grid, **kwargs)
        else:
            return self._regrid_array(data, backend, grid, **kwargs)

    def _regrid_array(self, data, backend, grid, **kwargs):
        import earthkit.data

        ds = data
        assert grid is not None

        # TODO: refactor this when this limitation is removed
        from earthkit.geo.regrid.gridspec import GridSpec

        out_grid = GridSpec.from_dict(grid)
        if not out_grid.is_regular_ll():
            raise ValueError(
                (f"Fieldlists can only be regridded to global regular lat-lon target grids. Target grid is {out_grid}")
            )

        fields = []
        for i, f in enumerate(ds):
            vv = f.to_numpy(flatten=True)

            in_grid = self.input_gridspec(f, i)

            v_res, out_grid = backend.regrid(
                vv,
                in_grid,
                out_grid,
                **kwargs,
            )

            fields.append(f.set({"values": v_res, "geography.grid_spec": out_grid}))

        r = earthkit.data.create_fieldlist(fields)

        return r

    def _regrid_grib(self, values, backend, grid, **kwargs):
        # TODO: remove this when ecCodes supports setting the gridSpec on a GRIB handle
        from earthkit.data.encoders.grib import GribEncoder
        from earthkit.data.field.grib.create import create_grib_field_from_buffer

        assert hasattr(backend, "regrid_grib")

        ds = values
        assert grid is not None

        encoder = GribEncoder()

        r = []
        for i, f in enumerate(ds):
            if f._get_grib():
                d = encoder.encode(f)
                message = d.to_bytes()
                if message is None:
                    raise ValueError(f"field at index={i} type={type(f)} is not supported in regrid!")

            v_res = backend.regrid_grib(message, grid, **kwargs)
            r.append(create_grib_field_from_buffer(v_res))

        return ds.from_fields(r)


class FieldDataHandler(DataHandler):
    @staticmethod
    def match(data):
        from earthkit.geo.utils import is_module_loaded

        if not is_module_loaded("earthkit.data"):
            return False

        try:
            import earthkit.data

            return isinstance(data, earthkit.data.Field) and hasattr(data, "get")
        except Exception:
            return False

    def regrid(self, data, **kwargs):
        from earthkit.data import FieldList

        ds = FieldList.from_fields([data])
        return FieldListDataHandler().regrid(ds, **kwargs)[0]


handler = [FieldListDataHandler, FieldDataHandler]
