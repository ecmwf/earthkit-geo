# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import json
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

        fields = []
        for i, field in enumerate(data):
            r = self._regrid_field(field, i, backend, grid, **kwargs)
            fields.append(r)

        from earthkit.data import FieldList

        return FieldList.from_fields(fields)

    def _regrid_field(self, field, index, backend, grid, **kwargs):
        message = None
        if field._get_grib():
            message = field.message()

        # GRIB data
        if message is not None:
            # TODO: remove this check for the final 1.0.0 release
            try:
                from earthkit.data.field.grib.create import create_grib_field_from_message
            except ImportError:
                raise RuntimeError("Minimum version of earthkit-data >= 1.0.0rc4 is required for GRIB regridding")

            # currently this if branch means the mir backend
            if hasattr(backend, "regrid_grib"):
                res_message = backend.regrid_grib(message, grid, **kwargs)

            # currently this if branch means the precomputed backend
            else:
                from earthkit.geo.regrid.gridspec import GridSpec

                vv = field.to_numpy(flatten=True)

                in_grid = self.input_gridspec(field, index)

                # for precomputed backend we need to build a special gridspec object
                # to match the matrix inventory items
                # TODO: remove this limitation
                in_grid = GridSpec.from_any(in_grid)
                out_grid = GridSpec.from_any(grid)

                v_res, out_grid = backend.regrid(
                    vv,
                    in_grid,
                    out_grid,
                    **kwargs,
                )
                from eckit.geo import Grid

                out_grid = Grid(out_grid)
                # convert gridspec to string as this is what ecCodes expects in the "gridSpec" key w
                out_spec = json.dumps(out_grid.spec)

                # when the field has an associated GRIB message, we
                # use a GribEncoder to encode the resulting data and preserve
                # the original GRIB metadata as much as possible (bar the grid)

                # TODO: avoid using a GribEncoder once the grid(spec) handling in the
                # earthkit-data field is improved
                from earthkit.data.encoders.grib import GribEncoder

                encoder = GribEncoder()
                d = encoder.encode(template=message, values=v_res, gridSpec=out_spec)
                res_message = d.to_bytes()

            return create_grib_field_from_message(res_message, template_field=field)
        else:
            vv = field.to_numpy(flatten=True)
            in_grid = self.input_gridspec(field, index)

            # currently this is the precomputed backend
            if not hasattr(backend, "regrid_grib"):
                # for precomputed backend we need to build a special gridspec object
                # to match the matrix inventory items
                # TODO: remove this limitation
                from earthkit.geo.regrid.gridspec import GridSpec

                in_grid = GridSpec.from_any(in_grid)
                out_grid = GridSpec.from_any(grid)
            else:
                out_grid = grid

            v_res, out_grid = backend.regrid(
                vv,
                in_grid,
                out_grid,
                **kwargs,
            )
            from eckit.geo import Grid

            out_grid = Grid(out_grid)

            return field.set({"values": v_res, "geography.grid_spec": out_grid})

        return None


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
