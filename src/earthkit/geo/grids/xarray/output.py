# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import logging

from earthkit.geo.grids._regrid.data.utils import create_grid_object

LOG = logging.getLogger(__name__)

# The names the CF "reduced_gaussian" encoding of a regridded output writes its
# grid under, as in the reference file of a reduced Gaussian NetCDF: the point
# dimension and its dimension coordinate of the same name (decision D2), the
# global latitude lines the points per latitude are counted along, the table of
# those counts, and the scalar variable carrying the grid mapping.
CF_POINT_DIM = "reduced_gaussian_index"
CF_LATITUDE_DIM = "lat"
CF_PL_NAME = "pl"
CF_MAPPING_NAME = "reduced_gaussian"

# the names of the tables the encoding writes, the point dimension aside
CF_TABLE_NAMES = (CF_LATITUDE_DIM, CF_PL_NAME, CF_MAPPING_NAME)


# the CF attributes of the variables the encoding writes, keyed by their role
_CF_ATTRS = {
    "index": {"standard_name": "reduced_gaussian_index", "units": "1"},
    "latitude": {"standard_name": "latitude", "units": "degrees_north"},
    "longitude": {"standard_name": "longitude", "units": "degrees_east"},
    "pl": {"long_name": "number of points per latitude", "units": "1"},
    "mapping": {"grid_mapping_name": "reduced_gaussian"},
}


def _cf_attrs(role: str, **extra) -> dict:
    """The CF attributes of a variable playing ``role``, a fresh dictionary each time.

    xarray keeps the dictionary it is handed, so two variables must not share one;
    the attributes of a grid mapping also name the table and the latitude dimension
    of *that* grid.
    """
    return dict(_CF_ATTRS[role], **extra)


# This is a temporary wrapper for the Grid interface
class GridWrapper:
    """Thin wrapper around an ``eckit.geo.Grid`` used to build output geography.

    Normalises access to a grid built from a grid spec (or an existing
    ``Grid`` instance) and adds helpers to extract flat or distinct
    lat/lon arrays for a given field shape.
    """

    def __init__(self, grid):
        """Initialise the wrapper.

        Parameters
        ----------
        grid_spec : Any
            A grid spec (dict/str) or an existing ``eckit.geo.Grid`` instance.
        """
        self._grid = create_grid_object(grid)
        self._latlons = None

    def __getattr__(self, name):
        """Delegate unknown attribute access to the wrapped ``Grid``."""
        return getattr(self._grid, name)

    def to_latlons(self):
        """Get the flat latitude and longitude arrays for the grid.

        The arrays are computed by the grid once and then cached: the
        geography builder asks for them while deciding the output dimensions
        and again while building the coordinates, and computing them can be
        expensive on a large grid.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            The latitude and longitude arrays.
        """
        if self._latlons is None:
            import numpy as np

            lat, lon = self._grid.to_latlons()
            self._latlons = (np.array(lat), np.array(lon))
        return self._latlons

    def is_spectral(self):
        """bool: Whether the grid is spectral (always False here)."""
        return False

    def to_distinct_latlons(self, field_shape):
        """Get the distinct (1D) latitude and longitude arrays for a 2D field.

        The grid is a regular mesh when its latitudes are constant along
        every row of the field and its longitudes constant along every
        column (which a rotated grid is not, however few points its rows
        and columns have); the distinct arrays are then the latitude of
        each row and the longitude of each column.

        Parameters
        ----------
        field_shape : Tuple[int, int]
            The shape of the field the grid is used for.

        Returns
        -------
        Tuple[Optional[np.ndarray], Optional[np.ndarray]]
            The distinct latitude and longitude arrays, or ``(None, None)``
            if the grid is not a regular mesh matching ``field_shape``.
        """
        import numpy as np

        if len(self._grid.shape) == 2:
            lat, lon = self.to_latlons()
            lat = lat.reshape(self._grid.shape)
            lon = lon.reshape(self._grid.shape)
            if (
                np.allclose(lat, lat[:, :1])
                and np.allclose(lon, lon[:1, :])
                and lat.shape[0] == field_shape[0]
                and lon.shape[1] == field_shape[1]
            ):
                return lat[:, 0], lon[0, :]

        return None, None


class XarrayOutputGeographyBuilder:
    """Builds output geography (dims/coords) for a regridded xarray variable.

    Wraps an output grid spec and derives the dimension names, coordinate
    arrays and coordinate-to-dimension mapping to attach to the regridded
    result.
    """

    def __init__(self, grid):
        """Initialise the builder with the given output grid.

        Parameters
        ----------
        grid : Any
            The output grid spec (dict/str) or an ``eckit.geo.Grid`` instance.
        """
        self.grid = GridWrapper(grid)
        self._distinct = None
        self._cf = None
        self._cf_checked = False

    @property
    def shape(self):
        """Tuple[int, ...]: The shape of the output grid."""
        return self.grid.shape

    @property
    def cf_reduced_gaussian(self):
        """Optional[ReducedGaussianEncoding]: How to write the grid as a CF reduced Gaussian.

        A plain octahedral ``O<n>`` or classical ``N<n>`` grid, global or with an
        area, is encoded with the CF ``reduced_gaussian`` grid mapping: the index
        of every point on the global grid, the global latitude lines and the table
        of points per latitude say everything about it (see
        :func:`~.reduced_gaussian.cf_reduced_gaussian_encoding`). Any other grid —
        a regular one, a HEALPix one, a reduced Gaussian one with a rotation or a
        table of its own — has no encoding, and its output keeps the per-point
        coordinates it always had (decision A10).

        The encoding is computed when it is first asked for, and cached: it reads
        the points of the grid, which a large grid does not take a moment to do.
        """
        if not self._cf_checked:
            from earthkit.geo.grids.xarray.reduced_gaussian import cf_reduced_gaussian_encoding

            self._cf_checked = True
            try:
                self._cf = cf_reduced_gaussian_encoding(self.grid)
            except Exception as e:
                LOG.debug("Could not encode the output grid as a CF reduced Gaussian: %s", e)
                self._cf = None
        return self._cf

    def _distinct_latlons(self):
        """Decide whether the grid carries distinct 1D latitudes and longitudes.

        This is the decision :meth:`coords` makes between latitude/longitude
        coordinates on their own dimensions and 2D ones on ``(y, x)``, cached
        so that :meth:`geo_dims` and :meth:`coords` share it, as well as the
        latitudes and longitudes the grid computes for it.

        Returns
        -------
        Tuple[Optional[np.ndarray], Optional[np.ndarray]]
            The distinct latitude and longitude arrays when the grid is a
            regular mesh matching its shape, or ``(None, None)``.
        """
        if self._distinct is None:
            lat = lon = None
            if len(self.shape) == 2:
                try:
                    d_lat, d_lon = self.grid.to_distinct_latlons(self.shape)
                except Exception as e:
                    LOG.debug("Could not extract distinct lat/lon coordinates: %s", e)
                else:
                    if (
                        d_lat is not None
                        and d_lon is not None
                        and len(d_lat) == self.shape[0]
                        and len(d_lon) == self.shape[1]
                    ):
                        lat, lon = d_lat, d_lon
            self._distinct = (lat, lon)
        return self._distinct

    def geo_dims(self):
        """Determine the geographical dimensions of the dataset.

        The decision is the one :meth:`coords` makes, shared through
        :meth:`_distinct_latlons`: a 2D grid has distinct 1D latitudes and
        longitudes (constant along every row and column, which a rotated
        grid is not) only when its coordinates are built on the ``latitude``
        and ``longitude`` dimensions; any other 2D grid carries 2D
        coordinates on ``(y, x)``. A 1D grid has a single point dimension,
        named after the index which numbers its points when the grid is
        encoded as a CF reduced Gaussian one (decision D2).

        Returns
        -------
        List[str]
            ``["latitude", "longitude"]`` for a 2D grid with distinct 1D
            latitudes and longitudes, ``["y", "x"]`` for any other 2D grid
            (e.g. a rotated one), and ``["reduced_gaussian_index"]`` or
            ``["values"]`` for a 1D grid.

        Raises
        ------
        ValueError
            If the grid has no geography this builder supports.
        """
        num = len(self.shape)
        if num >= 2:
            lat, lon = self._distinct_latlons()
            if lat is not None and lon is not None:
                return ["latitude", "longitude"]
            return ["y", "x"]
        if num == 1:
            return [self.point_dim]

        raise ValueError("Geography is not supported.")

    @property
    def point_dim(self):
        """Hashable: The name of the dimension the points of a 1D grid are numbered by.

        ``reduced_gaussian_index`` for a grid the CF encoding writes, the point
        dimension of the reference file (decision D2); ``values`` for any other
        1D grid, as it has always been.
        """
        return CF_POINT_DIM if self.cf_reduced_gaussian is not None else "values"

    def cf_indices(self):
        """np.ndarray: The global index of every point of the grid, in the order a CF file holds them.

        That is the increasing order the
        :attr:`~.reduced_gaussian.ReducedGaussianEncoding.order` of
        :attr:`cf_reduced_gaussian` sorts the points of the grid into (decision
        A3), the order :meth:`add_cf_reduced_gaussian` gathers the values into.
        """
        import numpy as np

        encoding = self.cf_reduced_gaussian
        if encoding is None:
            raise ValueError("The grid is not a plain reduced Gaussian one, and has no CF index.")

        indices = np.asarray(encoding.indices, dtype=np.int64).reshape(-1)
        if encoding.order is not None:
            indices = indices[np.asarray(encoding.order, dtype=np.int64).reshape(-1)]

        # the reference file numbers the points of its grid with 32 bit integers,
        # and an index which cannot be held by them is beyond any grid a machine
        # holds in memory: widened, never wrapped
        if indices.size and int(indices[-1]) < np.iinfo(np.int32).max:
            indices = indices.astype(np.int32)

        return indices

    def cf_coords(self, point_dim=None, as_dataset=True) -> dict:
        """The CF coordinates which locate a grid on the points it numbers, as ``xr.Variable`` objects.

        Parameters
        ----------
        point_dim : Optional[Hashable]
            The name of the point dimension, :attr:`point_dim` by default. A
            caller who named the output dimensions names them too.
        as_dataset : bool
            Whether the object to encode is a Dataset: its points are located by
            the global latitude lines the ``pl`` table counts (decision D3), the
            points of a DataArray by their own latitude and longitude (decision
            A4), which is all a DataArray can hold of a grid mapping.

        Returns
        -------
        Dict[Hashable, xr.Variable]
            The index coordinate, holding the increasing global index of every
            point, and either the latitude lines of the global grid or the
            latitude and longitude of every point of this one.
        """
        import numpy as np
        import xarray as xr

        point_dim = self.point_dim if point_dim is None else point_dim

        coords = {point_dim: xr.Variable((point_dim,), self.cf_indices(), _cf_attrs("index"))}
        if as_dataset:
            coords[CF_LATITUDE_DIM] = xr.Variable(
                (CF_LATITUDE_DIM,),
                np.asarray(self.cf_reduced_gaussian.latitudes, dtype=np.float64),
                _cf_attrs("latitude"),
            )
            return coords

        lat, lon = self.grid.to_latlons()
        if lat is None or lon is None:
            raise ValueError(f"The grid '{self.grid.spec}' gives no latitude and longitude for its points.")

        order = self.cf_reduced_gaussian.order
        if order is not None:
            order = np.asarray(order, dtype=np.int64).reshape(-1)
            lat, lon = np.asarray(lat).reshape(-1)[order], np.asarray(lon).reshape(-1)[order]

        coords["latitude"] = xr.Variable((point_dim,), np.asarray(lat), _cf_attrs("latitude"))
        coords["longitude"] = xr.Variable((point_dim,), np.asarray(lon), _cf_attrs("longitude"))
        return coords

    def cf_data_vars(self) -> dict:
        """The CF data variables which describe a reduced Gaussian grid: its ``pl`` table and its mapping.

        Returns
        -------
        Dict[Hashable, xr.Variable]
            The number of points on every latitude line of the *global* grid,
            along the latitude dimension the mapping names, and the scalar grid
            mapping naming the grid, its subtype and those two variables — the
            tables a Dataset holds to place its points (decision D3). Empty for a
            grid with no CF encoding.
        """
        import numpy as np
        import xarray as xr

        encoding = self.cf_reduced_gaussian
        if encoding is None:
            return {}

        return {
            CF_PL_NAME: xr.Variable(
                (CF_LATITUDE_DIM,),
                np.asarray(encoding.pl, dtype=np.int32),
                _cf_attrs("pl"),
            ),
            CF_MAPPING_NAME: xr.Variable(
                (),
                np.int64(0),
                _cf_attrs(
                    "mapping",
                    grid_subtype=encoding.subtype,
                    latitude_dimension=CF_LATITUDE_DIM,
                    points_per_latitude=CF_PL_NAME,
                ),
            ),
        }

    def add_cf_reduced_gaussian(self, obj, point_dim=None, as_dataset=True):
        """Gather a regridded object onto the CF point order, and encode its grid in it.

        The values of ``obj`` number the points of the grid in the order of the
        grid itself; a CF file must hold them in the increasing order of their
        global index, so they are gathered with that permutation first (decision
        A3) — an indexing which leaves a chunked object as lazy as it was.

        A Dataset is then encoded in full (decision D3): the index coordinate
        naming every point, the global latitude lines as the ``lat`` coordinate,
        the ``pl`` table and the grid mapping as data variables, and a
        ``grid_mapping`` attribute on each of the regridded variables. It holds no
        per-point latitude or longitude: the mapping and the index say where its
        points are. A DataArray can hold neither the table nor the mapping a
        variable refers to, so it is encoded partially (decisions D4 and A4): its
        index coordinate, and the latitude and longitude of every point.

        Parameters
        ----------
        obj : xr.Dataset or xr.DataArray
            The regridded object, its values on the points of the grid.
        point_dim : Optional[Hashable]
            The name of the point dimension of ``obj``, which a caller who named
            its ``out_dims`` has chosen over :attr:`point_dim`.
        as_dataset : bool
            Whether ``obj`` is a Dataset, and can therefore hold the tables.

        Returns
        -------
        xr.Dataset or xr.DataArray
            ``obj``, gathered onto the CF order and encoded.

        Raises
        ------
        ValueError
            If the grid has no CF encoding, or a name the encoding writes is
            already taken by the regridded data (decision A9).
        """
        import numpy as np

        encoding = self.cf_reduced_gaussian
        if encoding is None:
            raise ValueError("The grid is not a plain reduced Gaussian one, and has no CF encoding to write.")

        point_dim = self.point_dim if point_dim is None else point_dim
        clash = _cf_name_clash(obj, point_dim, as_dataset)
        if clash:
            raise ValueError(
                f"Cannot encode the reduced Gaussian grid '{encoding.grid}' in the regridded data: "
                f"the encoding writes {', '.join(clash)}. Name the points with an 'out_dims' "
                "dimension which leaves those names free, or rename the data which takes them."
            )

        if encoding.order is not None:
            obj = obj.isel({point_dim: np.asarray(encoding.order, dtype=np.int64).reshape(-1)})

        if not as_dataset:
            return obj.assign_coords(self.cf_coords(point_dim, as_dataset=False))

        names = list(obj.data_vars)
        obj = obj.assign_coords(self.cf_coords(point_dim))
        for name, var in self.cf_data_vars().items():
            obj[name] = var
        for name in names:
            obj[name].attrs["grid_mapping"] = CF_MAPPING_NAME

        return obj

    def coords(self):
        """Build the output coordinate arrays for the grid.

        Returns
        -------
        Tuple[Dict[str, int], Dict[str, np.ndarray], Dict[str, Tuple[str, ...]]]
            The output dimension sizes, the coordinate arrays (e.g.
            ``latitude``/``longitude``), and the dimensions each coordinate
            is defined on.

        Notes
        -----
        A grid the CF reduced Gaussian encoding writes is located by the index of
        its points and the global latitude lines, not by the latitude and
        longitude of every point (decision D3); :meth:`cf_coords` hands those out
        with the attributes a CF file requires of them.
        """
        import math

        field_shape = self.grid.shape

        coords = {}
        dims = {}
        coords_dim = {}

        if self.grid.is_spectral():
            if len(field_shape) == 1:
                dims["values"] = field_shape[0]
        else:
            if len(field_shape) == 1:
                point_dim = self.point_dim
                if self.cf_reduced_gaussian is not None:
                    for name, var in self.cf_coords().items():
                        coords[name] = var.data
                        coords_dim[name] = var.dims
                        for d in var.dims:
                            dims[d] = var.sizes[d]
                else:
                    dims[point_dim] = field_shape[0]
                    try:
                        lat, lon = self.grid.to_latlons()
                        if lat is not None and lon is not None:
                            coords["latitude"] = lat
                            coords["longitude"] = lon
                            coords_dim = {k: (point_dim,) for k in coords}
                    except Exception:
                        pass
            elif len(field_shape) == 2:
                lat, lon = self._distinct_latlons()
                if lat is not None and lon is not None:
                    coords["latitude"] = lat
                    coords["longitude"] = lon
                    coords_dim["latitude"] = ("latitude",)
                    coords_dim["longitude"] = ("longitude",)
                    dims["latitude"] = lat.size
                    dims["longitude"] = lon.size

                if not coords or not dims:
                    lat, lon = self.grid.to_latlons()
                    # print("to_latlons:", type(lat), type(lon))
                    if lat is not None and lon is not None:
                        lat = lat.reshape(field_shape)
                        lon = lon.reshape(field_shape)
                        coords["latitude"] = lat
                        coords["longitude"] = lon
                        coords_dim = {k: ("y", "x") for k in coords}
                        dims["y"] = field_shape[0]
                        dims["x"] = field_shape[1]
                        # print("field_shape:", field_shape, lat.shape, lon.shape)
                        assert coords["latitude"].shape == field_shape
                        assert coords["longitude"].shape == field_shape
            else:
                raise ValueError(f"Unsupported field shape {field_shape}")

        for k, v in coords.items():
            assert k in coords_dim, f"{k=}, {coords_dim=}"
            assert all(x in dims for x in coords_dim[k]), f"{k=}, {coords_dim=} {dims=}"
            assert v.size == math.prod([dims[x] for x in coords_dim[k]])

        return dims, coords, coords_dim


def _cf_name_clash(obj, point_dim, as_dataset=True) -> list:
    """The names of the CF reduced Gaussian encoding a regridded object has taken already.

    Parameters
    ----------
    obj : xr.Dataset or xr.DataArray
        The regridded data the encoding is about to be written into.
    point_dim : Hashable
        The dimension the points are numbered by, which takes a coordinate of its
        own name.
    as_dataset : bool
        Whether ``obj`` is a Dataset, which writes the tables next to the points.

    Returns
    -------
    List[str]
        One phrase per name in the way, empty when the encoding fits. The
        dimension of the points is not in the way of its own coordinate, so a
        variable which merely sits on it is no clash; a variable or coordinate
        NAMED like it is, because xarray has made the dimension coordinate of
        that name out of it, or will refuse the result.

    Notes
    -----
    Writing a table over data which holds one of those names would silently drop
    the regridded variable, or place the grid's latitude lines on the points of
    the grid: the clash is refused instead (decision A9). A regridded variable
    named like the point dimension has already been turned into its dimension
    coordinate when the result was assembled, and the index coordinate the
    encoding writes would overwrite it.
    """
    if not as_dataset:
        # a DataArray writes no tables, but the index coordinate of its point
        # dimension goes into it all the same, over any name it finds there
        if obj.name == point_dim or point_dim in obj.coords:
            return [f"a variable or coordinate named '{point_dim}' taking the name of its point dimension"]
        return []

    taken = (set(obj.variables) | set(obj.dims)) - {point_dim}

    clash = [f"a variable, coordinate or dimension named '{name}'" for name in CF_TABLE_NAMES if name in taken]
    if point_dim in obj.variables:
        clash.append(f"a variable or coordinate named '{point_dim}' taking the name of its point dimension")
    if point_dim in CF_TABLE_NAMES:
        clash.append(f"the name of a table as the name of its point dimension ('{point_dim}')")

    return clash
