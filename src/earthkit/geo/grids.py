# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


from eckit.geo import Grid as ECKIT_Grid


class Grid:
    """
    Geographical grid representation with support for various grid types.

    The Grid class provides a unified interface for working with different types
    of geographical grids used in numerical weather prediction, climate modeling,
    and geospatial applications. It wraps the eckit.geo Grid implementation and
    provides Pythonic access to grid properties and operations.

    The class uses a factory pattern: when you call `Grid.from_string()` or
    `Grid.from_dict()` with a grid code or spec, it automatically returns the appropriate
    specialised subclass (HEALPix, Octahedral, ReducedGaussian, etc.).

    Supported Grid Types
    --------------------
    - **Equirectangular** (list spec [dlon, dlat]): Regular lat-lon grid
    - **HEALPix** (code prefix 'H'): Hierarchical Equal Area isoLatitude Pixelization
    - **Octahedral** (code prefix 'O'): Octahedral reduced Gaussian grid
    - **ReducedGaussian** (code prefix 'N'): Original reduced Gaussian grid
    - **ICON** grids: Special grids like 'ICON_CH1_EPS', 'ICON_CH2_EPS'
    - **Generic**: Any other grid code supported by eckit.geo

    Parameters
    ----------
    **kwargs : dict
        Grid specification as keyword arguments. Typically includes:
        - grid : str
            The grid code (e.g., 'H32', 'O1280', 'N320')
        - Additional grid-specific parameters

    Examples
    --------
    >>> # Factory method automatically returns the right subclass
    >>> grid = Grid.from_string("H32")
    >>> isinstance(grid, HEALPix)
    True
    >>> grid.nside
    32

    >>> # Different grid types
    >>> octahedral = Grid.from_string("O1280")
    >>> reduced_gaussian = Grid.from_string("N320")

    >>> # From dictionary grid specification
    >>> grid = Grid.from_dict({"grid": "H64", "order": "nested"})
    >>> grid.nside
    64
    >>> grid.order
    'nested'

    >>> # Regular lat-lon grid from list specification
    >>> grid = Grid.from_dict({"grid": [1, 1]})
    >>> isinstance(grid, Equirectangular)
    True
    >>> grid.dlon
    1.0

    >>> # Get grid properties
    >>> grid = Grid.from_string("O640")
    >>> grid.size()  # Total number of points
    >>> grid.shape  # Grid shape
    >>> latlon = grid.to_latlon()  # Get lat/lon arrays
    """

    NAME = None

    def __init__(self, **kwargs):
        self._grid_spec = kwargs
        self._grid = ECKIT_Grid(self.grid_spec)

    @classmethod
    def _get_grid_class_and_params(cls, code):
        """
        Determine the appropriate Grid subclass and constructor parameters
        based on a grid code string or specification.

        Parameters
        ----------
        code : str or list
            The grid code string or grid specification. Can be:
            - String codes like "H32", "O1280", "N320"
            - List of [dlon, dlat] for regular lat-lon grids

        Returns
        -------
        tuple
            (grid_class, params_dict) where grid_class is the appropriate
            Grid subclass and params_dict contains the constructor parameters
        """
        # Handle list-based grid specifications (e.g., [5, 5] for regular lat-lon)
        if isinstance(code, (list, tuple)):
            if len(code) == 2:
                try:
                    dlon, dlat = float(code[0]), float(code[1])
                    return Equirectangular, {"dlon": dlon, "dlat": dlat}
                except (ValueError, TypeError):
                    pass  # Fall through to generic Grid

        if not isinstance(code, str):
            return Grid, {"grid": code}

        # Parse the code prefix to determine the grid type
        if code.startswith("H"):
            try:
                nside = int(code[1:])
                return HEALPix, {"nside": nside}
            except ValueError:
                pass  # Fall through to generic Grid
        elif code.startswith("O"):
            try:
                nlat = int(code[1:])
                return Octahedral, {"nlat": nlat}
            except ValueError:
                pass  # Fall through to generic Grid
        elif code.startswith("N"):
            try:
                nlat = int(code[1:])
                return ReducedGaussian, {"nlat": nlat}
            except ValueError:
                pass  # Fall through to generic Grid
        elif code == "ICON_CH1_EPS":
            return ICON_CH1_EPS, {}
        elif code == "ICON_CH2_EPS":
            return ICON_CH2_EPS, {}

        # Fallback to generic Grid for unknown or non-standard codes
        return Grid, {"grid": code}

    @classmethod
    def from_string(cls, code, **kwargs):
        """
        Create a Grid instance from a grid code string.

        When called on the base Grid class, this method automatically determines
        the appropriate subclass based on the code prefix:
        - "H" prefix: HEALPix grid (e.g., "H32")
        - "O" prefix: Octahedral grid (e.g., "O1280")
        - "N" prefix: Reduced Gaussian grid (e.g., "N320")
        - Special codes: "ICON_CH1_EPS", "ICON_CH2_EPS"

        When called on a subclass, it creates an instance of that subclass.

        Parameters
        ----------
        code : str
            The grid code string
        **kwargs
            Additional keyword arguments passed to the grid constructor

        Returns
        -------
        Grid or subclass
            An instance of the appropriate Grid subclass

        Examples
        --------
        >>> grid = Grid.from_string("H32")
        >>> isinstance(grid, HEALPix)
        True
        >>> grid.nside
        32
        """
        # Determine the appropriate class and base parameters
        grid_class, params = cls._get_grid_class_and_params(code)

        # If called on a subclass, verify the code matches that subclass
        if cls is not Grid:
            if grid_class is not cls:
                raise TypeError(
                    f"Code '{code}' resolves to {grid_class.__name__}, " f"but called on {cls.__name__}"
                )

        # Merge with any additional kwargs
        params.update(kwargs)

        return grid_class(**params)

    @classmethod
    def from_dict(cls, spec):
        """
        Create a Grid instance from a grid specification dictionary.

        When called on the base Grid class, this method automatically determines
        the appropriate subclass based on the "grid" key in the specification.
        The "grid" value can be:
        - A string code (e.g., "H32", "O1280", "N320")
        - A list of [dlon, dlat] for regular latitude-longitude grids

        When called on a subclass, it creates an instance of that subclass.

        Parameters
        ----------
        spec : dict
            Dictionary containing grid specification. Should contain a "grid" key
            with either a grid code string or a list specification.

        Returns
        -------
        Grid or subclass
            An instance of the appropriate Grid subclass

        Examples
        --------
        >>> # String-based grid codes
        >>> grid = Grid.from_dict({"grid": "H32", "order": "nested"})
        >>> isinstance(grid, HEALPix)
        True
        >>> grid.nside
        32
        >>> grid.order
        'nested'

        >>> # List-based grid specification for regular lat-lon
        >>> grid = Grid.from_dict({"grid": [5, 5]})
        >>> isinstance(grid, Equirectangular)
        True
        >>> grid.dlon
        5.0
        >>> grid.dlat
        5.0
        """
        if "grid" in spec:
            code = spec["grid"]
            grid_class, params = cls._get_grid_class_and_params(code)

            # If called on a subclass, verify the spec matches that subclass
            if cls is not Grid:
                if grid_class is not cls:
                    raise TypeError(
                        f"Grid spec {spec} resolves to {grid_class.__name__}, "
                        f"but called on {cls.__name__}"
                    )

            # For subclasses (not base Grid), extract only the extra params from spec
            # (excluding "grid" which the subclass constructors build themselves)
            if grid_class is not Grid:
                # Start with parsed params (e.g., nside=32 for HEALPix)
                merged_params = params.copy()
                # Add any extra params from spec, excluding "grid"
                for key, value in spec.items():
                    if key != "grid":
                        merged_params[key] = value
            else:
                # For base Grid class, pass everything including "grid"
                merged_params = spec

            return grid_class(**merged_params)

        # Fallback to generic Grid if no "grid" key
        return Grid(**spec)

    @property
    def name(self):
        """Human-readable name of the grid type."""
        if self.NAME is not None:
            return self.NAME
        return str(self)

    @property
    def grid_spec(self):
        """Grid specification dictionary."""
        return self._grid_spec

    def __repr__(self):
        return super().__repr__()

    def __str__(self):
        return str(self.grid_spec)

    def size(self):
        """Total number of grid points."""
        return self._grid.size()

    @property
    def shape(self):
        """Shape of the grid."""
        return self._grid.shape

    def bounding_box(self):
        """Bounding box of the grid."""
        return self._grid.bounding_box()

    def to_latlon(self):
        """Get latitude and longitude arrays for all grid points."""
        import numpy as np

        lats, lons = self._grid.to_latlons()
        return {"lat": np.array(lats), "lon": np.array(lons)}

    @property
    def code(self):
        """Grid code string."""
        return self.grid_spec["grid"]

    def plot(self, *args, **kwargs):
        """
        Visualise the grid on a map.

        Parameters
        ----------
        *args : tuple
            Additional arguments passed to the `earthkit.plots.Map` constructor.
        **kwargs : dict
            Additional keyword arguments passed to the `earthkit.plots.Map` constructor.
        """
        try:
            import earthkit.plots
        except ImportError:
            raise ImportError(
                "earthkit.plots is required to plot grids. Please install it with `pip install earthkit-plots`."
            )
        import cartopy.crs as ccrs

        chart = earthkit.plots.Map(*args, size=(7, 4), **kwargs)
        latlon = self.to_latlon()
        lats = latlon["lat"]
        lons = latlon["lon"]

        # Get the underlying cartopy geoaxes and its projection
        ax = chart.ax
        projection = ax.projection

        # Get the current extent of the geoaxes
        extent = ax.get_extent()

        # Transform lat/lon points to the projection's coordinate system
        # Source CRS is PlateCarree (lat/lon), target is the map's projection
        transformed = projection.transform_points(ccrs.PlateCarree(), lons, lats)
        x = transformed[..., 0]
        y = transformed[..., 1]

        # Filter points that fall within the extent
        mask = (x >= extent[0]) & (x <= extent[1]) & (y >= extent[2]) & (y <= extent[3])
        filtered_lons = lons[mask]
        filtered_lats = lats[mask]

        chart.coastlines(color="#333", linewidth=0.4)
        chart.borders(color="#333", linewidth=0.4, linestyle="--")
        chart.gridlines(draw_labels=False, color="#ddd")
        chart.scatter(x=filtered_lons, y=filtered_lats, marker=".", c="red", edgecolors="none", zorder=0, s=4)

        title_string = f"{self.NAME} grid ({self.code})"
        if "domain" in kwargs:
            title_string += "\n[Shown over {domain}]"

        chart.title(title_string, fontsize=10, color="grey")

        chart.show()


class Equirectangular(Grid):
    """
    Regular latitude-longitude grid (also known as equirectangular or plate carrée grid).

    A regular latitude-longitude grid has uniform spacing in both latitude and longitude
    directions. Grid cells are defined by their angular size in degrees. This is one of
    the simplest and most common grid types used in climate and meteorological data.

    Parameters
    ----------
    dlon : float
        Grid spacing in the longitude direction (degrees). For example, 1.0 creates
        grid cells that are 1° wide in longitude.
    dlat : float
        Grid spacing in the latitude direction (degrees). For example, 1.0 creates
        grid cells that are 1° tall in latitude.

    Examples
    --------
    Create a 5° × 5° grid:

    >>> grid = Equirectangular(dlon=5, dlat=5)
    >>> grid.dlon
    5.0
    >>> grid.dlat
    5.0

    Using the factory method:

    >>> grid = Grid.from_dict({"grid": [5, 5]})
    >>> isinstance(grid, Equirectangular)
    True
    >>> grid.dlon
    5.0
    >>> grid.dlat
    5.0

    Visualise the grid:

    .. plot::
       :include-source:

       from earthkit.geo.grids import Equirectangular

       grid = Equirectangular(dlon=5, dlat=5)
       grid.plot()
    """

    NAME = "Equirectangular"

    def __init__(self, dlon, dlat):
        super().__init__(grid=[dlon, dlat])

    @property
    def dlon(self):
        """Grid spacing in longitude direction (degrees)."""
        return self.grid_spec["grid"][0]

    @property
    def dlat(self):
        """Grid spacing in latitude direction (degrees)."""
        return self.grid_spec["grid"][1]

    @property
    def code(self):
        """Grid code string."""
        return f"{self.dlon}° × {self.dlat}°"


class HEALPix(Grid):
    """
    HEALPix (Hierarchical Equal Area isoLatitude Pixelization) grid.

    A HEALPix grid divides the sphere into equal-area pixels arranged in
    iso-latitude rings. It is commonly used in astrophysics and climate science
    for its equal-area property and efficient hierarchical structure.

    Parameters
    ----------
    nside : int
        The HEALPix resolution parameter. Total number of pixels is :math:`12 \\times nside^2`.
        Common values: 32, 64, 128, 256, 512, 1024, 2048, etc.
    order : str, optional
        Pixel ordering scheme. Default is 'ring'.
        - 'ring': Pixels are ordered along iso-latitude rings
        - 'nested': Hierarchical tree-based pixel ordering

    Examples
    --------
    Create a HEALPix grid:

    >>> grid = HEALPix(nside=8)
    >>> grid.nside
    8
    >>> grid.order
    'ring'

    Using the factory method:

    >>> grid = Grid.from_dict({"grid": "H16"})
    >>> isinstance(grid, HEALPix)
    True
    >>> grid.nside
    16

    Visualise the grid structure:

    .. plot::
       :include-source:

       from earthkit.geo.grids import HEALPix

       grid = HEALPix(nside=32, order='ring')
       grid.plot()
    """

    NAME = "HEALPix"

    def __init__(self, nside, order="ring"):
        super().__init__(grid=f"H{nside}", order=order)

    @property
    def nside(self):
        """The HEALPix resolution parameter. Total number of pixels is :math:`12 \\times nside^2`."""
        code = self.grid_spec["grid"]
        return int(code[1:])

    @property
    def order(self):
        """HEALPix ordering scheme ('ring' or 'nested')."""
        return self.grid_spec.get("order", "ring")


class ReducedGaussian(Grid):
    """
    Reduced Gaussian grid.

    A reduced Gaussian grid uses Gaussian latitudes but reduces the number of
    longitude points at higher latitudes where meridians converge. This reduces
    computational cost while maintaining numerical accuracy. Widely used in
    numerical weather prediction and climate modeling.

    The grid code follows the pattern 'NXXX' where XXX is the number of latitude
    lines between pole and equator.

    Parameters
    ----------
    nlat : int
        Number of latitude lines between pole and equator.
        Common values: 32, 48, 64, 80, 128, 160, 200, 256, 320, 400, 640, 1280.

    Examples
    --------
    Create a reduced Gaussian grid:

    >>> grid = ReducedGaussian(nlat=32)
    >>> grid.nlat
    32

    Using the factory method:

    >>> grid = Grid.from_dict({"grid": "N48"})
    >>> isinstance(grid, ReducedGaussian)
    True
    >>> grid.nlat
    48

    Visualise the grid structure:

    .. plot::
       :include-source:

       from earthkit.geo.grids import ReducedGaussian

       grid = ReducedGaussian(64)
       grid.plot()

    See Also
    --------
    Octahedral : Octahedral variant of reduced Gaussian grid.
    """

    NAME = "Reduced Gaussian"

    def __init__(self, nlat):
        super().__init__(grid=f"N{nlat}")

    @property
    def nlat(self):
        """Number of latitude lines between pole and equator."""
        code = self.grid_spec["grid"]
        return int(code[1:])


class Octahedral(Grid):
    """
    Octahedral reduced Gaussian grid.

    An octahedral grid is a variant of the reduced Gaussian grid that applies
    a new rule for computing the number of longitude points per latitude circle.
    Inspired by the Collignon projection of a sphere onto an octahedron, it has
    20 longitude points at the latitude circle closest to the poles, with the
    number of points increasing by 4 at each latitude line toward the equator.

    The octahedral grid has approximately 20% fewer grid points than the equivalent
    original reduced Gaussian grid (e.g., O1280 vs N1280), reducing computational
    cost while improving the calculation of local derivatives in grid-point space.

    Total number of grid points: :math:`4 \\times nlat \\times (nlat + 9)`.

    The grid code follows the pattern 'OXXX' where XXX is the number of latitude
    lines between pole and equator.

    Parameters
    ----------
    nlat : int
        Number of latitude lines between pole and equator.
        Common values: 32, 48, 64, 80, 128, 160, 200, 256, 320, 400, 640, 1280.

    Examples
    --------
    Create an octahedral grid:

    >>> grid = Octahedral(nlat=32)
    >>> grid.nlat
    32

    Using the factory method:

    >>> grid = Grid.from_dict({"grid": "O48"})
    >>> isinstance(grid, Octahedral)
    True
    >>> grid.nlat
    48

    Visualise the grid structure:

    .. plot::
       :include-source:

       from earthkit.geo.grids import Octahedral

       grid = Octahedral(32)
       grid.plot()

    See Also
    --------
    ReducedGaussian : Original reduced Gaussian grid implementation.
    """

    NAME = "Octahedral"

    def __init__(self, nlat):
        super().__init__(grid=f"O{nlat}")

    @property
    def nlat(self):
        """Number of latitude lines between pole and equator."""
        code = self.grid_spec["grid"]
        return int(code[1:])


class ICON_CH1_EPS(Grid):
    """
    ICON-CH1-EPS grid for Switzerland.

    A predefined unstructured grid used by the ICON (ICOsahedral Nonhydrostatic)
    ensemble prediction system for the CH1 (Switzerland) domain. This is a fixed
    grid specification with no configurable parameters.

    Examples
    --------
    >>> grid = ICON_CH1_EPS()
    >>> grid.code
    'ICON_CH1_EPS'

    Using the factory method:

    >>> grid = Grid.from_dict({"grid": "ICON_CH1_EPS"})
    >>> isinstance(grid, ICON_CH1_EPS)
    True

    Visualise the grid:

    .. plot::
       :include-source:

       from earthkit.geo.grids import ICON_CH1_EPS

       grid = ICON_CH1_EPS()
       grid.plot()

    See Also
    --------
    ICON_CH2_EPS : ICON-CH2-EPS grid for extended Switzerland domain
    """

    NAME = "ICON_CH1_EPS"

    def __init__(self):
        super().__init__(grid="ICON_CH1_EPS")


class ICON_CH2_EPS(Grid):
    """
    ICON-CH2-EPS grid for extended Switzerland domain.

    A predefined unstructured grid used by the ICON (ICOsahedral Nonhydrostatic)
    ensemble prediction system for the CH2 (extended Switzerland) domain. This is
    a fixed grid specification with no configurable parameters.

    Examples
    --------
    >>> grid = ICON_CH2_EPS()
    >>> grid.code
    'ICON_CH2_EPS'

    Using the factory method:

    >>> grid = Grid.from_dict({"grid": "ICON_CH2_EPS"})
    >>> isinstance(grid, ICON_CH2_EPS)
    True

    Visualise the grid:

    .. plot::
       :include-source:

       from earthkit.geo.grids import ICON_CH2_EPS

       grid = ICON_CH2_EPS()
       grid.plot()

    See Also
    --------
    ICON_CH1_EPS : ICON-CH1-EPS grid for Switzerland domain
    """

    NAME = "ICON_CH2_EPS"

    def __init__(self):
        super().__init__(grid="ICON_CH2_EPS")
