# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Backend discovery and dispatch for the regrid point-to-point interpolation.

A :class:`Backend` implements the actual numeric interpolation of gridded
values from one grid to another (e.g. via MIR, or precomputed sparse
matrices). :class:`BackendMaker` discovers the built-in backend classes
(any module in this package exposing a ``backend`` attribute) and caches
constructed backend instances, keyed by their construction arguments so a
given backend configuration is only instantiated once.
"""

from __future__ import annotations

import logging
import os
from abc import ABCMeta, abstractmethod

# from collections import namedtuple
from importlib import import_module
from types import ModuleType
from typing import TYPE_CHECKING, Any, Hashable, Tuple, Type, Union

from earthkit.geo.grids._regrid.data.utils import normalise_grid_spec as normalise_grid_spec

if TYPE_CHECKING:
    from eckit.geo import Grid as EckitGeoGrid
    from numpy.typing import NDArray

    from earthkit.geo.grids import Grid as EarthkitGeoGrid

    # A grid spec: a dict, a (possibly JSON-encoded) string, or a Grid
    # instance. ``earthkit.geo.grids.Grid`` currently re-exports
    # ``eckit.geo.Grid``, but the two are expected to diverge in the future,
    # so both are listed explicitly here.
    GridSpec = Union[dict, str, EckitGeoGrid, EarthkitGeoGrid]

LOG = logging.getLogger(__name__)


class Backend(metaclass=ABCMeta):
    """Abstract base class for point-to-point interpolation backends.

    Attributes
    ----------
    name : str or None
        The name the backend is registered and looked up under (e.g. via
        :func:`get_backend`).
    """

    name: str | None = None

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        pass

    @abstractmethod
    def regrid(
        self, values: NDArray, in_grid: GridSpec, out_grid: GridSpec, method: Any, **kwargs: Any
    ) -> Tuple[NDArray, Any]:
        """Interpolate ``values`` from ``in_grid`` onto ``out_grid``.

        Parameters
        ----------
        values : NDArray
            The values to interpolate, defined on ``in_grid``.
        in_grid : dict, str, eckit.geo.Grid or earthkit.geo.grids.Grid
            The input grid spec: a dict, a (possibly JSON-encoded) string, or
            a ``Grid`` instance. ``earthkit.geo.grids.Grid`` currently
            re-exports ``eckit.geo.Grid``, but the two are expected to
            diverge in the future.
        out_grid : dict, str, eckit.geo.Grid or earthkit.geo.grids.Grid
            The output grid spec, in the same forms as ``in_grid``.
        method : Any
            The interpolation method (e.g. ``"linear"``, ``"nearest-neighbour"``).
        **kwargs : Any
            Backend-specific extra arguments.

        Returns
        -------
        Tuple[NDArray, Any]
            The interpolated values and the (possibly backend-adjusted)
            output grid spec.
        """
        pass


class BackendLoader:
    """Loads a backend class from a module or a plugin entry point."""

    kind = "backend"

    def load_module(self, module: str) -> Type[Backend]:
        """Load the ``backend`` attribute from a submodule of this package.

        Parameters
        ----------
        module : str
            The (possibly relative) module name to import.

        Returns
        -------
        Type[Backend]
            The backend class exposed by the module as ``backend``.
        """
        return import_module(module, package=__name__).backend

    def load_entry(self, entry: Any) -> Type[Backend]:
        """Load a backend class from a plugin entry point.

        Parameters
        ----------
        entry : importlib.metadata.EntryPoint
            The entry point to load.

        Returns
        -------
        Type[Backend]
            The loaded backend class, either the entry point's target itself
            (if callable) or its ``backend`` attribute.
        """
        entry = entry.load()
        if callable(entry):
            return entry
        return entry.backend

    def load_remote(self, name: str) -> None:
        """Load a remote backend by name. Not implemented."""
        return None


class BackendMaker:
    """Discovers and constructs (and caches) built-in :class:`Backend` instances."""

    BACKENDS: dict[Hashable, Type[Backend]] = {}
    BACKEND_OBJECTS: dict[Hashable, Backend] = {}

    def __init__(self) -> None:
        self.BACKENDS = self._builtins()

    def _make_key(self, name: str, *args: Any, **kwargs: Any) -> Hashable:
        """Build a hashable cache key identifying a backend and its construction arguments.

        Parameters
        ----------
        name : str
            The backend name.
        *args : Any
            Positional arguments passed to the backend constructor.
        **kwargs : Any
            Keyword arguments passed to the backend constructor.

        Returns
        -------
        Hashable
            ``name`` alone when there are no extra arguments, otherwise a
            tuple combining ``name`` with ``args`` and the sorted ``kwargs``
            items.
        """
        if args or kwargs:
            key = [name, *args, *list(kwargs.items())]
            return tuple(key)
        else:
            return name

    def __call__(self, name: str, *args: Any, **kwargs: Any) -> Backend:
        """Get (constructing and caching if needed) a backend instance.

        Parameters
        ----------
        name : str
            The registered backend name.
        *args : Any
            Positional arguments passed to the backend constructor.
        **kwargs : Any
            Keyword arguments passed to the backend constructor.

        Returns
        -------
        Backend
            The (possibly cached) backend instance.

        Raises
        ------
        ValueError
            If ``name`` is not a known backend.
        """
        key = self._make_key(name, *args, **kwargs)
        if key in self.BACKEND_OBJECTS:
            return self.BACKEND_OBJECTS[key]

        if name in self.BACKENDS:
            klass = self.BACKENDS[name]
        else:
            # TODO: implement a plugin loader
            raise ValueError(f"Unknown backend: {name}")

        backend = klass(*args, **kwargs)
        self.BACKEND_OBJECTS[key] = backend

        return backend

    def _builtins(self) -> dict[Hashable, Type[Backend]]:
        """Scan for built-in backend classes.

        Imports every sibling module (or subpackage) of this package and, for
        each one exposing a ``backend`` attribute, registers it under the
        module's name, or, if ``backend`` is a dict, registers each of its
        items under their own key.

        Returns
        -------
        dict[Hashable, Type[Backend]]
            The discovered backend classes keyed by name.
        """
        r: dict[Hashable, Type[Backend]] = {}
        here = os.path.dirname(__file__)
        for path in sorted(os.listdir(here)):
            if path[0] in ("_", "."):
                continue

            if path.endswith(".py") or os.path.isdir(os.path.join(here, path)):
                name, _ = os.path.splitext(path)
                try:
                    module: ModuleType = import_module(f".{name}", package=__name__)
                    if hasattr(module, "backend"):
                        w = getattr(module, "backend")
                        if isinstance(w, dict):
                            for k, v in w.items():
                                r[k] = v
                        else:
                            r[name] = w
                except Exception:
                    LOG.exception("Error loading backend %s", name)

        LOG.debug(f"built-in backend classes: {r}")
        return r


MAKER = BackendMaker()


def get_backend(name: str, *args: Any, **kwargs: Any) -> Backend:
    """Get a backend by name.

    Parameters
    ----------
    name : str
        The registered backend name.
    *args : Any
        Positional arguments passed to the backend constructor.
    **kwargs : Any
        Keyword arguments passed to the backend constructor.

    Returns
    -------
    Backend
        The (possibly cached) backend instance.
    """
    return MAKER(name, *args, **kwargs)
