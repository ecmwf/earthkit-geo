# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Precomputed interpolation matrix inventory for the "precomputed" regrid backend.

Defines :class:`MatrixDb`, which looks up a precomputed sparse
interpolation matrix for a given input/output grid pair and interpolation
method. The inventory is described by an ``index.json`` file (loaded into a
:class:`MatrixIndex`) plus one ``.npz`` sparse-matrix file per entry, either
stored locally (:class:`LocalAccessor`) or fetched and cached from a remote
URL (:class:`UrlAccessor`/:class:`SystemAccessor`). :data:`SYS_DB` is the
built-in system inventory.
"""

import json
import logging
import os
from abc import ABCMeta, abstractmethod

from scipy.sparse import load_npz

from earthkit.geo.grids._regrid.backends.precomputed.gridspec import _GridWrapper
from earthkit.geo.utils import no_progress_bar
from earthkit.geo.utils.download import download_and_cache
from earthkit.geo.utils.url import join_url_path

LOG = logging.getLogger(__name__)

VERSION = 1

_SYSTEM_URL = "https://sites.ecmwf.int/repository/earthkit/regrid/db/1/"
_INDEX_FILENAME = "index.json"
_INDEX_SHA_FILENAME = "index.json.sha256"
_INDEX_GZ_FILENAME = "index.json.gz"
_METHOD_ALIAS = {"nearest-neighbour": ("nn", "nearest-neighbor")}

_GRIDBOX_DEFAULT = {
    "type": "grid-box-average",
    "nonLinear": [{"type": "missing-if-heaviest-missing"}],
    "solver": {"type": "multiply"},
    "cropping": False,
    "lsmWeightAdjustment": 0.2,
    "pruneEpsilon": 1e-10,
    "poleDisplacement": 0,
}


def is_gridbox_default(inter):
    """Check if the interpolation method is the default grid-box-average.

    In this case it should be just the string "grid-box-average" but now it
    is a dictionary. Until it is fixed in MIR we need this check.
    """
    method = inter["method"]
    if isinstance(method, dict):
        return method == _GRIDBOX_DEFAULT
    elif isinstance(method, str):
        return method == "grid-box-average"
    else:
        return False


def make_sha(data):
    """Compute the SHA-256 hex digest of ``data``.

    Parameters
    ----------
    data : str or Any
        The data to hash. A str is hashed as-is (UTF-8 encoded); anything
        else is first JSON-serialised (with sorted keys).

    Returns
    -------
    str
        The hex digest.
    """
    import hashlib

    m = hashlib.sha256()
    if isinstance(data, str):
        m.update(data.encode("utf-8"))
    else:
        m.update(json.dumps(data, sort_keys=True).encode("utf-8"))
    return m.hexdigest()


class MatrixAccessor(metaclass=ABCMeta):
    """Abstract base class giving access to the files of a matrix inventory.

    A ``MatrixAccessor`` knows where the inventory's ``index.json`` file and
    ``.npz`` matrix files live (locally or remotely) and how to fetch/reload
    them.
    """

    @abstractmethod
    def path(self):
        """Return the accessor's location (a local path or URL)."""
        pass

    @abstractmethod
    def is_local(self):
        """Return True if the inventory is stored locally."""
        pass

    @abstractmethod
    def index_path(self):
        """Return the local filesystem path of the inventory's index file."""
        pass

    @abstractmethod
    def matrix_path(self, name):
        """Return the local filesystem path of the named matrix file.

        Parameters
        ----------
        name : str
            The matrix file's relative path within the inventory, as given
            by :meth:`MatrixIndex.matrix_path`.

        Returns
        -------
        str
            The local filesystem path of the file.
        """
        pass

    @abstractmethod
    def reload(self, strict=False):
        """Force the index file to be (re-)fetched.

        Parameters
        ----------
        strict : bool, default=False
            Accessor-specific flag controlling how strictly the reload
            behaves.
        """
        pass

    def checked_remote(self):
        """Return True if the remote index has already been checked for updates."""
        return False

    @abstractmethod
    def reset(self):
        """Clear any cached state so the index is looked up afresh."""
        pass


class UrlAccessor(MatrixAccessor):
    """:class:`MatrixAccessor` fetching and caching the inventory from a URL."""

    def __init__(self, url):
        """Initialise the accessor with the inventory's base URL.

        Parameters
        ----------
        url : str
            The base URL the index and matrix files are served from.
        """
        self._url = url
        self._index_path = None
        self._checked_remote = False

    def path(self):
        """str: The base URL of the inventory."""
        return self._url

    def is_local(self):
        """bool: Always False."""
        False

    def checked_remote(self):
        """bool: Whether the remote index has already been checked for updates in this session."""
        return self._checked_remote

    def reset(self):
        """Clear the cached index path and remote-checked flag."""
        self._index_path = None
        self._checked_remote = False

    def reload(self, force=False):
        """Re-fetch the index file, checking the remote for updates.

        Parameters
        ----------
        force : bool, default=False
            When True, force re-downloading the remote checksum and index
            file even if a cached copy exists.
        """
        self._index_path = self._get_index(check_remote=True, force=force)

    def index_path(self):
        """Return the local (downloaded and cached) path of the index file.

        Returns
        -------
        str
            The local path of the index file, downloading and caching it
            first if not already available.
        """
        if self._index_path is None or not os.path.exists(self._index_path):
            self._index_path = self._get_index()
        return self._index_path

    def _get_index(self, check_remote=False, force=False):
        """Download (if needed) and return the local path of the uncompressed index file.

        Parameters
        ----------
        check_remote : bool, default=False
            If True, compare the local cached checksum against the remote
            one and re-download the index file if they differ.
        force : bool, default=False
            If True, unconditionally re-download the remote checksum and
            index file.

        Returns
        -------
        str
            The local path of the uncompressed index file.
        """
        from earthkit.geo.utils.caching import cache_file

        url = join_url_path(self._url, _INDEX_FILENAME)

        def _compare_sha(args, path, owner_data):
            """Decide if the index file should be downloaded and cached again."""
            LOG.info("UrlAccessor: compare local and remote index file checksums")
            LOG.info(f"UrlAccessor: cached entry {owner_data=}")
            if owner_data is None:
                return True
            local_sha = owner_data.get("sha256", None)
            LOG.info(f"UrlAccessor: local (cached) checksum={local_sha}")
            if local_sha is None:
                return True

            remote_sha = self._remote_sha()
            self._checked_remote = True

            if local_sha != remote_sha:
                LOG.info(
                    (
                        f"UrlAccessor: remote checksum={remote_sha} differs from "
                        "local (cached) checksum. Downloading new index file."
                    )
                )
                return True
            else:
                LOG.info(
                    (
                        f"UrlAccessor: remote checksum={remote_sha} is the same as "
                        "local (cached)  checksum. Use cached index file."
                    )
                )
                return False

        def _force_download(args, path, owner_data):
            """Decide if the index file should be downloaded and cached again."""
            LOG.info("UrlAccessor: forcefully download remote checksum and new index file")
            self._remote_sha()
            self._checked_remote = True
            return True

        def _create(target, args):
            """Download and cache gzipped index file, uncompress it and generate
            checksum from contents
            """
            path_gz = self._gzip_file()
            LOG.info(f"UrlAccessor: uncompress gzipped index file={path_gz}")
            import gzip

            with gzip.open(path_gz, "rb") as f:
                data = f.read()
                with open(target, "wb") as f_out:
                    f_out.write(data)

            with open(target, "r") as f:
                data = f.read()
                sha = make_sha(data)

            # the returned data will be stored in the cache in the entry's owner_data
            LOG.info(f"UrlAccessor: index file checksum={sha}")
            return {"sha256": sha}

        if force:
            force = _force_download
        elif check_remote:
            force = _compare_sha
        else:
            LOG.info("UrlAccessor: check if cached index file is available")
            force = None

        path = cache_file(
            "regrid",
            _create,
            (url,),
            force=force,
            extension=".cache",
        )

        LOG.info(f"UrlAccessor: index file={path}")
        return path

    def _remote_sha(self):
        """Download and return the remote index file's expected SHA-256 checksum.

        Returns
        -------
        str
            The checksum read from the remote ``.sha256`` file.

        Raises
        ------
        Exception
            If the checksum file cannot be downloaded (logged before being
            re-raised).
        """
        try:
            url = join_url_path(self._url, _INDEX_SHA_FILENAME)
            path = download_and_cache(
                url,
                owner="url",
                verify=True,
                force=True,
                chunk_size=1024 * 1024,
                http_headers=None,
                update_if_out_of_date=True,
                progress_bar=no_progress_bar,
                maximum_retries=0,
                retry_after=10,
            )
        except Exception:
            LOG.error(f"UrlAccessor: could not download index checksum file={url}")
            raise

        with open(path, "r") as f:
            sha = f.read().strip()
        return sha

    def _gzip_file(self):
        """Download and return the local path of the gzipped remote index file.

        Returns
        -------
        str
            The local (cached) path of the downloaded ``.json.gz`` file.

        Raises
        ------
        Exception
            If the file cannot be downloaded (logged before being re-raised).
        """
        try:
            url = join_url_path(self._url, _INDEX_GZ_FILENAME)
            LOG.info(f"Download gzipped index file={url}")
            path = download_and_cache(
                url,
                owner="url",
                verify=True,
                force=True,
                chunk_size=1024 * 1024,
                http_headers=None,
                update_if_out_of_date=True,
                progress_bar=no_progress_bar,
                maximum_retries=0,
                retry_after=10,
            )
        except Exception:
            LOG.error(f"Could not download index file={url}")
            raise

        return path

    def matrix_path(self, name):
        """Download (if needed) and return the local path of a matrix file.

        Parameters
        ----------
        name : str
            The matrix file's relative path within the inventory.

        Returns
        -------
        str
            The local (cached) path of the downloaded matrix file.

        Raises
        ------
        Exception
            If the file cannot be downloaded (logged before being re-raised).
        """
        try:
            url = join_url_path(self._url, name)
            path = download_and_cache(
                url,
                owner="url",
                verify=True,
                force=None,
                chunk_size=1024 * 1024,
                http_headers=None,
                update_if_out_of_date=False,
                maximum_retries=5,
                retry_after=10,
            )
        except Exception:
            LOG.error(f"Could not download matrix file={url}")
            raise

        return path


class LocalAccessor(MatrixAccessor):
    """:class:`MatrixAccessor` reading the inventory from a local directory."""

    def __init__(self, path):
        """Initialise the accessor with the inventory's local directory.

        Parameters
        ----------
        path : str
            The local directory the index and matrix files live under.
        """
        self._path = path

    def path(self):
        """str: The local directory of the inventory."""
        return self._path

    def is_local(self):
        """bool: Always True."""
        True

    def index_path(self):
        """str: The local path of the index file."""
        return os.path.join(self._path, _INDEX_FILENAME)

    def matrix_path(self, name):
        """Return the local path of a matrix file.

        Parameters
        ----------
        name : str
            The matrix file's relative path within the inventory.

        Returns
        -------
        str
            The local path of the matrix file.
        """
        return os.path.join(self._path, name)

    def reload(self, strict=False):
        """No-op: a local inventory never needs reloading."""
        pass

    def reset(self):
        """No-op: a local inventory has no cached state to reset."""
        pass


class SystemAccessor(UrlAccessor):
    """:class:`UrlAccessor` for the built-in system matrix inventory."""

    def __init__(self):
        """Initialise the accessor pointing at the system inventory URL."""
        super().__init__(_SYSTEM_URL)


class MatrixIndex(dict):
    """In-memory representation of the matrix inventory index file."""

    def load(self, path):
        """Load and parse an ``index.json`` inventory file into this dict.

        Parameters
        ----------
        path : str
            Local path of the index file to load.

        Raises
        ------
        ValueError
            If the index file's ``"version"`` does not match the expected
            :data:`VERSION`.
        """
        with open(path, "r") as f:
            index = json.load(f)
            version = index.get("version", None)
            if version != VERSION:
                raise ValueError(f"Invalid index file version: expected {VERSION}, got {version}")
            for name, entry in index["matrix"].items():
                # it is possible that the inventory is already updated with a new
                # gridspecs type, but a given earthkit-geo version is not
                # yet supporting it. In this case loading the index should not crash.
                try:
                    in_gs = _GridWrapper.from_dict(entry["input"])
                    out_gs = _GridWrapper.from_dict(entry["output"])
                    raw = entry
                    entry = dict(**entry)
                    entry["input"] = in_gs
                    entry["output"] = out_gs
                    entry["_name"] = name
                    entry["_raw"] = raw
                    self[name] = entry
                except Exception as e:
                    LOG.exception(f"Could not load matrix inventory entry {name}: {e}")

    @staticmethod
    def interpolation_method_name(item):
        """Return the interpolation method name of an inventory entry.

        Parameters
        ----------
        item : dict
            A matrix inventory entry.

        Returns
        -------
        str
            The method name: ``item["interpolation"]["method"]`` if it is a
            string, or its ``"type"`` if it is a dict.

        Raises
        ------
        ValueError
            If the method is neither a string nor a dict.
        """
        inter = item["interpolation"]
        method = inter["method"]
        if isinstance(method, str):
            return method
        if isinstance(method, dict):
            return method["type"]

        raise ValueError(f"Invalid interpolation method: {method}")

    @staticmethod
    def interpolation_method(item):
        """Return the raw ``"method"`` value (str or dict) of an inventory entry."""
        return item["interpolation"]["method"]

    @staticmethod
    def make_interpolation_uid(item):
        """Compute a unique id for an inventory entry's interpolation options.

        Parameters
        ----------
        item : dict
            A matrix inventory entry.

        Returns
        -------
        str
            The method name itself when the interpolation options are the
            default ``"grid-box-average"`` options or contain only
            ``"method"``/``"engine"``/``"version"``, otherwise a SHA-256
            digest of the interpolation options (see :func:`make_sha`).
        """
        inter = item["interpolation"]
        method = MatrixIndex.interpolation_method_name(item)
        # TODO: remove this when MIR is fixed
        if method == "grid-box-average" and is_gridbox_default(inter):
            uid = method
        elif isinstance(MatrixIndex.interpolation_method(item), dict):
            uid = make_sha(inter)
        elif set(inter.keys()) == {"method", "engine", "version"}:
            uid = method
        else:
            uid = make_sha(inter)
        return uid

    @staticmethod
    def matrix_dir_name(item):
        """Return the subdirectory name a matrix inventory entry's file is stored under.

        Parameters
        ----------
        item : dict
            A matrix inventory entry.

        Returns
        -------
        str
            ``"{engine}_{version}_{method_name}"``.
        """
        # TODO: review this logic when non-default interpolation options will
        # be available for a given method
        inter = item["interpolation"]
        engine = inter["engine"]
        version = inter["version"]
        method_name = MatrixIndex.interpolation_method_name(item)
        # uid =  inter.get("_uid", method_name)
        # uid = inter.get("_uid", inter["method"])
        return f"{engine}_{version}_{method_name}"

    @staticmethod
    def matrix_path(item):
        """Return an inventory entry's matrix file path, relative to the inventory root.

        Parameters
        ----------
        item : dict
            A matrix inventory entry.

        Returns
        -------
        str
            ``"<matrix_dir_name>/<entry name>.npz"``.
        """
        return os.path.join(MatrixIndex.matrix_dir_name(item), item["_name"] + ".npz")

    def find(self, gridspec_in, gridspec_out, method):
        """Find the inventory entry matching a grid pair and interpolation method.

        Parameters
        ----------
        gridspec_in : Any
            The input grid spec, in a form accepted by
            :meth:`~.gridspec._GridWrapper.from_any`.
        gridspec_out : Any
            The output grid spec, in the same form as ``gridspec_in``.
        method : str
            The interpolation method name.

        Returns
        -------
        dict or None
            The matching entry, or None if ``gridspec_in``/``gridspec_out``
            could not be parsed or no entry matches.
        """
        gridspec_in = _GridWrapper.from_any(gridspec_in)
        gridspec_out = _GridWrapper.from_any(gridspec_out)

        if gridspec_in is None or gridspec_out is None:
            return None

        for _, entry in self.items():
            if MatrixIndex.match(entry, gridspec_in, gridspec_out, method):
                return entry
        return None

    @staticmethod
    def match(item, gs_in, gs_out, method):
        """Check whether an inventory entry matches a grid pair and method.

        Parameters
        ----------
        item : dict
            A matrix inventory entry.
        gs_in : _GridWrapper
            The input grid to match against ``item["input"]``.
        gs_out : _GridWrapper
            The output grid to match against ``item["output"]``.
        method : str
            The interpolation method name to match.

        Returns
        -------
        bool
            True if the entry's method, input grid and output grid all
            match, False otherwise.
        """
        if (
            MatrixIndex.interpolation_method_name(item) == method
            and item["input"] == gs_in
            and item["output"] == gs_out
        ):
            return True
        return False

    @staticmethod
    def matrix_filename(item):
        """Return an inventory entry's matrix filename (without directory)."""
        return item["_name"] + ".npz"

    def subset(self, filters, fail_on_missing=True, raw=False):
        """Build a new :class:`MatrixIndex` with only the entries matching ``filters``.

        Parameters
        ----------
        filters : List[dict]
            A list of filter dicts, each with ``"input"``, ``"output"`` grid
            specs and an optional ``"method"`` (default ``"linear"``).
        fail_on_missing : bool, default=True
            If True, raise when a filter matches no entry; otherwise skip it
            with a warning.
        raw : bool, default=False
            If True, store each matched entry's raw (un-parsed) form instead
            of its parsed form.

        Returns
        -------
        MatrixIndex
            The subset index containing the matched entries.

        Raises
        ------
        ValueError
            If ``fail_on_missing`` is True and a filter matches no entry.
        """
        res = MatrixIndex()

        for i, item in enumerate(filters):
            gs_in = item["input"]
            gs_out = item["output"]
            method = item.get("method", "linear")
            LOG.info(f"ITEM[{i}]: {gs_in=} {gs_out=} {method=}")
            entry = self.find(gs_in, gs_out, method)
            if entry is not None:
                LOG.info("  found DB entry:" + entry["_name"])
                if raw:
                    res[entry["_name"]] = entry["_raw"]
                else:
                    res[entry["_name"]] = entry
            else:
                if fail_on_missing:
                    raise ValueError("No DB entry found!")
                else:
                    LOG.warning("  no DB entry found!")
        return res

    def to_raw(self):
        """Convert this index back to the raw ``index.json``-shaped dict.

        Returns
        -------
        dict
            ``{"version": VERSION, "matrix": {name: raw_entry, ...}}``.
        """
        res = dict(version=VERSION, matrix={})
        for _, entry in self.items():
            res["matrix"][entry["_name"]] = entry["_raw"]
        return res


class MatrixDb:
    """Interpolation matrix inventory for the "precomputed" regrid backend.

    Wraps a :class:`MatrixAccessor` (local directory or remote URL) and its
    parsed :class:`MatrixIndex`, and finds/loads the precomputed sparse
    interpolation matrix for a given input/output grid pair and
    interpolation method (see :meth:`find`).
    """

    # Parsing the index file builds a real eckit.geo.Grid per entry, which is
    # expensive (thousands of entries). Cache the parsed MatrixIndex per
    # (path, mtime, size), so re-loading an unchanged index file (e.g. after
    # _clear_index()) is free. Keyed the same way test_remote_index.py itself
    # detects whether an index file was actually re-downloaded.
    _INDEX_CACHE: dict = {}

    def __init__(self, accessor):
        """Initialise the matrix inventory with the given accessor.

        Parameters
        ----------
        accessor : MatrixAccessor
            Gives access to the inventory's index and matrix files.
        """
        self._index = None
        self._accessor = accessor

    @property
    def index(self):
        """MatrixIndex: The parsed inventory index, loaded lazily on first access."""
        if self._index is None:
            self._load_index()
        return self._index

    def _load_index(self):
        """Load (using the process-wide cache when possible) the accessor's index file."""
        path = self._accessor.index_path()

        st = os.stat(path)
        key = (path, st.st_mtime_ns, st.st_size)

        index = MatrixDb._INDEX_CACHE.get(key)
        if index is None:
            index = MatrixIndex()
            index.load(path)
            MatrixDb._INDEX_CACHE[key] = index

        self._index = index

    def _method_alias(self, method):
        """Return the canonical interpolation method name for a known alias.

        Parameters
        ----------
        method : str
            The interpolation method name, possibly an alias (see
            :data:`_METHOD_ALIAS`).

        Returns
        -------
        str
            The canonical method name, or ``method`` unchanged if it is not
            a known alias.
        """
        for k, v in _METHOD_ALIAS.items():
            if method in v:
                return k
        return method

    def find(
        self,
        gridspec_in,
        gridspec_out,
        method,
        **kwargs,
    ):
        """Find (loading and caching in memory if needed) a precomputed matrix.

        Parameters
        ----------
        gridspec_in : Any
            The input grid spec, in a form accepted by
            :meth:`~.gridspec._GridWrapper.from_any`.
        gridspec_out : Any
            The output grid spec, in the same form as ``gridspec_in``.
        method : str
            The interpolation method name (or a known alias, see
            :data:`_METHOD_ALIAS`).
        **kwargs : dict
            Additional keyword arguments forwarded to the in-memory cache
            lookup (``earthkit.geo.grids.utils.memcache.MEMORY_CACHE``).

        Returns
        -------
        Tuple[scipy.sparse.spmatrix or None, Tuple[int, ...] or None]
            The interpolation matrix and the output grid's shape, or
            ``(None, None)`` if the grid specs could not be parsed or no
            matching entry was found.
        """
        try:
            gridspec_in = _GridWrapper.from_any(gridspec_in)
            gridspec_out = _GridWrapper.from_any(gridspec_out)
        except Exception as e:
            LOG.warning(f"Cannot parse gridspecs with eckit.geo: {e}")
            gridspec_in = None
            gridspec_out = None

        if gridspec_in is None or gridspec_out is None:
            return None, None

        from earthkit.geo.grids.utils.memcache import MEMORY_CACHE

        return MEMORY_CACHE.get(
            gridspec_in,
            gridspec_out,
            method,
            create=self._create_matrix,
            find_entry=self.find_entry,
            create_from_entry=self._create_matrix_from_entry,
            **kwargs,
        )

    def _create_matrix(self, gridspec_in, gridspec_out, method):
        """Find the matching inventory entry and load its matrix.

        Parameters
        ----------
        gridspec_in : Any
            The input grid spec.
        gridspec_out : Any
            The output grid spec.
        method : str
            The interpolation method name.

        Returns
        -------
        Tuple[scipy.sparse.spmatrix or None, Tuple[int, ...] or None]
            See :meth:`_create_matrix_from_entry`.
        """
        return self._create_matrix_from_entry(self.find_entry(gridspec_in, gridspec_out, method))

    def _create_matrix_from_entry(self, entry):
        """Load an inventory entry's matrix and its output shape.

        Parameters
        ----------
        entry : dict or None
            A matrix inventory entry, or None.

        Returns
        -------
        Tuple[scipy.sparse.spmatrix or None, Tuple[int, ...] or None]
            The loaded matrix and the output grid's shape, or
            ``(None, None)`` if ``entry`` is None.
        """
        if entry is not None:
            z = self.load_matrix(entry)
            return z, entry["output"].shape
        return None, None

    def find_entry(self, gridspec_in, gridspec_out, method):
        """Find the inventory entry for a grid pair and method, reloading the index if needed.

        If no entry is found and the accessor is remote and has not yet
        checked for updates, forces a reload of the remote index and
        retries once.

        Parameters
        ----------
        gridspec_in : Any
            The input grid spec, already wrapped (or wrappable) via
            :meth:`~.gridspec._GridWrapper.from_any`.
        gridspec_out : Any
            The output grid spec, in the same form as ``gridspec_in``.
        method : str
            The interpolation method name (or a known alias).

        Returns
        -------
        dict or None
            The matching entry, or None if none is found.
        """
        method = self._method_alias(method)
        entry = self.index.find(gridspec_in, gridspec_out, method)
        if entry is None and not self._accessor.is_local() and not self._accessor.checked_remote():
            LOG.info(f"Matrix not found in DB for {gridspec_in=} {gridspec_out=} {method=}")
            LOG.info("Try to fetch remote index file to check for updates")
            self._accessor.reload()
            self._load_index()
            entry = self.index.find(gridspec_in, gridspec_out, method)

        return entry

    def load_matrix(self, entry):
        """Load an inventory entry's sparse interpolation matrix from disk.

        Parameters
        ----------
        entry : dict
            A matrix inventory entry.

        Returns
        -------
        scipy.sparse.spmatrix
            The loaded matrix.
        """
        path = self._matrix_fs_path(entry)
        z = load_npz(path)
        return z

    def _matrix_index_filename(self, entry):
        """Return an inventory entry's matrix filename (without directory)."""
        return self.index.matrix_filename(entry)

    def _matrix_index_path(self, entry):
        """Return an inventory entry's matrix path, relative to the inventory root."""
        return self.index.matrix_path(entry)

    def _matrix_fs_path(self, entry):
        """Return an inventory entry's local filesystem matrix path."""
        return self._accessor.matrix_path(self._matrix_index_path(entry))

    def subset_index(self, filters, **kwargs):
        """Build a subset :class:`MatrixIndex` matching the given filters.

        Parameters
        ----------
        filters : List[dict]
            See :meth:`MatrixIndex.subset`.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :meth:`MatrixIndex.subset`.

        Returns
        -------
        MatrixIndex
            The subset index.
        """
        return self.index.subset(filters, **kwargs)

    def copy_matrix_file(self, entry, out_dir, exist_ok=False, dry_run=False):
        """Copy an inventory entry's matrix file into a local directory.

        Parameters
        ----------
        entry : dict
            A matrix inventory entry.
        out_dir : str
            The destination directory root (the entry's subdirectory
            structure is preserved under it).
        exist_ok : bool, default=False
            If False, raise when the target file already exists.
        dry_run : bool, default=False
            If True, do not actually copy the file (only compute and return
            the target path, logging a warning instead of raising if it
            already exists).

        Returns
        -------
        str
            The target file path.

        Raises
        ------
        FileExistsError
            If the target file already exists, ``exist_ok`` is False and
            ``dry_run`` is False.
        """
        import shutil

        matrix_index_path = self._matrix_index_path(entry)
        src_file = self._matrix_fs_path(entry)
        target_file = os.path.join(out_dir, matrix_index_path)

        if not exist_ok and os.path.exists(target_file):
            if not dry_run:
                raise FileExistsError(f"target file already exists! {target_file}")
            else:
                LOG.warning("target file already exists! {target_file}")

        target_dir = os.path.dirname(target_file)

        if not dry_run:
            os.makedirs(target_dir, exist_ok=True)
            shutil.copyfile(src_file, target_file)

        return target_file

    def index_file_path(self):
        """str: The local path of the inventory's index file."""
        return self._accessor.index_path()

    def matrix_source(self):
        """str: The inventory's location (a local path or URL)."""
        return self._accessor.path()

    @staticmethod
    def from_path(path):
        """Build a :class:`MatrixDb` backed by a local directory.

        Parameters
        ----------
        path : str
            The local directory the inventory is stored in.

        Returns
        -------
        MatrixDb
        """
        return MatrixDb(LocalAccessor(path))

    @staticmethod
    def from_url(url):
        """Build a :class:`MatrixDb` backed by a remote URL.

        Parameters
        ----------
        url : str
            The base URL the inventory is served from.

        Returns
        -------
        MatrixDb
        """
        return MatrixDb(UrlAccessor(url))

    def __len__(self):
        """Return the number of entries in the inventory index."""
        return len(self.index)

    def _clear_index(self):
        """For testing only."""
        self._index = None

    def _reset(self):
        """For testing only."""
        self._index = None
        self._accessor.reset()


SYS_DB = MatrixDb(SystemAccessor())
