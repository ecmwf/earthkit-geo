.. _caching:

Disk-based caching of regridding precomputed weights caching
===============================================================

*New in version 1.0.0.*


.. note::

    At the moment this configuration is only related to the :ref:`precomputed <precomputed-backend>`
    backends in :func:`regrid`.


Purpose
-------

earthkit-geo uses a dedicated **directory** to store interpolation matrices and the related index file downloaded from the remote inventory. By default this directory serves a **cache** and is **managed** (its size is checked/limited).  It means if we run :func:`regrid` again with the same input and output grid it will load the matrix from the cache instead of downloading it again. Additionally, caching offers **monitoring and disk space management**. When the cache is full, cached data is deleted according to the configuration (i.e. oldest data is deleted first). The cache is implemented by using a sqlite database running in a separate thread.

Please note that the earthkit-geo cache configuration is managed through the :doc:`config`.

.. warning::

    The earthkit-geo cache is intended to be used by a single user.
    Sharing cache with **multiple users is not recommended**.
    Downloading a local copy of data on a shared disk to have multiple
    users working is a different use case and should be supported
    through using mirrors.

.. _cache_location:
.. _cache_policies:

Cache policies
------------------------------

The primary config option to control the cache is ``cache-policy``, which can take the following values:

  - :ref:`user <user_cache_policy>`  (default)
  - :ref:`temporary <temporary_cache_policy>`
  - :ref:`off <off_cache_policy>`

The cache location can be read and modified with Python (see the details below).

.. tip::

   See the :ref:`/how-tos/misc/cache.ipynb` notebook for examples.

.. note::

  It is recommended to restart your Jupyter kernels after changing
  the cache policy or location.

.. _user_cache_policy:

User cache policy
+++++++++++++++++++

When the ``cache-policy`` is "user" the **cache will be active** and created in a **managed directory** defined by the ``user-cache-directory`` config option. This is the **default** value.

.. note::

    The default location of the user cache directory is ``"~/.cache/earthkit-geo"`` and its maximum size is 5 GB.

The user cache directory is **not cleaned up on exit**. So next time you start earthkit-geo it will be there again unless it is deleted manually or it is set in way that on each startup a different path is assigned to it. Also, when you run multiple sessions of earthkit-geo under the same user they will share the same cache.

We can query the directory path via the :doc:`config` and also by calling the :meth:`~earthkit.geo.utils.caching.Cache.directory` :ref:`cache method <cache_methods>`.

.. code-block:: python

  >>> from earthkit.geo import cache, config
  >>> config.set("cache-policy", "user")
  >>> config.get("user-cache-directory")
  '/Users/username/.cache/earthkit-geo'
  >>> cache.directory()
  '/Users/username/.cache/earthkit-geo'


The following code shows how to change the ``user-cache-directory`` config option:

.. code:: python

  >>> from earthkit.geo import config
  >>> config.get("user-cache-directory")  # Find the current cache directory
  '/Users/username/.cache/earthkit-geo'
  >>> # Change the value of the setting
  >>> config.set("user-cache-directory", "/big-disk/earthkit-geo-cache")

  # Python kernel restarted

  >>> from earthkit.geo import config
  >>> config.get("user-cache-directory")  # Cache directory has been modified
  '/big-disk/earthkit-geo-cache'

More generally, the earthkit-geo config options can be read, modified, reset
to their default values from Python,
see the :doc:`Configs documentation <config>`.

.. _temporary_cache_policy:

Temporary cache policy
++++++++++++++++++++++++

When the ``cache-policy`` is "temporary" the **cache will be active and located in a managed** temporary directory created by ``tempfile.TemporaryDirectory``. This directory will be unique for each earthkit-geo session. When the directory object goes out of scope (at the latest on exit) the cache is **cleaned up**.

Due to the temporary nature of this directory path it cannot be queried via the :doc:`config`, but we need to call the :meth:`~earthkit.geo.utils.caching.Cache.directory` :ref:`cache method <cache_methods>`.

.. code-block:: python

  >>> from earthkit.geo import cache, config
  >>> config.set("cache-policy", "temporary")
  >>> cache.directory()
  '/var/folders/ng/g0zkhc2s42xbslpsywwp_26m0000gn/T/tmp_5bf5kq8'

We can specify the parent directory for the the temporary cache by using the ``temporary-cache-directory-root`` config option. By default it is set to Non  e (no parent directory specified).

.. code-block:: python

  >>> from earthkit.geo import cache, config
  >>> s = {
  ...     "cache-policy": "temporary",
  ...     "temporary-cache-directory-root": "~/my_demo_cache",
  ... }
  >>> config.set(s)
  >>> cache.directory()
  '~/my_demo_cache/tmp0iiuvsz5'


.. _off_cache_policy:

Off cache policy
++++++++++++++++++++++++

When the ``cache-policy`` is "off" no disk-based caching is available. In this case all files are downloaded into an **unmanaged** temporary directory created by ``tempfile.TemporaryDirectory``. Since caching is disabled, all repeated calls to :func:`regrid` will download the interpolation matrix again! This temporary directory will be unique for each earthkit-geo session. When the directory object goes out of scope (at the latest on exit) the directory will be **cleaned up**.

Due to the temporary nature of this directory path it cannot be queried via the :doc:`config`, but we need to call the :meth:`~earthkit.geo.utils.caching.Cache.directory` :ref:`cache method <cache_methods>`.

.. code-block:: python

  >>> from earthkit.geo import cache, config
  >>> config.set("cache-policy", "off")
  >>> cache.directory()
  '/var/folders/ng/g0zkhc2s42xbslpsywwp_26m0000gn/T/tmp_5bf5kq8'

We can specify the parent directory for the the temporary directory by using the ``temporary-directory-root`` config. By default it is set to None (no parent directory specified).

.. code-block:: python

  >>> from earthkit.geo import cache, config
  >>> s = {
  ...     "cache-policy": "off",
  ...     "temporary-directory-root": "~/my_demo_tmp",
  ... }
  >>> config.set(s)
  >>> cache.directory()
  '~/my_demo_tmp/tmp0iiuvsz5'


.. _cache_object:
.. _cache_methods:

Cache methods
-------------------------

The cache is controlled by a global object, which we can access as ``earthkit.geo.cache``.

.. code:: python

  >>> from earthkit.geo import cache
  >>> cache
  <earthkit.geo.utils.caching.Cache object at 0x117be7040>


When ``cache-policy`` is :ref:`user <user_cache_policy>` or :ref:`temporary <temporary_cache_policy>`
there are a set of methods available on this object to manage and interact with the cache.

.. list-table:: Methods/properties of the cache object
   :header-rows: 1

   * - Methods
     - Description

   * - :attr:`~earthkit.geo.utils.caching.Cache.policy`
     - Get the current cache policy object.
   * - :meth:`~earthkit.geo.caching.Cache.directory`
     - Return the path to the current cache directory
   * - :meth:`~earthkit.geo.utils.caching.Cache.size`
     - Return the total number of bytes stored in the cache
   * - :meth:`~earthkit.geo.utils.caching.Cache.check_size`
     - Check the cache size and trim it down when needed.
   * - :meth:`~earthkit.geo.utils.caching.Cache.entries`
     - Dump the entries stored in the cache
   * - :meth:`~earthkit.geo.utils.caching.Cache.summary_dump_database`
     - Return the number of items and total size of the cache
   * - :meth:`~earthkit.geo.utils.caching.Cache.purge`
     - Delete entries from the cache

.. warning::

    :meth:`~earthkit.geo.utils.caching.Cache.check_size` automatically runs when a new
    entry is added to the cache or any of the :ref:`cache_config` changes.

Examples:

.. code:: python

      >>> from earthkit.geo import cache
      >>> cache.policy.name
      'user'
      >>> cache.directory()
      '/Users/username/.cache/earthkit-geo/''
      >>> cache.size()
      846785699
      >>> cache.summary_dump_database()
      (40, 846785699)
      >>> d = cache.entries()
      >>> len(d)
      40
      >>> d[0].get("creation_date")
      '2023-10-30 14:48:31.320322'


Cache limits
------------

.. warning::

  These config options do not work when ``cache-policy`` is :ref:`off <off_cache_policy>`.



Maximum-cache-size
  The ``maximum-cache-size`` setting ensures that earthkit-geo does not
  use to much disk space.  Its value sets
  the maximum disk space used by earthkit-geo cache.  When earthkit-geo cache disk
  usage goes above this limit, earthkit-geo triggers its cache cleaning mechanism  before
  downloading additional data.  The value of cache-maximum-size is
  absolute (such as "10G", "10M", "1K"). To disable it use None.

Maximum-cache-disk-usage
  The ``maximum-cache-disk-usage`` setting ensures that earthkit-geo
  does not fill your disk. It specifies the maximum disk usage (as a percentage) on the filesystem
  containing the cache directory. When the total disk usage (so this is not the cache usage alone) goes above
  this limit, earthkit-geo triggers its cache cleaning mechanism to free up space before
  downloading additional data.
  The value of maximum-cache-disk-usage is relative (such as "90%" or "100%").
  To disable it use None.

.. warning::
    If your disk is filled by another application, earthkit-geo will happily
    delete its cached data to make room for the other application as soon
    as it has a chance.


.. .. note::
..     When tweaking the cache config, it is recommended to set the
..     ``maximum-cache-size`` to a value below the user disk quota (if applicable)
..     and ``maximum-cache-disk-usage`` to ``None``.


.. _cache_config:

Cache config parameters
-------------------------------

.. module-output:: generate_config_rst cache-policy maximum-cache-disk-usage maximum-cache-size temporary-cache-directory-root user-cache-directory

Other earthkit-geo config options can be found :ref:`here <config_table>`.
