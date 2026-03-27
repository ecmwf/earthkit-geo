
Earthkit-geo's documentation
======================================================

|Static Badge| |image1| |License: Apache 2.0| |Latest Release|

.. |Static Badge| image:: https://github.com/ecmwf/codex/raw/refs/heads/main/ESEE/foundation_badge.svg
   :target: https://github.com/ecmwf/codex/raw/refs/heads/main/ESEE
.. |image1| image:: https://github.com/ecmwf/codex/raw/refs/heads/main/Project%20Maturity/emerging_badge.svg
   :target: https://github.com/ecmwf/codex/raw/refs/heads/main/Project%20Maturity
.. |License: Apache 2.0| image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
   :target: https://opensource.org/licenses/apache-2-0
.. |Latest Release| image:: https://img.shields.io/github/v/release/ecmwf/earthkit-geo?color=blue&label=Release&style=flat-square
   :target: https://github.com/ecmwf/earthkit-geo/releases


.. important::

    This software is **Emerging** and subject to ECMWF's guidelines on `Software Maturity <https://github.com/ecmwf/codex/raw/refs/heads/main/Project%20Maturity>`_.


**earthkit-geo** is a Python package providing geospatial computations. It is part of the :xref:`earthkit` ecosystem.

The API features the :func:`regrid` function taking inputs of Numpy arrays, earthkit-data GRIB :xref:`field` or :xref:`fieldlist` objects, and  Xarray DataArrays or Datasets. It is implemented with various backends, the :ref:`default backend <mir-backend>` uses ECMWF's **MIR (Meteorological Interpolation and Regridding)** library to perform the regridding.



.. grid:: 1
   :gutter: 2

   .. grid-item-card:: Installation and Getting Started
      :img-top: _static/rocket.svg
      :link: getting-started
      :link-type: doc
      :class-card: sd-shadow-sm

      New to earthkit-geo? Start here with installation and a quick overview.

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item-card:: Tutorials
      :img-top: _static/book.svg
      :link: tutorials/index
      :link-type: doc
      :class-card: sd-shadow-sm

      Step-by-step guides to learn earthkit-geo.

   .. grid-item-card:: How-tos
      :img-top: _static/tool.svg
      :link: how-tos/index
      :link-type: doc
      :class-card: sd-shadow-sm

      Practical recipes for common tasks.

   .. grid-item-card:: Concepts and Explanations
      :img-top: _static/bulb.svg
      :link: explanations/index
      :link-type: doc
      :class-card: sd-shadow-sm

      Understand the core ideas behind earthkit-geo.

   .. grid-item-card:: API Reference Guide
      :img-top: _static/brackets-contain.svg
      :link: autoapi/earthkit/geo/index
      :link-type: doc
      :class-card: sd-shadow-sm

      Detailed documentation of all functions and classes.


.. toctree::
   :caption: User guide
   :maxdepth: 2
   :hidden:

   getting-started
   installation
   tutorials/index
   how-tos/index
   explanations/index
   API Reference guide <autoapi/earthkit/geo/index.rst>

.. toctree::
   :caption: Developer guide
   :maxdepth: 2
   :hidden:

   development


.. toctree::
   :maxdepth: 2
   :caption: Extras
   :hidden:

   release-notes/index
   licence
   genindex

