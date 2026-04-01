.. _install:

Installation
============

Installing from PyPI
------------------------------------

Default installation
+++++++++++++++++++++++++

Install **earthkit-geo** with python3 (>= 3.10) and ``pip`` as follows:

.. code-block:: bash

    python3 -m pip install earthkit-geo>=1.0.0rc0

This will not install the optional dependencies.

Installing all the optional packages
++++++++++++++++++++++++++++++++++++++++

You can install **earthkit-geo** with all the optional packages in one go by using:

.. code-block:: bash

    python3 -m pip install earthkit-geo[all]>=1.0.0rc0

Please note in **zsh** you need to use quotes around the square brackets:

.. code-block:: bash

    python3 -m pip install "earthkit-geo[all]>=1.0.0rc0"


Installing individual optional packages
+++++++++++++++++++++++++++++++++++++++++

Alternatively, you can install the following components (on top of the default installation) individually:

    - cartography: enables the usage of the :py:meth:`earthkit.geo.cartography.country_polygons` method

Usage:

.. code-block:: bash

    python3 -m pip install "earthkit-geo[cartography]>=1.0.0rc0"
