poliastro - Astrodynamics in Python
===================================

.. image:: _static/logo_text.png
   :width: 675px
   :align: center

**poliastro** is an open source (MIT) collection of Python functions useful
in Astrodynamics and Orbital Mechanics, focusing on interplanetary applications.
It provides a simple and intuitive API and handles physical quantities with
units.

View `source code`_ of poliastro!

.. _`source code`: https://github.com/poliastro/poliastro

Some of its awesome features are:

* Analytical and numerical orbit propagation
* Conversion between position and velocity vectors and classical orbital
  elements
* Coordinate frame transformations
* Hohmann and bi-elliptic maneuvers computation
* Trajectory plotting: static and interactive
* Initial orbit determination (Lambert's problem)
* Planetary ephemerides (using SPICE kernels via Astropy)
* Computation of Near-Earth Objects (NEOs)
* Atmospheric models
* Mission design tools: porkchops generator

And more to come!

poliastro is developed by an open, international community. Release
announcements and general discussion take place on our `mailing list`_
and `chat`_.

.. _`mailing list`: https://groups.io/g/poliastro-dev
.. _`chat`: https://riot.im/app/#/room/#poliastro:matrix.org

.. include:: form.rst

The `source code`_, `issue tracker`_ and `wiki`_ are hosted on GitHub, and all
contributions and feedback are more than welcome. You can test poliastro in your
browser using binder, a cloud Jupyter notebook server:

.. image:: https://img.shields.io/badge/launch-binder-e66581.svg?style=flat-square
   :target: https://mybinder.org/v2/gh/poliastro/poliastro/master?filepath=index.ipynb

.. _`source code`: https://github.com/poliastro/poliastro
.. _`issue tracker`: https://github.com/poliastro/poliastro/issues
.. _`wiki`: https://github.com/poliastro/poliastro/wiki/

See `benchmarks`_ for the performance analysis of poliastro.

.. _`benchmarks`: https://blog.poliastro.space/poliastro-benchmarks/

poliastro works on recent versions of Python and is released under
the MIT license, hence allowing commercial use of the library.

One the greatest features of poliastro apart from its accuracy, is its set of
plotting utilities: static or interactive 3D plotters, porkchop generator...
Check out a simple transfer orbit between Earth and Mars. For more information
on how to use this package, refer to the examples gallery.

.. raw:: html
    :file: auto_examples/images/sphx_glr_plot_going_to_mars_with_python_using_poliastro_001.html

.. include:: success.rst

Contents
--------

.. toctree::
    :maxdepth: 2

    about
    getting_started
    user_guide
    auto_examples/index
    api/safe/safe_index
    api/core/core_index
    changelog
    contributing
    references


.. note::
    Older versions of poliastro relied on some Fortran subroutines written by David A. Vallado for
    his book "Fundamentals of Astrodynamics and Applications" and available on
    the Internet as the `companion software of the book`__.
    The author explicitly gave permission to redistribute these subroutines
    in this project under a permissive license.

.. __: http://celestrak.com/software/vallado-sw.asp
