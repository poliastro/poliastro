User guide
==========

Defining the orbit: :code:`State` objects
-----------------------------------------

The core of poliastro are the :py:class:`~poliastro.twobody.State` objects
inside the :py:mod:`poliastro.twobody` module. They store all the required
information to define an orbit:

* The body acting as the central body of the orbit, for example the Earth.
* The position and velocity vectors or the orbital elements.
* The time at which the orbit is defined.

First of all, we have to import the relevant modules and classes:

.. code-block:: python

    import numpy as np
    from astropy import units as u
    
    from poliastro.bodies import Earth, Sun
    from poliastro.twobody import State

From position and velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several methods available to create
:py:class:`~poliastro.twobody.State` objects. For example, if we have the
position and velocity vectors we can use
:py:meth:`~poliastro.twobody.State.from_vectors`:

.. code-block:: python

    # Data from Curtis, example 4.3
    r = [-6045, -3490, 2500] * u.km
    v = [-3.457, 6.618, 2.533] * u.km / u.s
    
    ss = State.from_vectors(Earth, r, v)

And that's it! Notice a couple of things:

* Defining vectorial physical quantities using Astropy units is terribly easy.
  The list is automatically converted to a :code:`Quantity`, which is actually
  a subclass of NumPy arrays.
* If no time is specified, then a default value is assigned::

    >>> ss.epoch
    <Time object: scale='utc' format='jyear_str' value=J2000.000>
    >>> ss.epoch.iso
    '2000-01-01 12:00:00.000'

.. figure:: _static/curtis.png
   :align: right
   :figwidth: 350
   :alt: Plot of the orbit

If we're working on interactive mode (for example, using the wonderful IPython
notebook) we can immediately plot the current state::

    from poliastro.plotting import plot
    plot(ss)

This plot is made in the so called *perifocal frame*, which means:

* we're visualizing the plane of the orbit itself,
* the \\(x\\) axis points to the pericenter, and
* the \\(y\\) axis is turned \\(90 \\mathrm{^\\circ}\\) in the
  direction of the orbit.

The dotted line represents the *osculating orbit*:
the instantaneous Keplerian orbit at that point. This is relevant in the
context of perturbations, when the object shall deviate from its Keplerian
orbit.

From classical orbital elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can also define a :py:class:`~poliastro.twobody.State` using a set of
six parameters called orbital elements. Although there are several of
these element sets, each one with its advantages and drawbacks, right now
poliastro supports the *classical orbital elements*:

* Semimajor axis \\(a\\).
* Eccentricity \\(e\\).
* Inclination \\(i\\).
* Right ascension of the ascending node \\(\\Omega\\).
* Argument of pericenter \\(\\omega\\).
* True anomaly \\(\\nu\\).

In this case, we'd use the method
:py:meth:`~poliastro.twobody.State.from_classical`:

.. code-block:: python

    # Data for Mars at J2000 from JPL HORIZONS
    a = 1.523679 * u.AU
    ecc = 0.093315 * u.one
    inc = 1.85 * u.deg
    raan = 49.562 * u.deg
    argp = 286.537 * u.deg
    nu = 23.33 * u.deg
    
    ss = State.from_classical(Sun, a, ecc, inc, raan, argp, nu)

Notice that whether we create a ``State`` from \\(r\\) and \\(v\\) or from
elements we can access many mathematical properties individually::

    >>> ss.period.to(u.day)
    <Quantity 686.9713888628166 d>
    >>> ss.v
    <Quantity [  1.16420211, 26.29603612,  0.52229379] km / s>

To see a complete list of properties, check out the
:py:class:`poliastro.twobody.State` class on the API reference.

Changing the orbit: :code:`Maneuver` objects
--------------------------------------------

poliastro helps us define several in-plane and general out-of-plane
maneuvers with the :py:class:`~poliastro.maneuver.Maneuver` class inside the
:py:mod:`poliastro.maneuver` module.

Each ``Maneuver`` consists on a list of impulses \\(\\Delta v_i\\)
(changes in velocity) each one applied at a certain instant \\(t_i\\). The
simplest maneuver is a single change of velocity without delay: you can
recreate it either using the :py:meth:`~poliastro.maneuver.Maneuver.impulse`
method or instantiating it directly.

.. code-block:: python

    dv = [5, 0, 0] * u.m / u.s
    
    man = Maneuver.impulse(dv)
    man = Maneuver((0 * u.s, dv))  # Equivalent

There are other useful methods you can use to compute common in-plane
maneuvers, notably :py:meth:`~poliastro.maneuver.Maneuver.hohmann` and
:py:meth:`~poliastro.maneuver.Maneuver.bielliptic` for `Hohmann`_ and
`bielliptic`_ transfers respectively. Both return the corresponding
``Maneuver`` object, which in turn you can use to calculate the total cost
in terms of velocity change (\\(\\sum \|\\Delta v_i|\\)) and the transfer
time::

    >>> ss_i = State.circular(Earth, alt=700 * u.km)
    >>> hoh = Maneuver.hohmann(ss_i, 36000 * u.km)
    >>> hoh.get_total_cost()
    <Quantity 3.6173981270031357 km / s>
    >>> hoh.get_total_time()
    <Quantity 15729.741535747102 s>

You can also retrieve the individual vectorial impulses::

    >>> hoh.impulses[0]
    (<Quantity 0 s>, <Quantity [ 0.        , 2.19739818, 0.        ] km / s>)
    >>> hoh[0]  # Equivalent
    (<Quantity 0 s>, <Quantity [ 0.        , 2.19739818, 0.        ] km / s>)
    >>> tuple(_.decompose([u.km, u.s]) for _ in hoh[1])
    (<Quantity 15729.741535747102 s>, <Quantity [ 0.        , 1.41999995, 0.        ] km / s>)

.. _Hohmann: http://en.wikipedia.org/wiki/Hohmann_transfer_orbit
.. _bielliptic: http://en.wikipedia.org/wiki/Bi-elliptic_transfer

To actually retrieve the resulting ``State`` after performing a maneuver, use
the method :py:meth:`apply_maneuver`::

    >>> ss_f = ss_i.apply_maneuver(hoh)
    >>> ss_f.rv()
    (<Quantity [ -3.60000000e+04, -7.05890200e-11, -0.00000000e+00] km>, <Quantity [ -8.97717523e-16, -3.32749489e+00, -0.00000000e+00] km / s>)

More advanced plotting: :code:`OrbitPlotter` objects
----------------------------------------------------

We previously saw the :py:func:`poliastro.plotting.plot` function to easily
plot orbits. Now we'd like to plot several orbits in one graph (for example,
the maneuver me computed in the previous section). For this purpose, we
have :py:class:`~poliastro.plotting.OrbitPlotter` objects in the
:py:mod:`~poliastro.plotting` module.

These objects hold the perifocal plane of the first ``State`` we plot in
them, projecting any further trajectories on this plane. This allows to
easily visualize in two dimensions:

.. code-block:: python

    from poliastro.plotting import OrbitPlotter
    
    op = OrbitPlotter()
    ss_a, ss_f = ss_i.apply_maneuver(hoh, intermediate=True)
    op.plot(ss_i, label="Initial orbit")
    op.plot(ss_a, label="Transfer orbit")
    op.plot(ss_f, label="Final orbit")

Which produces this beautiful plot:

.. figure:: _static/hohmann.png
   :align: center
   :alt: Hohmann transfer
   
   Plot of a Hohmann transfer.

Where are the planets? Computing ephemerides
--------------------------------------------

.. versionadded:: 0.3.0

Thanks to the awesome jplephem package, poliastro can now read Satellite
Planet Kernel (SPK) files, part of NASA's SPICE toolkit. This means that
we can query the position and velocity of the planets of the Solar System.

The first time we import :py:mod:`poliastro.ephem` we will get a warning
indicating that no SPK files are present::

    >>> import poliastro.ephem
    No SPICE kernels found under ~/.poliastro. Please download them manually or using

      poliastro download-spk [-d NAME]

    to provide a default kernel, else pass a custom one as an argument to `planet_ephem`.

This is because poliastro does not download any data when installed: SPK files
weight several MiB and that would slow the download process. Instead, we are
requested to download them from NASA website or use the builtin command-line
utility::

    $ poliastro download-spk --name de421
    No SPICE kernels found under ~/.poliastro. Please download them manually or using

      poliastro download-spk [-d NAME]

    to provide a default kernel, else pass a custom one as an argument to `planet_ephem`.
    Downloading de421.bsp from http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/, please wait...
    Not Found
    Downloading de421.bsp from http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/, please wait...

If no ``--name`` argument is provided, de430 will be downloaded.
Alternatively, we can use :py:func:`poliastro.ephem.download_kernel` from a
Python session::

    >>> from poliastro import ephem
    >>> ephem.download_kernel("de421")
    File de421.bsp already exists under /home/juanlu/.poliastro
    >>>

In this case, the ``name`` argument is required.

Once we have downloaded an SPK file we can already compute the position and
velocity vectors of the planets with the
:py:func:`poliastro.ephem.planet_ephem` function. All we need is the body
we are querying and an ``astropy.time.Time`` scalar or vector variable::

    >>> from astropy import time
    >>> epoch = time.Time("2015-05-09 10:43")
    >>> from poliastro import ephem
    >>> r, v = ephem.planet_ephem(ephem.EARTH, epoch)
    >>> r
    <Quantity [ -9.99802065e+07, -1.03447226e+08, -4.48696791e+07] km>
    >>> v
    <Quantity [ 1880007.6848216 ,-1579126.15900176, -684591.24441181] km / s>

.. note:: The position and velocity vectors are usually given with respect to the
    Solar System Barycenter in the **International Celestial Reference Frame**
    (ICRF), which means approximately equatorial coordinates.

Traveling through space: solving the Lambert problem
----------------------------------------------------

The determination of an orbit given two position vectors and the time of
flight is known in celestial mechanics as **Lambert's problem**, also
known as two point boundary value problem. This contrasts with Kepler's
problem or propagation, which is rather an initial value problem.

The module :py:mod:`poliastro.iod` allows as to solve Lambert's problem,
provided the main attractor's gravitational constant, the two position
vectors and the time of flight. As you can imagine, being able to compute
the positions of the planets as we saw in the previous section is the
perfect complement to this feature!

For instance, this is a simplified version of the example
`Going to Mars with Python using poliastro`_, where the orbit of the
Mars Science Laboratory mission (rover Curiosity) is determined::

    >>> from astropy import time
    >>> date_launch = time.Time('2011-11-26 15:02', scale='utc')
    >>> date_arrival = time.Time('2012-08-06 05:17', scale='utc')
    >>> tof = date_arrival - date_launch
    >>> from poliastro import ephem
    >>> r0, _ = ephem.planet_ephem(ephem.EARTH, date_launch)
    >>> r, _ = ephem.planet_ephem(ephem.MARS, date_arrival)
    >>> from poliastro import iod
    >>> from poliastro.bodies import Sun
    >>> (v0, v), = iod.lambert(Sun.k, r0, r, tof)
    >>> v0
    <Quantity [-29.29150998, 14.53326521,  5.41691336] km / s>
    >>> v
    <Quantity [ 17.6154992 ,-10.99830723, -4.20796062] km / s>


.. figure:: _static/msl.png
   :align: center
   :alt: MSL orbit

   Mars Science Laboratory orbit.

.. _`Going to Mars with Python using poliastro`: http://nbviewer.ipython.org/github/poliastro/poliastro/blob/master/examples/Going%20to%20Mars%20with%20Python%20using%20poliastro.ipynb

*Per Python ad astra* ;)
