:noindex:
:orphan:

.. _api-reference:

API reference
=============

This page holds poliastro's API documentation, which might be helpful for final
users or developers to create their own poliastro-based utilities. Among the
different sub-packages and modules, we might differentiate two big categories:
core utilities and high-level ones.

* **Core API:** These routines are located in :py:mod:`poliastro.core`
  and can be considered raw Orbital Mechanics algorithms.
  These kind of routines are meant for developers and advanced users,
  since they provide low-level utilities.

* **High-level API:** All logical files included in this part define the main
  poliastro objects such us :py:class:`poliastro.twobody.Orbit`,
  :py:class:`poliastro.maneuver.Maneuver`, and their interaction with each other.
  In their definitions, they make use of previously presented core routines
  and extend them in such a way that
  it is easy for users to solve classical Orbital Mechanics problems.

.. toctree::
   :maxdepth: 5

   {% for page in pages %}
   {% if page.top_level_object and page.display %}
   {{ page.include_path }}
   {% endif %}
   {% endfor %}
