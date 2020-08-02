:noindex:
:orphan:

API reference
=============

This page holds poliastro's API documentation, which might be helpful for final
users or developers to create their own poliastro-based utilities. Among the
different sub-packages and modules, we might differentiate two big categories:
core utilities and high-level ones.


* **Core API:** This routines are locate within the `poliastro.core` sub-package
  and could be understood as raw orbital mechanics algorithms. These kind of
  routines are not supposed to be used by poliastro's users but rather by
  developers, since they provide low-level utilities.


* **High-level API:** All logical files included in this part define main
  `poliastro` objects such us `poliastro.twobody.Orbit`,
  `poliastro.maneuver.Maneuver` and how do they interact with each other. In
  their definitions, they make use of previously presented core routines and
  extend them in such a way that is easy for users to solve classical orbital
  mechanics problems.


.. toctree::
   :maxdepth: 5

   {% for page in pages %}
   {% if page.top_level_object and page.display %}
   {{ page.include_path }}
   {% endif %}
   {% endfor %}

