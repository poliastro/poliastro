# Related software

These are some projects which share similarities with poliastro or which
served as inspiration:

- [astropy](http://www.astropy.org/): According to its website, \"The
  Astropy Project is a community effort to develop a single core
  package for Astronomy in Python and foster interoperability between
  Python astronomy packages\". Not only does it provide important core
  features for poliastro like [time](https://docs.astropy.org/en/stable/time/) and physical [units](https://docs.astropy.org/en/stable/units/) handling, but
  also sets a [high bar](https://docs.astropy.org/en/stable/index.html) for code quality and documentation standards. A
  truly inspiring project.
- [Skyfield](https://rhodesmill.org/skyfield/): Another Astronomy
  Python package focused on computing observations of planetary bodies
  and Earth satellites written by [Brandon Rhodes](https://rhodesmill.org/brandon/). It is the successor
  of [pyephem](https://rhodesmill.org/pyephem/), also written by him, but skyfield is a pure Python
  package and provides a much cleaner API.
- [Plyades](https://plyades.readthedocs.io/): A pioneering
  astrodynamics library written in Python by [Helgee Eichhorn](https://helgeeichhorn.de/). Its
  clean and user-friendly API inspired me to completely refactor
  poliastro 0.2 so it could be much easier to use. It is now deprecated by the author, with [Astrodynamics.jl](https://juliaastrodynamics.github.io/) being its successor (poliastro, too!)
- [orbital](https://pythonhosted.org/OrbitalPy/): Yet another orbital
  mechanics Python library written by [Frazer McLean](https://www.frazermclean.co.uk/). It is very
  similar to poliastro (orbital plotting module was inspired by mine)
  but its internal structure is way smarter. It is more focused in
  plotting and it even provides 3D plots and animations.
- [orekit-python-wrapper](https://www.orekit.org/forge/projects/orekit-python-wrapper/wiki):
  According to its website, \"The Orekit python wrapper enables to use
  Orekit within a normal python environment\", using [JCC](https://lucene.apache.org/pylucene/jcc/index.html). [Orekit](https://www.orekit.org/) is a
  well-stablished, mature open source library for Astrodynamics
  written in Java strongly supported by several space agencies. The
  Python wrapper is developed by the [Swedish Space Corporation](https://sscspace.com/).
- [beyond](https://github.com/galactics/beyond/): A young flight
  dynamics library written in Python with a focus on developing \"a
  simple API for space observations\". Some parts overlap with
  poliastro, but it also introduces many interesting features, and the
  examples look promising. Worth checking!
- [SpiceyPy](https://github.com/andrewannex/SpiceyPy): This Python
  library wraps the [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html), a huge software collection
  developed by NASA which offers advanced astrodynamics functionality.
  Among all the wrappers available on the Internet, at the time of
  writing this is the most advanced and well-maintained one, although
  there are others.