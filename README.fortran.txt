These are some Fortran subroutines required by the program. They come from
the companion Fortran subroutines of Vallado, available from
http://celestrak.com/software/vallado-sw.asp, but were slightly modified
due to some errors and inconsistencies.

To create the required shared libraries, these components must be present on
the system:

* `f2py`
* A Fortran compiler (`gfortran` assumed by default)

The `f2py` and Fortran compiler binaries as well as the flags for the compiler
can be set up in the Makefile, lines 2 to 6. Then, just open a terminal and
run::

  $ make

The tests are in the same folder instead of a separate one for the sake of
simplicity. To compile them, just run::

  $ make test

To remove the leftover files::

  $ make clean
