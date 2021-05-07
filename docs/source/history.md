# History

## Creation

I started poliastro as a wrapper of some MATLAB and Fortran algorithms
that I needed for a University project: having good performance was a
must, so pure Python was not an option. As a three language project, it
was only known to work in my computer, and I had to fight against [oct2py](https://pypi.org/project/oct2py/)
and [f2py](https://numpy.org/doc/stable/f2py/) for long hours.

Later on, I enhanced the plotting capabilities of poliastro to serve me in
additional university tasks. I removed the MATLAB (Octave) code and kept
only the Fortran algorithms. Finally, when numba was mature enough, I
implemented everything in pure Python and poliastro 0.3 was born.

## Future ideas

These are some things that I would love to implement in poliastro to
expand its capabilities:

- 3D plotting of orbits
- Continuous thrust maneuvers
- Tisserand graphs
- Porkchop plots

## Acknowledgement from the original author

I am Juan Luis Cano Rodríguez (two names and two surnames, it\'s the
Spanish way!), an Aerospace Engineer with a passion for Astrodynamics
and the Open Source world. Before poliastro started to be a truly
community project, I started it when I was an Erasmus student at
Politecnico di Milano, an important technical university in Italy which
deeply influenced my life and ambitions and gave name to the library
itself. It is and always will be my tiny tribute to a country that will
always be in my heart and to people that never ceased to inspire me.
*Grazie mille!*
