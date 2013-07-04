#!/usr/bin/env python

from os.path import join


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration("poliastro", parent_package, top_path)

    config.add_library('ast2body',
                       sources=[join('poliastro', 'fortran', '*.for')])
    config.add_library('astiod',
                       sources=[join('poliastro', 'fortran', '*.for')])

    config.add_extension('_ast2body',
                         sources=['poliastro/ast2body.pyf'],
                         libraries=['ast2body'])
    config.add_extension('_astiod',
                         sources=['poliastro/astiod.pyf'],
                         libraries=['astiod'])

    config.add_data_dir(('tests', 'poliastro/tests'))

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    data_files = ['poliastro/octave/{}'.format(fn) for fn in
                  ('angl.m', 'constmath.m', 'mag.m',
                   'newtonm.m', 'newtonnu.m',
                   'rv2coe.m', 'uplanet_2013.m')]
    setup(version='0.1.0',
          description="poliastro - Utilities and Python wrappers for "
                      "Orbital Mechanics",
          author="Juan Luis Cano",
          author_email="juanlu001@gmail.com",
          url="https://github.com/Pybonacci/poliastro",
          license="BSD",
          keywords=[
              "aero", "aerospace", "engineering",
              "astrodynamics", "orbits", "kepler", "orbital mechanics"
          ],
          requires=["numpy", "scipy"],
          data_files=[('poliastro/octave', data_files)],
          packages=['poliastro'],
          classifiers=[
              "Development Status :: 2 - Pre-Alpha",
              "Intended Audience :: Education",
              "Intended Audience :: Science/Research",
              "License :: OSI Approved :: BSD License",
              "Operating System :: OS Independent",
              "Programming Language :: Python",
              "Programming Language :: Python :: 3",
              "Programming Language :: Python :: Implementation :: CPython",
              "Topic :: Scientific/Engineering",
              "Topic :: Scientific/Engineering :: Physics"
          ],
          long_description=open('README.rst').read(),
          configuration=configuration)
