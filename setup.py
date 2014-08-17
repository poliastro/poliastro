#!/usr/bin/env python

from os.path import join

from numpy.distutils.core import setup

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration("poliastro", parent_package, top_path)

    config.add_library('ast2body',
                       sources=[join('poliastro', 'src', '*.for')])
    config.add_library('astiod',
                       sources=[join('poliastro', 'src', '*.for')])

    config.add_extension('twobody._ast2body',
                         sources=[join('poliastro', 'src', 'ast2body.pyf')],
                         libraries=['ast2body'])
    config.add_extension('twobody._astiod',
                         sources=[join('poliastro', 'src', 'astiod.pyf')],
                         libraries=['astiod'])

    config.add_data_dir(('tests', 'poliastro/tests'))

    return config

if __name__ == '__main__':
    setup(version='0.3.0-dev',
          description="poliastro - Utilities and Python wrappers for "
                      "Orbital Mechanics",
          author="Juan Luis Cano",
          author_email="juanlu001@gmail.com",
          url="http://pybonacci.github.io/poliastro/",
          download_url="https://github.com/Pybonacci/poliastro",
          license="BSD",
          keywords=[
              "aero", "aerospace", "engineering",
              "astrodynamics", "orbits", "kepler", "orbital mechanics"
          ],
          requires=["numpy", "astropy", "pytest"],
          packages=['poliastro', 'poliastro.twobody'],
          classifiers=[
              "Development Status :: 3 - Pre-Alpha",
              "Intended Audience :: Education",
              "Intended Audience :: Science/Research",
              "License :: OSI Approved :: BSD License",
              "Operating System :: OS Independent",
              "Programming Language :: Python",
              "Programming Language :: Python :: 3",
              "Programming Language :: Python :: Implementation :: CPython",
              "Topic :: Scientific/Engineering",
              "Topic :: Scientific/Engineering :: Physics",
              "Topic :: Scientific/Engineering :: Astronomy",
          ],
          long_description=open('README').read(),
          configuration=configuration)
