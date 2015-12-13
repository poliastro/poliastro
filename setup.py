#!/usr/bin/env python
# coding: utf-8

# http://stackoverflow.com/a/10975371/554319
import io
from setuptools import setup, find_packages


# http://blog.ionelmc.ro/2014/05/25/python-packaging/
setup(
    name="poliastro",
    version='0.4.0.dev0',
    description="poliastro - Python package for Orbital Mechanics",
    author="Juan Luis Cano",
    author_email="juanlu001@gmail.com",
    url="http://poliastro.github.io/",
    download_url="https://github.com/poliastro/poliastro",
    license="MIT",
    keywords=[
      "aero", "aerospace", "engineering",
      "astrodynamics", "orbits", "kepler", "orbital mechanics"
    ],
    requires=["numpy",
            "numba",
            "astropy",
            "matplotlib",
            "jplephem",
            "scipy",
            "pytest"],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points={
      'console_scripts': [
          'poliastro = poliastro.cli:main'
      ]
    },
    classifiers=[
      "Development Status :: 3 - Alpha",
      "Intended Audience :: Education",
      "Intended Audience :: Science/Research",
      "License :: OSI Approved :: MIT License",
      "Operating System :: OS Independent",
      "Programming Language :: Python",
      "Programming Language :: Python :: 2",
      "Programming Language :: Python :: 2.7",
      "Programming Language :: Python :: 3",
      "Programming Language :: Python :: 3.3",
      "Programming Language :: Python :: 3.4",
      "Programming Language :: Python :: Implementation :: CPython",
      "Topic :: Scientific/Engineering",
      "Topic :: Scientific/Engineering :: Physics",
      "Topic :: Scientific/Engineering :: Astronomy",
    ],
    long_description=io.open('README', encoding='utf-8').read(),
    package_data={"poliastro": ['tests/*.py']}
)
