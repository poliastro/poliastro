#!/usr/bin/env python
# coding: utf-8

# http://stackoverflow.com/a/10975371/554319
import io
from setuptools import setup, find_packages


# http://blog.ionelmc.ro/2014/05/25/python-packaging/
setup(
    name="poliastro",
    version='0.5.dev0',
    description="Python package for Orbital Mechanics",
    author="Juan Luis Cano",
    author_email="hello@juanlu.space",
    url="http://poliastro.github.io/",
    download_url="https://github.com/poliastro/poliastro",
    license="MIT",
    keywords=[
      "aero", "aerospace", "engineering",
      "astrodynamics", "orbits", "kepler", "orbital mechanics"
    ],
    install_requires=[
        "numpy",
        "numba>=0.23",
        "astropy>=1.0",
        "matplotlib",
        "jplephem",
        "scipy",
    ],
    tests_require=[
        "pytest"
    ],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points={
      'console_scripts': [
          'poliastro = poliastro.cli:main'
      ]
    },
    classifiers=[
      "Development Status :: 4 - Beta",
      "Intended Audience :: Education",
      "Intended Audience :: Science/Research",
      "License :: OSI Approved :: MIT License",
      "Operating System :: OS Independent",
      "Programming Language :: Python",
      "Programming Language :: Python :: 3",
      "Programming Language :: Python :: 3.4",
      "Programming Language :: Python :: 3.5",
      "Programming Language :: Python :: Implementation :: CPython",
      "Topic :: Scientific/Engineering",
      "Topic :: Scientific/Engineering :: Physics",
      "Topic :: Scientific/Engineering :: Astronomy",
    ],
    long_description=io.open('README', encoding='utf-8').read(),
    package_data={"poliastro": ['tests/*.py']}
)
