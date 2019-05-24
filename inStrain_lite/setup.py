#!/usr/bin/env python

import os
from setuptools import setup

from inStrain_lite._version import __version__

setup(name='inStrain_lite',
      version=__version__,
      description='Calculation of strain-level metrics',
      url='https://github.com/alexcritschristoph/soil_popgen',
      author='Matt Olm and Alex Crits-Christoph',
      author_email='crits-christoph@berkeley.edu',
      license='MIT',
      package_data={'inStrain_lite': ['helper_files/combined_null1000000.txt']},
      include_package_data=True,
      packages=['inStrain_lite'],
      scripts=['bin/inStrain_lite'],
      python_requires='>=3.4.0',
      install_requires=[
          'numpy',
          'pandas',
          'seaborn',
          'matplotlib',
          'biopython',
          'scikit-learn',
          'pytest',
          'tqdm',
          'pysam',
          'networkx'
      ],
      zip_safe=False)
