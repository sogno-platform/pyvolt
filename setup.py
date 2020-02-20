#!/usr/bin/env python

from setuptools import setup

setup(name='pyvolt',
      version='0.2',
      description='State Estimation',
      author='Marco Pau, Markus Mirz, Jan Dinkelbach',
      author_email='acs-software@eonerc.rwth-aachen.de',
      packages = [ 'pyvolt' ],
      install_requires = [
        'matplotlib',
        'numpy',
        'cimpy'        
      ]
)