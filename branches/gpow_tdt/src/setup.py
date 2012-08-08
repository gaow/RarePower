#!/usr/bin/env python

from distutils.core import setup
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

import sys, os
setup(name = 'gpow_ext',
    version = '0.0.1',
    py_modules = [
        'pyProgressBar.__init__',
        'pyProgressBar.terminal',
        'pyProgressBar.progressbar'
    ],
    scripts = ['gpow-ext'],
    cmdclass = {'build_py': build_py },
    package_dir = {'pyProgressBar': 'pyProgressBar'}
    )
