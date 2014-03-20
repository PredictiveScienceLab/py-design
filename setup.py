#!/usr/bin/env python


from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import os
import glob


setup(name='Design points for random experiments',
      author='Ilias Bilionis',
      version='0.0',
      ext_modules=[Extension('design._design',
                            glob.glob(os.path.join('src', '*.f90')),
                             libraries=['python2.7', 'pthread'],
                             extra_link_args=['-shared'])],
      packages=['design'])
