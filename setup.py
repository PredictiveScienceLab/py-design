#!/usr/bin/env python


from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import os
import glob


setup(name='py-design',
      version='2.0',
      description='Design points for random experiments',
      author='Ilias Bilionis',
      author_email='ibilion@purdue.edu',
      url='https://github.com/ebilionis/py-design',
      download_url='https://github.com/ebilionis/py-design/tarball/1.0',
      keywards=['design', 'random experiment', 'latin hypercube',
                'sparse grid', 'computer experiments'],
      ext_modules=[Extension('design._design',
                            glob.glob(os.path.join('src', '*.f90')),
                            extra_link_args=['-undefined dynamic_lookup'])],
      packages=['design'])
