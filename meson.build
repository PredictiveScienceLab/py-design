project('design', 'c', 'fortran',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

fortran_src_files = [
  'src/faure.f90',
  'src/halton.f90',
  'src/hammersley.f90',
  'src/ihs.f90',
  'src/lambert.f90',
  'src/latin_center.f90',
  'src/latin_cover.f90',
  'src/latin_edge.f90',
  'src/latin_random.f90',
  'src/latinize.f90',
  'src/sandia_sparse.f90',
  'src/sobol.f90'
]

designfortran_source = custom_target('_designfortranmodule.c',
  input : fortran_src_files,
  output : ['_designfortranmodule.c', '_designfortran-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_designfortran', '--lower']
)

inc_np = include_directories(incdir_numpy, incdir_f2py)

py.extension_module('_designfortran',
  fortran_src_files,
  designfortran_source,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true
)

py.install_sources(
  'design/__init__.py',
  subdir: 'design'
)