from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from sys import platform as _platform
import os
import numpy as np

#openmp_arg = '-fopenmp'
#if _platform == "win32":
#    openmp_arg = '-openmp'

extensions = [
  Extension(
    'multicut', ['multicut.pyx', 'src/nl-lmp.cxx'],
    language="c++",
    include_dirs=[np.get_include(), '.', 'include', 'src'],
    extra_compile_args=['-std=c++11','-O3', '-DHAVE_CPP11_INITIALIZER_LISTS'],
    extra_link_args=['-std=c++11', '-L./']
  )
] 

setup(
    name = 'multicut',
    ext_modules = cythonize(extensions)
)
