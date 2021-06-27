#!/usr/bin/env python
# coding=utf-8

#from distutils.core import setup, Extension
from Cython.Build import cythonize
try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension


extra_compile_args = []

extra_compile_args = extra_compile_args + ['-O3', '-Wno-error=declaration-after-statement']

extensions = Extension("mantaAss_module", 
        sources = ["mantaAssembler.cpp"], 
        depends = ['mantaAssembler.hpp'],
        extra_compile_args=extra_compile_args)

setup(name="mantaAss_module",
      license='MIT',
      description='Python implementation of local assembly based on de Bruijn graphg',
      author='ydliu',
      ext_modules=[extensions]
    )
