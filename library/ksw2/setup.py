#!/usr/bin/env python
# coding=utf-8

#from distutils.core import setup, Extension
from Cython.Build import cythonize
try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension

import sys, platform

extra_compile_args = []

if platform.machine() in ["aarch64", "arm64"]:
	extra_compile_args.extend(['-ftree-vectorize', '-DKSW_SSE2_ONLY', '-D__SSE2__'])
else:
	extra_compile_args.append('-msse4.1') # WARNING: ancient x86_64 CPUs don't have SSE4

extra_compile_args = extra_compile_args + ['-O3', '-Wno-error=declaration-after-statement']

extensions = Extension("ksw2_module", 
        sources = ["ksw2.c", "ksw2_extz2_sse.c", "ksw2_extd2_sse.c", "ksw2_dispatch.c", 'ksw2_ll_sse.c'], 
        depends = ['ksw2.h'],
        libraries = ['z', 'm'],
        extra_compile_args=extra_compile_args)

setup(name="ksw2_module",
      license='MIT',
      url = 'https://www.github.com/lh3/ksw2',
      description='Python implementation of ksw2 written by LiHeng',
      long_description="Python implementation of ksw2 (Used for local and global alignment). Please visit the github pages for more details.",
      author='Heng Li',
      ext_modules=[extensions]
    )
