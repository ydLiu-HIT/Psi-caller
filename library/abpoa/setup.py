#!/usr/bin/env python
# coding=utf-8
import os
from distutils.core import setup, Extension
from Cython.Build import cythonize


simd_flag = "-msse2"

SSE41 = "__SSE4_1__"
AVX2 = "__AVX2__"
AVX512F = "__AVX512F__"
AVX512BW = "__AVX512BW__"

FLAG_SSE2 = "-msse2"
FLAG_SSE41 = "-msse4.1"
FLAG_AVX2 = "-mavx2"
FLAG_AVX512F = "-mavx512f"
FLAG_AVX512BW = "-mavx512bw"

cmd = "./simd_check 2>/dev/null"
sflag = os.popen(cmd).read().strip()
print(sflag)

if sflag == AVX512BW:
    simd_flag = FLAG_AVX512BW
elif sflag == AVX512F:
    simd_flag = FLAG_AVX512F
elif sflag == AVX2:
    simd_flag = FLAG_AVX2
elif sflag == SSE41:
    simd_flag = FLAG_SSE41
print(simd_flag)

extra_compile_args = ['-O3', '-Wno-error=declaration-after-statement', simd_flag]

# define the extension module
extensions = cythonize([
    Extension("poa_module",
        sources = ['pyPOA.c', 'abpoa_align.c', 'abpoa_graph.c', 'simd_abpoa_align.c', 'utils.c', 'simd_check.c', 'abpoa_dot_plot.c', 'abpoa_plot.c', 'agglo_hier_clu.c'],
        depends = ['abpoa.h', 'abpoa_align.h', 'abpoa_graph.h', 'simd_abpoa_align.h', 'utils.h', 'simd_instruction.h', 'ksort.h', 'kseq.h', 'seq.h', 'kdq.h', 'align.h', 'agglo_hier_clu.h'],
        libraries = ['z', 'm', 'pthread'],
        extra_compile_args=extra_compile_args)
])

setup(name="POA",
      version='0.1',
      license='MIT',
      description='Python implementation of abPOA',
      long_description="Multiple sequence alignment. Please visit the github pages for more details.",
      author='ydliu, yangao007',
      author_email="ydliu@hit.edu.cn",
      ext_modules=extensions
    )
