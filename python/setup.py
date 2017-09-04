#!/usr/bin/env python3
from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(
    name="kmkm",
    ext_modules=cythonize(Extension(
        "kmkm._kmkm",
        sources=["kmkm/_kmkm.pyx", ],
        include_dirs=["../src", "../src/ext"],
        extra_compile_args=['-std=c++14', ],
        libraries=['boost_serialization', 'boost_iostreams', 'z'],
        language="c++",)),
)
