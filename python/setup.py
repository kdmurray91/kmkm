#!/usr/bin/env python3
from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(
    name="kmkm",
    ext_modules=cythonize(Extension(
        "kmkm/_kmkm",
        sources=["kmkm/_kmkm.pyx", ],
        include_dirs=["../src", "../src/ext"],
        language="c++",)),
)
