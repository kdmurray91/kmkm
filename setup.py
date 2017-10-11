#!/usr/bin/env python3
from setuptools import setup, Extension
from Cython.Build import cythonize

inst_deps = [
    'zarr',
    'numpy',
]


setup(
    name="kmkm",
    ext_modules=cythonize(Extension(
        "kmkm._kmkm",
        sources=["kmkm/_kmkm.pyx", ],
        include_dirs=["src", "src/ext"],
        extra_compile_args=['-std=c++14', ],
        libraries=['boost_serialization', 'boost_system', 'boost_filesystem',
                   'boost_iostreams', 'z'],
        language="c++",)),
    install_requires=inst_deps,
)
