#!/usr/bin/env python3
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

inst_deps = [
    'zarr',
    'numpy',
]


setup(
    name="kmkm",
    ext_modules=cythonize(Extension(
        "kmkm._kmkm",
        sources=["kmkm/_kmkm.pyx", ],
        include_dirs=["src", "src/ext", np.get_include()],
        extra_compile_args=['-std=c++14', ],
        libraries=['boost_serialization', 'boost_system', 'boost_filesystem',
                   'boost_iostreams', 'z'],
        language="c++",)),
    install_requires=inst_deps,
    setup_requires=[
        'numpy',
    ],
)
