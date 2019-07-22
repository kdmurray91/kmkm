#!/usr/bin/env python3
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

inst_deps = [
    'zarr',
    'zict',
    'lmdb',
    'numpy',
]


setup(
    name="kmkm",
    packages=['kmkm', ],
    ext_modules=cythonize(Extension(
        "kmkm._kmkm",
        sources=["kmkm/_kmkm.pyx", ],
        include_dirs=["src", "src/ext", np.get_include()],
        extra_compile_args=['-std=c++14', ],
        libraries=['boost_serialization', 'boost_system', 'boost_filesystem',
                   'boost_iostreams', 'z'],
        language="c++",)),
    install_requires=inst_deps,
    entry_points="""
        [console_scripts]
        kmkm=kmkm.main:main
    """,
    setup_requires=[
        'numpy',
        'cython',
    ],
)
