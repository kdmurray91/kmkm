# distutils: language = c++
#
# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdint cimport uint8_t
cimport numpy as np
import numpy as np
import sys

cdef extern from "kmkm.hh" namespace "kmkm":
    cdef cppclass KmerCounter[T]:
        KmerCounter()
        KmerCounter(const string &filename)
        KmerCounter(int ksize, size_t cvsize, bool canonical, size_t cbf_tables) except +
        void consume_from(const string &filename) except +
        void consume(const string &sequence) except +
        vector[T] & counts() except +
        const T* data() except +
        int k() except +

ctypedef KmerCounter[uint8_t] KmerCounterU8

cdef class PyKmerCounter:
    cdef readonly int ksize
    cdef readonly size_t cvsize
    cdef KmerCounterU8 *ctr

    def __init__(self, int ksize = 0, int cvsize = 0, bool canonical=True,
                 int cbf_tables=0, str filename=None):
        if filename is None:
            self.ksize = ksize
            self.cvsize = cvsize
            self.ctr = new KmerCounterU8(self.ksize, self.cvsize, canonical,
                                        cbf_tables)
        else:
            self.ctr = new KmerCounterU8(filename.encode('utf-8'))
            self.ksize = self.ctr.k()
            self.cvsize = self.ctr.counts().size()


    def count_file(self, str filename):
        self.ctr.consume_from(filename.encode("utf-8"))

    def counts(self):
        n = np.asarray(<np.uint8_t[:self.cvsize]>self.ctr.data())
        return n.reshape((1, self.cvsize))

    def __dealloc__(self):
        if self.ctr is not NULL:
            del self.ctr
            self.ctr = NULL
