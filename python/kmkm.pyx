# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdint cimport uint8_t
cimport numpy as np
import numpy as np

cdef extern from "kmkm.hh" namespace "kmkm":
    cdef cppclass KmerCounter[T]:
        KmerCounter (int ksize, size_t cvsize, bool canonical, size_t cbf_tables) except +
        void consume_from(const string &filename) except +
        void consume(const string &sequence) except +
        vector[T] counts() except +

ctypedef KmerCounter[uint8_t] KmerCounterU8



cdef class KmerCounterSet:
    cdef int ksize
    cdef size_t cvsize
    cdef bool canonical
    cdef size_t cbf_tables
    cdef public list samples
    cdef public np.ndarray global_array

    def __init__(self, int ksize, int cvsize, bool canonical=True, int cbf_tables=0):
        self.ksize = ksize
        self.cvsize = cvsize
        self.canonical = canonical
        self.cbf_tables = cbf_tables
        self.global_array = np.zeros((0, cvsize), dtype=np.uint8)
        self.samples = []


    def count_file(self, str filename):
        cdef KmerCounterU8 *ctr = new KmerCounterU8(self.ksize, self.cvsize, self.canonical, self.cbf_tables)
        ctr.consume_from(filename.encode("ascii"))
        c = ctr.counts()
        n = np.asarray(<np.uint8_t[:c.size()]>c.data()).reshape((1, self.cvsize))
        self.global_array = np.concatenate((self.global_array, n), axis=0)
        self.samples.append(filename)


