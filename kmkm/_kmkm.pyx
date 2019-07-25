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
from cython.operator import dereference as deref
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
        void save(const string &filename) except +
        void load(const string &filename) except +
        vector[T] & counts() except +
        const T* data() except +
        int k() except +
        size_t nnz() except +

cdef extern from "kmseq.hh" namespace "kmseq":
    cdef cppclass KSeq:
        KSeq()
        string name, seq, qual

    cdef cppclass KSeqReader:
        KSeqReader()
        KSeqReader(const string &filename) except +
        void open(const string &filename) except +
        bool next_read(KSeq &ks) except +
        size_t next_chunk(vector[KSeq] &sequences, size_t max) except +

ctypedef KmerCounter[uint8_t] KmerCounterU8

@cython.auto_pickle(True)
cdef class PySeq:
    cdef KSeq *kseq

    def __init__(self, **kwargs):
        self.kseq = new KSeq()
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def __repr__(self):
        if not self.name and not self.seq:
            return ""
        if self.qual:
            return "@{}\n{}\n+\n{}\n".format(
                self.name.decode('utf-8'),
                self.seq.decode('utf-8'),
                self.qual.decode('utf-8'),
            )
        return ">{}\n{}\n".format(
            self.name.decode('utf-8'),
            self.seq.decode('utf-8'),
        )

    property seq:
        def __get__(self):
            return self.kseq.seq
        def __set__(self, seq):
            if isinstance(seq, str):
                seq = seq.encode("ascii")
            if isinstance(seq, bytes):
                self.kseq.seq = seq
            else:
                raise TypeError("seq must be str or bytes")

    property qual:
        def __get__(self):
            return self.kseq.qual
        def __set__(self, qual):
            if isinstance(qual, str):
                qual = qual.encode("ascii")
            if isinstance(qual, bytes):
                self.kseq.qual = qual
            else:
                raise TypeError("qual must be str or bytes")

    property name:
        def __get__(self):
            return self.kseq.name
        def __set__(self, name):
            if isinstance(name, str):
                name = name.encode("ascii")
            if isinstance(name, bytes):
                self.kseq.name = name
            else:
                raise TypeError("name must be str or bytes")

    def __dealloc__(self):
        if self.kseq is not NULL:
            del self.kseq
            self.kseq = NULL


@cython.auto_pickle(True)
cdef class PySeqReader:
    cdef KSeqReader *rdr

    def __init__(self, str filename):
        self.rdr = new KSeqReader(filename.encode('utf-8'))

    def __iter__(self):
        return self

    def __next__(self):
        seq = PySeq()
        ok = self.rdr.next_read(deref(seq.kseq))
        if not ok:
            raise StopIteration
        return seq

    def next_chunk(self, int max = 1000000):
        seqs = []
        for i in range(max):
            seq = PySeq()
            if not self.rdr.next_read(deref(seq.kseq)):
                break
            seqs.append(seq)
        return seqs

    def __dealloc__(self):
        if self.rdr is not NULL:
            del self.rdr
            self.rdr = NULL


@cython.auto_pickle(True)
cdef class PyKmerCounter:
    cdef readonly int ksize
    cdef readonly size_t cvsize
    cdef KmerCounterU8 *ctr

    def __init__(self, int ksize = 21, int cvsize = 1000000, bool canonical=True,
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


    def count_sequences(self, list sequences):
        for seq in sequences:
            assert isinstance(seq, PySeq)
            self.ctr.consume(seq.seq)

    def count_file(self, str filename):
        self.ctr.consume_from(filename.encode("utf-8"))

    def counts(self):
        n = np.asarray(<np.uint8_t[:self.cvsize]>self.ctr.data())
        return n.reshape((1, self.cvsize))

    def save(self, str filename):
        self.ctr.save(filename.encode("utf-8"))

    def __dealloc__(self):
        if self.ctr is not NULL:
            del self.ctr
            self.ctr = NULL

    property nnz:
        def __get__(self):
            return self.ctr.nnz()
