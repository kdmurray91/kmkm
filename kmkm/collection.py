# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import dask
import zarr
from zict import LMDB
import numpy as np

from os.path import basename

from ._kmkm import PyKmerCounter as KmerCounter
from .logger import LOGGER as LOG, enable_logging


class KmerCollection(object):
    def __init__(self, outprefix, mode='a', ksize=0, cvsize=0, cbf_tables=0):
        self.outprefix = outprefix
        self.store = LMDB(outprefix + '.zarr')
        self.array = zarr.open_array(store=self.store, mode=mode, shape=(0, 0),
                                     chunks=(1, 2**20), dtype='u1',
                                     compressor=zarr.Blosc('zstd', clevel=3))
        self.array.attrs['samples'] = self.array.attrs.get("samples", [])
        LOG.debug("Opened array at %s. Size %r, samples %r", outprefix,
                  self.array.shape, self.samples)
        self.ksize = ksize
        self.cvsize = cvsize
        self.cbf_tables = cbf_tables

    def _fname_to_sample(self, fname):
        bn = basename(fname)
        for ext in [".gz", ".kmr", ".fastq", ".fasta"]:
            if bn.endswith(ext):
                bn = bn[:-len(ext)]
        return bn

    @property
    def samples(self):
        return self.array.attrs["samples"]

    def set_countparams(self, ksize=31, cvsize=2**20, cbf_tables=3):
        self.ksize = ksize
        self.cvsize = cvsize
        self.cbf_tables = cbf_tables

    def count_readfile(self, filename):
        LOG.info("Counting %s", filename)
        kmr = KmerCounter(self.ksize, self.cvsize, cbf_tables=self.cbf_tables)
        kmr.count_file(filename)
        self._add_counter(kmr, filename)

    def add_file(self, filename):
        LOG.info("Adding %s", filename)
        kmr = KmerCounter(filename=filename)
        self._add_counter(kmr, filename)

    def _add_counter(self, kmr, filename):
        nrow, ncol = self.array.shape
        if ncol == 0:
            ncol = kmr.cvsize
            self.array.resize(0, ncol)
        elif ncol != kmr.cvsize:
            raise ValueError("Kmer counter size mismatch at " + filename)
        self.array.attrs['samples'].append(self._fname_to_sample(filename))
        LOG.debug("shape before adding %s: %r", filename, self.array.shape)
        LOG.debug("shape of %s:  %r", filename,  kmr.counts().shape)
        self.array.append(kmr.counts(), axis=0)
        LOG.debug("shape after %r", self.array.shape)

    def add_files(self, filenames):
        for f in filenames:
            self.add_file(f)
