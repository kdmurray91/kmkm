# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import zarr
from zict import LMDB
import numpy as np
from numcodecs import Blosc

from os.path import basename, exists, isdir
import shutil

from ._kmkm import PyKmerCounter as KmerCounter
from .logger import LOGGER as LOG, enable_logging


class KmerCollection(object):
    def __init__(self, outfile, mode='a', ksize=0, cvsize=0, cbf_tables=0):
        self.outfile = outfile
        if mode.startswith("w"):
            if exists(outfile):
                LOG.debug("Removing existing array at %s", outfile)
                shutil.rmtree(outfile)
        if not exists(outfile):
            LOG.debug("Opening new array at %s", outfile)
            self.array = zarr.open(
                outfile, mode=mode, shape=(0, 0), chunks=(1, 2**24),
                dtype='u1', compressor=Blosc('zstd', clevel=3)
            )
        else:
            LOG.debug("Opening existing array at %s", outfile)
            self.array = zarr.open(outfile, mode=mode)
        self.array.attrs['samples'] = self.array.attrs.get("samples", [])

        LOG.debug("Opened array at %s. Size %r, samples %r", outfile,
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
    @samples.setter
    def samples(self, samps):
        self.array.attrs["samples"] = samps


    def set_countparams(self, ksize=21, cvsize=2**20, cbf_tables=1):
        self.ksize = ksize
        self.cvsize = cvsize
        self.cbf_tables = cbf_tables

    def count_seqfile(self, filename):
        LOG.info("Counting %s", filename)
        kmr = KmerCounter(self.ksize, self.cvsize, cbf_tables=self.cbf_tables)
        kmr.count_file(filename)
        self._add_counter(kmr, filename)

    def add_file(self, filename):
        LOG.info("Adding %s", filename)
        if os.path.isdir(filename):
            raise ValueError("%s must be a file", filename)
        kmr = KmerCounter(filename=filename)
        self._add_counter(kmr, filename)

    def _add_counter(self, kmr, filename):
        nrow, ncol = self.array.shape
        if ncol == 0:
            ncol = kmr.cvsize
            self.array.resize(0, ncol)
        elif ncol != kmr.cvsize:
            raise ValueError("Kmer counter size mismatch at " + filename)
        sname = self._fname_to_sample(filename)
        if sname in self.samples:
            raise ValueError("Duplcate sample name! (" + sname + ")")
        self.samples = self.samples + [sname]
        LOG.debug("shape before adding %s: %r", filename, self.array.shape)
        LOG.debug("shape of %s:  %r", filename,  kmr.counts().shape)
        self.array.append(kmr.counts(), axis=0)
        LOG.debug("shape after %r", self.array.shape)
        LOG.debug("samples after %r", self.samples)

    def add_files(self, filenames):
        for f in filenames:
            self.add_file(f)
