# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import zarr
from zict import LMDB
import numpy as np
from numcodecs import Blosc


import multiprocessing
from concurrent.futures import as_completed, ThreadPoolExecutor

from os.path import basename, exists, isdir
import shutil

from ._kmkm import PyKmerCounter as KmerCounter
from .logger import LOGGER as LOG, enable_logging



class KmerCollection(object):
    def __init__(self, outfile, mode='a', ksize=0, cvsize=0, cbf_tables=0, chunks=(2**3, 2**17)):
        self.outfile = outfile
        # sync disabled due to new threading code
        self.sync = None #zarr.ProcessSynchronizer(outfile+'.sync')
        if mode.startswith("w"):
            if exists(outfile):
                LOG.debug("Removing existing array at %s", outfile)
                shutil.rmtree(outfile)
        if not exists(outfile):
            LOG.debug("Opening new array at %s", outfile)
            self.array = zarr.open(
                outfile, mode=mode, shape=(0, 0), chunks=chunks,
                dtype='u1', compressor=Blosc('zstd', clevel=8),
                synchronizer=self.sync
            )
        else:
            LOG.debug("Opening existing array at %s", outfile)
            self.array = zarr.open(outfile, mode=mode, synchronizer=self.sync)
        self.array.attrs['samples'] = self.array.attrs.get("samples", [])

        LOG.debug("Opened array at %s. Size %r, samples %r", outfile,
                  self.array.shape, self.samples)
        self.ksize = ksize
        self.cvsize = cvsize
        self.cbf_tables = cbf_tables

    def _fname_to_sample(self, fname):
        bn = basename(fname)
        for ext in [".gz", ".bz2", ".kmr", ".fastq", ".fasta"]:
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

    def _make_counter(self):
        return KmerCounter(self.ksize, self.cvsize, cbf_tables=self.cbf_tables)

    def count_seqfile(self, filename, stats=True):
        kmr = self._make_counter()
        kmr.count_file(filename)
        self.add_counter(kmr, filename)
        if stats:
            return filename, kmr.nnz

    def add_file(self, filename):
        LOG.info("Adding %s", filename)
        if os.path.isdir(filename):
            raise ValueError("%s must be a file", filename)
        kmr = KmerCounter(filename=filename)
        self.add_counter(kmr, filename)

    def add_counter(self, kmr, samplename):
        self.add_counts(kmr.counts(), samplename)

    def add_counts(self, countvec, samplename):
        nrow, ncol = self.array.shape
        if ncol == 0:
            ncol = countvec.shape[1]
            self.array.resize(0, ncol)
        elif ncol != countvec.shape[1]:
            raise ValueError("Kmer counter size mismatch at " + samplename)
        sname = self._fname_to_sample(samplename)
        if sname in self.samples:
            raise ValueError("Duplcate sample name! (" + sname + ")")
        self.samples = self.samples + [sname]
        self.array.append(countvec, axis=0)

    def count_seqfiles(self, filenames, njobs=1, callback=print):
        with ThreadPoolExecutor(max_workers=njobs) as executor:
            for res in executor.map(self.count_seqfile, filenames):
                try:
                    callback(*res)
                except Exception as exc:
                    LOG.error(str(exc))
