# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from ._kmkm import PyKmerCounter as KmerCounter
import zarr
import numpy as np


class KmerCounterSet(object):

    def __init__(self, ksize, cvsize, canonical=True, cbf_tables=0):
        self.ksize = ksize
        self.cvsize = cvsize
        self.canonical = canonical
        self.cbf_tables = cbf_tables
        self.global_array = zarr.zeros((0, cvsize))
        self.samples = []

    def count_file(self, filename):
        ctr = KmerCounter(self.ksize, self.cvsize, self.canonical,
                          self.cbf_tables)
        ctr.consume_from(filename)
        self.global_array = np.concatenate((self.global_array, ctr.counts()), axis=0)
        self.samples.append(filename)

    def load_file(self, filename):
        pass
