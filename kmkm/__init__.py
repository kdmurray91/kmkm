# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from ._kmkm import PyKmerCounter as KmerCounter
from .logger import LOGGER as LOG, enable_logging
from .collection import KmerCollection

__all__ = [
    "KmerCounter",
    "KmerCollection",
    "enable_logging",
]
