# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import logging as L

__all__ = [
    'enable_logging',
    'LOGGER',
]

LOGGER = L.getLogger("kmkm")
LOGGER.setLevel(L.DEBUG)
LOGGER.addHandler(L.NullHandler())


def enable_logging(level=L.INFO):
    s = L.StreamHandler()
    s.setFormatter(L.Formatter('%(asctime)s: %(message)s'))
    s.setLevel(level)
    LOGGER.addHandler(s)
