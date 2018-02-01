# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


def partition(items, nparts):
    """ Split ``items`` into ``nparts`` partitions.

    :param items: Items to be partitions
    :param nparts: Number of partitions
    :returns: list of lists containing each partition

    Examples
    --------

    >>> partition(range(10), 2)
    [[0, 2, 4, 6, 8], [1, 3, 5, 7, 9]]
    >>> partition(range(11), 2)
    [[0, 2, 4, 6, 8, 10], [1, 3, 5, 7, 9]]
    >>> partition(range(10), 3)
    [[0, 3, 6, 9], [1, 4, 7], [2, 5, 8]]
    >>> partition(range(3), 4)
    [[0], [1], [2], []]
    """
    parts = [list() for _ in range(nparts)]
    for i, item in enumerate(items):
        parts[i % nparts].append(item)
    return parts


