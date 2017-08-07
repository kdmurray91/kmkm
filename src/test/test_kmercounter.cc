// Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


TEST_CASE("KmerCounter init", "[KmerCounter]") {
    const size_t cvsize = 10000;
    const int k = 4;
    KmerCounter<uint8_t> ctr(k, cvsize);

    SECTION("Counter init") {
        const auto &cv = ctr.counts();
        REQUIRE(cv.size() == cvsize);
        bool all_zero = true;
        for (size_t i = 0; i < cvsize; i++) {
            all_zero &= (cv[i] == 0);
        }
        REQUIRE(all_zero);

        KmerIterator i("AAAA", 20);
        REQUIRE(i.finished());
    }

    SECTION("Collision rate") {
        ctr.consume("AAAA");
        REQUIRE(ctr.nnz() == 1);
        REQUIRE(ctr.collision_rate() == 1.0/double(cvsize));

        // Check there is no collision with self
        ctr.consume("AAAA");
        REQUIRE(ctr.nnz() == 1);
        REQUIRE(ctr.collision_rate() == 1.0/double(cvsize));
    }
}



// vim:set et sw=4 ts=4:
