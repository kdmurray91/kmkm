// Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



TEST_CASE("KmerIterator basics", "[KmerIterator]") {
    SECTION("Too small") {
        KmerIterator i("AAAA", 20);
        REQUIRE(i.size() == 0);
        REQUIRE(i.finished());
    }

    SECTION("Length correct") {
        const string seq = "AATTAATTAATT";
        const int k = 4;
        KmerIterator i(seq, k);
        size_t len = 0;
        REQUIRE(i.size() == seq.size() - k + 1);
        REQUIRE(!i.finished());
        while (!i.finished()) {
            len++;
            i.next();
        }
        REQUIRE(len == seq.size() - k + 1);
    }
}


void _value_test(KmerIterator &ki, vector<uint64_t> &expected)
{
    size_t i = 0;
    while(!ki.finished()) {
        auto got = ki.next();
        INFO("Testing at " << i << " " << got);
        REQUIRE(got == expected[i++]);
    }
    REQUIRE(i == ki.size());
}

TEST_CASE("kmer values", "[KmerIterator]") {
    SECTION("1mer") {
        KmerIterator ki("ACGTACGT", 1, false);
        vector<uint64_t> expected {0, 1, 2, 3, 0, 1, 2, 3};
        size_t i = 0;
        while(!ki.finished()) {
            auto got = ki.next();
            INFO("Testing at " << i << " " << got);
            REQUIRE(got == expected[i++]);
        }
        REQUIRE(i == ki.size());
    }

    SECTION("1mer canoncial") {
        KmerIterator ki("ACGTACGT", 1, true);
        vector<uint64_t> expected {0, 1, 1, 0, 0, 1, 1, 0};
        size_t i = 0;
        while(!ki.finished()) {
            auto got = ki.next();
            INFO("Testing at " << i << " " << got);
            REQUIRE(got == expected[i++]);
        }
        REQUIRE(i == ki.size());
    }

}


// vim:set et sw=4 ts=4:
