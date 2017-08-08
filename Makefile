# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CXXFLAGS += -std=c++14 -O3 -Wall -g
CPPFLAGS += -I src -isystem src/ext
LIBS += -lboost_program_options -lboost_serialization -lboost_iostreams -lz

prefix ?= /usr/local
PREFIX ?= $(prefix)

lib_srcs := 
lib_headers := src/kmkm.hh
test_srcs := src/test/main.cc
test_prog := bin/kmkm_tests
count_prog := bin/kmkm_count

.PHONY: all
all: $(test_prog) $(count_prog)

$(test_prog): $(test_srcs) $(lib_headers) $(wildcard src/test/*.cc)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -o $@ $<

$(count_prog): ./src/kmer_counter.cc $(lib_srcs) $(lib_headers)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -o $@ $<

.PHONY: test
test: $(test_prog)
	./$(test_prog) -s -r compact

.PHONY: install
install:
	mkdir -p $(PREFIX)/include
	cp src/kmkm.hh $(PREFIX)/include/

.PHONY: clean
clean:
	rm -f $(test_prog) $(count_prog)
