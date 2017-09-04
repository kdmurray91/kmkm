# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CXXFLAGS += -std=c++14 -O3 -Wall -g
CPPFLAGS += -I src -isystem src/ext -fopenmp
LIBS += -lboost_program_options -lboost_serialization -lboost_iostreams -lz -larmadillo -lhdf5

prefix ?= /usr/local
PREFIX ?= $(prefix)

lib_srcs := 
lib_headers := src/kmkm.hh
test_srcs := src/test/main.cc
test_prog := bin/kmkm_tests
count_prog := bin/kmkm_count
blup_prog := bin/kblup
PROGS = $(blup_prog) $(count_prog)

.PHONY: all
all: $(test_prog) $(PROGS)

bin/%: ./src/%.cc $(lib_srcs) $(lib_headers)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -o $@ $<

$(test_prog): ./src/test/main.cc $(lib_srcs) $(lib_headers)
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
	rm -f $(test_prog) $(PROGS)
