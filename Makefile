# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CXXFLAGS += -std=c++14 -O3 -Wall
CPPFLAGS += -I src

prefix ?= /usr/local
PREFIX ?= $(prefix)

lib_srcs := 
lib_headers := src/kmkm.hh
test_srcs := src/test/main.cc

tests: $(lib_srcs) $(test_srcs) | $(lib_headers) $(wildcard src/test/*.cc)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^

.PHONY: test
test: tests
	./tests -s -r compact


.PHONY: install
install:
	mkdir -p $(PREFIX)/include
	cp src/kmkm.hh $(PREFIX)/include/

.PHONY: clean
clean:
	rm -f tests
