# Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CXXFLAGS += -std=c++14 -O3 -Wall -g $(shell pkg-config --cflags zlib)
CPPFLAGS += -I. -isystem ext -fopenmp # -fsanitize=address
LIBS += -lboost_filesystem -lboost_system  -lboost_serialization -lboost_iostreams $(shell pkg-config --libs zlib)

prefix ?= /usr/local
PREFIX ?= $(prefix)

lib_srcs := 
lib_headers := $(wildcard *.hh)
test_srcs := test/main.cc $(wildcard test/test_*.cc)
test_prog := bin/kmkm_tests

.PHONY: all
all: $(test_prog) $(PROGS)

$(test_prog): $(test_srcs) $(lib_srcs) $(lib_headers)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -o $@ $<

.PHONY: test
test: $(test_prog)
	./$(test_prog) -s -r compact

.PHONY: install
install:
	mkdir -p $(PREFIX)/include
	cp $(lib_headers) $(PREFIX)/include/

.PHONY: clean
clean:
	rm -f $(test_prog) $(PROGS)
