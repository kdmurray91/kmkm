// Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <numeric>

#include "kmkm.hh"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace kmkm;
using namespace kmseq;

struct CountOpts
{
    int ksize;
    int cvsize;
    size_t cbf_tables;
    string outfile;
    vector<string> readfiles;
};

int count_main(CountOpts &opt)
{
    KmerCounter<> ctr(opt.ksize, UINT64_C(1)<<opt.cvsize);
    size_t nread = 0;
    for (auto readfile: opt.readfiles) {
        cerr << "  - " << readfile << endl;
        nread += ctr.consume_from(readfile);
    }
    cerr << "n_reads: " << nread << endl;
    cerr << "distinct_kmers: " << ctr.nnz() << endl;
    cerr << "total_kmers: " << accumulate(ctr.counts().begin(), ctr.counts().end(), 0) << endl;
    ctr.save(opt.outfile);
    return 0;
}

int parse_args(CountOpts &opt, int argc, char **argv)
{
    namespace po = boost::program_options;
    po::options_description flags("Options");
    flags.add_options()
        ("help,h",
            "Print this help")
        ("ksize,k", po::value<int>(&opt.ksize)->default_value(21),
            "Kmer size")
        ("cvsize,z", po::value<int>(&opt.cvsize)->default_value(25),
            "Counter size (size is 2^N bytes)")
        ("tables,t", po::value<size_t>(&opt.cbf_tables)->default_value(0),
            "Number of counting bloom filter tables. (0 disables the CBF, uses simple count vector)")
        ("outfile,o", po::value<string>(&opt.outfile)->required(),
            "Output filename")
        ("readfiles", po::value<vector<string>>(&opt.readfiles)->required(),
            "Input FASTX files");
    po::positional_options_description pos;
    pos.add("readfiles", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
                .options(flags).positional(pos).run(), vm);

    if (vm.count("help")) {
        cout << flags << endl;
        return 0;
    }
    try {
        po::notify(vm);
    } catch (po::error &e) {
        cerr << e.what() << endl << endl;
        cerr << flags << endl;
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    CountOpts opt;

    int ret = parse_args(opt, argc, argv);
    if (ret != 0) return ret;

    cerr << "ksize: " << opt.ksize << endl;
    cerr << "cvsize: " << opt.cvsize << endl;
    cerr << "readfiles: " << endl;
    return count_main(opt);
}
