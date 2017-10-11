// Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <numeric>
#include <omp.h>

#include "kmkm.hh"
#include "kmercollection.hh"
#include "minheap.hh"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#define ARMA_USE_HDF5
#include <armadillo>

using namespace std;
using namespace kmkm;
using namespace kmercollection;
using namespace arma;
namespace bfs = boost::filesystem;

struct BlupOpts
{
    string outfile;
    vector<string> countfiles;
    size_t top_n;
    vector<string> samplenames;
};

int blup_main(BlupOpts &opt)
{
    fmat freq;

    {
        KmerCollection samps(opt.top_n);
        if (samps.load(opt.outfile, opt.countfiles)) {
            cerr << "Using saved counts... ";
            for (const auto &sn: samps.names()) {
                cout << "  - " << sn << endl;
            }
        } else {
            samps.add_samples(opt.countfiles);
            samps.save(opt.outfile);
        }

        cerr << "Normalise, scale and centre counts... ";
         // L1 norm of counts, dim0 (cols, i.e. samples) have unit length
        freq = normalise(samps.counts(), 1, 0);
    }
    auto ukmer = mean(freq, 1); // Row mean
    auto sdkmer = stddev(freq, 0, 1); // row stddev with N-1 norm
    freq.each_col() -= ukmer;
    freq.each_col() /= sdkmer;
    // Remove NANs
    freq.for_each([](fmat::elem_type& v) { if (!is_finite(v)) v = 0.; });
    cerr << "and save them... ";
    freq.save(opt.outfile + ".scaledcounts", hdf5_binary_trans);
    cerr << "done" << endl;

    cerr << "Calc covar...";
    fmat covar = cov(freq);
    cerr << "and save it...";
    covar.save(opt.outfile + ".covar", csv_ascii);
    cerr << " done" << endl;


    cerr << "Calc cor...";
    fmat correl = cor(freq);
    cerr << "and save it...";
    correl.save(opt.outfile + ".cor", csv_ascii);
    cerr << " done" << endl;

    return 0;
}

int parse_args(BlupOpts &opt, int argc, char **argv)
{
    namespace po = boost::program_options;
    int num_threads;
    po::options_description flags("Options");
    flags.add_options()
        ("help,h",
            "Print this help")
        ("threads,t", po::value<int>(&num_threads)->default_value(16),
            "Number of parallel threads")
        ("top-n,n", po::value<size_t>(&opt.top_n)->default_value(1000000),
            "Number of count vector entries to consider")
        ("outfile,o", po::value<string>(&opt.outfile)->required(),
            "Output filename")
        ("countfiles", po::value<vector<string>>(&opt.countfiles)->required(),
            "Input count vector files (.kmr/.kmr.gz)");
    po::positional_options_description pos;
    pos.add("countfiles", -1);

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
    omp_set_num_threads(num_threads);
    return 0;
}

int main(int argc, char *argv[])
{
    BlupOpts opt;

    int ret = parse_args(opt, argc, argv);
    if (ret != 0) return ret;

    cout << "top_n: " << opt.top_n << endl;
    cout << "samples: " << endl;
    return blup_main(opt);
}

