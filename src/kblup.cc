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

#include "kmkm.hh"


#include <boost/program_options.hpp>
#define ARMA_USE_HDF5
#include <armadillo>

using namespace std;
using namespace kmkm;

struct BlupOpts
{
    string outfile;
    vector<string> countfiles;
    size_t top_n;
    int n_threads;
};

int blup_main(BlupOpts &opt)
{
    using namespace arma;
    const size_t nsamp = opt.countfiles.size();
    map<size_t, string> names;

    fmat counts (opt.top_n, nsamp);

    cerr << "Loading individual countfiles into count matrix..." << endl;
    #pragma omp parallel for schedule(dynamic) num_threads(opt.n_threads) shared(names)
    for (size_t j = 0; j < nsamp; j++) {
        const auto & countfile = opt.countfiles[j];
        KmerCounter<> ctr(countfile);
        const auto &cv = ctr.counts();
        const size_t N = opt.top_n == 0 ? cv.size() : opt.top_n; // -n 0 = use all bins
        for (size_t i = 0; i < n; i++) {
            counts(i, j) = cv[i];
        }
        #pragma omp critical
        {
            cout << "  - " << countfile << endl;
            names[j] = countfile;
        }
    }

    cerr << "Save counts...";
    counts.save(opt.outfile + ".counts.h5", hdf5_binary_trans);
    cerr << " done" << endl;

    cerr << "Calc covar...";
    fmat freq = counts;//normalise(counts, 1, 0);
    freq.save(opt.outfile + ".freq.h5", hdf5_binary_trans);
    auto ukmer = mean(freq, 1); // Row mean
    auto sdkmer = stddev(freq, 0, 1); // row stddev with N-1 norm
    freq.each_col() -= ukmer;
    freq.each_col() /= sdkmer;
    freq.save(opt.outfile + ".normfreq.h5", hdf5_binary_trans);
    fmat covar = cov(freq);
    cerr << " done" << endl;

    cerr << "Save covar...";
    covar.save(opt.outfile + ".covar.asc", csv_ascii);
    cerr << " done" << endl;

    return 0;
}

int parse_args(BlupOpts &opt, int argc, char **argv)
{
    namespace po = boost::program_options;
    po::options_description flags("Options");
    flags.add_options()
        ("help,h",
            "Print this help")
        ("threads,t", po::value<int>(&opt.n_threads)->default_value(16),
            "Number of parallel threads")
        ("top-n,n", po::value<size_t>(&opt.top_n)->default_value(1000000),
            "Number of count vector entries to consider")
        ("outfile,o", po::value<string>(&opt.outfile)->required(),
            "Output filename")
        ("countfiles", po::value<vector<string>>(&opt.countfiles)->required(),
            "Input FASTX files");
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
    return 0;
}

int main(int argc, char *argv[])
{
    BlupOpts opt;

    int ret = parse_args(opt, argc, argv);
    if (ret != 0) return ret;

    cout << "top_n: " << opt.top_n << endl;
    cout << "countfiles: " << endl;
    return blup_main(opt);
}

