// Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KMERCOLLECTION_HH_PQKNVYSJ
#define KMERCOLLECTION_HH_PQKNVYSJ

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <exception>

#include <boost/filesystem.hpp>
#define ARMA_USE_HDF5
#include <armadillo>

#include "kmkm.hh"


namespace kmercollection
{

using namespace kmkm;
using namespace arma;
using namespace std;
namespace bfs = boost::filesystem;


class KmerCollection
{
private:
    size_t _top_n;
    vector<string> _samplenames;
    fmat _counts;

public:

    KmerCollection(size_t top_n)
        : _top_n(top_n)
    {}

    inline const fmat &counts()
    {
        return _counts;
    }

    inline const vector<string> &names()
    {
        return _samplenames;
    }

    inline string normalise_samplename(const string&filepath)
    {
        bfs::path name = bfs::path(filepath).filename();
        for (const string &ext: {".gz", ".kmr"}) {
            if (name.extension() == ext) {
                name = name.stem();
            }
        }
        return name.string();
    }

    void add_samples(const vector<string> &samples)
    {
        const size_t nsamp = samples.size();

        _samplenames.clear();
        _samplenames.assign(nsamp, "");
        _counts.set_size(_top_n, nsamp);

        cerr << "Loading individual countfiles into count matrix..." << endl;
        #pragma omp parallel for schedule(dynamic) default(shared)
        for (size_t j = 0; j < nsamp; j++) {
            const auto &countfile = samples[j];
            KmerCounter<> ctr(countfile);
            const auto &cv = ctr.counts();
            const size_t N = _top_n == 0 ? cv.size() : _top_n;
            for (size_t i = 0; i < N; i++) {
                _counts(i, j) = cv[i];
            }

            const string name = normalise_samplename(countfile);
            _samplenames[j] = name;
            #pragma omp critical
            {
                cout << "  - " << name << endl;
            }
        }
    }

    void save(const string &basename)
    {
        const string &bn = basename;
        _counts.save(bn + ".counts", hdf5_binary_trans);
        
        ofstream nf(bn + ".samples");
        for (const auto &s: _samplenames) {
            nf << s << "\n";
        }

    }

    /*! \brief Attempts to load a kmer collection from disc
     *
     * Loads count matrix, and list of samples. Checks the shape of the count
     * matrix, and (if samples is given) check the list of sample identifiers
     * against those saved on disc. 
     *
     * \param basename Prefix of count and sample files on disc
     * \param samples list of sample names or count files
     * \return true on success, else false
     */
    bool load(const string &basename, vector<string> samples={})
    {
        const string &bn = basename;

        if (!bfs::exists(bn + ".counts") || !bfs::exists(bn + ".samples")) return false;

        _counts.load(bn + ".counts", hdf5_binary_trans);

        if (_counts.n_rows != _top_n) return false;
        if (samples.size() > 0 && _counts.n_cols != samples.size()) return false;

        for (auto &s: samples) s = normalise_samplename(s);

        
        _samplenames.clear();
        ifstream nf(bn + ".samples");
        for (string line; getline(nf, line); ) {
            _samplenames.push_back(line);
        }
        if (samples.size() > 0 && samples != _samplenames) {
            return false;
        }
        if (_counts.n_cols != _samplenames.size()) {
            return false;
        }
        return true;
    }

    
#if 0
    bool add_sample(string &filename)
    {
    }

    bool insert_sample(string &filename, size_t idx)
    {
    }
#endif
};

    
} /* kmercollection */ 

#endif /* end of include guard: KMERCOLLECTION_HH_PQKNVYSJ */
