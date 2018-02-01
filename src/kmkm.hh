// Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KMKM_HH_WY0MH8N1
#define KMKM_HH_WY0MH8N1

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <functional>

//#include <boost/multi_array.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include "kmseq.hh"


using namespace std;

namespace kmkm {

/**
* @brief: 64 bit integer hash function
*
* @param: uint64_t x
*
* @return: uint64_t
*/
static inline uint64_t inthash64(uint64_t x)
{
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
}

static inline uint64_t kmer_revcomp(uint64_t x, int k)
{
    // Reverse 2-bit sections
    x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>  2;
    x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4;
    x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8;
    x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16;
    x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32;
    // Complement
    x = ~x;
    // Shift back within bitmask
    x >>= 64 - 2 * k;
    return x;
}

/*! \class KmerIterator
 *  \brief Iterator over the kmers of a DNA Sequence
 *
 *  This ought be used like:
 *
 *  KmerIterator i("ACGTAGAG", 2);
 *  while (!i.finished()) {
 *      nt = i.next();
 *  }
 */
class KmerIterator
{
public:
    KmerIterator (const string &sequence, int k, bool canonical=true)
        : _k(k) , _seq(sequence) , _len(sequence.size()) , _pos(0) , _canonical(canonical)
        , _last_nthash(0) , _mask((UINT64_C(1) << (2*k)) - 1)
    {
    }


    /*! \brief Returns the next k-mer in sequence
     *
     *  K-mers are returned as 2-bit encoded nucleotide sequences, such that
     *  AAAA... = 0, and TTTT... = 4^k.
     *
     * \return 2-bit encoded kmer as uint64_t
     */
    inline uint64_t next()
    {
        ssize_t skip = 0;
        do {
            if (this->finished()) {
                /* Prevent an out-of-bounds read */
                return 0;
            }
            if (skip > 0) {
                skip--;
            }
            char nucl = _seq[_pos++] & 0x5f; // Force uppercase
            uint64_t n = 0; // Numeric nucleotide repr
            switch (nucl) {
                case 'A':
                    n = 0;
                    break;
                case 'C':
                    n = 1;
                    break;
                case 'G':
                    n = 2;
                    break;
                case 'T':
                    n = 3;
                    break;
                default:
                    skip = _k;
                    break;
            }
            _last_nthash = ((_last_nthash << 2) | n) & _mask;
        } while (skip > 0 || _pos < _k);
        if (_canonical) {
            return min(_last_nthash, kmer_revcomp(_last_nthash, _k));
        } else {
            return _last_nthash;
        }
    }


    /*! \brief Returns the next k-mer in sequence, after hashing
     *
     *  K-mers are obtained from this->next(), and passed through an integer
     *  hash function to scramble lexographic order.
     *
     * \return Hashed kmer as uint64_t
     */
    uint64_t next_hashed()
    {
        return inthash64(this->next());
        // std::hash not guranteed to be identical between implementations/runs
        // of program. So use our own.
        //return std::hash<uint64_t>{}(this->next());
    }

    /*! \brief Check if iterator has compeleted iteration
     *
     * \return true when finished
     */
    inline bool finished()
        
    {
        return _len < _k || _pos >= _len;
    }

    /*! \brief Resets iterator to start position
     */
    void rewind()
    {
        _pos = 0;
        _last_nthash = 0;
    }

    /*! @brief Gives the (maximum) number of k-mers
     *
     * Note: This is length - ksize + 1. The actual number of k-mers obtained
     * can be lower if some kmers are invalid (e.g. because the sequence has
     * ambiguous nucleotides).
     *
     * @return The maximum number of k-mers
     */
    inline size_t size() const
    {
        return max(ssize_t(_len) - _k + 1, ssize_t(0));
    }


private:
    const unsigned int _k;
    const string &_seq;
    const size_t _len;
    size_t _pos;
    bool _canonical;
    uint64_t _last_nthash;
    uint64_t _mask;
};


/**********************************************************************
*                            KmerCounter                             *
**********************************************************************/



/*! \class KmerCounter
 *  \brief Counting Bloom Filter-based k-mer counter
 *
 */
template <typename ElType = uint8_t>
class KmerCounter
{
public:
    KmerCounter()
        : _k(0)
        , _cbf_tables(0)
        , _canonical(false)
    { }

    KmerCounter (int k, size_t vecsize, bool canonical=true, size_t cbf_tables=0)
        : _k(k)
        , _cbf_tables(cbf_tables)
        , _canonical(canonical)
        , _counts(vecsize, 0)
        , _cbf((vecsize/2) * _cbf_tables, 0)
    {
    }

    KmerCounter(KmerCounter&& x)
        : _k(x._k)
        , _cbf_tables(x._cbf_tables)
        , _canonical(x._canonical)
        , _counts(std::move(x._counts))
        , _cbf(std::move(x._cbf))
    { }

    KmerCounter(const string &filename)
        : _k(0)
        , _cbf_tables(0)
        , _canonical(false)
    {
        this->load(filename);
    }

    inline KmerCounter & operator=(KmerCounter&& x)
    {
        if (this != &x) {
            _k = x._k;
            _cbf_tables = x._cbf_tables;
            _canonical = x._canonical;
            _counts = std::move(x._counts);
            _cbf = std::move(x._cbf);
        }
    }

    inline void count(uint64_t hashed_kmer)
    {
        const size_t cvidx = hashed_kmer % _counts.size();
        if (_counts.size() == 0) {
            throw "CBF not initialised";
        }


        ElType current;
        if (_cbf_tables > 0) {
            // If using CBF, do a count-min operation to determine current
            // count.
            const size_t cbfsize = _counts.size() / 2;
            const size_t cbfidx = hashed_kmer % cbfsize;
            current = numeric_limits<ElType>::max();
            for (size_t t = 0; t < _cbf_tables; t++) {
                ElType tval = _cbf[t * cbfsize + cbfidx];
                if (tval < current) {
                    current = tval;
                }
            }
        } else {
            // Otherwise, read current count directly from the CV
            current = _counts[cvidx];
        }
        _counts[cvidx] = current + 1;
    }

    inline void consume(const string &sequence)
    {
        KmerIterator ki(sequence, _k, _canonical);
        while (!ki.finished()) {
            this->count(ki.next_hashed());
        }
    }

    void consume(const vector<string> &sequences)
    {
        for (auto seq: sequences) this->consume(seq);
    }

    inline const vector<ElType>& counts() const
    {
        return _counts;
    }

    inline const ElType * data() const
    {
        return _counts.data();
    }

    inline size_t nnz() const
    {
        size_t nnz = 0;
        for (auto v: _counts) if (v > 0) nnz++;
        return nnz;
    }

    inline double collision_rate() const
    {
        if (_counts.size() == 0) return -1;
        return double(this->nnz()) / double(_counts.size());
    }

    inline int k() const
    {
        return _k;
    }

    void save(const string &filename)
    {
        using namespace boost::iostreams;
        ofstream fp(filename, ios_base::out | ios_base::binary);
        filtering_streambuf<output> out;
        if (boost::algorithm::ends_with(filename, ".gz")) {
            out.push(gzip_compressor());
        }
        out.push(fp);
        boost::archive::binary_oarchive ar(out);
        ar << *this;
    }

    void load(const string &filename)
    {
        using namespace boost::iostreams;
        ifstream fp(filename, ios_base::in | ios_base::binary);
        filtering_streambuf<input> in;
        if (boost::algorithm::ends_with(filename, ".gz")) {
            in.push(gzip_decompressor());
        }
        in.push(fp);

        boost::archive::binary_iarchive ar(in);
        ar >> *this;
    }

    size_t consume_from(const string &filename)
    {
        kmseq::KSeqReader seqs(filename);
        size_t n = 0;
        for (kmseq::KSeq seq; seqs.next_read(seq);) {
            this->consume(seq.seq);
            n++;
        }
        return n;
    }

protected:
    const unsigned int _k;
    const size_t _cbf_tables;
    const bool _canonical;
    vector<ElType> _counts;
    vector<ElType> _cbf;

    // Serialization
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & const_cast<unsigned int &>(_k);
        ar & const_cast<bool &>(_canonical);
        ar & _counts;
    }
};


} // end namespace kmercount
#endif /* end of include guard: KMKM_HH_WY0MH8N1 */

// vim:set et sw=4 ts=4:
