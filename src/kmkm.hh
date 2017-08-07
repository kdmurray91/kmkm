// Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
#ifndef KMKM_HH_WY0MH8N1
#define KMKM_HH_WY0MH8N1

#include <vector>
#include <string>
#include <iostream>
#include <boost/multi_array.hpp>

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
    KmerIterator (const string &sequence, int k, bool canoncial=true, int64_t seed=3301)
        : _k(k) , _seq(sequence) , _len(sequence.size()) , _pos(0) , _canonical(canoncial)
        , _last_nthash(0) , _mask((1 << (2*k)) - 1)
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
 *  TODO: Implement CBF counting
 */
template <typename T>
class KmerCounter
{
public:
    KmerCounter (int k, size_t vecsize, bool canoncial=true, bool use_cbf=true, int64_t seed=3301)
        : _k(k)
        , _use_cbf(use_cbf)
        , _canonical(canoncial)
        , _counts(vecsize, 0)
    {
    }

    inline void consume(const string &sequence)
    {
        KmerIterator ki(sequence, _k, _canonical);
        while (!ki.finished()) {
            auto mer = ki.next_hashed();
            _counts[mer % _counts.size()]++;
        }
    }

    void consume(const vector<string> &sequences)
    {
        for (auto seq: sequences) this->consume(seq);
    }


    const vector<T>& counts() const
    {
        return _counts;
    }

    size_t nnz() const
    {
        size_t nnz = 0;
        for (auto v: _counts) if (v > 0) nnz++;
        return nnz;
    }

    double collision_rate() const
    {
        return double(this->nnz()) / double(_counts.size());
    }
    
#if 0
    void consume_from(const string &filename);
    void save(const string &filename) const;
    static KmerCounter load(const string &filename);
#endif

protected:
    const unsigned int _k;
    const bool _use_cbf;
    const bool _canonical;

    vector<T> _counts;
    vector<T> _cbf;

};

//typedef KmerCounter<uint8_t> KmerCounter;


} // end namespace kmercount
#endif /* end of include guard: KMKM_HH_WY0MH8N1 */

// vim:set et sw=4 ts=4:
