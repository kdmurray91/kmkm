// A thin C++ wrapper around kseq.h
#ifndef KMSEQ_HH_T0ZK1CKV
#define KMSEQ_HH_T0ZK1CKV

#include <zlib.h>
#include <string>
#include "kseq.h"

namespace kmkm
{
namespace kmseq
{

using namespace std;
 
KSEQ_INIT(gzFile, gzread)  
  
struct KSeq
{
    string sequence;
    bool valid;
};



class KSeqReader
{
public:
    KSeqReader(const string &filename)
    {
        _fp = gzopen(filename.c_str(), "r");
        if (_fp == NULL) {
            throw "Could not open file";
        }
        _seq = kseq_init(_fp);
    }

    ~KSeqReader()
    {
        kseq_destroy(_seq);
        gzclose(_fp);
    }

    bool next_read(string &sequence)
    {
        int l = kseq_read(_seq);
        if (l < 1) return false;
        sequence.assign(_seq->seq.s, _seq->seq.l);
        return true;
    }

    size_t next_chunk(vector<string> &sequences, size_t size)
    {
        size_t count = 0;
        for (auto seq: sequences) {
            if (count >= size) break;
            if (!this->next_read(seq)) break;
        }
        sequences.resize(count);
        return count;
    }

protected:
    gzFile _fp;  
    kseq_t *_seq;  
};
    
} /* end namespace kmseq */ 
} /* end namespace kmkm */ 

#endif /* end of include guard: KMSEQ_HH_T0ZK1CKV */

// vim:set et sw=4 ts=4:
