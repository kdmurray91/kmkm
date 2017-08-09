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
    for (auto readfile: opt.readfiles) {
	const size_t n = ctr.consume_from(readfile);
        uint64_t sum = accumulate(ctr.counts().begin(), ctr.counts().end(), 0);
        cerr << ctr.k() << " " << readfile << " " << n << " " << ctr.nnz() <<  " " << sum << " " << endl;
    }
    ctr.save(opt.outfile);
    KmerCounter<> loaded = KmerCounter<>::load(opt.outfile);
    uint64_t sum = accumulate(loaded.counts().begin(), loaded.counts().end(), 0);
    cerr << loaded.k() << " " << loaded.nnz() <<  " " << sum << " " << endl;
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
        ("outfile", po::value<string>(&opt.outfile)->required(),
            "Output filename")
        ("readfiles", po::value<vector<string>>(&opt.readfiles)->required(),
            "Input FASTX files");
    po::positional_options_description pos;
    pos.add("outfile", 1).add("readfiles", -1);

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

    cerr << opt.ksize << endl;
    cerr << opt.cvsize << endl;
    return count_main(opt);
}
