//
// Created by wulong on 8/3/16.
//


#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <ctime>
#include <iomanip>
#include <boost/program_options.hpp>

#include "../librarymsms/PeakList.h"
#include "../librarymsms/SWATHmzXMLReader.h"
//#include <thread>

//#include "mzXMLReader.h"

using namespace std;


void extractMzXMLToMatrix(string inputfile) {
    mzXMLReader mr("2");
    vector<PeakList *> x = mr.ReadmzXMLToPeakLists(inputfile);
    cout << "[Info] Spectra found in " << inputfile << " : " << x.size() << endl;
    mr.exportToCSV(x, inputfile + ".csv");
    releaseVectorPeakListPtr(x);
}


int main(int argc, char *argv[]) {
    SimpleTimer st("mzXMLReader(converter)");
    // set usage by boost program options
    namespace po = boost::program_options;
    po::options_description mzXMLReader_po("usage");
    // declare two options
    mzXMLReader_po.add_options()
            ("help,h", po::bool_switch()->default_value(false), "print help information")
            ("inputfile,i", po::value<std::vector<string> >(), "input mzXML")

            ("ms1", po::bool_switch()->default_value(false), "only calculate dot product for ms1");
//            ("ms2in1", po::bool_switch()->default_value(false), "add all the ms2 into one")
//            ("ms1+ms2", po::bool_switch()->default_value(false), "only calculate dot product for ms1 and ms2")
//            ("ms2", po::bool_switch()->default_value(false), "only calculate dot product for ms2")
//            ("similarity", po::value<string>()->default_value("1"), "1, Dot product; 2, Lp norm (p=2 by default); 3, Shared peak count; 4, Pearson Correlation Coefficient")
//            ("pnorm", po::value<double>()->default_value(2), "Configuation of p for Lp norm. p = 2 by default")
//            ("mergeMethod", po::value<string>()->default_value("Arithmetic"), "Arithmetic, Harmonic, Geometric")
//            ("filters", po::value<string>()->default_value("Dummy"), "Dummy; PeakSquareRoot; PeakToProb; AboveMean; AboveMedian; Top50")
//            ("ms2filters", po::value<string>()->default_value("Top20+NormalizeByMax"), "Dummy; NormalizeByMax; Top20")
//            ("verbose", po::bool_switch()->default_value(false), "print more progress information");

    // I would like to check how the merge works in git tool.

    // parse parameters from argc, and argv.

    po::variables_map vm;
//    po::store(po::parse_command_line(argc,argv,mzXMLReader_po),vm);
//    po::notify(vm);

    po::positional_options_description p;
    p.add("inputfile", -1);

    po::store(po::command_line_parser(argc, argv).options(mzXMLReader_po).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("inputfile")) {
        vector<string> filelist = vm["inputfile"].as<vector<string> >();
        for (int i = 0; i < filelist.size(); ++i) {
            extractMzXMLToMatrix(filelist[i]);
        }
    }


    return 0;
}