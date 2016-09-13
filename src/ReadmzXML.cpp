/*
Read mzXML files by the code from SpectraST Ramp Module.
Attention:
    a. zlib.h and  glob.h are used;
    b. Do not compile by Visual Studio compiler
    c. This project is upload to MEGA disk


Author: L. WU
Date: July 20

*/


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


void TestPeakList();

void GetSwathMS2Rand(const string &inputfile1, const string &inputfile2, string &outputbinaryfile,
                     SWATHmzXMLReader &smr, string mergeMethod);

void GetSwathMS2InOne(const string &inputfile1, const string &inputfile2, string &outputbinaryfile,
                      SWATHmzXMLReader &smr, string mergeMethod);

void GetSwathMS1DotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile,
                           SWATHmzXMLReader &smr, string mergeMethod);

// I hope this function is added.
void GetSwathMS2DotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile,
                           SWATHmzXMLReader &smr, string mergeMethod);


void
GetSwathDotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile, SWATHmzXMLReader &smr,
                   string mergeMethod);

void CalculateMSnDotProduct(const string &inputfile1, const string &inputfile2, const string &outputbinaryfile,
                            const string mslevel);

void getRTMatches(const vector<vector<long>> &Matches, const vector<PeakList *> &c, const vector<PeakList *> &d,
                  vector<vector<double>> &RTMatches);

void outputRTMatches(const vector<vector<double>> &RTMatches, const string &RTMatchesFile);


void RunDTWRTMatch(const vector<PeakList *> &c, const vector<PeakList *> &d, string &outputbinaryfile,
                   Matrix &dpMatrix);

vector<vector<long> > DTW(Matrix &dp) {
    //cout << "[Info] Runnint DTW" << endl;
    long r = dp.getM_row(), c = dp.getM_col();
    Matrix score(r, c);

    // update upper-left conner
    score.set(0, 0, dp.get(0, 0));

    // update first column
    for (int i = 1; i < r; ++i) {
        double val = dp.get(i, 0) < score.get(i - 1, 0) ? score.get(i - 1, 0) : dp.get(i, 0);
        score.set(i, 0, val);
    }
    //score.Print();
    // update first row
    for (int i = 1; i < c; ++i) {
        double val = dp.get(0, i) < score.get(0, i - 1) ? score.get(0, i - 1) : dp.get(0, i);
        score.set(0, i, val);
    }

    //score.Print();

    // update matrix
    for (int i = 1; i < r; ++i) {
        for (int j = 1; j < c; ++j) {
            double s1 = score.get(i - 1, j - 1);
            double s2 = score.get(i - 1, j);
            double s3 = score.get(i, j - 1);
            double s4 = s1 + dp.get(i, j);
            if (s4 >= s3 && s4 >= s2) // s4 > s2, s3
            {
                score.set(i, j, s4);
            }
            else if (s3 >= s2) // s3 [s4] s2 [s4]
            {
                score.set(i, j, s3);
            }
            else {
                score.set(i, j, s2);
            }
            //score.Print();
        }

    }

    //score.Print();
    //cout << "[Info] Start tracing back" << endl;
    vector<vector<long> > Matches;
    // trace back
    int i = r - 1, j = c - 1;
    while (i >= 1 && j >= 1) {
        double s1 = score.get(i - 1, j - 1);
        double s2 = score.get(i - 1, j);
        double s3 = score.get(i, j - 1);
        double s4 = s1 + dp.get(i, j);
        if (s4 >= s3 && s4 >= s2) {
            vector<long> match;
            match.push_back(i);
            match.push_back(j);
            Matches.push_back(match);
            i--;
            j--;
        }
        else if (s3 >= s2) {
            j--;
        }
        else {
            i--;
        }
    }


    //cout << "[Info] Finish DTW" << endl;
    return Matches;
}

void print_usage() {
    cout << "usage" << endl;
}

int main(int argc, char *argv[]) {
    SimpleTimer st("mzXMLReader(CalculateDotProduct)");
    // set usage by boost program options
    namespace po = boost::program_options;
    po::options_description mzXMLReader_po("usage");
    // declare two options
    mzXMLReader_po.add_options()
            ("help,h", po::bool_switch()->default_value(false), "print help information")
            ("inputfile,i", po::value<std::vector<string> >(), "input mzXML")

            ("ms1", po::bool_switch()->default_value(false), "only calculate dot product for ms1")
            ("ms2in1", po::bool_switch()->default_value(false), "add all the ms2 into one")
            ("ms1+ms2", po::bool_switch()->default_value(false), "only calculate dot product for ms1 and ms2")
            ("ms2", po::bool_switch()->default_value(false), "only calculate dot product for ms2")
            ("ms2rand", po::bool_switch()->default_value(false),
             "Calculate dot product of product vector of MS2 matrices with a random vector")
            ("similarity", po::value<string>()->default_value("1"),
             "1, Dot product; 2, Lp norm (p=2 by default); 3, Shared peak count; 4, Pearson Correlation Coefficient")
            ("pnorm", po::value<double>()->default_value(2), "Configuation of p for Lp norm. p = 2 by default")
            ("mergeMethod", po::value<string>()->default_value("Arithmetic"), "Arithmetic, Harmonic, Geometric")
            ("filters", po::value<string>()->default_value("Dummy"),
             "Dummy; PeakSquareRoot; PeakToProb; AboveMean; AboveMedian; Top50")
            ("ms2filters", po::value<string>()->default_value("Top20+NormalizeByMax"), "Dummy; NormalizeByMax; Top20")
            ("verbose", po::bool_switch()->default_value(false), "print more progress information");

    // I would like to check how the merge works in git tool.

    // parse parameters from argc, and argv.

    po::variables_map vm;
//    po::store(po::parse_command_line(argc,argv,mzXMLReader_po),vm);
//    po::notify(vm);

    po::positional_options_description p;
    p.add("inputfile", -1);

    po::store(po::command_line_parser(argc, argv).options(mzXMLReader_po).positional(p).run(), vm);
    po::notify(vm);

    if (vm["help"].as<bool>()) {
        cout << "Build: " << __DATE__ << " " << __TIME__ << endl;
        cout << "Read mzXML and calculate dot product..." << endl;
        cout << mzXMLReader_po << endl;
        return 0;
    }
    if (vm["ms1"].as<bool>()) {
        cout << "only output ms1" << endl;
    }

    vector<string> filelist;

    if (vm.count("inputfile")) {

        filelist = vm["inputfile"].as<vector<string> >();
//        cout << "Input " << filelist.size() << " mzXML files" << endl;
//        for (int i = 0; i < filelist.size(); ++i) {
//            cout << filelist[i] << endl;
//        }

    }
    else {
        cout << "Please input two mzXML files." << endl;
        cout << "Example:" << endl;
        cout << argv[0] << " inputmzXML1 inputmzXL2 " << endl;
        return 0;
    }
    //TestPeakList();
    //exit(0);
    //TestDTW();

    try {
        //TestMatrix();
        string SLASH = string("/");
        if (argc == 1) {

            cout << argv[0] << " input1.mzXML[string] input2.mzXML[string] msLevel[int]" << endl;
            exit(0);
        }
        string inputpath1, inputname1;
        string inputfile1 = filelist[0], inputfile2 = filelist[1];
        int found = inputfile2.find_last_of(SLASH);
        string MSLEVEL = "1";
        if (argc == 4) {
            MSLEVEL = argv[3];
        }
        else if (argc > 4) {
            MSLEVEL = "z";
        }

        if (string::npos == found) {
            inputpath1 = "./";
            inputname1 = inputfile1;
        }
        else {
            inputpath1 = inputfile1.substr(0, found + 1);
            inputname1 = inputfile1.substr(found + 1);
        }
        //string outputfile = inputfile2+"_"+inputname1+"_MS"+MSLEVEL+".csv";
        // We do not need this .
        string outputbinaryfile =
                inputfile2.substr(0, inputfile2.length() - 6) + "_" + inputname1.substr(0, inputname1.length() - 6) +
                "_MS" + MSLEVEL + "_binary.dat";


        //CalculateMSnDotProduct(inputfile1, inputfile2, outputbinaryfile,MSLEVEL);


        outputbinaryfile =
                inputfile2.substr(0, inputfile2.length() - 6) + "_" + inputname1.substr(0, inputname1.length() - 6) +
                ".dat";
        // Attention
        SWATHmzXMLReader mr;
        SimlarityMetric *simMetric = NULL;
        if (vm.count("similarity")) {
//            cout << "found similarity" << endl;
        }
        if (vm.count("pnorm")) {
//            cout << "found pnorm" << endl;
        }
        string sim = "";
        if (vm["similarity"].as<string>() == "1") {
            simMetric = new DotProduct();
            sim = "Dot product";
        }
        else if (vm["similarity"].as<string>() == "2") {
            double p = vm["pnorm"].as<double>();
            cout << "p" << p << endl;
            simMetric = new PNorm(p);
//            cout << "pnorm" << endl;
            sim = "Lp-Norm (p=" + to_string(p) + ")";
        }
        else if (vm["similarity"].as<string>() == "3") {
            simMetric = new SPC();
            sim = "Shared peak count";
//            cout << "spc" <<endl;
        }
        else if (vm["similarity"].as<string>() == "4") {
            simMetric = new PCC();
            sim = "Pearson Correlation Coefficient";
        }
        else if (vm["similarity"].as<string>() == "5") {
            simMetric = new JSD();
            sim = "Jensen-Shannon Divergence";
        }

        else {
//            cout << "wrong" << endl;

            simMetric = new DotProduct();
            sim = "Default(DotProduct)";

        }
        cout << "[Info] Using " << sim << " to calculate similarity of two spectra." << endl;
        string mergeMethod = "Arithmetic";
        if (vm.count("mergeMethod")) {
            mergeMethod = vm["mergeMethod"].as<string>();
        }
        cout << "[Info] Using " << mergeMethod << " method to merge meultiple matrix." << endl;
        mr.setSimlarityMetric(simMetric);
        if (vm.count("filters")) {
            mr.setPeakFilters(vm["filters"].as<string>());
        }

        if (vm.count("ms2filters")) {
            mr.setMs2PeakListFilters(vm["ms2filters"].as<string>());
        }

//        cout << "error" << endl;
//        cout << MSLEVEL << endl;
        //cout << "We are doing denoise now" << endl;
        if (MSLEVEL == "1" || vm["ms1"].as<bool>()) {
            cout << "[Info] MS1 based DTW" << endl;
            GetSwathMS1DotProduct(inputfile1, inputfile2, outputbinaryfile, mr, mergeMethod);
        }

        if (MSLEVEL == "x" || vm["ms1+ms2"].as<bool>()) {
//            cout << "MS@=-----------------" << endl;
            cout << "[Info] MS1+MS2 based DTW" << endl;
            GetSwathDotProduct(inputfile1, inputfile2, outputbinaryfile, mr, mergeMethod);
        }

        // I hope those lines are added into the new version.
        if (MSLEVEL == "2" || vm["ms2"].as<bool>()) {
//            cout << "Ms2 aline" << endl;
            cout << "[Info] MS2 based DTW" << endl;
            GetSwathMS2DotProduct(inputfile1, inputfile2, outputbinaryfile, mr, mergeMethod);

        }

        if (vm["ms2in1"].as<bool>()) {
            cout << "[Info] Add MS2 spectra together" << endl;
            GetSwathMS2InOne(inputfile1, inputfile2, outputbinaryfile, mr, mergeMethod);
        }

        if (vm["ms2rand"].as<bool>()) {
            cout << "[Info] Random Vector Multiplication of MS2 matrix " << endl;
            GetSwathMS2Rand(inputfile1, inputfile2, outputbinaryfile, mr, mergeMethod);
        }







        //string rtmatchesfile = outputbinaryfile + "1.dp.rt.matches";
        //SWATHDenoise(inputfile1, inputfile2,rtmatchesfile);




    }
    catch (const char *s) {
        cout << "[Error] " << s << endl;
    }
    catch (...) {
        cout << "[Error] Unexpected error occurred! Program will terminate." << endl;
    }
    return 0;
}

void GetSwathMS2Rand(const string &inputfile1, const string &inputfile2, string &outputbinaryfile,
                     SWATHmzXMLReader &smr, string mergeMethod) {
    mzXMLReader mr("2");
    int cycle1 = smr.GetSwathCycleNum(inputfile1);
    int cycle2 = smr.GetSwathCycleNum(inputfile2);

    if (cycle1 != cycle2) {
        cout << "[Error] SWATH cycle in two input file must be the same." << endl;
        throw "Inconsistent swath cycle setting";

    }

    vector<PeakList *> ms2_all_a = mr.ReadmzXMLToPeakLists(inputfile1);
    vector<PeakList *> ms2_all_b = mr.ReadmzXMLToPeakLists(inputfile2);
    string curOutputfile = outputbinaryfile + "_ms2rand.dp";

    int threads = getProperThreads();
//    cout << "[Info] Calculating dot product of swath " << i << endl;
    smr.RandMultiThreadDotProduct(ms2_all_a, ms2_all_b, curOutputfile, threads, cycle1);

    Matrix dp(curOutputfile + "");
    vector<PeakList *> ms1_a = smr.ReadmzXMLToPeakLists(inputfile1);
    vector<PeakList *> ms1_b = smr.ReadmzXMLToPeakLists(inputfile2);
    RunDTWRTMatch(ms1_a, ms1_b, curOutputfile, dp);

    dp.outputAsText(curOutputfile + ".txt");

    releaseVectorPeakListPtr(ms2_all_a);
    releaseVectorPeakListPtr(ms2_all_b);

    releaseVectorPeakListPtr(ms1_a);
    releaseVectorPeakListPtr(ms1_b);

}

void GetSwathMS2InOne(const string &inputfile1, const string &inputfile2, string &outputbinaryfile,
                      SWATHmzXMLReader &smr, string mergeMethod) {
    //    SWATHmzXMLReader mr;
    int cycle1 = smr.GetSwathCycleNum(inputfile1);
    int cycle2 = smr.GetSwathCycleNum(inputfile2);
    if (cycle1 != cycle2) {
        cout << "[Error] SWATH cycle in two input file must be the same." << endl;
        throw "Inconsistent swath cycle setting";

    }
    // Reading first MS2
    int start_swath = 2;
    int topN = 1;
    cout << "[Info] Keep top " << topN << " peaks from each MS2..." << endl;
    vector<PeakList *> c = smr.ReadSWATHmzXMLToPeakLists(inputfile1, start_swath, cycle1);

    vector<PeakList *> d = smr.ReadSWATHmzXMLToPeakLists(inputfile2, start_swath, cycle1);
    for (int j = 0; j < c.size(); j++) {
        if (c[j] != NULL) {
            smr.FilterMS2PeakList(c[j]);
//            c[j]->KeepTopN(topN);
//            c[j]->NormalizedToSum();
        }
        if (d[j] != NULL) {
            smr.FilterMS2PeakList(d[j]);
//            d[j]->KeepTopN(topN);
//            d[j]->NormalizedToSum();
        }

    }

    for (int i = start_swath + 1; i <= cycle1; i += 1) {
        SimpleTimer st("Loading " + to_string(i) + "/" + to_string(cycle1) + " SWATH...");
        vector<PeakList *> a = smr.ReadSWATHmzXMLToPeakLists(inputfile1, i, cycle1);
        vector<PeakList *> b = smr.ReadSWATHmzXMLToPeakLists(inputfile2, i, cycle1);
        for (int j = 0; j < a.size(); ++j) {
            if (a[j] != NULL) {
                smr.FilterMS2PeakList(a[j]);
//                a[j]->KeepTopN(topN);
//                a[j]->NormalizedToSum();
                if (c[j] == NULL) {
                    c[j] = new PeakList(*a[j]);

                }
                else {
                    c[j]->addPeakList(*a[j]);
                }

            }
            if (b[j] != NULL) {
                smr.FilterMS2PeakList(b[j]);
//                b[j]->KeepTopN(topN);
//                b[j]->NormalizedToSum();
                if (d[j] == NULL) {
                    d[j] = new PeakList(*b[j]);
                }
                else {
                    d[j]->addPeakList(*b[j]);
                }

            }
        }
        releaseVectorPeakListPtr(b);
        releaseVectorPeakListPtr(a);


    }
    string curOutputfile = outputbinaryfile + "_ms2in1.dp";

    int threads = getProperThreads();
//    cout << "[Info] Calculating dot product of swath " << i << endl;
    smr.MultiThreadDotProduct(c, d, curOutputfile, threads);

    Matrix dp(curOutputfile + "");

    RunDTWRTMatch(c, d, curOutputfile, dp);
    dp.outputAsText(curOutputfile + ".txt");
    releaseVectorPeakListPtr(c);
    releaseVectorPeakListPtr(d);


}

// The following function is added by myself today , July 20. this line should be kept.
void GetSwathMS2DotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile,
                           SWATHmzXMLReader &smr, string mergeMethod) {

//    SWATHmzXMLReader mr;
    int cycle1 = smr.GetSwathCycleNum(inputfile1);
    int cycle2 = smr.GetSwathCycleNum(inputfile2);
    if (cycle1 != cycle2) {
        cout << "[Error] SWATH cycle in two input file must be the same." << endl;
        throw "Inconsist swath cycle setting";

    }

    // sum of several matrix
    vector<Matrix> dpMatrix;
    // double * res = NULL;

    // MS1: i = 1
    // MS2 SWATH: i = 2 to cycle=32
    int threads = getProperThreads();
    for (int i = 2; i <= cycle1; ++i) {
        SimpleTimer st("calculate dot product of cycle " + to_string(i));
        string curOutputfile = outputbinaryfile + to_string(i) + ".dp";

        vector<PeakList *> c = smr.ReadSWATHmzXMLToPeakLists(inputfile1, i, cycle1);
        vector<PeakList *> d = smr.ReadSWATHmzXMLToPeakLists(inputfile2, i, cycle1);
        cout << "[Info] Calculating dot product of swath " << i << endl;
        smr.MultiThreadDotProduct(c, d, curOutputfile, threads);

        Matrix dp(curOutputfile + "");

        RunDTWRTMatch(c, d, curOutputfile, dp);
        dp.outputAsText(curOutputfile + ".txt");
        releaseVectorPeakListPtr(c);
        releaseVectorPeakListPtr(d);

        dpMatrix.push_back(dp);

    }
    // to get the row and col, we have the read the file again.
    // How could we skip this step?
    long row = 0;
    long col = 0;
    vector<PeakList *> c = smr.ReadSWATHmzXMLToPeakLists(inputfile1, 1, cycle1);
    vector<PeakList *> d = smr.ReadSWATHmzXMLToPeakLists(inputfile2, 1, cycle1);
    row = c.size();
    col = d.size();

    Matrix sum_of_dpMatrix(row, col);

    sum_of_dpMatrix.mergeMatrix(dpMatrix, mergeMethod);

    outputbinaryfile += "_sum_ms2";

    RunDTWRTMatch(c, d, outputbinaryfile, sum_of_dpMatrix);

    sum_of_dpMatrix.outputBinary(outputbinaryfile + ".dp");
    sum_of_dpMatrix.outputAsText(outputbinaryfile + ".dp.txt");
    releaseVectorPeakListPtr(c);
    releaseVectorPeakListPtr(d);


}


void GetSwathMS1DotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile,
                           SWATHmzXMLReader &smr, string mergeMethod) {



//    SWATHmzXMLReader mr;
    int cycle1 = smr.GetSwathCycleNum(inputfile1);
    int cycle2 = smr.GetSwathCycleNum(inputfile2);
    if (cycle1 != cycle2) {
        cout << "[Error] SWATH cycle in two input file must be the same." << endl;

        // I hope the throw error is added,
        throw "Inconsisit swath cycle num";


    }

    // sum of several matrix
    vector<Matrix> dpMatrix;
    // double * res = NULL;

    // MS1: i = 1
    // MS2 SWATH: i = 2 to cycle
    int threads = getProperThreads();


    string curOutputfile = outputbinaryfile + to_string(1) + ".dp";

    // There is memroy leak here! please pay attention to this; fixed
    vector<PeakList *> c = smr.ReadSWATHmzXMLToPeakLists(inputfile1, 1, cycle1);
    vector<PeakList *> d = smr.ReadSWATHmzXMLToPeakLists(inputfile2, 1, cycle1);
    cout << "[Info] Calculating dot product of swath " << 1 << endl;
    smr.MultiThreadDotProduct(c, d, curOutputfile, threads);


    Matrix dp(curOutputfile + "");


    RunDTWRTMatch(c, d, curOutputfile, dp);
    dp.outputAsText(curOutputfile + ".txt");
    releaseVectorPeakListPtr(c);
    releaseVectorPeakListPtr(d);

}
// TODO:
// work here.
// load the spectra from mzXML
// merge all the ms2 with certain mz shift into one spectra M-MS2
// do the same workflow just as it is an MS1

//void GetSwathDotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile)
//{
//    // Load all the MS2 and connecting them
//    // inpput file1 & file2
//    SWATHmzXMLReader mr;
//    int cycle1 = mr.GetSwathCycleNum(inputfile1);
//    int cycle2 = mr.GetSwathCycleNum(inputfile2);
//    if (cycle1 != cycle2)
//    {
//        cout << "[Error] SWATH cycle in two input file must be the same." << endl;
//        //quit program
//    }
//
//    // double * res = NULL;
//
//    // MS1: i = 1
//    // MS2 SWATH: i = 2 to cycle
//    int threads = getProperThreads();
//
//
//    string curOutputfile = outputbinaryfile + to_string(1)+".dp";
//    cout << "checking file exist" << curOutputfile << endl;
//    // create vector of first ms2
//    vector<PeakList *> c = mr.ReadSWATHmzXMLToPeakLists(inputfile1, 2, cycle1);
//    vector<PeakList *> d = mr.ReadSWATHmzXMLToPeakLists(inputfile2, 2, cycle2);
//    for(int i = 3; i <= cycle1; ++i) {
//        //get next MS2 data, append them one by one
//        vector<PeakList *> c_temp = mr.ReadSWATHmzXMLToPeakLists(inputfile1, i, cycle1);
//        vector<PeakList *> d_temp = mr.ReadSWATHmzXMLToPeakLists(inputfile2, i, cycle1);
//        //push back PeakList* in c_temp / d_temp to c / d
//        //merge c_temp first
//        int size_c_temp = c_temp.size();
//        for (int j = 0; j < size_c_temp; ++j) {
//            PeakList *temp = c_temp[j];
//            c[j]->addPeakList(*temp, i - 2);
//        }
//        int size_d_temp = d_temp.size();
//        for (int k = 0; k < size_d_temp; ++k) {
//            PeakList *temp = d_temp[k];
//            d[k]->addPeakList(*temp, i - 2);
//        }
//    }
//    cout << "[Info] Calculating dot product of swath " << 1 << endl;
//    mr.MultiThreadDotProduct(c,d,curOutputfile,threads,cycle1);
//
//
//    Matrix dp(curOutputfile+"");
//    cout << "[Info] Checking file " << curOutputfile + ".txt" << endl;
//
//
//    RunDTWRTMatch(c, d, curOutputfile, dp);
//    dp.outputAsText(curOutputfile+".txt");
//
//}
//



void
GetSwathDotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile, SWATHmzXMLReader &smr,
                   string mergeMethod) {


//    SWATHmzXMLReader mr;
    int cycle1 = smr.GetSwathCycleNum(inputfile1);
    int cycle2 = smr.GetSwathCycleNum(inputfile2);
    if (cycle1 != cycle2) {
        cout << "[Error] SWATH cycle in two input file must be the same." << endl;
        throw "Inconsisit swath cycle num";
    }

    vector<Matrix> dpMatrix;
    // MS1: i = 1
    // MS2 SWATH: i = 2 to cycle=32
    int threads = getProperThreads();
    for (int i = 1; i <= cycle1; ++i) {
        SimpleTimer st("calculate dot product of cycle " + to_string(i));
        string curOutputfile = outputbinaryfile + to_string(i) + ".dp";

        // There is memroy leak here! please pay attention to this
        vector<PeakList *> c = smr.ReadSWATHmzXMLToPeakLists(inputfile1, i, cycle1);
        vector<PeakList *> d = smr.ReadSWATHmzXMLToPeakLists(inputfile2, i, cycle1);
        cout << "[Info] Calculating dot product of swath " << i << endl;
        smr.MultiThreadDotProduct(c, d, curOutputfile, threads, i);

        Matrix dp(curOutputfile + "");

        RunDTWRTMatch(c, d, curOutputfile, dp);
        dp.outputAsText(curOutputfile + ".txt");
        dpMatrix.push_back(dp);
        releaseVectorPeakListPtr(c);
        releaseVectorPeakListPtr(d);

    }
    long row = 0;
    long col = 0;
    vector<PeakList *> c = smr.ReadSWATHmzXMLToPeakLists(inputfile1, 1, cycle1);
    vector<PeakList *> d = smr.ReadSWATHmzXMLToPeakLists(inputfile2, 1, cycle1);
    row = c.size();
    col = d.size();

    Matrix sum_of_dpMatrix(row, col);
    sum_of_dpMatrix.mergeMatrix(dpMatrix, mergeMethod);


    RunDTWRTMatch(c, d, outputbinaryfile, sum_of_dpMatrix);

    sum_of_dpMatrix.outputBinary(outputbinaryfile + "_sum.dp");
    sum_of_dpMatrix.outputAsText(outputbinaryfile + "_sum.dp.txt");
    releaseVectorPeakListPtr(c);
    releaseVectorPeakListPtr(d);


}

void RunDTWRTMatch(const vector<PeakList *> &c, const vector<PeakList *> &d, string &outputbinaryfile,
                   Matrix &dpMatrix) {
    cout << "[Info] Run DTW ..." << endl;
    vector<vector<long> > Matches = DTW(dpMatrix);
    vector<vector<double> > RTMatches;
    getRTMatches(Matches, c, d, RTMatches);
    string RTMathesFile = outputbinaryfile + ".rt.matches";
    //cout << "output sum dp" << RTMathesFile << endl;
    outputRTMatches(RTMatches, RTMathesFile);
}

void outputRTMatches(const vector<vector<double>> &RTMatches, const string &RTMatchesFile) {
    cout << "[Info] Export DTW result to " << RTMatchesFile << endl;
    FILE *pfile = fopen(RTMatchesFile.c_str(), "w");
    for (int i = 0; i < RTMatches.size(); ++i) {
        fprintf(pfile, "%.5lf\t%.5lf\n", RTMatches[i][0], RTMatches[i][1]);
    }
    fclose(pfile);
}


void getRTMatches(const vector<vector<long>> &Matches, const vector<PeakList *> &c, const vector<PeakList *> &d,
                  vector<vector<double>> &RTMatches) {
    for (int i = 0; i < Matches.size(); ++i) {
        vector<double> rtmatch;

        rtmatch.push_back(c[Matches[i][0]]->getRTinSeconds());
        rtmatch.push_back(d[Matches[i][1]]->getRTinSeconds());

//        if(Matches[i][0] != Matches[i][1])
//            cout << rtmatch[0] << " "<< rtmatch[1] << endl;
        RTMatches.push_back(rtmatch);
    }
}


void CalculateMSnDotProduct(const string &inputfile1, const string &inputfile2, const string &outputbinaryfile,
                            const string mslevel) {
    mzXMLReader mr(mslevel);
    vector<PeakList *> a = mr.ReadmzXMLToPeakLists(inputfile1);
    vector<PeakList *> b = mr.ReadmzXMLToPeakLists(inputfile2);

    AnalysisMassDiff amd;
    vector<double> massdiff = amd.CalculateMassDifference(a, 0.5);
    int threads = 24;
    mr.MultiThreadDotProduct(a, b, outputbinaryfile, threads);

}

