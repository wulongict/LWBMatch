/*
Read mzXML files by the code from SpectraST Ramp Module.
Attention:
    a. zlib.h and  glob.h are used;
    b. Do not compile by Visual Studio compiler
    c. This project is upload to MEGA disk




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

void GetSwathMS1DotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile);
void GetSwathDotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile);

void CalculateMSnDotProduct(const string &inputfile1, const string &inputfile2, const string &outputbinaryfile,
                            const string mslevel);

void getRTMatches(const vector<vector<long>> &Matches, const vector<PeakList *> &c, const vector<PeakList *> &d,
                  vector<vector<double>> &RTMatches);

void outputRTMatches(const vector<vector<double>> &RTMatches, const string &RTMathesFile);


void RunDTWRTMatch(const vector<PeakList *> &c, const vector<PeakList *> &d, string &outputbinaryfile,
                   Matrix &dpMatrix);

vector<vector<long> >  DTW(Matrix & dp)
{
    cout << "[Info] Runnint DTW" << endl;
    long r = dp.getM_row(), c = dp.getM_col();
    Matrix score(r,c);

    // update upper-left conner
    score.set(0,0, dp.get(0,0));

    // update first column
    for (int i = 1; i < r; ++i) {
        double val = dp.get(i,0) < score.get(i-1,0) ? score.get(i-1,0) : dp.get(i,0);
        score.set(i,0,val);
    }
    //score.Print();
    // update first row
    for (int i = 1; i < c; ++i) {
        double val = dp.get(0,i) < score.get(0,i-1) ? score.get(0,i-1) : dp.get(0,i);
        score.set(0,i,val);
    }

    //score.Print();

    // update matrix
    for (int i = 1; i < r; ++i) {
        for (int j = 1; j < c; ++j) {
            double s1 = score.get(i-1,j-1);
            double s2 = score.get(i-1,j);
            double s3 = score.get(i,j-1);
            double s4 = s1 + dp.get(i,j);
            if (s4 >= s3 && s4 >= s2) // s4 > s2, s3
            {
                score.set(i,j,s4);
            }
            else if (s3 >= s2) // s3 [s4] s2 [s4]
            {
                score.set(i,j,s3);
            }
            else
            {
                score.set(i,j,s2);
            }
            //score.Print();
        }
    }

    //score.Print();
    cout << "[Info] Start tracing back" << endl;
    vector<vector<long> > Matches;
    // trace back
    int i = r-1, j = c-1;
    while (i>=1 && j >= 1)
    {
        double s1 = score.get(i-1,j-1);
        double s2 = score.get(i-1,j);
        double s3 = score.get(i,j-1);
        double s4 = s1 + dp.get(i,j);
        if (s4 >= s3 && s4 >= s2 )
        {
            vector<long> match;
            match.push_back(i);
            match.push_back(j);
            Matches.push_back(match);
            i --;
            j --;
        }
        else if(s3 >= s2)
        {
            j --;
        }
        else
        {
            i--;
        }
    }


    cout << "[Info] Finish DTW" << endl;
    return Matches;
}
void print_usage()
{
    cout << "usage" << endl;
}

int main(int argc, char * argv[])
{
    SimpleTimer st;
    // set usage by boost program options
    namespace po = boost::program_options;
    po::options_description mzXMLReader_po("usage");
    // declare two options
    mzXMLReader_po.add_options()
            ("help,h",po::bool_switch()->default_value(false),"print help information")
            ("inputfile,i",po::value< std::vector<string> >(),"input mzXML")
            ("ms1",po::bool_switch()->default_value(false),"only calculate dot product for ms1")

            ;
    // parse parameters from argc, and argv.

    po::variables_map vm;
//    po::store(po::parse_command_line(argc,argv,mzXMLReader_po),vm);
//    po::notify(vm);

    po::positional_options_description p;
    p.add("inputfile",-1);

    po::store(po::command_line_parser(argc,argv).options(mzXMLReader_po).positional(p).run(),vm);
    po::notify(vm);

    if (vm["help"].as<bool>())
    {
        cout << "Build: " << __DATE__ << " " << __TIME__ << endl;
        cout << "Read mzXML and calculate dot product..." << endl;
        cout << mzXMLReader_po << endl;
        exit(0);
    }
    if(vm["ms1"].as<bool>())
    {
        cout << "only output ms1" << endl;
    }

    vector<string> filelist;

    if(vm.count("inputfile"))
    {

        filelist = vm["inputfile"].as<vector<string> >();
        cout << "Input "<< filelist.size() <<" mzXML files" << endl;
        for (int i = 0; i < filelist.size(); ++i) {
            cout << filelist[i] << endl;
        }

    }
    else
    {
        cout << "Please input two mzXML files." << endl;
        cout << "Example:" << endl;
        cout << argv[0] << " inputmzXML1 inputmzXL2 " << endl;
        exit(0);
    }
    //TestPeakList();
    //exit(0);
    //TestDTW();

    try{
        //TestMatrix();
        string SLASH = string("/");
        if(argc == 1)
        {

            cout << argv[0] << " input1.mzXML[string] input2.mzXML[string] msLevel[int]" << endl;
            exit(0);
        }
        string inputpath1,inputname1;

        string inputfile1 = filelist[0], inputfile2 = filelist[1];
        string MSLEVEL = argv[3];
        int found = inputfile2.find_last_of(SLASH);
        if(string::npos==found)
        {
            inputpath1 = "./";
            inputname1 = inputfile1;
        }
        else
        {
            inputpath1 = inputfile1.substr(0,found+1);
            inputname1 = inputfile1.substr(found+1);
        }
        //string outputfile = inputfile2+"_"+inputname1+"_MS"+MSLEVEL+".csv";
        string outputbinaryfile = inputfile2.substr(0,inputfile2.length()-6)+"_"+inputname1.substr(0,inputname1.length()-6)+"_MS"+MSLEVEL+"_binary.dat";


        //CalculateMSnDotProduct(inputfile1, inputfile2, outputbinaryfile,MSLEVEL);


        outputbinaryfile = inputfile2.substr(0,inputfile2.length()-6)+"_"+inputname1.substr(0,inputname1.length()-6)+".dat";
        // Attention
        //cout << "We are doing denoise now" << endl;
        if (MSLEVEL == "1" || vm["ms1"].as<bool>())
        {
            GetSwathMS1DotProduct(inputfile1,inputfile2,outputbinaryfile);
        }
        else{
            GetSwathDotProduct(inputfile1, inputfile2, outputbinaryfile);
        }


        //string rtmatchesfile = outputbinaryfile + "1.dp.rt.matches";
        //SWATHDenoise(inputfile1, inputfile2,rtmatchesfile);




    }
    catch (const char * s)
    {
        cout << "[Error] " << s << endl;
    }
    catch (...)
    {
        cout << "[Error] Unexpected error occurred! Program will terminate." << endl;
    }
    return 0;
}



void GetSwathMS1DotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile) {


    SWATHmzXMLReader mr;
    int cycle1 = mr.GetSwathCycleNum(inputfile1);
    int cycle2 = mr.GetSwathCycleNum(inputfile2);
    if (cycle1 != cycle2)
    {
        cout << "[Error] SWATH cycle in two input file must be the same." << endl;

    }

    // sum of several matrix
    vector<Matrix> dpMatrix;
    // double * res = NULL;

    // MS1: i = 1
    // MS2 SWATH: i = 2 to cycle
    int threads = getProperThreads();

    SimpleTimer st("calculate dot product of cycle " + to_string(1));
    string curOutputfile = outputbinaryfile + to_string(1)+".dp";


    if (!isFileExist(curOutputfile))// never skip this step
    {
        // There is memroy leak here! please pay attention to this
        vector<PeakList*> c = mr.ReadSWATHmzXMLToPeakLists(inputfile1,1,cycle1);
        vector<PeakList*> d = mr.ReadSWATHmzXMLToPeakLists(inputfile2,1,cycle1);
        cout << "[Info] Calculating dot product of swath " << 1 << endl;
        mr.MultiThreadDotProduct(c,d,curOutputfile,threads);

    }
    Matrix dp(curOutputfile+"");
    if (! isFileExist(curOutputfile+".txt"))
    {
        vector<PeakList*> c = mr.ReadSWATHmzXMLToPeakLists(inputfile1,1,cycle1);
        vector<PeakList*> d = mr.ReadSWATHmzXMLToPeakLists(inputfile2,1,cycle1);
        RunDTWRTMatch(c, d, curOutputfile, dp);
        dp.outputAsText(curOutputfile+".txt");
    }







}

void GetSwathDotProduct(const string &inputfile1, const string &inputfile2, string &outputbinaryfile) {


    SWATHmzXMLReader mr;
    int cycle1 = mr.GetSwathCycleNum(inputfile1);
    int cycle2 = mr.GetSwathCycleNum(inputfile2);
    if (cycle1 != cycle2)
    {
        cout << "[Error] SWATH cycle in two input file must be the same." << endl;

    }

    // sum of several matrix
    vector<Matrix> dpMatrix;
    // double * res = NULL;

    // MS1: i = 1
    // MS2 SWATH: i = 2 to cycle
    int threads = getProperThreads();
    for (int i = 1; i <= cycle1; ++i) {
        SimpleTimer st("calculate dot product of cycle " + to_string(i));
        string curOutputfile = outputbinaryfile + to_string(i)+".dp";


        if (!isFileExist(curOutputfile))// never skip this step
        {
            // There is memroy leak here! please pay attention to this
            vector<PeakList*> c = mr.ReadSWATHmzXMLToPeakLists(inputfile1,i,cycle1);
            vector<PeakList*> d = mr.ReadSWATHmzXMLToPeakLists(inputfile2,i,cycle1);
            cout << "[Info] Calculating dot product of swath " << i << endl;
            mr.MultiThreadDotProduct(c,d,curOutputfile,threads);

        }
        Matrix dp(curOutputfile+"");
        if (! isFileExist(curOutputfile+".txt"))
        {
            vector<PeakList*> c = mr.ReadSWATHmzXMLToPeakLists(inputfile1,i,cycle1);
            vector<PeakList*> d = mr.ReadSWATHmzXMLToPeakLists(inputfile2,i,cycle1);
            RunDTWRTMatch(c, d, curOutputfile, dp);
            dp.outputAsText(curOutputfile+".txt");
        }


        dpMatrix.push_back(dp);


    }
    long row = 0;
    long col = 0;
    vector<PeakList*> c = mr.ReadSWATHmzXMLToPeakLists(inputfile1,1,cycle1);
    vector<PeakList*> d = mr.ReadSWATHmzXMLToPeakLists(inputfile2,1,cycle1);
    row = c.size();
    col = d.size();

    Matrix sum_of_dpMatrix(row, col);
    sum_of_dpMatrix.AddMatrix(dpMatrix);

    RunDTWRTMatch(c, d, outputbinaryfile, sum_of_dpMatrix);

    sum_of_dpMatrix.outputBinary(outputbinaryfile+"_sum.dp");
    sum_of_dpMatrix.outputAsText(outputbinaryfile+"_sum.dp.txt");
}

void RunDTWRTMatch(const vector<PeakList *> &c, const vector<PeakList *> &d, string &outputbinaryfile,
                   Matrix &dpMatrix) {
    cout << "[Info] Run DTW ..." << endl;
    vector<vector<long> > Matches = DTW(dpMatrix);
    vector<vector<double> > RTMatches;
    getRTMatches(Matches, c, d, RTMatches);
    string RTMathesFile = outputbinaryfile + ".rt.matches";
    cout << "output sum dp" << RTMathesFile << endl;
    outputRTMatches(RTMatches, RTMathesFile);
}

void outputRTMatches(const vector<vector<double>> &RTMatches, const string &RTMathesFile) {
    FILE *pfile = fopen(RTMathesFile.c_str(),"w");
    for (int i = 0; i < RTMatches.size(); ++i) {
        fprintf(pfile,"%.5lf\t%.5lf\n", RTMatches[i][0],RTMatches[i][1]);
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
    vector<PeakList*> a = mr.ReadmzXMLToPeakLists(inputfile1);
    vector<PeakList*> b = mr.ReadmzXMLToPeakLists(inputfile2);

    AnalysisMassDiff amd;
    vector<double> massdiff = amd.CalculateMassDifference(a,0.5);
    int threads = 24;
    mr.MultiThreadDotProduct(a,b,outputbinaryfile,threads);

}

