//
// Created by wulong on 7/12/15.
//

#ifndef PROJECT_MZXMLREADER_H
#define PROJECT_MZXMLREADER_H

#include <vector>
#include <iostream>
#include <glob.h>
#include <ctime>
#include <zlib.h>
#include <cmath>
#include "SpectraST_cramp.hpp"
#include "PeakList.h"
#include <fstream>
#include <thread>
#include "Util.h"
typedef string mzXMLFilename ;


class mzXMLReader
{
    int MaxScanDiff;
    string SLASH; // or
    string MSLEVEL;
public:

    mzXMLReader(string mslevel);

    vector<PeakList*> ReadmzXMLToPeakLists(mzXMLFilename f);

    void CalculateDotProduct(vector<PeakList*> &a, vector<PeakList*> &b, int start, int end, int thread_ID,double * res);

    void SingleThreadDotProduct(vector<PeakList*> &a, vector<PeakList*> &b, string outputfile);

    void CreateBinningList(vector<PeakList*> &a);

    void MultiThreadDotProduct(vector<PeakList*> &a, vector<PeakList*> &b, string outputbinaryfile, int threadNum);

    void ApplyMemoryforMatrix(const long aSize, const long bSize, long long int &nSize, double *&res);

    void CalcDotProduct(vector<PeakList *> &a, const vector<PeakList *> &b, double *res, int threadNum);



};
#endif //PROJECT_MZXMLREADER_H
