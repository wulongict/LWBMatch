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
#include <random>     /* srand, rand */

typedef string mzXMLFilename;


class mzXMLReader {
    int MaxScanDiff;
    string SLASH; // or
    string MSLEVEL;
    SimlarityMetric *m_SimCalculator;
    vector<VectorFilter *> m_PeakFilters;
    vector<PeakListFilter *> m_ms2PeakListFilters;
public:

    mzXMLReader(string mslevel, int max_scan_th = 10000);

    void setPeakFilters(string filters) {
        cout << "[Info] Seting filters: " << filters << endl;
        while (filters.length() > 0) {
            int found = filters.find_first_of("+");
            if (found != string::npos) {
                string param = filters.substr(0, found);
                filters = filters.substr(found + 1);
                m_PeakFilters.push_back(VectorFilterFactoryMethod(param));
                cout << "[Info] Add filter" << param << endl;
            }
            else {
                m_PeakFilters.push_back(VectorFilterFactoryMethod(filters));
                filters = "";
            }

        }

    }

    void setMs2PeakListFilters(string ms2PeakFilters) {
        cout << "[Info] Setting filters:  " << ms2PeakFilters << endl;
        while (ms2PeakFilters.length() > 0) {
            int found = ms2PeakFilters.find_first_of("+");
            if (found != string::npos) {
                string param = ms2PeakFilters.substr(0, found);
                ms2PeakFilters = ms2PeakFilters.substr(found + 1);
                m_ms2PeakListFilters.push_back(CreatePeakListFilterFactoryMethod(param));
            }
            else {
                m_ms2PeakListFilters.push_back(CreatePeakListFilterFactoryMethod(ms2PeakFilters));
                ms2PeakFilters = "";
            }
        }
    }

    void setSimlarityMetric(SimlarityMetric *simMetric) {
        if (m_SimCalculator != NULL) {
            delete m_SimCalculator;

        }
        m_SimCalculator = simMetric;
    }


    virtual ~mzXMLReader();

    vector<PeakList *> ReadmzXMLToPeakLists(mzXMLFilename f);

    void
    CalculateDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, int start, int end, int thread_ID, double *res);

    void SingleThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputfile);

    void SingleThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, double *res);

    //overload with vector of vector
    void SingleThreadDotProduct(vector<vector<double> > a, vector<vector<double> > b, double *res);

    //
    void
    RandomVectorMultiplication(vector<PeakList *> &a, vector<vector<double> > &product, vector<double> random_vector,
                               int num_RT, int cycle);

    void CreateBinningList(vector<PeakList *> &a);

    void MultiThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputbinaryfile, int threadNum,
                               int cycle = 1);

    //
    void RandMultiThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputbinaryfile, int threadNum,
                                   int cycle);

    void ApplyMemoryforMatrix(const long aSize, const long bSize, long long int &nSize, double *&res);

    void CalcDotProduct(vector<PeakList *> &a, const vector<PeakList *> &b, double *res, int threadNum);


    void setMaxScanDiff(int MaxScanDiff) {
        mzXMLReader::MaxScanDiff = MaxScanDiff;
    }

    void FilterMS2PeakList(PeakList *x) {
        for (int i = 0; i < m_ms2PeakListFilters.size(); i++) {
            m_ms2PeakListFilters[i]->run(x);
        }
    }

    void exportToCSV(vector<PeakList *> &spec, string csvFileName) {
        cout << "[Info] Start exporting spectra to csv file " << csvFileName << endl;
        CreateBinningList(spec);
        ofstream fout;
        fout.open(csvFileName.c_str(), ios::out);
        for (int i = 0; i < spec.size(); ++i) {
            if (i % 101 == 0 or i == spec.size() - 1)
                cout << "\rExporting: " << i << flush;
            if (spec[i] == NULL) {
                fout << "," << endl;
                continue;
            }
            BinningPeakList *b = spec[i]->getM_binningPeakList();
            if (b == NULL) {
                cout << "[Info] Empty binning list" << endl;
                throw "[Info] Invalid binning list";
            }
            for (int k = 0; k < b->GetBinNum(); k++) {
                fout << b->GetIntensity(k) << ",";
            }
            fout << endl;

        }
        fout.close();
        cout << endl << "[Info] Finish\n" << endl;
    }
};

vector<double> normalize(vector<vector<double> > a);

#endif //PROJECT_MZXMLREADER_H
