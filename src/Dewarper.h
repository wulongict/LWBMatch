//
// Created by wulong on 7/13/15.
//

#ifndef RT_ALIGNMENT_DEWARPER_H
#define RT_ALIGNMENT_DEWARPER_H


#include "featuremap.h"
#include "lowess.h"
#include "parameters.h"

#include <cstdlib>
#include <string.h>
#include <memory.h>
#include <iostream>

using namespace std;

#define  Set(a, b)  memset(a,b,sizeof(a))


struct LOWESS_node {
    double x, y;
};

class RTNodes {
public:
    RTNodes() {
        x = 0;
        y = 0;
    }

    RTNodes(const RTNodes &other) {
        x = other.x;
        y = other.y;
    }

    static bool comp(const RTNodes &a, const RTNodes &b) {
        if (a.x == b.x) return a.y < b.y;
        else return a.x < b.x;
    }

    double x;
    double y;
};


class DummyDewarper {


protected:
    string m_filename;
    string sample_file;
    string reference_file;
    parameters *pParam;
    string mzXMLPath;
    string sample_mzXML;
    string reference_mzXML;
    string cmdline;

    FeatureMap &m_sample; // what does it mean for a reference member?
    FeatureMap &m_reference;
    vector<double> rt_x;
    vector<double> rt_y;
    vector<double> delta_x_y;
public:
    DummyDewarper(FeatureMap &sample, FeatureMap &reference) : m_sample(sample), m_reference(reference) {
        sample_file = sample.features[0].FileList[0];
        reference_file = reference.features[0].FileList[0];
//    cout << sample_file << reference_file << endl;

        pParam = parameters::GetParam();
        //cout << pParam->getLwbmatch_bin_path();
        mzXMLPath = pParam->getPath(pParam->getFileListPath());

        sample_mzXML = mzXMLPath + sample_file.substr(0, sample_file.length() - 11) + ".mzXML";
        reference_mzXML = mzXMLPath + reference_file.substr(0, reference_file.length() - 11) + ".mzXML";
        cmdline = "";

        this->generateM_filename(sample_file, reference_mzXML);
    }

    virtual ~DummyDewarper();

    virtual void Run();

    virtual void GetRTMatches();

    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML);

    void setM_filename(const string &filename) {
        m_filename = filename;
    }

    const string getM_filename() const { return m_filename; }

    void alignScans();
};

int GetIntegerRT(double rt, const int RT_binSize);//{ return floor(rt / RT_binSize); }
bool LOWESS_comp(const LOWESS_node &a, const LOWESS_node &b);

class LowessDewarper : public DummyDewarper {
public:
    LowessDewarper(FeatureMap &sample, FeatureMap &reference);

    LOWESSData alignmentLC(FeatureMap *reference, FeatureMap *sample);

    int GetMaxValIndex(vector<int> V);

    int GetSum(vector<int> V);

    void GetRTMatches();

    bool isInWindow(double rt, double mz, double curRT, double curMZ, double max_rt, double max_mz);

    double searchForPotentialRTShifts(FeatureMap *reference, FeatureMap *sample, int &length, int max_rt,
                                      vector<LOWESS_node> &Nodes, double th_MZ);

    ~LowessDewarper();

    void GetLowessXY(int i, int length, int RT_binSize, LOWESSData &ld, vector<LOWESS_node> &Nodes) const;

    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML);

};

class DTWDewarper : public DummyDewarper {
private:

    vector<RTNodes> Nodes;

public:
    // ToDo: finish this class
    // get msn and calculate the DTW matrix
    DTWDewarper(FeatureMap &sample, FeatureMap &reference);

    ~DTWDewarper() {}

    void GetRTMatches();


    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML);

};

class DTWDewarperMSn : public DTWDewarper {
public:

    DTWDewarperMSn(FeatureMap &sample, FeatureMap &reference);

    virtual ~DTWDewarperMSn() {}

    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML);

//    void alignScans(const parameters *pParam, const string &sample_mzXML, const string &reference_mzXML) const;
};

class DTWDewarperMS2 : public DTWDewarper {
public:
    DTWDewarperMS2(FeatureMap &sample, FeatureMap &reference);

    virtual ~DTWDewarperMS2() {}

    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML);

//    void alignScans();
};

class DTWDewarperMS2InOne : public DTWDewarper {
public:
    DTWDewarperMS2InOne(FeatureMap &sample, FeatureMap &reference);

    virtual ~DTWDewarperMS2InOne() {}

    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML);

//    void alignScans();
};


class DTWDewarperMS2Rand : public DTWDewarper {
public:
    DTWDewarperMS2Rand(FeatureMap &sample, FeatureMap &reference);

    virtual ~DTWDewarperMS2Rand() {}

    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML);

//    void alignScans();
};

#endif //RT_ALIGNMENT_DEWARPER_H
