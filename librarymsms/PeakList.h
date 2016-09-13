#ifndef PEAKLIST_H_INCLUDED
#define PEAKLIST_H_INCLUDED

#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include <algorithm>

//#include "Util.h"
using namespace std;
const bool FLANKYBIN = false;
const double EPSILON = 1e-8;

// A new class for peaks transform/Filter
//
class BinningPeakList;

class VectorFilter {
    string m_Type;
public:
    VectorFilter() { m_Type = "Dummy"; }

    void setType(string type) { m_Type = type; }

    string getType() { return m_Type; }

    virtual ~VectorFilter() {}

    virtual void run(BinningPeakList *x);
};

class PeakSquareRoot : public VectorFilter {
public:
    PeakSquareRoot() { setType("PeakSquareRoot"); }

    void run(BinningPeakList *x);
};

class NormalizedToProb : public VectorFilter {
public:
    NormalizedToProb() { setType("PeakToProb"); }

    void run(BinningPeakList *x);
};

class RemovePeakLessThanMean : public VectorFilter {
public:
    // peak intenisyt - mean
    //
    RemovePeakLessThanMean() { setType("AboveMean"); }

    void run(BinningPeakList *x);
};

class RemovePeakLessThanMedian : public VectorFilter {
public:
    RemovePeakLessThanMedian() { setType("AboveMedian"); }

    void run(BinningPeakList *x);
};

class KeepTopNPeaks : public VectorFilter {
private:
    int m_topN;
public:
    KeepTopNPeaks() {
        m_topN = 50;
        setType("TopN");
    }

    void setTopN(int topN) {
        m_topN = topN;
    }

    void run(BinningPeakList *x);
};

// FactoryMethod
VectorFilter *VectorFilterFactoryMethod(string FilterType);

class BinningPeakList {
    vector<double> m_intensityList;
    int BinSize;
    double norm;


public:
    BinningPeakList();

    BinningPeakList(const vector<double> &x, const vector<double> &y);

    ~BinningPeakList();

    void setNorm(double x);

    double CalcNorm();// error here we will stop using this one
    double CalcDotProductByNonZeros(BinningPeakList &other);

    double CalcSimilarByPnorm(BinningPeakList &other, double p);

    int GetBinNum();

    double GetIntensity(int k) const;

    vector<double> GetIntensityList() const {
        return m_intensityList;
    }

    void setIntensity(int i, double intensity) {
        m_intensityList[i] = intensity;
    }

    void setIntensitylist(vector<double> intensitylist) {
        for (int i = 0; i < intensitylist.size(); ++i) {
            cout << intensitylist[i];
            m_intensityList.push_back(intensitylist[i]);
        }
    }

    BinningPeakList *filterPeaks(VectorFilter *vf) {
        vf->run(this);
        vector<int> tmp;
//        cout << "before " << tmp.size() << " " << nonzeros.size()<< endl;
        for (int i = 0; i < nonzeros.size(); i++) {
            if (fabs(m_intensityList[nonzeros[i]]) > EPSILON)
                tmp.push_back(nonzeros[i]);
        }
        nonzeros.swap(tmp);
//        cout << "after " << tmp.size() << " " << nonzeros.size()<< endl;
//        delete vf;
        return this;
    }

    void Print();

    double GetNorm();

public:
    vector<int> nonzeros;
};

class PeakList {
private:
private:
    double RTinSeconds;
    vector<double> m_mzList;
    vector<double> m_intensityList;
    BinningPeakList *m_binningPeakList; // release the resources in the destruction function
public:
    PeakList();

    PeakList(const PeakList &pl);

    ~PeakList();

    double getRTinSeconds() const;

    void setRTinSeconds(double RTinSeconds);

    void setM_mzList(const vector<double> &m_mzList);

    void setM_intensityList(const vector<double> &m_intensityList);

    BinningPeakList *getM_binningPeakList() const;

    const vector<double> getM_mzList() const;

    const vector<double> getM_intensityList() const;

    void InsertPeak(double mz, double intensity);

    void setBinningPeakList(BinningPeakList *binlist) {
        m_binningPeakList = binlist;
    }

    PeakList getPeaksWithin(double leftmz, double rightmz);

    void addPeakList(PeakList &pl, int ms2_number = 0);

    double CalcDotProduct(PeakList &other);

//    double CalcDotProduct(PeakList &other, SimlarityMetric * simcalculator)
    BinningPeakList *CreateBinningPeaks();

    void NormalizedToSum() {
        double sum = 0;
        for (int i = 0; i < m_intensityList.size(); ++i) {
            sum += m_intensityList[i];
        }
        for (int j = 0; j < m_intensityList.size(); ++j) {
            m_intensityList[j] /= sum;
        }
    }

    void KeepTopN(int N = 50) {
//        cout << "[Alert] The BinningPeakList is not updated." << endl;
        double threshold = 0;
        vector<double> v = m_intensityList;
        sort(v.begin(), v.end(), std::greater<double>());
        if (N > v.size()) {
            // do nothing, peaks less than topN

        }
        else {
            threshold = v[N - 1];
        }
        vector<double> tmp_mz, tmp_intensity;
        for (int i = 0; i < m_mzList.size(); ++i) {
            if (m_intensityList[i] >= threshold) {
                tmp_mz.push_back(m_mzList[i]);
                tmp_intensity.push_back(m_intensityList[i]);
            }
        }
//        cout << "[Info] After filter by top N #Peaks: " << tmp_mz.size() << endl;
        m_mzList = tmp_mz;
        m_intensityList = tmp_intensity;

    }
};

void releaseVectorPeakListPtr(vector<PeakList *> &vpl);


class PeakListFilter {
    string m_FilterType;
public:
    const string &getM_FilterType() const {
        return m_FilterType;
    }


    void setM_FilterType(const string &FilterType) {
        PeakListFilter::m_FilterType = FilterType;
    }

public:
    PeakListFilter() {
        m_FilterType = "Dummy";
    }

    virtual ~PeakListFilter() {}

    virtual void run(PeakList *x) {}

};


class PeakListKeepTopN : public PeakListFilter {
private:
    int m_topN;
public:
    PeakListKeepTopN() {
        m_topN = 20;
        setM_FilterType("KeepTopN");
    }

    void setTopN(int topN) {
        m_topN = topN;
    }

    void run(PeakList *x) {
        x->KeepTopN(m_topN);
    }
};

class PeakListNormalizedByMaximal : public PeakListFilter {
public:
    PeakListNormalizedByMaximal() {
        setM_FilterType("NormalizeByMax");
    }

    void run(PeakList *x) {
        x->NormalizedToSum();
    }

};

PeakListFilter *CreatePeakListFilterFactoryMethod(string param);

#endif // PEAKLIST_H_INCLUDED
