#ifndef PEAKLIST_H_INCLUDED
#define PEAKLIST_H_INCLUDED
#include <vector>
using namespace std;
const bool FLANKYBIN = false;
const double EPSILON = 1e-8;

class BinningPeakList
{
    vector<double> m_intensityList;
    int BinSize;
public:
    BinningPeakList(const vector<double> &x, const vector<double> &y);
    ~BinningPeakList();
    double CalcDotProduct(BinningPeakList & other);
    int GetBinNum();
    double GetIntensity(int k);
    void Print();
};

class PeakList
{
private:
private:
    double RTinSeconds;
    vector<double> m_mzList;
    vector<double> m_intensityList;
    BinningPeakList * m_binningPeakList; // release the resources in the destruction function
public:
    PeakList();
    PeakList(const PeakList & pl);
    ~PeakList();
    double getRTinSeconds() const;
    void setRTinSeconds(double RTinSeconds);
    void setM_mzList(const vector<double> &m_mzList);
    void setM_intensityList(const vector<double> &m_intensityList);
    BinningPeakList * getM_binningPeakList() const;
    const vector<double> getM_mzList() const;
    const vector<double> getM_intensityList() const;
    void InsertPeak(double mz, double intensity);
    PeakList  getPeaksWithin(double leftmz, double rightmz);
    void addPeakList(PeakList &pl);
    double CalcDotProduct(PeakList & other);
    BinningPeakList * CreateBinningPeaks();
};

void releaseVectorPeakListPtr(vector<PeakList*> &vpl);

#endif // PEAKLIST_H_INCLUDED
