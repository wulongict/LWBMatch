#include "PeakList.h"
#include <cmath>
#include <iostream>
using namespace std;
BinningPeakList::BinningPeakList(const vector<double> &x, const vector<double> &y)
{
    const double BinSize = 0.5;
    const double  MaxMZ = 2001;
    int BinNum = floor(MaxMZ / BinSize) + 1;
    m_intensityList.assign(BinNum,0);
    for(int i = 0; i < x.size(); i ++)
    {
        int idx = floor(x[i] / BinSize);
        m_intensityList[idx] += y[i];
        if (idx -1 >= 0 && FLANKYBIN) m_intensityList[idx-1] += y[i]/2;
        if (idx +1 < BinNum && FLANKYBIN) m_intensityList[idx+1] += y[i]/2;
    }
}
BinningPeakList::~BinningPeakList()
{

}
void BinningPeakList::Print()
{
    for(int i = 0; i < GetBinNum(); i ++)
    {
        cout << m_intensityList[i] << " ";
    }
    cout << endl;
}
int BinningPeakList::GetBinNum()
{
    return m_intensityList.size();
}

double BinningPeakList::GetIntensity(int k)
{
    return m_intensityList.at(k);
}

double BinningPeakList::CalcDotProduct(BinningPeakList &other)
{
    //Print();
    //other.Print();
    double dot = 0;
    for(int i = 0; i < GetBinNum(); i ++)
    {
        dot += (other.GetIntensity(i)*GetIntensity(i));
    }
    return dot;
}



PeakList::PeakList()
{
    m_binningPeakList = NULL;
    RTinSeconds = 0;
}
PeakList::~PeakList()
{
   // cout << "Release memory" << endl;
    if (m_binningPeakList != NULL)
        delete m_binningPeakList;
}
void PeakList::InsertPeak(double mz, double intensity)
{
    this->m_mzList.push_back(mz);
    this->m_intensityList.push_back(intensity);

}

BinningPeakList * PeakList::CreateBinningPeaks()
{
    if (m_binningPeakList == NULL)
        m_binningPeakList = new BinningPeakList(m_mzList, m_intensityList);
    return m_binningPeakList;
}

double PeakList::CalcDotProduct(PeakList &other)
{
    BinningPeakList* x = this->CreateBinningPeaks();
    BinningPeakList* y = other.CreateBinningPeaks();
    double dot = x->CalcDotProduct(*y);
    double normx = x->CalcDotProduct(*x);
    double normy = y->CalcDotProduct(*y);
    //cout << normy << " " << normx << " " << dot << endl;
    if(fabs(normy)<EPSILON || fabs(normx) < EPSILON)
    {
        dot = 0;
    }
    else
    {
        dot = dot/sqrt(normx*normy);
    }


    //delete x;
    //delete y;
    return dot;
}

void TestPeakList()
{
    cout << "Test" << endl;
    PeakList x;
    x.InsertPeak(1,2);
    x.InsertPeak(2,4);
    x.InsertPeak(3,6);
    PeakList y;
    y.InsertPeak(2,4);
    PeakList z = y;

    cout << "should be (2.4) " << z.getM_mzList()[0] << " " << z.getM_intensityList()[0] << endl;

    cout << "DotProduct" << x.CalcDotProduct(y) << endl;
    cout << "Test Done" << endl;
}

PeakList PeakList::getPeaksWithin(double leftmz, double rightmz) {
    // This function used to get the peaks within a small window
    PeakList pl;
    //ToDo: rewrite this funtion to speedup.

    if (m_mzList.size() > 0)
    {
        std::vector<double>::iterator lower, upper;
        lower = std::lower_bound(m_mzList.begin(),m_mzList.end(),leftmz);
        upper = std::lower_bound(m_mzList.begin(),m_mzList.end(),rightmz);
        int start = lower - m_mzList.begin();
        int end = upper - m_mzList.begin();

        for (int i = start; i < end; ++i) {

                pl.InsertPeak(m_mzList[i],m_intensityList[i]);


        }

    }



    return pl;
}

PeakList::PeakList(const PeakList &pl) {
    // ToDo. write a blog for this one.
    // cout << "call" << endl; this is very wiered for me to look? why the error can not be solve ?
    //setM_mzList(pl.getM_mzList());
    //setM_intensityList(pl.getM_intensityList());
    m_mzList = pl.getM_mzList();
    m_intensityList = pl.getM_intensityList();
    if(pl.getM_binningPeakList()==NULL)
    {
        m_binningPeakList = NULL;
    }
    else
    {
        m_binningPeakList = new BinningPeakList(m_mzList,m_intensityList);
    }
    RTinSeconds = pl.getRTinSeconds();

}

void PeakList::addPeakList(PeakList &pl) {
    m_mzList.insert(m_mzList.end(),pl.getM_mzList().begin(), pl.getM_mzList().end());
}

double PeakList::getRTinSeconds() const {
    return RTinSeconds;
}

void PeakList::setRTinSeconds(double RTinSeconds) {
    PeakList::RTinSeconds = RTinSeconds;
}

void PeakList::setM_mzList(const vector<double> &m_mzList) {
    PeakList::m_mzList = m_mzList;
}

void PeakList::setM_intensityList(const vector<double> &m_intensityList) {
    PeakList::m_intensityList = m_intensityList;
}

BinningPeakList *PeakList::getM_binningPeakList() const {
    return m_binningPeakList;
}

const vector<double> PeakList::getM_mzList() const {
    return m_mzList;
}

const vector<double> PeakList::getM_intensityList() const {
    return m_intensityList;
}

void releaseVectorPeakListPtr(vector<PeakList*> &vpl) {
    for (int i = 0; i < vpl.size(); ++i) {
        if(vpl[i] != NULL )
        {
            delete vpl[i];
            vpl[i]=NULL;
        }
    }
}