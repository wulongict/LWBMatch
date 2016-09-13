#include "PeakList.h"
#include "Util.h"
#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>


using namespace std;

BinningPeakList::BinningPeakList() {
    m_intensityList.clear();
}

BinningPeakList::BinningPeakList(const vector<double> &x, const vector<double> &y) {
    const double BinSize = 0.5;
    const double MaxMZ = 2001;// use the old value
    const float cutoff_intensity = 25;

    int BinNum = floor(MaxMZ / BinSize) + 1;
    bool spread = true;
    m_intensityList.assign(BinNum, 0);
    for (int i = 0; i < x.size(); i++) {
        //cout << i << " ";
        int idx = floor(x[i] / BinSize);
//<<<<<<< HEAD
        m_intensityList[idx] += y[i];

        if (idx - 1 >= 0 && FLANKYBIN)
            m_intensityList[idx - 1] += y[i] / 2;

        if (idx + 1 < BinNum && FLANKYBIN)
            m_intensityList[idx + 1] += y[i] / 2;
//=======
//        if(y[i] > cutoff_intensity) {
//            double temp = y[i];
//            m_intensityList[idx] += temp;
//        }

//            if (idx - 1 >= 0 && spread)
//                m_intensityList[idx - 1] += y[i] / 5;
//
//            if (idx + 1 < BinNum && spread)
//                m_intensityList[idx + 1] += y[i] / 5;
//>>>>>>> ms2_alignment
    }
    //sharpen peak with neighbour

//    double c = 5;
//    vector<double> temp = m_intensityList;
//    for(int i = 1; i < BinNum -2 ; i++){
//        double newpeak = (-c)*(temp[i-1]+temp[i+1])+(1+2*c)*temp[i];
//        m_intensityList[i] = newpeak > 0 ? newpeak : 0;
//    }


//    norm = this->CalcNorm(); // calculate norm
    norm = -1;
    // ToDo: remove this  calcNorm

    for (int j = 0; j < BinNum; ++j) {
        if (fabs(m_intensityList[j]) > EPSILON) nonzeros.push_back(j);
    }

}

BinningPeakList::~BinningPeakList() {

}

void BinningPeakList::Print() {
    for (int i = 0; i < GetBinNum(); i++) {
        cout << m_intensityList[i] << " ";
    }
    cout << endl;
}

int BinningPeakList::GetBinNum() {
    return m_intensityList.size();
}

double BinningPeakList::GetIntensity(int k) const {
    return m_intensityList.at(k);
}

double BinningPeakList::CalcDotProductByNonZeros(BinningPeakList &other) {
    //Print();
    //other.Print();
    double dot = 0;
    vector<int> v;
    int m = 0, n = 0;
    //double meanx = 0;
    // double meany = 0;
//    for(int i = 0; i < m_intensityList.size(); i++)
    while (m < nonzeros.size() && n < other.nonzeros.size()) {
        if (nonzeros[m] == other.nonzeros[n]) {
            v.push_back(nonzeros[m]);
            m++;
            n++;
        }
        else if (nonzeros[m] < other.nonzeros[n]) m++;
        else n++;
    }

    for (int i = 0; i < v.size(); i++) {
        double x = GetIntensity(v[i]);
        double y = other.GetIntensity(v[i]);
        // Richard: Tried Harmonic Mean
        // fixed by lwu
        dot += (x * y);// /(x+y);
    }
    double normx = GetNorm();
    double normy = other.GetNorm();
    //cout << normy << " " << normx << " " << dot << endl;
    double ret = 0;
    if (fabs(normy) > EPSILON && fabs(normx) > EPSILON) {
        ret = dot / (normx * normy);
        cout << "normx, normy and dot" << normx << " " << normy << " " << dot << endl;
    }
    //cout << "dot product by nonzeros " << dot << endl;
    return ret;
}

// TODO:now just copy dotproduct one.
// this function seems never called?
double BinningPeakList::CalcSimilarByPnorm(BinningPeakList &other, double p) {
    cout << "--------------------We should not be here----------------" << endl;
    double dot = 0;
    double normx_p = this->CalcNorm();
    double normy_p = other.CalcNorm();
    if (normx_p < EPSILON || normy_p < EPSILON) {
        dot = 1;
    } else {
        for (int i = 0; i < m_intensityList.size(); i++) {
            double x = GetIntensity(i);
            double y = other.GetIntensity(i);
            if (x == 0) {
                dot += y / normy_p;
            } else if (y == 0)
                dot += fabs(other.GetIntensity(i) / normy_p - GetIntensity(i) / normx_p);
        }
    }
//<<<<<<< HEAD
    double normx = GetNorm();
    double normy = other.GetNorm();
    //cout << normy << " " << normx << " " << dot << endl;
    double ret = 0;
    if (fabs(normy) > EPSILON && fabs(normx) > EPSILON) {
        ret = dot / (normx * normy);
        //cout << "normx, normy and dot" << normx << " " << normy << " " << dot << endl;
    }


    return ret;

}

// the norm here is normalized by power of 1/p
double BinningPeakList::CalcNorm() {
    double distance = 0;
    for (int i = 0; i < GetBinNum(); i++) {
        distance += (GetIntensity(i) * GetIntensity(i));
    }
    distance = sqrt(distance);
    return distance;
}

double BinningPeakList::GetNorm() {
    if (norm == -1) return CalcNorm();
    else return norm;
}


PeakList::PeakList() {
    m_binningPeakList = NULL;
    RTinSeconds = 0;
}

PeakList::~PeakList() {
    // cout << "Release memory" << endl;
    if (m_binningPeakList != NULL)
        delete m_binningPeakList;
}

void PeakList::InsertPeak(double mz, double intensity) {
    this->m_mzList.push_back(mz);
    this->m_intensityList.push_back(intensity);

}

BinningPeakList *PeakList::CreateBinningPeaks() {
    if (m_binningPeakList == NULL)
        m_binningPeakList = new BinningPeakList(m_mzList, m_intensityList);
    return m_binningPeakList;
}

//double PeakList::CalcDotProduct(PeakList &other, SimlarityMetric * simcalculator)
//{
//    BinningPeakList* x = this->CreateBinningPeaks();
//
//    BinningPeakList* y = other.CreateBinningPeaks();
//
////    double dot = x->CalcDotProductByNonZeros(*y);
//    double dot = simcalculator->calc(*x, *y);
//    return dot;
//}

double PeakList::CalcDotProduct(PeakList &other) {
    BinningPeakList *x = this->CreateBinningPeaks();

    BinningPeakList *y = other.CreateBinningPeaks();
//<<<<<<< HEAD
    // ToDo: is it possible that we have do the normalization twice? both in PeakList and Binning list?
    // Yes, fixed.
    double dot = x->CalcDotProductByNonZeros(*y);
////=======
//
//    // dot product score
//    double dot ;
//    // p-norm score
//    // calculate norm in advance, Here we only get it.
////    double normx = x->GetNorm();
////    double normy = y->GetNorm();
//    //cout << normy << " " << normx << " " << dot << endl;
//    if(fabs(normy)<EPSILON || fabs(normx) < EPSILON)
//    {
//        dot = 0;
//    }
//    else
//    {
//        dot = x->CalcDotProductByNonZeros(*y);
////        dot = dot/sqrt(normx*normy);
//    }
//    //TODO: other similarity approach
    // norm approach

//
//    double p = 1;
//    double dot = x->CalcSimilarByPnorm(*y, p);
//>>>>>>> ms2_alignment


    return dot;
}

void TestPeakList() {
    cout << "Test" << endl;
    PeakList x;
    x.InsertPeak(1, 2);
    x.InsertPeak(2, 4);
    x.InsertPeak(3, 6);
    PeakList y;
    y.InsertPeak(2, 4);
    PeakList z = y;

    cout << "should be (2.4) " << z.getM_mzList()[0] << " " << z.getM_intensityList()[0] << endl;

    cout << "DotProduct" << x.CalcDotProduct(y) << endl;
    cout << "Test Done" << endl;
}

PeakList PeakList::getPeaksWithin(double leftmz, double rightmz) {
    // This function used to get the peaks within a small window
    PeakList pl;
    //ToDo: rewrite this funtion to speedup.

    if (m_mzList.size() > 0) {
        std::vector<double>::iterator lower, upper;
        lower = std::lower_bound(m_mzList.begin(), m_mzList.end(), leftmz);
        upper = std::lower_bound(m_mzList.begin(), m_mzList.end(), rightmz);
        int start = lower - m_mzList.begin();
        int end = upper - m_mzList.begin();

        for (int i = start; i < end; ++i) {

            pl.InsertPeak(m_mzList[i], m_intensityList[i]);


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
    if (pl.getM_binningPeakList() == NULL) {
        m_binningPeakList = NULL;
    }
    else {
        m_binningPeakList = new BinningPeakList(m_mzList, m_intensityList);
    }
    RTinSeconds = pl.getRTinSeconds();

}

///<<<<<<< HEAD
//void PeakList::addPeakList(PeakList &pl) {
////This code is not correct. Richard will update this part.
//    m_mzList.insert(m_mzList.end(),pl.getM_mzList().begin(), pl.getM_mzList().end());
//=======
void PeakList::addPeakList(PeakList &pl, int ms2_number) {
    vector<double> temp_mz = pl.getM_mzList();

    //mz shift
    double mz_shift = 2000 * ms2_number;
    if (fabs(mz_shift) > EPSILON) {
        cout << "[Info] mzhshift=" << mz_shift << endl;
        //#pragma omp parallel for  schedule(dynamic)
        for (int i = 0; i < temp_mz.size(); i++) {
            temp_mz[i] += mz_shift;
        }

    }

    vector<double> temp_intensity = pl.getM_intensityList();

    m_mzList.insert(m_mzList.end(), temp_mz.begin(), temp_mz.end());
    m_intensityList.insert(m_intensityList.end(), temp_intensity.begin(), temp_intensity.end());
//>>>>>>> ms2_alignment
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

void releaseVectorPeakListPtr(vector<PeakList *> &vpl) {
    for (int i = 0; i < vpl.size(); ++i) {
        if (vpl[i] != NULL) {
            delete vpl[i];
            vpl[i] = NULL;
        }
    }
}

VectorFilter *VectorFilterFactoryMethod(string FilterType) {
    if (FilterType == "PeakSquareRoot")
        return new PeakSquareRoot();
    else if (FilterType == "PeakToProb")
        return new NormalizedToProb();
    else if (FilterType == "AboveMean") {
        return new RemovePeakLessThanMean();
    }
    else if (FilterType == "AboveMedian") {
        return new RemovePeakLessThanMedian();
    }
    else if (FilterType.length() > 3 && FilterType.substr(0, 3) == "Top") {
        KeepTopNPeaks *ret = new KeepTopNPeaks();
        int N = atoi(FilterType.substr(3).c_str());
        cout << "[Info] Top-" << N << " filter is created" << endl;
        ret->setTopN(N);
        return ret;
    }
    else {
        cout << "[Error] Invalid Filter type: " << FilterType << endl;
        cout << "[Alert] Using base class" << endl;
        return new VectorFilter();
    }
}

PeakListFilter *CreatePeakListFilterFactoryMethod(string param) {
    if (param == "Dummy") {
        return new PeakListFilter();
    }
    else if (param == "NormalizeByMax") {
        return new PeakListNormalizedByMaximal();
    }
    else if (param.length() > 3 && param.substr(0, 3) == "Top") {
        param = param.substr(3);
        int N = atoi(param.c_str());
        PeakListKeepTopN *ret = new PeakListKeepTopN();
        ret->setTopN(N);
        return ret;
    }
    else {
        cout << "[Info] Invalid MS2 PeakList filter: " << param << endl;
        throw "[Error] Invalid PeakList filter for MS2";
    }
}

void BinningPeakList::setNorm(double x) {
    this->norm = x;
}

void PeakSquareRoot::run(BinningPeakList *x) {
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], sqrt(intensity));
    }
}

void VectorFilter::run(BinningPeakList *x) {
    //cout << "[Info] Running base class" << endl;
}

void NormalizedToProb::run(BinningPeakList *x) {
    double sum = 0;
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        sum += intensity;
    }

    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], (intensity) / sum);
    }
}

void RemovePeakLessThanMean::run(BinningPeakList *x) {
    double mean = 0;
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        mean += intensity;
    }
    mean = mean / x->GetBinNum();

    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], (intensity - mean > 0) ? intensity - mean : 0);
    }
}

void RemovePeakLessThanMedian::run(BinningPeakList *x) {
    // how to ge the median?
    if (x->nonzeros.size() < 40) {

        return;
    }

    vector<double> p;
    for (int i = 0; i < x->nonzeros.size(); ++i) {
        p.push_back(x->GetIntensity(x->nonzeros[i]));
    }
    sort(p.begin(), p.end());
    double median = 0;
    int psize = p.size();
    if (psize < 3) {
        cout << "[Info] too small peaks" << endl;
    }
    else if (psize % 2 == 0) {
        median = (p[psize / 2 - 1] + p[psize / 2]) / 2;
    }
    else {
        median = p[psize / 2];
    }
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], (intensity - median > 0) ? intensity - median : 0);
    }

}

void KeepTopNPeaks::run(BinningPeakList *x) {
//    double topN = 50;
    if (m_topN > x->nonzeros.size())
        return;
    vector<double> p;
    for (int i = 0; i < x->nonzeros.size(); ++i) {
        p.push_back(x->GetIntensity(x->nonzeros[i]));
    }
    sort(p.begin(), p.end(), std::greater<double>());
//    cout << "test" <<  p.size() <<endl;
    double threshold = p[m_topN - 1];
    for (int i = 0; i < x->nonzeros.size(); i++) {
        double intensity = x->GetIntensity(x->nonzeros[i]);
        x->setIntensity(x->nonzeros[i], (intensity > threshold) ? intensity : 0);
    }


}
