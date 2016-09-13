//
// Created by wulong on 10/29/15.
//

#ifndef PROJECT_UTIL_H
#define PROJECT_UTIL_H
//#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <thread>
//#include <iostream>
#include <algorithm>
#include "PeakList.h"
#include <chrono>
#include <cfloat>

using namespace std;

namespace statistic {
    double calcmean(vector<double> v);

    double calcstd(vector<double> v);

    double calcGeneralizedMean(vector<double> v, double p);

//    double getNonzeroMin(const vector<double> &v);
//
//    void checkAndFixVector(vector<double> &v, double &min, double &max);

    double getNonzeroMin(const vector<double> &v);

    void checkAndFixVector(vector<double> &v, double &min, double &max);


    double calcHarmonicMean(vector<double> v);


    double calcGeometricMean(vector<double> v);

    double calcArithmeticMean(vector<double> v);

    double calcQuadraticMean(vector<double> v);

    double calcCubicMean(vector<double> v);

}


class SimpleTimer {
private:
    std::chrono::time_point<std::chrono::steady_clock> m_start;
//    int m_start;
    std::chrono::time_point<std::chrono::steady_clock> m_end;
//    int m_end;
//    double m_used;
    std::chrono::duration<double> m_used;
    string m_taskname;
public:
    SimpleTimer();

    SimpleTimer(string taskname);

    ~SimpleTimer();

};

int getProperThreads();

void getfilesize(const string &filename, long &filesize);

void TestMatrix();

bool isFileExist(const string &curOutputfile);

class Matrix {
private:
    long m_row;
    long m_col;
public:
    long getM_row() const;

    long getM_col() const;

private:
    long m_Size;
    double *m_entries;
public:

    Matrix();

    Matrix(const Matrix &other);

    Matrix(long r, long c);

    Matrix(long r, long c, const double *matrix);

    ~Matrix();


    void Print();

    void initialize(double val = 0, double r = -1, double c = -1);

    void set(long i, long j, double value);

    double get(long i, long j) const;


    Matrix(string filename);


    void outputBinary(string outputfilename);

    void outputAsText(string outputfilename);

    void mergeMatrix(vector<Matrix> &vMatrix, string mergeMethod);

//<<<<<<< HEAD
    void HomonicMatrix(vector<Matrix> &vMatrix);

//=======
    void gmOfMatrix(vector<Matrix> &vMatrix);
//>>>>>>> ms2_alignment

};

class AnalysisMassDiff {
public:

    double minDiff(double mz, vector<double> mzlist);

    std::vector<double> CalculateMassDifference(std::vector<PeakList *> &vpl, double tolerance);
};


class SimlarityMetric {
public:

    SimlarityMetric() {}

    virtual ~SimlarityMetric() {}

    virtual double calcDistance(vector<double> &a, vector<double> &b) {
        cout << "Base function called" << endl;
        return 0;
    }

    virtual double calc(BinningPeakList &a, BinningPeakList &b) {
        cout << "Base function called" << endl;
        return 0;
    }

    virtual double calcNorm(BinningPeakList &a) {
        cout << "Base function called" << endl;
        return 1;
    }

    double calcDistance(PeakList &a, PeakList &b) {
//        cout << "In base function peak a and b" << endl;
        // Here x and y are safe, because a and b will release them.
        BinningPeakList *x = a.CreateBinningPeaks();

        BinningPeakList *y = b.CreateBinningPeaks();


        double dot = this->calc(*x, *y);
        return dot;
    }

};

// Remember: This is cosine of two vector, not simply dot-product.
// Tested by gtest from Google
class DotProduct : public SimlarityMetric {
public:
    DotProduct() {}

    double calcDistance(vector<double> &a, vector<double> &b) {
        double XX = 0, YY = 0, XY = 0;
        for (int i = 0; i < a.size(); i++) {
            XX += a[i] * a[i];
            XY += a[i] * b[i];
            YY += b[i] * b[i];
        }
        double sim = 0;
//        printf("XX=%lf, YY = %lf, XY=%lf\n", XX, YY, XY);
        if (fabs(XX) > EPSILON && fabs(YY) > EPSILON)
            sim = XY / sqrt(XX * YY);
//        cout << "simlarity " << sim << endl;
        return sim;

    }

    double calcNorm(BinningPeakList &a) {
        double ret = 0;
        for (int i = 0; i < a.nonzeros.size(); i++) {
            double x = a.GetIntensity(a.nonzeros[i]);
            ret += (x * x);
        }
        ret = sqrt(ret);
        a.setNorm(ret);
        return ret;
    }

    double calc(BinningPeakList &a, BinningPeakList &b) {
//        cout << "In dot product class" << endl;
//        a.Print();
//        b.Print();
        double dot = 0;
        vector<int> v;
        int m = 0, n = 0;
        while (m < a.nonzeros.size() && n < b.nonzeros.size()) {
            if (a.nonzeros[m] == b.nonzeros[n]) {
                v.push_back(a.nonzeros[m]);
                m++;
                n++;
            }
            else if (a.nonzeros[m] < b.nonzeros[n]) m++;
            else n++;
        }

        for (int i = 0; i < v.size(); i++) {

            dot += (b.GetIntensity(v[i]) * a.GetIntensity(v[i]));
        }
        double normx = a.GetNorm();
        double normy = b.GetNorm();
//        cout <<"x, y ,dot " <<  normy << " " << normx << " " << dot << endl;
        double ret = 0;
        if (fabs(normy) > EPSILON && fabs(normx) > EPSILON) {
            ret = dot / (normx * normy);
        }
//        cout << "ret = " << ret << endl;

        if (ret < 0) {
            cout << "x, y ,dot " << normy << " " << normx << " " << dot << endl;
        }


        return ret;

    }
};

//
// This is the pearson's correlation coefficient
// ToDo;  to be tested.
class PCC : public SimlarityMetric {
public:
    PCC() {}

    double calcDistance(vector<double> &a, vector<double> &b) {
        double XX = 0, YY = 0, XY = 0;
        double meana = statistic::calcmean(a);
        double meanb = statistic::calcmean(b);
        for (int i = 0; i < a.size(); i++) {
            XX += (a[i] - meana) * (a[i] - meana);
            XY += (a[i] - meana) * (b[i] - meanb);
            YY += (b[i] - meanb) * (b[i] - meanb);
        }
        double sim = 0;
//        printf("XX=%lf, YY = %lf, XY=%lf\n", XX, YY, XY);
        if (fabs(XX) > EPSILON && fabs(YY) > EPSILON)
            sim = XY / sqrt(XX * YY);
//        cout << "simlarity " << sim << endl;
        return (sim + 1) / 2;

    }

    double calcMean(BinningPeakList &a) {
        double ret = 0;
        for (int i = 0; i < a.nonzeros.size(); ++i) {
            ret += a.GetIntensity(a.nonzeros[i]);
        }
        return ret / a.GetBinNum();
    }

    double calcNorm(BinningPeakList &a) {

        double ret = 0;
        double meana = calcMean(a);
        for (int i = 0; i < a.nonzeros.size(); i++) {
            double x = a.GetIntensity(a.nonzeros[i]);
            ret += ((x - meana) * (x - meana));
        }
        ret = sqrt(ret);
        a.setNorm(ret);
        return ret;
    }

    double calc(BinningPeakList &a, BinningPeakList &b) {
//        cout << "In dot product class" << endl;
//        a.Print();
//        b.Print();
        double dot = 0;
        double meana = calcMean(a);
        double meanb = calcMean(b);

        for (int i = 0; i < a.GetBinNum(); i++) {

            dot += ((b.GetIntensity(i) - meanb) * (a.GetIntensity(i) - meana));
        }
        double normx = a.GetNorm();
        double normy = b.GetNorm();
//        cout <<"x, y ,dot " <<  normy << " " << normx << " " << dot << endl;
        double ret = 0;
        if (normy > EPSILON && normx > EPSILON) {
            ret = dot / (normx * normy);
        }
//        ret = fabs(ret);
//        cout << "ret = " << ret << endl;
        ret = (ret + 1) / 2;

        if (ret < 0) {
            cout << "x, y ,dot " << normy << " " << normx << " " << dot << endl;
        }


        return ret;

    }
};

// Todo: to be tested by Gtest.
// ToDo: new scores to be added
// 1. Shared peak counts:
// 2. K-L divergence.
// 3. Pearson Correlaiton Coefficient: r

// Todo: This is not real Minkowski distance,
// We add 1/n in the formula
// (sum(xi^p)/n)^(1/p)
// This term is removed when we normalize the right-hand side, divided by sum of p-norm(a) and pnorm(b)
class PNorm : public SimlarityMetric {
    // If each feature/peak contains some information, then combining all those informations by p-norm might lead to a very powerful discriminator.
    // Here we use Lp-distance
    // The distance is bounded by Minkowski inequality
    // pNorm(f+g) <= pNorm(f) + pNorm(g)
    // A problem is: if f and g are of very different magnitude, then, f could be the dominate factor in pNorm(f+g)
    //
private:
    double m_p;
    // valid for
    // m_p >= 1; Regular p-Norm
    // m_p = 0;  Geometric distance
    // m_p = -1;  Harmonic distance
public:
    PNorm(double p = 2) : m_p(p) {
        if ((p < 1 && p > 0) || (p < 0 && p > -1) || (p < -1)) {
            cout << "[Error] Invalid p=" << p << " for calculating p-norm. Valid range [1, +Inf) & {0,-1}" << endl;
            throw "Error in calculate p-norm! Program will exit.";
        }

        // Todo: normalize for p ={0,-1}
        if (p == -1 || p == 0) {
            cout << "[Alert] L" << p << "-Norm is not normalized to [0,1] range" << endl;
        }

    }

    // Todo: check if a.size == b.size
    // if not, exit
    double calcDistance(vector<double> &a, vector<double> &b) {
        vector<double> distance(a.size(), 0), abs_a(a.size(), 0), abs_b(b.size(), 0);
        for (int i = 0; i < a.size(); ++i) {
            distance[i] = fabs(a[i] - b[i]);
            abs_a[i] = fabs(a[i]);
            abs_b[i] = fabs(b[i]);

        }
        double ret = 0, norm_a = 0, norm_b = 0;
        if (m_p > 0) {
            ret = statistic::calcGeneralizedMean(distance, m_p);
            // fixed: I would like to introduce the new metric by m_p, however, there are some problem here
            // How could we bound the distance,
            // should we do the normalization first?
            // We bound this by Minkowski inequality.
            // from this site, http://math.stackexchange.com/questions/581257/equality-in-minkowskis-theorem
            // we got the following conclusion.
            // The equality of pNorm(f+g) == pNorm(f) + pNorm(g) will stand when f = cg almost everywhere.
            norm_a = statistic::calcGeneralizedMean(abs_a, m_p);
            norm_b = statistic::calcGeneralizedMean(abs_b, m_p);
//            cout <<"[Info] p-norm(a-b)="  << ret << " p-norm(a)=" << norm_a << " p-norm(b)=" << norm_b << endl;
            ret = ret / (norm_a + norm_b);
        }

        else if (m_p == -1) {

            double min = *std::min_element(distance.begin(), distance.end());
            double max = *std::max_element(distance.begin(), distance.end());
            statistic::checkAndFixVector(distance, min, max);
            if (max > 0) {
                ret = statistic::calcGeneralizedMean(distance, m_p);
            }

        }
        else if (m_p == 0) {
            double min = *std::min_element(distance.begin(), distance.end());
            double max = *std::max_element(distance.begin(), distance.end());
            statistic::checkAndFixVector(distance, min, max);
            if (max > 0) {
                ret = statistic::calcGeometricMean(distance);
            }
        }
        return 1 - ret;
    }

    double calcNorm(BinningPeakList &a) {
        double ret = 0;
        for (int i = 0; i < a.nonzeros.size(); ++i) {
            ret += pow(fabs(a.GetIntensity(a.nonzeros[i])), m_p);
        }
        ret /= a.GetBinNum();
        ret = pow(ret, 1 / m_p);
        a.setNorm(ret);
        return ret;
    }

    double calc(BinningPeakList &a, BinningPeakList &b) {
        // Attention, this function is not optimized by the non-zeros.
//        cout << "In pNorm" << endl;
        // Todo: to be optimized by non-zeros
        double ret = 0;

        vector<double> v;
        vector<int> nonzeros_x;
        int m = 0, n = 0;
        while (m < a.nonzeros.size() && n < b.nonzeros.size()) {
            double x = 0, y = 0;
            if (a.nonzeros[m] < b.nonzeros[n]) {
                x = a.GetIntensity(a.nonzeros[m]);
                m++;
            }
            else if (a.nonzeros[m] > b.nonzeros[n]) {
                y = b.GetIntensity(b.nonzeros[n]);
                n++;
            }
            else {
                x = a.GetIntensity(a.nonzeros[m]);
                m++;
                y = b.GetIntensity(b.nonzeros[n]);
                n++;
            }
            ret += pow(fabs(x - y), m_p);
        }
        while (m < a.nonzeros.size()) {
            ret += pow(fabs(a.GetIntensity(a.nonzeros[m])), m_p);
            m++;
        }

        while (n < b.nonzeros.size()) {
            ret += pow(fabs(b.GetIntensity(b.nonzeros[n])), m_p);
            n++;
        }

        ret /= a.GetBinNum();
        ret = pow(ret, 1 / m_p);

//        for (int i = 0; i < a.GetBinNum(); ++i) {
//            x.push_back(a.GetIntensity(i));
//            y.push_back(b.GetIntensity(i));
//        }

        double normx = a.GetNorm();
        double normy = b.GetNorm();
        if (normx == -1 || normy == -1)
            throw "Invalid normx and normy!";
        if (normy + normx > EPSILON) { ret /= (normx + normy + EPSILON); }
        // if two of them are zero
        if (ret > 1) {
            cout << "[Info] Error ret > 1" << endl;
            cout << ret << " " << normx << " " << normy << endl;
            a.Print();
            b.Print();
            throw "Error, ret > 1";
        }
        return 1 - ret;
    }
};


class SPC : public SimlarityMetric {
public:
    SPC() {}

    double calcDistance(vector<double> &a, vector<double> &b) {
        double ret = 0, sharedcount = 0;
        double acount = 0, bcount = 0;
        for (int i = 0; i < a.size(); ++i) {
            if (fabs(a[i]) > EPSILON && fabs(b[i]) > EPSILON) {
                sharedcount += 1;
            }
            if (fabs(a[i]) > EPSILON) {
                acount += 1;
            }
            if (fabs(b[i]) > EPSILON) {
                bcount += 1;
            }

        }
        ret = sharedcount / (acount + bcount - sharedcount);

        return ret;
    }

    double calc(BinningPeakList &a, BinningPeakList &b) {
        double ret = 0, sharedcount = 0;
        double acount = a.nonzeros.size(), bcount = b.nonzeros.size();

        int m = 0, n = 0;
        while (m < acount && n < bcount) {
            if (a.nonzeros[m] == b.nonzeros[n]) {
                sharedcount += 1;
                m++;
                n++;
            }
            else if (a.nonzeros[m] < b.nonzeros[n]) {
                m++;
            }
            else { n++; }
        }
        ret = sharedcount / (acount + bcount - sharedcount);
        return ret;
    }

    double calcNorm(BinningPeakList &a) {
        return a.nonzeros.size();

    }

};

class JSD : public SimlarityMetric {
public:
    JSD() {}

    double calcDistance(vector<double> &a, vector<double> &b) {
        return 1;
    }

    double calcNorm(BinningPeakList &a) {
        return 1;
    }

    double calc(BinningPeakList &a, BinningPeakList &b) {
        //cout << "JSD" << endl;
        double div = 0;
        double sum_p = 0;
        double sum_q = 0;
        int m = 0;
        int n = 0;
        vector<double> p;
        vector<double> q;
        for (int i = 0; i < a.nonzeros.size(); i++) {
            double x = a.GetIntensity(a.nonzeros[i]);
            sum_p += x;
            p.push_back(x);
        }
        for (int i = 0; i < b.nonzeros.size(); i++) {
            double x = b.GetIntensity(b.nonzeros[i]);
            sum_q += x;
            q.push_back(x);
        }
        while (m < a.nonzeros.size() && n < b.nonzeros.size()) {
            if (a.nonzeros[m] < b.nonzeros[n]) {
                double p_m = p[m] / sum_p;
                div += p_m * log(2) / 2;
                m++;
            }
            else if (a.nonzeros[m] > b.nonzeros[n]) {
                double q_n = q[n] / sum_q;
                div += q_n * log(2) / 2;
                n++;
            }
            else {
                double p_m = p[m] / sum_p;
                double q_n = q[n] / sum_q;
                double mean = (q_n + p_m) / 2;
                div += (p_m * log(p_m / mean) + q_n * log(q_n / mean)) / 2;
                m++;
                n++;
            }
        }
        while (m < a.nonzeros.size()) {
            double p_m = p[m] / sum_p;
            div += p_m * log(2) / 2;
            m++;
        }

        while (n < b.nonzeros.size()) {
            double q_n = q[n] / sum_q;
            div += q_n * log(2) / 2;
            n++;
        }
        return 1 - div;

    }
};


#endif //PROJECT_UTIL_H
