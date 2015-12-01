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
#include "PeakList.h"
using namespace std;

namespace statistic{
    double calcmean(vector<double> v);
    double calcstd(vector<double> v);
}


class SimpleTimer{
private:
    int m_start;
    int m_end;
    double m_used;
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

class Matrix{
private:
    long m_row;
    long m_col;
public:
    long getM_row() const;

    long getM_col() const;

private:
    long m_Size;
    double * m_entries;
public:

    Matrix(const Matrix & other);

    Matrix(long r, long c);

    Matrix(long r, long c, const double * matrix);

    ~Matrix();


    void Print();

    void initialize();

    void set(long i, long j, double value);

    double get(long i, long j) const;


    Matrix(string filename);


    void outputBinary(string outputfilename);

    void outputAsText(string outputfilename);

    void AddMatrix(vector<Matrix> &vMatrix);


};

class AnalysisMassDiff
{
public:

    double minDiff(double mz, vector<double> mzlist);

    std::vector<double> CalculateMassDifference(std::vector<PeakList*> &vpl, double tolerance);
};

#endif //PROJECT_UTIL_H
