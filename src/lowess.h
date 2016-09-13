/*
    Copyright (c) <year> <copyright holders>

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use,
    copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following
    conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef LOWESS_H
#define LOWESS_H

#include <vector>


//#include "alignmentevaluation.h"
using namespace std;

//class LWOESSData;
class LOWESSData {
    // double lowess_x[MAXN], lowess_y[MAXN], lowess_rw[MAXN], lowess_ys[MAXN], lowess_res[MAXN];
    vector<double> m_lowess_x;
    vector<double> m_lowess_y;
    vector<double> m_lowess_rw;
    vector<double> m_lowess_ys;
    vector<double> m_lowess_res;
public:
    LOWESSData() {}

    ~LOWESSData() {}

    LOWESSData(const LOWESSData &other) {
        setM_lowess_x(other.getM_lowess_x());
        setM_lowess_ys(other.getM_lowess_ys());
        setM_lowess_y(other.getM_lowess_y());
        setM_lowess_rw(other.getM_lowess_rw());
        setM_lowess_res(other.getM_lowess_res());
    }


    const vector<double> &getM_lowess_x() const {
        return m_lowess_x;
    }

    void setM_lowess_x(const vector<double> &m_lowess_x) {
        LOWESSData::m_lowess_x = m_lowess_x;
    }

    const vector<double> &getM_lowess_y() const {
        return m_lowess_y;
    }

    void setM_lowess_y(const vector<double> &m_lowess_y) {
        LOWESSData::m_lowess_y = m_lowess_y;
    }

    const vector<double> &getM_lowess_rw() const {
        return m_lowess_rw;
    }

    void setM_lowess_rw(const vector<double> &m_lowess_rw) {
        LOWESSData::m_lowess_rw = m_lowess_rw;
    }

    const vector<double> &getM_lowess_ys() const {
        return m_lowess_ys;
    }

    void setM_lowess_ys(const vector<double> &m_lowess_ys) {
        LOWESSData::m_lowess_ys = m_lowess_ys;
    }

    const vector<double> &getM_lowess_res() const {
        return m_lowess_res;
    }

    void setM_lowess_res(const vector<double> &m_lowess_res) {
        LOWESSData::m_lowess_res = m_lowess_res;
    }
};


class LOWESS {
public:
    void lowess(int n, double f, int nsteps, double delta, LOWESSData &ld);
};

class DummyLowess {
public:
    int getM_n() const {
        return m_n;
    }

    void setM_n(int m_n) {
        DummyLowess::m_n = m_n;
    }

    double getM_f() const {
        return m_f;
    }

    void setM_f(double m_f) {
        DummyLowess::m_f = m_f;
    }

    int getM_nsteps() const {
        return m_nsteps;
    }

    void setM_nsteps(int m_nsteps) {
        DummyLowess::m_nsteps = m_nsteps;
    }

    double getM_delta() const {
        return m_delta;
    }

    void setM_delta(double m_delta) {
        DummyLowess::m_delta = m_delta;
    }

private:
    int m_n;
    double m_f;
    int m_nsteps;
    double m_delta;
public:
    DummyLowess(int n, double f, int nsteps, double delta) : m_n(n), m_delta(delta), m_f(f), m_nsteps(nsteps) {}

    ~DummyLowess() {}

    virtual void
    Run(vector<double> &x, vector<double> &y, vector<double> &ys, vector<double> &rw, vector<double> &res) {}
};


#endif // LOWESS_H
