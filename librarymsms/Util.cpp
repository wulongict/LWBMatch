//
// Created by wulong on 10/29/15.
//

#include <algorithm>
#include <vector>
#include "Util.h"
#include <cmath>
#include <iomanip>
//#include <vector>
#include <iostream>
//#include "PeakList.h"


using namespace std;


namespace statistic{

    

    double calcmean(vector<double> v) {
        double s = 0;
        for (int i = 0; i < v.size(); ++i) {
            s+= v[i];
        }
        return s/v.size();
    }



    double calcstd(vector<double> v) {
        double avg = calcmean(v);
        double sum_of_square = 0;
        for (int i = 0; i < v.size(); ++i) {
            sum_of_square += (v[i] - avg)*(v[i]-avg);
        }
        return sqrt(sum_of_square/v.size());
    }
}


std::vector<double> AnalysisMassDiff::CalculateMassDifference(std::vector<PeakList*> &vpl, double tolerance) {
    std::cout << "[Info] Calculate mass difference with tolerance " << tolerance << "Th" << endl;
    int Count = 0;
    int MaxDataPointsNum = 10000;
    vector<double> massdiff;
    for (int i = 0; i < vpl.size()-1; ++i) {
        PeakList * currentpl = vpl[i];
        cout << i << " / " << vpl.size() << "\r"  << flush;
        if (Count > MaxDataPointsNum) break;
        for (int j = 0; j <currentpl->getM_mzList().size(); ++j) {
            double mz = currentpl->getM_mzList()[j];
            PeakList tmp = vpl[i+1]->getPeaksWithin(mz-tolerance, mz + tolerance);
            if (tmp.getM_mzList().size() > 0)
            {
                massdiff.push_back(minDiff(mz,tmp.getM_mzList()));
                Count ++;
            }


        }
    }
    cout << endl;
    double AvgDiff = statistic::calcmean(massdiff);
    double stdDiff = statistic::calcstd(massdiff);
    cout << "[Resu] Mass deviation:" << endl;
    cout << "[Resu] mean " << AvgDiff << endl;
    cout << "[Resu] std " << stdDiff << endl;
    return massdiff;
}

double AnalysisMassDiff::minDiff(double mz, vector<double> mzlist) {
    double minMassDiff = 100;
    for (int i = 0; i < mzlist.size(); ++i) {
        double massdiff = mzlist[i] - mz;
        if( fabs(minMassDiff) > fabs(massdiff))
        {
            minMassDiff = massdiff;
        }
    }
    return minMassDiff;
}


long Matrix::getM_col() const {
    return m_col;
}

long Matrix::getM_row() const {
    return m_row;
}

Matrix::Matrix(const Matrix & other) {
    m_row = other.getM_row();
    m_col = other.getM_col();
    m_Size = m_row * m_col;
    m_entries = new double[m_Size];
    if (m_entries == NULL)
    {
        cout << "[Info] Fail to create matrix" << endl;
        throw "Fail to create matrix";
    }
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            set(i,j,other.get(i,j));
        }
    }
}

Matrix::Matrix(long r, long c) {
    m_row = r;
    m_col = c;
    m_Size = r * c;
    m_entries = new double [m_Size];
    if (m_entries == NULL)
    {
        cout << "[Info] Fail to create matrix" << endl;
        throw "Fail to create matrix";
    }
    initialize();
}

Matrix::Matrix(long r, long c, const double * matrix) {
    m_row = r;
    m_col = c;
    m_Size = r * c;
    m_entries = new double [m_Size];
    if (m_entries == NULL)
    {
        cout << "[Info] Fail to create matrix" << endl;
        throw "Fail to create matrix";
    }
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            set(i,j,matrix[i*m_col+j]);
        }
    }

}

Matrix::~Matrix() {
    if (m_entries!=NULL)
        delete [] m_entries;
}

void Matrix::Print() {
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            cout << get(i,j) << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Matrix::AddMatrix(vector<Matrix> &vMatrix) {
    cout << "[Info] Calculate sum of Matrix" << endl;
    initialize();// This is just add all of them directly together, without any tune of the weight.
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            for (int k = 0; k < vMatrix.size(); ++k) {
                double val = get(i,j) + vMatrix[k].get(i,j);
                set(i,j,val);
            }

        }
    }
}

void Matrix::outputAsText(string outputfilename) {
    cout << "[Info] Export Matrix to " << outputfilename << endl;
    FILE *pfile = fopen(outputfilename.c_str(),"w");
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            fprintf(pfile, "%.4lf\t", m_entries[i*m_col+j]);
        }
        fprintf(pfile,"\n");

    }
    fclose(pfile);
}

void Matrix::outputBinary(string outputfilename) {
    cout <<"[Info] Export Matrix to " << outputfilename << endl;
    FILE *pfile = fopen(outputfilename.c_str(),"wb");
    fwrite(&m_row,sizeof(long),1,pfile);
    fwrite(&m_col, sizeof(long),1,pfile);
    fwrite(m_entries,sizeof(double),m_Size,pfile);

    fclose(pfile);
}

Matrix::Matrix(string filename) {
    FILE *pfile = fopen(filename.c_str(),"rb");
    if (pfile == NULL)
    {
        cout << "can not open file" << endl;
        throw "Fail to open file" ;
    }

    int len = fread(&m_row, sizeof(long), 1, pfile);
    len = fread(&m_col, sizeof(long), 1, pfile);

    m_Size = m_row * m_col;

    m_entries = new double[m_Size];

    len = fread(m_entries,sizeof(double),m_Size,pfile);
    fclose(pfile);
}

double Matrix::get(long i, long j) const {
    return m_entries[i*m_col +j];
}

void Matrix::set(long i, long j, double value) {
    m_entries[i*m_col +j] = value;
}

void Matrix::initialize() {
    for (int i = 0; i < m_row; ++i) {
        for (int j = 0; j < m_col; ++j) {
            set(i,j,0);
        }
    }
}



void TestMatrix() {
    try {
        Matrix m(2,2);
        m.set(1,1, 2.0);
        m.set(0,0, 1.0);
        string filename = "testMatrix.binary";
        m.outputBinary(filename);
        m.Print();
        Matrix n(filename+"x");
        n.Print();
        cout << "[True] Testing Matrix" << endl;
    }
    catch (char const *s)
    {
        cout << "[Error] " << s << endl;
        cout << "[False] Testing Matrix" << endl;
    }





}


void getfilesize(const string &filename, long &filesize) {
    FILE *pfile = fopen(filename.c_str(),"rb");
    fseek(pfile, 0,SEEK_END);
    filesize= ftell(pfile);
    rewind(pfile);
    fclose(pfile);
}



bool isFileExist(const string &curOutputfile) {
    bool ret = false;
    FILE * pfile = fopen(curOutputfile.c_str(),"r");
    if (pfile!=NULL)
    {
        fclose(pfile);
        ret = true;
    }
    return ret;
}

SimpleTimer::SimpleTimer() {
    m_start = clock();
    m_end = m_start;
    m_used = 0;
    m_taskname = "";
    cout << "[Timer] " << flush;
    cout << "Timer is started..." << endl;

}

SimpleTimer::SimpleTimer(string taskname) {
    m_start = clock();
    m_end = m_start;
    m_used = 0;
    m_taskname = taskname;
    cout << "[Timer] " << flush;
    cout << m_taskname << " Timer is started..." << endl;
}

SimpleTimer::~SimpleTimer() {
    m_end = clock();
    m_used = (m_end -m_start)*1.0/CLOCKS_PER_SEC;
    cout << "[Timer] " << flush;
    std::cout << m_taskname << " Time used: " << std::fixed << setprecision(2)<< m_used << "s." << endl;
}

int getProperThreads() {
    int threads = std::thread::hardware_concurrency() -1;
    if (threads == 0)
    {
        threads += 1;
    }
    return threads;


}