#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <map>
#include <iostream>

using namespace std;
const long long MAXN = 200005l;

class Position {
public:
    int m_i;
    int m_j;

    Position(const Position &p) {
        m_i = p.m_i;
        m_j = p.m_j;
    }

    Position(int i, int j) : m_i(i), m_j(j) {}

    //Position & UpdatePos(int i, int j){m_i = i; m_j = j;}
    int GetX() const { return this->m_i; }

    int GetY() const { return this->m_j; }

    bool operator<(const Position &other) const {
        if (other.GetX() > this->m_i) {
            return true;
        }
        if (other.GetX() < this->m_i) {
            return false;
        }
        return this->m_j < other.GetY();
    }
};

class WeightMatrix {
public:
    WeightMatrix() {}

    virtual ~WeightMatrix() {}

    virtual bool InitMatrix(int featureNum) {
        return true;
    }

    virtual short int GetVal(Position p) {}

    virtual void PushVal(Position p, int value) {}

};


class SparseMatrix : public WeightMatrix {
public:
    SparseMatrix() {
        //m_maxrow = 0;
        //m_maxcol = 0;
        //m_Matrix = map<Position, short int>();
        // cout << "sizeof(int):" << sizeof(int) << endl
        // << "sizeof(short int):" << sizeof(short int) << endl
        //  << "sizeof(short): " << sizeof(short) << endl;
    }

    bool InitMatrix() {
        return true;
    }

    short int GetVal(Position p) {
        //cout << "MatrixSize: " << m_maxrow << "-by-" << m_maxcol << endl;

        if (m_Matrix.find(p) != m_Matrix.end()) {
            return m_Matrix[p];
        }
        else {
            return 0;
        }
    }

    void PushVal(Position p, int value) {
        //if(p.GetX() > m_maxrow) m_maxrow = p.GetX();
        //if(p.GetY() > m_maxcol) m_maxcol = p.GetY();
        if (value == 0) {
            return;
        }
        else {
            m_Matrix[p] = value;
        }
    }

private:
    map<Position, short int> m_Matrix;
    //int m_maxrow;
    //int m_maxcol;
};


// for basis class WeightMatrix, the sub-class SparseMatrix and origMatrix are different for memory allocation.
class origMatrix : public WeightMatrix {
public:
    origMatrix() {}

    ~origMatrix() {}

    bool InitMatrix(int featureNum) {
        freeMemory();
        return applyMemory(featureNum);
    }

    short int GetVal(Position p) {
        return m_Matrix[p.GetX()][p.GetY()];
    }

    void PushVal(Position p, int val) {
        m_Matrix[p.GetX()][p.GetY()] = val;
    }

protected:
    void freeMemory() {
        if (m_Matrix != NULL) {
            free(m_Matrix);
            m_Matrix = NULL;
        }
    }

    bool applyMemory(int featureNum) {
        m_Matrix = NULL;
        m_Matrix = (short int (*)[MAXN]) malloc(sizeof(short int) * (featureNum * MAXN));
        if (m_Matrix == NULL) {
            cout << "[Info] Memroy allocation fails" << endl;
            return true;
        }
        else {
            cout << "[Info] Memory allocation completes" << endl;
            return false;
        }

    }

private:
    short int (*m_Matrix)[MAXN];
};


#endif
