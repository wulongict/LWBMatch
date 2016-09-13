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

#ifndef LWBMATCH_H
#define LWBMATCH_H

#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <cstdlib>
#include <string.h>
#include <memory.h>
#include "lowess.h"
#include "featuremap.h"
#include "SparseMatrix.h"
#include "parameters.h"
#include "Dewarper.h"

using namespace std;


//par
const int INF = 1000000000;


//par
//class LOWESS_node
//{
//LOWESS_node(){
//    x = 0; y = 0;
//}
//double x;
//    double y;
//
//};



class LWBMatch {
public:
    LWBMatch(WeightMatrix *Mat) { M = Mat; }

    ~LWBMatch() {
        if (M != NULL) delete M;
        if (lowessModel != NULL)delete lowessModel;
    }

    FeatureMap runLWBMatch();

    bool hierarchicalClustering(int clusteringResult[][2]);

    int kuhn_munkras(int m, int n, int *match1, int *match2);

    //LOWESSData alignmentLC(FeatureMap *reference, FeatureMap * sample);
    void applyMemory(int featureNum);

    void createWeightMatrix(FeatureMap *reference, FeatureMap *sample);

    int makeConsensus(FeatureMap &reference, FeatureMap &sample);
    //double searchForPotentialRTShifts(FeatureMap *reference, FeatureMap *sample, int &length, int max_rt,vector<LOWESS_node> &Nodes);

    short int GetWeight(int i, int j) {
        //return weightMat[i][j];
        //Position p(0,0);
        return M->GetVal(Position(i, j));

    }

    void SetWeight(int i, int j, short int w) {
        //weightMat[i][j] = w;
        //Position p(0,0);
        M->PushVal(Position(i, j), w);
    }

    int sampleNum;
    vector<FeatureMap> samples;
    double th_RT, th_MZ, th_intensity;
    //short int (*weightMat)[MAXN]; // another method is use in WeightMatrix by Polymorphism
    int match1[MAXN], match2[MAXN];
    WeightMatrix *M;

    LOWESS *lowessModel;

    //void GetLowessXY(int i, int length, int RT_binSize) const;
    //void GetLowessXY(int i, int length, int RT_binSize, LOWESSData &ld, vector<LOWESS_node> &Nodes) const;
    DummyDewarper *WarpingfactoryMethod(FeatureMap &sample, FeatureMap &reference) const;
};


#endif // LWBMATCH_H
