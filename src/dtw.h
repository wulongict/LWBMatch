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

#ifndef DTW_H
#define DTW_H

#include "featuremap.h"

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

using namespace std;


class DTW {
public:
    enum ScoreFuntions {
        SPC, SIGMOD
    };

    FeatureMap runDTW(double th_RT, double th_MZ, double th_intensity);

    void applyMemory(int featureNum);

    int makeConsensus(FeatureMap &reference, FeatureMap &sample);

    double DynamicTimeWarpping(FeatureMap &reference, FeatureMap &sample);

    int score_func(Feature &reference_s, Feature &query_s, ScoreFuntions sf);

    int sharePeakCount(Feature &f1, Feature &f2);

    int sigmod(Feature &f1, Feature &f2);

    int getMapIDWithHighestFeatureNumbers();

    static const long long MAXN = 200005l;
    double th_RT, th_MZ, th_intensity;
    int sampleNum;
    vector<FeatureMap> samples;
    int (*weightMat)[MAXN];
    int match1[MAXN], match2[MAXN];


};

#endif // DTW_H
