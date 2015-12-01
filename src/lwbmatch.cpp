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

#include "lwbmatch.h"
#include "Dewarper.h"


//debug
template<class T>
void OUT(T x, int n) {
    for (int i = 0; i < n; ++i) cout << x[i] << ' ';
    cout << endl;
}

template<class T>
void OUT(T x, int n, int m) {
    for (int i = 0; i < n; ++i) OUT(x[i], m);
    cout << endl;
}

const double eps = 1e-8;
#define  LL long long
#define  out(x) (cout << #x << " = " << x << endl)

#define e 2.71828182846
#define inf 0xfffffff
#define _clr(x) memset(x,0xff,sizeof(int)*MAXN)
bool debugPar = false;

void SparsityCalculator(short *weightMat, int length);

//debug
// Attention: do not waste the stack memory.


//KM algorithm
// Attention, m <= n should be satisfied to make sure the program will terminated
int LWBMatch::kuhn_munkras(int m, int n, int *match1, int *match2) {
    cout << "[Info] Running KM algorithm..." << endl;
    //famous kuhn_munkras algorithm using DFS to find augment path
    // comment by lwu:
    // match1 stores the one matched with ith label in l1
    // match1 stores the one matched with jth label in l2
    //int ut = -inf;
    //cout << ut << endl;
    int *s = new int[MAXN];
    int *t = new int[MAXN];
    int *l1 = new int[MAXN];
    int *l2 = new int[MAXN];
    //int  l1[MAXN], l2[MAXN];
    int p, q, ret = 0, i, j, k;
    //cout << "program started" << endl;
    // commentby lwu:
    // The length of l1 is m, the size of sample.featurenum
    // The length of l2 is n, the size of reference.featurenum
    // Set label, l1 for X, l2 for Y
    // initial label l1 is the maximal of the row
    // intial label l2 is zero

    //cout << ut << endl;
    for (i = 0; i < m; i++) {
        for (l1[i] = -inf, j = 0; j < n; j++) {
            l1[i] = GetWeight(i, j) > l1[i] ? GetWeight(i, j) : l1[i];
        }
    }


    for (i = 0; i < n; l2[i++] = 0);

    // Start matching
    _clr(match1), _clr(match2);
    for (i = 0; i < m; i++) {
        // For each free node in X
        _clr(t);
        q = 0;
        //commnetby lwu: first point is the current point
        // if i is not the free point, next.
        // if there is no free one, break.

        // Push the free node into S
        s[0] = i;


        for (p = 0; p <= q && match1[i] < 0; p++) {
            // For each node in S
            // commentby lwu: s and t are the relationship or belongship
            // k point to the current free point.
            k = s[p];
            // Start BFS search from free node i = s[0],
            // find every node in the equal Graph
            for (j = 0; j < n && match1[i] < 0; j++) {
                // For each node in Y, check whether the sum of the label of node k and node j is equal weight.
                // if the two endpoint of an edge has a weight equal to the sum of l1+l2
                // node k and node j is linked!
                // node j is not used before;
                // (k,j) is in the Equal Graph
                //
                if (l1[k] + l2[j] == GetWeight(k, j) && t[j] < 0) {
                    // If So, put node j into T, Push new node in X which link to j into S
                    q++;
                    // next free point
                    // who is linked with j?
                    // another free node?
                    s[q] = match2[j];// s[q] record the original one which related to node j.
                    t[j] = k; // node j is used now! Update t[j], tells the j notde is linked to k.
                    // if match2[j] is free, means the one linked with match2 is free then
                    // If we find a second free node, Augment the path is done.
                    // from the second free node, route back to the first free node, update match1 and match2.
                    if (s[q] < 0) // if matched2[j] is negative, unused!
                    {
                        // Excited! We find a alternative (argumentation) path, from free node to free node.
                        // Attention: Once the program get into this loop. match1[i]<0 will no longer satisfied.
                        for (p = j; p >= 0; j = p)
                            // frist round p >= 0 is satisfied.
                            //recursively find the first free node, start from the second free node.
                            // once we find the first free point, p is reset to -1. the loop will exit.
                        {//commentby lwu: j linked to k
                            // recover k from t[j], who linked with j
                            k = t[j];
                            // for the first time, match1[k] is -1, k = s[0] = i,, match1[k]=match1[i]=-1
                            // so p will be -1 in the first round.
                            // Then break out, i ++
                            // If for every point, we find a equal edge as a match, then the solution is done.
                            // If for some point, i, we could not find the equal edge, then we will go to the next label update procedure
                            p = match1[k];// originally, who is linked with k? node p
                            match2[j] = k;// now ,j is matched
                            match1[k] = j;// now, node j is linked with match1[k], but node p is free now.
                        }
                    }
                }
            }
        }
        // commentby lwu: Update the label:
        // commentby lwu: alternative tree
        // negative means free point.
        // if i-th point is free
        // Now the alternative path is used and the graph is larger!
        if (match1[i] < 0) {
            // commentby lwu: as we have i++ in the loop, here i-- will keep the same i.
            i--;
            p = inf;
            for (k = 0; k <= q; k++) {
                // commentby lwu: s[0,...,q]
                // for each node in S
                for (j = 0; j < n; j++) {
                    // for each node not in T
                    // commentby lwu: for each Node in Y.
                    // If node j belongs to T, continue
                    // we want the node in Equal Graph but not in T
                    if (t[j] >= 0) {
                        continue;
                    }
                    // commentby lwu: find the minimal price
                    // find the minimual gap between l1+l2 and weightmat
                    if (l1[s[k]] + l2[j] - GetWeight(s[k], j) < p) {
                        p = l1[s[k]] + l2[j] - GetWeight(s[k], j);
                    }
                }
            }

            // commentby lwu: add p for each column
            // for each node in T & Y, add p to the label,
            for (j = 0; j < n; j++) {
                // if the point is not in T
                // do nothing
                // else increase the label by p
                if (t[j] < 0) {
                    continue;
                }
                else {
                    l2[j] += p;
                }

            }
            //commentby lwu: minus p for each row.
            // for each node in the augmentation path, label decrease.
            for (k = 0; k <= q; k++) {
                // for each node record in S, decrease the label by p
                l1[s[k]] -= p;
            }
            //finally, the weight decrease by p.
        }
    }
    // calculate the total weight.
    for (i = 0; i < m; i++)
        ret += GetWeight(i, match1[i]);
    delete[] s;
    delete[] t;
    delete[] l1;
    delete[] l2;
    return ret;
}


int node_num;






//int FindFirstIDXLT(const double query, const vector<int> &AsdFiveNum) {
//    int k = 0;
//    for (k = 0; k < 5; k++) { if (query <= AsdFiveNum[k]) break; }
//    return k;
//    // This is the old way, without the else break, now we have break;
//    // otherwise, always return 5
//    //return 5;
//}

//void UpdateFiveNum(int *fiveNum, int *fiveIndex, double currentCounter, int currentIndex) {
//    vector<int> vfiveNum(fiveNum, fiveNum + 5);
//    int k = FindFirstIDXLT(currentCounter, vfiveNum);
//    if (k != 0) {
//        for (int kk = 0; kk < k - 1; kk++) {
//            fiveNum[kk] = fiveNum[kk + 1];
//            fiveIndex[kk] = fiveIndex[kk + 1];
//        }
//        fiveNum[k - 1] = currentCounter;
//        fiveIndex[k - 1] = currentIndex;
//    }
//}
//
//void ResetFiveNum(int *Five) {
//    for (int j = 0; j < 5; j++) { *(Five + j) = 0; }
//}

//void GetTopFiveNumIDX(int *fiveNum, int *fiveIndex, int *counter, int arrayLength) {
//    ResetFiveNum(fiveIndex);
//    ResetFiveNum(fiveNum);
//    for (int j = 0; j < arrayLength; j++) {
//        if (counter[j] == 0) continue;
//        double current = counter[j];
//        UpdateFiveNum(fiveNum, fiveIndex, current, j);
//
//    }
//}




//// This function is bad, please validate it.
//void UpdateCounter(FeatureMap *sample, int *counter, double curRT, double curMZ, double max_rt, double max_mz,
//                   int RT_binSize, int middlePoint) {
//    for (int j = 0; j < sample->featureNum; j++) {
//        if (isInWindow(sample->features[j].rt, sample->features[j].mz, curRT, curMZ, max_rt, max_mz)) {
//            int RTBinNum = (curRT - sample->features[j].rt) / RT_binSize + middlePoint;
//            counter[RTBinNum]++;
//        }
//    }
//}







//weight function
double sigmod(double dist) {
    return 1 / (1 + pow(e, 18 * dist - 6));
}

short int calcScore2(Feature *f1, Feature *f2, double th_RT, double th_MZ, double th_intensity) {
    if (fabs(f1->rt - f2->rt) > th_RT || fabs(f1->mz - f2->mz) > th_MZ ||
        fabs(f1->intensity - f2->intensity) > th_intensity)
        return 0;
    double dist = fabs(f1->rt - f2->rt) * 1.0 / th_RT + fabs(f1->mz - f2->mz) * 1.0 / th_MZ +
                  fabs(f1->intensity - f2->intensity) * 1.0 / th_intensity;

    return floor(1000 * (sigmod(dist)));
}
//weight function

void LWBMatch::createWeightMatrix(FeatureMap *reference, FeatureMap *sample) {
    cout << "[Info] Creating Weight Matrix..." << endl;
    int i, j;
    int index = 0;

    for (i = 0; i < sample->featureNum; i++) {
        for (j = 0; j < reference->featureNum; j++) {
            //                  weightMat[i][j]=calcScore(reference.features[i],sample.features[j]);
            short int w = calcScore2(&(sample->features[i]), &(reference->features[j]), th_RT, th_MZ,
                                     th_intensity);// delete now
            //cout << "(i,j,w) = " << i << "," << j << "," << w << endl;
            SetWeight(i, j, w);
            if (debugPar) cout << "Be careful " << endl;
            if (debugPar) cout << GetWeight(i, j) << " "; // delete now
        }
        if (debugPar) cout << endl;
    }
}

void LWBMatch::applyMemory(int featureNum) {
    cout << "[Info] Applying memory..." << endl;
    M->InitMatrix(featureNum);
    // commentby lwu: Why use MAXN here? can not understand
    //weightMat = 0;
    //weightMat = (short int(*)[MAXN])malloc(sizeof(short int)*(featureNum*MAXN));// delete

    //if (weightMat == 0) cout << "Apply memory unsuccessfully" << endl;
    //else cout << "Apply memory successfully" << endl;

}


struct Edge {
    int from, to;
    double weight;
};

bool Edge_comp(Edge a, Edge b) {
    return a.weight < b.weight;
}

bool LWBMatch::hierarchicalClustering(int clusteringResult[][2]) {

    bool clusterFlag[1000];
    int length, index = 0;
    vector<Edge> edges;
    Edge edge;
    edges.clear();
    vector<LOWESS_node> Nodes;
    //calc difference between each pair
    for (int i = 0; i < sampleNum; i++) {
        clusterFlag[i] = 1;
        for (int j = i + 1; j < sampleNum; j++) {
            LowessDewarper *lwp = new LowessDewarper(samples[i], samples[j]);
            double deviation = lwp->searchForPotentialRTShifts(&samples[i], &samples[j], length, 2500, Nodes, th_MZ);
            if (samples[i].featureNum > samples[j].featureNum) edge.from = i, edge.to = j;
            else edge.from = j, edge.to = i;
            edge.weight = deviation;
            edges.push_back(edge);
        }
    }

    //sort by difference from smallest to largest
    stable_sort(edges.begin(), edges.end(), Edge_comp);

    //decide alignment order
    for (int i = 0; i < edges.size(); i++) {
        if (clusterFlag[edges[i].from] && clusterFlag[edges[i].to]) {
            clusteringResult[index][0] = edges[i].from;
            clusteringResult[index++][1] = edges[i].to;
            clusterFlag[edges[i].to] = 0;
        }
    }
    return true;
    ///
}

//void dewarp(FeatureMap &sample) {
//	// sample
//    // lowess_x
//    // node_num
//    // lowess_ys
//    // eps
//    // max()
//	for (int i = 0; i < sample.featureNum; i++) {
//		int p = lower_bound(lowess_x, lowess_x + node_num, sample.features[i].rt) - lowess_x;
//		double line_k, line_b, line_y;
//		if (p == node_num) line_y = lowess_ys[node_num - 1];
//		else if (p == 0) line_y = lowess_ys[0];
//		else {
//			if (fabs(lowess_x[p] - sample.features[i].rt) < eps)
//			line_y = lowess_ys[p];
//			else {
//				line_k = (lowess_ys[p] - lowess_ys[p - 1]) / (lowess_x[p] - lowess_x[p - 1]);
//				line_b = lowess_ys[p] - line_k*lowess_x[p];
//				line_y = line_k*sample.features[i].rt + line_b;
//			}
//		}
//		//         if(p<node_num) cout<<lowess_x[p]<<":"<<sample.features[i].rt<<"->";
//		sample.features[i].rt = max(0.0, sample.features[i].rt - line_y);
//		//         cout<<sample.features[i].rt<<endl;
//
//	}
//}

int LWBMatch::makeConsensus(FeatureMap &reference, FeatureMap &sample) {
    cout << "[Info] Making consensus..." << endl;
    int matches = 0;

    //combine common features in both maps
    for (int i = 0; i < reference.featureNum; i++) {
        if (match2[i] != -1) {
            if (GetWeight(match2[i], i) != 0) {
                matches++;
                double rf = reference.features[i].rt;
                double sf = sample.features[match2[i]].rt;
//                cout << rf << " " << sf << " " << rf - sf << GetWeight(match2[i], i) << " " << matches << endl;
                reference.features[i].rt = (reference.features[i].rt * reference.features[i].consensusNum +
                                            sample.features[match2[i]].rt * sample.features[match2[i]].consensusNum) /
                                           (reference.features[i].consensusNum +
                                            sample.features[match2[i]].consensusNum);
                reference.features[i].mz = (reference.features[i].mz * reference.features[i].consensusNum +
                                            sample.features[match2[i]].mz * sample.features[match2[i]].consensusNum) /
                                           (reference.features[i].consensusNum +
                                            sample.features[match2[i]].consensusNum);
                reference.features[i].intensity =
                        (reference.features[i].intensity * reference.features[i].consensusNum +
                         sample.features[match2[i]].intensity * sample.features[match2[i]].consensusNum) /
                        (reference.features[i].consensusNum + sample.features[match2[i]].consensusNum);
                reference.features[i].consensusNum += sample.features[match2[i]].consensusNum;
                for (int j = 0; j < sample.features[match2[i]].consensusNum; j++) {
                    reference.features[i].FileList.push_back(sample.features[match2[i]].FileList[j]);
                    reference.features[i].RTList.push_back(sample.features[match2[i]].RTList[j]);
                    reference.features[i].MZList.push_back(sample.features[match2[i]].MZList[j]);
                    reference.features[i].IntList.push_back(sample.features[match2[i]].IntList[j]);
                    //cout <<j << ": " <<  reference.features[i].rt << " " << sample.features[match2[i]].rt << endl;
                }

            }
        }
    }

    //add features only in sample to reference
    for (int i = 0; i < sample.featureNum; i++) {
        if (match1[i] == -1 || (match1[i] != -1 && (GetWeight(i, match1[i]) == 0))) {
            reference.features.push_back(sample.features[i]);
            reference.featureNum++;
        }
    }
    if (debugPar) {
        for (int i = 0; i < reference.featureNum; i++) {
            for (int k = 0; k < reference.features[i].consensusNum; k++) {

                printf("%-20.3lf%-20.3lf%-20.3lf", reference.features[i].RTList[k], reference.features[i].MZList[k],
                       reference.features[i].IntList[k]);

            }
            puts("");
        }
    }
    return matches;
}

// Lwu: Read this one, There are three functions. First:
FeatureMap  LWBMatch::runLWBMatch() {
    parameters *pParam = parameters::GetParam();
    this->th_RT = pParam->getTh_RT();
    this->th_MZ = pParam->getTh_MZ();
    this->th_intensity = pParam->getTh_intensity();
    int min_x, min_y;
    int clusteringResult[100][2];//commentby lwu: not more than 200 files input
    FeatureMap sample, reference;
    //hierarchicalClustering to decide multi-alignment order
    cout << "[Info] Running Hierarchical Clustering..." << endl;
    hierarchicalClustering(clusteringResult);
    cout << "[Info] Start alignment" << endl;
    for (int t = 0; t < sampleNum - 1; t++) {
        min_x = clusteringResult[t][0];
        min_y = clusteringResult[t][1];

        //the id for the map with more features should be put in min_y
        if (samples[min_x].featureNum > samples[min_y].featureNum) swap(min_x, min_y);

        sample = samples[min_x];
        reference = samples[min_y];

        cout << "[Info] Construct Dewarper" << endl;

        DummyDewarper *dw = WarpingfactoryMethod(sample, reference);
        dw->Run();

        applyMemory(sample.featureNum);
        memset(match1,0,sizeof(match1));
        memset(match2,0,sizeof(match2));

        createWeightMatrix(&reference, &sample);
        //SparsityCalculator((short int *)weightMat,sample.featureNum*MAXN);
        //run KM algorithm and get the matching result to establish the consensus map

        int value = kuhn_munkras(sample.featureNum, reference.featureNum, match1, match2);


        int matches = makeConsensus(reference, sample);

        cout << "[Resu] FeatureMap " << min_x << " vs FeatureMap " << min_y << ":" << "matches->" << matches << endl;
        samples[min_x] = reference;
        samples[min_y] = reference;

    }
    return reference;
}

DummyDewarper *LWBMatch::WarpingfactoryMethod(FeatureMap &sample, FeatureMap &reference) const {
    parameters *pParam = parameters::GetParam();
    DummyDewarper *dw = NULL;
    if (pParam->getWarper() == "0")
    {
        dw = new DummyDewarper(sample, reference);
    }

    else if (pParam->getWarper() == "1")
    {
        dw = new LowessDewarper(sample, reference);
    }

    else if (pParam->getWarper() == "2")
    {
        dw = new DTWDewarper(sample,reference);
    }

    else if (pParam->getWarper() == "3")
    {
        dw = new DTWDewarperMSn(sample,reference);
    }

    else
    {
        dw = NULL;
        cout << "[Error] Incorrect warping function, try [0 ~ 3] " << endl;
        throw "Incorrect warping funciton";
    }
    return dw;
}

void SparsityCalculator(short int *weightMat, int length) {
    double Sparsity = 0;
    for (int i = 0; i < length; i++) {
        if (weightMat[i] == 0) {
            Sparsity += 1;
        }
    }
    cout << "Sparsity: " << Sparsity / length << endl;
}

