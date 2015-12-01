//
// Created by wulong on 7/13/15.
//

#ifndef RT_ALIGNMENT_DEWARPER_H
#define RT_ALIGNMENT_DEWARPER_H


#include "featuremap.h"
#include "lowess.h"
#include "parameters.h"

#include <cstdlib>
#include <string.h>
#include <memory.h>
#include <iostream>
using namespace std;

#define  Set(a, b)  memset(a,b,sizeof(a))




struct LOWESS_node {
    double x, y;
};

class RTNodes{
public:
    RTNodes(){x = 0; y = 0;}
    RTNodes(const RTNodes& other)
    {
        x = other.x;
        y = other.y;
    }
    static bool comp(const RTNodes & a,const RTNodes & b)
    {
        if (a.x == b.x) return a.y < b.y;
        else return a.x < b.x;
    }
    double x;
    double y;
};



class DummyDewarper{
protected:
    FeatureMap &m_sample;
    FeatureMap &m_reference;
    vector<double> rt_x;
    vector<double> rt_y;
    vector<double> delta_x_y;
public:
    DummyDewarper(FeatureMap &sample, FeatureMap &reference):m_sample(sample),m_reference(reference){}
    virtual ~DummyDewarper();
    virtual void Run();
    virtual void GetRTMatches();

};
int GetIntegerRT(double rt, const int RT_binSize) ;//{ return floor(rt / RT_binSize); }
bool LOWESS_comp(const LOWESS_node &a, const LOWESS_node &b);

class LowessDewarper: public DummyDewarper {
public:
    LowessDewarper(FeatureMap &sample, FeatureMap &reference);
    LOWESSData alignmentLC(FeatureMap *reference, FeatureMap *sample) {
        //LOWESS_node Nodes[MAXN];
        vector<LOWESS_node> Nodes;

        cout << "[Info] LC alignment" << endl;

        int length;
        parameters * pParam = parameters::GetParam();
        double th_mz = pParam->getTh_MZ();
        //calc potential correct deviations
        searchForPotentialRTShifts(reference, sample, length, 2500, Nodes, th_mz);

        //apply lowess model
        stable_sort(Nodes.begin(), Nodes.begin() + length, LOWESS_comp);

        int i, RT_binSize = 50;

        double f = 0.66667, nsteps = 3, delta = (Nodes[length - 1].x - Nodes[0].x) * 0.01;
        //node_num = length;
        cout << "[Info] Get Lowess x and y" << endl;
        LOWESSData ld;

        GetLowessXY(i, length, RT_binSize, ld, Nodes);

        cout << "[Info] start LOWESS" << endl;
        LOWESS * lowessModel = new LOWESS();

        lowessModel->lowess(length, f, nsteps, delta, ld);

        delete lowessModel;


        cout << "[Info] finish LOWESS" << endl;
        return ld;


    }
    int GetMaxValIndex(vector<int> V) {
        int maxValIndex = V.size() / 2;
        for (int i = 0; i < V.size(); i++) {
            if (V[maxValIndex] < V[i]) { maxValIndex = i; }
        }
        return maxValIndex;
    }

    int GetSum(vector<int> V) {
        int sum = 0;
        for (int i = 0; i < V.size(); i++) { sum += V[i]; }
        return sum;
    }

    void GetRTMatches();

    bool isInWindow(double rt, double mz, double curRT, double curMZ, double max_rt, double max_mz) {
        return fabs(rt - curRT) <= max_rt && fabs(curMZ - mz) < max_mz;
    }
    // Question: How to understand this function
    double searchForPotentialRTShifts(FeatureMap *reference, FeatureMap *sample, int &length, int max_rt,
                                  vector<LOWESS_node> &Nodes, double th_MZ) {
        int i, j, q, index = 0, sameScan = 0;
        const int RT_binSize = 50;// This number is too large I think
        const int arrayLength = 400, middlePoint = 200;

        //int fiveIndex[5], fiveNum[5];
        int counter[arrayLength] = {0};// non-zero entry
        reference->sortByRT();

        int th_sameScan = floor(0.005 * reference->featureNum);// what does this mean?

        // Nodes are used for lowess.
        Nodes.push_back(LOWESS_node());
        Nodes[0].x = GetIntegerRT(reference->features[0].rt, RT_binSize);
        for (i = 0; i < reference->featureNum; i++) {
            q = GetIntegerRT(reference->features[i].rt, RT_binSize);
            //cout << reference->features[0].FileList[0] << endl;

            if (Nodes[index].x != q) {
                if (sameScan > th_sameScan) { // Will start a new hist gram
                    int max_index = GetMaxValIndex(vector<int>(counter, counter + arrayLength));
                    int sum = GetSum(vector<int>(counter, counter + arrayLength));
                    //GetTopFiveNumIDX(fiveNum, fiveIndex, counter, arrayLength);

                    Nodes[index].y = (max_index - middlePoint) * RT_binSize; // rough shift of RT
                    //cout << "Nodes" <<Nodes[index].y << endl;

                    if (sum >= 10) {
                        index++;
                        Nodes.push_back(LOWESS_node());
                    }
                }
                Set(counter, 0);
                Nodes[index].x = q;
                sameScan = 0;
            }
            // else Node[index].x == q
            sameScan++;
            double curMZ = reference->features[i].mz;
            double curRT = reference->features[i].rt;
            for (j = 0; j < sample->featureNum; j++)// histgram of 50 seconds bin around the reference rt
            {
                if (isInWindow(sample->features[j].rt, sample->features[j].mz, curRT, curMZ, max_rt, th_MZ)) {
                    int RTBinNum = (curRT - sample->features[j].rt) / RT_binSize + middlePoint;
                    counter[RTBinNum]++;
                }
            }
        }
        length = index;


        //calc sum of deviations
        double deviation_sum = 0;
        for (i = 0; i < index; i++) {
            deviation_sum += fabs(Nodes[i].y);
        }

        if (index != 0) return deviation_sum / index; //commentby lwu: return the summation of dRT/delta_RT
        else return 0;
    }
    ~LowessDewarper();

    void GetLowessXY(int i, int length, int RT_binSize, LOWESSData &ld, vector<LOWESS_node> &Nodes) const {
        vector<double> lowess_x, lowess_y;
        lowess_x.assign(length, 0);
        lowess_y.assign(length, 0);
        for (i = 0; i < length; i++) {
            lowess_x[i] = Nodes[i].x * RT_binSize;
            lowess_y[i] = Nodes[i].y;
            cout << lowess_x[i] << " " << lowess_y[i] << endl;
        }
        ld.setM_lowess_x(lowess_x);
        ld.setM_lowess_y(lowess_y);
    }


};

class DTWDewarper: public DummyDewarper
{
private:

    vector<RTNodes> Nodes;
public:
    void setM_filename(const string &m_filename) {
        DTWDewarper::m_filename = m_filename;
    }

private:
    string m_filename;
public:
    // ToDo: finish this class
    // get msn and calculate the DTW matrix
    DTWDewarper(FeatureMap &sample, FeatureMap &reference);

    ~DTWDewarper(){}
    void GetRTMatches();


    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML);
};

class DTWDewarperMSn: public DTWDewarper
{
public:

    DTWDewarperMSn(FeatureMap &sample, FeatureMap &reference) : DTWDewarper(sample,
                                                                            reference) { }

    virtual ~DTWDewarperMSn() { }

    virtual void generateM_filename(const string &sample_file, const string &reference_mzXML) ;
};


#endif //RT_ALIGNMENT_DEWARPER_H
