//
// Created by wulong on 7/13/15.
//
#include <algorithm>
#include "Dewarper.h"
#include "parameters.h"
#include <fstream>
using namespace std;

bool LOWESS_comp(const LOWESS_node &a, const LOWESS_node &b) {
if (a.x == b.x) return a.y < b.y;
else return a.x < b.x;
}

int GetIntegerRT(double rt, const int RT_binSize) { return floor(rt / RT_binSize); }

LowessDewarper::LowessDewarper(FeatureMap &sample, FeatureMap &reference):DummyDewarper(sample,reference){




}

LowessDewarper::~LowessDewarper() {

}



DummyDewarper::~DummyDewarper() {

}

void DummyDewarper::Run() {
    this->GetRTMatches();
    cout << "[Info] Calibrate RT..." << endl;
    int number_of_nodes = rt_x.size();
    if (number_of_nodes == 0) return ;
    double eps = 1e-8;
    for (int i = 0; i < m_sample.featureNum; i++) {
        double current_feature_rt = m_sample.features[i].rt;
        int p = lower_bound(rt_x.begin(), rt_x.begin() + number_of_nodes, current_feature_rt) - rt_x.begin();
        double line_y;
        if (p == number_of_nodes) line_y = delta_x_y[number_of_nodes - 1];
        else if (p == 0) line_y = delta_x_y[0];
        else {
            if (fabs(rt_x[p] - current_feature_rt) < eps)
                line_y = delta_x_y[p];
            else {
                double lambda = (current_feature_rt - rt_x[p - 1]) / (rt_x[p] - rt_x[p - 1]);
                line_y = lambda * delta_x_y[p] + (1 - lambda) * delta_x_y[p - 1];
            }
        }


        m_sample.features[i].rt = max(0.0, current_feature_rt - line_y);
    }
}


DTWDewarper::DTWDewarper(FeatureMap &sample, FeatureMap &reference):DummyDewarper(sample,reference) {
//
    string sample_file = sample.features[0].FileList[0];
    string reference_file = reference.features[0].FileList[0];
    cout << sample_file << reference_file << endl;
    parameters * pParam = parameters::GetParam();
    //cout << pParam->getLwbmatch_bin_path();
    string mzXMLPath = pParam->getPath(pParam->getFileListPath());

    string sample_mzXML = mzXMLPath + sample_file.substr(0, sample_file.length()-11) + ".mzXML";
    string reference_mzXML = mzXMLPath + reference_file.substr(0, reference_file.length()-11) + ".mzXML";
    string cmdline = pParam->getLwbmatch_bin_path() + "mzXMLReader "+ sample_mzXML + " " + reference_mzXML + "  1";
    cout << "[Info] >" << cmdline << endl;
    system(cmdline.c_str());

    generateM_filename(sample_file, reference_mzXML);

}

// DTW results on only MS1
void DTWDewarper::generateM_filename(const string &sample_file, const string &reference_mzXML)
{
    m_filename = reference_mzXML.substr(0, reference_mzXML.length()-6)+ "_" + sample_file.substr(0,sample_file.length()-11)  + ".dat1.dp.rt.matches";
}

void DTWDewarper::GetRTMatches() {
    ifstream fin;
    double x, y;
    fin.open(m_filename.c_str(),ios::in);
    if (!fin.good())
    {
        cout << "[Info] can not open RT match file: "  << m_filename<< endl;
        throw "can not open file" ;
    }
    while(fin>> x >> y)
    {
        RTNodes tmp;
        tmp.x = x;
        tmp.y = y;
        Nodes.push_back(tmp);
        //cout <<"method2 " << x << " " << x-y << endl;
    }
    fin.close();
    stable_sort(Nodes.begin(),Nodes.end(),RTNodes::comp);

    for (int i = 0; i < Nodes.size(); ++i) {
        rt_x.push_back(Nodes[i].x);
        rt_y.push_back(Nodes[i].y);
    }

    for (int j = 0; j < rt_x.size(); ++j) {
        delta_x_y.push_back(rt_x[j] - rt_y[j]);
    }

}

// retention pairs generated from result from DTW on sum of all the dotproduct,
void DTWDewarperMSn::generateM_filename(const string &sample_file, const string &reference_mzXML) {
    string filename = reference_mzXML.substr(0, reference_mzXML.length()-6)+ "_" + sample_file.substr(0,sample_file.length()-11)  + ".dat.rt.matches";
    setM_filename(filename);
}

void LowessDewarper::GetRTMatches() {
    LOWESSData ld = alignmentLC(&m_sample, &m_reference);
    rt_x = ld.getM_lowess_x();
    delta_x_y = ld.getM_lowess_ys();
}

void DummyDewarper::GetRTMatches() {

}
