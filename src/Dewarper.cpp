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

LowessDewarper::LowessDewarper(FeatureMap &sample, FeatureMap &reference) : DummyDewarper(sample, reference) {

    this->generateM_filename(sample_file, reference_mzXML);


}

LowessDewarper::~LowessDewarper() {

}


DummyDewarper::~DummyDewarper() {

}

void DummyDewarper::Run() {
    this->alignScans();

    this->GetRTMatches();
    cout << "[Info] Calibrate RT..." << endl;
    string predictfile = getM_filename() + ".predicted";
    ofstream fout;
    fout.open(predictfile.c_str(), ios::out);
    cout << "[Info] Output RT -- RTshift to file " << predictfile << endl;
    int number_of_nodes = rt_x.size();
    if (number_of_nodes == 0) return;
    double eps = 1e-8;
    for (int i = 0; i < m_sample.m_featureNum; i++) {
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
        double current_rt_shift = current_feature_rt > line_y ? line_y : current_feature_rt;
        fout << current_feature_rt << " " << current_rt_shift << endl;

    }
    fout.close();
    cout << "[Info] Output done" << endl;
}


DTWDewarper::DTWDewarper(FeatureMap &sample, FeatureMap &reference) : DummyDewarper(sample, reference) {
//    cout << "construct DTW dewarper" << endl;
//
    cmdline = pParam->getLwbmatch_bin_path() + "mzXMLReader  " + sample_mzXML + " " + reference_mzXML + " --ms1";

    string similarity = pParam->getSimilarity();
    double pnorm = pParam->getPnorm();
    string mergeMethod = pParam->getMergeMethod();
    string DTWOptions =
            " --similarity " + similarity + " --pnorm " + to_string(pnorm) + " --mergeMethod " + mergeMethod +
            " --filters " + pParam->getFilter() + " --ms2filters " + pParam->getMs2filter();;
    cmdline = cmdline + DTWOptions;
//    alignScans();

    this->generateM_filename(sample_file, reference_mzXML);

}


void DTWDewarper::GetRTMatches() {
    ifstream fin;
    double x, y;
    fin.open(getM_filename().c_str(), ios::in);
    if (!fin.good()) {
        cout << "[Info] can not open RT match file: " << getM_filename() << endl;
        throw "can not open file";
    }
    cout << "[Info] Reading RT->RTshifts from file " << getM_filename() << endl;
    while (fin >> x >> y) {
        RTNodes tmp;
        tmp.x = x;
        tmp.y = y;
        Nodes.push_back(tmp);
        //cout <<"method2 " << x << " " << x-y << endl;
    }
    fin.close();
    stable_sort(Nodes.begin(), Nodes.end(), RTNodes::comp);

    for (int i = 0; i < Nodes.size(); ++i) {
        rt_x.push_back(Nodes[i].x);
        rt_y.push_back(Nodes[i].y);
    }

    for (int j = 0; j < rt_x.size(); ++j) {
        delta_x_y.push_back(rt_x[j] - rt_y[j]);
    }

}


void LowessDewarper::GetRTMatches() {
//    cout << "********************************" << endl;
    LOWESSData ld = alignmentLC(&m_sample, &m_reference);
    rt_x = ld.getM_lowess_x();
    delta_x_y = ld.getM_lowess_ys();

    // output lowess results
    cout << "[Info] Output lowess results into " << getM_filename() << endl;
    ofstream fout;
    fout.open(getM_filename().c_str(), ios::out);
    for (int i = 0; i < rt_x.size(); i++) {
        fout << rt_x[i] << " " << delta_x_y[i] << endl;
    }
    fout.close();
    cout << "[Info] Output complete" << endl;
}

void DummyDewarper::GetRTMatches() {

}

LOWESSData LowessDewarper::alignmentLC(FeatureMap *reference, FeatureMap *sample) {
    //LOWESS_node Nodes[MAXN];
    vector<LOWESS_node> Nodes;

    cout << "[Info] LC alignment" << endl;

    int length;
    parameters *pParam = parameters::GetParam();
    double th_mz = pParam->getTh_MZ();
    //calc potential correct deviations
    searchForPotentialRTShifts(reference, sample, length, 2500, Nodes, th_mz);

    //apply lowess model
    stable_sort(Nodes.begin(), Nodes.begin() + length, LOWESS_comp);

    int i, RT_binSize = 50;

    double f = 0.66667, nsteps = 3, delta = (Nodes[length - 1].x - Nodes[0].x) * 0.01;
    cout << "[Info] Parameters for LOWESS: (f, nsteps, delta) = (" << f << ", " << nsteps << ", " << delta << ")"
         << endl;
    //node_num = length;
    cout << "[Info] Loading points (x, y) ..." << endl;
    LOWESSData ld;

    GetLowessXY(i, length, RT_binSize, ld, Nodes);

    cout << "[Info] Start LOWESS..." << endl;
    LOWESS *lowessModel = new LOWESS();

    lowessModel->lowess(length, f, nsteps, delta, ld);

    delete lowessModel;


    cout << "[Info] Finish LOWESS..." << endl;
    return ld;


}

int LowessDewarper::GetMaxValIndex(vector<int> V) {
    // ToDo: Why initialized as v.size/2
    int maxValIndex = V.size() / 2;
    for (int i = 0; i < V.size(); i++) {
        if (V[maxValIndex] < V[i]) { maxValIndex = i; }
    }
    return maxValIndex;
}

int LowessDewarper::GetSum(vector<int> V) {
    int sum = 0;
    for (int i = 0; i < V.size(); i++) { sum += V[i]; }
    return sum;
}

bool LowessDewarper::isInWindow(double rt, double mz, double curRT, double curMZ, double max_rt, double max_mz) {
    return fabs(rt - curRT) <= max_rt && fabs(curMZ - mz) < max_mz;
}

double LowessDewarper::searchForPotentialRTShifts(FeatureMap *reference, FeatureMap *sample, int &length, int max_rt,
                                                  vector<LOWESS_node> &Nodes, double th_MZ) {

    int i, j, q, index = 0, feature_num_counter = 0;
    const int RT_binSize = 50;// This number is too large I think
    const int arrayLength = 400, middlePoint = 200;

    //int fiveIndex[5], fiveNum[5];
    int counter[arrayLength] = {0};// non-zero entry
    reference->sortByRT();
    sample->sortByRT();

    int th_min_feature_num = floor(0.005 * reference->m_featureNum);// what does this mean?

    // Nodes are used for lowess.
    Nodes.push_back(LOWESS_node());
    cout << "search for potential rt shifts 0" << endl;
    Nodes[0].x = GetIntegerRT(reference->features[0].rt, RT_binSize);
    cout << "search for potential rt shifts 1 " << endl;
    for (i = 0; i < reference->m_featureNum; i++) {
        if (i % 101 == 0)
            cout << "search for potential rt shifts i = " << i << endl;
        // for each feauture, calculate bin num
        q = GetIntegerRT(reference->features[i].rt, RT_binSize);


        if (Nodes[index].x != q) {
            //current reference feature comes from a new bin, different from the bin we are working
            if (feature_num_counter > th_min_feature_num) { // Will start a new hist gram
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
            feature_num_counter = 0;
        }
        // in the old bin
        // else Node[index].x == q
        feature_num_counter++;
        double curMZ = reference->features[i].mz;
        double curRT = reference->features[i].rt;
        for (j = 0; j < sample->m_featureNum; j++)// histgram of 50 seconds bin around the reference rt
        {
//            if (curRT - sample->features[j].rt < )
            int RTBinNum = (curRT - sample->features[j].rt) / RT_binSize + middlePoint;
            if (RTBinNum < 0) {
                continue;
            }
            else if (RTBinNum >= arrayLength) {
                break;
            }
            else if (isInWindow(sample->features[j].rt, sample->features[j].mz, curRT, curMZ, max_rt, th_MZ)) {

                counter[RTBinNum]++; // counter contains 400 bins,working as a histgram, or bar chart
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

void LowessDewarper::GetLowessXY(int i, int length, int RT_binSize, LOWESSData &ld, vector<LOWESS_node> &Nodes) const {
    vector<double> lowess_x, lowess_y;
    lowess_x.assign(length, 0);
    lowess_y.assign(length, 0);

    for (i = 0; i < length; i++) {
        lowess_x[i] = Nodes[i].x * RT_binSize;
        lowess_y[i] = Nodes[i].y;


    }
    ld.setM_lowess_x(lowess_x);
    ld.setM_lowess_y(lowess_y);

    string lowessinputstr = getM_filename() + ".input";
    cout << "[Info] Output the input of LOWESS to " << lowessinputstr << endl;
    ofstream fout;
    fout.open(lowessinputstr.c_str(), ios::out);
    for (i = 0; i < length; i++) {
        fout << lowess_x[i] << " " << lowess_y[i] << endl;
    }
    fout.close();
    cout << "[Info] Output done" << endl;
}

void DummyDewarper::generateM_filename(const string &sample_file, const string &reference_mzXML) {
//    cout << sample_file << " " << reference_mzXML << endl;
//    cout << "DTW dummy" << endl;
    string filename = reference_mzXML.substr(0, reference_mzXML.length() - 6) + "_" +
                      sample_file.substr(0, sample_file.length() - 11) + ".dummy.rt.matches";
    setM_filename(filename);
}

void LowessDewarper::generateM_filename(const string &sample_file, const string &reference_mzXML) {
//    cout << sample_file << " " << reference_mzXML << endl;

    string filename = reference_mzXML.substr(0, reference_mzXML.length() - 6) + "_" +
                      sample_file.substr(0, sample_file.length() - 11) + ".lowess.rt.matches";
    setM_filename(filename);
    cout << "[Info] construct output file name: " << filename << endl;
}

// DTW results on only MS1
void DTWDewarper::generateM_filename(const string &sample_file, const string &reference_mzXML) {
//    cout << sample_file << " " << reference_mzXML << endl;
//    cout << "filename is" << getM_filename().c_str() << endl;
//    cout << "DTWMS1" << endl;
    string filename = reference_mzXML.substr(0, reference_mzXML.length() - 6) + "_" +
                      sample_file.substr(0, sample_file.length() - 11) + ".dat1.dp.rt.matches";
    this->setM_filename(filename);
}

// retention pairs generated from result from DTW on sum of all the dotproduct,
void DTWDewarperMSn::generateM_filename(const string &sample_file, const string &reference_mzXML) {
//    cout << "DTWMSn" << getM_filename().c_str() << endl;
//    cout << sample_file << " " << reference_mzXML << endl;
    string filename = reference_mzXML.substr(0, reference_mzXML.length() - 6) + "_" +
                      sample_file.substr(0, sample_file.length() - 11) + ".dat.rt.matches";
    setM_filename(filename);
}

DTWDewarperMSn::DTWDewarperMSn(FeatureMap &sample, FeatureMap &reference) : DTWDewarper(sample,
                                                                                        reference) {
    cmdline = pParam->getLwbmatch_bin_path() + "mzXMLReader " + sample_mzXML + " " + reference_mzXML + "  --ms1+ms2";
    string similarity = pParam->getSimilarity();
    double pnorm = pParam->getPnorm();
    string mergeMethod = pParam->getMergeMethod();
    string DTWOptions =
            " --similarity " + similarity + " --pnorm " + to_string(pnorm) + " --mergeMethod " + mergeMethod +
            " --filters " + pParam->getFilter() + " --ms2filters " + pParam->getMs2filter();
//    string DTWOptions = " --similarity " + similarity + " --pnorm " + to_string(pnorm);
    cmdline = cmdline + DTWOptions;
//    alignScans(pParam, sample_mzXML, reference_mzXML);

//    generateM_filename(sample_file, reference_mzXML);

}

//void DTWDewarperMSn::alignScans(const parameters *pParam, const string &sample_mzXML, const string &reference_mzXML) const {
//
//    cout << "[Info] >" << cmdline << endl;
//    system(cmdline.c_str());
//}

DTWDewarperMS2::DTWDewarperMS2(FeatureMap &sample, FeatureMap &reference) : DTWDewarper(sample, reference) {


//    alignScans(pParam, sample_mzXML, reference_mzXML);

    cmdline = pParam->getLwbmatch_bin_path() + "mzXMLReader " + sample_mzXML + " " + reference_mzXML + "  --ms2";
    string similarity = pParam->getSimilarity();
    double pnorm = pParam->getPnorm();
//    string DTWOptions = " --similarity " + similarity + " --pnorm " + to_string(pnorm);
    string mergeMethod = pParam->getMergeMethod();
//    string filter = pParam->getFilter();
    string DTWOptions =
            " --similarity " + similarity + " --pnorm " + to_string(pnorm) + " --mergeMethod " + mergeMethod +
            " --filters " + pParam->getFilter() + " --ms2filters " + pParam->getMs2filter();
    cmdline = cmdline + DTWOptions;
    //cout << "cmdline" << endl;
    this->generateM_filename(sample_file, reference_mzXML);


    return;
}

//void DTWDewarperMS2::alignScans() {
//    cout << "[Info] >" << cmdline << endl;
//
//    system(cmdline.c_str());
//}

void DTWDewarperMS2::generateM_filename(const string &sample_file, const string &reference_mzXML) {
    string filename = reference_mzXML.substr(0, reference_mzXML.length() - 6) + "_" +
                      sample_file.substr(0, sample_file.length() - 11) + ".dat_sum_ms2.rt.matches";
    setM_filename(filename);
//    cout << "ms2 file " << filename << endl;
}

void DummyDewarper::alignScans() {
//    cout << "should not be dummy" << endl;
    cout << "[Info] > " << cmdline << endl;
    system(cmdline.c_str());
}

DTWDewarperMS2InOne::DTWDewarperMS2InOne(FeatureMap &sample, FeatureMap &reference) : DTWDewarper(sample, reference) {
    cmdline = pParam->getLwbmatch_bin_path() + "mzXMLReader " + sample_mzXML + " " + reference_mzXML + "  --ms2in1";
    string similarity = pParam->getSimilarity();
    double pnorm = pParam->getPnorm();
//    string DTWOptions = " --similarity " + similarity + " --pnorm " + to_string(pnorm);
    string mergeMethod = pParam->getMergeMethod();
//    string filter = pParam->getFilter();
    string DTWOptions =
            " --similarity " + similarity + " --pnorm " + to_string(pnorm) + " --mergeMethod " + mergeMethod +
            " --filters " + pParam->getFilter() + " --ms2filters " + pParam->getMs2filter();
    cmdline = cmdline + DTWOptions;
    //cout << "cmdline" << endl;
    this->generateM_filename(sample_file, reference_mzXML);
}

void DTWDewarperMS2InOne::generateM_filename(const string &sample_file, const string &reference_mzXML) {
    string filename = reference_mzXML.substr(0, reference_mzXML.length() - 6) + "_" +
                      sample_file.substr(0, sample_file.length() - 11) + ".dat_ms2in1.dp.rt.matches";
    setM_filename(filename);
}

DTWDewarperMS2Rand::DTWDewarperMS2Rand(FeatureMap &sample, FeatureMap &reference) : DTWDewarper(sample, reference) {
    cmdline = pParam->getLwbmatch_bin_path() + "mzXMLReader " + sample_mzXML + " " + reference_mzXML + "  --ms2rand";
    string similarity = pParam->getSimilarity();
    double pnorm = pParam->getPnorm();
//    string DTWOptions = " --similarity " + similarity + " --pnorm " + to_string(pnorm);
    string mergeMethod = pParam->getMergeMethod();
//    string filter = pParam->getFilter();
    string DTWOptions =
            " --similarity " + similarity + " --pnorm " + to_string(pnorm) + " --mergeMethod " + mergeMethod +
            " --filters " + pParam->getFilter() + " --ms2filters " + pParam->getMs2filter();
    cmdline = cmdline + DTWOptions;
    //cout << "cmdline" << endl;
    this->generateM_filename(sample_file, reference_mzXML);
}

void DTWDewarperMS2Rand::generateM_filename(const string &sample_file, const string &reference_mzXML) {
    string filename = reference_mzXML.substr(0, reference_mzXML.length() - 6) + "_" +
                      sample_file.substr(0, sample_file.length() - 11) + ".dat_ms2rand.dp.rt.matches";
    setM_filename(filename);
}
