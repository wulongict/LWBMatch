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

#include "alignmentevaluation.h"
#include "lwbmatch.h"
#include "dtw.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <string.h>

using namespace std;
parameters *parameters::pParam = NULL;

const string &parameters::getSimilarity() const {
    return similarity;
}

double parameters::getPnorm() const {
    return pnorm;
}

const string &parameters::getMergeMethod() const {
    return mergeMethod;
}

const string &parameters::getFilter() const {
    return filter;
}

const string &parameters::getMs2filter() const {
    return ms2filter;
}

void parameters::setMs2filter(const string &ms2filter) {
    parameters::ms2filter = ms2filter;
}

AlignmentEvaluation::AlignmentEvaluation() {
    // ToDo: initialize the vector;
    samples = vector<FeatureMap>();
}

void AlignmentEvaluation::readFeatureMaps(const char *fileListPath) {
    freopen(fileListPath, "r", stdin);
    int i, j;

    //get file name from the absolute path
    for (i = 0; i < fileNum; i++) {
        scanf("%s", fileName[i]);
        string temp(fileName[i]), fileName;
        for (j = 0; j < temp.length(); j++) {
            if (temp[j] == '/' || temp[j] == '\\') fileName.clear();// attention windows
            else {
                fileName += temp[j];
            }
        }
        strcpy(shortFileName[i], fileName.c_str());
    }
    fclose(stdin);

    //read feature map one by one
    for (i = 0; i < fileNum; i++) {
        // Here we open the file as stdin
        // I think there minght be some problem here.
        freopen(fileName[i], "r", stdin);
        cout << "[Info] loading features from " << i << "th file: " << fileName[i] << endl;
        samples[i].readFeatureMap(shortFileName[i]);
        // The stdin is closed.
        fclose(stdin);
        // ToDo: test if we could input something.
    }
}

#ifdef _FREOPEN_FILE_
void AlignmentEvaluation::writeResults(FeatureMap *consensusMap, const char * outputPath)
{
    freopen(outputPath, "w", stdout);
    double ZERO = 0;
    printf("%d %d\n", consensusMap->m_featureNum, fileNum);
    for (int i = 0; i < fileNum; i++) {
        for (int j = 0; j < strlen(shortFileName[i]); j++) {
            if (shortFileName[i][j] == '.') {
                shortFileName[i][j] = 0;
                printf("%-60s", shortFileName[i]);
                shortFileName[i][j] = '.';
            }
        }

    }
    puts("");
    for (int i = 0; i < consensusMap->m_featureNum; i++) {
        for (int j = 0; j < fileNum; j++) {
            bool suc = 0;
            for (int k = 0; k < consensusMap->features[i].consensusNum; k++) {
                if (strcmp(consensusMap->features[i].FileList[k].c_str(), shortFileName[j]) == 0) {
                    suc = 1;
                    printf("%-20.3lf%-20.3lf%-20.3lf", consensusMap->features[i].RTList[k], consensusMap->features[i].MZList[k], consensusMap->features[i].IntList[k]);
                    break;
                }
            }
            if (!suc) {
                printf("%-20.3lf%-20.3lf%-20.3lf", ZERO, ZERO, ZERO);
            }
        }
        puts("");
    }
    //freopen("/dev/tty","w",stdout); //for Linux
    freopen("CON", "w", stdout);//for Windows
}
#else

void AlignmentEvaluation::writeResults(FeatureMap *consensusMap, const char *outputPath) {
    FILE *pfile = fopen(outputPath, "w");
    //freopen(outputPath, "w", stdout);

    double ZERO = 0;

    fprintf(pfile, "%d %d\n", consensusMap->m_featureNum, fileNum);
    //printf("%d %d\n", consensusMap->m_featureNum, fileNum);
    for (int i = 0; i < fileNum; i++) {
        for (int j = 0; j < strlen(shortFileName[i]); j++) {
            if (shortFileName[i][j] == '.') {
                shortFileName[i][j] = 0;
                fprintf(pfile, "%-60s", shortFileName[i]);
                //printf("%-60s", shortFileName[i]);
                shortFileName[i][j] = '.';
            }
        }

    }
    //fputs("",pfile);
    fprintf(pfile, "\n");
    //puts("");
    for (int i = 0; i < consensusMap->m_featureNum; i++) {
        for (int j = 0; j < fileNum; j++) {
            bool suc = 0;
            for (int k = 0; k < consensusMap->features[i].consensusNum; k++) {
                if (strcmp(consensusMap->features[i].FileList[k].c_str(), shortFileName[j]) == 0) {
                    suc = 1;
                    fprintf(pfile, "%-20.3lf%-20.3lf%-20.3lf", consensusMap->features[i].RTList[k],
                            consensusMap->features[i].MZList[k], consensusMap->features[i].IntList[k]);
                    //printf("%-20.3lf%-20.3lf%-20.3lf", consensusMap->features[i].RTList[k], consensusMap->features[i].MZList[k], consensusMap->features[i].IntList[k]);
                    break;
                }
            }
            if (!suc) {
                fprintf(pfile, "%-20.3lf%-20.3lf%-20.3lf", ZERO, ZERO, ZERO);
                //printf("%-20.3lf%-20.3lf%-20.3lf", ZERO, ZERO, ZERO);
            }
        }
        fprintf(pfile, "\n");
        //fputs("", pfile);
        //puts("");
    }
    //freopen("/dev/tty","w",stdout); //for Linux
    //freopen("CON", "w", stdout);//for Windows
    fclose(pfile);
}

#endif

bool isClose(Feature &a, Feature &b, double th_RT, double th_MZ, double th_intensity) {
    if (fabs(a.mz - b.mz) < th_MZ && fabs(a.rt - b.rt) < th_RT && fabs(a.intensity - b.intensity) < th_intensity) {
        //cout << a.mz << " " << b.mz << " " << a.rt << " " << b.rt << " " << a.intensity << " " << b.intensity << endl;// Why intensity of b is zero.
        return 1;
    }
    else
        return 0;

}

void
AlignmentEvaluation::validateResult(int consensusFeatureNum, int consensusMapFileNum, int gtFeatureNum, int gtFileNum,
                                    double th_RT, double th_MZ, double th_intensity) {
    cout << "Threshold: " << th_RT << " " << th_MZ << " " << th_intensity << endl;
    int i, j, k;
    for (i = 0; i < consensusFeatureNum; i++) {// for each consensus feature,line
        for (j = 0; j < consensusMapFileNum; j++) { // for each file
            if (fabs(consensus[i][j].mz - 0) > 1e-8) {  // if mz not empty
                int index = consensus2gt[j];   // the file column of groundtruth.
                for (k = 0; k < gtFeatureNum; k++) {
                    if (fabs(gt[k][index].mz - 0) > 1e-8) {  // if mz not empty
                        if (isClose(consensus[i][j], gt[k][index], th_RT, th_MZ,
                                    th_intensity)) { // whether close feature.
                            break;
                        }
                    }

                }
                if (k == gtFeatureNum) {// if not found! any close feature in groundtruth, set consensus as zero.
                    consensus[i][j].mz = 0;
                    consensus[i][j].rt = 0;
                }
            }
        }
    }

    for (i = 0; i < gtFeatureNum; i++) { // for each feature in groundtruth
        if (gtHasConsensus[i] == 0) continue; // if groundtruth has corresponding feature in consensus.
        for (j = 0; j < gtFileNum; j++) {  // for each gt file
            if (fabs(gt[i][j].mz - 0) > 1e-8) {
                int index = gt2consensus[j];  // get the corresponding consensus index
                for (k = 0; k < consensusFeatureNum; k++) { // for each conssus feature. from file index
                    if (fabs(consensus[k][index].mz - 0) > 1e-8) {
                        if (isClose(gt[i][j], consensus[k][index], th_RT, th_MZ,
                                    th_intensity)) { // if there is a close feature
                            break;
                        }
                    }

                }
                if (k == consensusFeatureNum) { // if not found close feature
                    gt[i][j].mz = 0;
                    gt[i][j].rt = 0;
                }
            }
        }
    }
}

struct Pair {
    Feature a, b;
} pair_temp;

void AlignmentEvaluation::calcPrecisionAndRecall(int consensusFeatureNum, int consensusMapFileNum, int gtFeatureNum,
                                                 int gtFileNum, double th_RT, double th_MZ, double th_intensity,
                                                 int separatePoint) {
    cout << "Start calculating Recall and Precision" << endl;

    int i, j, k, sum = 0, correct_sum = 0, found_sum = 0, wrong_sum = 0, size, pair_num = 0, index, fNum = 0, p, q;
    double mz, precision = 0, recall = 0;

    vector<Pair> pairs[100][100];

    validateResult(consensusFeatureNum, consensusMapFileNum, gtFeatureNum, gtFileNum, th_RT, th_MZ, th_intensity);
    // set some feature as zero feature in validateResult
    cout << "groundtruth feature num " << gtFeatureNum << endl;
    for (i = 0; i < gtFeatureNum; i++) {   // for each line in groundtruth
        if (separatePoint == 0) {  // if separatePoint is ZERO
            for (j = 0; j < gtFileNum; j++) {  // for each gt File
                //cout << "j gtHasConsensus[j]" << j <<  " " << gtHasConsensus[j]<< endl;// This is the problem , zero and continue!
                if (gtHasConsensus[j] == 0)
                    continue;  // If there is no such file in groundtruth? How could this happen?

                if (fabs(gt[i][j].mz - 0) > 1e-8) { // if line i , and filename j with non-zero
                    for (k = j + 1; k < gtFileNum; k++) { // from next filename
                        if (gtHasConsensus[k] == 0) continue; // if there is information of file k
                        if (fabs(gt[i][k].mz - 0) > 1e-8) { // If line with non-zero mz
                            pair_temp.a = gt[i][j];
                            pair_temp.b = gt[i][k];
                            pairs[j][k].push_back(pair_temp); // Get a feature pair
                            sum++;
                        }
                    }
                }
            }
        }
        else { // if separatePoint is not zero
            for (j = 0; j < separatePoint; j++) {  // for each file before separatePoint

                if (gtHasConsensus[j] == 0) {
                    cout << "gtHasConsensus[j]==0" << gtHasConsensus[j] << endl;
                    //exit(0);
                    continue;
                }// if this file is in groundtruth
                if (fabs(gt[i][j].mz - 0) > 1e-8) { // if mz is non-zero
                    for (k = separatePoint; k < gtFileNum; k++) { // for each file after separatePoint
                        if (gtHasConsensus[k] == 0) continue;
                        if (fabs(gt[i][k].mz - 0) > 1e-8) {
                            pair_temp.a = gt[i][j];
                            pair_temp.b = gt[i][k];
                            pairs[j][k].push_back(pair_temp); // pushed as a feature pair
                            sum++;
                        }
                    }
                }
            }
        }
        // loop for this feature line
    }
    int matchedpair = 0;
    cout << "Get paired features: " << sum << endl;
    memset(matched, sizeof(matched), 0);
    bool suc = 0, have = 0;
    for (i = 0; i < consensusFeatureNum; i++) { // for each row, the alignment.
        if (separatePoint == 0) {
            for (j = 0; j < consensusMapFileNum; j++) { // consensusMapFileNum, for each file j
                if (fabs(consensus[i][j].mz - 0) > 1e-8) { // if the mz is not zero
                    for (k = j + 1; k < consensusMapFileNum; k++) {// for every file after this file
                        if (fabs(consensus[i][k].mz - 0) > 1e-8) { // if the mz is also not zero
                            suc = 0;
                            have = 0;
                            int indexj = consensus2gt[j], indexk = consensus2gt[k];// corresponding gt file is indexk, and the other is indexj

                            for (p = 0; p <
                                        pairs[indexj][indexk].size(); p++) { // for all those pairs of feature file j and k,
                                // Attention:
                                // Today, I comapred my evaluation medula with jijie's contourpart. I found the difference is that Jijie only check whether the consensus is exist in groundtruth.
                                // That is good, I think if you want to calculate the precision, that is perfect to calculate.
                                // However, When you want to get the recall, the question is whether the groundtruth feature can be found in consensus, you have to check every groundtruth.
                                // Unfortunately, if there are co-eluted peptides in groundtruth, you should not break out of the check when you find the first peptide from groundtruth,
                                // otherwise, you will never check the second peptide which is coeluted with the the first one.
                                //if (floor(consensus[i][j].rt) == 3102)
                                //{
                                //	double tmz = pairs[0][1][p].a.mz;
                                //	double trt = pairs[0][1][p].a.rt;
                                //	cout << tmz << " " << trt << "ground truth " << p << endl;;
                                //	tmz = consensus[i][0].mz;
                                //	trt = consensus[i][0].rt;
                                //	cout << tmz << " " << trt << "consensus " << i << endl;
                                //	//exit(0);
                                //}
                                have = 1;
                                // if the two features are paired and close with each other.
                                if (isClose(consensus[i][j], pairs[indexj][indexk][p].a, th_RT, th_MZ, th_intensity) &&
                                    isClose(consensus[i][k], pairs[indexj][indexk][p].b, th_RT, th_MZ, th_intensity)) {
                                    matchedpair++;
                                    //cout <<"checked :" <<  indexj << " " << indexk<< " " << p << endl;
                                    if (matched[indexj][indexk][p] == 0) { //
                                        matched[indexj][indexk][p] = 1;
                                        correct_sum++;
                                    }
                                    suc = 1;
                                    break;
                                }
                            }

                            //cout << "suc have !suc&&have" << suc << " " << have << " " << (!suc&&have) << endl;
                            if (!suc && have) { // if fail but have, so we make one wrong match.
                                // 				cout<<consensus[i][j].rt<<" "<<consensus[i][j].mz<<" "<<consensus[i][k].rt<<" "<<consensus[i][k].mz<<endl;
                                wrong_sum++;
                            }
                        }
                    }
                }
            }
        }
        else {
            for (j = 0; j < separatePoint; j++) {// for every file before the separatePoint
                if (fabs(consensus[i][j].mz - 0) > 1e-8) {
                    for (k = separatePoint; k < consensusMapFileNum; k++) { // for every file after the separatePoint
                        if (fabs(consensus[i][k].mz - 0) > 1e-8) {
                            suc = 0;
                            have = 0;
                            int indexj = consensus2gt[j], indexk = consensus2gt[k];

                            for (p = 0; p < pairs[indexj][indexk].size(); p++) {
                                have = 1;
                                if (isClose(consensus[i][j], pairs[indexj][indexk][p].a, th_RT, th_MZ, th_intensity) &&
                                    isClose(consensus[i][k], pairs[indexj][indexk][p].b, th_RT, th_MZ, th_intensity)) {
                                    if (matched[indexj][indexk][p] == 0) {
                                        matched[indexj][indexk][p] = 1;
                                        correct_sum++;
                                    }
                                    suc = 1;
                                    break;
                                }
                            }
                            if (!suc && have) {
                                // 				cout<<consensus[i][j].rt<<" "<<consensus[i][j].mz<<" "<<consensus[i][k].rt<<" "<<consensus[i][k].mz<<endl;
                                wrong_sum++;
                            }
                        }
                    }
                }
            }
        }

    }

    cout << "matched pair in consensus: " << matchedpair << endl;

    cout << "The number of all pairs in the groundtruth:" << sum << endl;
    found_sum = wrong_sum + correct_sum;
    cout << "The number of all pairs found by tool:" << found_sum << endl;
    cout << "The number of correct pairs found by tool:" << correct_sum << endl;
    precision = correct_sum * 1.0 / found_sum;
    recall = correct_sum * 1.0 / sum;
    cout << "Precision:" << precision << endl;
    cout << "Recall:" << recall << endl;
}

void
AlignmentEvaluation::evaluateConsensusMap(const char *consensusMapPath, const char *gtPath, double th_RT, double th_MZ,
                                          double th_intensity, int separatePoint) {
    int consensusFeatureNum, consensusMapFileNum, gtFeatureNum, gtFileNum;
    char cShortFileName[100][500], gtShortFileName[100][500], sequence[500], temp[500];
    //read consensusMap
    cout << "Start Reading consencus Map " << endl;
    cout << "Reading file " << consensusMapPath << endl;
    freopen(consensusMapPath, "r", stdin);
    scanf("%d%d", &consensusFeatureNum, &consensusMapFileNum);

    cout << "Reading consensus file name" << endl;
    for (int i = 0; i < consensusMapFileNum; i++) scanf("%s", cShortFileName[i]);

    consensus.reserve(consensusFeatureNum);
    consensus.resize(consensusFeatureNum);

    for (int i = 0; i < consensusFeatureNum; i++) {
        consensus[i].reserve(consensusMapFileNum);
        consensus[i].resize(consensusMapFileNum);
        for (int j = 0; j < consensusMapFileNum; j++) {
            scanf("%lf%lf%lf", &consensus[i][j].rt, &consensus[i][j].mz, &consensus[i][j].intensity);
        }

    }
    fclose(stdin);
    cout << "Read consensus done" << endl;

    //
    // apply memory for gtHasConsensus
    memset(gtHasConsensus, sizeof(gtHasConsensus), 0);
    cout << "Reading Ground truth" << endl;
    freopen(gtPath, "r", stdin);
    scanf("%d%d", &gtFeatureNum, &gtFileNum);
    scanf("%s%s", temp, temp);
    for (int i = 0; i < gtFileNum; i++) {
        // Read file names
        scanf("%s", gtShortFileName[i]);
        cout << gtShortFileName[i] << "---> filename " << endl;
        for (int j = 0; j < consensusMapFileNum; j++) {
            cout << cShortFileName[j] << "--> consensus file name " << endl;
            if (strcmp(gtShortFileName[i], cShortFileName[j]) == 0) {// What is cShortFilename[j]
                gt2consensus[i] = j; // grountTruth File i  == consensus Map fil j
                consensus2gt[j] = i; // vise versa
                gtHasConsensus[i] = 1; // GroundTruth Colum i Can Be Found in Consensus.
            }
        }
    }

    //read groudtruth
    gt.reserve(gtFeatureNum);// gt is vector<vector<feature>>
    gt.resize(gtFeatureNum);
    for (int i = 0; i < gtFeatureNum; i++) { // for each line in ground truth file
        gt[i].reserve(gtFileNum);
        gt[i].resize(gtFileNum);
        scanf("%s", sequence);// read the sequence
        for (int j = 0; j < gtFileNum; j++) { // for each file in gtfile
            scanf("%s", temp);// Read a string
            if (strcmp(temp, "NA") == 0) gt[i][j].rt = 0; // If it is NA, then rt is zero
            else sscanf(temp, "%lf", &gt[i][j].rt);

            scanf("%s", temp); // If it is NA, then mz is zero
            if (strcmp(temp, "NA") == 0) gt[i][j].mz = 0;
            else sscanf(temp, "%lf", &gt[i][j].mz);
            string stemp(sequence);
            gt[i][j].peptideSequence = stemp; // store sequence
        }
    }
    fclose(stdin);
    calcPrecisionAndRecall(consensusFeatureNum, consensusMapFileNum, gtFeatureNum, gtFileNum, th_RT, th_MZ,
                           th_intensity, separatePoint);
}

void AlignmentEvaluation::runEvaluation() {
    parameters *pParam = parameters::GetParam();
    int i, j;

    // Read from featureXML files get the RT, MZ, Intensity triple
    cout << "[Info] Loading featureXML..." << endl;
    this->fileNum = pParam->getFileNum();// ->fileNum;
    samples.reserve(fileNum);
    samples.resize(fileNum);
    readFeatureMaps(pParam->getFileListPath().c_str()/*.fileListPath.c_str()*/ );


    //////////////////////////////////////////////LWBMatch
    if (pParam->getMethodID() == 1) {
        WeightMatrix *WM;
        if (pParam->getFaster() == "1") WM = new origMatrix();
        else WM = new SparseMatrix();
        LWBMatch *lwbMatch = new LWBMatch(WM);
        lwbMatch->sampleNum = fileNum;
        lwbMatch->samples = samples;

        cout << "[Info] Running LWBMatch..." << endl;
        consensusMap = lwbMatch->runLWBMatch();


        cout << "[Info] Writing consensus..." << endl;
        writeResults(&consensusMap, pParam->getOutputPath().c_str());


        delete lwbMatch;
    }


    /////////////////////////////////////////////DTW
    if (pParam->getMethodID() == 2) {
        DTW *dtw = new DTW();
        dtw->sampleNum = fileNum;
        dtw->samples = samples;

        puts("DTW starts running");
        consensusMap = dtw->runDTW(pParam->getTh_RT(), pParam->getTh_MZ(), pParam->getTh_intensity());
        puts("DTW finishs running");

        puts("start writing consensuss...");
        writeResults(&consensusMap, pParam->getOutputPath().c_str());
        puts("finish writing consensuss");

        delete dtw;
    }

    ////////////////////////////////////////////
    if (pParam->getGtPath().length() != 0) {
        puts("start evaluation...");
        evaluateConsensusMap(pParam->getOutputPath().c_str(), pParam->getGtPath().c_str(), 50, 0.5, INF,
                             pParam->getSeparatePoint());              //for Mix datasets
        //     evaluateConsensusMap(outputPath,gtPath,50,0.2,INF,separatePoint);          //for U2OS datasets
        puts("finish evaluation...");
    }
}

void AlignmentEvaluation::runEvaluation(int methodID, double th_RT, double th_MZ, double th_intensity, int fileNum,
                                        char *fileListPath, char *outputPath, char *gtPath, int separatePoint) {
    int i, j;

    // Read from featureXML files get the RT, MZ, Intensity triple
    cout << "[Info] Reading featureXML..." << endl;
    this->fileNum = fileNum;
    samples.reserve(fileNum);
    samples.resize(fileNum);
    readFeatureMaps(fileListPath);

    //////////////////////////////////////////////LWBMatch
    if (methodID == 1) {
        LWBMatch *lwbMatch = new LWBMatch(new origMatrix);
        lwbMatch->sampleNum = fileNum;
        lwbMatch->samples = samples;

        cout << "[Info] Running LWBMatch..." << endl;
        consensusMap = lwbMatch->runLWBMatch();


        cout << "[Info] Writing consensus..." << endl;
        writeResults(&consensusMap, outputPath);

        delete lwbMatch;
    }


    /////////////////////////////////////////////DTW
    if (methodID == 2) {
        DTW *dtw = new DTW();
        dtw->sampleNum = fileNum;
        dtw->samples = samples;

        puts("DTW starts running");
        consensusMap = dtw->runDTW(th_RT, th_MZ, th_intensity);
        puts("DTW finishs running");

        puts("start writing consensuss...");
        writeResults(&consensusMap, outputPath);
        puts("finish writing consensuss");

        delete dtw;
    }

    ////////////////////////////////////////////
    if (strlen(gtPath) != 0) {
        puts("start evaluation...");
        evaluateConsensusMap(outputPath, gtPath, 50, 0.5, INF, separatePoint);              //for Mix datasets
        //     evaluateConsensusMap(outputPath,gtPath,50,0.2,INF,separatePoint);          //for U2OS datasets
        puts("finish evaluation...");
    }
}
