#include "alignmentevaluation.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include "parameters.h"
#include "../librarymsms/Util.h"

using namespace std;

/*! \mainpage LWBmatch
 LWBMatch is a retention time (RT) alignment tool. It use a graph-based approach to align the features from several replicates. It first use LOWESS algorithm to calculate the warping function. This function is subtracted from the original RT curve to calibrate the RT shift on large extent. Then it sets up a bipartite map for a pair of features and find the maximal weighted match. The weight of the edge from feature i to feature j is calculate by some sigmoid function of linear combination of the absolute deviation of m/z, intention and rt between i and j.
 */

void printUsage() {
    cout << "Build on " << __TIME__ << " " << __DATE__ << endl;
    cout << endl;
    cout << "Usage 1:" << endl;
    cout << "lwbmatch [options] " << endl;
    cout << "options:" << endl;
    cout << "\t-m methodID\t (default 1; 1 for LWBMatch, 2 for DTW )" << endl
         << "\t-r th_RT\t (default 500; RT tolerance; seconds)" << endl
         << "\t-z th_MZ\t (default 0.2; MZ tolerance; Th)" << endl
         << "\t-i th_intensity\t (default 100000000; intensity tolerance)" << endl
         << "\t-n fileNum\t (default 2; the number of sample files)" << endl
         << "\t-l fileListPath\t (default \"\"; the path of the file containing all the sample files)" << endl
         << "\t-o outpath\t (default \"\"; the result file)" << endl
         << "\t-g groundtruth\t (default \"\"; the path of the groundtruth)" << endl
         << "\t-s separate\t (default 0; separate point for combined datasets. e.g. A1 A2 A3 B1 B2 B3 B4, separatePoint = 3)"
         << endl
         << "\t-f faster\t (default 1; set as 1, apply more than 4GB memory for the weight matrix; set as 0, utilizing the sparsity of the weight matrix, ~3 times slower)"
         << endl
         << "\t-w warper\t (default 2; 0, no warping; 1, LOWESS warping; 2, MS1-based DTW; 3, MS1+MS2-based DTW; 4, MS2-based DTW Matrix superposition; 5, MS2 spectra superpososition, DTW); 6, Product of MS2 spectra and a gaussian random vector"
         << endl
         << "\t-similarity \t (default 1; 1, dot product; 2, Lp-Norm (p=2 by default); 3, Shared peak count; 4, Pearson Correlation Coefficient); 5, Jensen-Shannon Divergence"
         << endl
         << "\t-pnorm \t (default 2.0; p >= 1)" << endl
         << "\t-filter \t (default Dummy; Dummy; PeakSquareRoot; PeakToProb; AboveMean; AboveMedian; Top50)" << endl
         << "\t-ms2filter \t (default Dummy; Dummy; NormalizeByMax; Top50)" << endl
         << "\t-mergeMethod \t (default Arithmetic; Arithmetic; Harmonic; Geometric;)" << endl;
    cout << "Example:" << endl;
    cout << "./lwbmatch -n 2 -l /path/to/your/FeatureXMLList.txt" << endl;
    cout << endl;
    cout << "Usage 2:" << endl;
    cout << "\tlwbmatch methodID th_RT th_MZ th_intensity fileNum fileListPath outputPath gtPath separatePoint" << endl;
    cout << "or" << endl;
    cout << "\tlwbmatch methodID th_RT th_MZ th_intensity fileNum fileListPath outputPath" << endl;
    cout << "Example:" << endl;
    cout
            << "./lwbmatch 1 500 0.2 100000000 20 /path/to/your/FeatureXMLList.txt /path/to/your/resutfile /path/to/the/groundtruth 10"
            << endl;
    cout << endl;

    /*cout<<endl;
    cout<<"*****************************Usage***********************************"<<endl;
    cout<<"lwbmatch methodID th_RT th_MZ th_intensity fileNum fileListPath outputPath gtPath separatePoint"<<endl;
    cout<<"Or"<<endl;
    cout<<"lwbmatch methodID th_RT th_MZ th_intensity fileNum fileListPath outputPath"<<endl;
    cout<<endl;
    cout<<"Example:"<<endl;
    cout<<"/LWBMatch/build/lwbmatch 1 500 0.2 1000000000 20 /evaluation/dataset/Standard_Protein_Mix_Database/Mix_4+Mix_7/FeatureXMLList /evaluation/result/Standard_Protein_Mix_Database/Mix_4+Mix_7/Mix_4+Mix_7_1.resu /evaluation/dataset/Standard_Protein_Mix_Database/Mix_4+Mix_7/Mix_4+Mix_7_GroundTruth.resu 10"<<endl;
    cout<<"Or"<<endl;
    cout<<"/LWBMatch/build/lwbmatch 1 500 0.2 1000000000 20 /evaluation/dataset/Standard_Protein_Mix_Database/Mix_4+Mix_7/FeatureXMLList /evaluation/result/Standard_Protein_Mix_Database/Mix_4+Mix_7/Mix_4+Mix_7_1.resu"<<endl;
    cout<<endl;
    cout<<"methodID: 			1 for LWBMatch, 2 for DTW"<<endl;
    cout<<"th_RT,th_MZ,th_intensity: 	tolerance"<<endl;
    cout<<"fileNum: 			the number of sample files"<<endl;
    cout<<"fileListPath: 			the path of the file containing all the sample files' paths"<<endl;
    cout<<"outputPath: 			the path of the file to store the alignment results"<<endl;
    cout<<"gtPath: 			the path of the file containing ground truth"<<endl;
    cout<<"separatePoint: 			separate point for combined datasets. For example, a combined dataset consists of 3 samples from the first dataset and 4 samples from second datasets,like A1 A2 A3 B1 B2 B3 B4, so the separate point is 3. This parameter is for heterogeneous alignment evaluation. " <<endl;
    cout<<"*********************************************************************"<<endl;
    cout<<endl;*/
}


int main(int argc, char **argv) {
    parameters *pParam = parameters::GetParam();
    bool debug = false;
    char **nargv = new char *[7];
    for (int i = 0; i < 7; i++) {
        nargv[i] = new char[400];
    }
    strcpy(nargv[0], "lwbmatch");
    strcpy(nargv[1], "-l");
    strcpy(nargv[2], "F:\\SpecBackUp\\wulongspec\\working\\ms2Align\\LongSwath_UPS1_1ug_rep2.featureXML_List");
    strcpy(nargv[3], "-g");
    strcpy(nargv[4], "F:\\SpecBackUp\\wulongspec\\LWBMatch\\LongSwath_UPS1_1ug_rep1_2_groundtruth.txt.txt");
    strcpy(nargv[5], "-o");
    strcpy(nargv[6], "F:\\SpecBackUp\\wulongspec\\working\\ms2Align\\LongSwath_UPS1_1ug_rep2_1.consensus.resu");



//	clock_t start = clock();
    if (argc == 1 && debug == false) {
        printUsage();
        return 0;
    }


    /*double th_RT,th_MZ,th_intensity;
    int fileNum,separatePoint,methodID;
    char fileListPath[500],outputPath[500],gtPath[500];
    if(argc == 8 || argc ==10) {*/
    try {
        if (debug) {
            cout << "Running the default param" << endl;
            pParam->ParseCMD(7, (char **) nargv);
        }
        else {
            pParam->ParseCMD(argc, argv);
        }
        pParam->print();
        /*///working function: 1 for LWBMatch, 2 for DTW
        sscanf(argv[1],"%d",&methodID);
        ///parameters for LWBMatch
        sscanf(argv[2],"%lf",&th_RT);
        sscanf(argv[3],"%lf",&th_MZ);
        sscanf(argv[4],"%lf",&th_intensity);

        ///samples
        sscanf(argv[5],"%d",&fileNum);
        sscanf(argv[6],"%s",fileListPath);
        sscanf(argv[7],"%s",outputPath);
        if(argc == 10)sscanf(argv[8],"%s",gtPath);
        else gtPath[0]=0;

        ///separate point for combined datasets. For example, a combined dataset consists of 3 samples from the first dataset and
        ///4 samples from second datasets,like A1 A2 A3 B1 B2 B3 B4, so the separate point is 3. This parameter is for heterogeneous alignment evaluation.
        if(argc == 10 )sscanf(argv[9],"%d",&separatePoint);
        else separatePoint=0;*/
    }
    catch (exception e) {
        printUsage();
        return -1;
    }
    SimpleTimer st("LWBMatch");
    AlignmentEvaluation *evaluation = new AlignmentEvaluation();
    evaluation->runEvaluation();
    //evaluation->runEvaluation(methodID,th_RT,th_MZ,th_intensity,fileNum,fileListPath,outputPath,gtPath,separatePoint);
    delete evaluation;


//    clock_t end = clock();
//    cout << "[Info] --------- Time Used: " << (end - start)*1.0/CLOCKS_PER_SEC << " s. ---------" << endl;
    return 0;
}
