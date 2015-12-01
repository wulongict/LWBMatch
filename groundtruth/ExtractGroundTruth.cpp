/*
written:  Jijie Wang
modified: Landon
 */
#include <iostream>
#include <vector>

#include <algorithm>

#include "../librarymsms/FileReader.h"
using namespace FileReader;
using namespace std;

void PrintUsage()
{
	cout << "Usage:" << endl;
	cout << "./ExtractGroundTruth mzXMLfileNum mzXMLList pepXMLPath threshold groundTruthPath percent" << endl;
	cout << "\n" << endl;
}
int main(int argc, char * argv[])
{
	//Spectrum_query  *querys = new Spectrum_query[maxn];
	if(argc==1){PrintUsage(); return 0;}
	//int i,j,k,t,p,q,T,x,y,sum,index;
	double threshold,percent=0;
	int fileNum,featureNum,queryNum;
	string mzXMLFileList,pepXMLfile,groundTruthPath;
	if (argc == 7)
	{
		sscanf(argv[1], "%d", &fileNum);
		mzXMLFileList = argv[2];
		pepXMLfile = argv[3];

		sscanf(argv[4], "%lf", &threshold);
		groundTruthPath = argv[5];
		sscanf(argv[6], "%lf", &percent);
	}
	else
	{
        PrintUsage();
        return 0;
		//fileNum = 2;
		//fileListPath = "F:\\SpecBackUp\\wulongspec\\LWBMatch\\mzXMLs.txt";
		//pepXmlPath = "F:\\SpecBackUp\\wulongspec\\LWBMatch\\interact-combine.ipro.pep.xml";
		//threshold = 0.9;
		//scanNum2RTPath = "F:\\SpecBackUp\\wulongspec\\LWBMatch\\LongSwath_UPS1_1ug_rep1_2_ScanNum2RT.resu";
		//groundTruthPath = "F:\\SpecBackUp\\wulongspec\\LWBMatch\\LongSwath_UPS1_1ug_rep1_2_groundtruth.txt";

		//percent = 0.1;
	}

	// First, Read Scan2RT
	//cerr << "1. ReadScan 2 RT" << endl;
	//Reader * scan2rt = new Scan2RTReader(scanNum2RTPath);
	//scan2rt->parse();


	// Second, Read mzXML file List
	//cerr << "2. Read mzXML file List " << endl;
	cerr << "Step 1/2, Read mzXML file list" << endl;
	mzXMLListReader* mzxml = new mzXMLListReader(mzXMLFileList);
	mzxml->parse();

	// Third, Read PepXML file
	cerr << "Step 2/2. Read ipro.pep.xml file" << endl;
	pepxmlReader *pepxml = new pepxmlReader(pepXMLfile);
	pepxml->parse();
	pepxml->GetGroundTruth( mzxml->getshortnames() );
	pepxml->OutputGroundTruth(groundTruthPath, mzxml->getshortnames());

	// Fourth, Read the groundtruth file
	//groundTruthPath = "F:\\SpecBackUp\\wulongspec\\LWBMatch\\groundtruth_2item.txt";
	//FileReader::groundtruthReader * groundtruth = new FileReader::groundtruthReader(groundTruthPath );
	//groundtruth->parse();
	//groundtruth->print();

	// Fifth, Read the consensus file
	//string consensusfile =   "F:\\SpecBackUp\\wulongspec\\working\\ms2Align\\LongSwath_UPS1_1ug_rep2_1.consensus.resu";
	//string consensusfile = "F:\\SpecBackUp\\wulongspec\\LWBMatch\\consensus_2item.resu";
	//FileReader::ConsensusReader * consensus = new FileReader::ConsensusReader(consensusfile);
	//consensus->parse();
	//consensus->print();
	//cout << *(consensus->GetPeakFeature(0, 0)) << endl;

	// Sixth, Get recall  and precision
	//groundtruth->ValidateConsensus(consensus);
	//groundtruth->CalculateRecallPrecision(consensus);


	return 0;
}


