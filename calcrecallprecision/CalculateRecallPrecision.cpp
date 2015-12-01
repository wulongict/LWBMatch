/*
ID: jijie.w1
PROG:
LANG: C++
 */
#include <iostream>
#include <vector>

#include <algorithm>

#include "../librarymsms/FileReader.h"
#include <iostream>


using namespace std;

void PrintUsage()
{
	cout << "Usage:" << endl;
	cout << "./CalculateRecallPrecision groundtruthfile consensusfile mz_th[default:0.1][recommend:0.1~0.5] rt_th[default:5][recommend:5~50]" << endl;
}
int main(int argc, char * argv[])
{
	if(argc==1){PrintUsage(); return 0;}

	string groundtruthfile, consensusfile;
	double mz_th = 0.1;
	double rt_th = 5;
	if (argc >= 3)
	{
		groundtruthfile = argv[1];// input ground truth
		consensusfile = argv[2];
		mz_th = 0.1;
	}
	if(argc >= 4)
	{
		mz_th = stod(argv[3]);
	}
	if(argc >= 5)
	{
		rt_th = stod(argv[4]);
	}
	if (argc < 3)
	{
		PrintUsage();
		return 0;
	}

	
	FileReader::groundtruthReader * groundtruth = new FileReader::groundtruthReader(groundtruthfile);
	groundtruth->parse();
	groundtruth->Set_mz_threshold(mz_th);
	groundtruth->Set_rt_threshold(rt_th);
	groundtruth->CalcMassDiff();
	

	FileReader::ConsensusReader * consensus = new FileReader::ConsensusReader(consensusfile);
	consensus->parse();
	consensus->Set_mz_threshold(mz_th);
	consensus->Set_rt_threshold(rt_th);
	consensus->CalcMassDiff();


	groundtruth->ValidateConsensus(consensus);
	//groundtruth->CalculateRecallPrecision(consensus);
	consensus->ValidateConsensus(groundtruth);
	groundtruth->CalculateRecallPrecision(consensus);

	return 0;
}


