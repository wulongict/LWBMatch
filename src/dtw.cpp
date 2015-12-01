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

#include "dtw.h"

int DTW::getMapIDWithHighestFeatureNumbers(){
	int maxNum = 0, maxIndex = 0;
	for (int i = 0; i < sampleNum; i++){
		if (samples[i].featureNum > maxNum) maxNum = samples[i].featureNum, maxIndex = i;
	}
	return maxIndex;
}

void DTW::applyMemory(int featureNum){
	weightMat = (int(*)[MAXN])malloc(sizeof(int)*(featureNum*MAXN));
}


int DTW::sharePeakCount(Feature & f1, Feature &f2)
{
	if (fabs(f1.mz - f2.mz) < th_MZ&&fabs(f1.intensity - f2.intensity) < th_intensity) return 1;
	//     if(fabs(f1.rt-f2.rt)<th_RT&&fabs(f1.mz-f2.mz)<th_MZ&&fabs(f1.intensity-f2.intensity)<th_intensity) return 1; 
	else return 0;
}

int DTW::sigmod(Feature & f1, Feature &f2)
{
	if (fabs(f1.mz - f2.mz) > th_MZ) return 0;
	double dist = fabs(f1.mz - f2.mz)*1.0 / th_MZ;
	double sigmod_ = 1 / (1 + pow(2.71828, 18 * dist - 6));
	return floor(1000 * (sigmod_));
}

int DTW::score_func(Feature &f1, Feature &f2, ScoreFuntions sf)
{
	switch (sf)
	{
	case SPC:
		return sharePeakCount(f1, f2);
	case SIGMOD:
		return sigmod(f1, f2);
	default:
		return 0;
	}
}

double DTW::DynamicTimeWarpping(FeatureMap &reference, FeatureMap & sample){
	//
	//                    1. weightMat[i-1][j]
	//weightMat[i][j]=max 2. weightMat[i-1][j-1] + score_func(sample.features[i],reference.features[j],sf)
	//                    3. weightMat[i][j-1]
	//
	reference.sortByRT();
	sample.sortByRT();
	int i, j;
	//ScoreFuntions sf = SPC;
	ScoreFuntions sf = SIGMOD;

    // initialize:
	for (i = 0; i < sample.featureNum; i++){
		for (j = 0; j < reference.featureNum; j++)
			weightMat[i][j] = 0;
	}
	weightMat[0][0] = score_func(sample.features[0], reference.features[0], sf);

	for (i = 1; i < reference.featureNum; i++)
	{
        weightMat[0][i] = score_func(sample.features[0], reference.features[i], sf);
	}
	for (i = 1; i < sample.featureNum; i++)
	{
        weightMat[i][0] = score_func(sample.features[i], reference.features[0], sf);
	}

    // There is a bug here
    //ToDo:
    // Why the first row and first column is not calculate using the dp formula?
	for (i = 1; i < sample.featureNum; i++)
		for (j = 1; j < reference.featureNum; j++)
		{
			// The layout of s1, s2, s3, s4
			// s3[i-1,j-1]	s1[i-1, j]
			// s2[i,j-1]	s4[i,j]
			int s4 = score_func(sample.features[i], reference.features[j], sf);
			int s1 = (weightMat[i - 1][j]);
			int s2 = (weightMat[i][j - 1]);
			int s3 = (weightMat[i - 1][j - 1] + s4); // What ? Why s3 + s4 -> s3
			if (s1 > s2 && s1 > s3)
			{
				weightMat[i][j] = weightMat[i - 1][j];
			}
			else if (s2 > s1 && s2 > s3)
			{
				weightMat[i][j] = weightMat[i][j - 1];
			}
			else
			{
				weightMat[i][j] = weightMat[i - 1][j - 1] + s4;
			}
		}
	i = sample.featureNum - 1; j = reference.featureNum - 1;


	//trace back to find the matching
	while (i != 0 && j != 0){
		int s4, s1, s2, s3;
		s4 = score_func(sample.features[i], reference.features[j], sf);
		s1 = (weightMat[i - 1][j]);
		s2 = (weightMat[i][j - 1]);
		s3 = (weightMat[i - 1][j - 1] + s4);
		if (s1 > s2 && s1 > s3)
		{
			i = i - 1;
		}
		else if (s2 > s1 && s2 > s3)
		{
			j = j - 1;
		}
		else
		{
			if (s4 != 0){
				// 		sum++; 
				match1[i] = j;
				match2[j] = i;
				// 		cout<<reference[i].rt<<" "<<query[j].rt<<endl;
			}
			i = i - 1; j = j - 1;
		}
	}
	if (score_func(sample.features[i], reference.features[j], sf)){
		match1[i] = j;
		match2[j] = i;
	}
	//     freopen("/extra2/jwang/weightMat.txt","w",stdout);
	//     for (i = 1; i < sample.featureNum; i++){
	//         for (j = 1; j < reference.featureNum; j++){
	// 	    cout<<weightMat[i][j]<<" ";
	// 	}
	// 	cout<<endl;
	//     }
	//     freopen("/dev/tty", "w",stdout);
	return weightMat[sample.featureNum - 1][reference.featureNum - 1];
}

int DTW::makeConsensus(FeatureMap &reference, FeatureMap &sample){
	int i, matches = 0;
	for (i = 0; i < reference.featureNum; i++) {
		if (match2[i] != -1) {
			matches++;
			reference.features[i].rt = (reference.features[i].rt*reference.features[i].consensusNum + sample.features[match2[i]].rt) / (reference.features[i].consensusNum + 1);
			reference.features[i].mz = (reference.features[i].mz*reference.features[i].consensusNum + sample.features[match2[i]].mz) / (reference.features[i].consensusNum + 1);
			reference.features[i].intensity = (reference.features[i].intensity*reference.features[i].consensusNum + sample.features[match2[i]].intensity) / (reference.features[i].consensusNum + 1);
			reference.features[i].consensusNum++;
			reference.features[i].FileList.push_back(sample.features[match2[i]].FileList[0]);
			reference.features[i].RTList.push_back(sample.features[match2[i]].rt);
			reference.features[i].MZList.push_back(sample.features[match2[i]].mz);
			reference.features[i].IntList.push_back(sample.features[match2[i]].intensity);
		}
	}
	for (i = 0; i < sample.featureNum; i++) {
		if (match1[i] == -1) {
			reference.features.push_back(sample.features[i]);
			reference.featureNum++;
		}
		// 	    out(reference.featureNum);
	}
	return matches;
}

FeatureMap  DTW::runDTW(double th_RT, double th_MZ, double th_intensity){
	this->th_RT = th_RT;
	this->th_MZ = th_MZ;
	this->th_intensity = th_intensity;

	int i, j, k, t, T, x, y, sum, index, min_x, min_y;
	FeatureMap sample, reference;
	//find the map with highest number of features
	min_y = getMapIDWithHighestFeatureNumbers();

	for (t = 0; t < sampleNum; t++){
		if (min_y == t) continue;
		min_x = t;

		sample = samples[min_x];
		reference = samples[min_y];

		applyMemory(sample.featureNum);

		for (i = 0; i < MAXN; i++) match1[i] = match2[i] = -1;

		double value = DynamicTimeWarpping(reference, sample);

		int matches = makeConsensus(reference, sample);

		cout << "matches:" << matches << endl;

		samples[min_y] = reference;

		free(weightMat);
	}

	return reference;
}

