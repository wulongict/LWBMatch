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

#include "featuremap.h"
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;

bool mzComp(Feature a, Feature b)
{
    if (a.mz==b.mz) return a.rt<b.rt;
    else return a.mz < b.mz;
}

FeatureMap::FeatureMap()
{
	features = vector<Feature>();
	featureNum = 0;
}

void FeatureMap::sortByMZ()
{
    stable_sort(features.begin(),features.end(),mzComp);
}

bool rtComp(Feature a, Feature b)
{
    if (a.rt==b.rt) return a.mz<b.mz;
    else return a.rt < b.rt;
}

void FeatureMap::sortByRT()
{
    stable_sort(features.begin(),features.end(),rtComp);
}


char line[50005],in1[200],in2[200],in3[200];
bool contains(char *line,const char *word) {
    int i,j,lenx=strlen(line),leny=strlen(word);
    for (i=0;i<lenx;i++) {
        for (j=0;j<leny;j++) {
            if (line[i+j]!=word[j]) break;
        }
        if (j==leny) return true;
    }
    return false;
}
void FeatureMap::readFeatureMap(char * shortFileName)
{
    int i,j,k,t,T,index=0;
    double rt=-1,mz=-1,intensity;
    string sshortFileName(shortFileName);
    featureNum=0;
    features.clear();
    bool unique = false;
    int subfeature=0;
    while (fgets(line,sizeof(line),stdin)) {
        if (contains(line,"<subordinate>")) {
            subfeature ++;
        }
        else if (contains(line,"</subordinate>")) {
            subfeature --;
        }


        if (subfeature!=0) continue;
        if (contains(line,"count=")) {
            for (i=0;i<strlen(line);i++) {
                if (line[i]=='"') line[i]=' ';
            }
            sscanf(line,"%s%s%d%s",in1,in2,&featureNum,in3);
            features.reserve(featureNum);
            features.resize(featureNum);
            cout<<"[Info] "<<shortFileName<<"  #Features: "<<featureNum<<endl;
        }
        else if (contains(line,"<position dim")) {

            for (i=0;i<strlen(line);i++) {
                if (line[i]=='<'||line[i]=='>') line[i]=' ';
            }
            if (rt==-1) sscanf(line,"%s%s%lf%s",in1,in2,&rt,in3);
            else if (mz==-1) sscanf(line,"%s%s%lf%s",in1,in2,&mz,in3);
            if (rt!=-1&&mz!=-1) {
				// Check for redudency!
                for (i=0;i<index;i++) {
                    if (fabs(features[i].rt-rt)<1 && fabs(features[i].mz-mz) < 0.01) break;
                }
                if (i<index) unique = false;
                else unique = true;
                features[index].rt=rt;
                features[index].mz=mz;

				// Reset the rt and mz to -1
				rt=-1,mz=-1;
            }
        }
        else if (contains(line,"<intensity>")) {
            for (i=0;i<strlen(line);i++) {
                if (line[i]=='<'||line[i]=='>') line[i]=' ';
            }
            sscanf(line,"%s%lf%s",in1,&intensity,in3);
			// Question: So the unique variable is useless here
            if (1) {//if(unique)
                features[index].consensusNum=1;
                features[index].FileList.push_back(sshortFileName);
                features[index].RTList.push_back(features[index].rt);
                features[index].MZList.push_back(features[index].mz);
                features[index].IntList.push_back(intensity);
                features[index++].intensity=intensity;
                unique = false;
            }
        }

    }
    featureNum = index;

}
