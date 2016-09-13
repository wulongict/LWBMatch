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

#ifndef FEATUREMAP_H
#define FEATUREMAP_H

#include <string>
#include <vector>

using namespace std;

class Feature {
    // Question: Why each feature could be both one single feature, rt, mz, intensity
    // and could also be a list of feature RTlist, mzlist, intlist?
    // This is used for making consensus set.
    // Inherited from Jijie's code.

public:
    double rt, mz, intensity;
    int consensusNum;
    string peptideSequence;
    vector<string> FileList;
    vector<double> RTList;
    vector<double> MZList;
    vector<double> IntList;
};

class FeatureMap {
public:
    //ToDo: It's crash when  we try to use features.clear()
    //I thought the bug is : use the vector before it is initialized.
    FeatureMap();

    vector<Feature> features;
    int m_featureNum;

    void sortByMZ();

    void sortByRT();

    void readFeatureMap(char *FileName);
//    {
//        if (contains(FileName, "features.tsv"))
//        {
//            readFeatureMapDinosaur(FileName);
//        }
//        else if (contains(FileName, "featureXML"))
//        {
//            readFeatureMapFeatureXML(FileName);
//        }
//        else
//        {
//            throw "[Error] Invalid feature file type.";
//        }
//    }

    void readFeatureMapFeatureXML1(char *FileName);

    vector<double> splitAsDouble(char *line);

    void readFeatureMapDinosaur(char *FileName);
};

// Now we add another FeatureMap Reader for Dinosaur feautre


#endif // FEATUREMAP_H
