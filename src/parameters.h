//
// Created by wulong on 7/16/15.
//

#ifndef RT_ALIGNMENT_PARAMETERS_H
#define RT_ALIGNMENT_PARAMETERS_H
#include <string>
#include <iostream>
#include <sstream>
//#include <bits/stl_map.h>
#include <map>

using namespace std;

class ParamOption
{
public:
    ParamOption(string shortOpt, string longOpt):shortOption(shortOpt), longOption(longOpt){}
    bool operator== (string inputOption)
    {
        return shortOption == inputOption || longOption == inputOption;
    }
    string shortOption;
    string longOption;
};

class dummyParam{
    map<ParamOption,string> ParamMap;
public:
    dummyParam(){

    }
};

class parameters{
private:
    string lwbmatch_bin_path;
    static parameters* pParam;
    parameters()
    {
        lwbmatch_bin_path = "";

        methodID = 1;
        th_RT = 500;
        th_MZ = 0.2;
        th_intensity = 100000000;

        fileNum = 2;
        separatePoint = 0;
        fileListPath = "";
        outputPath = "";
        gtPath = "";

        faster = "1";

        Warper = "2";
    }
public:
    static parameters * GetParam()
    {
        if (pParam == NULL)
            pParam = new parameters();
        return pParam;
    }


    void print()
    {
        cout << "Parameters: " << methodID << " " << th_RT << " " << th_MZ << " " << th_intensity << " " << fileNum << " "
        << fileListPath << " " << outputPath << " " << gtPath << " " << separatePoint << endl;

    }
    string getPath(string file_path)
    {
        int found =  file_path.find_last_of("\\");
        if (found == string::npos)
        {
            found = file_path.find_last_of("/");
        }
        if (found == string::npos)
            return "./";
        else
            return file_path.substr(0, found+1);
    }

    const string &getLwbmatch_bin_path() const {
        return lwbmatch_bin_path;
    }

    void setLwbmatch_bin_path(const string &lwbmatch_bin_path) {
        parameters::lwbmatch_bin_path = lwbmatch_bin_path;
    }

    void SetOption(string option, string value)
    {
        istringstream iss(value);

        if (option == "m" || option == "method")
        {
            iss >> methodID;
        }
        else if (option == "r" || option == "rt")
        {
            iss >> th_RT;
        }
        else if (option == "z" || option == "mz")
        {
            iss >> th_MZ;
        }
        else if (option == "i" || option == "intensity")
        {
            iss >> th_intensity;
        }
        else if (option == "n" || option == "filenum")
        {
            iss >> fileNum;

        }
        else if (option == "g" || option == "gtpath")
        {
            iss >> gtPath;
        }
        else if (option == "o" || option == "outputpath")
        {
            iss >> outputPath;
        }
        else if (option == "l" || option == "filelist")
        {
            iss >> fileListPath;
        }
        else if (option == "s" || option == "separatepoint")
        {
            iss >> separatePoint;
        }
        else if (option == "f" || option == "faster")
        {
            iss >> faster;
        }
        else if (option == "w" || option == "warper")
        {
            iss >> Warper;
        }


        else
        {
            cout << "Bad option value pair: <" << option << ", " << value << ">." << endl;
            throw "Program fails to parse the option!";
        }
    }
    void ParseCMD(int argc, char ** argv)
    {
        lwbmatch_bin_path = getPath(argv[0]);
        if (argv[1][0] == '-' )
        {
            ParseCMDWithOption(argc, argv);
            checkparam();
        }
        else
        {
            ParseCMDOldVersion(argc, argv);
        }
    }
    void checkparam()
    {
        if (outputPath == "") // User does not define the output path.
        {
            int found = fileListPath.find("featureXMLList");
            if(found == string::npos)
            {
                cout << "[Error] could not find \"featureXMLList in input filename\"" << endl;
                outputPath = fileListPath + ".resu";
            }
            else
            {
                outputPath = fileListPath.substr(0,found) + "resu";
            }
        }


    }
    void ParseCMDOldVersion(int argc, char **argv)
    {
        if (argc < 8)
        {
            cout << "Please provide seven parameters at least!" << endl;
            throw "Too less parameters!";
        }
        else if (argc >= 8)
        {
            methodID = atoi(argv[1]);
            th_RT = atof(argv[2]);
            th_MZ = atof(argv[3]);
            th_intensity = atof(argv[4]);

            fileNum = atoi(argv[5]);
            fileListPath = argv[6];
            outputPath = argv[7];

            if (argc == 10)
            {
                gtPath = argv[8];
                separatePoint = atoi(argv[9]);

            }
        }
        else
        {
            cout << "Too more parameters or incorrect parameters!" << endl;
            throw "Please Check the parameters: no more than 9 parameters!";
        }
    }
    void ParseCMDWithOption(int argc, char ** argv)
    {
        int i = 1;
        string option = "";
        while (i < argc)
        {
            if (argv[i][0] == '-' && option == "")
            {
                argv[i]++;
                option = argv[i];
            }
            else {
                SetOption(option, argv[i]);
                option = "";
            }
            i++;
        }
    }
private:
    double th_RT;
    double th_MZ;
    double th_intensity;

    int methodID;
public:
    double getTh_RT() const {
        return th_RT;
    }

    void setTh_RT(double th_RT) {
        parameters::th_RT = th_RT;
    }

    double getTh_MZ() const {
        return th_MZ;
    }

    void setTh_MZ(double th_MZ) {
        parameters::th_MZ = th_MZ;
    }

    double getTh_intensity() const {
        return th_intensity;
    }

    void setTh_intensity(double th_intensity) {
        parameters::th_intensity = th_intensity;
    }

    int getMethodID() const {
        return methodID;
    }

    void setMethodID(int methodID) {
        parameters::methodID = methodID;
    }

    int getFileNum() const {
        return fileNum;
    }

    void setFileNum(int fileNum) {
        parameters::fileNum = fileNum;
    }

    int getSeparatePoint() const {
        return separatePoint;
    }

    void setSeparatePoint(int separatePoint) {
        parameters::separatePoint = separatePoint;
    }

    const string &getFileListPath() const {
        return fileListPath;
    }

    void setFileListPath(const string &fileListPath) {
        parameters::fileListPath = fileListPath;
    }

    const string &getOutputPath() const {
        return outputPath;
    }

    void setOutputPath(const string &outputPath) {
        parameters::outputPath = outputPath;
    }

    const string &getGtPath() const {
        return gtPath;
    }

    void setGtPath(const string &gtPath) {
        parameters::gtPath = gtPath;
    }

    const string &getFaster() const {
        return faster;
    }

    void setFaster(const string &faster) {
        parameters::faster = faster;
    }

private:
    int fileNum;
    int separatePoint;

    string fileListPath;
    string outputPath;
    string gtPath;

    string faster;

    string Warper;

public:
    const string &getWarper() const {
        return Warper;
    }

    void setWarper(const string &Warper) {
        parameters::Warper = Warper;
    }
};


#endif //RT_ALIGNMENT_PARAMETERS_H
