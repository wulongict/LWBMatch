#ifndef _FILEREADER_H_
#define _FILEREADER_H_

#ifdef _WIN32
#define PATH_SLASH '\\'
#else
#define PATH_SLASH '/'
#endif
#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <string.h>
#include <memory.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;
namespace FileReader
{
    double calcmean(vector<double> d)
    {
        if (d.size() == 0)
        {
            cout << "Can not calculate mean of empty vector! Program exits!" << endl;
            exit(0);
        }
        double s = 0;
        for (unsigned int i = 0; i < d.size(); i++)
        {
            s += d[i];
        }
        return s / d.size();
    }
    double calcSTD(vector<double> d)
    {
        double dmean = calcmean(d);
        double sum_of_square = 0;// square(d)
        for (unsigned int i = 0; i < d.size(); i++)
        {
            sum_of_square += (d[i] - dmean)*(d[i] - dmean);
        }
        return sqrt(sum_of_square/d.size());
    }
class parameters
{
    map<string, string> m_param;
    int argc;
    char ** argv;
public:
    parameters() {}
    virtual void InitParaMap(string options)
    {
        // options is split by space.
        string option;
        istringstream iss(options);
        while (iss >> option)
        {
            if (option=="")
            {
                cout << "Empty option is invalid!" << endl;
                break;
            }
            m_param[option] = "";
        }

    }//put something into the map
    virtual void Help()
    {
        cout << "Usage:" << endl;
        // ToDo:
        // Add usage function here.
        cout << "***.exe ";
        for (map<string, string>::iterator it = m_param.begin(); it != m_param.end(); it ++)
        {
            cout << " -" << it->first << " [value]";
        }
        cout << endl;
    }
    virtual void Print()
    {

        for (map<string, string>::iterator it = m_param.begin(); it != m_param.end(); it++)
        {
            cout <<"[Param] " << it->first << " =  " << it->second << endl;;
        }

    }

    virtual ~parameters() {}

    virtual void ParseCMD(int argc, char **argv) {}

    virtual string Get(string option)
    {
        if (m_param.find(option)!=m_param.end())
        {
            return m_param[option];
        }
        else
        {
            cout << "No value for option: " << option << endl;
            exit(0);
        }
    }

    virtual void ParseCMDWithOption(int argc, char ** argv)
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
            else
            {
                Set(option, argv[i]);
                option = "";
            }
            i++;
        }
    }
private:
    virtual void Set(string option, string value)
    {
        if (m_param.find(option) != m_param.end())
        {
            m_param[option] = value;
        }
        else
        {
            cout << "Unexpected option: " << option << endl;
            exit(0);
        }
    }




};


class LWBMatchParam : public parameters
{
public:
    double th_RT;
    double th_MZ;
    double th_intensity;

    int methodID;
    int fileNum;
    int separatePoint;

    string fileListPath;
    string outputPath;
    string gtPath;
    LWBMatchParam()
    {
        methodID = 1;
        th_RT = 500;
        th_MZ = 0.2;
        th_intensity = 100000000;

        fileNum = 2;
        separatePoint = 0;
        fileListPath = "";
        outputPath = "";
        gtPath = "";
    }
    void print()
    {
        cout << "Parameters: " << methodID << " " << th_RT << " " << th_MZ << " " << th_intensity << " " << fileNum << " "
             << fileListPath << " " << outputPath << " " << gtPath << " " << separatePoint << endl;

    }
    void ParseCMD(int argc, char ** argv)
    {
        if (argv[1][0] == '-')
        {
            ParseCMDWithOption(argc, argv);
        }
        else
        {
            ParseCMDOldVersion(argc, argv);
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

    void Set(string option, string value)
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
        else
        {
            cout << "Bad option value pair: <" << option << ", " << value << ">." << endl;
            throw "Program fails to parse the option!";
        }
    }
};

vector<string> SplitByFirst(string input, char firstchar)
{
    vector<string> ret;
    int found = input.find_first_of(firstchar);
    if (found != string::npos)
    {
        ret.push_back(input.substr(0, found));
        ret.push_back(input.substr(found + 1));
    }
    else
    {
        cout << "Can not find " << firstchar << " in " << input << ". Please check!" << endl;
        exit(0);
    }
    return ret;
}

class PSM
{
public:
    PSM()
    {
        pmass = 1.007276;
        precursor_neutral_mass= 0;
        calc_neutral_pep_mass = 0;
        retention_time_sec = 0;
        probability = 0;
        iprobability = 0;
        mz = 0;
        isDecoy = false;
        spectrum = "";
        peptide = "";
        scan = "";
        assumed_charge = 0;
    }
    static bool asdsortbyPeptide(const PSM & x, const PSM & y)
    {
        return x.peptide.compare(y.peptide) < 0;
    }
    friend ostream & operator<<(ostream &os, PSM &psm)
    {
        os <<psm.spectrum << " " << psm.scan  << " " << psm.assumed_charge << " " << psm.peptide << " "
           << psm.precursor_neutral_mass << " " << psm.mz << " " << psm.retention_time_sec << " " << psm.GetPepMZ()
           << " " << psm.iprobability << " " << psm.probability << " " << psm.isDecoy << endl;
        return os;
    }
    double GetPepMZ()
    {
        return calc_neutral_pep_mass / assumed_charge+pmass;
    }
    double GetPrecursorMZ()
    {
        return precursor_neutral_mass / assumed_charge + pmass;
    }
    double pmass;
    string spectrum;
    double precursor_neutral_mass;
    double calc_neutral_pep_mass;
    int assumed_charge;
    double retention_time_sec;
    string peptide;
    double probability;
    double iprobability;// iproprophet probability.
    double mz;
    string scan;
    bool isDecoy;
};


class PeakFeature
{
public:
    PeakFeature(double mz = 0, double rt = 0, double inten = 0)
    {
        SetMZIntenRT(mz, rt, inten);
    }
    PeakFeature(const PeakFeature &a)
    {
        SetMZIntenRT(a.GetMZ(), a.GetRT(), a.GetInten());
    }
    virtual ~PeakFeature() {}
    virtual void PrintID()
    {
        cout << "PeakFeature" << endl;
    }
    // Get
    double GetMZ() const
    {
        return m_mz;
    }
    double GetInten() const
    {
        return m_inten;
    }
    double GetRT() const
    {
        return m_rt;
    }
    virtual void ReadPeak(istream & ism)
    {}




    // Set
    void SetMZ(double mz)
    {
        m_mz = mz;
    }
    void SetInten(double inten)
    {
        m_inten = inten;
    }
    void SetRT(double rt)
    {
        m_rt = rt;
    }
    void SetMZIntenRT(double mz, double rt, double inten)
    {
        m_mz = mz;
        m_rt = rt;
        m_inten = inten;
    }

    bool isNonZero()
    {
        if (m_mz < 0.1)
        {
            return false;
        }
        return true;
    }

    // friend
    friend istream & operator >>(istream & ism, PeakFeature & pf)
    {
        pf.ReadPeak(ism);
        return ism;
    }
    friend ostream & operator <<(ostream & osm, PeakFeature & pf)
    {
        osm << pf.GetRT() << " " << pf.GetMZ() << " " << pf.GetInten();
        return osm;
    }

private:
    double m_mz;
    double m_rt;
    double m_inten;
};

class PeakFeatureGT : public PeakFeature
{
public:
    PeakFeatureGT(double mz = 0, double rt = 0, double inten = 0) :PeakFeature(mz, rt, inten) {}
    void  ReadPeak(istream & ism)
    {
        double mz, rt;
        string linemz, linert;
        ism >> linert >> linemz;
        if (linemz == "NA")
            mz = 0;
        else
            mz = stod(linemz);

        if (linert == "NA")
            rt = 0;
        else
            rt = stod(linert);
        //ism >> rt >> mz;
        SetMZIntenRT(mz, rt, 0);
    }
    void PrintID()
    {
        cout << "PeakFeature GT" << endl;
    }

};
class PeakFeatureConsensus : public PeakFeature
{
public:
    void  ReadPeak(istream & ism)
    {
        double mz, rt, inten;
        string linemz, linert, lineinten;
        ism >> linert >> linemz >> lineinten;
        if (linemz == "NA")	mz = 0;
        else mz = stod(linemz);

        if (linert == "NA")	rt = 0;
        else rt = stod(linert);

        if (lineinten == "NA")	inten = 0;
        else inten = stod(lineinten);
        //ism >> rt >> mz >> inten;
        SetMZIntenRT(mz, rt, inten);
    }
    void PrintID()
    {
        cout << "PeakFeature consensus " << endl;
    }
};

class AlignedFeatures
{
    vector<PeakFeature> m_Peaks;
    int m_FileNum;
public:
    AlignedFeatures(int filenum)
    {
        m_FileNum = filenum;
    }
    AlignedFeatures(const AlignedFeatures & af)
    {
        m_FileNum = af.GetFileNum();
        for (int i = 0; i < m_FileNum; i++)
        {
            PeakFeature pf = af.GetPeakFeature(i);
            m_Peaks.push_back(pf);
        }

    }

    ~AlignedFeatures() {}
    int GetFileNum() const
    {
        return m_FileNum;
    }
    PeakFeature  GetPeakFeature(int i) const
    {
        return m_Peaks[i];
    }
    void AddFeature(double mz, double rt, double inten)
    {
        m_Peaks.push_back(PeakFeature(mz, rt, inten));
    }
    void AddFeature(PeakFeature pf)
    {
        m_Peaks.push_back(pf);
    }
    friend istream & operator >>(istream & ism, AlignedFeatures & af)
    {

        for (int i = 0; i < af.GetFileNum(); i++)
        {
            PeakFeature pf;
            ism >> pf;
            af.AddFeature(pf);
        }
        return ism;
    }
    friend ostream & operator <<(ostream & osm, AlignedFeatures &  af)
    {
        //PeakFeature pf;
        for (int i = 0; i < af.GetFileNum(); i++)
        {
            osm << " ";
            PeakFeature pf = af.GetPeakFeature(i);
            osm << pf;
        }
        return osm;
    }
};

class FeatureMatrix
{
public:
    // constructor
    FeatureMatrix(int AlignmentNum, int FileNum) :m_FileNum(FileNum), m_AlignmentNum(AlignmentNum)
    {
        // create line;
        m_fm = new PeakFeature**[AlignmentNum];
        // create each entry
        for (int i = 0; i < AlignmentNum; i++)
        {
            m_fm[i] = new PeakFeature*[FileNum];
        }
    }
    // destructor
    virtual ~FeatureMatrix()
    {
        for (int i = 0; i < m_AlignmentNum; i++)
        {
            for (int j = 0; j < m_FileNum; j++)
            {
                if (m_fm[i][j] != NULL)
                {
                    delete m_fm[i][j];
                    m_fm[i][j] = NULL;
                }
            }
            if (m_fm[i] != NULL)
            {
                delete[] m_fm[i];
                m_fm[i] = NULL;
            }
        }
        delete[] m_fm;


    }// be careful , make virtual destruction function.


    // get
    int GetAlignmentNum()
    {
        return m_AlignmentNum;
    }
    int GetFileNum()
    {
        return m_FileNum;
    }
    PeakFeature * GetPeakFeature(int i, int j)
    {
        return m_fm[i][j];
    }

    // set
    void SetPeakFeature(int i, int j, PeakFeature *a)
    {
        m_fm[i][j] = a;
    }
private:
    // Interface
    virtual void PrintID()
    {
        cout << "Feature Marix" << endl;
    }
    virtual void WritePrepare(ostream &osm, int AlignmentIndex) {}
    virtual void Readprepare(istream &ism) {}
    virtual PeakFeature * Create()
    {
        return NULL;
    }

public:
    // friend IO
    friend istream & operator >> (istream & ism, FeatureMatrix & fm)
    {
        for (int i = 0; i < fm.GetAlignmentNum(); i++)
        {
            fm.Readprepare(ism);
            for (int j = 0; j < fm.GetFileNum(); j++)
            {
                PeakFeature * pf = fm.Create();
                ism >> *pf;
                fm.SetPeakFeature(i, j, pf);
            }
        }
        return ism;
    }
    friend ostream & operator << (ostream & osm, FeatureMatrix & fm)
    {
        //AlignedFeatures af;
        for (int i = 0; i < fm.GetAlignmentNum(); i++)
        {
            fm.WritePrepare(osm, i);
            for (int j = 0; j < fm.GetFileNum(); j++)
            {
                osm << *(fm.GetPeakFeature(i, j)) << " ";
            }
            osm << endl;
        }
        return osm;
    }
private:
    int m_AlignmentNum; // row number
    int m_FileNum; // column number
    PeakFeature *** m_fm; // feature matrix
};

class GroundTruth : public FeatureMatrix
{
    vector<string> m_Peptides;
public:
    // consttructor
    GroundTruth(int AlignmentNum, int FileNum) : FeatureMatrix(AlignmentNum, FileNum) { }
    // destructor
    ~GroundTruth() {}

    // Identify yourself
    void PrintID()
    {
        cout << "Feature Marix ground truth" << endl;
    }
private:
    void Readprepare(istream & ism)
    {
        string pep;
        ism >> pep;
        AddPeptide(pep);
    }
    void WritePrepare(ostream & osm, int AlignmentIndex)
    {
        osm << GetPeptide(AlignmentIndex) << " ";
    }
    PeakFeature * Create()
    {
        return new PeakFeatureGT;
    }
    // Get from vector
    string GetPeptide(int i)
    {
        return m_Peptides[i];
    }
    // Add
    void AddPeptide(string pepSQ)
    {
        m_Peptides.push_back(pepSQ);

    }

};

class FeaturePair
{
public:
    static double m_th_rt;
    static double m_th_mz;
    static double m_th_inten;
public:
    // cosntructor
    FeaturePair(PeakFeature *a, PeakFeature *b) :x(a), y(b) {}
    // destructor
    ~FeaturePair() {}
    // equal function
    bool operator ==(const FeaturePair & other)
    {
        if (fabs(other.x->GetMZ() - x->GetMZ()) > FeaturePair::m_th_mz)
        {
            return false;
        }
        else if (fabs(other.y->GetMZ() - y->GetMZ()) > FeaturePair::m_th_mz)
        {
            return false;
        }
        else if (fabs(other.x->GetRT() - x->GetRT()) > FeaturePair::m_th_rt)
        {
            return false;
        }
        else if (fabs(other.y->GetRT() - y->GetRT()) > FeaturePair::m_th_rt)
        {
            return false;
        }
        else if (fabs(other.x->GetInten() - x->GetInten()) > FeaturePair::m_th_inten)
        {
            return false;
        }
        else if (fabs(other.y->GetInten() - y->GetInten()) > FeaturePair::m_th_inten)
        {
            return false;
        }
        return true;

    }
    PeakFeature * x;
    PeakFeature * y;
};
double FeaturePair::m_th_rt = 5;
double FeaturePair::m_th_mz = 0.1;
double FeaturePair::m_th_inten = 1e20;

class ValidPairs
{
    double eps;
public:
    vector<FeaturePair> m_pairs;
    ValidPairs()
    {
        eps = 1e-8;
    }
    virtual ~ValidPairs() {}
    void AddFeaturePair(PeakFeature *x, PeakFeature *y)
    {
        if (fabs(x->GetMZ()) < eps || fabs(y->GetMZ()) < eps)
        {
            return;
        }
        m_pairs.push_back(FeaturePair(x, y));
    }

};


class Consensus : public FeatureMatrix
{

public:
    Consensus(int AlignmentNum, int FileNum) : FeatureMatrix(AlignmentNum, FileNum) {}
    ~Consensus() {}
    void PrintID()
    {
        cout << "Feature Marix Consensus" << endl;
    }
private:
    // factory method
    PeakFeature * Create()
    {
        return  new PeakFeatureConsensus;
    }
};

class Reader
{
    string m_filename;
public:
    Reader(string input)
    {
        m_filename = input;
    }
    virtual ~Reader() {}
    virtual void parse() {}
    string GetFileName()
    {
        return m_filename;
    }
    virtual void Print() {}
};

class FeatureMatrixReader : public Reader
{
    int m_AlignmentNum;
    int m_FileNum;
    FeatureMatrix *m_featurematrix;
    vector<string> m_FileNames;

public:
    double getM_th_rt() const {
        return m_th_rt;
    }

    double getM_th_mz() const {
        return m_th_mz;
    }

    double getM_th_inten() const {
        return m_th_inten;
    }

private:
// three threholds for validation. is feature A  equal to features B?
    double m_th_rt;
    double m_th_mz;
    double m_th_inten;
public:
    //Get
    int GetAlignmentNum()
    {
        return m_AlignmentNum;
    }
    int GetFileNum()
    {
        return m_FileNum;
    }
    // constructor
    FeatureMatrixReader(string input) : Reader(input)
    {
        m_th_rt = 50;
        m_th_mz = 0.1; // This does matter.// Why 0.5? // How to change this? by the user?
        m_th_inten = 1e9;// This does not matter
    }
	void Set_mz_threshold(double new_mz_threshold=0.5)
	{
		// 0.5 is the default threshold in jijie's code
		m_th_mz = new_mz_threshold;
	}
    void Set_rt_threshold(double new_rt_threshold=50)
    {
        m_th_rt = new_rt_threshold;
    }
    // destructor
    virtual ~FeatureMatrixReader()
    {
        if (m_featurematrix != NULL) delete m_featurematrix;
    }
    virtual FeatureMatrix *CreateFeatureMatrix(int AlignmentNum, int FileNUm)
    {
        return NULL;
    }
    virtual void ReadTwoStr(istream & ism) {}
    virtual void ValidateConsensus(FeatureMatrixReader * fmr)
    {

        for (int i = 0; i < fmr->GetAlignmentNum(); i++)
        {
            // for each line
            for (int j = 0; j < fmr->GetFileNum(); j++)
            {
                PeakFeature * consensuspf = fmr->GetPeakFeature(i, j);
                if (!isValid(consensuspf, j))
                {
                    consensuspf->SetMZIntenRT(0, 0, 0);
                }
            }
        }
    }
    virtual bool isValid(PeakFeature *pf, int FileID)
    {

        for (int i = 0; i < GetAlignmentNum(); i++)
        {
            PeakFeature * curpf = GetPeakFeature(i, FileID);

            if (fabs(curpf->GetMZ() - pf->GetMZ()) < m_th_mz && fabs(curpf->GetRT() - pf->GetRT()) < m_th_rt && fabs(curpf->GetInten() - pf->GetInten()) < m_th_inten)
            {
                return true;

            }
        }
        return false;
    }
    void CheckFeaturePair(FeatureMatrixReader * fmr, int & fmrFeaturePairSum, int &fmrFeaturePairChecked)
    {
        FeaturePair::m_th_mz = m_th_mz;
        FeaturePair::m_th_rt = m_th_rt;
        int found = 0;
        fmrFeaturePairChecked = 0, fmrFeaturePairSum = 0;

        for (int i = 0; i < fmr->GetAlignmentNum(); i++)
        {
            for (int j = 0; j < fmr->GetFileNum() - 1; j++)
            {
                for (int k = j + 1; k < fmr->GetFileNum(); k++)
                {
                    PeakFeature *fmrpfx = fmr->GetPeakFeature(i, j);
                    PeakFeature *fmrpfy = fmr->GetPeakFeature(i, k);

                    if (fmrpfx->isNonZero() && fmrpfy->isNonZero())
                        fmrFeaturePairSum += 1;
                    else
                    {
                        continue;
                    }

                    found = 0;
                    for (int m = 0; m < GetAlignmentNum(); m++)
                    {
                        PeakFeature *curpfx = GetPeakFeature(m, j);
                        PeakFeature *curpfy = GetPeakFeature(m, k);

                        if (FeaturePair(curpfx, curpfy) == FeaturePair(fmrpfx, fmrpfy))
                        {
                            found = 1;
                            break;
                        }
                    }
                    fmrFeaturePairChecked += found;
                    if (found == 0)
                    {
                        cout << *fmrpfx << " - " << *fmrpfy << endl;
                    }
                }
            }
        }
    }
    void parse()
    {
        string line;
        ifstream fin;
        fin.open(GetFileName(), ios::in);
        fin >> m_AlignmentNum >> m_FileNum;
        m_featurematrix = CreateFeatureMatrix(m_AlignmentNum, m_FileNum);
        ReadTwoStr(fin);
        for (int i = 0; i < m_FileNum; i++)
        {
            fin >> line;
            m_FileNames.push_back(line);

        }
        cout << "Start Reading " << GetFileName() << endl;

        fin >> *m_featurematrix;
        fin.close();
    }
    void print()
    {
        cout << *m_featurematrix << endl;
    }
    PeakFeature * GetPeakFeature(int i, int j)
    {
        return m_featurematrix->GetPeakFeature(i, j);
    }


    void CalcMassDiff()
    {
        bool output_flag = false;
        vector<double> massDiff;
        for (int i = 0; i < GetAlignmentNum(); ++i) {
            PeakFeature * x = GetPeakFeature(i,0);
            PeakFeature * y = GetPeakFeature(i,1);
            double massdiff = x->GetMZ() - y->GetMZ();
            if(fabs(massdiff)< 0.5)
            {
                massDiff.push_back(massdiff);
            }
            else if (output_flag)
            {
                cout << *x << endl;
                cout << *y << endl;
                cout << "--" << endl;
            }

        }
        double avg = calcmean(massDiff);
        double massSTD = calcSTD(massDiff);
        cout << "Mass Deviation:" << endl;
        cout << "Mean " << avg << endl;
        cout << "STD " << massSTD << endl;
    }
};

class ConsensusReader : public FeatureMatrixReader
{
public:
    // constructor
    ConsensusReader(string input) : FeatureMatrixReader(input) {	}
    // destructor
    ~ConsensusReader() {  }
private:
    // factory method
    FeatureMatrix * CreateFeatureMatrix(int AlignmentNum, int FileNum)
    {
        return new Consensus(AlignmentNum, FileNum);
    }


};

class groundtruthReader : public FeatureMatrixReader
{
public:
    groundtruthReader(string input) : FeatureMatrixReader(input)
    {
    }
    ~groundtruthReader() { }
    FeatureMatrix *CreateFeatureMatrix(int AlignmentNum, int FileNum)
    {
        return new GroundTruth(AlignmentNum, FileNum);
    }
    void ReadTwoStr(istream &ism)
    {
        string line;
        ism >> line >> line;
    }


    void CalculateRecallPrecision(ConsensusReader *consensus)
    {
        int checked = 0, sum = 0;
        CheckFeaturePair(consensus, sum, checked);
        ofstream fout;
        string recall_precision_file = consensus->GetFileName() + ".recall_precision.txt";
        fout.open(recall_precision_file.c_str(),ios::out | ios::app);
        //fout << recall_precision_file << endl;
        // we need get too.
        fout << getM_th_mz() << "," << getM_th_rt() << ",";

        //fout <<"threhold of mz and rt: " << getM_th_mz() << " " << getM_th_rt() << endl;
        fout << checked << "," << sum << "," << checked *1.0 / sum << ",";
        //fout << "Precision: " << checked << "/" << sum << " = " << checked *1.0 / sum << endl;
        cout << "Precision: " << checked << "/" << sum << " = " << checked *1.0 / sum << endl;
        consensus->CheckFeaturePair(this, sum, checked);
        fout << checked << "," << sum << "," << checked *1.0 / sum << endl;
        //fout << "Recall: " << checked << "/" << sum << " = " << checked *1.0 / sum << endl;
        cout << "Recall: " << checked << "/" << sum << " = " << checked*1.0 / sum << endl;
        fout.close();

    }

};

class pepxmlReader : public Reader
{
    vector<PSM> m_PSMs;
    map<string, vector<double> > GroundTruth;
    double threshold;
public:
    pepxmlReader(string input) : Reader(input)
    {
        if(GetFileName().find("ipro.pep.xml") == string::npos )
        {
            cout << "ATTENTION: Input file is not ipro.pep.xml format, process of this file may be incorrect!" << endl;
        }
    }
    ~pepxmlReader() {}
    void parse()
    {
        string line, preline = "";
        ifstream fin;
        fin.open(GetFileName(), ios::in);
        int i = 0;
        while (fin >> line)
        {
            while (line.find("<error_point") != string::npos)
            {
                // If find error point line
                fin >> line; //error =
                double error = GetVal(line);
                fin >> line; // min_prob=
                double min_prob = GetVal(line);
                fin >> line; // numcorr =
                double num_corr = GetVal(line);
                fin >> line; // numincorr=
                double num_incorr = GetVal(line);
                fin>> line; // next error_point
                if (fabs(error - 0.01) < 1e-4)
                {
                    threshold = min_prob;
                    cout << "threshold=" << threshold << " error=" << error << " num_corr=" << num_corr << " num_incorr=" << num_incorr<< endl;
                }
            }
            if (line.find("<spectrum_query") != string::npos)
            {
                i++;
                m_PSMs.push_back(PSM());
            }
            else if (line.find("spectrum=") != string::npos)
            {
                line = SplitByFirst(line, '\"')[1];
                line = SplitByFirst(line, '.')[0];
                m_PSMs[i - 1].spectrum = line;
            }
            else if (line.find("precursor_neutral_mass=") != string::npos)
            {
                line = SplitByFirst(line, '\"')[1];
                line = SplitByFirst(line, '\"')[0];
                m_PSMs[i - 1].precursor_neutral_mass = stod(line);

            }
            else if (line.find("calc_neutral_pep_mass=") != string::npos)
            {
                m_PSMs[i-1].calc_neutral_pep_mass = GetVal(line);
            }
            else if (line.find("assumed_charge=") != string::npos)
            {
                line = SplitByFirst(line, '\"')[1];
                line = SplitByFirst(line, '\"')[0];
                m_PSMs[i - 1].assumed_charge = stod(line);
                m_PSMs[i - 1].mz = m_PSMs[i - 1].precursor_neutral_mass / m_PSMs[i - 1].assumed_charge + 1;
            }
            else if (line.find("retention_time_sec=") != string::npos)// Original do noting!
            {
                line = SplitByFirst(line, '\"')[1];
                line = SplitByFirst(line, '\"')[0];
                m_PSMs[i - 1].retention_time_sec = stod(line);
            }
            else if (line.find("peptide=") != string::npos)// ToDo
            {
                // in case that we will find modified_peptide=
                if (line.find("peptide=") < 1)
                {
                    line = SplitByFirst(line, '\"')[1];
                    line = SplitByFirst(line, '\"')[0];
                    m_PSMs[i - 1].peptide = line + to_string(m_PSMs[i - 1].assumed_charge);

                }

            }
            else if (line.find("protein_descr=") != string::npos)// ToDo
            {
                //cout << line << endl;
                line = SplitByFirst(line, '\"')[1];
                //if (line == "Decoy")
                //cout << line << endl;
                if(line.find("Decoy") != string::npos)
                {
                    //cout <<"found:" <<  line << endl;
                    m_PSMs[i-1].isDecoy = true;
                }
                else
                {
                    m_PSMs[i-1].isDecoy = false;
                }
                //protein_descr
            }
            else if (line.find("start_scan=") != string::npos)// new entry
            {
                line = SplitByFirst(line, '\"')[1];
                line = SplitByFirst(line, '\"')[0];
                m_PSMs[i - 1].scan = line; // Get RT by scan? They are using something without rt
                // Let's use the correct RT

            }
            else if (line.find("probability=") != string::npos && preline.find("<peptideprophet_result") != string::npos)// ToDo
            {
                line = SplitByFirst(line, '\"')[1];
                line = SplitByFirst(line, '\"')[0];
                m_PSMs[i - 1].probability = stod(line);// ToDo, use threshold. i--
                if (m_PSMs[i - 1].probability < threshold && GetFileName().find("pep.xml") != string::npos && GetFileName().find("ipro.pep.xml")==string::npos)
                {
                    // Do nothing here
                    m_PSMs.pop_back();
                    i--;
                }
            }
            else if (line.find("probability=") != string::npos && preline.find("<interprophet_result") != string::npos)
            {
                //cout << "i-1" << i-1 << endl;
                m_PSMs[i-1].iprobability = GetVal(line);
                //if(m_PSMs[i-1].isDecoy == true && m_PSMs[i - 1].iprobability >= threshold) cout << i-1 << " Decoy: " << m_PSMs[i-1].iprobability << endl;
                if (m_PSMs[i - 1].iprobability < threshold)
                {
                    m_PSMs.pop_back();

                    i--;
                }

            }
            //cout << line << endl;
            preline = line;
        }
        cout << "querynum:" << m_PSMs.size() << endl;
        cout << "decoynum:" << CountDecoy() << endl;
        fin.close();
    }
    void Print()
    {
        for (int i = 0; i < m_PSMs.size(); i ++)
        {
            cout << m_PSMs[i];
        }
    }

    void CalcMassDiff()
    {
        vector<double> massdiff;
        for (int i = 0; i < m_PSMs.size() ; ++i) {
            //cout << m_PSMs[i].precursor_neutral_mass - m_PSMs[i].calc_neutral_pep_mass << endl;
            massdiff.push_back((m_PSMs[i].precursor_neutral_mass-m_PSMs[i].calc_neutral_pep_mass)/m_PSMs[i].assumed_charge);
        }
        double avg = calcmean(massdiff);
        double massSTD = calcSTD(massdiff);
        cout << "Mass Deviation:" << endl;
        cout << "mean " << avg << endl;
        cout << "STD " << massSTD << endl;
    }
    void OutputToTSV()
    {
        string filename = GetFileName() + ".tsv";
        ofstream fout;
        fout.open(filename.c_str(), ios::out);
        //std::setprecision(15);
        for (int i = 0; i < m_PSMs.size(); ++i) {
            fout << m_PSMs[i].spectrum << "\t" << m_PSMs[i].scan << "\t" << std::setprecision(5)
                 << m_PSMs[i].precursor_neutral_mass << "\t" << m_PSMs[i].calc_neutral_pep_mass << "\t"
                 << m_PSMs[i].precursor_neutral_mass - m_PSMs[i].calc_neutral_pep_mass << "\t"
                 << m_PSMs[i].assumed_charge << "\t" << m_PSMs[i].mz << "\t"
                 << m_PSMs[i].peptide << "\t" << m_PSMs[i].probability << "\t"
                  << m_PSMs[i].iprobability << std::fixed << "\t" << m_PSMs[i].isDecoy << "\t" << endl;
        }
        fout.close();
    }

private:
    int CountDecoy()
    {
        int decoynum = 0;
        for (int i = 0; i < m_PSMs.size(); i ++)
        {
            if (m_PSMs[i].isDecoy == true)
            {
                decoynum ++;
            }
        }
        return decoynum;
    }
    double GetVal(string line)
    {
        line = SplitByFirst(line, '\"')[1];
        line = SplitByFirst(line, '\"')[0];
        return stod(line);
    }


    void groupmzrt(double &mzmean, double &rtmean, vector<double> mz, vector<double> rt)
    {
        //cout << "start group" << endl;
        if (0 == mz.size())
        {
            mzmean = 0;
            rtmean = 0;
            return;
        }
        else if (1 == mz.size())
        {
            mzmean = mz[0];
            rtmean = rt[0];
            return;
        }
        rtmean = calcmean(rt);
        mzmean = calcmean(mz);
        double rtSTD = calcSTD(rt);

        bool valid = true;
        if (rtSTD > 100) valid = false;
        for (unsigned int i = 0; i < rt.size(); i++)
        {
            if (rt[i] - rtmean > 2 * rtSTD)
            {
                valid = false;
            }
        }

        if (!valid)
        {
            mzmean = 0;
            rtmean = 0;
        }
        //cout << "end group" << endl;
    }
public:
    void GetGroundTruth(vector<string> filenames)
    {
        cout << "Looking for gound truth" << endl;
        sort(m_PSMs.begin(), m_PSMs.end(), PSM::asdsortbyPeptide);
        for (unsigned int i = 0; i < m_PSMs.size(); i++)
        {
            // for a peptide
            string peptide_charge = m_PSMs[i].peptide;
            //cout << peptide_charge << endl;
            vector<double> rt_mz_list;
            int k = i;
            // for a filename
            for (unsigned int j = 0; j < filenames.size(); j++)
            {
                string filename = filenames[j] + "_Q1";
                // for each PSM with correct peptide and correct filename
                k = i;
                vector<double> mz, rt;
                while (k < m_PSMs.size() && 0 == m_PSMs[k].peptide.compare(peptide_charge))
                {
                    if (0 == m_PSMs[k].spectrum.compare(filename))// if file name correct
                    {
                        // add mz, rt
                        mz.push_back(m_PSMs[k].mz);//
                        rt.push_back(m_PSMs[k].retention_time_sec);
                    }
                    k++;
                }
                double mzmean, rtmean;
                groupmzrt(mzmean, rtmean, mz, rt);
                rt_mz_list.push_back(rtmean);

                // Different way of update mz
                // use the peptide mass if it is not zero
                if(fabs(mzmean-m_PSMs[i].GetPepMZ()) < 3) // GetPrecursorMz(): modyfy here today. there is no difference
                {
                    rt_mz_list.push_back(m_PSMs[i].GetPepMZ());//.GetPepMZ()); // modify here today: wu long 06/07/2015
                }
                else
                {
                    rt_mz_list.push_back(mzmean);
                }
                // use the mzmean anyway, jiji's code
                //rt_mz_list.push_back(mzmean);

            }
            int countzero = 0;
            for (unsigned int idx = 0; idx < rt_mz_list.size(); idx += 2)
            {
                if (fabs(rt_mz_list[idx]) < 0.001)
                    countzero++;
            }
            if (countzero <= filenames.size() * 0.1)
                GroundTruth[peptide_charge] = rt_mz_list;
            i = k - 1;


        }
    }
    void OutputGroundTruth(string GroundTruthPath, vector<string> filenames)
    {
        cout << "Start output ground truth" << endl;
        ofstream fout;
        fout.open(GroundTruthPath, ios::out);
        fout << GroundTruth.size() << " " << filenames.size() << endl;
        fout << "Peptide Sequence";
        for (unsigned int i = 0; i < filenames.size(); i++)
        {
            fout << " " << filenames[i];
        }
        fout << endl;

        // output groundtruth
        for (map<string, vector<double> >::iterator it = GroundTruth.begin(); it != GroundTruth.end(); it++)
        {
            fout << it->first;
            for (unsigned int i = 0; i < it->second.size(); i++)
            {
                fout << " " << it->second[i];
            }
            fout << endl;
        }
        fout.close();
        cout << "end output ground truth" << endl;
    }

};


class mzXMLListReader : public Reader
{
    vector<string> m_shortFileName;
    map<string, int> m_File2ID;
public:
    mzXMLListReader(string input) : Reader(input) {}
    ~mzXMLListReader() {}
    void parse()
    {
        string filename;
        ifstream fin;
        fin.open(GetFileName(), ios::in);
        int i = 0;
        while (fin >> filename)
        {
            int found = filename.find_last_of(PATH_SLASH);
            if (found != string::npos)
            {
                filename = filename.substr(found + 1);
            }
            found = filename.find_last_of('.');
            if (found != string::npos)
            {
                filename = filename.substr(0, found);
            }
            m_shortFileName.push_back(filename);
            m_File2ID[filename] = i + 1;
            cout << filename << endl;
        }
    }
    // get
    vector<string> getshortnames()
    {
        return m_shortFileName;
    }
};

class Scan2RTReader : public Reader
{
    map<string, double> m_Scan2RT;
public:
    // constructor
    Scan2RTReader(string input) :Reader(input)	{}
    // destructor
    ~Scan2RTReader() {}
    void parse()
    {
        int sum = 0;
        string label, scan;
        double rt;
        ifstream fin;
        fin.open(GetFileName(), ios::in);
        while (fin >> label >> scan >> rt)
        {
            sum++;
            m_Scan2RT[label + scan] = rt;
        }
        fin.close();
        cout << "Scan2RT pair number " << sum << endl;
    }
};
}

#endif // !_FILEREADER_H_

