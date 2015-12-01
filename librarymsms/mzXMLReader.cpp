#include "mzXMLReader.h"

mzXMLReader::mzXMLReader(string mslevel) {
    MaxScanDiff = 1000; // Every two scan should calculate the Dot Product.
    //cout << "[Info] Dot product bandwidth = " << MaxScanDiff << endl;
    SLASH = "/";
    // or "\\"
    MSLEVEL = mslevel;
}

vector<PeakList*> mzXMLReader::ReadmzXMLToPeakLists(mzXMLFilename f) {
    cout << "[Info] Loading " << f << endl;
    //f = "D:\\Projects\\mzXMLReader\\GZ_1_1.mzXML";
    //f = "D:\\Swath\\DiaumpireData\\UPS\\LongSwath_UPS1_1ug_rep1_SWATH2.mzXML";
    cRamp* cramp = new cRamp(f.c_str());
    if (!cramp->OK()) {
        cout << "Cannot open file \"" << f << "\". File skipped." << endl;
        delete (cramp);
    }
    rampRunInfo* runInfo = cramp->getRunInfo();

    if (!runInfo) {
        cout << "Cannot open file \"" << f << "\". File skipped." << endl;
        exit(0);

    }
    rampInstrumentInfo* instr = cramp->getInstrumentInfo();
    vector<PeakList*> vpl;
    for(int k = 1; k <= cramp->getLastScan(); k ++)
    {
        rampScanInfo* scanInfo = cramp->getScanHeaderInfo(k);
        PeakList * p = new PeakList();
        // If the scan is NULL, skip this scan
        if (!scanInfo) {continue; }
        p->setRTinSeconds(scanInfo->getRetentionTimeSeconds());

        // If the msLevel is not equal to MSLEVEL, skip this scan
        if(scanInfo->m_data.msLevel!=atoi(MSLEVEL.c_str()))  continue;

        rampPeakList* peaks = cramp->getPeakList(scanInfo->m_data.acquisitionNum);



        //cout << scanInfo->m_data.precursorScanNum<< " " << scanInfo->m_data.msLevel << endl;
        if (!peaks)
        {
            continue;
        }
        double sumint = 0;

        for(int j = 0; j < peaks->getPeakCount(); j ++)
        {
            double mz = peaks->getPeak(j)->mz;
            float intensity = (float)(peaks->getPeak(j)->intensity);
            p->InsertPeak(mz,intensity);
        }
        vpl.push_back(p);
    }
    delete cramp;
    return vpl;
}

void mzXMLReader::CalculateDotProduct(vector<PeakList*> &a, vector<PeakList*> &b, int start, int end, int thread_ID,double * res) {
    end = end < a.size()? end: a.size();
    //cout << start << " --> " << end << endl;
    int LastPercent = 0;
    long bSize = b.size();
    long long  k = start*bSize;
    for (int i = start; i < end ; i ++)
    {
        // Debug: only when start == 0
        //if(start != 0)
        // break;
        for (int j = 0; j < b.size(); ++j)
        {
            double dot = 0;
            if(i - MaxScanDiff > j || i + MaxScanDiff <j)
            {
                dot = 0;
            }
            else if(b[j] == NULL || a[i] == NULL){
                dot = 0;

            }
            else
            {
                dot = a[i]->CalcDotProduct(*(b[j]));
            }
            res[k] = dot;
            // cout << dot << endl;
            k++;
        }
        int percent = (100*(i+1-start))/(end-start);
        if(percent%50==0 && LastPercent != percent)
        {
            LastPercent = percent;
            cout <<"." << flush;
        }

    }

}

void mzXMLReader::CreateBinningList(vector<PeakList*> &a) {
    for (int i = 0; i < a.size(); ++i) {
        if (a[i]!=NULL)
        {
            BinningPeakList * x = a[i]->CreateBinningPeaks();
        }

    }
}

void mzXMLReader::MultiThreadDotProduct(vector<PeakList*> &a, vector<PeakList*> &b, string outputbinaryfile, int threadNum) {
    cout << "[Info] Calculating Dot Product.." << endl;

    long long int nSize;
    double *matrix;
    ApplyMemoryforMatrix(a.size(), b.size(), nSize, matrix);

    if(matrix == NULL)
    {
        cout << "[Error] Fail to apply memory" << endl;
        exit(0);
    }

    CreateBinningList(a);
    CreateBinningList(b);
    CalcDotProduct(a, b, matrix, threadNum);

    Matrix m(a.size(),b.size(),matrix);

    m.outputBinary(outputbinaryfile);



    //outputMatrix(outputbinaryfile, nSize, matrix);
    delete [] matrix;

    //exit(0);
}

void mzXMLReader::CalcDotProduct(vector<PeakList *> &a, const vector<PeakList *> &b, double *res, int threadNum) {
    int MinSize = a.size()/threadNum+1;
    vector<thread> tasks;
    int k = 0;
    cout << "[Info] Start dot product" << endl;
    for (int i = 0; i < a.size(); i+=MinSize)
    {
        tasks.push_back(std::thread(bind(&mzXMLReader::CalculateDotProduct,mzXMLReader(MSLEVEL),a,b,i,i+MinSize,k++,res)));
    }
    for (int i = 0; i < threadNum; ++i)
    {
        tasks[i].join();
    }
    cout << endl;
    cout << "[Info] End dot product" << endl;
}

void mzXMLReader::SingleThreadDotProduct(vector<PeakList*> &a, vector<PeakList*> &b, string outputfile) {
    //ofstream fout;
    FILE * pfile = fopen(outputfile.c_str(),"w");
    //fout.open("D:\\Swath\\DiaumpireData\\UPS\\LongSwath_UPS1_1ug_rep1_2.mzXML.dot.tsv",ios::out);
    for(int i = 0; i < a.size() ; i ++)
    {
        for(int j = 0; j < b.size() ; j ++)
        {
            double dot = 0;
            if(i - MaxScanDiff > j || i + MaxScanDiff <j)
                dot = 0;
            else
                dot = a[i]->CalcDotProduct(*(b[j]));
            fprintf(pfile,"%.2lf,",dot);
        }
        fprintf(pfile,"\n");
        cout << "Progress: " << i << " /" << a.size() << "\r" << flush;
        //fout << endl;
    }
    //fout.close();
    fclose(pfile);
}

void mzXMLReader::ApplyMemoryforMatrix(const long aSize, const long bSize, long long int &nSize, double *&res) {
    nSize= aSize * bSize;
    res= new double [nSize];
}