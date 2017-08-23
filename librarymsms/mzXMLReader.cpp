#include "mzXMLReader.h"
#include "math.h"

mzXMLReader::mzXMLReader(string mslevel, int max_scan_th) {
    MaxScanDiff = max_scan_th; // Every two scan should calculate the Dot Product.
    //cout << "[Info] Dot product bandwidth = " << MaxScanDiff << endl;
    SLASH = "/";
    // or "\\"
    MSLEVEL = mslevel;

    m_SimCalculator = new DotProduct();

}

vector<PeakList *> mzXMLReader::ReadmzXMLToPeakLists(mzXMLFilename f) {
    cout << "[Info] Loading " << f << endl;
    int emptySpecNum = 0;
    //f = "D:\\Projects\\mzXMLReader\\GZ_1_1.mzXML";
    //f = "D:\\Swath\\DiaumpireData\\UPS\\LongSwath_UPS1_1ug_rep1_SWATH2.mzXML";
    cRamp *cramp = new cRamp(f.c_str());
    if (!cramp->OK()) {
        cout << "Cannot open file \"" << f << "\". File skipped." << endl;
        delete (cramp);
    }
    rampRunInfo *runInfo = cramp->getRunInfo();

    if (!runInfo) {
        cout << "Cannot open file \"" << f << "\". File skipped." << endl;
        exit(0);

    }
    rampInstrumentInfo *instr = cramp->getInstrumentInfo();
    vector<PeakList *> vpl;
    for (int k = 1; k <= cramp->getLastScan(); k++) {
        rampScanInfo *scanInfo = cramp->getScanHeaderInfo(k);
        PeakList *p = new PeakList();
        // If the scan is NULL, skip this scan
        if (!scanInfo) { continue; }
        p->setRTinSeconds(scanInfo->getRetentionTimeSeconds());

        // If the msLevel is not equal to MSLEVEL, skip this scan
        if (scanInfo->m_data.msLevel != atoi(MSLEVEL.c_str())) continue;

        rampPeakList *peaks = cramp->getPeakList(scanInfo->m_data.acquisitionNum);



        //cout << scanInfo->m_data.precursorScanNum<< " " << scanInfo->m_data.msLevel << endl;
        if (!peaks) {
            vpl.push_back(NULL);
            emptySpecNum++;
            continue;
        }
        double sumint = 0;

        for (int j = 0; j < peaks->getPeakCount(); j++) {
            double mz = peaks->getPeak(j)->mz;
            float intensity = (float) (peaks->getPeak(j)->intensity);
            p->InsertPeak(mz, intensity);
        }
        vpl.push_back(p);
    }
    delete cramp;
    cout << "[Info] Empty spectra: " << emptySpecNum << endl;
    return vpl;
}

void mzXMLReader::CalculateDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, int start, int end, int thread_ID,
                                      double *res) {
    end = end < a.size() ? end : a.size();
    //cout << start << " --> " << end << endl;
    int LastPercent = 0;
    long bSize = b.size();
    long long k = start * bSize;
    for (int i = start; i < end; i++) {
        // Debug: only when start == 0
        //if(start != 0)
        // break;
        for (int j = 0; j < b.size(); ++j) {
            double dot = 0;
            if (i - MaxScanDiff > j || i + MaxScanDiff < j) {
                dot = 0;
            }
            else if (b[j] == NULL || a[i] == NULL) {
                dot = 0;

            }
            else {
                dot = a[i]->CalcDotProduct(*(b[j]));
            }
            res[k] = dot;
            // cout << dot << endl;
            k++;
        }
        int percent = (100 * (i + 1 - start)) / (end - start);
        if (percent % 50 == 0 && LastPercent != percent) {
            LastPercent = percent;
            cout << "." << flush;
        }

    }

}

void mzXMLReader::CreateBinningList(vector<PeakList *> &a) {

    for (int i = 0; i < a.size(); ++i) {
        if (a[i] != NULL) {
            BinningPeakList *x = a[i]->CreateBinningPeaks();
            // We could filter it here
            for (int j = 0; j < m_PeakFilters.size(); ++j) {
                x = x->filterPeaks(m_PeakFilters[j]);
            }
            if (m_SimCalculator == NULL) {
                cout << "[Error] Invalid simMetric. please call " << endl;

            }

            double norm = m_SimCalculator->calcNorm(*x);


        }

    }

}

//<<<<<<< HEAD
//void mzXMLReader::MultiThreadDotProduct(vector<PeakList*> &a, vector<PeakList*> &b, string outputbinaryfile, int threadNum) {
//    cout << "[Info] Calculating Dot Product.." << endl;
//    SimpleTimer st("calculate dot product of cycle " + to_string(1));
//=======
void
mzXMLReader::MultiThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputbinaryfile, int threadNum,
                                   int cycle) {
    //  cout << "[Info] Calculating Dot Product.." << endl;
    // SimpleTimer st("calculate dot product of cycle " + to_string(1));
//>>>>>>> ms2_alignment

    long long int nSize;
    double *matrix;
    ApplyMemoryforMatrix(a.size(), b.size(), nSize, matrix);

    if (matrix == NULL) {
        cout << "[Error] Fail to apply memory" << endl;
        throw "Not enough memory for dot product matrix.";
    }
//    m_SimCalculator = new DotProduct();
//    m_SimCalculator = new PNorm(1);
//    m_SimCalculator = new SPC();
    // ToDo: this declaration could be earlier
    cout << "[Info] Filters used:";
    for (int k = 0; k < m_PeakFilters.size(); ++k) {
        cout << " " << m_PeakFilters[k]->getType();
    }
    cout << endl;
    CreateBinningList(a);
    CreateBinningList(b);

//    CalcNorm(a, b, matrix, threadNum);
    SingleThreadDotProduct(a, b, matrix);

    Matrix m(a.size(), b.size(), matrix);

    m.outputBinary(outputbinaryfile);



    //outputMatrix(outputbinaryfile, nSize, matrix);
    delete[] matrix;
//    delete m_SimCalculator;
//    m_SimCalculator = NULL;

    //exit(0);
}

void mzXMLReader::CalcDotProduct(vector<PeakList *> &a, const vector<PeakList *> &b, double *res, int threadNum) {
    int MinSize = a.size() / threadNum + 1;
    vector<thread> tasks;
    int k = 0;
    cout << "[Info] Start calculating dot product by multi-thread..." << endl;
//    int max_scan_diff = 1000;
    SimpleTimer st("MultiThread DotProduct");
//    cout << "[Info] Maximal scan deviation: " << max_scan_diff << endl;
    for (int i = 0; i < a.size(); i += MinSize) {
        tasks.push_back(std::thread(
                ::std::bind(&mzXMLReader::CalculateDotProduct, mzXMLReader(MSLEVEL), a, b, i, i + MinSize, k++, res)));
    }
    for (int i = 0; i < threadNum; ++i) {
        tasks[i].join();
    }
    cout << endl;
    cout << "[Info] Finish dot-product calculating" << endl;
}

void mzXMLReader::SingleThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, double *res) {


    SimpleTimer st("OpenMP DotProduct");
    // This is how we do this work!

//    m_SimCalculator = new PNorm();
//>>>>>>> ms2_alignment
#pragma omp parallel for  schedule(dynamic) collapse(2)
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < b.size(); ++j) {

            double similarity = 0;
            double similarity_another = 0;
            // dot product

            if (b[j] == NULL || a[i] == NULL) {
                similarity = 0;

            } else {

//                similarity = a[i]->CalcDotProduct(*(b[j]));
                // try another way to calculate simlarity
                similarity = m_SimCalculator->calcDistance(*a[i], *b[j]);

//                cout << similarity << " " << similarity_another << endl;
//                if (fabs(similarity-similarity_another) > EPSILON)
//                    throw "Error inconsistence similarity";
            }
            res[i * b.size() + j] = similarity;

        }
    }


}

void mzXMLReader::SingleThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputfile) {
    //ofstream fout;
    FILE *pfile = fopen(outputfile.c_str(), "w");
    //fout.open("D:\\Swath\\DiaumpireData\\UPS\\LongSwath_UPS1_1ug_rep1_2.mzXML.dot.tsv",ios::out);
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < b.size(); j++) {
            double dot = 0;
            if (i - MaxScanDiff > j || i + MaxScanDiff < j)
                dot = 0;
            else
                dot = a[i]->CalcDotProduct(*(b[j]));
            fprintf(pfile, "%.2lf,", dot);
        }
        fprintf(pfile, "\n");
        cout << "Progress: " << i << " /" << a.size() << "\r" << flush;
        //fout << endl;
    }
    //fout.close();
    fclose(pfile);
}

void mzXMLReader::ApplyMemoryforMatrix(const long aSize, const long bSize, long long int &nSize, double *&res) {
    nSize = aSize * bSize;
    res = new double[nSize];
}

mzXMLReader::~mzXMLReader() {
    if (m_SimCalculator != NULL) {
        delete m_SimCalculator;
        m_SimCalculator = NULL;
    }

    for (int i = 0; i < m_PeakFilters.size(); ++i) {
        if (m_PeakFilters[i] != NULL) {
            delete m_PeakFilters[i];
            m_PeakFilters[i] = NULL;
        }

    }
}


void mzXMLReader::RandMultiThreadDotProduct(vector<PeakList *> &a, vector<PeakList *> &b, string outputbinaryfile,
                                            int threadNum, int cycle) {
    long long int nSize;
    double *matrix;
    int num_RT_a = a.size() / (cycle - 1);
    int num_RT_b = b.size() / (cycle - 1);
    ApplyMemoryforMatrix(num_RT_a, num_RT_b, nSize, matrix);

    if (matrix == NULL) {
        cout << "[Error] Fail to apply memory" << endl;
        throw "Not enough memory for dot product matrix.";
    }


    CreateBinningList(a);
    CreateBinningList(b);

    vector<vector<double> > product_a;
    vector<vector<double> > product_b;

    vector<double> random_vector;
    default_random_engine generator;
    normal_distribution<double> distribution(0, 1);
    for (int i = 0; i < a[0]->getM_binningPeakList()->GetBinNum(); i++) {
        double number = distribution(generator);
        random_vector.push_back(number);
    }


    RandomVectorMultiplication(a, product_a, random_vector, num_RT_a, cycle);
    RandomVectorMultiplication(b, product_b, random_vector, num_RT_b, cycle);

//    CalcNorm(a, b, matrix, threadNum);
    SingleThreadDotProduct(product_a, product_b, matrix);

    Matrix m(num_RT_a, num_RT_b, matrix);

    m.outputBinary(outputbinaryfile);


    //outputMatrix(outputbinaryfile, nSize, matrix);
    delete[] matrix;
//    delete m_SimCalculator;
//    m_SimCalculator = NULL;

    //exit(0);
}

//calculate the product vector of each MS2 matrices
void mzXMLReader::RandomVectorMultiplication(vector<PeakList *> &a, vector<vector<double> > &product,
                                             vector<double> random_vector, int num_RT, int cycle) {

//#pragma omp parallel for  schedule(dynamic) collapse(2)
    for (int i = 0; i < num_RT; ++i) {
        vector<double> random_product;
        for (int j = 0; j < cycle - 1; ++j) {
            double p = 0;
            if (a[i * (cycle - 1) + j] != NULL) {
                vector<double> ms2spectrum = a[i * (cycle - 1) + j]->getM_binningPeakList()->GetIntensityList();
                for (int k = 0; k < random_vector.size(); ++k) {
                    p += ms2spectrum[k] * random_vector[k];
                }
            }
            random_product.push_back(p);
        }
        product.push_back(random_product);
    }
}

//overloaded function for vector of vector instead of vector of peaklist
void mzXMLReader::SingleThreadDotProduct(vector<vector<double> > a, vector<vector<double> > b, double *res) {


    SimpleTimer st("OpenMP DotProduct");
    // This is how we do this work!


//#pragma omp parallel for  schedule(dynamic) collapse(2)
    for (int i = 0; i < a.size(); ++i) {
        for (int j = 0; j < b.size(); ++j) {
            double similarity = 0;
            // dot product
            similarity = m_SimCalculator->calcDistance(a[i], b[j]);

            res[i * b.size() + j] = similarity;
        }
    }
}

vector<double> normalize(vector<vector<double> > a) {
    vector<double> length_vector;
    //RT dimension
    for (int i = 0; i < a.size(); i++) {
        double length = 0;
        for (int j = 0; j < a[i].size(); j++) {
            length += pow(a[i][j], 2);
        }
        length = sqrt(length);
        length_vector.push_back(length);
    }
    return length_vector;
}
