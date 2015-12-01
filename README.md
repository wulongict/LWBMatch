---------------------------------------------------------

# User Manual for LWBMatch

A profile-based and feature-based hybrid algorithm for retention time alignment of SWATH MS data

Date: Dec. 1, 2015  
Authors: Long Wu, Henry Lam

If you have any question about this manual, feel free to ask lwuai @ connect [dot] ust [dot] hk for help. 

----------------------------------------------------------


## Three tools are provided in this package.

1. Retention time alignment tool: **lwbmatch**
2. Generte ground truth from identification results: **groundtruth**
3. Recall and precision calculator: **calcrecallprecision**


## Build / Compile from source code
### Pre-requirment, 
- Boost (Version >= 1.52)
- CMake (Version >= 3.2)

### Commandlines for build:

```
cd /path/to/sourcecode
cmake .
make
```
Then we get **_lwbmatch, groundtruth_** and **_calcrecallprecision_** in /path/to/source/src, /path/to/source/groundtruth  and /path/to/source/calcrecallprecision.


----------------------------------------------

## Data Analysis: Run Test data set step by step

### 1. Download and unzip dataset: 
```
wget http://nesvilab.org/tsouc/DIA_Umpire_SampleDataUPS.zip
unzip DIA_Umpire_SampleDataUPS.zip
```
### 2. Run featurefinder
#### 2.1 convert each mzXML to mzML
```
FileConverter -threads 18 -in_type mzXML -in  [inputmzXML]  -out [outputname].tmp.mzML -out_type mzML
FileFilter -threads 18 -in_type mzML -in [outputname].tmp.mzML -out [outputname].mzML -out_type mzML
```
#### 2.2 Extract featureXML
```
FeatureFinderCentroided -threads 18 -in [outputname].mzML -out [outputname].featureXML
```
### 3. Run lwbmatch
#### 3.1 Create input featureXMLList
FileName: test.featureXML  
Write the two lines into this file, replace the /path/to/outputname with your own filename 
```
/path/to/outputname1.featureXML
/path/to/outputname2.featureXML
```
#### 3.2 Run lwbmatch with DTW
```
./lwbmatch -l /path/to/test.featureXML -w 2
```
#### 3.3 Run lwbmatch with LOWESS (old version)
```
./lwbmatch -l /path/to/test.featureXML -w 1
```
### 4. Create ground truth and calculate recall and precision
#### 4.1 Run DIA-Umpire, get [inputfile_Q1.mgf]
```
java -jar -Xmx16G /path/to/DIA_Umpire_SE.jar inputmzXML  [diaumpire_configurefile]
```
There will be three mzXML with Q1 ~ Q3 in name, we use one with Q1. 

**Instructions for downloading DIA-Umpire and its default parameter [diaumpire_configurefile]** 

- Find `diaumpire_configurefile` from [this page](http://sourceforge.net/projects/diaumpire/files/Parameter%20files/)  
	Fileanme: diaumpire.se_params

-  Find DiA-Umpire from [this page](http://sourceforge.net/projects/diaumpire/files/JAR%20executables/)  
	We use: DIA-Umpire_v1_284.zip

- Unzip DIA-Umpire_v1_284.zip and use DIA_Umpire_SE.jar

#### 4.2 Convert mgf to mzXML, [inputfile_Q1.mgf] -> [inputfile_Q1.mzXML]  (Proteowizard required)
```
msconvert --mzXML -o [outputpath] [inputfile_Q1.mgf]
```
#### 4.3 Run Comet, search [inputfile1_Q1.mzXML] and [inputfile2_Q1.mzXML] together
```
/path/to/comet -P[parameterfile] [inputfile1_Q1.mzXML] [inputfile2_Q1.mzXML]
```

- Find Comet from [this page](http://sourceforge.net/projects/comet-ms/files/)  
  
	We use comet version 2013.01 rev. 0

#### 4.4 Filter by xinteract , to get [*.ipro.pep.xml] (TPP required)
input: search result file [\*.pep.xml]  
output: [\*.ipro.pep.xml]  
For more information, please refer to [this page](http://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20(TPP)/)

#### 4.5 Create ground truth from [.ipro.pep.xml]
```
/path/to/groundtruth [filename:2] [mzXMLListtxtName_contains_two_mzXML_filename]  [ipropepxml]  [threshold_for_filter_default_0.9]  [groundtruth_filename]  [dummy_parameter:1]
```
#### 4.6 Calculate recall and precision
```
/path/to/calcrecallprecision [groundtruth_filename]  [RTAlignResult_with_extension_resu]  [mz_tolerance:0.1] [rt_tolerance: 1 ~ 25]
```

------------------------------------------------
