---------------------------------------------------------

# User Manual for feature-based label-free quantification workflow

Update: Sep. 15, 2016 
Authors: Long Wu, Yak Chak Li, Henry Lam

If you have any question about this manual, feel free to ask lwuai @ connect [dot] ust [dot] hk for help.

----------------------------------------------------------


## Requirements
### System requirements
- Operating System: Ubuntu 14.04  
We will transplant the workflow for Windows platform later.
- Memory: 32 GB or above  
This requirement is to make sure we can run large mzXML (> 4GB) files successfully.


### Python and several packages used
- Python 2.7.* (Version 2.7.6 or higher, but not Version 3.*) 
- Scipy (Tested Version 0.16.1)
- Numpy (Tested Version 1.10.1)
- Matplotlib (Tested Version 1.5.0)

### Third party software tools
- DIA-Umpire [Download](http://sourceforge.net/projects/diaumpire/files/Parameter%20files/)  
To extract pseudo-MS2 spectra from SWATH data. The parameter should also be provided.


- Comet [Download](http://sourceforge.net/projects/comet-ms/files/)  
To search pseudo-MS2 spectra against a sequence database.

- Xinteract (from TPP) [Download](http://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20(TPP)/)  
To control FDR and report the identified peptides/proteins above a FDR threshold


- OpenMS [Download](https://github.com/OpenMS/OpenMS)  
To extract the features from given SWATH MS data file. The features are used for retention time alignment and quantification. There are three tools used in our workflow including FeatureFinderCentroided, FileConverter, and FileFilter



## How to use the workflow
### Configuration
The first step to run the workflow is to set up several parameters. The parameter file is RTQuant.ini, which is in the **quantifyscripts** folder. The content of this file is as following:
```
[basic settings]
tpppath = /path/to/tpp/bin/
msconvertpath = /path/to/msconvert/bin/
extname = 
cometpath = /path/to/comet/binary/
recallprecision_binpath = /path/to/lwbmatch/bin/
binpath = /path/to/lwbmatch/bin/ 
lwbmatchpath = /path/to/lwbmatch/bin/
openmspath = /path/to/openms/bin/
mkdecoydbperl = /path+fullname/to/decoy-creating/script/creatDecoy.pl
diaumpirepath = /path/to/DIA-Umpire/bin/
systemname = Linux
cometname = filename-of-comet
```
There are several parameters in the above example are set as **/path/to/lwbmatch/bin/**. This path is /path/to/source/bin/, which contains **lwbmatch**, **groundtruth**, and **recallprecision**. Later we will consider to remove those redundant path settings.


### Run the workflow
```
python /path/to/quantifyscripts/RTQuant.py input1.mzXML input2.mzXML
```


## The workflow contains several different parts.
1. LWBMatch.py  
Retention time alignment: Scan- and Feature- based retention time alignment.

2. DIAumpireFlow.py  
Purpose: For get the identified PSMs and get the ground truth map.  
2.1 Pseudo-MS2 file extraction;   
2.2 Using Comet to search pseudo-MS2 spectra  
2.3 Using Xinteract to control FDR of identified PSMs.

FDR: false discovery rate  
PSM: peptide spectrum match  


----------------------------------------------