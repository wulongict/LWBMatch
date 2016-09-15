#!/usr/bin/python
# Python 2.7 Required
# Comet search. This script can not be launched by itself.
# Use another script to call this one.
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------


"""Run comet search, with 0.3Th tolerance?"""

from PipelineConfig import *
import glob


# wp = PipelineConfig()


class comet:
    def __init__(self, wfc, wp):
        """
        :type wfc: WorkflowController
        """
        print '[Info] Running comet...'
        self.wp = wp
        self.wfc = wfc
        if wp.extname == '':
            # on linux system ToDo:
            self.cometbinary = wp.cometname  # ToDo: I changed my code back to the old comet
            # self.cometbinary = "comet.2015012.linux.exe"
            print "*" * 50
            print "[Info] We are using new comet now"
            print "*" * 50
        else:
            self.cometbinary = 'comet.2015012.win64.exe'

    def convert_mgf_to_mzxml(self, MGFList):
        mzXMLList = []
        for eachMgf in MGFList:
            tmpfile = glob.glob(eachMgf[0:-4] + '.mzXML')
            mzXMLPath = os.path.split(eachMgf)[0]
            mzXMLList.append(eachMgf[0:-4] + '.mzXML')
            if len(tmpfile) == 0 or self.wfc.run_ms_convert:
                if self.wp.extname == '':  # running on linux
                    cmdline = "\"" + os.path.join(self.wp.msconvertpath,
                                                  'msconvert' + self.wp.extname) + "\"" + ' ' + eachMgf + ' --mzXML ' + ' -o ' + mzXMLPath
                else:
                    cmdline = "\"" + os.path.join(self.wp.msconvertpath,
                                                  'msconvert' + self.wp.extname) + "\"" + ' ' + eachMgf + ' --mzXML --outfile ' + eachMgf[
                                                                                                                                  0:-4] + '.mzXML' + ' -o ' + mzXMLPath
                os.system(cmdline)
        # end
        return mzXMLList

    def Run(self, MGFspectra, database):
        print '[Info] Running comet...'
        mzXMLList = self.convert_mgf_to_mzxml(MGFspectra)

        cometConfigFile = self.write_comet_config_file(database)
        pepxmllist = []
        cometCMD = os.path.join(self.wp.cometpath, self.cometbinary) + ' -P' + cometConfigFile + ' '
        for i in range(len(mzXMLList)):
            pepxmllist.append(mzXMLList[i][:-5] + 'pep.xml')
            print pepxmllist
            if len(glob.glob(pepxmllist[-1])) == 0 or self.wfc.run_comet:
                print '[Info] >' + cometCMD + ' ' + mzXMLList[i]
                os.system(cometCMD + ' ' + mzXMLList[i])
        # end

        return pepxmllist

    def write_comet_config_file(self, database):
        # read params from template file
        cometparamfile = "comet.params"
        if self.cometbinary == "comet.2015012.linux.exe":
            # linux and windows are using the same name.
            cometparamfile = "comet.params.openswath"
        configtemplate = os.path.join(self.wp.binpath, cometparamfile)
        fid = open(configtemplate, 'r')
        lines = fid.readlines()
        fid.close()

        print '[Info] Comet parameter template...'
        print ''.join(lines)

        # output params file
        cometconfigfile = os.path.join(os.path.split(database)[0], "comet.params.openswath")
        print '[Info] Start writing new parameter file to ' + cometparamfile

        fid = open(cometconfigfile, 'w')
        for i in range(len(lines)):
            # print lines[i]
            if 'database_name' in lines[i]:
                fid.write('database_name = %s\n' % database)
                print '[Info] database = ' + database
            elif 'decoy_search' in lines[i]:
                fid.write('decoy_search = 0\n')
            # ToDo: In the future, we need to change the search tolerance.
            # elif 'peptide_mass_tolerance' in lines[i]:
            #     fid.write('peptide_mass_tolerance = %d\n' %wp.search_parent_ion_tolerance_Th)
            else:
                fid.write(lines[i])
                # end
        # end
        fid.close()
        # pause
        return cometconfigfile
