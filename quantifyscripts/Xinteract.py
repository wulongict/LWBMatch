#!/usr/bin/python
# Python 2.7 Required
# Xinteract module. This script can not run by itself. Use another script to call it.
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------




from PipelineConfig import *

import glob





class Xinteract:
    def __init__(self, wfc, wp):
        """
        :type wfc: WorkflowController
        :param wfc: work flow controller
        :type wp: PipelineConfig
        """
        print '[Info] Runnint Xinteract...'
        self.wfc = wfc
        self.wp = wp
        self.params = ' -dDECOY -p0 -ip -OPd '
        self.params = ' -dDECOY -p0 -ip -OPdAt '

    def GetOutputPath(self, pepxml):
        """
        :type pepxml: str
        """
        outputpath = os.path.split(pepxml)[0]

        if self.wp.extname != '':  # Running linux
            outputpath = os.path.join(self.wp.TPPpath, r'../wwwroot')
        # end if
        if outputpath == '':
            outputpath = './'
        return outputpath

    def GetPepXMLFrommzXML(self, inputmzXMLpair):
        pepxmlpair = []
        for eachmzXML in inputmzXMLpair:
            pepxmlname = eachmzXML[:-6] + '_Q1.' + 'pep.xml'
            pepxmlpair.append(pepxmlname)
        return pepxmlpair

    def RunXinteractForPepXMLPair(self, inputmzXMLpair):
        if len(inputmzXMLpair) == 0:
            print '[Info] No mzXML found..'
            exit(0)
        xinteractoutputpath = self.GetOutputPath(inputmzXMLpair[0])
        print '[Info] Running Xinteract for pepxml pair'

        os.chdir(xinteractoutputpath)
        XinteractCMD = "\"" + os.path.join(self.wp.TPPpath, 'xinteract' + self.wp.extname) + "\"" + self.params + ' -N'

        pepxml_pair = self.GetPepXMLFrommzXML(inputmzXMLpair)

        pepxmlname = os.path.split(pepxml_pair[0])[1][:-8] + '_' + os.path.split(pepxml_pair[1])[1][:-8]

        output_filename = "interact-" + pepxmlname + '.ipro.pep.xml'

        ipro_for_pair_mzXML = os.path.join(xinteractoutputpath, output_filename)

        if os.path.isfile(ipro_for_pair_mzXML) and not self.wfc.run_xinteract:
            # never skip this
            return ipro_for_pair_mzXML

        cmdline = XinteractCMD + pepxmlname + " " + pepxml_pair[0] + " " + pepxml_pair[1]
        print '[Info] >' + cmdline
        os.system(cmdline)
        return ipro_for_pair_mzXML
        # last function!

    def RunAll(self, pepxmllist):
        ipro_pepxml_list = []
        if len(pepxmllist) == 0:
            print '[Info] No input *pep.xml file!'
            exit(0)

        xinteractoutputpath = self.GetOutputPath(pepxmllist[0])
        os.chdir(xinteractoutputpath)
        XinteractCMD = "\"" + os.path.join(self.wp.TPPpath, 'xinteract' + self.wp.extname) + "\"" + self.params + ' -N'
        pepxmlname = ''
        for each_pepxml in pepxmllist:
            pepxmlname = pepxmlname + os.path.split(each_pepxml)[1][0:-8]

        cmdline = XinteractCMD + pepxmlname
        outputfilename = "interact-" + pepxmlname + '.ipro.pep.xml'
        ipro_pepxml_list.append(outputfilename)
        if len(glob.glob(outputfilename)) != 0 and not self.wfc.run_xinteract:
            return ipro_pepxml_list

        for each_pepxml in pepxmllist:
            cmdline = cmdline + ' ' + each_pepxml

        print '[Info] >' + cmdline
        os.system(cmdline)

        return ipro_pepxml_list

    def Run(self, pepxmllist):
        ipropepxmllist = []
        if len(pepxmllist) == 0:
            print '[Info] No input *pep.xml file!'
            exit(0)

        xinteract_output_path = self.GetOutputPath(pepxmllist[0])

        print xinteract_output_path
        os.chdir(xinteract_output_path)
        XinteractCMD = "\"" + os.path.join(self.wp.TPPpath,
                                           'xinteract' + self.wp.extname) + "\"" + ' -dDECOY -p0 -ip -OPd -N'
        for eachpepxml in pepxmllist:
            pepxmlname = os.path.split(eachpepxml)[1][0:-8]
            output_filename = "interact-" + pepxmlname + '.ipro.pep.xml'
            ipropepxmllist.append(output_filename)
            if os.path.isfile(output_filename) and not self.wfc.run_xinteract:
                continue

            if not self.wfc.run_xinteractQ1 and 'Q1.pep.xml' in eachpepxml:
                continue

            if not self.wfc.run_xinteractQ2 and 'Q2.pep.xml' in eachpepxml:
                continue

            if not self.wfc.run_xinteractQ3 and 'Q3.pep.xml' in eachpepxml:
                continue

            cmdline = XinteractCMD + pepxmlname + " " + eachpepxml
            print '[Info] >' + cmdline
            os.system(cmdline)

        # end
        return ipropepxmllist
        # print 'xinteract ', XinteractCMD
        # os.system(XinteractCMD)


class High_Resolution_Xinteract(Xinteract):
    def __init__(self):
        pass
