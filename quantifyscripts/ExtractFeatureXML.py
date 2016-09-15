#!/usr/bin/python
# Python 2.7 Required
# Extract features from the input mzXML file
# Input:
# mzXML file
#
#
# Output:
# FeatureXML file
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------


from PipelineConfig import *
import sys
from WorkflowController import *
import glob


class ExportFeatureXML:
    def __init__(self, wfc, wp):
        print '[Info] Exporting featureXML file...'
        self.wfc = wfc
        self.wp = wp

    def Run(self, inputmzXML):
        filename = inputmzXML[0:-6]
        featureXML = filename + '.featureXML'

        if len(glob.glob(featureXML)) != 0 and not self.wfc.run_feature_finder:
            return featureXML
        # print filename
        print '[Info] Running Featurefinder, mzXML (MS1)-> featureXML '

        cmdline = '\"' + os.path.join(self.wp.OpenMSpath,
                                      'FileConverter' + self.wp.extname) + '\"' + ' -threads 18 -in_type mzXML -in ' + inputmzXML + ' -out ' + filename + '.tmp.mzML -out_type mzML'
        print '[Info] >' + cmdline
        os.system(cmdline)

        cmdline = '\"' + os.path.join(self.wp.OpenMSpath,
                                      'FileFilter' + self.wp.extname) + '\"' + ' -threads 18 -in_type mzML -in ' + filename + '.tmp.mzML' + ' -out ' + filename + '.mzML -out_type mzML'

        print '[Info] >' + cmdline
        os.system(cmdline)

        cmdline = '\"' + os.path.join(self.wp.OpenMSpath,
                                      'FeatureFinderCentroided' + self.wp.extname) + '\"' + ' -threads 18 -in ' + filename + '.mzML' + ' -out ' + filename + '.featureXML'
        print '[Info] >' + cmdline
        os.system(cmdline)
        featureXML = filename + '.featureXML'
        return featureXML


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print '[Info] extract the feature xml file from mzXML'
        print 'Usage:'
        print '%s  inputmzXML' % sys.argv[0]

    else:
        code_path = os.path.split(sys.argv[0])[0]
        wp = PipelineConfig()
        wp.load_from_config_file(os.path.join(code_path, 'RTQuant.ini'))
        wfc = WorkflowController()
        efxml = ExportFeatureXML(wfc,wp)
        efxml.Run(sys.argv[1])
