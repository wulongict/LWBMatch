#!/usr/bin/python
# Python 2.7 Required
# Recall and precision reader, This program is used to check the index of mzXML files.
# dependency: indexmzXML from TPP
# input:
# mzXML file
#
# output:
# updated mzXML file
#
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------


from PipelineConfig import *
import shutil
import sys
import os
from WorkflowController import *


class Converter:
    """ return nothing. """

    def __init__(self, wfc,wp):
        """     no return.
        :type wfc: WorkflowController
        @type wfc: WorkflowController
        @type wp: PipelineConfig
        """
        print '[Info] Starting a converter...'
        self.wfc = wfc
        self.wp = wp

    def check_mzxml(self, inputmzxml):
        """update index
        @type inputmzxml: string
        :param inputmzxml: string
        """
        print '[Info] Starting check'
        # wp = PipelineConfig()
        cmdline = os.path.join(self.wp.TPPpath, 'indexmzXML' + self.wp.extname) + ' ' + inputmzxml

        print '[Info] >' + cmdline
        if self.wfc.run_indexmzXML:  # this is a very important bug, fix it on Nov 30. 2015. check mzXML anyway.
            os.system(cmdline)
        possible_output = inputmzxml + '.new'

        if os.path.isfile(possible_output):
            print '[Info] Updating mzXML file...'
            shutil.move(possible_output, inputmzxml)
        print '[Info] Checking completes ...'


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print "Usage:\n%s inputmzXML" % (sys.argv[0])
    else:
        print sys.argv
        code_path = os.path.split(sys.argv[0])[0]
        wp = PipelineConfig()
        wp.load_from_config_file(os.path.join(code_path, 'RTQuant.ini'))
        wfc = WorkflowController()
        cnvt = Converter(wfc, wp)
        cnvt.check_mzxml(sys.argv[1])
