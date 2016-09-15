#!/usr/bin/python
# Python 2.7 Required
# Run LWBMatch on input FeatureXMLList file
# Input:
# FeatureXMLList file
#
# Output:
# Resu file, feature based alignment result
# dp.txt file, dot product results
# rt.matches, retention time warping curve by DTW (Dynamic time warping)
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------


from PipelineConfig import *
import plotQuiver
from WorkflowController import *
import glob
import sys


class LWBMatch:
    def __init__(self, wfc, wp):
        """
        :param wfc: work flow controller
        @type wfc: WorkflowControler
        :param wp: workpath configurations
        @type wp: PipelineConfig
        """
        print '[Info] LWBMatch is running...'
        self.wfc = wfc
        self.wp = wp

    def Run(self, fileFeatureXMLList):
        RTAlignResult = fileFeatureXMLList[0:-14] + 'resu'
        cmdline = os.path.join(self.wp.lwbmatchpath,
                               'lwbmatch' + self.wp.extname) + '  -l ' + fileFeatureXMLList + ' -o ' + RTAlignResult + ' -w 2'
        print '[Info] >' + cmdline
        if not os.path.isfile(RTAlignResult) or self.wfc.run_lwbmatch:
            # always run True
            os.system(cmdline)
        # self.plot_figures_by_resu(RTAlignResult)
        return RTAlignResult

    def plot_figures_by_resu(self, RTAlignResult):
        print '[Info] Plot result of LWBmatch plot results of lwbmatch  =============='
        # Load the results
        print RTAlignResult
        if len(glob.glob(RTAlignResult + '_QuiverDotPlot.png')) != 0 and not self.wfc.run_plot_figure:
            print '[Info] figures exists! exit!'
            return

        plotQuiver.DrawFiguresOfResuFile(RTAlignResult)
        # reader = plotQuiver.LWBMatchResuReader(RTAlignResult)
        # reader.Read()
        # reader.FilterZeros()
        # reader.DrawIntensityAlignCurve()
        # reader.DrawMZAlignCurve()
        # reader.DrawRTAlignCurve()
        # reader.DrawRTShiftDistribution()
        # reader.DrawQuiverPlot()

    def calculate_recall_precision(self, ipropepxml, RTAlignResult, mzXMLpair):
        print '[Info] calculate recall & precision ==========='
        groundtruth = RTAlignResult + '.gt.resu'
        if self.wfc.run_extract_ground_truth or not os.path.isfile(groundtruth):
            ExtractGroundtruth(ipropepxml, mzXMLpair, groundtruth, self.wp)
        if self.wfc.run_recall_precision or True:
            CalculateRecallPrecision(groundtruth, RTAlignResult, self.wp)


def OutputmzXMLpair(mzXMLpair, mzXMLListtxtName):
    print '[Info] Write to mzXML.txt file'
    if not os.path.isfile(mzXMLListtxtName):
        fid = open(mzXMLListtxtName, 'w')
        fid.write(mzXMLpair[0] + "\n")
        fid.write(mzXMLpair[1] + "\n")
        fid.close()
    return


def ExtractGroundtruth(ipropepxml, mzXMLpair, groundtruth, wp):
    print 'Extracting Ground truth'
    print mzXMLpair
    # write mzXML to txt
    mzXMLListtxtName = groundtruth + '_mzXMLList.txt'
    OutputmzXMLpair(mzXMLpair, mzXMLListtxtName)

    cmdline = os.path.join(wp.binpath,
                           'groundtruth' + wp.extname) + " 2 " + mzXMLListtxtName + " " + ipropepxml + " 0.9 " + groundtruth + " 1"
    print '[Info] >' + cmdline
    # if not os.path.isfile(groundtruth):
    os.system(cmdline)
    # exit(0)

    return


def CalculateRecallPrecision(groundtruth, RTAlignResult, wp, recordsNum=50):
    # N = recordsNum
    # max_rt = 25.0
    # rt_step = max_rt / N
    # mz = [0.1 + i*0 for i in range(N)]
    # rt = [1 + i * rt_step for i in range(N)]
    #
    # for i in range(N):
    #     print 'recall-precision-mz-rt'
    #     cmdline = os.path.join(wp.recallprecision_binpath, 'recallprecision' + wp.extname) + " " + groundtruth + " " + RTAlignResult + " " + str(mz[i]) + " " + str(rt[i])
    #     print '[Info] >' + cmdline
    #     os.system(cmdline)

    print 'Calculating recall precisoin with 0.1 Th mz_threshold'
    cmdline = os.path.join(wp.recallprecision_binpath,
                           'recallprecision' + wp.extname) + " " + groundtruth + " " + RTAlignResult + " 0.05 0.15"
    print '[Info] >' + cmdline
    os.system(cmdline)

    print 'Calculating recall precisoin with 0.5 Th mz_threshold'
    cmdline = os.path.join(wp.recallprecision_binpath,
                           'recallprecision' + wp.extname) + " " + groundtruth + " " + RTAlignResult + " 50 150"
    print '[Info] >' + cmdline
    os.system(cmdline)
    print 'recall precision done'


def RunLWBMatchWorkFlow(fileFeatureXMlList, wp):
    print '[Info] Running LWBMatch'
    wfc.run_plot_figure = True
    wfc.run_extract_ground_truth = True
    wfc.run_recall_precision = True
    lwbmatch = LWBMatch(wfc, wp)
    RTAlignResult = lwbmatch.Run(fileFeatureXMlList)
    lwbmatch.plot_figures_by_resu(RTAlignResult)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print 'usage: %s featurexmllist ' % (sys.argv[0])
    else:
        code_path = os.path.split(sys.argv[0])[0]
        wp = PipelineConfig()
        wp.load_from_config_file(os.path.join(code_path, 'RTQuant.ini'))
        wfc = WorkflowController()
        RunLWBMatchWorkFlow(sys.argv[1], wp)
