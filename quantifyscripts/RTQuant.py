#!/usr/bin/python
# Python 2.7 Required
# Python interface of swath data retention time alignment and quantification tool
# Input:
# SWATH data
# Technical or biological replicates
#
# Output:
# Scan base alignment
# Feature based alignment
# Intensity based quantification
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------


import argparse
import SwathmzXML
from DIAUmpireFlow import *
from ExtractFeatureXML import *
from LWBMatch import *
from ReadRecallPrecision import *
from WorkflowController import *


def RunDTW(filenamebase, wfc, SwathCyclePair, mzXMLpair, wp):
    print '[Info] Run DTW of ms1 alone'
    cmdline = os.path.join(wp.lwbmatchpath, 'mzXMLReader' + wp.extname)
    cmdline = cmdline + " "
    cmdline = cmdline + mzXMLpair[0] + " "
    cmdline = cmdline + mzXMLpair[1] + ' x'
    print '[Info] >' + cmdline

    filename = filenamebase + ".dat1.dp.txt.png"
    if not os.path.isfile(filename) or wfc.run_dtw:
        print 'filename:', filename
        print wfc.run_dtw
        os.system(cmdline)

    print '[Info] Run DTW of ms2 alone'
    cmdline = os.path.join(wp.lwbmatchpath, 'mzXMLReader' + wp.extname)
    cmdline = cmdline + " "
    cmdline = cmdline + mzXMLpair[0] + " "
    cmdline = cmdline + mzXMLpair[1] + ' --ms2'
    print '[Info] >' + cmdline

    filename = filenamebase + ".dat_sum_ms2.dp.txt.png"
    if not os.path.isfile(filename) or wfc.run_dtw:
        print 'filename:', filename
        print wfc.run_dtw
        os.system(cmdline)

    filename = filenamebase + ".dat1.dp.txt"
    import plotMatrix
    if not os.path.isfile(filename + '.png') or wfc.run_dtw:
        print '[Info] plot dot product matrix of ms1 and ms2 pairs'
        plotMatrix.plot_msn_dp_with_multi_thread(filename, SwathCyclePair[0])

    filename = filenamebase + '.dat_sum.dp.txt'
    if not os.path.isfile(filename + '.png') or wfc.run_dtw:
        print '[Info] plot dot product matrix of sumation of ms1 and ms2 together'
        plotMatrix.ReadPlotSumFile(filename)

    filename = filenamebase + '.dat_sum_ms2.dp.txt'
    if not os.path.isfile(filename + '.png') or wfc.run_dtw:
        print '[Info] plot dot product matrix of summation of ms2'
        plotMatrix.ReadPlotSumFile(filename)

    filename = filenamebase + '.dat.rt.matches'

    if not os.path.isfile(filename + '.png') or wfc.run_dtw:
        print filename
        print wfc.run_dtw
        print '[Info] plot rt matches curves individually'
        plotMatrix.plot_each_rt_matches_curve_with_multi_thread(filename, SwathCyclePair[0])

    print '[Info] plot rt matches curves into one figure'
    plotMatrix.plot_all_rt_matches_curve_into_single_file(filename, SwathCyclePair[0])


def create_featurexml_file_pair(mzxml_pair):
    """construct the featurexml_pair and writes it to the list file
    :param mzxml_pair: a list of two mzXML file
    @type mzxml_pair: list
    :rtype: str
    """
    featurexml_pair = [mzxml_pair[i][:-6] + '.featureXML' for i in range(len(mzxml_pair))]
    # featurexml_pair.append(mzxml_pair[0][:-6] + '.featureXML')
    # featurexml_pair.append(mzxml_pair[1][:-6] + '.featureXML')

    featurexml_list_file = featurexml_pair[1][0:-11] + '_' + os.path.split(featurexml_pair[0])[1] + 'List'

    fid = open(featurexml_list_file, 'w')
    fid.write(featurexml_pair[0] + '\n')
    fid.write(featurexml_pair[1])
    fid.close()
    print featurexml_list_file
    return featurexml_list_file


def calculateSWATHCycle(mzXMLpair):
    SwathCyclePair = []
    maxswath = SwathmzXML.CalculateSwathCycle(mzXMLpair[0])
    SwathCyclePair.append(maxswath)
    maxswath = SwathmzXML.CalculateSwathCycle(mzXMLpair[1])
    SwathCyclePair.append(maxswath)
    print '[Info] SWATH Cycle ', SwathCyclePair[0]
    return SwathCyclePair


def create_mzxmlpair(inputfile1, inputfile2):
    """return a pair of mzXML file
    :param inputfile1: one file
    :param inputfile2: another file
    :type inputfile1: str
    :type inputfile2: str
    :rtype: list
    """
    if os.path.split(inputfile1)[0] == "":
        inputfile1 = os.path.join("./", inputfile1)
    if os.path.split(inputfile2)[0] == "":
        inputfile2 = os.path.join("./", inputfile2)

    return [inputfile1, inputfile2]


def run_dia_umpire_comet_xinteract(inputmzXML, wfc, pipelinecfg):
    """
    :param wfc: work flow controller
    @type wfc: WorkflowController
    """
    workflowname = 'DIAUmpire -->-- Comet -->-- Xinteract Workflow'
    print '[Info] ' + workflowname + ' is running...'

    diaumpire = DIAUmpire(wfc, pipelinecfg)
    MGFList = diaumpire.Run(pipelinecfg, inputmzXML)
    # self.RunDIAUmpire()
    cometse = comet(wfc, pipelinecfg)
    fastapath = os.path.split(inputmzXML)[0]
    db = DB(wfc)
    pepxmllist = cometse.Run(MGFList, db.Run(fastapath, pipelinecfg))

    xint = Xinteract(wfc, pipelinecfg)
    xint.Run(pepxmllist)
    # xint.RunAll(pepxmllist)
    # RunXinteract()
    print '[Info] DIAUmpire -->-- Comet -->-- Xinteract Workflow finished'


def lwbmatch_workflow_groundtruth(mzXMLpair, wfc, wp):
    """run lwbmatch on two mzXML and calculate recall and precision on the ground truth"""
    # exchange mslevel MS1 <-> MS2
    print '[Info] DIAUmpire -->-- comet --> Xinteract -->-- FeatureFinder } step_1'
    print '[Info] DIAUmpire -->-- comet --> Xinteract -->-- FeatureFinder } step_2'
    print '[Info] {step_1 step_2 } -->-- Xinteract --> DTW-->-- LWBMatch -->-- groundtruth -->-- RecallPrecision -->-- Plot Figures'
    swath_cycle_pair = calculateSWATHCycle(mzXMLpair)

    run_dia_umpire_comet_xinteract(mzXMLpair[0], wfc, wp)

    efxml = ExportFeatureXML(wfc, wp)

    # efxml = ExportFeatureXML(wfc)
    efxml.Run(mzXMLpair[0])

    run_dia_umpire_comet_xinteract(mzXMLpair[1], wfc, wp)
    efxml.Run(mzXMLpair[1])

    featurexml_list_file = create_featurexml_file_pair(mzXMLpair)

    xint = Xinteract(wfc, wp)
    print '[Info] run xinteract to get the search results of the two file filtered together...'

    ipropepxml = xint.RunXinteractForPepXMLPair(mzXMLpair)

    filenamebase = featurexml_list_file[:-15]
    # filename = filenamebase + ".dat1.dp.txt.png"
    # if not os.path.isfile(filename):
    RunDTW(filenamebase, wfc, swath_cycle_pair, mzXMLpair, wp)

    lwbm = LWBMatch(wfc, wp)
    rt_alignment_result = lwbm.Run(featurexml_list_file)

    # self.RunXinteractForPepXMLPair()

    lwbm.calculate_recall_precision(ipropepxml, rt_alignment_result, mzXMLpair)
    # ExtractGroundtruth(ipropepxml, self.mzXMLpair,groundtruth)
    #
    # self.RunExtractGroundtruth()
    # self.calculate_recall_precision()

    lwbm.plot_figures_by_resu(rt_alignment_result)

    print 'Workflow 3 is done'
    # # self.align_each_swath()
    #
    # # continue run the workflow
    # # Here is where denoise is working
    # print '*' * 20
    # print '**************       Becareful ... We just skip the denoise part automatically *************'
    # print '*' * 20
    # return
    #
    # import DeNoise
    # mzML = DeNoise.RunSwathNoise(self.mzXMLpair[0], self.mzXMLpair[1], runanyway=True)
    # DeNoise.RunNoiseFilter(mzML)
    # self.run_dia_umpire_comet_xinteract(mzML[:-5] + '_denoised.mzXML')
    #
    # # ToDo: rt matches  should only calculate once.
    # mzML = DeNoise.RunSwathNoise(self.mzXMLpair[1], self.mzXMLpair[0], runanyway=True)
    # DeNoise.RunNoiseFilter(mzML)
    # self.run_dia_umpire_comet_xinteract(mzML[:-5] + '_denoised.mzXML')


def GenerateReport(mzXMLPath):
    import generateReport
    generateReport.GetReport(mzXMLPath)


def run(args):
    # wf = Workflow1.Workflow(args.inputfile1)
    wfc = WorkflowController()
    print wfc.__dict__
    wp = PipelineConfig()
    # wp.output_as_configure_file(args.Configfile)
    if os.path.isfile(args.Configfile):
        wp.load_from_config_file(args.Configfile)
    else:
        print '[Error] Configuration file %s does not exist.' % args.Configfile
        exit(0)
    # exit()
    mzXMLpair = create_mzxmlpair(args.inputfile1, args.inputfile2)
    lwbmatch_workflow_groundtruth(mzXMLpair, wfc, wp)
    mzXMLPath = os.path.split(mzXMLpair[0])[0]
    GenerateReport(mzXMLPath)


def parse_options():
    code_path = os.path.split(sys.argv[0])[0]
    parser = argparse.ArgumentParser(
        description='Align the retention time between two replicates of the same sample; Quantification of features based on intensity ')
    parser.add_argument("inputfile1", type=str, help="first replicate of SWATH data")
    parser.add_argument("inputfile2", type=str, help="second replicate of SWATH data")
    # parser.add_argument("-A", "--AlignmentMethod", type=str, default="MS1-DTW",
    #                     choices={"MS1-DTW", "MS2-DTW", "MS1+MS2-DTW", "LOWESS"},
    #                     help="methods of retention time alignment")
    parser.add_argument("-C", "--Configfile", default=os.path.join(code_path, "RTQuant.ini"), type=str,
                        help="configurations of software tools involved")
    # parser.add_argument("-O", "--OutDir", default="./", type=str, help="path for output")
    # parser.add_argument("-N", "--OutFileName", default="test", type=str, help="name / identifier of output files")
    args = parser.parse_args()
    # print args.help
    # print args.AlignmentMethod

    run(args)


if __name__ == "__main__":
    parse_options()
