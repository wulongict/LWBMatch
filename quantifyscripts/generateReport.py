#!/usr/bin/python
# Python 2.7 Required
# Generating several summary files for png and jpg figures
# Input:
# path of all the figures
#
#
# Output:
# Several HTML files, each contains some figures about the alignment result
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------



import os
import glob
import sys

from abc import abstractmethod


class htmlReport:
    def __init__(self, inputpath_name):
        self.filename = inputpath_name
        self.fid = open(self.filename, 'w')

    def CreateHead(self):
        self.fid.write(
            '<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\"  \"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\"> '
            '<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\"><head><meta http-equiv=\"Content-Type\" '
            'content=\"text/html;charset=utf-8\" />')
        # self.fid.write('<link rel=\"stylesheet\" href=\"jemdoc.css\" type=\"text/css\" />')
        self.fid.write('<title>Report</title>')
        self.fid.write('  <style type="text/css">')
        self.fid.write("""
                        .figureFeatured {
                            display:table;
                            width:1px;
                            text-align:center;
                        }
                        .figureFeatured img {
                            display:block;
                            float:none;
                        }
                        .figcaption {
                            display:table;
                            clear:both
                        }
                        </style>"""
                       )
        self.fid.write('</head><body>')

    def CreateTail(self):
        self.fid.write('</body></html>')

    @abstractmethod
    def __exit__(self, exc_type, exc_val, exc_tb):
        print 'exit html report'

    def Close(self):
        print 'closing file'
        self.fid.close()
        return self.filename

    def InsertPic(self, pic):
        pictemplate = '<figure><img src="imgfile" alt="Report" width="600px" title="emptytitle"/> <figcaption>dummyname</figcaption></figure>'
        picStr = pictemplate.replace('imgfile', pic)
        picStr = picStr.replace("dummyname", os.path.split(pic)[1])
        picStr = picStr.replace("emptytitle", os.path.split(pic)[1])
        self.fid.write(picStr)

    def InsertPicList(self, piclist, figcaption=""):

        self.fid.write('<figure>')
        for eachpic in piclist:
            self.InsertPic(eachpic)
        self.fid.write('<figcaption>%s</figcaption></figure>' % figcaption)

    def InsertHyperLink(self, htmldictionary={}):
        for eachkey in htmldictionary.keys():
            self.fid.write("<p><a href=\"%s\">%s</a></p>" % (htmldictionary[eachkey], eachkey))


class SimpleFigureReport(htmlReport):
    def __init__(self, inputpath):
        htmlReport.__init__(self, os.path.join(inputpath, 'simplefigurereport.html'))
        self.inputpath = inputpath

    def __enter__(self):
        self.CreateHead()
        print self.inputpath
        picturelist = glob.glob(os.path.join(self.inputpath, '*.png'))
        print picturelist
        self.InsertPicList(picturelist, "")

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.CreateTail()
        self.Close()


def CreateDTWFigsAsHTML(filepath):
    hr = htmlReport(os.path.join(filepath, 'DTWFigs.html'))
    hr.CreateHead()
    filelist = glob.glob(os.path.join(filepath, '*same*.jpg'))
    print filelist
    hr.InsertPicList(filelist, "DTW aligned RT on similarities Matrix of MS1, MS2 and sum of all the MSn")

    filelist = glob.glob(os.path.join(filepath, '*same*.png'))
    print filelist
    hr.InsertPicList(filelist, "DTW aligned RT on similarities Matrix of MS1, MS2 and sum of all the MSn")

    filelist = glob.glob(os.path.join(filepath, '*imshow.jpg'))
    print filelist
    hr.InsertPicList(filelist, "DTW aligned RT on similarities Matrix of MS1, MS2 and sum of all the MSn")

    filelist = glob.glob(os.path.join(filepath, '*imshow.png'))
    print filelist
    hr.InsertPicList(filelist, "DTW aligned RT on similarities Matrix of MS1, MS2 and sum of all the MSn")

    hr.CreateTail()
    return hr.Close()


def CreateRTShiftsAsHTML(filepath):
    hr = htmlReport(os.path.join(filepath, 'RTShifts.html'))
    hr.CreateHead()
    filelist = glob.glob(os.path.join(filepath, '*resu*RTshift.png'))
    print filelist[0]
    hr.InsertPicList(filelist, "Real RT shift vs DTW Curve")
    hr.CreateTail()
    return hr.Close()


def CreateDotProductMatrixImage(filepath):
    hr = htmlReport(os.path.join(filepath, 'DotProduct.html'))
    hr.CreateHead()
    filelist = glob.glob(os.path.join(filepath, '*dp.txt.png'))
    print filelist[0]
    hr.InsertPicList(filelist, "Dot Product Images")
    hr.CreateTail()
    return hr.Close()


def CreateAlignmentReport(alignmentresult):
    hr = htmlReport(alignmentresult + 'alignment.html')
    hr.CreateHead()
    filelist = glob.glob(alignmentresult + "*.png")
    print filelist
    hr.InsertPicList(filelist, "LWBMatch alignment result")
    hr.CreateTail()
    return hr.Close()


def CreateAlignmentReportForPath(inputpath):
    hr = htmlReport(os.path.join(inputpath, "alignment_all.html"))
    hr.CreateHead()
    filelist = glob.glob(os.path.join(inputpath, "*.resu.gt.resu"))

    allfilelist = []
    for eachfile in filelist:
        allfilelist.extend(glob.glob(eachfile[:-8] + "*.png"))

    print allfilelist

    hr.InsertPicList(allfilelist, "Alignment of All")
    hr.CreateTail()
    return hr.Close()


def CreateDenoiseReport(inputpath):
    hr = htmlReport(os.path.join(inputpath, 'denoise.html'))
    hr.CreateHead()

    filelist = sorted(glob.glob(os.path.join(inputpath, '*mgfdenoise.png')))

    print filelist
    hr.InsertPicList(filelist, "performance of noise reduction")

    hr.CreateTail()
    return hr.Close()


def CreateXICReportZoomIn(inputpath):
    hr = htmlReport(os.path.join(inputpath, 'XICSummaryZoomIn.html'))
    hr.CreateHead()

    filelist = sorted(glob.glob(os.path.join(inputpath, '*ZoomIn*.png')))

    print filelist
    hr.InsertPicList(filelist, "Performance of RT warping for XIC")

    hr.CreateTail()
    return hr.Close()


def CreateXICReportZoomOut(inputpath):
    hr = htmlReport(os.path.join(inputpath, 'XICSummaryZoomOut.html'))
    hr.CreateHead()

    filelist = sorted(glob.glob(os.path.join(inputpath, '*ZoomOut.png')))

    print filelist
    hr.InsertPicList(filelist, "Performance of RT warping for XIC")

    hr.CreateTail()
    return hr.Close()


def create_scatterplot_spikedin_peptides_highlited(inputpath):
    hr = htmlReport(os.path.join(inputpath, 'SpikedInPeptidesHighlighted.html'))
    hr.CreateHead()

    filelist = sorted(glob.glob(os.path.join(inputpath, '*pep.intensity_scatter.png')))

    print filelist
    hr.InsertPicList(filelist, "SpikedInPeptide high lighted")

    hr.CreateTail()
    return hr.Close()


def CreateSummary(inputpath, HTMLfileDictionary):
    hr = htmlReport(os.path.join(inputpath, 'Datasetsummary.html'))
    hr.CreateHead()
    hr.InsertHyperLink(HTMLfileDictionary)
    hr.CreateTail()
    return hr.Close()


def CreateDTW_LOWESS_Comparison(inputpath):
    hr = htmlReport(os.path.join(inputpath, 'dtw_lowess.html'))
    hr.CreateHead()

    filelist = sorted(glob.glob(os.path.join(inputpath, '*dtw_lowess.png')))

    print filelist
    hr.InsertPicList(filelist, "DTW LOWESS Comparision")

    hr.CreateTail()
    return hr.Close()


def Create_DTW_ERROR_Boxplot_Report(inputpath):
    hr = htmlReport(os.path.join(inputpath, 'dtw_error_boxplot.html'))
    hr.CreateHead()

    filelist = sorted(glob.glob(os.path.join(inputpath, '*list_boxplot_error*.png')))

    print filelist
    hr.InsertPicList(filelist, "DTW Error Boxplot")

    hr.CreateTail()
    return hr.Close()


def GetReport(filepath='./'):
    if filepath == "":
        filepath = './'
    rtshiftshtml = CreateRTShiftsAsHTML(filepath)
    dtwhtml = CreateDTWFigsAsHTML(filepath)
    dphtml = CreateDotProductMatrixImage(filepath)
    denoisehtml = CreateDenoiseReport(filepath)
    xic_summary_zoom_in = CreateXICReportZoomIn(filepath)
    xic_summary_zoom_out = CreateXICReportZoomOut(filepath)
    spikedinpep = create_scatterplot_spikedin_peptides_highlited(filepath)
    dtw_lowess = CreateDTW_LOWESS_Comparison(filepath)
    dtw_error = Create_DTW_ERROR_Boxplot_Report(filepath)
    alignmenthtml = CreateAlignmentReportForPath(filepath)

    HTMLfileDictionary = {'DTW Error Boxplot': dtw_error, 'DTW vs LOWESS': dtw_lowess, \
                          'Spiked in peptide intensity scatter plot (OpenSWATH data only)': spikedinpep, \
                          'XIC mass trace summary (zoom out)': xic_summary_zoom_out, \
                          'XIC mass trace summary zoom in': xic_summary_zoom_in, 'Retention time shifts': rtshiftshtml, \
                          'DTW warping function': dtwhtml, 'Dot product matrix': dphtml, 'Noise reduction': denoisehtml, \
                          'Alignment': alignmenthtml}
    CreateSummary(filepath, HTMLfileDictionary)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'Usage: *.py figure_path'

    else:
        GetReport(sys.argv[1])
        # os.system('google-chrome ' + htmlfile)
