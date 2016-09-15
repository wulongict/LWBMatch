# -*- coding: utf-8 -*-
# !/usr/bin/python
# Python 2.7 Required
# plot figures based on the input alignment result
# Input:
# resu file by LWBMatch
# Output:
# Several different figures related to the alignment result (resu file)
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------

from abc import abstractmethod
from abc import ABCMeta
import sys
import numpy as np
import math



import matplotlib.pyplot as plt




class HistgramPloter:
    def __init__(self, datalist, outputfilename, xlabel, ylabel='Probability', facecolor='green', number_of_bins=100):
        print '[Info] Plotting histgram ' + outputfilename
        self.outputfilename = outputfilename
        self.ylabel = ylabel
        self.xlabel = xlabel
        self.facecolor = facecolor
        self.datalist = datalist
        self.number_of_bins = number_of_bins

    def plot(self):
        [maxdata, mindata] = [max(self.datalist), min(self.datalist)]
        step = (maxdata - mindata) / (self.number_of_bins - 1)

        # This is a left_close_right_open interval
        left_boundaries = [mindata + i * step for i in range(0, self.number_of_bins)]
        n, bins, patches = plt.hist(self.datalist, left_boundaries, normed=1, facecolor=self.facecolor,
                                    edgecolor=self.facecolor, alpha=0.5)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.xlim([mindata - step, maxdata + 2 * step])
        plt.savefig(self.outputfilename, dpi=180)
        # plt.savefig(self.inputfile+'_MZDeviationDistribution.png',dpi=180)
        plt.close()


class ResuReader:
    __metaclass__ = ABCMeta

    def __init__(self, filename):
        self.inputfile = filename
        self.data = []
        self.filenames = []
        self.type_name = ""

    def Print(self):
        for eachrow in self.data:
            print eachrow

    @abstractmethod
    def GetSelfType(self):
        pass

    @abstractmethod
    def GetRTList(self):
        pass

    def DrawRTAlignCurve(self):
        RT1, RT2 = self.GetRTList()
        if len(RT1) == 0 or len(RT2) == 0:
            print "Empty RT, check the data please"
            return

        plt.plot(RT1, RT2, 'g<', markersize=0.8)
        plt.xlabel('RT1 (sec)')
        plt.ylabel('RT2 (sec)')
        # plt.show()
        plt.savefig(self.inputfile + '_RTAlgin.png', dpi=360)
        plt.close()

    def DrawRTshift_MS1(self):
        name_of_font = 'Arial'
        size_of_font = 12
        RT1, RT2 = self.GetRTList()
        if len(RT1) == 0 or len(RT2) == 0:
            print "Empty RT, check the data please"
            return

        rtshift = []
        for i in range(len(RT1)):
            rtshift.append(float(RT1[i]) - float(RT2[i]))
        # end


        plt.plot(RT1, rtshift, 'g.', markersize=3.0, label=self.type_name)
        plt.xlabel('$RT_1$ (sec)')
        plt.ylabel('$\Delta$ $RT$ = $RT_1$ - $RT_2$ (sec)')

        x, z = self.get_warping_delta()
        std_times = 3.0
        zmin, zmax = np.mean(z) - np.std(z) * std_times, np.mean(z) + np.std(z) * std_times
        plt.plot(x, z, 'b-', markersize=1, label='DTW on MS1')
        plt.legend(loc='best')
        plt.ylim([zmin, zmax])
        plt.savefig(self.inputfile + "_RTshiftonlyMS1.png", dpi=90)
        plt.close()

        # Figure 2 (a,b) in the paper # boxplot overlaid picture is somewhere else (figure b)
        # using grayscale here
        plt.plot(RT1, rtshift, markersize=3.0, label=self.type_name, alpha=0.7, linestyle='None', marker='.',
                 markerfacecolor='#808080', markeredgecolor='#808080')
        plt.xlabel('$RT_1$ (sec)', fontname=name_of_font, fontsize=size_of_font)
        plt.ylabel('$\Delta$ $RT$ = $RT_1$ - $RT_2$ (sec)', fontname=name_of_font, fontsize=size_of_font)

        x, z = self.get_warping_delta()
        std_times = 3.0
        zmin, zmax = np.mean(z) - np.std(z) * std_times, np.mean(z) + np.std(z) * std_times
        plt.plot(x, z, 'k-', markersize=1, label='DTW on MS1', linewidth=2.0)
        plt.legend(loc='best')
        plt.ylim([zmin, zmax])
        plt.savefig(self.inputfile + "_grayscale_RTshiftonlyMS1.png", dpi=90)
        plt.close()

    def get_warping_delta(self):  # only ms1
        commonfilename = self.inputfile.replace(".resu.gt.resu", "")
        commonfilename = commonfilename.replace(".resu", "")
        print commonfilename
        RTMatchesFile = commonfilename + ".dat1.dp.rt.matches"
        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]
        return x, z

    def get_rt_shift_estimated_by_DTW_on_MS1_MS2(self):
        commonfilename = self.inputfile.replace(".resu.gt.resu", "")
        commonfilename = commonfilename.replace(".resu", "")
        # print commonfilename

        RTMatchesFile = commonfilename + ".dat.rt.matches"

        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]

        return x, z

    def get_rt_shift_estimated_by_DTW_on_MS1(self):
        commonfilename = self.inputfile.replace(".resu.gt.resu", "")
        commonfilename = commonfilename.replace(".resu", "")
        # print commonfilename

        RTMatchesFile = commonfilename + ".dat1.dp.rt.matches"

        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]

        return x, z

    def get_rt_shift_estimation_SD_on_MS1(self):
        x, z = self.get_rt_shift_estimated_by_DTW_on_MS1()
        RT1, RT2 = self.GetRTList()
        if len(RT1) == 0 or len(RT2) == 0:
            print "Empty RT, check the data please"
            return
        rtshift = []
        for i in range(len(RT1)):
            rtshift.append(float(RT1[i]) - float(RT2[i]))
        # end

        estimation_errors = []
        print '[Info] Calculating error ...'
        for i in range(len(RT1)):
            eachrt = float(RT1[i])
            Real_rt_shift = rtshift[i]

            # rt_distance = 100000
            # estimated_rt_shift = 0
            # matched_rt = 0
            # for j in range(len(x)):
            #     currentrt = float(x[j])
            #
            #     if abs(currentrt - eachrt) < rt_distance:
            #         matched_rt = currentrt
            #         rt_distance = abs(currentrt - eachrt)
            #         estimated_rt_shift = float(z[j])
            #     #end
            # print 'pair or rt: ', matched_rt, eachrt
            estimated_rt_shift = z[np.abs(np.array(x) - eachrt).argmin()]
            # estimated_rt_shift = z[min(range(len(x)), key=lambda i: abs(x[i]-eachrt))]
            # print estimated_rt_shift
            # pause
            # abs_rt_shift = [rt_item - eachrt for rt_item in x]
            # a,b = min(abs_rt_shift)
            # print a, b
            # pause
            error_of_estimation = estimated_rt_shift - Real_rt_shift
            # print 'pair of error:', estimated_rt_shift, Real_rt_shift
            estimation_errors.append(error_of_estimation)
        print '[Info] Standard deviation of error of estimation (MS1): ', np.std(estimation_errors)
        sd_estimation = np.std(estimation_errors)
        return sd_estimation, estimation_errors, rtshift

    def get_rt_shift_estimation_SD_on_MS1_MS2(self):
        x, z = self.get_rt_shift_estimated_by_DTW_on_MS1_MS2()
        RT1, RT2 = self.GetRTList()
        if len(RT1) == 0 or len(RT2) == 0:
            print "Empty RT, check the data please"
            return
        rtshift = []
        for i in range(len(RT1)):
            rtshift.append(float(RT1[i]) - float(RT2[i]))
        # end

        estimation_errors = []
        print '[Info] Calculating error ...'
        for i in range(len(RT1)):
            eachrt = float(RT1[i])
            Real_rt_shift = rtshift[i]

            # rt_distance = 100000
            # estimated_rt_shift = 0
            # matched_rt = 0
            # for j in range(len(x)):
            #     currentrt = float(x[j])
            #
            #     if abs(currentrt - eachrt) < rt_distance:
            #         matched_rt = currentrt
            #         rt_distance = abs(currentrt - eachrt)
            #         estimated_rt_shift = float(z[j])
            #     #end
            # # print 'pair or rt: ', matched_rt, eachrt
            estimated_rt_shift = z[np.abs(np.array(x) - eachrt).argmin()]
            error_of_estimation = estimated_rt_shift - Real_rt_shift
            estimation_errors.append(error_of_estimation)
        print '[Info] Standard deviation of error of estimation (MS1 & MS2): ', np.std(estimation_errors)
        sd_estimation = np.std(estimation_errors)
        return sd_estimation, estimation_errors, rtshift

    def DrawRTShiftDistribution(self):
        RT1, RT2 = self.GetRTList()
        if len(RT1) == 0 or len(RT2) == 0:
            print "Empty RT, check the data please"
            return

        rtshift = []
        for i in range(len(RT1)):
            rtshift.append(float(RT1[i]) - float(RT2[i]))
        # end

        plt.plot(RT1, rtshift, 'g.', markersize=3.0, label=self.type_name)
        plt.xlabel('$RT_1$ (sec)')
        plt.ylabel('$\Delta$ $RT$ = $RT_1$ - $RT_2$ (sec)')

        # plt.ylim([-500,500]) # No y-limit any more
        commonfilename = self.inputfile.replace(".resu.gt.resu", "")
        commonfilename = commonfilename.replace(".resu", "")
        print commonfilename

        RTMatchesFile = commonfilename + ".dat.rt.matches"

        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]
        plt.plot(x, z, 'r-', markersize=2, label='DTW on MS1+MS2')
        std_times = 3.0
        zmin, zmax = np.mean(z) - np.std(z) * std_times, np.mean(z) + np.std(z) * std_times

        # -----------------------added for ms2 alone version -------------------
        RTMatchesFile = commonfilename + ".dat_sum_ms2.rt.matches"
        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]
        zmin2, zmax2 = np.mean(z) - np.std(z) * std_times, np.mean(z) + np.std(z) * std_times

        zmin, zmax = min([zmin, zmin2]), max([zmax, zmax2])
        plt.plot(x, z, 'k->', markersize=1, label='DTW on MS2')
        # ------------------------end-------------------------------------------


        RTMatchesFile = commonfilename + ".dat1.dp.rt.matches"
        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]
        zmin1, zmax1 = np.mean(z) - np.std(z) * std_times, np.mean(z) + np.std(z) * std_times
        plt.plot(x, z, 'b-', markersize=1, label='DTW on MS1')

        plt.legend(loc='best')

        plt.ylim(min([zmin, zmin1]), max([zmax, zmax1]))

        plt.savefig(self.inputfile + "_RTshift.png", dpi=90)
        plt.close()
        # example data
        # mu = 100 # mean of distribution
        # sigma = 15 # standard deviation of distribution
        # x = mu + sigma * np.random.randn(10000)

        print 'max and min of rtshift are: ', min(rtshift), max(rtshift)

        [rtshiftmin, rtshiftmax] = [min(rtshift), max(rtshift)]
        number_of_bins = 100
        step = (rtshiftmax - rtshiftmin) / (number_of_bins - 1)
        bin_centers = [rtshiftmin + i * step for i in range(0, 100)]
        # the histogram of the data
        n, bins, patches = plt.hist(rtshift, bin_centers, normed=1, facecolor='green', alpha=0.5)
        # add a 'best fit' line
        # y = mlab.normpdf(bins, mu, sigma)
        # plt.plot(bins, y, 'r--')
        plt.xlabel('RT shifts (sec)')
        plt.ylabel('Probability')
        plt.title(r'Distribution of RT shifts')

        plt.savefig(self.inputfile + '_RTshiftsDistribution.png', dpi=360)
        plt.close()

    def ExtractColumn(self, colnum):
        col = []
        for eachitem in self.data:
            col.append(eachitem[colnum])
        return col

    @abstractmethod
    def GetMzList(self):
        pass

    def DrawMZAlignCurve(self):
        MZ1, MZ2 = self.GetMzList()
        diffmz = []
        for i in range(len(MZ1)):
            diffmzitem = float(MZ1[i]) - float(MZ2[i])
            diffmzitem = 1000000 * diffmzitem
            diffmzitem = diffmzitem / float(MZ1[i])

            diffmz.append(diffmzitem)
        plt.plot(MZ1, MZ2, 'go', markersize=0.5)
        plt.xlabel('MZ1 (Th)')
        plt.ylabel('MZ2 (Th)')
        # plt.show()
        plt.savefig(self.inputfile + '_MZAlgin.png', dpi=360)
        plt.close()
        plt.plot(MZ1, diffmz, 'go', markersize=1.0)
        plt.xlabel('MZ1 (Th)')
        plt.ylabel('MZ deviation (ppm)')
        plt.savefig(self.inputfile + '_mzdiff.png', dpi=180)
        plt.close()

        hp = HistgramPloter(diffmz, self.inputfile + '_MZDeviationDistribution.png', 'MZ deviation (ppm)')
        hp.plot()
        # [maxmz, minmz] = [max(diffmz), min(diffmz)]
        # number_of_bins = 100
        # step = (maxmz - minmz) / (number_of_bins -1)
        # bin_centers = [minmz + i*step for i in range(0,number_of_bins)]
        # n, bins, patches = plt.hist(diffmz, bin_centers, normed=1, facecolor='green', alpha=0.5)
        # plt.xlabel('MZ deviation (ppm)')
        # plt.ylabel('Probability')
        # plt.savefig(self.inputfile+'_MZDeviationDistribution.png',dpi=180)
        # plt.close()

    def DrawQuiverPlot(self):
        RT1, RT2 = self.GetRTList()
        MZ1, MZ2 = self.GetMzList()
        # RT1 = self.ExtractColumn(0);
        # MZ1 = self.ExtractColumn(1)
        # #Int1 = self.ExtractColumn(2)
        # RT2 = self.ExtractColumn(3);
        # MZ2 = self.ExtractColumn(4)
        # #Int2 = self.ExtractColumn(5)

        plt.plot(RT1, MZ1, 'b.')

        plt.xlabel('$RT_1$ (sec)')
        plt.ylabel('$MZ_1$ (Th)')
        plt.savefig(self.inputfile + '_DotPlotFirst.png', dpi=360)
        plt.close()

        plt.plot(RT2, MZ2, 'g.')
        plt.xlabel('$RT_2$ (sec)')
        plt.ylabel('$MZ_2$ (Th)')
        plt.savefig(self.inputfile + '_DotPlotSecond.png', dpi=360)
        plt.close()

        plt.plot(RT1, MZ1, 'b.', RT2, MZ2, 'g.')
        plt.xlabel('RT (sec)')
        plt.ylabel('MZ (Th)')

        plt.savefig(self.inputfile + '_DotPlot.png', dpi=360)
        plt.close()

        plt.plot(RT1, MZ1, 'b.', RT2, MZ2, 'g.')
        plt.xlabel('RT (sec)')
        plt.ylabel('MZ (Th)')
        for i in range(len(RT1)):
            plt.plot([RT1[i], RT2[i]], [MZ1[i], MZ2[i]], 'r-')

        plt.savefig(self.inputfile + '_QuiverDotPlot.png', dpi=360)
        plt.close()


def loadRTMatches(RTMatchesFile):
    print '[Info] Loading ' + RTMatchesFile
    fid = open(RTMatchesFile, "r")
    lines = fid.readlines()
    fid.close()
    lines = [eachline.strip() for eachline in lines]
    lines = [eachline.split() for eachline in lines]
    lines = [[float(item) for item in eachline] for eachline in lines]
    # print lines[0]
    x = [eachline[0] for eachline in lines]
    y = [eachline[1] for eachline in lines]
    return x, y


class GroundTruthReader(ResuReader):
    def __init__(self, filename):
        ResuReader.__init__(self, filename)
        self.type_name = "ground truth map"
        self.peptideseq = []

    def Read(self):
        fid = open(self.inputfile, 'r')
        lines = fid.readlines()
        fid.close()

        lines = [eachline.strip() for eachline in lines]
        lines = [eachline.split() for eachline in lines]
        self.peptideseq = [eachline[0] for eachline in lines]
        self.peptideseq = self.peptideseq[2:]
        lines = [eachline[1:] for eachline in lines]
        self.filenames = lines[1]
        print '[Info] filenames line in ground truth ', self.filenames
        self.data = lines[2:]
        if self.data == []:
            print "No Ground truth found, check the search result please. There might be no identifications."

    def GetRTList(self):
        RTcol = 0;
        RT1 = self.ExtractColumn(RTcol)
        RT2 = self.ExtractColumn(RTcol + 2);
        return RT1, RT2

    def GetMzList(self):
        MZcol = 1;
        MZ1 = self.ExtractColumn(MZcol)
        MZ2 = self.ExtractColumn(MZcol + 2);
        return MZ1, MZ2

    def GetSelfType(self):
        return self.type_name


class LWBMatchResuReader(ResuReader):
    def __init__(self, filename):
        # why this is incorrect?
        ResuReader.__init__(self, filename)
        self.type_name = "predicted map"

    def Read(self):
        fid = open(self.inputfile, 'r')
        lines = fid.readlines()
        fid.close()
        lines = [eachline.strip() for eachline in lines]
        lines = [eachline.split() for eachline in lines]
        self.filenames = lines[1]
        self.data = lines[2:]

    def FilterZeros(self):
        tmp = []
        for eachitem in self.data:
            if '0.000' in eachitem:
                continue
            else:
                tmp.append(eachitem)
        self.data = tmp

    def GetFilenames(self):
        return self.filenames

    def GetRTList(self):
        RTcol = 0;
        RT1 = self.ExtractColumn(RTcol)
        RT2 = self.ExtractColumn(RTcol + 3);
        if len(RT1) != len(RT2):
            print "Error in RT list"
        return RT1, RT2

    def GetMzList(self):
        MZcol = 1;
        MZ1 = self.ExtractColumn(MZcol)
        MZ2 = self.ExtractColumn(MZcol + 3);
        return MZ1, MZ2

    def GetIntenList(self):
        Intencol = 2;
        Inten1 = self.ExtractColumn(Intencol)
        Inten2 = self.ExtractColumn(Intencol + 3);
        return Inten1, Inten2

    def DrawIntensityAlignCurve(self):
        # Draw the scatter plot for alignment of intensity
        # Figure 1, RAW intensity
        col = 2;
        Int1 = self.ExtractColumn(col)
        Int2 = self.ExtractColumn(col + 3);
        plt.plot(Int1, Int2, 'g.')
        plt.xlabel('Intensity1')
        plt.ylabel('Intensity2')

        # plt.show()
        plt.savefig(self.inputfile + '_IntAlgin.png', dpi=360)
        plt.close()

        # Figure 2, log scaled intensity

        logint1 = []
        logint2 = []
        for i in range(len(Int1)):
            logint1.append(math.log(float(Int1[i]), 2.0))
            logint2.append(math.log(float(Int2[i]), 2.0))
        plt.plot(logint1, logint2, 'g.', markersize=1)

        log2_shift = [logint2[i] - logint1[i] for i in range(len(logint2))];
        median_log_shift = np.median(log2_shift)
        mean_log_shift = np.mean(log2_shift)

        xitem = np.arange(min(logint1), max(logint2))
        yitem = [eachxitem + mean_log_shift for eachxitem in xitem]
        plt.plot(xitem, yitem, 'r-.')
        plt.legend(["predicted feature pairs", "estimated regression line"])
        plt.xlabel('$log_2$(Intensity1)')
        plt.ylabel('$log_2$(Intensity2)')

        # plt.show()
        plt.savefig(self.inputfile + '_logIntAlgin.png', dpi=360)
        plt.close()

        # Figure 3, hist 2d, a heat map of the results
        plt.hist2d(logint1, logint2, bins=100)
        plt.colorbar()
        plt.savefig(self.inputfile + "_logIntAlign_hist2d.png", dpi=360)
        plt.close()

        # Figure 4, histgram, the difference of log scaled intensity difference.
        plt.hist(log2_shift, bins=500)
        plt.ylabel("frequency")
        plt.xlabel("$Intensity Ratio (log_2 scale)$")
        plt.title(r"$Mean: \bar{x}$ = %.2lf  Median: $\~x$ = %.2lf" % (mean_log_shift, median_log_shift))
        plt.savefig(self.inputfile + "_logIntShiftsAlign_hist.png", dpi=360)
        plt.close()
        from scatter_hist import scatter_hist_plot
        # scatter_hist_plot(logint1,logint2)

    def split_rt_shit_based_on_rt(self, rtshift, rt, interval=100):
        # split the rtshift data into several parts, based on the rt windows.
        # each 1000s as one part
        # interval = 100
        rtshift_parts = []
        rt_position = []
        for i in range(int(max(rt) / interval) + 1):
            rtshift_parts.append(
                [rtshift[j] for j in range(len(rt)) if rt[j] >= i * interval and rt[j] < (i + 1) * interval])
            rt_position.append((i + 1 / 2.0) * interval)
        # print rtshift_parts
        return rtshift_parts, rt_position

    def plot_boxplots_on_RT(self):  # Figure 2b
        # plot a new figure for the delta rt
        # fixed the x axis tick overlap problem
        rt1, rt2 = self.GetRTList()
        rt1 = map(float, rt1)
        rt2 = map(float, rt2)

        rtshift = [a - b for a, b in zip(rt1, rt2)]

        xlim_window = self.plot_boxplot_for_rt_delta_rt(rt1, rtshift)  # I think this is it.

        mz_ranges = [[400, 600], [600, 800], [800, 1000], [1000, 1200]]
        mz1, mz2 = self.GetMzList()
        mz1 = map(float, mz1)
        mz2 = map(float, mz2)

        rt_index = [[i for i in range(len(mz1)) if mz1[i] > interval[0] and mz1[i] <= interval[1]] for interval in
                    mz_ranges]
        rt1 = [[rt1[i] for i in rt_index[j]] for j in range(len(mz_ranges))]
        rtshift = [[rtshift[i] for i in rt_index[j]] for j in range(len(mz_ranges))]

        for i in range(len(mz_ranges)):
            self.plot_boxplot_for_rt_delta_rt(rt1[i], rtshift[i], mz_ranges[i], str(mz_ranges[i][1]), xlim_window)

            # exit(0)
            # pause

    def plot_boxplot_for_rt_delta_rt(self, rt1, rtshift, mz_range=[], mass_range_right="", xlim_windows=[0, 0]):
        name_of_font = "Arial"
        size_of_font = 12

        meanrtshift = np.mean(rtshift)
        sd_rtshift = np.std(rtshift)
        ylim_interval = [meanrtshift - 3 * sd_rtshift, meanrtshift + 3 * sd_rtshift]

        interval = 1000
        rt_shift_parts, rt_position = self.split_rt_shit_based_on_rt(rtshift, rt1, interval)
        # plt.plot(rt1, rtshift, 'g.', alpha=0.4)
        # plt.show()
        # plt.plot(rt1, rtshift, 'g.', alpha=0.4, label = self.type_name)

        commonfilename = self.inputfile.replace(".resu.gt.resu", "")
        commonfilename = commonfilename.replace(".resu", "")
        print commonfilename

        RTMatchesFile = commonfilename + ".dat.rt.matches"
        RTMatchesFile = commonfilename + ".dat1.dp.rt.matches"
        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]

        plt.plot(rt1, rtshift, alpha=0.7, label=self.type_name, marker='.', markersize=3.0, markerfacecolor='#808080',
                 markeredgecolor='#808080', linestyle='None')
        plt.plot(x, z, 'k-', markersize=1, label='DTW on MS1', linewidth=2)
        plt.legend(loc='best')
        bp = plt.boxplot(rt_shift_parts, positions=rt_position, widths=[int(interval / 1)] * len(rt_position),
                         showfliers=False, whis=[10, 90])
        plt.setp(bp['boxes'], color='black')
        plt.setp(bp['whiskers'], color='black')
        plt.setp(bp['medians'], color='black', alpha=0.7)

        if xlim_windows[1] == 0:
            xlim_windows = [0, rt_position[-1] + interval]
        plt.xlim(xlim_windows)
        xticks = range(0, 9001, 1000)

        plt.xticks(xticks, map(str, xticks), fontname=name_of_font, fontsize=size_of_font)
        ax = plt.gca()
        plt.setp(ax.get_xticklabels(), fontsize=size_of_font)
        # plt.tick_params(labelsize=size_of_font)
        # plt.xticks([0] + rt_position, ['0'] + map(int, map(int, rt_position)), fontname=name_of_font, fontsize=size_of_font)
        if len(mz_range) == 2:
            plt.title('$mz \in [%d, %d]$' % (mz_range[0], mz_range[1]))
        plt.xlabel('$RT_1$ (sec)')
        plt.ylabel('$\Delta RT = RT_1 - RT_2$ (sec)')
        plt.savefig(self.inputfile + mass_range_right + "_rtshift_boxplotZoomOut.png", dpi=360)
        zmax, zmax1, zmin, zmin1 = self.calculate_ylim_from_rt_matches()
        # plt.ylim(min([zmin, zmin1]), max([zmax, zmax1]))
        plt.ylim(ylim_interval)

        # plt.show()
        plt.savefig(self.inputfile + mass_range_right + "_rtshift_boxplotZoomIn.png", dpi=360)
        plt.close()
        return xlim_windows

    def calculate_ylim_from_rt_matches(self):
        commonfilename = self.inputfile.replace(".resu.gt.resu", "")
        commonfilename = commonfilename.replace(".resu", "")
        print commonfilename
        RTMatchesFile = commonfilename + ".dat.rt.matches"
        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]
        # plt.plot(x,z,'r-',markersize=2,label='DTW on MS1+MS2')
        std_times = 3.0
        zmin, zmax = np.mean(z) - np.std(z) * std_times, np.mean(z) + np.std(z) * std_times
        RTMatchesFile = commonfilename + ".dat1.dp.rt.matches"
        x, y = loadRTMatches(RTMatchesFile)
        z = [x[i] - y[i] for i in range(len(x))]
        zmin1, zmax1 = np.mean(z) - np.std(z) * std_times, np.mean(z) + np.std(z) * std_times
        return zmax, zmax1, zmin, zmin1

    def plot_scatterplot_for_intensity(self):
        # Todo: This function is not finished yet, I will implement this part in another function
        # plot several figures based on the alignment result
        # load data
        inten1, inten2 = self.GetIntenList()
        inten1 = map(float, inten1)
        inten2 = map(float, inten2)

        loginten1 = map(np.log2, inten1)
        loginten2 = map(np.log2, inten2)

    def GetSelfType(self):
        return self.type_name


def TestDraw():
    x = range(9);
    y = x;
    plt.plot(x, y);
    plt.show();


def GenerateReport(filename):
    import generateReport
    import os
    generateReport.CreateAlignmentReport(filename)
    inputpath = os.path.split(filename)[0]
    generateReport.GetReport(inputpath)


def DrawFiguresOfResuFile(inputfile):  # only lwbmatch file (.resu) not groundtruth file (.gt.resu)
    print '[Info] Plot figures of the alignment result'
    reader = LWBMatchResuReader(inputfile)
    reader.Read()
    reader.FilterZeros()
    reader.get_rt_shift_estimation_SD_on_MS1()
    reader.get_rt_shift_estimation_SD_on_MS1_MS2()
    reader.plot_boxplots_on_RT()  # Figure 2 b
    reader.DrawRTAlignCurve()
    reader.DrawMZAlignCurve()
    reader.DrawIntensityAlignCurve()
    # print 'Debug'
    reader.DrawRTShiftDistribution()
    reader.DrawRTshift_MS1()
    reader.DrawQuiverPlot()
    reader = GroundTruthReader(inputfile + '.gt.resu')
    reader.Read()
    reader.get_rt_shift_estimation_SD_on_MS1()
    reader.get_rt_shift_estimation_SD_on_MS1_MS2()
    # reader.FilterZeros()
    # reader.DrawRTAlignCurve()
    # reader.DrawMZAlignCurve()
    # reader.DrawIntensityAlignCurve()
    print 'Debug again------------------------'
    # exit(0)
    reader.DrawRTShiftDistribution()
    print 'Debug again and again'
    reader.DrawRTshift_MS1()
    # reader.DrawQuiverPlot()
    # reader.Print();

    GenerateReport(inputfile)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'Usage:'
        print sys.argv[0] + " lwbmatch_results"

        # hp = HistgramPloter([1.0,2,3,5,6.0,100],'mydata.png','test data')
        # hp.plot()

    else:
        for each_lwbmatch_resu_file in sys.argv[1:]:
            DrawFiguresOfResuFile(each_lwbmatch_resu_file)
