#!/usr/bin/python
# Python 2.7 Required
# plot figures of the dot product matrix
# Input:
# dp.txt files contains all the dot product between a pair of input files
# Output:
# Figures for each dp.txt file
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------



import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import copy
import multiprocessing


class PlotDTWRTAlignCurve:
    def __init__(self):
        pass


class PlotDotProductImage:
    def __init__(self):
        pass


def print_usage():
    print 'Usage:'
    print 'plotMatrix.py inputmatrix[tsv]'



def plot_matrix_image(inputdata, outputname):
    """
    :type l: multiprocessing.Lock()
    :param l: a lock
    """

    plt.imshow(inputdata, interpolation='nearest')
    plt.colorbar()

    plt.savefig(outputname, dpi=180)

    plt.close()

    print outputname


def load_matrix_to_list(inputfile):
    print "[Info] loading data from " + inputfile
    fid = open(inputfile, 'r')
    lines = fid.readlines()
    fid.close()
    lines = [lines[i].strip() for i in range(len(lines))]
    lines = [lines[i].split('\t') for i in range(len(lines))]

    lines = [[float(lines[i][j]) for j in range(len(lines[i]))] for i in range(len(lines))]

    # for i in range(5):
    #     print lines[i][1:10]
    return lines


def ReadPlotSumFile(filename):
    lines = load_matrix_to_list(filename)
    dotproduct = np.array(lines)

    plot_matrix_image(dotproduct, filename + '.png')


def plot_msn_dp_with_multi_thread(ms1_dp_txt_filename, cyclenum):
    lines = load_matrix_to_list(ms1_dp_txt_filename)
    dotproduct = np.array(lines)

    # plot MS1 dp figure
    lock = multiprocessing.Lock()
    plot_matrix_image(dotproduct, ms1_dp_txt_filename + ".png")

    # initialize the sum of dp matrix
    sum_of_dp = copy.deepcopy(dotproduct)
    sum_of_dp = sum_of_dp * 2.0 / cyclenum;

    # mythreads = []
    # l = multiprocessing.Lock()
    for i in range(2, cyclenum):
        print '[Info] plot %s-th SWATH ' % str(i - 1)
        # for j in range(5):
        #     print sum_of_dp[j][0:10]

        msi_dp_txt_filename = ms1_dp_txt_filename.replace("dat1", "dat" + str(i));
        lines = load_matrix_to_list(msi_dp_txt_filename)
        current_dp = np.array(lines)
        sum_of_dp = sum_of_dp + current_dp / cyclenum
        # t = multiprocessing.Process(target=plot_matrix_image, args=(l,current_dp,ms1_dp_txt_filename+'.png',))
        # t.start()
        # mythreads.append(t)
        if os.path.isfile(msi_dp_txt_filename + '.png'):
            continue
        plot_matrix_image(current_dp, msi_dp_txt_filename + ".png")

    # for i in range(2):
    #     mythreads[i].start()


    # for j in range(5):
    #     print sum_of_dp[j][0:10]

    print '[Info] plot the sum of ms1 and ms2 dot product matrix...'
    if not os.path.isfile(ms1_dp_txt_filename + 'same.png'):
        plot_matrix_image(sum_of_dp, ms1_dp_txt_filename + "sum.png")


def plot_dtw_rt_matches(filename):
    """
    :type filename: str
    """
    if os.path.isfile(filename + '.shift.png'):
        print '[Info] figure exist!'
        return

    print '[Info] plot rt matches curve for file %s ' % filename


    x, y, z = load_rt_matches_data(filename)


    plt.plot(x, z, 'g.', markersize=1)
    plt.xlabel('RT1(sec)')
    plt.ylabel('RT2 - RT1 by DTW (sec)')
    plt.savefig(filename + '.shift.png', dpi=180)
    plt.close()


def load_rt_matches_data(filename):
    """
    :type filename: str
    """
    print '[Info] loading rt matches file: %s' % filename
    fid = open(filename, "r")
    lines = fid.readlines()
    fid.close()
    lines = [eachline.strip() for eachline in lines]
    lines = [eachline.split() for eachline in lines]
    lines = [[float(item) for item in eachline] for eachline in lines]
    # print lines[0]
    x = [eachline[0] for eachline in lines]
    y = [eachline[1] for eachline in lines]
    z = [x[i] - y[i] for i in range(len(x))]
    return x, y, z


def plot_each_rt_matches_curve_with_multi_thread(filename, cyclenum):
    commonbase = filename.replace(".rt.matches", "")
    print '[Info] plot the rt matches of summation of ms1 and ms2'
    plot_dtw_rt_matches(filename)

    for i in range(cyclenum):
        print '[Info] plot the rt matches of each single rt matches file'
        inputfile = commonbase + str(i + 1) + ".dp.rt.matches"
        # t = threading.Thread(target=plot_dtw_rt_matches, args=(inputfile,))
        # t.start()
        plot_dtw_rt_matches(commonbase + str(i + 1) + ".dp.rt.matches")


def plot_all_rt_matches_curve_into_single_file(filename, cyclenum):
    """
    :type filename: str
    :type cyclenum: int
    """
    if os.path.isfile(filename + '.same.jpg'):
        print '[Info] file exists!'
        return

    commonbase = filename.replace(".rt.matches", "")

    # plt.savefig(filename+'.shift.png',dpi = 720)
    number = cyclenum
    cmap = plt.get_cmap('gnuplot')  # ('rainbow')#('gnuplot')
    colors = [cmap(i) for i in np.linspace(0, 1, number)]
    for i in range(1, cyclenum):
        curfilename = commonbase + str(i + 1) + ".dp.rt.matches"
        x, y, z = load_rt_matches_data(curfilename)
        plt.plot(x, z, color=colors[i], markersize=1, linewidth=0.8, label='${i}$'.format(i=i + 1))

    # end

    curfilename = commonbase + str(0 + 1) + ".dp.rt.matches"
    x, y, z = load_rt_matches_data(curfilename)
    plt.plot(x, z, markersize=1, linewidth=1, label='$ 1$')

    x, y, z = load_rt_matches_data(filename)
    plt.plot(x, z, 'g-', markersize=1, linewidth=1, label='$ all$')

    ax = plt.gca()
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.25,
                     box.width, box.height * 0.75])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
              fancybox=True, shadow=True, ncol=9, prop={'size': 8})
    # plt.legend(frameon=False)

    # plt.legend(prop={'size':6},loc='best')

    plt.xlabel('$RT_1$ (sec)')
    plt.ylabel('$\Delta RT = RT_1-RT_2$ (sec)')
    # plt.savefig(filename + ".same.png", dpi=180)
    plt.savefig(filename + ".same.jpg", dpi=180)

    plt.close()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print_usage()

    else:
        filename = sys.argv[1]
        plot_each_rt_matches_curve_with_multi_thread(filename, 35)
        plot_all_rt_matches_curve_into_single_file(filename, 35)

