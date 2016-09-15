#!/usr/bin/python
# Python 2.7 Required
# Recall and precision reader,
# Input:
# Folder contains recall and precision file.
# Output:
# A CSV file contains all the recall precision values
#
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------

""""""
import os
import sys
import glob


class RecallPrecionReader:
    def __init__(self):
        print '[Info] Reading recall precision file...'

    def run(self, inputpath):
        file_list = self.get_file_list(inputpath)
        recall_precision = self.read_files(file_list)
        self.output(inputpath, recall_precision)

    def output(self, inputpath, recall_precision):
        outputfile = os.path.join(inputpath, 'rp.csv')
        fid = open(outputfile, 'w')
        for eachitem in recall_precision:
            fid.write(','.join(eachitem))
            fid.write('\n')

        fid.close()
        outputfile = os.path.join(inputpath, 'rp.tsv')
        fid = open(outputfile, 'w')
        for eachitem in recall_precision:
            fid.write('\t'.join(eachitem))
            fid.write('\n')

        fid.close()

    def get_file_list(self, inputpath):
        file_list = glob.glob(os.path.join(inputpath, '*recall*precision*'))
        print '[Info] find file '
        print  file_list
        return file_list

    def read_files(self, file_list):
        recall_precision = []
        for eachfile in file_list:
            data_list = self.read_one_file(eachfile)
            recall_precision.extend(data_list)
        # end
        return recall_precision

    def read_one_file(self, recall_precision_file):
        fid = open(recall_precision_file, 'r')
        all_lines = fid.readlines()

        fid.close()
        striped_lines = [eachline.strip() for eachline in all_lines]

        recall_precision_data = []
        for i in range(0, len(all_lines), 3):
            print all_lines[i:i + 3]
            filename = os.path.split(all_lines[i])[1][:-27]
            precision = all_lines[i + 1].split(':')[1].strip()
            recall = all_lines[i + 2].split(':')[1].strip()
            recall_precision_data.append([filename, recall, precision])

        print recall_precision_data
        return recall_precision_data


if __name__ == '__main__':
    rpReader = RecallPrecionReader()
    if len(sys.argv) == 1:
        print 'Usage:\n%s inputpath' % sys.argv[0]

    else:
        inputpath = sys.argv[1]
        rpReader.run(inputpath)
