#!/usr/bin/python
# Python 2.7 Required
# Remove proteins if it has a name started with ">rev_"
# Input:
# Fasta file
#
#
# Output:
# Rewrite the same file. Remove the protein entries with ">rev_" in their name.
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------


import sys


# comment: get the input filename and rewrite with no decoy.

class FastaReader:
    def __init__(self, inputfile):
        self.filename = inputfile
        self.AC = []
        self.DE = []
        self.SQ = []
        self.isDecoy = []

    def printbyindex(self, idx):
        if idx < len(self.AC):
            print 'AC: %s\nDE: %s\nSQ: %s\n' % (self.AC[idx], self.DE[idx], self.SQ[idx])

    def lookforPep(self, peptideqry):
        pos = -1
        for i in range(len(self.AC)):
            pos = self.SQ[i].find(peptideqry)
            if pos == -1:
                continue
            else:
                print 'Found peptide %s in %d-th protein (0 ~ num-1)' % (peptideqry, i)
                print 'Protein name: %s' % (self.AC[i])
                break
                # end
        # end

        if pos == -1:
            print 'Fail to find peptide %s ' % (peptideqry)

        return pos

    def loadseq(self):
        # do not forget to load the data first.
        fid = open(self.filename, 'r')
        alllines = fid.readlines()
        fid.close()

        i = 0

        while i < len(alllines):
            curline = alllines[i].strip()
            # print curline
            if curline[0] == '>':
                items = curline.split()
                self.AC.append(items[0])
                self.DE.append(' '.join(items[1:]))
                self.SQ.append('')
                i = i + 1
            else:

                self.SQ[len(self.AC) - 1] = self.SQ[len(self.AC) - 1] + curline
                i = i + 1

                # end
        # end
        print 'Load %d proteins' % (len(self.AC))
        # self.printbyindex(len(self.AC)-1)

    def rewrite(self):

        fid = open(self.filename, 'r')
        alllines = fid.readlines()
        fid.close()

        fidout = open(self.filename, 'w')
        for line in alllines:
            if '>rev_' in line:  # find the reverse start, flip the flag
                reverse_sq = True
            elif '>' in line:  # find the non-reverse, flip the flag
                reverse_sq = False
                fidout.write(line)
            elif reverse_sq:  # if we are in the reverst protein, skip
                continue
            else:  # we input.
                fidout.write(line)
        # end for

        fidout.close()


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print sys.argv[0] + ' input.fasta'
        print 'This program could filter out the sequence started with >rev_'
    else:
        fr = FastaReader(sys.argv[1])
        fr.rewrite()
