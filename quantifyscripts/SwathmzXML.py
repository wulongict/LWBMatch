#!/usr/bin/python
# Read and split mzXML file into different files.
import sys


def CalculateSwathCycle(filename):
    Cycle = 0
    fid = open(filename, 'r')
    line = fid.readline()
    found = False
    while line != '':
        if 'msLevel=\"1\"' in line:
            if found:
                break
            else:
                Cycle += 1
                found = True
        elif 'msLevel=\"2\"' in line:
            Cycle += 1
        line = fid.readline()
    # end while
    fid.close()
    return Cycle


def SplitSwath(filename, i):
    Cycle = CalculateSwathCycle(filename)
    if i >= Cycle:
        return;
    print 'Total Cycle number ' + str(Cycle) + '. Spliting SWATH ' + str(i)
    fid = open(filename, 'r');
    fod = open(filename[0:-6] + '_SWATH' + str(i) + ".mzXML", "w")
    line = fid.readline()
    while line != '':
        if 'scanCount=' in line:  # scancount line
            scanCountstr = line.split('\"')[1]
            line = line.replace(scanCountstr, str(int(scanCountstr) * 2 / Cycle))
            fod.write(line)
            line = fid.readline()
        elif '<scan num=' in line:  # start line
            scanNumstr = line.split('\"')[1]
            if int(scanNumstr) % Cycle == 1 or (int(scanNumstr) - (i + 1)) % Cycle == 0:
                fod.write(line)
                line = fid.readline()
                while line != '' and '</scan>' not in line:
                    fod.write(line)
                    line = fid.readline()
                fod.write(line)  # write </scan>
                line = fid.readline()  # Read new line
            else:
                # Not the scan we want, skip all those scans, from <scan> to </scan>
                line = fid.readline()
                while line != '' and '</scan>' not in line:
                    line = fid.readline()
                line = fid.readline()  # Read new line

        else:  # other lines
            fod.write(line)
            line = fid.readline()  # Read new line

    fid.close()
    fod.close()


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'split SWATH mzXML file into several different file, each comes from one isolation window'
        print "Usage: %s SWATHData.mzXML startSwath endSwath\n e.g. \n %s abc.mzXML 1 1" % (sys.argv[0], sys.argv[0])
    else:
        filename = sys.argv[1]
        startSwath = int(sys.argv[2])
        endSwath = startSwath
        if len(sys.argv) == 4:
            endSwath = int(sys.argv[3])

        for i in range(startSwath, endSwath + 1):
            SplitSwath(filename, i)
