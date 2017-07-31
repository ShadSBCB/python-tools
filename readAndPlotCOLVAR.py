########################################################################################################################
#                                                                                                                      #
#   Script to plot from a COLVAR file (output of plumed run).                                                          #
#                                                                                                                      #
#   Usage:                                                                                                             #
#                                                                                                                      #
#   python readAndPlotCOLVAR.py /path/to/COLVAR read                                                                   #
#                                                                                                                      #
#   returns all the fields available for plotting.                                                                     #
#                                                                                                                      #
#   python readAndPlotCOLVAR.py /path/to/COLVAR plot field1 field2                                                     #
#                                                                                                                      #
#   Saves plots under the name field1_vs_field2_COLVAR.svg                                                             #
#                                                                                                                      #
#   python readAndPlotCOLVAR.py /path/to/COLVAR save field1 field2                                                     #
#                                                                                                                      #
#   Saves data to csv file.                                                                                            #
#                                                                                                                      #
#   Script developed by Naushad Velgy. Please send questions to naushad.velgy@dtc.ox.ac.uk                             #
#                                                                                                                      #
########################################################################################################################

import csv
import sys
import numpy
import os
from matplotlib import pyplot

with open(sys.argv[1], 'r') as f:
    d = csv.reader(f, delimiter=' ')
    for line in d:
        if line[0] == '#!':
            entries = line[2:]
            data = {entries[i]:[] for i in range(len(entries))}
            continue
        else:
            for i,key in enumerate(entries):
                data[key].append(line[i+1])

if str(sys.argv[2]).lower() == 'read':
    print 'Select from the following list:\n'
    for entry in data:
        print entry
elif str(sys.argv[2]).lower() == 'save':
    with open('%s_vs_%s_%s.csv' % (str(sys.argv[3]), str(sys.argv[4]), str(sys.argv[1])), 'w') as f:
        for entries in zip(data[str(sys.argv[3])], data[str(sys.argv[4])]):
            f.write('%s,%s\n' % (entries[0], entries[1]))
elif str(sys.argv[2]).lower() == 'plot':
    print 'Plotting...'
    try:
        pyplot.plot(data[str(sys.argv[3])], data[str(sys.argv[4])], 'k.')
        pyplot.xlabel(str(sys.argv[3]))
        pyplot.ylabel(str(sys.argv[4]))
        plottitle = '%s_vs_%s.svg' % (str(sys.argv[3]), str(sys.argv[4]))
        pyplot.savefig('%s_vs_%s_%s.svg' % (str(sys.argv[3]), str(sys.argv[4]), str(sys.argv[1])))
    except:
        print 'Variable not in COLVAR file.\n'
        print 'Select from the following list:\n'
        for entry in data:
            print entry
else:
    print 'Option not recognised.'
    print 'Usage:'
    print '/path/to/python /path/to/script.py /path/to/COLVAR read'
    print 'shows options for plotting.'
    print '/path/to/python /path/to/script.py /path/to/COLVAR plot field1 field2'
    print 'plots field1 vs field2 in svg format.'
    print '/path/to/python /path/to/script.py /path/to/COLVAR save field1 field2'
    print 'saves field1 vs field2 to csv file.'
