########################################################################################################################
#                                                                                                                      #
#   Script to plot from a COLVAR file (output of plumed run). It also pickles the COLVAR file for re-use.              #
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
#   Script developed by Naushad Velgy. Please send questions to naushad.velgy@dtc.ox.ac.uk                             #
#                                                                                                                      #
########################################################################################################################

import csv
import sys
import numpy
import os
from matplotlib import pyplot
import pickle

if os.path.exists('%s.p' % sys.argv[1]) == True:
    data = pickle.load(open('%s.p' % sys.argv[1], 'r'))
else:
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

if os.path.exists('%s.p' % sys.argv[1]) == False:
    pickle.dump(data, open('%s.p' % sys.argv[1], 'w'))

if str(sys.argv[2]).lower() == 'read':
    print 'Select from the following list:\n'
    for entry in data:
        print entry
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
    print 'plots field1 vs field2'

