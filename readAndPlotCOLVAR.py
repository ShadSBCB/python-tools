########################################################################################################################
#                                                                                                                      #
#   Script to plot from a COLVAR file (output of plumed run).                                                          #
#                                                                                                                      #
#   Usage:                                                                                                             #
#                                                                                                                      #
#   python readAndPlotCOLVAR.py -cv /path/to/COLVAR --read                                                             #
#                                                                                                                      #
#   returns all the fields available for plotting.                                                                     #
#                                                                                                                      #
#   python readAndPlotCOLVAR.py -cv /path/to/COLVAR --plot --field1 --field2                                           #
#                                                                                                                      #
#   Saves plots under the name field1_vs_field2_COLVAR.svg                                                             #
#                                                                                                                      #
#   python readAndPlotCOLVAR.py -cv /path/to/COLVAR --save --field1 --field2                                           #
#                                                                                                                      #
#   Saves data to csv file.                                                                                            #
#                                                                                                                      #
#   Script developed by Naushad Velgy. Please send questions to naushad.velgy@dtc.ox.ac.uk                             #
#                                                                                                                      #
########################################################################################################################

from argparse import ArgumentParser
import csv
import sys
import numpy
import os
from matplotlib import pyplot

################################################################################

parser = ArgumentParser(description="Script to parse and plot data from COLVAR files. Requires Matplotlib.")
parser.add_argument('-cv', '--colvar', dest='colvar', default=None, help="Path to COLVAR file. Make sure the file "
                                                                       "has only one '#!' header.")
parser.add_argument('--field1', dest='field1', default=None, help="1st field to plot or save.")
parser.add_argument('--field2', dest='field2', default=None, help="2nd field to plot or save.")
parser.add_argument('--read', dest='read', action='store_true',
                    help="Reads COLVAR file and prints fields to plot or save.")
parser.add_argument('--plot', dest='plot', action='store_true', help="Reads COLVAR file and saves field1_vs_field2.svg.")
parser.add_argument('--save', dest='save', action='store_true', help="Reads COLVAR file and saves field1_vs_field2.csv.")

args = parser.parse_args()
################################################################################

#Sanity check 1

if args.read == True and (args.field1 != None and args.field2 != None):
    print "The point of '--read' is to advise on which fields can be plotted."
    print "Giving me a '--read' command and a set of fields to plot defeats that purpose."
    print "Please, either use the '--plot' and/or '--save' commands to plot fields."
    sys.exit()
elif args.read == True and args.save == True:
    print "The point of '--read' is to advise on which fields can be plotted."
    print "Giving me a '--read' command and the '--save' command implies you know which fields you want to save."
    print "Please, either use the '--plot' and/or '--save' commands to plot fields."
    sys.exit()
elif args.read == True and args.plot == True:
    print "The point of '--read' is to advise on which fields can be plotted."
    print "Giving me a '--read' command and the '--plot' command implies you know which fields you want to plot."
    print "Please, either use the '--plot' and/or '--save' commands to plot fields."
    sys.exit()

#Sanity check 2
with open(args.colvar, 'r') as f:
    seen = set()
    for line in f:
        line_lower = line.lower()
        if line_lower in seen:
            print "It would seem that your header file has more than one header (starting with #!)."
            print "Ensure that it only has one and try again."
            print "This is important because the script will only get data after then last #! line"
            print "(which should be the first line of the file."
            sys.exit()
with open(args.colvar, 'r') as f:
    d = csv.reader(f, delimiter=' ')
    for line in d:
        if line[0] == '#!':
            entries = line[2:]
            data = {entries[i]:[] for i in range(len(entries))}
            continue
        else:
            for i,key in enumerate(entries):
                data[key].append(line[i+1])

if args.read == True:
    print "COLVAR file contains the following fields:\n"
    for entry in data:
        print entry
    sys.exit()

#Check if fields are in data

if (str(args.field1) not in data.keys()) and (str(args.field2) not in data.keys()):
    print "Neither field provided is not in the COLVAR file."
    print "Use /path/to/python /path/to/script.py -cv /path/to/COLVAR --read"
    print "to get a list of fields in the COLVAR file."
    sys.exit()
elif str(args.field1) not in data.keys():
    print "Field 1 is not in the COLVAR file."
    print "Use /path/to/python /path/to/script.py -cv /path/to/COLVAR --read"
    print "to get a list of fields in the COLVAR file."
    sys.exit()
elif str(args.field2) not in data.keys():
    print "Field 2 is not in the COLVAR file."
    print "Use /path/to/python /path/to/script.py -cv /path/to/COLVAR --read"
    print "to get a list of fields in the COLVAR file."
    sys.exit()

elif args.save == True:
    if (args.field1 != None and args.field2 != None):
        with open('%s_vs_%s_%s.csv' % (str(args.field1), str(args.field2), str(args.colvar)), 'w') as f:
            for entries in zip(data[str(args.field1)], data[str(args.field2)]):
                f.write('%s,%s\n' % (entries[0], entries[1]))
    else:
        print "Can't save fields if I have none to save."
        print "/path/to/python /path/to/script.py -cv /path/to/COLVAR --save --field1 field1 --field2 field2"
        print "is what you are looking for."

elif args.plot == True:
    if (args.field1 != None and args.field2 != None):
        print 'Plotting...'
        try:
            pyplot.plot(data[str(args.field1)], data[str(args.field2)], 'k.')
            pyplot.xlabel(str(args.field1))
            pyplot.ylabel(str(args.field2))
            plottitle = '%s_vs_%s.svg' % (str(args.field1), str(args.field2))
            pyplot.savefig('%s_vs_%s_%s.svg' % (str(args.field1), str(args.field2), str(args.colvar)))
        except:
            print 'Variable not in COLVAR file.\n'
            print 'Select from the following list:\n'
            for entry in data:
                print entry

else:
    print 'Option not recognised.'
    print 'Usage:'
    print '/path/to/python /path/to/script.py -cv /path/to/COLVAR --read'
    print 'shows options for plotting.'
    print '/path/to/python /path/to/script.py -cv /path/to/COLVAR --plot --field1 field1 --field2 field2'
    print 'plots field1 vs field2 in svg format.'
    print '/path/to/python /path/to/script.py -cv /path/to/COLVAR --save --field1 field1 --field2 field2'
    print 'saves field1 vs field2 to csv file.'
