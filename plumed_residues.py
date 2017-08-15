########################################################################################################################
#                                                                                                                      #
#   Script to return a text file with a list of residue-by-residue contact, a file to visualise the residues in        #
#   PyMOL, and a file to visualise the residues in VMD. Use the same protein in 2 different states.                    #
#   Ensure the protein has the same number of atoms and the that all each atom on one state maps to the same atom      #
#   on the other state. Ensure all files are in GROMACS format (.gro).                                                 #
#                                                                                                                      #
#   Usage:                                                                                                             #
#                                                                                                                      #
#   python plumed.py --input1 /path/to/inactive/state.gro --input2 /path/to/active/state.gro --gpcr BW_numbers.txt     #
#                                                                                                                      #
#   The last file is optional and only works with GPCR structures.                                                     #                                                             #
#                                                                                                                      #
#   Script developed by Jan Domanski and Naushad Velgy. Please send questions to naushad.velgy@dtc.ox.ac.uk            #
#                                                                                                                      #
########################################################################################################################

from MDAnalysis import Universe
from MDAnalysis.lib import distances as cdists
import numpy as np
import csv
import re
import os
import sys

#Map to nearest non-symmetrical atom (often C, also N and O and S)

selection = """
protein and not name H* and (
   (resname ALA and (backbone or name CB))
or (resname ARG and (backbone or name CB or name CG or name NE))
or (resname ASN and (backbone or name CG))
or (resname ASP and (backbone or name CG))
or (resname CYS and (backbone or name SG))
or (resname GLN and (backbone or name CB or name CD))
or (resname GLU and (backbone or name CB or name CD))
or (resname GLY)
or (resname HIS and (backbone or name CE1 or name NE2))
or (resname ILE and (backbone or name CB))
or (resname LEU and (backbone or name CG))
or (resname LYS and (backbone or name CB or name CZ or name NZ))
or (resname MET and (backbone or name CB or name SD))
or (resname SER and (backbone or name OG))
or (resname PHE and (backbone or name CG or name CZ))
or (resname PRO and (backbone or name CD))
or (resname THR and (backbone or name CB))
or (resname TRP and (name CE3 or name NE1))
or (resname TYR and (backbone or name CG or name CZ))
or (resname VAL and (backbone or name CB))
)
"""

################################################################################

parser = ArgumentParser(description="Generates a text file with a list of unique residues, as well as  files for "
                                    "visualisation. Requires MDAnalysis and numpy.")
parser.add_argument('--input1', dest=struct1, default=None, help="Input structure 1, in GRO format.")
parser.add_argument('--input2', dest=struct2, default=None, help="Input structure 2, in GRO format.")
parser.add_argument('--scaling_factor', dest=SF, default=1.8, help="Scaling factor. Default is 1.8.")
parser.add_argument('--gpcr', dest=gpcr, default=None, help="Input text with Ballesteros-Weinstein numbers."
                                                            "See notes for more details."
                                                            "If your protein is not a GPCR, do not use this flag.")
parser.add_argument('--pymol', dest=pymol, action='store_true', help="Outputs a PML file to be used in PyMOL."
                                                                     "See notes for more details.")
parser.add_argument('--vmd', dest=vmd, action='store_true', help="Outputs a TCL file to be used in VMD."
                                                                     "See notes for more details.")
parser.add_argument('--help', dest=help, action='store_true', help="Print this page and exit.")

args = parser.parse_args()
################################################################################

if args.help == True:
    print "Usage:"
    print "/path/to/python /path/to/plumed_residues.py --input1 /path/to/state1.gro --input2 /path/to/state2.gro"
    print "--input1 and --input2 are 2 receptor files in different conformations. \n" \
          "Make sure atoms are in the same order."
    print "--scaling_factor (optional) determines the scaling factor. The default is 1.8."
    print "--gpcr (optional) signals script to output residue list with Ballesteros-Weinstein numbers. See notes."
    print "--pymol (optional) signals script to output a PML file for visualisation of contacts in PyMOL. See notes."
    print "--vmd (optional) signals script to output a TCL file for visualisation of contacts in VMD. See notes."
    print "--help prints this page and exits."
    sys.exit()

try:
    inactive = Universe(args.struct1)
    active = Universe(args.struct2)
    #Sanity check 1
    if inactive.filename.split('.')[1] != 'gro' and active.filename.split('.')[1] != 'gro':
        print "This scripts works better if both inputs are in the GRO format."
        print "This is because the nomenclature of GRO files is consistent, which is important."
        print "Consider converting your files to GRO files using, for example,"
        print "gmx editconf -f input.notgro -o input.gro"
        sys.exit()
except IndexError:
    print "Usage:"
    print "python plumed.py --input1 /path/to/inactive/state.gro --input2 /path/to/active/state.gro"
    sys.exit()

#Sanity check 2 (very important)
for i, atoms in enumerate(zip(inactive.select_atoms('protein').atoms.names,
                              active.select_atoms('protein').atoms.names)):
    if atoms[0] == atoms[1]:
        continue
    else:
        print "Atom {} is different in the two structures.".format(i)

SCALING_FACTOR = args.SF

def generate_contacts(
        universe,
        selection,
        upper_cutoff=6.0,
        lower_cutoff=0.0,
        resid_separation=3
    ):

    gr = universe.select_atoms(selection)
    distances = cdists.distance_array(gr.positions, gr.positions)

    triu = np.triu(distances)
    mask = np.argwhere((triu < upper_cutoff) & (triu > lower_cutoff))

    contacts = [(atom1, atom2) for atom1, atom2 in mask \
        if not abs(gr.atoms[atom1].resid - gr.atoms[atom2].resid) <= resid_separation]
    return gr, triu, set([ frozenset(pair) for pair in contacts ])

def get_residues(gr, contacts):

    return [ [gr[key[0]].resid, gr[key[1]].resid] for key in contacts.keys() ]

def get_BW_numbers(text):
    BW_numbers = {re.compile("([a-zA-Z]+)([0-9]+)").match(line[2].strip()).group(2): line[1].strip() for line in
                  csv.reader(open(text, 'r'), delimiter='\t') }
    return BW_numbers

atoms_active, triu_active, contacts_active = generate_contacts(active, selection)
atoms_inactive, triu_inactive, contacts_inactive = generate_contacts(inactive, selection)

unique_active_contacts = contacts_active - contacts_inactive
unique_inactive_contacts = contacts_inactive - contacts_active

contacts_active = {}
for idA, idB in unique_active_contacts:
    da = triu_active[idA, idB] or triu_active[idB, idA]
    di = triu_inactive[idA, idB] or triu_inactive[idB, idA]
    if min(da, di) * SCALING_FACTOR > max(da, di): continue
    #print("active", da, di)
    contacts_active[(idA, idB)] = da

contacts_inactive = {}
for idA, idB in unique_inactive_contacts:
    da = triu_active[idA, idB] or triu_active[idB, idA]
    di = triu_inactive[idA, idB] or triu_inactive[idB, idA]
    if min(da, di) * SCALING_FACTOR > max(da, di): continue
    #print("inactive", da, di)
    contacts_inactive[(idA, idB)] = di

all_res_inactive_list = get_residues(atoms_inactive, contacts_inactive)
inactive_res_list = []
a = [ inactive_res_list.append(elem) for elem in all_res_inactive_list if elem not in inactive_res_list ]

all_res_active_list = get_residues(atoms_active, contacts_active)
active_res_list = []
b = [ active_res_list.append(elem) for elem in all_res_active_list if elem not in active_res_list ]

if args.gpcr == True:
    try:
        BW_numbers = get_BW_numbers(args.gpcr)
        with open('Unique_contact_residues.dat', 'w') as f:
            writer = csv.writer(f)
            f.write('Inactive\n')
            f.write('Resid1 \t BW_number \t Resid2 \t BW_number \n')
                for entry in inactive_res_list:
                    try:
                        f.write('%s \t %s \t %s \t %s \n' % (entry[0], BW_numbers[str(entry[0])], entry[1],
                                                             BW_numbers[str(entry[1])]))
                except KeyError:
                    continue

            f.write('Active\n')
            f.write('Resid1 \t Resid2\n')
            for entry in active_res_list:
                try:
                    f.write('%s \t %s \t %s \t %s \n' % (entry[0], BW_numbers[str(entry[0])], entry[1],
                                                         BW_numbers[str(entry[1])]))
                except KeyError:
                    continue

    except IndexError:
        print 'BW numbers file not given, or doesn\'t have the correct format.'
        with open('Unique_contact_residues.dat', 'w') as f:
            writer = csv.writer(f)
            f.write('Inactive\n')
            f.write('Resid1 \t Resid2\n')
            for entry in inactive_res_list:
                f.write('%s \t %s \n' % (entry[0], entry[1]))
            f.write('Active\n')
            f.write('Resid1 \t Resid2\n')
            for entry in active_res_list:
                f.write('%s \t %s \n' % (entry[0], entry[1]))

if args.pymol == True:
    with open('PyMOL_representation.pml', 'w') as f:
        for entry in get_residues(atoms_inactive, contacts_inactive):
            f.write('show sticks, resid %s or resid %s and not name H*\n' % (entry[0], entry[1]))
        for entry in get_residues(atoms_active, contacts_active):
            f.write('show sticks, resid %s or resid %s and not name H*\n' % (entry[0], entry[1]))
        f.write('hide (hydro)')

if args.vmd == True:
    with open('VMD_representation.tcl', 'w') as f:
        for entry in get_residues(atoms_inactive, contacts_inactive):
            f.write('mol representation Licorice\n')
            f.write('mol selection "noh(resid %s or resid %s)"\n' % (entry[0], entry[1]))
            f.write('mol addrep top\n')
        for entry in get_residues(atoms_active, contacts_active):
            f.write('mol representation Licorice\n')
            f.write('mol selection "noh(resid %s or resid %s)"\n' % (entry[0], entry[1]))
            f.write('mol addrep top\n')
