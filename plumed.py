########################################################################################################################
#                                                                                                                      #
#   Script to return a plumed.dat file with an all-atom contact map. Use the same protein in 2 different states.       #
#   Ensure the protein has the same number of atoms and the that all each atom on one state maps to the same atom      #
#   on the other state. Ensure all files are in GROMACS format (.gro).                                                 #
#                                                                                                                      #
#   Usage:                                                                                                             #
#                                                                                                                      #
#   python plumed.py /path/to/inactive/state.gro /path/to/active/state.gro                                             #
#                                                                                                                      #
#                                                                                                                      #
#   Script developed by Jan Domanski and Naushad Velgy. Please send questions to naushad.velgy@dtc.ox.ac.uk            #
#                                                                                                                      #
########################################################################################################################

#libraries to be used
from MDAnalysis import Universe as u
from MDAnalysis.lib import distances as cdists
import numpy as np
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

#Load receptors into universe
try:
    inactive = u(sys.argv[1])
    active = u(sys.argv[2])
except IndexError:
    print "Usage:"
    print "python plumed.py /path/to/inactive/state.gro /path/to/active/state.gro"
    sys.exit()

SCALING_FACTOR = 1.8

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

def generate_lines(gr, contacts):

    return [ "ATOMS%s=%s,%s SWITCH%s={EXP R_0=0.02000 D_0=%.4f} WEIGHT%s=%.6f\n" \
             % (i+1, gr[idA].id+1, gr[idB].id+1, i+1, distance * 0.1 * SCALING_FACTOR, i+1, 1.0/len(contacts)) \
            for i, ((idA, idB), distance) in enumerate(contacts.items())
    ]

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

with open('plumed.dat', 'w') as f:

    f.write("CONTACTMAP ...\n")
    f.write("".join(generate_lines(atoms_inactive, contacts_inactive)))
    f.write("LABEL=cmap_inactive\n")
    f.write("SUM\n")
    f.write("... CONTACTMAP\n\n")

    f.write("CONTACTMAP ...\n")
    f.write("".join(generate_lines(atoms_active, contacts_active)))
    f.write("LABEL=cmap_active\n")
    f.write("SUM\n")
    f.write("... CONTACTMAP\n")

    f.write("cmap: COMBINE ARG=cmap_inactive,cmap_active COEFFICIENTS=1,-1 PERIODIC=NO\n")

    f.write("PRINT ARG=* FILE=COLVAR")
