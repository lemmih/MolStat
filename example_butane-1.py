"""
  Example of calculating the energy
  of butane with a random dihedral angle
"""

import matplotlib.pyplot as plt
import molecule as mol
import numpy as np


#Molecule Specific
no_atoms = 4

#Its important to call this function before the others as it generates the global variable mol
mol.generate_chain(no_atoms)
# Number of dihedral angles corresponding
# to the molecule alkane chain
no_dihedral = 1

# Create a list of random dihedral angles
# between 0.0 and 360.0 degrees
dihedral_list = np.random.uniform(0.0, 360.0, no_dihedral)

mol.set_dihedral(dihedral_list)

energy = mol.get_energy()


mol.save_molecule('butane.xyz')

print "The energy for \omega = ", dihedral_list, "is", energy


