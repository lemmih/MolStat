import matplotlib.pyplot as plt
import molecule as mol
import numpy as np
import math

no_generations = 100

def greedy_optimize(alkane_size):
    mol.generate_chain(alkane_size)
    no_dihedral = alkane_size-3
    dihedral_list = np.random.uniform(0.0, 360.0, no_dihedral)
    mol.set_dihedral(dihedral_list)
    energy = mol.get_energy()
    for _ in range(no_generations):
        new_list = np.random.uniform(0.0, 360.0, no_dihedral)
        mol.set_dihedral(new_list)
        new_energy = mol.get_energy()
        if new_energy < energy:
            dihedral_list = new_list
            energy = new_energy
    print 'Optimized configuration:', energy, dihedral_list

greedy_optimize(4)
greedy_optimize(10)
greedy_optimize(20)
greedy_optimize(40)
