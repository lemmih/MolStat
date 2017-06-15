import matplotlib.pyplot as plt
import molecule as mol
import numpy as np

no_atoms = 4
mol.generate_chain(no_atoms)

xaxis = range(0,180)
points = []
for i in xaxis:
    mol.set_dihedral([i])
    energy = mol.get_energy()
    # print i/10.0, energy
    points.append(energy)

plt.plot(xaxis, points, 'r')
plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Angle (degrees)')
plt.title('Energy of butane at different dihedral angles')
# plt.show()
plt.savefig('butane_energy.png')
