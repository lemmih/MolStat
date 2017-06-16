import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import molecule as mol
import numpy as np
import math
import sys

def reset():
    global no_generations
    global no_parents
    global temperature
    global mutation_rate

    no_generations = 100
    no_parents = 20
    temperature = 1.0
    mutation_rate = 0.1

reset()

def configuration_energy(dihedral_list):
    mol.set_dihedral(dihedral_list)
    return mol.get_energy()




def greedy_optimize(alkane_size):
    mol.generate_chain(alkane_size)
    no_dihedral = alkane_size-3
    dihedral_list = np.random.uniform(0.0, 360.0, no_dihedral)
    mol.set_dihedral(dihedral_list)
    energy = mol.get_energy()
    lows = []
    for _ in range(no_generations):
        new_list = np.random.uniform(0.0, 360.0, no_dihedral)
        mol.set_dihedral(new_list)
        new_energy = mol.get_energy()
        if new_energy < energy:
            dihedral_list = new_list
            energy = new_energy
        lows.append(energy)
    return lows





def mate_mutate(parent_alpha, parent_beta):
    n = len(parent_alpha)
    m = np.random.randint(0,n)
    left = np.concatenate([parent_alpha[0:m], parent_beta[m:n]])
    right = np.concatenate([parent_beta[0:m], parent_alpha[m:n]])
    if np.random.random() < mutation_rate:
        left[np.random.randint(0,n)] = np.random.uniform(0.0, 360.0)
    if np.random.random() < mutation_rate:
        right[np.random.randint(0,n)] = np.random.uniform(0.0, 360.0)
    return (left, right)

def mate_mutate_improved(parent_alpha, parent_beta):
    n = len(parent_alpha)
    left = []
    right = []
    for i in range(n):
        if np.random.random() < 0.5:
            left.append(parent_alpha[i])
            right.append(parent_beta[i])
        else:
            right.append(parent_alpha[i])
            left.append(parent_beta[i])

    if np.random.random() < mutation_rate:
        left[np.random.randint(0,n)] = np.random.uniform(0.0, 360.0)
    if np.random.random() < mutation_rate:
        right[np.random.randint(0,n)] = np.random.uniform(0.0, 360.0)
    return (left, right)


def newVector(no_dihedral):
    return np.random.uniform(0.0, 360.0, no_dihedral)

def parenticide(child_energy, parent_energy):
    delta = abs(child_energy-parent_energy)
    return np.random.random() < math.exp(-delta/temperature)

def genetic_optimize(alkane_size):
    mol.generate_chain(alkane_size)
    no_dihedral = alkane_size-3

    parents = [ newVector(no_dihedral) for _ in range(no_parents) ]
    energy = [ configuration_energy(parent) for parent in parents ]

    lows = []
    means = []

    converged = no_generations

    for gen in range(no_generations):
        print 'Gen:', gen, 'size', alkane_size, 'temp', temperature, 'parents', no_parents, 'mut', mutation_rate
        # Each parent mates exactly once
        for alpha in range(0,no_parents,2):
            # mate alpha with random other.
            # then move random other to alpha+1
            beta = np.random.randint(alpha+1,no_parents)
            (left, right) = mate_mutate(parents[alpha], parents[beta])
            left_energy = configuration_energy(left)
            right_energy = configuration_energy(right)

            if left_energy < energy[alpha] or parenticide(left_energy, energy[alpha]):
                energy[alpha] = left_energy
                parents[alpha] = left
                # print 'New energy', left_energy
            elif left_energy < energy[beta] or parenticide(left_energy, energy[beta]):
                energy[beta] = left_energy
                parents[beta] = left
                # print 'New energy', left_energy
            if right_energy < energy[alpha] or parenticide(right_energy, energy[alpha]):
                energy[alpha] = right_energy
                parents[alpha] = right
                # print 'New energy', right_energy
            elif right_energy < energy[beta] or parenticide(right_energy, energy[beta]):
                energy[beta] = right_energy
                parents[beta] = right
                # print 'New energy', right_energy

            parents[alpha+1], parents[beta] = parents[beta], parents[alpha+1]
            energy[alpha+1], energy[beta] = energy[beta], energy[alpha+1]

        lows.append(min(energy))
        means.append(np.mean(energy))
        if min(energy) / np.mean(energy) > 0.95 and converged > gen:
            converged = gen
    return (lows, means, converged)

def plot_gens(alkane_size):
    reset()
    print 'Plot by gen'

    (lows, means, converged) = genetic_optimize(alkane_size)
    plt.plot(range(no_generations), lows, 'r', label='Minimum')
    plt.plot(range(no_generations), means, 'k', label='Mean')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('Generation')
    plt.yscale('log')
    # plt.text(50, max(means)/2, 'Final configuration energy:',
    #     verticalalignment='top', horizontalalignment='center',
    #     color='green', fontsize=15)
    plt.annotate('Final configuration energy: '+str(lows[-1]), xy=(0.3, 0.5), xytext=(0.3, 0.5), xycoords='figure fraction', textcoords='figure fraction')
    plt.title('C'+str(alkane_size)+'H'+str(alkane_size*2+2) + ' energy')
    plt.suptitle('temp=' + str(temperature) + ', mut_rate=' + str(mutation_rate) + ', parents=' + str(no_parents))
    low_patch = mpatches.Patch(color='red', label='Minimum')
    mean_patch = mpatches.Patch(color='black', label='Mean')
    plt.legend(handles=[low_patch, mean_patch])
    # plt.show()
    plt.savefig('by_gen_n_'+str(alkane_size)+'_t_'+str(temperature)+'_m_'+str(mutation_rate)+'_p_'+str(no_parents)+'.png')
    plt.clf()

def plot_temp(alkane_size):
    global temperature
    reset()

    print 'Plot temperature', alkane_size

    reset()
    xaxis = [i/10.0 for i in range(5,105,5)]
    lows = []
    means = []
    for t in xaxis:
        temperature = t
        (l, m, converged) = genetic_optimize(alkane_size)
        lows.append(l[-1])
        means.append(m[-1])
    plt.plot(xaxis, lows, 'r', label='Minimum')
    plt.plot(xaxis, means, 'k', label='Mean')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('Temperature')
    plt.title('C'+str(alkane_size)+'H'+str(alkane_size*2+2) + ' energy')
    plt.suptitle('mut_rate=' + str(mutation_rate) + ', parents=' + str(no_parents))
    low_patch = mpatches.Patch(color='red', label='Minimum')
    mean_patch = mpatches.Patch(color='black', label='Mean')
    plt.legend(handles=[low_patch, mean_patch])
    # plt.show()
    plt.savefig('by_temp_n_'+str(alkane_size)+'_m_'+str(mutation_rate)+'_p_'+str(no_parents)+'.png')
    plt.clf()

def plot_mut(alkane_size):
    global mutation_rate

    print 'Plot mutation rate', alkane_size

    reset()
    xaxis = [i/100.0 for i in range(0,55,5)]
    lows = []
    means = []
    for rate in xaxis:
        mutation_rate = rate
        (l, m, converged) = genetic_optimize(alkane_size)
        lows.append(l[-1])
        means.append(m[-1])
    plt.plot(xaxis, lows, 'r', label='Minimum')
    plt.plot(xaxis, means, 'k', label='Mean')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('Mutation rate')
    plt.title('C'+str(alkane_size)+'H'+str(alkane_size*2+2) + ' energy')
    plt.suptitle('temp=' + str(temperature) + ', parents=' + str(no_parents))
    low_patch = mpatches.Patch(color='red', label='Minimum')
    mean_patch = mpatches.Patch(color='black', label='Mean')
    plt.legend(handles=[low_patch, mean_patch])
    # plt.show()
    plt.savefig('by_mut_n_'+str(alkane_size)+'_t_'+str(temperature)+'_p_'+str(no_parents)+'.png')
    plt.clf()

def plot_greedy(alkane_size):
    print 'Plot vs greedy'

    reset()

    greedy_lows = greedy_optimize(alkane_size)

    (lows, means, converged) = genetic_optimize(alkane_size)
    plt.plot(range(no_generations), lows, 'r', label='Minimum')
    plt.plot(range(no_generations), means, 'k', label='Mean')
    plt.plot(range(no_generations), greedy_lows, 'g', label='Greedy')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('Generation')
    plt.yscale('log')
    plt.title('C'+str(alkane_size)+'H'+str(alkane_size*2+2) + ' energy')
    plt.suptitle('temp=' + str(temperature) + ', mut_rate=' + str(mutation_rate) + ', parents=' + str(no_parents))
    low_patch = mpatches.Patch(color='red', label='Minimum')
    mean_patch = mpatches.Patch(color='black', label='Mean')
    greedy_patch = mpatches.Patch(color='green', label='Greedy')
    plt.legend(handles=[low_patch, mean_patch, greedy_patch])
    # plt.show()
    plt.savefig('by_greedy_n_'+str(alkane_size)+'_t_'+str(temperature)+'_m_'+str(mutation_rate)+'_p_'+str(no_parents)+'.png')
    plt.clf()

def plot_parents(alkane_size):
    global no_parents

    reset()

    print 'Plot parents', alkane_size

    reset()
    xaxis = range(4,100,10)
    converges = []
    for p in xaxis:
        no_parents = p
        (l, m, converged) = genetic_optimize(alkane_size)
        print converged
        converges.append(converged)
    plt.plot(xaxis, converges, 'r')
    plt.ylabel('Generations to converge')
    plt.xlabel('Number of parents')
    plt.title('C'+str(alkane_size)+'H'+str(alkane_size*2+2))
    plt.suptitle('temp=' + str(temperature) + ', mut_rate=' + str(mutation_rate))

    plt.savefig('by_parents_n_'+str(alkane_size)+'_t_'+str(temperature)+'_m_'+str(mutation_rate)+'_p_'+'.png')
    plt.clf()

def plot_gens_improved(alkane_size):
    global mate_mutate
    reset()
    tmp = mate_mutate
    mate_mutate = mate_mutate_improved

    print 'Plot by gen improved'

    (lows, means, converged) = genetic_optimize(alkane_size)
    mate_mutate = tmp

    plt.plot(range(no_generations), lows, 'r', label='Minimum')
    plt.plot(range(no_generations), means, 'k', label='Mean')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('Generation')
    plt.yscale('log')
    plt.annotate('Final configuration energy: '+str(lows[-1]), xy=(0.3, 0.5), xytext=(0.3, 0.5), xycoords='figure fraction', textcoords='figure fraction')
    plt.title('C'+str(alkane_size)+'H'+str(alkane_size*2+2) + ' energy')
    plt.suptitle('temp=' + str(temperature) + ', mut_rate=' + str(mutation_rate) + ', parents=' + str(no_parents))
    low_patch = mpatches.Patch(color='red', label='Minimum')
    mean_patch = mpatches.Patch(color='black', label='Mean')
    plt.legend(handles=[low_patch, mean_patch])
    # plt.show()
    plt.savefig('by_gen2_n_'+str(alkane_size)+'_t_'+str(temperature)+'_m_'+str(mutation_rate)+'_p_'+str(no_parents)+'.png')
    plt.clf()

sims = {'1': plot_gens
       ,'2': plot_temp
       ,'3': plot_mut
       ,'4': plot_greedy
       ,'5': plot_parents
       ,'6': plot_gens_improved }


if len(sys.argv) == 1:
    for i in sims:
        sims[i](10)
        sims[i](20)
        sims[i](40)
elif len(sys.argv) == 2:
    sims[sys.argv[1]](10)
    sims[sys.argv[1]](20)
    sims[sys.argv[1]](40)
elif len(sys.argv) == 3:
    if int(sys.argv[2]) == 1:
        sims[sys.argv[1]](10)
    elif int(sys.argv[2]) == 2:
        sims[sys.argv[1]](20)
    elif int(sys.argv[2]) == 3:
        sims[sys.argv[1]](40)
else:
    print 'What?'
