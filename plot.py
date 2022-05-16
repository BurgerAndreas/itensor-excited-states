import numpy as np
import matplotlib.pyplot as plt
import csv

log_folder = './logs/'
plot_folder = './plots/'


# plot energy in file ./log_folder/fname_energies.csv
def plot_energies(fname='dmrg'):
    # readout energies
    log_file = log_folder + fname + '_energies.csv'
    j = []
    energies = []
    with open(log_file, 'r') as ofile:
        log = csv.reader(ofile, delimiter=',')
        for row in log:
            if row[0].startswith('#'):
                continue  # ignore
            else:
                j.append(float(row[0]))
                energies.append(row[1:])
    energies = np.array(energies, dtype=float)
    ofile.close()

    # plot energies
    for e in range(np.shape(energies)[1]):
        plt.plot(j, energies[:, e], marker='x', lw=1, ms=4, label='E' + str(e))
    # plt.grid(axis='both')
    plt.title('Energies')
    plt.ylabel('energy')
    plt.xlabel('J')
    plt.legend(loc="upper right")
    plt.savefig(plot_folder + fname + '_energies.png', dpi=400)
    plt.show()
    plt.close()
    plt.cla()
    plt.clf()

    # plot energie differences
    for e in range(np.shape(energies)[1] - 1):
        diff = np.abs(energies[:, 0] - energies[:, e + 1])
        plt.plot(j, diff, marker='x', lw=1, ms=4, label='|E0-E' + str(e+1) + '|')
    # plt.grid(axis='both')
    plt.title('Energy gap')
    plt.ylabel('energy gap')
    plt.xlabel('J')
    plt.legend(loc="upper right")
    plt.savefig(plot_folder + fname + '_gap.png', dpi=400)
    plt.show()
    plt.close()
    plt.cla()
    plt.clf()


plot_energies(fname='ising')
