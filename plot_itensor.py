import numpy as np
import matplotlib.pyplot as plt
import csv


# plot energy in file ./log_folder/fname_energies.csv
def plot_energies(model='ising', N=16, spin='Half', compare=False):
    fname = model + str(N)
    # readout energies
    log_folder = './logs_spin' + spin + '/'
    plot_folder = './plots_spin' + spin + '/'
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

    # readout comparison energies
    if compare:
        file_comp = './exact-diagonalisation/' + model + "_" + str(N) + '.csv'
        j_comp = []
        energies_comp = []
        with open(file_comp, 'r') as compfile:
            log = csv.reader(compfile, delimiter=',')
            for row in log:
                if row[0].startswith('#'):
                    continue  # ignore
                else:
                    j_comp.append(float(row[0]))
                    energies_comp.append(row[1:])
        energies_comp = np.array(energies_comp, dtype=float)
        compfile.close()

    # sort energies, in case GS converged to an excited state
    for cnt, row in enumerate(energies):
        energies[cnt] = np.sort(energies[cnt])

    # plot energies
    #for e in range(np.shape(energies)[1]):
    for e in range(3):
        plt.plot(j, energies[:, e], marker='x', lw=1, ms=4, label='MPS' + str(e))
    # plot comparison energies
    if compare:
        for e in range(3):
            plt.plot(j_comp, energies_comp[:, e], marker='o', lw=1, ms=4, label='ED' + str(e))
        plt.xlim([0.1, 2.1])
    # plt.grid(axis='both')
    plt.title('ITensor Energies | ' + fname + ' spin' + spin)
    plt.ylabel('energy')
    plt.xlabel('J')
    plt.legend(loc="upper right")
    plt.savefig(plot_folder + fname + '_energies.png', dpi=400)
    plt.show()
    plt.close()
    plt.cla()
    plt.clf()

    # plot difference MPS to ED
    if compare:
        for e in range(3):
            plt.plot(j_comp, np.abs(energies[20:, e]-energies_comp[:, e]), marker='x', lw=1, ms=4, label='|MPS-ED|' + str(e))
        plt.xlim([0.1, 2.1])
        # plt.grid(axis='both')
        plt.title('ITensor vs ED | ' + fname  + ' spin' + spin)
        plt.ylabel('energy')
        plt.xlabel('J')
        plt.legend(loc="upper right")
        plt.savefig(plot_folder + fname + '_mps_vs_ed.png', dpi=400)
        plt.show()
        plt.close()
        plt.cla()
        plt.clf()

    # plot energie differences
    for e in range(np.shape(energies)[1] - 1):
        diff = np.abs(energies[:, 0] - energies[:, e + 1])
        plt.plot(j, diff, marker='x', lw=1, ms=4, label='|E0-E' + str(e + 1) + '|')
    # plt.grid(axis='both')
    plt.title('ITensor Energy gap | ' + fname + ' spin' + spin)
    plt.ylabel('energy gap')
    plt.xlabel('J')
    plt.legend(loc="upper right")
    plt.savefig(plot_folder + fname + '_gap.png', dpi=400)
    plt.show()
    plt.close()
    plt.cla()
    plt.clf()


# ising16, ising32, heisenberg16, heisenberg32
#plot_energies(model='ising', N=4, spin='Half', compare=True)
#plot_energies(model='ising', N=8, spin='Half', compare=True)
#plot_energies(model='ising', N=16, spin='Half', compare=True)
#plot_energies(model='ising', N=32, spin='Half', compare=False)
plot_energies(model='ising', N=64, spin='Half', compare=False)
plot_energies(model='ising', N=128, spin='Half', compare=False)

#edlf.table_to_files()