import numpy as np
import matplotlib.pyplot as plt
import csv


def table_to_files():
    j = []
    energies = []
    with open('table.csv', 'r') as ofile:
        log = csv.reader(ofile, delimiter=',')
        sites = 4
        for row_num, row in enumerate(log):
            if row[0].startswith('#'):
                continue  # ignore
            else:
                if row_num == 0:
                    continue
                if row_num == 1:
                    j.append(row[2:])
                    j = np.array(j, dtype=float)
                    #j = np.expand_dims(j, axis=0)
                    j = j.T
                    print('j shape', np.shape(j))
                    print(j, '\n')
                    continue

                energies.append(row[2:])
                print(row[2:])

                if row_num % 3 == 1:
                    print('\nsites', sites)
                    energies = np.array(energies, dtype=float)
                    energies = energies.T
                    print('energies', np.shape(energies))
                    print(energies)

                    # save in cleaned format
                    j_energies = np.append(j, energies, axis=1)

                    # save to file
                    np.savetxt('./0_exact-diagonalisation_archive-h/ising_' + str(sites) + '.csv', j_energies, delimiter=',')

                    # reset
                    sites += 1
                    energies = []

    energies = np.array(energies, dtype=float)
    ofile.close()
