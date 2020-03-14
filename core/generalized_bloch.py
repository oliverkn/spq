import numpy as np
import qutip as qu


def gell_mann_matrices(n):
    l_list = []

    for k in range(0, n):
        for j in range(0, k):
            X = np.zeros((n, n), dtype=complex)
            X[j, k] = 1
            X[k, j] = 1
            l_list.append(X)

            Y = np.zeros((n, n), dtype=complex)
            Y[j, k] = -1j
            Y[k, j] = 1j
            l_list.append(Y)

    for l in range(1, n):
        Z = np.zeros((n, n))
        Z[l, l] = -l
        for j in range(1, l + 1):
            Z[j - 1, j - 1] = 1
        Z = np.sqrt(2 / (l * (l + 1))) * Z
        l_list.append(Z)

    return l_list


for l in gell_mann_matrices(3):
    print(l)

def dm_from_bloch(r, n):
    rho = np.zeros((n,n))

    rho += 1/n * np.identity(n)
