import numpy as np
import qutip as qu


def get_dm(z_arr, eig_vals, n):
    assert abs(np.sum(eig_vals) - 1) < 1e-6

    U_n = qu.identity(n)

    for j in range(2, n + 1):
        zj = np.array(z_arr[j - 1])
        A_j = compute_V(n, j, zj)
        U_n = A_j * U_n

    D = np.diag(eig_vals)
    D = qu.Qobj(D)

    rho = U_n * D * U_n.dag()

    assert abs(rho.tr() - 1) < 1e-6
    assert rho.isherm

    return rho


def get_qubit_dm(r, theta, phi):
    rho = qu.identity(2)
    rho = rho + r * np.sin(theta) * np.cos(phi) * qu.sigmax()
    rho = rho + r * np.sin(theta) * np.sin(phi) * qu.sigmay()
    rho = rho + r * np.cos(theta) * qu.sigmaz()
    rho = 0.5 * rho
    return rho


def apply_choi(rho_in, rho_choi, dim_a):
    return dim_a * qu.ptrace(qu.tensor(rho_in.trans(), qu.identity(2)) * rho_choi, 1)


def compute_V(n, j, z):
    z_u = qu.Qobj(z).unit()
    r = np.linalg.norm(z, ord=2)
    c = np.cos(r)
    s = np.sin(r)

    W = qu.identity(j - 1) - (1 - c) * z_u * z_u.dag()

    V = np.zeros((n, n))
    V[0:j - 1, 0:j - 1] = W
    V[0:j - 1, j - 1] = s * z_u.trans()
    V[j - 1, 0:j - 1] = -s * z_u.conj().trans()
    V[j - 1, j - 1] = c

    V[j:, j:] = np.identity(n - j)

    return qu.Qobj(V)


def get_hermitian(x, n):
    assert len(x) == n ** 2

    H = np.zeros((n, n), dtype=np.complex_)

    i = 0  # index for parameter

    for r in range(0, n):
        H[r, r] = x[i]
        i += 1
        for c in range(r + 1, n):
            H[r, c] = x[i] + 1.j * x[i + 1]
            H[c, r] = np.conjugate(H[r, c])
            i += 2

    return H


def get_unitary(p, n):
    H = get_hermitian(p, n)
    H = qu.Qobj(H)

    return (1.j * H).expm()


# %%
print(get_hermitian([1, 2, 3, 4], 2))

U = get_unitary([1, 2, 3, 4], 2)
print(U.dag() * U)

# %%
z_arr = [[], [1.], [1., 1.], [1., 100., 1.]]
eig_vals = [0.5, 0.1, 0.1, 0.3]

rho = get_dm(z_arr, eig_vals, 4)

print(rho)
print(rho.tr(), rho.isherm)
