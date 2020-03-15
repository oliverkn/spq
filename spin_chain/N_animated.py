import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import qutip as q


# %%
def couple(A, B, i, N):
    '''
    :return: id_0 x ... x id_i-1 x C x id_i+2 x ... x id_N-1
    '''

    if i == 0:
        return q.tensor(A, B, q.identity([2] * (N - 2)))
    elif i == N - 2:
        return q.tensor(q.identity([2] * (N - 2)), A, B)
    elif i == N - 1:
        return q.tensor(B, q.identity([2] * (N - 2)), A)
    elif i >= N:
        raise Exception('invalid index')
    else:
        return q.tensor(q.identity([2] * i), A, B, q.identity([2] * (N - 2 - i)))


# %%

N = 10
dims = [2] * N

H = q.qzero(dimensions=dims)
H = H + N * q.identity(dims)

sx, sy, sz = q.sigmax(), q.sigmay(), q.sigmaz()

for i in range(0, N - 1):
    H = H - couple(sx, sx, i, N)
    H = H - couple(sy, sy, i, N)
    H = H - couple(sz, sz, i, N)

# print(H)
# print(H.eigenenergies())
# print(H.eigenstates())

# %%

N_e = 1

psi_0 = [q.basis(2, 0), ] * (N - N_e)
psi_0.extend([q.basis(2, 1), ] * N_e)

# psi_0 = [q.basis(2, 0), ] * (N)

# print(psi_0)

psi_0 = q.tensor(*psi_0)
# print(psi_0)


# %%
T = 10
T_steps = 1000
t_arr = np.linspace(0, T, T_steps)
res = q.sesolve(H, psi_0, t_arr)
# print(res.states)

# %%

psi = res.states[0]
X = []
Y = [0] * N
U, V = [], []
for i in range(0, N):
    X.append(i)
    U.append(q.expect(sx, psi.ptrace(i)))
    V.append(q.expect(sz, psi.ptrace(i)))

fig, ax = plt.subplots(1, 1)
Q = ax.quiver(X, Y, U, V, pivot='mid', color='r', scale=5)

ax.set_xlim(-1, N)
ax.set_ylim(-1, 1)

ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03])
slider = Slider(ax_slider, 't', 0, T_steps - 1, valinit=0, valstep=1)


def update_quiver(t):
    psi = res.states[int(t)]
    U, V = [], []
    for i in range(0, N):
        U.append(q.expect(sx, psi.ptrace(i)))
        V.append(q.expect(sz, psi.ptrace(i)))

    Q.set_UVC(U, V)


slider.on_changed(update_quiver)

fig.tight_layout()
plt.show()

update_quiver(0)
# %%
psi = res.states[0]

X = []
Y = [0] * N
U, V, C = [], [], []
for i in range(0, N):
    amp = psi[2 ** i, 0]
    X.append(i)
    U.append(0)
    V.append(np.abs(amp))
    C.append(np.angle(amp))

C[0] = -np.pi
C[N - 1] = np.pi

fig, ax = plt.subplots(1, 1)
Q = ax.quiver(X, Y, U, V, C, cmap='hsv', pivot='mid', scale=5)

ax.set_xlim(-1, N)
ax.set_ylim(-1, 1)

ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03])
slider = Slider(ax_slider, 't', 0, T_steps - 1, valinit=0, valstep=1)


def update_quiver(t):
    psi = res.states[int(t)]
    U, V, C = [], [], []
    for i in range(0, N):
        amp = psi[2 ** i, 0]
        # U.append(np.imag(amp))
        # V.append(np.real(amp))
        U.append(0)
        V.append(np.abs(amp))
        C.append(np.angle(amp))

    Q.set_UVC(U, V, C)


slider.on_changed(update_quiver)
fig.tight_layout()
plt.show()
