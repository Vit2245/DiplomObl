import sys
from datetime import datetime
from typing import Union

import matplotlib.pyplot as plt
import numpy as np

import scipy.integrate as integrate
from sympy import Symbol, pi, sin, cos, symbols, diff, Matrix, lambdify, S, integrate, latex
from numpy import linalg as la

sys.setrecursionlimit(10 ** 6)
Num = Union[int, float]

def print_latex(*args):
    print(*[latex(func) + '\n' for func in args])

start_time = datetime.now()

n = 1
N = n ** 2

x = Symbol('x')
y = Symbol('y')
q = Symbol('q')
i = Symbol('i')
aa = Symbol('aa')
bb = Symbol('bb')
kx = Symbol('kx')
ky = Symbol('ky')
E1 = Symbol('E1')
E2 = Symbol('E2')
k = Symbol('k')
r = Symbol('r')
z = Symbol('z')
mu12 = Symbol('mu12')
mu21 = Symbol('mu21')
h = Symbol('h')
G12 = Symbol('G12')
G13 = Symbol('G13')
G23 = Symbol('G23')
A = Symbol('A')
B = Symbol('B')

def create_functional(n):
    f = 6 * (1 / 4 - i ** 2 / h ** 2)

    X1 = sin(2 * i * pi * x / aa)
    X2 = sin((2 * i - 1) * pi * x / aa)
    X3 = sin((2 * i - 1) * pi * x / aa)
    X4 = cos((2 * i - 1) * pi * x / aa)
    X5 = sin((2 * i - 1) * pi * x / aa)
    Y1 = sin((2 * i - 1) * pi * y / bb)
    Y2 = sin(2 * i * pi * y / bb)
    Y3 = sin((2 * i - 1) * pi * y / bb)
    Y4 = sin((2 * i - 1) * pi * y / bb)
    Y5 = cos((2 * i - 1) * pi * y / bb)

    # print_latex(X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5)

    U = 0
    V = 0
    W = 0
    Psix = 0
    Psiy = 0

    u = [[None] * n for _ in range(n)]
    v = [[None] * n for _ in range(n)]
    w = [[None] * n for _ in range(n)]
    psix = [[None] * n for _ in range(n)]
    psiy = [[None] * n for _ in range(n)]

    for m in range(n):
        for k in range(n):
            u[m][k] = (symbols(f'u{m + 1}{k + 1}'))
            v[m][k] = (symbols(f'v{m + 1}{k + 1}'))
            w[m][k] = (symbols(f'w{m + 1}{k + 1}'))
            psix[m][k] = (symbols(f'psix{m + 1}{k + 1}'))
            psiy[m][k] = (symbols(f'psiy{m + 1}{k + 1}'))

    SN = []
    for line in u + v + w + psix + psiy:
        SN += line

    for m in range(n):
        for k in range(n):
            U += u[m][k] * X1.subs(i, m + 1) * Y1.subs(i, k + 1)
            V += v[m][k] * X2.subs(i, m + 1) * Y2.subs(i, k + 1)
            W += w[m][k] * X3.subs(i, m + 1) * Y3.subs(i, k + 1)
            Psix += psix[m][k] * X4.subs(i, m + 1) * Y4.subs(i, k + 1)
            Psiy += psiy[m][k] * X5.subs(i, m + 1) * Y5.subs(i, k + 1)

    # print_latex(U, V, W, Psix, Psiy)

    Theta1 = -(diff(W, x)) / A - kx * U

    Theta2 = -(diff(W, y)) / B - ky * V

    ex = (diff(U, x)) / A + (diff(A, y)) * V / (A * B) - kx * W + (S(1) / 2) * Theta1 ** 2

    ey = (diff(V, y)) / B + (diff(B, x)) * U / (A * B) - ky * W + (S(1) / 2) * Theta2 ** 2

    gxy = (diff(V, x)) / A + (diff(U, y)) / B - (diff(A, y)) * U / (A * B) - (diff(B, x) * V / (
            A * B) + Theta1 * Theta2)

    gxz = k * (f.subs(i, z)) * (Psix - Theta1)
    gyz = k * (f.subs(i, z)) * (Psiy - Theta2)

    varkappa1 = (diff(Psix, x)) / A + (diff(A, y)) * Psiy / (A * B)
    varkappa2 = (diff(Psiy, y)) / B + (diff(B, x)) * Psix / (A * B)
    varkappa12 = S(1) / 2 * (
            (diff(Psiy, x)) / A + (diff(Psix, y)) / B - ((diff(A, y)) * Psix + (diff(B, x)) * Psiy) / (
            A * B))

    print("#Mx")

    Mx = (S(1) / 12) * E1 * h ** 3 * (mu21 * varkappa2 + varkappa1) / (1 - mu12 * mu21)

    print("#My")

    My = (S(1) / 12) * E2 * h ** 3 * (mu12 * varkappa1 + varkappa2) / (1 - mu12 * mu21)

    print("#Mxy")
    Mxy = (S(1) / 6) * G12 * h ** 3 * varkappa12

    print("#Myx")
    Myx = (S(1) / 6) * G12 * h ** 3 * varkappa12

    print("#Nx")
    Nx = (E1 * h / (1 - mu12 * mu21)) * (ex + mu21 * ey)

    print("#Ny")
    Ny = (E2 * h / (1 - mu12 * mu21)) * (ey + mu12 * ex)

    Nxy = G12 * h * gxy
    print("Nxy")

    Nyx = G12 * h * gxy
    print("Nyx")

    Px = 0
    Py = 0
    Qx = G13 * k * h * (Psix - Theta1)
    print("Qx")

    Qy = G23 * k * h * (Psiy - Theta2)

    print("Qy")

    Epp1 = Nx * ex + Ny * ey
    Epp3 = Epp1 + S(1) / 2 * (Nxy + Nyx) * gxy
    Epp4 = Epp3 + Mx * varkappa1 + My * varkappa2
    Epp6 = Epp4 + (Mxy + Myx) * varkappa12
    Epp7 = Epp6 + Qx * (Psix - Theta1)
    Epp8 = Epp7 + Qy * (Psiy - Theta2)
    EP = S(1) / 2 * integrate(integrate(Epp8 * A * B, (y, 0, aa)), (x, 0, bb))
    AA = integrate(integrate((Px * U + Py * V + W * q) * A * B, (y, 0, aa)), (x, 0, bb))

    Es = EP - AA
    print_latex(Es)
    return Es, SN, W


Es, SN, W = create_functional(n)

print("functional created", datetime.now() - start_time)

print("Es", Es)

Jacobi2 = np.array([0] * 5 * N)

Jacobi = []

for i in SN:
    Jacobi.append(Es.diff(i))

Deter = []

for dpU in Jacobi:
    lineOfHessian = []
    for symb in SN:
        lineOfHessian.append(dpU.diff(symb))
    Deter.append(lineOfHessian)

Jacobi = Matrix(Jacobi)
Deter = Matrix(Deter)

Coef = np.zeros(len(SN), dtype=np.float)

XkPred = np.array(Coef)

MasRes = []
res2 = []

WC = []
WCC = []
wcWW = []
WC2 = []

BufV = np.zeros((5 * N), dtype=float)
Buf = np.zeros((5 * N), dtype=float)

dict_coef = dict(zip(SN, list(Coef)))
dict_coef.update({q: 0.})

constants = x, y, aa, bb, kx, ky, E1, E2, k, r, z, mu12, mu21, h, G12, G13, G23, A, B

lambda_deter = lambdify([*dict_coef.keys(), *constants], Deter)
print("Deter", Deter)
lambda_jacobi = lambdify([*dict_coef.keys(), *constants], Jacobi)
print("Jacobi", Jacobi)

print("preparations is done", datetime.now() - start_time)

# Computing is beginning
# ─────────────────────────────────────────────────────
h = 0.09
r = 225 * h
A = 1
B = 1
aa1 = 0
aa = round(60 * h, 2)
bb = round(60 * h, 2)
h = 0.09
E1 = 2.1 * 10 ** 5
E2 = 2.1 * 10 ** 5
kx = 1 / r
ky = 1 / r
z = -(1 / 2) * h
r = 225 * h
k = 5 / 6
mu12 = 0.3
mu21 = 0.3
G12 = 0.33 * 10 ** 5
G13 = 0.33 * 10 ** 5
G23 = 0.33 * 10 ** 5
x_center = (aa + aa1) / 2
x_quarter = (aa + aa1) / 4
y_center = bb / 2
y_quarter = bb / 4
epsillon = 1 * 10 ** (-5)
delq = 0.1
MAX = 33
x = x_center
y = y_center

numbers = x, y, aa, bb, kx, ky, E1, E2, k, r, z, mu12, mu21, h, G12, G13, G23, A, B

Q_y = []

for qi in range(0, MAX + 1):
    qq = round(delq * qi, 2)  # Увеличиваем нагрузку
    dict_coef.update({q: qq})
    print('Увеличиваем нагрузку qq={: f}'.format(qq), " коэффициенты: ", end="")
    delta = 1
    kol_iter = 0
    print(dict_coef)
    while delta > epsillon or kol_iter <= 15:
        dict_coef.update(zip(SN, list(Coef)))

        dict_values = dict_coef.values()
        Deter1 = lambda_deter(*dict_values, *numbers)
        Jacobi1 = lambda_jacobi(*dict_values, *numbers)

        Rans = np.dot(np.array(la.inv(Deter1)), Jacobi1).reshape(Coef.shape)
        tmp = Coef - Rans
        Coef = np.array(tmp)  # Находим решение методом Ньютона

        delta = np.sum(np.abs(Coef - XkPred)) / len(Coef)  # ??

        XkPred = np.array(Coef)

        kol_iter += 1

    print("kol_iter=", kol_iter, "delta=", delta)
    wc1 = W
    Xk_new = list(Coef)
    # масив значений функции W c подставленными коэф. с в завимости от q
    wc11 = wc1
    wc = wc11
    WC.append(wc)
    # Q_y.append(qq)
    # wc2 = W
    # Xk_new = list(Coef)
    # for wi in range(2 * N, 3 * N):
    #     wc2 = wc2.subs(SN[wi], Coef[wi])
    # # масив значений функции W c подставленными коэф. с в завимости от q
    # wc22 = wc2.subs(x, x_quarter)
    # wc23 = wc22.subs(y, y_quarter)
    # WC2.append(wc23)

print("answer calculated", datetime.now() - start_time)

fig = plt.figure(num=1, figsize=(8, 6))
plt.plot(WC, Q_y, color='r', linestyle='--', marker='o', markersize=3, label='W((a+a1)/2,b/2)')
# plt.plot(WC2, Q_y, color='b', linestyle='--', marker='o', markersize=3, label='W((a+a1)/4,b/4)')
plt.legend(loc='upper left')
grid1 = plt.grid(True)
plt.xlabel("W,м")
plt.ylabel("q,МПа")
plt.title('График прогиба W')
plt.show()

print(datetime.now() - start_time)
