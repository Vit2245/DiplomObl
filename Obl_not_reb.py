import sys
from datetime import datetime
from os import replace
from typing import Union

import matplotlib.pyplot as plt
import numpy as np

import scipy.integrate as integrate
from symengine import expand, lambdify as lambdify_se
from sympy import Symbol, pi, sin, cos, symbols, diff, Matrix, S, integrate, latex, Integral, init_printing, \
    Derivative, cse, lambdify, Poly, separatevars, collect, expand_power_base, Mul, Add
from numpy import linalg as la
from sympy.integrals.trigonometry import trigintegrate

sys.setrecursionlimit(10 ** 6)
Num = Union[int, float]

init_printing()


class StopWatch:
    def __init__(self):
        self.start_time = datetime.now()

    def lap(self):
        return datetime.now() - self.start_time
    
    def remark(self, comment: str):
        print(comment + ':', self.lap())


def print_latex(*args):
    print(*[latex(func) + '\n' for func in args])


def output_latex(path: str, *args):
    with open(path, 'w') as file:
        file.write(
            '\\documentclass{article}\n\\usepackage{breqn}\n\\usepackage[margin=0.1in]{geometry}\n\\begin{document}\n\\begin{dmath}\n')
        file.write(str(*[latex(func) + '\n ' for func in args]).replace('psi', '\\Psi^'))
        file.write('\\end{dmath}\n\\end{document}')


# def scipy_integrate(functional, definitions):
#     # like: scipy_integrate(Es, (x, 0, aa))
#     return integrate.quad(functional, definitions[1], definitions[2])


stop_watch = StopWatch()
remark = stop_watch.remark

n = 1
N = n ** 2

x = Symbol('x')
y = Symbol('y')
q = Symbol('q')
i = Symbol('i')
aa = Symbol('aa')
bb = Symbol('bb')
k_x = Symbol('kx')
k_y = Symbol('ky')
E_1 = Symbol('E1')
E_2 = Symbol('E2')
k = Symbol('k')
r = Symbol('r')
z = Symbol('z')
mu_12 = Symbol('mu12')
mu_21 = Symbol('mu21')
h = Symbol('h')
G_12 = Symbol('G12')
G_13 = Symbol('G13')
G_23 = Symbol('G23')
A = Symbol('A')
B = Symbol('B')


def create_functional(n):
    f = 6 * (1 / 4 - i ** 2 / h ** 2)

    X_1 = sin(2 * i * pi * x / aa)
    X_2 = sin((2 * i - 1) * pi * x / aa)
    X_3 = sin((2 * i - 1) * pi * x / aa)
    X_4 = cos((2 * i - 1) * pi * x / aa)
    X_5 = sin((2 * i - 1) * pi * x / aa)
    Y_1 = sin((2 * i - 1) * pi * y / bb)
    Y_2 = sin(2 * i * pi * y / bb)
    Y_3 = sin((2 * i - 1) * pi * y / bb)
    Y_4 = sin((2 * i - 1) * pi * y / bb)
    Y_5 = cos((2 * i - 1) * pi * y / bb)

    # print_latex(X_1, X_2, X_3, X_4, X_5, Y_1, Y_2, Y_3, Y_4, Y_5)

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
        for g in range(n):
            u[m][g] = (symbols(f'u{m + 1}{g + 1}'))
            v[m][g] = (symbols(f'v{m + 1}{g + 1}'))
            w[m][g] = (symbols(f'w{m + 1}{g + 1}'))
            psix[m][g] = (symbols(f'psix{m + 1}{g + 1}'))
            psiy[m][g] = (symbols(f'psiy{m + 1}{g + 1}'))

    SN = []
    for line in u + v + w + psix + psiy:
        SN += line

    for m in range(n):
        for g in range(n):
            U += u[m][g] * X_1.subs(i, m + 1) * Y_1.subs(i, g + 1)
            V += v[m][g] * X_2.subs(i, m + 1) * Y_2.subs(i, g + 1)
            W += w[m][g] * X_3.subs(i, m + 1) * Y_3.subs(i, g + 1)
            Psix += psix[m][g] * X_4.subs(i, m + 1) * Y_4.subs(i, g + 1)
            Psiy += psiy[m][g] * X_5.subs(i, m + 1) * Y_5.subs(i, g + 1)

    # print_latex(U, V, W, Psix, Psiy)

    Theta1 = -(diff(W, x)) / A - k_x * U

    Theta2 = -(diff(W, y)) / B - k_y * V

    ex = (diff(U, x)) / A + (diff(A, y)) * V / (A * B) - k_x * W + (S(1) / 2) * Theta1 ** 2

    ey = (diff(V, y)) / B + (diff(B, x)) * U / (A * B) - k_y * W + (S(1) / 2) * Theta2 ** 2

    gxy = (diff(V, x)) / A + (diff(U, y)) / B - (diff(A, y)) * U / (A * B) - (diff(B, x) * V / (
            A * B) + Theta1 * Theta2)

    gxz = k * (f.subs(i, z)) * (Psix - Theta1)
    gyz = k * (f.subs(i, z)) * (Psiy - Theta2)

    kappa1 = (diff(Psix, x)) / A + (diff(A, y)) * Psiy / (A * B)
    kappa2 = (diff(Psiy, y)) / B + (diff(B, x)) * Psix / (A * B)
    kappa12 = S(1) / 2 * (
            (diff(Psiy, x)) / A + (diff(Psix, y)) / B - ((diff(A, y)) * Psix + (diff(B, x)) * Psiy) / (
            A * B))

    print("#Mx")

    Mx = (S(1) / 12) * E_1 * h ** 3 * (mu_21 * kappa2 + kappa1) / (1 - mu_12 * mu_21)

    print("#My")

    My = (S(1) / 12) * E_2 * h ** 3 * (mu_12 * kappa1 + kappa2) / (1 - mu_12 * mu_21)

    print("#Mxy")
    Mxy = (S(1) / 6) * G_12 * h ** 3 * kappa12

    print("#Myx")
    Myx = (S(1) / 6) * G_12 * h ** 3 * kappa12

    print("#Nx")
    Nx = (E_1 * h / (1 - mu_12 * mu_21)) * (ex + mu_21 * ey)

    print("#Ny")
    Ny = (E_2 * h / (1 - mu_12 * mu_21)) * (ey + mu_12 * ex)

    Nxy = G_12 * h * gxy
    print("Nxy")

    Nyx = G_12 * h * gxy
    print("Nyx")

    Px = 0
    Py = 0
    Qx = G_13 * k * h * (Psix - Theta1)
    print("Qx")

    Qy = G_23 * k * h * (Psiy - Theta2)

    print("Qy")

    Epp1 = Nx * ex + Ny * ey
    Epp3 = Epp1 + S(1) / 2 * (Nxy + Nyx) * gxy
    Epp4 = Epp3 + Mx * kappa1 + My * kappa2
    Epp6 = Epp4 + (Mxy + Myx) * kappa12
    Epp7 = Epp6 + Qx * (Psix - Theta1)
    Epp8 = Epp7 + Qy * (Psiy - Theta2)
    EP = S(1) / 2 * Epp8
    AA = Px * U + Py * V + W * q

    Es = EP - AA
    return Es, SN, W


Es, SN, W = create_functional(n)

values = {
    aa: round(60 * 0.09, 2),
    bb: round(60 * 0.09, 2),
    k_x: 1 / 225 * 0.09,
    k_y: 1 / 225 * 0.09,
    E_1: 2.1 * 10 ** 5,
    E_2: 2.1 * 10 ** 5,
    k: 5 / 6,
    r: 225 * 0.09,
    z: -(1 / 2) * 0.09,
    mu_12: 0.3,
    mu_21: 0.3,
    h: 0.09,
    G_12: 0.33 * 10 ** 5,
    G_13: 0.33 * 10 ** 5,
    G_23: 0.33 * 10 ** 5,
    A: 1.,
    B: 1.
}

Es_lambda_values = lambdify(values.keys(), Es, 'sympy')
Es = Es_lambda_values(*values.values())

remark('replacing is done')

Es = expand(Es)

remark('expanding is done')

Es_diff = [Mul(*arg.args[:-2]) for arg in Es.args]
Es_int = [Mul(*arg.args[-2:]) for arg in Es.args]

remark('Derivatives and Integrals have been separated')

Es_int_x = [Mul(*[nested_arg for nested_arg in arg.args if nested_arg.has(x)]) for arg in Es_int]
Es_int_y = [Mul(*[nested_arg for nested_arg in arg.args if nested_arg.has(y)]) for arg in Es_int]

remark('functional created')

Int_lambdas = [lambdify(term, x) for term in Es_int_x]
Jacobi = []
for i in SN:
    Jacobi.append(diff(Es_integrated, i))

Hessian = []
for dpU in Jacobi:
    lineOfHessian = []
    for symb in SN:
        lineOfHessian.append(diff(dpU, symb))
    Hessian.append(lineOfHessian)

remark('Jacobi and Hessian created (without integrating)')

Jacobi = Matrix(Jacobi)
Hessian = Matrix(Hessian)

remark('matrices created')

Coef = np.zeros(len(SN), dtype=np.float)
XkPred = np.array(Coef)

MasRes = []
res2 = []
WC = []
WCC = []
wcWW = []
WC2 = []
cBufV = np.zeros((5 * N), dtype=float)
Buf = np.zeros((5 * N), dtype=float)
dict_coef = dict(zip(SN, list(Coef)))
dict_coef.update({q: 0.})

print(Jacobi)

lambda_hessian = lambdify(dict_coef.keys(), Hessian, 'sympy')
lambda_jacobi = lambdify(dict_coef.keys(), Jacobi, 'sympy')

remark('preparations is done')

# # Computing is beginning
# # ─────────────────────────────────────────────────────
aa1 = 0
x_center = (values[aa] + aa1) / 2
x_quarter = (values[aa] + aa1) / 4
y_center = values[bb] / 2
y_quarter = values[bb] / 4
epsilon = 1 * 10 ** (-5)
delta_q = 0.1
MAX = 33

Q_y = []

for qi in range(0, MAX + 1):
    qq = round(delta_q * qi, 2)  # Увеличиваем нагрузку
    dict_coef.update({q: qq})
    print('Увеличиваем нагрузку qq={: f}'.format(qq), " коэффициенты: ", end="")
    delta = 1
    kol_iter = 0
    print(dict_coef)
    while delta > epsilon and kol_iter <= 15:
        dict_coef.update(zip(SN, list(Coef)))

        dict_values = dict_coef.values()
        Hessian = lambda_hessian(*dict_values)
        Jacobi = lambda_jacobi(*dict_values)

        print(Jacobi)

        Rans = np.dot(np.array(la.inv(Hessian)), Jacobi).reshape(Coef.shape)
        tmp = Coef - Rans
        Coef = np.array(tmp)  # Находим решение методом Ньютона

        delta = np.sum(np.abs(Coef - XkPred)) / len(Coef)  # ??

        XkPred = np.array(Coef)

        kol_iter += 1

    print("kol_iter=", kol_iter, "delta=", delta)
    wc1 = W
    Xk_new = list(Coef)
    for wi in range(2 * N, 3 * N):
        wc1 = wc1.subs(SN[wi], Coef[wi])
    # масив значений функции W c подставленными коэф. с в завимости от q
    wc11 = wc1.subs(x, x_center)
    wc = wc11.subs(y, y_center)
    WC.append(wc)
    Q_y.append(qq)
    wc2 = W
    Xk_new = list(Coef)
    for wi in range(2 * N, 3 * N):
        wc2 = wc2.subs(SN[wi], Coef[wi])
    # масив значений функции W c подставленными коэф. с в завимости от q
    wc22 = wc2.subs(x, x_quarter)
    wc23 = wc22.subs(y, y_quarter)
    WC2.append(wc23)

remark('answer calculated')

fig = plt.figure(num=1, figsize=(8, 6))
plt.plot(WC, Q_y, color='r', linestyle='--', marker='o', markersize=3, label='W((a+a1)/2,b/2)')
# plt.plot(WC2, Q_y, color='b', linestyle='--', marker='o', markersize=3, label='W((a+a1)/4,b/4)')
plt.legend(loc='upper left')
grid1 = plt.grid(True)
plt.xlabel("W,м")
plt.ylabel("q,МПа")
plt.title('График прогиба W')
plt.show()

remark('full time')
