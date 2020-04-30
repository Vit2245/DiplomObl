import sys
from datetime import datetime
from typing import Union

import cupy
import matplotlib.pyplot as plt
import numpy
import numpy as np
import scipy
import scipy.integrate as integrate
import sympy
from scipy.linalg import inv
from symengine import expand
from sympy import Symbol, pi, sin, cos, symbols, Derivative, S, latex, init_printing, \
    lambdify, Mul, Add, Basic

sys.setrecursionlimit(10 ** 6)
Num = Union[int, float]
Expr = Union[Basic, None]

init_printing()


class StopWatch:
    def __init__(self):
        self.start_time = datetime.now()

    def lap(self):
        return datetime.now() - self.start_time

    def remark(self, comment: str, *args, **kwargs):
        print(comment, ':', self.lap(), *args, **kwargs)


def print_latex(*args):
    print(*[latex(func) + '\n' for func in args])


def output_latex(path: str, *args):
    with open(path, 'w') as file:
        file.write(
            '\\documentclass{article}\n\\usepackage{breqn}\n\\usepackage[margin=0.1in]{geometry}\n\\begin{'
            'document}\n\\begin{dmath}\n')
        file.write(str(*[latex(func) + '\n ' for func in args]).replace('psi', '\\psi_^'))
        file.write('\\end{dmath}\n\\end{document}')


stop_watch = StopWatch()
remark = stop_watch.remark

n = 1
N = n ** 2

x = Symbol('x')
y = Symbol('y')
q = Symbol('q')
i = Symbol('i')
upper_limit_x = Symbol('aa')
upper_limit_y = Symbol('bb')
principal_curvature_x = Symbol('kx')
principal_curvature_y = Symbol('ky')
young_modulus_1 = Symbol('E1')
young_modulus_2 = Symbol('E2')
k = Symbol('k')
r = Symbol('r')
z = Symbol('z')
Poisson_coefficient_12 = Symbol('mu12')
Poisson_coefficient_21 = Symbol('mu21')
h = Symbol('h')
Shear_modulus_12 = Symbol('G12')
Shear_modulus_13 = Symbol('G13')
Shear_modulus_23 = Symbol('G23')
Lame_A = Symbol('A')
Lame_B = Symbol('B')


def create_functional(n):
    f = 6 * (1 / 4 - i ** 2 / h ** 2)

    approximate = {x: {}, y: {}}

    approximate[x][1] = sin(2 * i * pi * x / upper_limit_x)
    approximate[x][2] = sin((2 * i - 1) * pi * x / upper_limit_x)
    approximate[x][3] = sin((2 * i - 1) * pi * x / upper_limit_x)
    approximate[x][4] = cos((2 * i - 1) * pi * x / upper_limit_x)
    approximate[x][5] = sin((2 * i - 1) * pi * x / upper_limit_x)
    approximate[y][1] = sin((2 * i - 1) * pi * y / upper_limit_y)
    approximate[y][2] = sin(2 * i * pi * y / upper_limit_y)
    approximate[y][3] = sin((2 * i - 1) * pi * y / upper_limit_y)
    approximate[y][4] = sin((2 * i - 1) * pi * y / upper_limit_y)
    approximate[y][5] = cos((2 * i - 1) * pi * y / upper_limit_y)

    U = 0
    V = 0
    W = 0
    psi_x = 0
    psi_y = 0

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

    def fill(dim_1, dim_2, symbols_list, group, param) -> list:
        return symbols_list[dim_1][dim_2] \
               * approximate[x][group].subs(param, dim_1 + 1) \
               * approximate[y][group].subs(param, dim_2 + 1)

    for m in range(n):
        for g in range(n):
            U += fill(m, g, u, 1, i)
            V += fill(m, g, v, 2, i)
            W += fill(m, g, w, 3, i)
            psi_x += fill(m, g, psix, 4, i)
            psi_y += fill(m, g, psiy, 5, i)

    Theta1 = -(Derivative(W, x)) / Lame_A - principal_curvature_x * U
    Theta2 = -(Derivative(W, y)) / Lame_B - principal_curvature_y * V

    ex = (Derivative(U, x)) / Lame_A + (Derivative(Lame_A, y)) * V / (Lame_A * Lame_B) - principal_curvature_x * W + (
            S(1) / 2) * Theta1 ** 2

    ey = (Derivative(V, y)) / Lame_B + (Derivative(Lame_B, x)) * U / (Lame_A * Lame_B) - principal_curvature_y * W + (
            S(1) / 2) * Theta2 ** 2

    gxy = (Derivative(V, x)) / Lame_A + (Derivative(U, y)) / Lame_B - (Derivative(Lame_A, y)) * U / (
            Lame_A * Lame_B) - (
                  Derivative(Lame_B, x) * V / (
                  Lame_A * Lame_B) + Theta1 * Theta2)

    gxz = k * (f.subs(i, z)) * (psi_x - Theta1)
    gyz = k * (f.subs(i, z)) * (psi_y - Theta2)

    kappa1 = (Derivative(psi_x, x)) / Lame_A + (Derivative(Lame_A, y)) * psi_y / (Lame_A * Lame_B)
    kappa2 = (Derivative(psi_y, y)) / Lame_B + (Derivative(Lame_B, x)) * psi_x / (Lame_A * Lame_B)
    kappa12 = S(1) / 2 * (
            (Derivative(psi_y, x)) / Lame_A + (Derivative(psi_x, y)) / Lame_B - (
            (Derivative(Lame_A, y)) * psi_x + (Derivative(Lame_B, x)) * psi_y) / (
                    Lame_A * Lame_B))

    moment_x = (S(1) / 12) * young_modulus_1 * h ** 3 * (Poisson_coefficient_21 * kappa2 + kappa1) / (
            1 - Poisson_coefficient_12 * Poisson_coefficient_21)
    moment_y = (S(1) / 12) * young_modulus_2 * h ** 3 * (Poisson_coefficient_12 * kappa1 + kappa2) / (
            1 - Poisson_coefficient_12 * Poisson_coefficient_21)
    moment_xy = (S(1) / 6) * Shear_modulus_12 * h ** 3 * kappa12
    moment_yx = (S(1) / 6) * Shear_modulus_12 * h ** 3 * kappa12
    exertion_x = (young_modulus_1 * h / (1 - Poisson_coefficient_12 * Poisson_coefficient_21)) * (
            ex + Poisson_coefficient_21 * ey)
    exertion_y = (young_modulus_2 * h / (1 - Poisson_coefficient_12 * Poisson_coefficient_21)) * (
            ey + Poisson_coefficient_12 * ex)
    exertion_xy = Shear_modulus_12 * h * gxy
    exertion_yx = Shear_modulus_12 * h * gxy

    Px = 0
    Py = 0

    Qx = Shear_modulus_13 * k * h * (psi_x - Theta1)
    Qy = Shear_modulus_23 * k * h * (psi_y - Theta2)

    # TODO: refactor strings below using +=
    potential_energy1 = exertion_x * ex + exertion_y * ey
    potential_energy3 = potential_energy1 + S(1) / 2 * (exertion_xy + exertion_yx) * gxy
    potential_energy4 = potential_energy3 + moment_x * kappa1 + moment_y * kappa2
    potential_energy6 = potential_energy4 + (moment_xy + moment_yx) * kappa12
    potential_energy7 = potential_energy6 + Qx * (psi_x - Theta1)
    potential_energy8 = potential_energy7 + Qy * (psi_y - Theta2)

    EP = S(1) / 2 * potential_energy8
    AA = Px * U + Py * V + W * q

    Es = EP - AA

    return Es, SN, W


Es, SN, W = create_functional(n)


# function for for searching all ways that walk to the specified symbol
def recursive_cleaning_internal(expr, symbol):
    stack = []
    if expr == symbol:
        return expr
    for arg in expr.args:
        branch = recursive_cleaning_internal(arg, symbol)
        if branch:
            stack.append(branch)
    if stack:
        return expr.func(*stack)
    return


# Type hinting wrapper
def recursive_cleaning(expr: Basic, symbol: Symbol) -> Expr:
    return recursive_cleaning_internal(expr, symbol)


Ex = recursive_cleaning(Es, x)

# TODO: substitute the approximate functions to functional
# approximate[x][1] = sin(2 * i * pi * x / upper_limit_x)
# approximate[x][2] = sin((2 * i - 1) * pi * x / upper_limit_x)
# approximate[x][3] = sin((2 * i - 1) * pi * x / upper_limit_x)
# approximate[x][4] = cos((2 * i - 1) * pi * x / upper_limit_x)
# approximate[x][5] = sin((2 * i - 1) * pi * x / upper_limit_x)
# approximate[y][1] = sin((2 * i - 1) * pi * y / upper_limit_y)
# approximate[y][2] = sin(2 * i * pi * y / upper_limit_y)
# approximate[y][3] = sin((2 * i - 1) * pi * y / upper_limit_y)
# approximate[y][4] = sin((2 * i - 1) * pi * y / upper_limit_y)
# approximate[y][5] = cos((2 * i - 1) * pi * y / upper_limit_y)
# TODO: check the "7000-th" element of the functional. Does it have only four factors too?

values = {
    upper_limit_x: round(60 * 0.09, 2),
    upper_limit_y: round(60 * 0.09, 2),
    principal_curvature_x: 1 / (225 * 0.09),
    principal_curvature_y: 1 / (225 * 0.09),
    young_modulus_1: 2.1 * 10 ** 5,
    young_modulus_2: 2.1 * 10 ** 5,
    k: 5 / 6,
    r: 225 * 0.09,
    z: -(1 / 2) * 0.09,
    Poisson_coefficient_12: 0.3,
    Poisson_coefficient_21: 0.3,
    h: 0.09,
    Shear_modulus_12: 0.33 * 10 ** 5,
    Shear_modulus_13: 0.33 * 10 ** 5,
    Shear_modulus_23: 0.33 * 10 ** 5,
    Lame_A: 1.,
    Lame_B: 1.
}

W_lambda_values = lambdify(values.keys(), W, sympy)
W = W_lambda_values(*values.values())

Es_lambda_values = lambdify(values.keys(), Es, sympy)
Es = Es_lambda_values(*values.values())

remark('replacing is done')

Es = expand(Es)

remark('expanding is done')

Es_Derivative = [Mul(*[nested_arg for nested_arg in arg.args if not nested_arg.has(x) and not nested_arg.has(y)]) for
                 arg in Es.args]

remark('Derivatives have been separated')

Es_int_x = [Mul(*[nested_arg for nested_arg in arg.args if nested_arg.has(x)]) for arg in Es.args]
Es_int_y = [Mul(*[nested_arg for nested_arg in arg.args if nested_arg.has(y)]) for arg in Es.args]

remark('Integrals have been separated')

Int_lambdas_x = [lambdify(x, term * values[Lame_A], modules=[{'Mul': cupy.prod}, cupy]) for term in Es_int_x]
Int_lambdas_y = [lambdify(y, term * values[Lame_B], modules=[{'Mul': cupy.prod}, cupy]) for term in Es_int_y]

remark('lambdas for integrals have been created')

Integral_x = [scipy.integrate.quad(int_lambda, 0, values[upper_limit_x])[0] for int_lambda in Int_lambdas_x]
Integral_y = [scipy.integrate.quad(int_lambda, 0, values[upper_limit_y])[0] for int_lambda in Int_lambdas_y]

remark('integrals have been created')

Diff_sum = Add(*[Mul(term, Integral_x[i], Integral_y[i]) for i, term in enumerate(Es_Derivative)])

remark('General Derivativeerential function created')

Jacobi = []
for i in SN:
    Jacobi.append(Derivative(Diff_sum, i))

remark("Jacobi's list have been created")

Hessian = []
for dpU in Jacobi:
    lineOfHessian = []
    for symb in SN:
        lineOfHessian.append(Derivative(dpU, symb))
    Hessian.append(lineOfHessian)

remark('Hessian have been created')

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

lambda_hessian = lambdify(dict_coef.keys(), Hessian, numpy)
lambda_jacobi = lambdify(dict_coef.keys(), Jacobi, numpy)

remark('preparations were done')

# # Computing is beginning
# # ─────────────────────────────────────────────────────
aa1 = 0
x_center = (values[upper_limit_x] + aa1) / 2
x_quarter = (values[upper_limit_x] + aa1) / 4
y_center = values[upper_limit_y] / 2
y_quarter = values[upper_limit_y] / 4
epsilon = 1 * 10 ** (-5)
delta_q = 0.1
MAX = 40
MAX_ITERATION_COUNT = 150

Q_y = []

for qi in range(0, MAX + 1):
    qq = round(delta_q * qi, 2)  # Увеличиваем нагрузку
    dict_coef.update({q: qq})
    delta = 1
    kol_iter = 0

    remark(f'Increases the load q={qq}', end='\t')

    while delta > epsilon and kol_iter <= MAX_ITERATION_COUNT:
        dict_coef.update(zip(SN, list(Coef)))

        dict_values = dict_coef.values()
        Hessian = lambda_hessian(*dict_values)
        Jacobi = lambda_jacobi(*dict_values)

        Rans = np.dot(inv(np.array(Hessian)), Jacobi).reshape(Coef.shape)
        tmp = Coef - Rans
        Coef = np.array(tmp)  # Находим решение методом Ньютона

        delta = np.sum(np.abs(Coef - XkPred)) / len(Coef)  # ??

        XkPred = np.array(Coef)

        kol_iter += 1

    print("iterations count =", kol_iter, "delta =", delta)
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
plt.plot(WC2, Q_y, color='b', linestyle='--', marker='o', markersize=3, label='W((a+a1)/4,b/4)')
plt.legend(loc='upper left')
grid1 = plt.grid(True)
plt.xlabel("W,м")
plt.ylabel("q,МПа")
plt.title('График прогиба W')
plt.show()

remark('full time')
