from math import *
import sys
from datetime import datetime
from typing import Union

import matplotlib.pyplot as plt
import numpy as np

import scipy.integrate as integrate
import symengine as sm
import sympy as sp
from numpy import linalg as la

sys.setrecursionlimit(10 ** 6)
Num = Union[int, float]

start_time = datetime.now()

h = 0.09

n = 1
N = n ** 2

aa = round(60 * h, 2)
aa1 = 0
bb = round(60 * h, 2)

E1 = 2.1 * 10 ** 5
E2 = 2.1 * 10 ** 5
mu12 = 0.3
mu21 = 0.3

z = -(1 / 2) * h

r = 225 * h

k = 5 / 6
f = lambda i: 6 * (1 / 4 - i ** 2 / h ** 2)

A = 1
B = 1

kx = 1 / r
ky = 1 / r

G12 = 0.33 * 10 ** 5
G13 = 0.33 * 10 ** 5
G23 = 0.33 * 10 ** 5
x = sp.Symbol('x')
y = sp.Symbol('y')
q = sp.Symbol('q')
print(pi)

X1 = lambda i: sp.sin(2 * i * round(pi, 5) * x / aa)
X2 = lambda i: sp.sin((2 * i - 1) * round(pi, 5) * x / aa)
X3 = lambda i: sp.sin((2 * i - 1) * round(pi, 5) * x / aa)
X4 = lambda i: sp.cos((2 * i - 1) * round(pi, 5) * x / aa)
X5 = lambda i: sp.sin((2 * i - 1) * round(pi, 5) * x / aa)
Y1 = lambda i: sp.sin((2 * i - 1) * round(pi, 5) * y / bb)
Y2 = lambda i: sp.sin(2 * i * round(pi, 5) * y / bb)
Y3 = lambda i: sp.sin((2 * i - 1) * round(pi, 5) * y / bb)
Y4 = lambda i: sp.sin((2 * i - 1) * round(pi, 5) * y / bb)
Y5 = lambda i: sp.cos((2 * i - 1) * round(pi, 5) * y / bb)

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

for i in range(1, n + 1):
    for j in range(1, n + 1):
        u[i - 1][j - 1] = (sp.symbols(f'u{i}{j}'))
        v[i - 1][j - 1] = (sp.symbols(f'v{i}{j}'))
        w[i - 1][j - 1] = (sp.symbols(f'w{i}{j}'))
        psix[i - 1][j - 1] = (sp.symbols(f'psix{i}{j}'))
        psiy[i - 1][j - 1] = (sp.symbols(f'psiy{i}{j}'))

SN = []
for line in u + v + w + psix + psiy:
    SN += line

for i in range(1, n + 1):
    for j in range(1, n + 1):
        U += u[i - 1][j - 1] * X1(i) * Y1(j)
        V += v[i - 1][j - 1] * X2(i) * Y2(j)
        W += w[i - 1][j - 1] * X3(i) * Y3(j)
        Psix += psix[i - 1][j - 1] * X4(i) * Y4(j)
        Psiy += psiy[i - 1][j - 1] * X5(i) * Y5(j)

Theta1 = sm.expand(-(sp.diff(W, x)) / A - kx * U)

Theta2 = sm.expand(-(sp.diff(W, y)) / B - ky * V)

ex = sm.expand((sp.diff(U, x)) / A + (sp.diff(A, y)) * V / (A * B) - kx * W + (1 / 2) * Theta1 ** 2)

ey = sm.expand((sp.diff(V, y)) / B + (sp.diff(B, x)) * U / (A * B) - ky * W + (1 / 2) * Theta2 ** 2)

gxy = sm.expand((sp.diff(V, x)) / A + (sp.diff(U, y)) / B - (sp.diff(A, y)) * U / (A * B) - (sp.diff(B, x)) * V / (
        A * B) + Theta1 * Theta2)

gxz = sm.expand(k * (f(z)) * (Psix - Theta1))
gyz = sm.expand(k * (f(z)) * (Psiy - Theta2))

varkappa1 = sm.expand((sp.diff(Psix, x)) / A + (sp.diff(A, y)) * Psiy / (A * B))
varkappa2 = sm.expand((sp.diff(Psiy, y)) / B + (sp.diff(B, x)) * Psix / (A * B))
varkappa12 = sm.expand(1 / 2 * (
        (sp.diff(Psiy, x)) / A + (sp.diff(Psix, y)) / B - ((sp.diff(A, y)) * Psix + (sp.diff(B, x)) * Psiy) / (
        A * B)))

print("#Mx")

Mx = sm.expand((1 / 12) * E1 * h ** 3 * (mu21 * varkappa2 + varkappa1) / (1 - mu12 * mu21))

print("#My")

My = sm.expand((1 / 12) * E2 * h ** 3 * (mu12 * varkappa1 + varkappa2) / (1 - mu12 * mu21))

print("#Mxy")
Mxy = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)

print("#Myx")
Myx = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)

print("#Nx")
Nx = sm.expand((E1 * h / (1 - mu12 * mu21)) * (ex + mu21 * ey))

print("#Ny")
Ny = sm.expand((E2 * h / (1 - mu12 * mu21)) * (ey + mu12 * ex))

Nxy = sm.expand(G12 * h * gxy)
print("Nxy")

Nyx = sm.expand(G12 * h * gxy)
print("Nyx")

Px = 0
Py = 0
Qx = sm.expand(G13 * k * h * (Psix - Theta1))
print("Qx")

Qy = sm.expand(G23 * k * h * (Psiy - Theta2))

print("Qy")

Epp1 = Nx * ex + Ny * ey
print("Razbienie1")

del (Nx, ex, Ny, ey)
Epp3 = Epp1 + 1 / 2 * (Nxy + Nyx) * gxy
del (Nxy, Nyx, gxy)
print("Razbienie3")
Epp4 = Epp3 + Mx * varkappa1 + My * varkappa2
print("Razbienie4")
del (Mx, varkappa1, My, varkappa2)
Epp6 = Epp4 + (Mxy + Myx) * varkappa12
print("Razbienie6")
del (Mxy, Myx, varkappa12)
Epp7 = Epp6 + Qx * (Psix - Theta1)
print("Razbienie7")
del (Qx, Psix, Theta1)
Epp8 = Epp7 + Qy * (Psiy - Theta2)
print("Razbienie8")
del (Qy, Psiy, Theta2)

AllEpp = Epp8

del (Epp1, Epp3, Epp4, Epp6, Epp7, Epp8)

EPp = sm.expand(AllEpp)

print("Время раскрытия скобок")
print(datetime.now() - start_time)

Epp = str(EPp).split('+')

del (Epp[0])


def create_variables(n: Num, symbol: sp.Symbol, limit: Num) -> list:
    variables = []
    for i in range(1, n + 1):
        for st in range(4, 0, -1):
            variables.append(str(sp.sin(2 * i * round(pi, 5) * symbol / limit) ** st))
            variables.append(str(sp.cos(2 * i * round(pi, 5) * symbol / limit) ** st))
            variables.append(str(sp.sin((2 * i - 1) * round(pi, 5) * symbol / limit) ** st))
            variables.append(str(sp.cos((2 * i - 1) * round(pi, 5) * symbol / limit) ** st))
    return variables


variable_x = create_variables(n, x, aa)
variable_y = create_variables(n, y, bb)

EP = []

for xc in Epp:
    Epp = xc.split('-')
    if len(Epp) > 1:

        EP.append(Epp[0])
        del (Epp[0])
        for el in Epp:
            EP.append('-' + el)
    else:
        EP.append(Epp[0])

print(len(EP))
for i in EP:
    print(i)

my_dict = {}
dict_x = {}
dict_y = {}

int_x = open('out_x.txt')
with int_x as inp:
    for i in inp.readlines():
        key, val = i.strip().split(':')
        val = val.strip(' ')
        dict_x[key] = val
int_x.close()

int_y = open('out_y.txt')
with int_y as inp:
    for i in inp.readlines():
        key, val = i.strip().split(':')
        dict_y[key] = val.strip(' ')
int_y.close()


def replace_by_dict(ep: [str], variables: [str], dictionary: dict, a: Num, b: Num, x: sp.Symbol, lame: Num,
                    operator: str = '') -> None:
    for index, elem in enumerate(ep):
        zamena = ''
        for zz in variables:
            i = elem.find(zz)
            if i >= 0:
                elem = elem[:i] + '' + elem[i + len(zz) + 1:]
                zamena += '*' + zz
        zamena = zamena[1:]
        if zamena in dictionary:
            result = dictionary.get(zamena)
        else:
            zam_x = sm.expand(zamena)
            result = integrate.quad(lambda xx: (zam_x * lame).subs(x, xx), a, b)[0]
            dictionary.update({zamena: result})

        ep[index] = elem + operator + str(result)


replace_by_dict(EP, variable_x, dict_x, 0, 5.4, x, A)

print("Время раскрытия скобок")
print(datetime.now() - start_time)

with open('out_x.txt', 'w') as out:
    for key, val in dict_x.items():
        out.write('{}:{}\n'.format(key, val))

replace_by_dict(EP, variable_y, dict_y, 0, 5.4, y, B, '*')

print("Время раскрытия скобок")
print(datetime.now() - start_time)

with open('out_x.txt', 'w') as out:
    for key, val in dict_x.items():
        out.write('{}:{}\n'.format(key, val))
with open('out_y.txt', 'w') as out:
    for key, val in dict_y.items():
        out.write('{}:{}\n'.format(key, val))

int_x.close()
int_y.close()
number = 0
print("Resul")
for el in EP:
    print(sm.expand(el))

print("Время раскрытия скобок")
print(datetime.now() - start_time)
print("Resul1")

with open('агт.txt', 'w') as out:
    for el in EP:
        out.write(str(sm.expand(el)) + "\n")

allin = 0
for el in EP:
    allin += 1 / 2 * sm.expand(el)

print("Allin")
allin = sm.expand(allin)

AA = sp.integrate(sp.integrate((Px * U + Py * V + W * q) * A * B, (y, 0, 5.4)), (x, 0, 5.4))

AA = sm.expand(AA)
Es = allin - AA

Es = sm.expand(Es)
print("Es")
print(Es)

Jacobi2 = np.array([0] * 5 * N)
print(Jacobi2)

k = 0

Jacobi = []

for i in SN:
    Jacobi.append(sm.expand(Es.diff(i)))

Deter = []

for dpU in Jacobi:
    lineOfHessian = []
    for symb in SN:
        lineOfHessian.append(sm.expand(dpU.diff(symb)))
    Deter.append(lineOfHessian)

Jacobi = sp.Matrix(Jacobi)
Deter = sp.Matrix(Deter)

print('Начальный нулевой вектор ... ', end='')

epsillon = 1 * 10 ** (-5)

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

delq = 0.1
MAX = 33
Q_y = []

dict_coef = dict(zip(SN, list(Coef)))
dict_coef.update({q: 0.})

lambda_deter = sp.lambdify(dict_coef.keys(), Deter)
lambda_jacobi = sp.lambdify(dict_coef.keys(), Jacobi)

for qi in range(0, MAX + 1):
    qq = round(delq * qi, 2)  # Увеличиваем нагрузку
    dict_coef.update({q: qq})
    print('Увеличиваем нагрузку qq={: f}'.format(qq), " коэффициенты: ", end="")
    delta = 1
    kol_iter = 0
    print(dict_coef)
    while delta > epsillon:

        dict_coef.update(zip(SN, list(Coef)))

        dict_values = dict_coef.values()
        Deter1 = lambda_deter(*dict_values)
        Jacobi1 = lambda_jacobi(*dict_values)

        Rans = np.dot(np.array(la.inv(Deter1)), Jacobi1).reshape(Coef.shape)
        tmp = Coef - Rans
        Coef = np.array(tmp)  # Находим решение методом Ньютона

        delta = np.sum(np.abs(Coef - XkPred)) / len(Coef)  # ??

        XkPred = np.array(Coef)

        kol_iter += 1
        if kol_iter > 15:
            delta = 0

    print("kol_iter=", kol_iter, "delta=", delta)
    wc1 = W
    Xk_new = list(Coef)
    for wi in range(2 * N, 3 * N):
        wc1 = wc1.subs(SN[wi], Coef[wi])
    # масив значений функции W c подставленными коэф. с в завимости от q
    wc11 = wc1.subs(x, (aa + aa1) / 2)
    wc = wc11.subs(y, bb / 2)
    WC.append(wc)
    Q_y.append(qq)
    wc2 = W
    Xk_new = list(Coef)
    for wi in range(2 * N, 3 * N):
        wc2 = wc2.subs(SN[wi], Coef[wi])
    # масив значений функции W c подставленными коэф. с в завимости от q
    wc22 = wc2.subs(x, (aa + aa1) / 4)
    wc23 = wc22.subs(y, bb / 4)
    WC2.append(wc23)

print(datetime.now() - start_time)
fig = plt.figure(num=1, figsize=(8, 6))
plt.plot(WC, Q_y, color='r', linestyle='--', marker='o', markersize=3, label='W((a+a1)/2,b/2)')
plt.plot(WC2, Q_y, color='b', linestyle='--', marker='o', markersize=3, label='W((a+a1)/4,b/4)')
plt.legend(loc='upper left')
grid1 = plt.grid(True)
plt.xlabel("W,м")
plt.ylabel("q,МПа")
plt.title('График прогиба W')
plt.show()

print(datetime.now() - start_time)
