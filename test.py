# import time
from math import *
import sys
from datetime import datetime
# from threading import Thread
from typing import Union
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import symengine as sm
import sympy as sp

import numba as nb
from numba import njit, jit
from numpy import linalg as la
# from sympy.parsing.sympy_parser import parse_expr
from numba import jit, numba
from time import gmtime, strftime
import pickle as pk
class StopWatch:
    def __init__(self):
        self.start_time = datetime.now()

    def lap(self):
        return datetime.now() - self.start_time

    def remark(self, comment: str):
        print(comment + ':', self.lap())
start_all=datetime.now()
stop_watch = StopWatch()
remark = stop_watch.remark
sys.setrecursionlimit(10 ** 9)


Num = Union[int, float]


n = 2
N = np.power(n, 2)
##### Пологая
# h = 0.09
# aa = round(60 * h, 2)
# aa1 = 0
# bb = round(60 * h, 2)
# E1 = 2.1 * 10 ** 5
# E2 = 2.1 * 10 ** 5
# mu12 = 0.3
# mu21 = 0.3
# z = -(1 / 2) * h
# r = 225 * h
# k = 5 / 6
# f = lambda i: 6 * (1 / 4 - i ** 2 / h ** 2)
# A = 1
# B = 1
#
# kx = 1 / r
# ky = 1 / r
##### Пологая
x = sp.symbols('x')
y = sp.symbols('y')
q = sp.symbols('q')





#####Тор
h = 0.01
aa = round(1/2*round(pi, 5), 5)
aa1 = 0
bb = round(1/2*round(pi, 5), 5)
E1 = 2.1 * 10 ** 5
E2 = 2.1 * 10 ** 5
mu12 = 0.3
mu21 = 0.3
d=2
alpha=0
z = -(1 / 2) * h
r = 13
k = 5 / 6
f = lambda i: 6 * (1 / 4 - i ** 2 / h ** 2)
A = r
B = d+r*(sp.sin(x-alpha)+sp.sin(alpha))

kx = 1 / r
ky =sp.sin(x-alpha)/(d+r*(sp.sin(x-alpha)+sp.sin(alpha)))
#####Тор

G12 = 0.33 * 10 ** 5
G13 = 0.33 * 10 ** 5
G23 = 0.33 * 10 ** 5

remark('Начальные условия')


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

u = []
for l in range(0, n):
    u.append([0] * n)
v = []
for l in range(0, n):
    v.append([0] * n)
w = []
for l in range(0, n):
    w.append([0] * n)
psix = []
for l in range(0, n):
    psix.append([0] * n)
psiy = []
for l in range(0, n):
    psiy.append([0] * n)

for i in range(1, int(np.sqrt(N) + 1)):
    for j in range(1, int(np.sqrt(N) + 1)):
        u[i - 1][j - 1] = (sp.symbols('u' + str(i) + '' + str(j)))
        v[i - 1][j - 1] = (sp.symbols('v' + str(i) + '' + str(j)))
        w[i - 1][j - 1] = (sp.symbols('w' + str(i) + '' + str(j)))
        psix[i - 1][j - 1] = (sp.symbols('psix' + str(i) + '' + str(j)))
        psiy[i - 1][j - 1] = (sp.symbols('psiy' + str(i) + '' + str(j)))

SN = []

for i in range(1, int(np.sqrt(N) + 1)):
    for j in range(1, int(np.sqrt(N) + 1)):
        SN.append(sp.symbols('u' + str(i) + '' + str(j)))
for i in range(1, int(np.sqrt(N) + 1)):
    for j in range(1, int(np.sqrt(N) + 1)):
        SN.append(sp.symbols('v' + str(i) + '' + str(j)))
for i in range(1, int(np.sqrt(N) + 1)):
    for j in range(1, int(np.sqrt(N) + 1)):
        SN.append(sp.symbols('w' + str(i) + '' + str(j)))
for i in range(1, int(np.sqrt(N) + 1)):
    for j in range(1, int(np.sqrt(N) + 1)):
        SN.append(sp.symbols('psix' + str(i) + '' + str(j)))
for i in range(1, int(np.sqrt(N) + 1)):
    for j in range(1, int(np.sqrt(N) + 1)):
        SN.append(sp.symbols('psiy' + str(i) + '' + str(j)))

for i in range(1, int(np.sqrt(N) + 1)):
    for j in range(1, int(np.sqrt(N) + 1)):
        U = U + u[i - 1][j - 1] * X1(i) * Y1(j)
        V = V + v[i - 1][j - 1] * X2(i) * Y2(j)
        W = W + w[i - 1][j - 1] * X3(i) * Y3(j)
        Psix = Psix + psix[i - 1][j - 1] * X4(i) * Y4(j)
        Psiy = Psiy + psiy[i - 1][j - 1] * X5(i) * Y5(j)
remark('UVWPSIXPSIY')
# print(U)
# print(V)
# print(W)
# print(Psix)
# print(Psiy)

Theta1 = sm.expand(-(sp.diff(W, x)) / A - kx * U)
# print("Theta1")
# print(Theta1)


Theta2 = sm.expand(-(sp.diff(W, y)) / B - ky * V)
# print("Theta2")
# print(Theta2)
ex = sm.expand((sp.diff(U, x)) / A + (sp.diff(A, y)) * V / (A * B) - kx * W + (1 / 2) * Theta1 ** 2)
# print("#ex")
# print(ex)
ey = sm.expand((sp.diff(V, y)) / B + (sp.diff(B, x)) * U / (A * B) - ky * W + (1 / 2) * Theta2 ** 2)
# print("#ey")
# print(ey)
gxy = sm.expand((sp.diff(V, x)) / A + (sp.diff(U, y)) / B - (sp.diff(A, y)) * U / (A * B) - (sp.diff(B, x)) * V / (
        A * B) + Theta1 * Theta2)
# print("#gxy")
# print(gxy)
gxz = sm.expand(k * (f(z)) * (Psix - Theta1))
gyz = sm.expand(k * (f(z)) * (Psiy - Theta2))
# print("#gxz")
# print(gxz)
# print("#gyz")
# print(gyz)
varkappa1 = sm.expand((sp.diff(Psix, x)) / A + (sp.diff(A, y)) * Psiy / (A * B))
varkappa2 = sm.expand((sp.diff(Psiy, y)) / B + (sp.diff(B, x)) * Psix / (A * B))
varkappa12 = sm.expand(1 / 2 * (
        (sp.diff(Psiy, x)) / A + (sp.diff(Psix, y)) / B - ((sp.diff(A, y)) * Psix + (sp.diff(B, x)) * Psiy) / (
        A * B)))
# print("#varkappa1")

# print(varkappa1)
# print("#varkappa2")
# print(varkappa2)
# print("#varkappa12")
# print(varkappa12)

# print("#Mx")

Mx = sm.expand((1 / 12) * E1 * h ** 3 * (mu21 * varkappa2 + varkappa1) / (1 - mu12 * mu21))
# print(Mx)


# print("#My")

My = sm.expand((1 / 12) * E2 * h ** 3 * (mu12 * varkappa1 + varkappa2) / (1 - mu12 * mu21))
# print(My)
# print("#Mxy")
Mxy = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)

# print(Mxy)
# print("#Myx")
Myx = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)
# print(Myx)
# print("#Nx")
Nx = sm.expand((E1 * h / (1 - mu12 * mu21)) * (ex + mu21 * ey))
# print(Nx)

# print("#Ny")
Ny = sm.expand((E2 * h / (1 - mu12 * mu21)) * (ey + mu12 * ex))
# print(Ny)


Nxy = sm.expand(G12 * h * gxy)
# print("Nxy")
# print(Nxy)

Nyx = sm.expand(G12 * h * gxy)
# print("Nyx")
# print(Nyx)

Px = 0
Py = 0
Qx = sm.expand(G13 * k * h * (Psix - Theta1))
# print("Qx")
# print(Qx)
Qy = sm.expand(G23 * k * h * (Psiy - Theta2))
# print("Qy")
# print(Qy)

Epp1 = Nx * ex + Ny * ey
# print("Razbienie1")
# print(Epp1)
del (Nx, ex, Ny, ey)
Epp3 = Epp1 + 1 / 2 * (Nxy + Nyx) * gxy
del (Nxy, Nyx, gxy)
# print("Razbienie3")
Epp4 = Epp3 + Mx * varkappa1 + My * varkappa2
# print("Razbienie4")
del (Mx, varkappa1, My, varkappa2)
Epp6 = Epp4 + (Mxy + Myx) * varkappa12
# print("Razbienie6")
del (Mxy, Myx, varkappa12)
Epp7 = Epp6 + Qx * (Psix - Theta1)
# print("Razbienie7")
del (Qx, Psix, Theta1)
Epp8 = Epp7 + Qy * (Psiy - Theta2)
# print("Razbienie8")
del (Qy, Psiy, Theta2)

AllEpp = Epp8*A*B

# print(sm.expand(AllEpp))
del (Epp1, Epp3, Epp4, Epp6, Epp7, Epp8)

start=datetime.now()
EPp = sm.expand(AllEpp)

Epp=EPp.args
remark('Разбиение на слагемые')


Num = Union[int, float]

input = open('data_x.pkl', 'rb')
try:
    dict_x = pk.load(input)
except:
    dict_x={}
input.close()

input = open('data_y.pkl', 'rb')
try:
    dict_y = pk.load(input)
except:
    dict_y={}
input.close()
# print(dict_x)
# print(dict_y)

def replace_by_dict(ep, dictionary1: dict, dictionary2: dict, a: Num, b: Num,):
    massznach_x=[]
    massznach_y = []
    massznach=[]
    for arg in ep:
        # print(arg)
        mass_x = []
        mass_y = []
        mass=[]
        for nested_arg in arg.args:
            if not nested_arg.has(x) and not nested_arg.has(y):
                mass.append(nested_arg)
            if nested_arg.has(sp.sin(x)):
                mass_x.append(nested_arg)
            if nested_arg.has(y):
                mass_y.append(nested_arg)
        tx= sp.Mul(*mass_x)
        ty= sp.Mul(*mass_y)
        t=sp.Mul(*mass)
        if tx in dictionary1:
            massznach_x.append(dictionary1.get(tx))
        else:
            zam_x = sm.expand(tx)
            result = integrate.quad(lambda xx: (zam_x).subs(x, xx), a, b)[0]
            dictionary1.update({tx: result})
            massznach_x.append(dictionary1.get(tx))
        if ty in dictionary2:
            massznach_y.append(dictionary2.get(ty))
        else:
            zam_y = sm.expand(ty)
            result = integrate.quad(lambda yy: (zam_y).subs(y, yy), a, b)[0]
            dictionary2.update({ty: result})
            massznach_y.append(result)
        massznach.append(t)
    return massznach_x, dictionary1,massznach_y, dictionary2,massznach