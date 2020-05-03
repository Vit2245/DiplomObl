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
import pickle as pk
start_all=datetime.now()

sys.setrecursionlimit(10 ** 9)

start_time = datetime.now()
Num = Union[int, float]


n = 1
N = np.power(n, 2)
h = 0.09
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

print('Всего прошло времени')
print(datetime.now() - start_all)
x = sp.symbols('x')
y = sp.symbols('y')
q = sp.symbols('q')


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
print('Всего прошло времени')
print(datetime.now() - start_all)
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

print("#Mx")

Mx = sm.expand((1 / 12) * E1 * h ** 3 * (mu21 * varkappa2 + varkappa1) / (1 - mu12 * mu21))
# print(Mx)


print("#My")

My = sm.expand((1 / 12) * E2 * h ** 3 * (mu12 * varkappa1 + varkappa2) / (1 - mu12 * mu21))
# print(My)
print("#Mxy")
Mxy = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)

# print(Mxy)
print("#Myx")
Myx = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)
# print(Myx)
print("#Nx")
Nx = sm.expand((E1 * h / (1 - mu12 * mu21)) * (ex + mu21 * ey))
# print(Nx)

print("#Ny")
Ny = sm.expand((E2 * h / (1 - mu12 * mu21)) * (ey + mu12 * ex))
# print(Ny)


Nxy = sm.expand(G12 * h * gxy)
print("Nxy")
# print(Nxy)

Nyx = sm.expand(G12 * h * gxy)
print("Nyx")
# print(Nyx)

Px = 0
Py = 0
Qx = sm.expand(G13 * k * h * (Psix - Theta1))
print("Qx")
# print(Qx)

Qy = sm.expand(G23 * k * h * (Psiy - Theta2))

print("Qy")
# print(Qy)

# Epp=Nx * ex+Ny * ey+1 / 2 * (Nxy + Nyx) * gxy +Mx * varkappa1+My * varkappa2 +(Mxy + Myx) * varkappa12 +Qx * (Psix - Theta1)+Qy * (Psiy - Theta2)
# EPp=sm.expand(Epp)


Epp1 = Nx * ex + Ny * ey
print("Razbienie1")
# print(Epp1)
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

# print(sm.expand(AllEpp))
del (Epp1, Epp3, Epp4, Epp6, Epp7, Epp8)


start=datetime.now()
EPp = sm.expand(AllEpp)
print('Всего прошло времени')
print(datetime.now() - start)

start=datetime.now()
# EPp = sp.expand(EPp)
print('Всего прошло времени')
print(datetime.now() - start)

# print(EPp)
print("Время раскрытия скобок")
print(datetime.now() - start_time)
#
Epp=EPp.args
# EP = []
#
# for xc in Epp:
#     EP.append(str(sm.expand(xc)))

# for i in EP:
#     print(i)
print('Всего прошло времени')
print(datetime.now() - start_all)







# int_x = open('out_x.txt')
# with int_x as inp:
#     for i in inp.readlines():
#         key, val = i.strip().split(':')
#         val = val.strip(' ')
#         dict_x[key] = val
# int_x.close()
#
# int_y = open('out_y.txt')
# with int_y as inp:
#     for i in inp.readlines():
#         key, val = i.strip().split(':')
#         val = val.strip(' ')
#         dict_y[key] = val
# int_y.close()
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
print(dict_x)
print(dict_y)
def replace_by_dict(ep, dictionary: dict, a: Num, b: Num, s:sp.symbols):
    massznach=[]
    for arg in ep:
        # print(arg)
        mass=[]
        for nested_arg in arg.args:
            if nested_arg.has(s):
                mass.append(nested_arg)
        t= sp.Mul(*mass)
        if t in dictionary:
            # print("1")
            massznach.append(dictionary.get(t))
        else:
            # print("2")
            zam_x = sm.expand(t)
            if s==x:
                result = integrate.quad(lambda xx: (zam_x).subs(x, xx), a, b)[0]
                dictionary.update({t: result})
            else:
                result = integrate.quad(lambda yy: (zam_x).subs(y, yy), a, b)[0]
                dictionary.update({t: result})
            massznach.append(result)

    return massznach, dictionary



zxc,dict_x=replace_by_dict(Epp, dict_x, 0, 5.4, x)
print('Всего прошло времени1')
print(datetime.now() - start)


output = open('data_x.pkl', 'wb')
pk.dump(dict_x, output, 2)
output.close()


zxc2,dict_y=replace_by_dict(Epp, dict_y, 0, 5.4, y)
print('Всего прошло времени2')
print(datetime.now() - start)

output = open('data_y.pkl', 'wb')
pk.dump(dict_y, output, 2)
output.close()
#
# Es_diff = [sp.Mul(*[nested_arg for nested_arg in arg.args if not nested_arg.has(x) and not nested_arg.has(y)]) for arg in Epp]
# print('Всего прошло времени')
# print(datetime.now() - start)
# funcional=[]
# for i,el in enumerate(Es_diff):
#     funcional.append(sm.expand(1 / 2 * sm.expand(el)*zxc[i]*zxc2[i]))
# print('Всего прошло времени')
# print(datetime.now() - start)
#
# dict_all={}
#
# print("*****************")
# for i in funcional:
#     mas=i.args
#     if len(mas)!=0:
#         b=1
#         for ii in mas[1:]:
#             b*=ii
#
#         # print(b,mas[0])
#         # dict_all.update({b},mas[0])
#         if b in dict_all:
#             # print("нашел")
#             # print(zamena[1:])
#             result = dict_all.get(b)
#             result+=mas[0]
#             dict_all.update({b: result})
#         else:
#             dict_all.update({b:mas[0]})

# int_f = open('out_Fun.txt','w')
# with int_f as out:
#     for key, val in dict_all.items():
#         out.write('{}:{}\n'.format(key, val))
#
# int_f.close()

from functools import reduce
# itt=[]
# start_timeS = datetime.now()
# for key in dict_all:
#
#     itt.append(key*dict_all[key])

# print(itt)
#
#
#
# for arg in Epp:
#     print(arg)
#     t = 1
#     for nested_arg in arg.args:
#         if nested_arg.has(x):
#             t*=nested_arg
#     print(t)
#     arg.subs()
# def replace_by_dict(ep, variables: [str], dictionary: dict, a: Num, b: Num, s=sp.symbols) -> None:
#
#     for index, elem in enumerate(ep):
#         # start = datetime.now()
#         zamena = ''
#         for zz in variables:
#             if elem.has(s):
#                 i = elem.find(zz)
#                 if i < 0:
#                     continue
#                 if elem[:i]=='*':
#                     elem = elem[:i-1] + '' + elem[i-1 + len(zz) +1:]
#                 else:
#                     elem = elem[:i] + '' + elem[i  + len(zz) + 1:]
#                 zamena += '*' + zz
#                 ep[index] = elem
#             continue
#         if zamena[1:] in dictionary:
# #             result = dictionary.get(zamena[1:])
# #
# #             if elem[-1] == '*':
# #                 ep[index] = elem + '' + str(result)
# #                 # print("Нашел",datetime.now() - start)
# #             else:
# #                 ep[index] = elem + '*' + str(result)
# #                 # print("Нашел",datetime.now() - start)
# #         else:
# #             zam_x = sm.expand(zamena[1:])
# #             if s == 'x':
# #                 result = integrate.quad(lambda xx: (zam_x * A).subs(x, xx), a, b)[0]
# #
# #             else:
# #                 result = integrate.quad(lambda yy: (zam_x * B).subs(y, yy), a, b)[0]
# #             dictionary.update({zamena[1:]: result})
# #             if elem[-1] == '*':
# #                 ep[index] = elem +'' + str(result)
# #                 # print("Не Нашел",datetime.now() - start)
# #             else:
# #                 ep[index] = elem + '*' + str(result)
# #                 # print("Не Нашел",datetime.now() - start)
#

