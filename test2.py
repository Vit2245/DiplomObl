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
x = sp.symbols('x')
y = sp.symbols('y')
q = sp.symbols('q')
#### Пологая
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
#### Пологая






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
print(len(Epp))
remark('Разбиение на слагемые')

int_ff = open('aa1.txt','w')
with int_ff as out:
    for key in Epp:
        out.write(str(key)+'\n')
int_ff.close()
Num = Union[int, float]

try:
    input = open('data_x.pkl', 'rb')
    dict_x = pk.load(input)
    input.close()
except:
    dict_x={}



try:
    input = open('data_y.pkl', 'rb')
    dict_y = pk.load(input)
    input.close()
except:
    dict_y={}

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
            need=str(nested_arg)
            # print(need)
            # print(need.find('x'))
            # print(need.find('y'))
            # print(need.find('x'))
            if need.find('x')==-1 and need.find('y')==-1 or need.find('psi')!=-1:
                mass.append(nested_arg)
                continue
            if need.find('x')!=-1:
                mass_x.append(nested_arg)
                continue
            if need.find('y')!=-1:
                mass_y.append(nested_arg)
                continue
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


zxc,dict_x,zxc2,dict_y,mrg=replace_by_dict(Epp, dict_x,dict_y, 0, aa)
remark('Начальные условия')
output = open('data_x.pkl', 'wb')
pk.dump(dict_x, output, 2)
output.close()
output = open('data_y.pkl', 'wb')
pk.dump(dict_y, output, 2)
output.close()
remark('Начальные условия')
funcional=[]
for i in range(len(zxc)):
    funcional.append(sm.expand(1 / 2 * sm.expand(mrg[i])*zxc[i]*zxc2[i]))
remark('Начальные условия')
dict_all={}
print("*****************")
for i in funcional:
    mas=i.args
    if len(mas)!=0:
        b=1
        for ii in mas[1:]:
            b*=ii
        if b in dict_all:
            result = dict_all.get(b)
            result+=mas[0]
            dict_all.update({b: result})
        else:
            dict_all.update({b:mas[0]})

int_f = open('out_Fun.txt','w')
with int_f as out:
    for key, val in dict_all.items():
        out.write('{}:{}\n'.format(key, val))
int_f.close()

from functools import reduce
itt=[]
start_timeS = datetime.now()
for key in dict_all:
    itt.append(key*dict_all[key])


remark('Начальные условия')

nnn=len(itt)//16

product=0


start_timeS = datetime.now()
product += reduce((lambda x, y: x + y), itt[:nnn])
product += reduce((lambda x, y: x + y), itt[nnn:2*nnn])
product += reduce((lambda x, y: x + y), itt[2*nnn:3*nnn])
product += reduce((lambda x, y: x + y), itt[3*nnn:4*nnn])
product += reduce((lambda x, y: x + y), itt[4*nnn:5*nnn])
product += reduce((lambda x, y: x + y), itt[5*nnn:6*nnn])
product += reduce((lambda x, y: x + y), itt[6*nnn:7*nnn])
product += reduce((lambda x, y: x + y), itt[7*nnn:8*nnn])
product += reduce((lambda x, y: x + y), itt[8*nnn:9*nnn])
product += reduce((lambda x, y: x + y), itt[9*nnn:10*nnn])
product += reduce((lambda x, y: x + y), itt[10*nnn:11*nnn])
product += reduce((lambda x, y: x + y), itt[11*nnn:12*nnn])
product += reduce((lambda x, y: x + y), itt[12*nnn:13*nnn])
product += reduce((lambda x, y: x + y), itt[13*nnn:14*nnn])
product += reduce((lambda x, y: x + y), itt[14*nnn:15*nnn])
product += reduce((lambda x, y: x + y), itt[15*nnn:16*nnn])
product += reduce((lambda x, y: x + y), itt[16*nnn:])

product=sm.expand(product)

print("Время сумма2")
print(datetime.now() - start_timeS)

print("Время Allin")
start_time = datetime.now()
AA = sp.integrate(sp.integrate((Px * U + Py * V + W * q) * A * B, (y, 0, bb)), (x, 0, aa))
# print(AA)

AA = sm.expand(AA)
Es = product-AA
Es = sm.expand(Es)
print("Es")
print(datetime.now() - start_time)
# print(Es)
remark('Начальные условия')
# print(Coef)


start_time1 = datetime.now()
Jacobi = []
for i in SN:
    Jacobi.append(sm.expand(sm.diff(Es, i)))

print("Время первая производная")
print(datetime.now() - start_time1)
remark('Начальные условия')
Deter = []
start_time2 = datetime.now()
for dpU in Jacobi:
    lineOfHessian = []
    for symb in SN:
        lineOfHessian.append(sm.diff(dpU, symb))
    Deter.append(lineOfHessian)
del(itt)
del(dict_y)
del(dict_x)
del(product)
del(AA)
del(Es)
print("Время вторая производная")
print(datetime.now() - start_time2)




epsillon = 1 * 10 ** (-5)

print('Начальный нулевой вектор ... ')
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
MAX = 40
Q_y = []

dict_coef = dict(zip(SN, list(Coef)))
dict_coef.update({q: 0.})

start_time = datetime.now()
lambda_deter = sm.lambdify([*dict_coef.keys()], Deter)

print("Время матрицы Dict Deter")
print(datetime.now() - start_time)

start_time = datetime.now()
lambda_jacobi= sm.lambdify([*dict_coef.keys()], Jacobi)

print("Время матрицы dict Jacobi")
print(datetime.now() - start_time)

remark('Начальные условия')
start_time2 = datetime.now()
for qi in range(0, MAX + 1):
    qq = round(delq * qi, 2)  # Увеличиваем нагрузку
    dict_coef.update({q: qq})
    # print('Увеличиваем нагрузку qq={: f}'.format(qq), " коэффициенты: ", end="")
    delta = 1
    kol_iter = 0
    # print(dict_coef)
    while delta > epsillon:
        dict_coef.update(zip(SN, list(Coef)))
        # print(dict_coef)
        dict_values = dict_coef.values()
        # print(dict_values)
        Deter2 = lambda_deter(*dict_values)
        Jacobi2 = lambda_jacobi(*dict_values)
        Rans = np.dot(np.array(la.inv(Deter2)), Jacobi2).reshape(Coef.shape)
        tmp = Coef - Rans
        Coef = np.array(tmp)  # Находим решение методом Ньютона
        delta = np.sum(np.abs(Coef - XkPred)) / len(Coef)  # ??
        XkPred = np.array(Coef)
        kol_iter = kol_iter + 1
        if kol_iter > 100:
            delta = 0
    # print("kol_iter=", kol_iter, "delta=", delta)
    res2.append(Coef)
    wc1 = W
    Xk_new = list(Coef)
    for wi in range(2 * N, 3 * N):
        wc1 = wc1.subs(SN[wi], Coef[wi])
    wc11 = wc1.subs(x, (aa + aa1) / 2)
    wc = wc11.subs(y, bb / 2)
    WC.append(wc)
    Q_y.append(qq)
    wc2 = W
    Xk_new = list(Coef)
    for wi in range(2 * N, 3 * N):
        wc2 = wc2.subs(SN[wi], Coef[wi])

    wc22 = wc2.subs(x, (aa + aa1) / 4)
    wc23 = wc22.subs(y, bb / 4)
    WC2.append(wc23)

# print(datetime.now() - start_time2)
# fig,ax = plt.figure(num=1, figsize=(8, 6))
# ax.plot(WC, Q_y, color='r', linestyle='--', marker='o', markersize=3, label='W((a+a1)/2,b/2)')
# ax.plot(WC2, Q_y, color='b', linestyle='--', marker='o', markersize=3, label='W((a+a1)/4,b/4)')
# ax.legend(loc='upper left')
# grid1 = plt.grid(True)
# ax.xlabel("W,м")
# ax.ylabel("q,МПа")
# c= strftime("%H_%M_%S", gmtime())
# ax.title('График прогиба W n = '+str(n))
# fig.savefig('График прогиба W n = '+str(n)+' '+str(c))
# ax.show()



fig,ax = plt.subplots(num=1, figsize=(8, 6))
ax.plot(WC, Q_y, color='green', linestyle='-.', marker='o', markersize=3, label='W((a+a1)/2,b/2) N= '+str(N))
ax.plot(WC2, Q_y, color='purple', linestyle='-.', marker='o', markersize=3, label='W((a+a1)/4,b/4) N= '+str(N))
ax.legend(loc='lower center')
grid1 = plt.grid(True)
ax.set_xlabel("W,м")
ax.set_ylabel("q,МПа")
# ax.set_ylim([0, 5.1])
# ax.set_xlim([0, 0.45])
c= strftime("%H_%M_%S", gmtime())
ax.set_title('График прогиба W n = '+str(n))
# fig.savefig('График прогиба W n = '+str(n)+' '+str(c))
fig.show()



print(datetime.now() - start_time2)
remark('Начальные условия')

# fig1, ax1 = plt.subplots(num=1, figsize=(8, 6))
# ax1.plot(WC, Q_y, color='green', linestyle='-.', marker='o', markersize=3, label='W((a+a1)/2,b/2) N= '+str(N))
# ax1.plot(WC2, Q_y, color='purple', linestyle='-.', marker='o', markersize=3, label='W((a+a1)/4,b/4) N= '+str(N))
# ax1.legend(loc='upper left')
# grid1 = plt.grid(True)
# ax1.set_ylim([2.7, 3.4])
# ax1.set_xlim([0, 0.45])
# ax1.set_xlabel("W,м")
# ax1.set_ylabel("q,МПа")
#
# ax1.set_title('График прогиба W n = '+str(n))
# fig1.savefig('График прогиба приближен W n = '+str(n)+' '+str(c))
# fig1.show()


print(datetime.now() - start_time2)
remark('Начальные условия')
