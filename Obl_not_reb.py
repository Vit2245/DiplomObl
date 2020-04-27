import time
from math import *
import sys
from datetime import datetime
from threading import Thread
from typing import Union
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import symengine as sm
import sympy as sp
from numpy import linalg as la
from sympy.parsing.sympy_parser import parse_expr

start_all=datetime.now()


sys.setrecursionlimit(10 ** 6)

start_time = datetime.now()
Num = Union[int, float]


n =4
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

# print(AllEpp)
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
EP = []

for xc in Epp:
    EP.append(str(sm.expand(xc)))

# for i in EP:
#     print(i)
print('Всего прошло времени')
print(datetime.now() - start_all)


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
#
Num = Union[int, float]


def replace_by_dict(ep: [str], variables: [str], dictionary: dict, a: Num, b: Num, s=sp.symbols) -> None:
    for index, elem in enumerate(ep):
        zamena = ''
        # print(index)
        # put the polynom's members in the right order
        for zz in variables:
            if elem.find(zz):
                i = elem.find(zz)
                if i < 0:
                    continue
                if elem[:i]=='*':
                    elem = elem[:i-1] + '' + elem[i-1 + len(zz) +1:]
                else:
                    elem = elem[:i] + '' + elem[i  + len(zz) + 1:]
                zamena += '*' + zz
                ep[index] = elem
            continue
        if zamena[1:] in dictionary:
            # print("нашел")
            # print(zamena[1:])
            result = dictionary.get(zamena[1:])
            # print(result)
            # print(result)
            # ep[index] = elem[:-1] + operator + str(result)
            if elem[-1] == '*':
                ep[index] = elem + '' + str(result)
            else:
                ep[index] = elem + '*' + str(result)
        else:
            # print("не нашел")
            zam_x = sm.expand(zamena[1:])
            # print(zam_x)
            # result=Quadrature.simpson(lambda xx: (zam_x*A).subs(x,xx), 0, 5.4, rtol=1e-10)
            if s == 'x':
                result = integrate.quad(lambda xx: (zam_x * A).subs(x, xx), a, b)[0]
                if result > -0.00001 and result < 0.00001:
                    result = 0
            else:
                result = integrate.quad(lambda yy: (zam_x * B).subs(y, yy), a, b)[0]
                if result > -0.00001 and result < 0.00001:
                    result = 0

            dictionary.update({zamena[1:]: result})

            if elem[-1] == '*':
                ep[index] = elem +'' + str(result)
            else:
                ep[index] = elem + '*' + str(result)

replace_by_dict(EP, variable_x, dict_x, 0, 5.4, 'x')

int_x = open('out_x.txt','w')
with int_x as out:
    for key, val in dict_x.items():
        out.write('{}:{}\n'.format(key, val))

int_x.close()
replace_by_dict(EP, variable_y, dict_y, 0, 5.4, 'y')




print("Время раскрытия скобок3")
print(datetime.now() - start_time)
print('Всего прошло времени')
print(datetime.now() - start_all)
int_y = open('out_y.txt','w')
with int_y as out:
    for key, val in dict_y.items():
        out.write('{}:{}\n'.format(key, val))


int_y.close()
number = 0
zxc=[]
with open('агт.txt', 'w') as out:
    for el in EP:
        out.write(str(sm.expand(el)) + "\n")


for el in EP:
     zxc.append(sm.expand(1/2*sm.expand(el)))
print('Всего прошло времени')
print(datetime.now() - start_all)
dict_all={}

print("*****************")
for i in zxc:
    mas=i.args
    if len(mas)!=0:
        b=1
        for ii in mas[1:]:
            b*=ii

        # print(b,mas[0])
        # dict_all.update({b},mas[0])
        if b in dict_all:
            # print("нашел")
            # print(zamena[1:])
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


print('Всего прошло времени')
print(datetime.now() - start_all)

nnn=len(itt)//16

product=0
start_timeS = datetime.now()



# def Sum1(el,eln,nn,mass):
#     el += reduce((lambda x, y: x + y), mass[(eln-1)*nn:eln * nn])
# thread1 = Thread(target=Sum1, args=(product,2,nnn,itt,))
# thread2 = Thread(target=Sum1, args=(product,3,nnn,itt,))
# thread3 = Thread(target=Sum1, args=(product,4,nnn,itt,))
# thread4 = Thread(target=Sum1, args=(product,5,nnn,itt,))
# thread5 = Thread(target=Sum1, args=(product,6,nnn,itt,))
# thread6 = Thread(target=Sum1, args=(product,7,nnn,itt,))
# thread7 = Thread(target=Sum1, args=(product,8,nnn,itt,))
# thread8 = Thread(target=Sum1, args=(product,9,nnn,itt,))
# thread9 = Thread(target=Sum1, args=(product,10,nnn,itt,))
# thread10 = Thread(target=Sum1, args=(product,11,nnn,itt,))
# thread11 = Thread(target=Sum1, args=(product,12,nnn,itt,))
# thread12 = Thread(target=Sum1, args=(product,13,nnn,itt,))
# thread13 = Thread(target=Sum1, args=(product,14,nnn,itt,))
# thread14 = Thread(target=Sum1, args=(product,15,nnn,itt,))
# thread15 = Thread(target=Sum1, args=(product,16,nnn,itt,))

# thread1.start()
# thread2.start()
# thread3.start()
# thread4.start()
# thread5.start()
# thread6.start()
# thread7.start()
# thread8.start()
# thread9.start()
# thread10.start()
# thread11.start()
# thread12.start()
# thread13.start()
# thread14.start()
# thread15.start()
# print(product)
# thread1.join()
# thread2.join()
# thread3.join()
# thread4.join()
# thread5.join()
# thread6.join()
# thread7.join()
# thread8.join()
# thread9.join()
# thread10.join()
# thread11.join()
# thread12.join()
# thread13.join()
# thread14.join()
# thread15.join()
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
print("Время сумма2")
print(datetime.now() - start_timeS)



product=sm.expand(product)




print("Время Allin")
start_time = datetime.now()
AA = sp.integrate(sp.integrate((Px * U + Py * V + W * q) * A * B, (y, 0, aa)), (x, 0, bb))
# print(AA)

AA = sm.expand(AA)
Es = product - AA
Es = sm.expand(Es)
print("Es")
print(datetime.now() - start_time)
# print(Es)
print('Всего прошло времени')
print(datetime.now() - start_all)
# print(Coef)
start_time1 = datetime.now()
Jacobi = []
for i in SN:
    Jacobi.append(sm.expand(sm.diff(Es, i)))
print(type(Jacobi[0]))
print("Время первая производная")
print(datetime.now() - start_time1)
print('Всего прошло времени')
print(datetime.now() - start_all)
Deter = []
start_time2 = datetime.now()
for dpU in Jacobi:
    lineOfHessian = []
    for symb in SN:
        lineOfHessian.append(sm.diff(dpU, symb))
    Deter.append(lineOfHessian)

print(type(Deter[0][0]))
print("Время вторая производная")
print(datetime.now() - start_time2)


print('Всего прошло времени')
print(datetime.now() - start_all)

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
MAX = 34
Q_y = []


dict_coef = dict(zip(SN, list(Coef)))
dict_coef.update({q: 0.})
print(dict_coef.keys())
start_time = datetime.now()
SN.append(q)

# start_time = datetime.now()
# lambda_deter = sm.lambdify([*dict_coef.keys()], Deter)
#
# print("Время матрицы dict")
# print(datetime.now() - start_time)
# start_time = datetime.now()
# lambda_jacobi = sm.lambdify([*dict_coef.keys()], Jacobi)
#
# print("Время матрицы dict")
# print(datetime.now() - start_time)



start_time = datetime.now()
lambda_deter = sm.lambdify(SN, Deter)

print("Время матрицы SN")
print(datetime.now() - start_time)
start_time = datetime.now()
lambda_jacobi= sm.lambdify(SN, Jacobi)

print("Время матрицы SN")
print(datetime.now() - start_time)


print('Всего прошло времени')
print(datetime.now() - start_all)
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
        dict_values = dict_coef.values()
        Deter2 = lambda_deter(*dict_values)
        Jacobi2 = lambda_jacobi(*dict_values)
        Rans = np.dot(np.array(la.inv(Deter2)), Jacobi2).reshape(Coef.shape)
        tmp = Coef - Rans
        Coef = np.array(tmp)  # Находим решение методом Ньютона
        delta = np.sum(np.abs(Coef - XkPred)) / len(Coef)  # ??
        XkPred = np.array(Coef)
        kol_iter = kol_iter + 1
        if kol_iter > 22:
            delta = 0

    # print("kol_iter=", kol_iter, "delta=", delta)
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

print(datetime.now() - start_time2)
fig = plt.figure(num=1, figsize=(8, 6))
plt.plot(WC, Q_y, color='r', linestyle='--', marker='o', markersize=3, label='W((a+a1)/2,b/2)')
plt.plot(WC2, Q_y, color='b', linestyle='--', marker='o', markersize=3, label='W((a+a1)/4,b/4)')
plt.legend(loc='upper left')
grid1 = plt.grid(True)
plt.xlabel("W,м")
plt.ylabel("q,МПа")
plt.title('График прогиба W n = '+str(n))
print(datetime.now() - start_time2)
print('Всего прошло времени')
print(datetime.now() - start_all)
plt.show()