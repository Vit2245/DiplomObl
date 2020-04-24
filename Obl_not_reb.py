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
start_all=datetime.now()


sys.setrecursionlimit(10 ** 6)

start_time = datetime.now()
Num = Union[int, float]




h = 0.09

n =3
N = np.power(n, 2)

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
EPp = sm.expand(AllEpp)
# print(EPp)
print("Время раскрытия скобок")
print(datetime.now() - start_time)

Epp=EPp.args
EP = []
strr='0'
for xc in Epp:
    EP.append(str(sm.expand(xc)))

for i in EP:
    print(i)



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


def replace_by_dict(ep: [str], variables: [str], dictionary: dict, a: Num, b: Num, s=sp.symbols,
                    operator: str = '') -> None:
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
xc=EP[ 0:19]
print(len(xc))
for i in xc:
    print(i)
replace_by_dict(EP, variable_x, dict_x, 0, 5.4, 'x')
# print("-----------")
# for i in EP:
#     print(sm.expand(i))
# print("Время раскрытия скобок")
# print(datetime.now() - start_time)
#

int_x = open('out_x.txt','w')
with int_x as out:
    for key, val in dict_x.items():
        out.write('{}:{}\n'.format(key, val))

int_x.close()
replace_by_dict(EP, variable_y, dict_y, 0, 5.4, y, '*')
# print("-----------")
# for i in EP:
#     print(sm.expand(i))

# thread1 = Thread(target=Fun1, args=(EP,dict_x))
# thread2 = Thread(target=Fun2, args=(EP,dict_y))
#
#
#
# thread1.start()
# time.sleep(0.10)
# thread2.start()
# thread1.join()
# thread2.join()


print("Время раскрытия скобок")
print(datetime.now() - start_time)

int_y = open('out_y.txt','w')
with int_y as out:
    for key, val in dict_y.items():
        out.write('{}:{}\n'.format(key, val))


int_y.close()
number = 0
for el in EP:
     print(sm.expand(el))


EPP1=[]
start_time = datetime.now()
for i, el in enumerate(EP):
    EPP1.append(str( 1 / 2 * sm.expand(el)))


print(EPP1)

print(len(EPP1))
def Sum(massiv):
    aa = '+'.join(massiv)
    a=aa.replace('+-','-')
    otvet=sm.expand(a)
    return otvet

# bbbb='+'.join(EPP1)
# bbbb=bbbb.replace('+-','-')

print("Время Allin")
start_time = datetime.now()
allin = sm.expand(Sum(EPP1))
print(allin)
all=sm.expand(allin)
print(datetime.now() - start_time)

AA = sp.integrate(sp.integrate((Px * U + Py * V + W * q) * A * B, (y, 0, aa)), (x, 0, bb))
print(AA)
start_time = datetime.now()
AA = sm.expand(AA)
Es = all - AA
Es = sm.expand(Es)
print("Es")
print(datetime.now() - start_time)
print(Es)

# print(Coef)
start_time1 = datetime.now()
Jacobi = []
for i in SN:
    Jacobi.append(sm.diff(Es, i))
print(Jacobi)
print("Время первая производная")
print(datetime.now() - start_time1)

Deter = []
start_time2 = datetime.now()
for dpU in Jacobi:
    lineOfHessian = []
    for symb in SN:
        lineOfHessian.append(sp.diff(dpU, symb))
    Deter.append(lineOfHessian)


print("Время вторая производная")
print(datetime.now() - start_time2)

# start_time = datetime.now()
# Jacobi1 = sp.Matrix(Jacobi)
# Deter1 = sp.Matrix(Deter)
# print("Время матрицы")
# print(datetime.now() - start_time)

print('1')
epsillon = 1 * 10 ** (-5)
print('2')
print('Начальный нулевой вектор ... ')
Coef = np.zeros(len(SN), dtype=np.float)
# print(Coef)
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

start_time = datetime.now()
lambda_deter = sp.lambdify(dict_coef.keys(), Deter)
print("Время матрицы")
print(datetime.now() - start_time)
start_time = datetime.now()
lambda_jacobi = sp.lambdify(dict_coef.keys(), Jacobi)
print("Время матрицы")
print(datetime.now() - start_time)

start_time2 = datetime.now()
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
        Deter2 = lambda_deter(*dict_values)
        Jacobi2 = lambda_jacobi(*dict_values)
        Rans = np.dot(np.array(la.inv(Deter2)), Jacobi2).reshape(Coef.shape)
        tmp = Coef - Rans
        Coef = np.array(tmp)  # Находим решение методом Ньютона
        delta = np.sum(np.abs(Coef - XkPred)) / len(Coef)  # ??
        XkPred = np.array(Coef)
        kol_iter = kol_iter + 1
        if kol_iter > 16:
            delta = 0

    print("kol_iter=", kol_iter, "delta=", delta)
    wc1 = W
    Xk_new = list(Coef)
    for wi in range(2 * N, 3 * N):
        wc1 = wc1.subs(SN[wi], Coef[wi])
    # wc1=wc1
    # wcWW.append(wc1)  # масив значений функции W c подставленными коэф. с в завимости от q
    wc11 = wc1.subs(x, (aa + aa1) / 2)
    wc = wc11.subs(y, bb / 2)
    WC.append(wc)
    Q_y.append(qq)
    wc2 = W
    Xk_new = list(Coef)
    for wi in range(2 * N, 3 * N):
        wc2 = wc2.subs(SN[wi], Coef[wi])
    # wc1=wc1
    # wcWW.append(wc1)  # масив значений функции W c подставленными коэф. с в завимости от q
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
plt.title('График прогиба W')
print(datetime.now() - start_time)
print(datetime.now() - start_all)
plt.show()
#

