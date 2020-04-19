import math as mth
from math import *
import numpy as np
import sympy as sp
import symengine as sm
from datetime import datetime
import time
#import test as ts
from threading import Thread
import scipy.integrate as integrate
import scipy.special as special
import matplotlib.pyplot as plt
start_time = datetime.now()
class Quadrature:
    """Базовые определения для квадратурных формул"""
    __sum = 0.0
    __nseg = 1  # число отрезков разбиения
    __ncalls = 0 # считает число вызовов интегрируемой функции

    def __restart(func, x0, x1, nseg0, reset_calls = True):
        """Обнуление всех счётчиков и аккумуляторов.
           Возвращает интеграл методом трапеций на начальном разбиении"""
        if reset_calls:
            Quadrature.__ncalls = 0
        Quadrature.__nseg = nseg0
        # вычисление суммы для метода трапеций с начальным числом отрезков разбиения nseg0
        Quadrature.__sum = 0.5 * (func(x0) + func(x1))
        dx = 1.0 * (x1 - x0) / nseg0
        for i in range(1, nseg0):
            Quadrature.__sum += func(x0 + i * dx)

        Quadrature.__ncalls += 1 + nseg0
        return Quadrature.__sum * dx

    def __double_nseg(func, x0, x1):
        """Вдвое измельчает разбиение.
           Возвращает интеграл методом трапеций на новом разбиении"""
        nseg = Quadrature.__nseg
        dx = (x1 - x0) / nseg
        x = x0 + 0.5 * dx
        i = 0
        AddedSum = 0.0
        for i in range(nseg):
            AddedSum += func(x + i * dx)

        Quadrature.__sum += AddedSum
        Quadrature.__nseg *= 2
        Quadrature.__ncalls += nseg
        return Quadrature.__sum * 0.5 * dx

    def simpson(func, x0, x1, rtol = 1.0e-10, nseg0 = 1):
        start_time_simp = datetime.now()
        """Интегрирование методом парабол с заданной точностью.
           rtol - относительная точность,
           nseg0 - число отрезков начального разбиения"""
        old_trapez_sum = Quadrature.__restart(func, x0, x1, nseg0)
        new_trapez_sum = Quadrature.__double_nseg(func, x0, x1)
        ans = (4 * new_trapez_sum - old_trapez_sum) / 3
        old_ans = 0.0
        i=0
        err_est = max(1, abs(ans))
        while (err_est > abs(rtol * ans)):
            i=i+1
            if Quadrature.__ncalls>5000:
                # print("Total function calls: " + str(Quadrature.__ncalls))
                return ans
            else:
                old_ans = ans
                old_trapez_sum = new_trapez_sum
                new_trapez_sum = Quadrature.__double_nseg(func, x0, x1)
                ans = (4 * new_trapez_sum - old_trapez_sum) / 3
                err_est = abs(old_ans - ans)

        # print("Total function calls: " + str(Quadrature.__ncalls))
        # print("Время взятия интегралла")
        # print(datetime.now() - start_time_simp)
        return ans

def Fun1(zxccxz,dict_xx):
    zamena_x = ''
    number = 0
    for elem in zxccxz:
        # print(number)
        for zz in variable_x:
            if elem.find("*"+zz):
                i = elem.find("*"+zz)
                if i < 0:
                    continue
                elem = elem[:i] + '' + elem[i + len(zz) + 1:]
                zamena_x = zamena_x + '*' + zz
                EP[number] = elem.strip(' ')
            continue
        if zamena_x[1:] in dict_xx:
            resul_x = dict_xx.get(zamena_x[1:])
            EP[number] = elem.strip() + "*"+str(resul_x)
            zamena_x = ""
            number = number + 1
            int_x.close()
            continue
        else:

            zam_x = sm.expand(zamena_x[1:])
            resul_x = integrate.quad(lambda xx: (zam_x * A).subs(x, xx), 0, 5.4)[0]
            dict_xx.update({zamena_x[1:]: resul_x})
            EP[number] = elem.strip() +"*"+ str(resul_x)
            zamena_x = ""
            number = number + 1
    return EP

def Fun2(zxccxz,dict_yy):
    zamena_y = ''
    number = 0
    for elem in zxccxz:
        # print(number)
        for zz in variable_y:
            if elem.find("*"+zz):
                i = elem.find("*"+zz)
                if i < 0:
                    continue
                elem = elem[:i] + '' + elem[i + len(zz) + 1:]
                zamena_y = zamena_y + '*' + zz

                EP[number] = elem.strip()
            continue
        if zamena_y[1:] in dict_yy:
            # print("нашел")
            resul_y = dict_yy.get(zamena_y[1:])
            EP[number] = elem.strip() + '*' + str(resul_y)
            zamena_y = ""
            number = number + 1
            int_y.close()
            continue
        else:
            # print("не нашел")
            zam_y = sm.expand(zamena_y[1:])

            resul_y = integrate.quad(lambda yy: (zam_y * B).subs(y, yy), 0, 5.4)[0]

            dict_yy.update({zamena_y[1:]: resul_y})
            EP[number] = elem.strip() + '*' + str(resul_y)
            zamena_y = ""
            number = number + 1
    return EP


h = 0.09

n =2
N = np.power(n, 2)

aa =round(60 * h,2)
aa1 = 0
bb = round(60 * h,2)

E1 =2.1 * 10 ** 5
E2 =2.1 * 10 ** 5
mu12= 0.3
mu21=0.3



z = -(1 / 2) * h

r = 225 * h

k = 5 / 6
f = lambda i: 6 * (1 / 4 - i ** 2 / h ** 2)

A = 1
B = 1

kx =1 / r
ky =1 / r

G12 = 0.33 * 10 ** 5
G13 = 0.33 * 10 ** 5
G23 = 0.33 * 10 ** 5
x = sp.symbols('x')
y = sp.symbols('y')
q = sp.symbols('q')
print(mth.pi)

X1 = lambda i: sp.sin(2 * i * round(mth.pi, 5) * x / aa)
X2 = lambda i: sp.sin((2 * i - 1) * round(mth.pi, 5) * x / aa)
X3 = lambda i: sp.sin((2 * i - 1) * round(mth.pi, 5)*  x / aa)
X4 = lambda i: sp.cos((2 * i - 1) * round(mth.pi, 5) * x / aa)
X5 = lambda i: sp.sin((2 * i - 1) * round(mth.pi, 5) * x / aa)
Y1 = lambda i: sp.sin((2 * i - 1) * round(mth.pi, 5) * y / bb)
Y2 = lambda i: sp.sin(2 * i * round(mth.pi, 5) * y / bb)
Y3 = lambda i: sp.sin((2 * i - 1) * round(mth.pi, 5) * y / bb)
Y4 = lambda i: sp.sin((2 * i - 1) * round(mth.pi, 5) * y / bb)
Y5 = lambda i: sp.cos((2 * i - 1) * round(mth.pi, 5) * y / bb)

U = 0
V = 0
W = 0
Psix = 0
Psiy = 0

u=[]
for l in range(0, n):
    u.append([0]*n)
v=[]
for l in range(0, n):
    v.append([0]*n)
w=[]
for l in range(0, n):
    w.append([0]*n)
psix=[]
for l in range(0, n):
    psix.append([0]*n)
psiy=[]
for l in range(0, n):
    psiy.append([0]*n)


for i in range(1, int(np.sqrt(N)+1)):
    for j in range(1, int(np.sqrt(N)+1)):
        u[i-1][j-1]=(sp.symbols('u' + str(i) + '' + str(j)))
        v[i-1][j-1]=(sp.symbols('v' + str(i) + '' + str(j)))
        w[i-1][j-1]=(sp.symbols('w' + str(i) + '' + str(j)))
        psix[i-1][j-1]=(sp.symbols('psix' + str(i) + '' + str(j)))
        psiy[i-1][j-1]=(sp.symbols('psiy' + str(i) + '' + str(j)))




SN=[]



for i in range(1, int(np.sqrt(N)+1)):
    for j in range(1, int(np.sqrt(N)+1)):
        SN.append(sp.symbols('u' + str(i) + '' + str(j)))
for i in range(1, int(np.sqrt(N)+1)):
    for j in range(1, int(np.sqrt(N)+1)):
        SN.append(sp.symbols('v' + str(i) + '' + str(j)))
for i in range(1, int(np.sqrt(N)+1)):
    for j in range(1, int(np.sqrt(N)+1)):
        SN.append(sp.symbols('w' + str(i) + '' + str(j)))
for i in range(1, int(np.sqrt(N)+1)):
    for j in range(1, int(np.sqrt(N)+1)):
        SN.append(sp.symbols('psix' + str(i) + '' + str(j)))
for i in range(1, int(np.sqrt(N)+1)):
    for j in range(1, int(np.sqrt(N)+1)):
        SN.append(sp.symbols('psiy' + str(i) + '' + str(j)))

for i in range(1, int(np.sqrt(N)+1)):
    for j in range(1, int(np.sqrt(N)+1)):
        U = U + u[i-1][j-1] * X1(i) * Y1(j)
        V = V + v[i-1][j-1] * X2(i) * Y2(j)
        W = W + w[i-1][j-1] * X3(i) * Y3(j)
        Psix = Psix + psix[i-1][j-1] * X4(i) * Y4(j)
        Psiy = Psiy + psiy[i-1][j-1] * X5(i) * Y5(j)

#print(U)
#print(V)
#print(W)
#print(Psix)
#print(Psiy)

Theta1 = sm.expand(-(sp.diff(W, x)) / A - kx * U)
#print("Theta1")
#print(Theta1)



Theta2 =sm.expand( -(sp.diff(W, y)) / B - ky * V)
#print("Theta2")
#print(Theta2)
ex =sm.expand( (sp.diff(U, x)) / A + (sp.diff(A, y)) * V / (A * B) - kx * W + (1 / 2) * Theta1 ** 2)
#print("#ex")
#print(ex)
ey = sm.expand((sp.diff(V, y)) / B + (sp.diff(B, x)) * U / (A * B) - ky * W + (1 / 2) * Theta2 ** 2)
#print("#ey")
#print(ey)
gxy =sm.expand( (sp.diff(V, x)) / A + (sp.diff(U, y)) / B - (sp.diff(A, y)) * U / (A * B) - (sp.diff(B, x)) * V / (A * B) + Theta1 * Theta2)
#print("#gxy")
#print(gxy)
gxz =sm.expand( k * (f(z)) * (Psix - Theta1))
gyz =sm.expand( k * (f(z)) * (Psiy - Theta2))
#print("#gxz")
#print(gxz)
#print("#gyz")
#print(gyz)
varkappa1 =sm.expand( (sp.diff(Psix, x)) / A + (sp.diff(A, y)) * Psiy / (A * B))
varkappa2 =sm.expand( (sp.diff(Psiy, y)) / B + (sp.diff(B, x)) * Psix / (A * B))
varkappa12 =sm.expand( 1 / 2 * ((sp.diff(Psiy, x)) / A + (sp.diff(Psix, y)) / B - ((sp.diff(A, y)) * Psix + (sp.diff(B, x)) * Psiy) / (A * B)))
#print("#varkappa1")

#print(varkappa1)
#print("#varkappa2")
#print(varkappa2)
#print("#varkappa12")
#print(varkappa12)


print("#Mx")


Mx = sm.expand((1 / 12) * E1 * h ** 3 * (mu21 * varkappa2 + varkappa1) / (1-mu12 * mu21))
# print(Mx)


print("#My")


My = sm.expand((1 / 12) * E2 * h ** 3 * (mu12 * varkappa1 + varkappa2) / (1-mu12 * mu21))
# print(My)
print("#Mxy")
Mxy = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)

# print(Mxy)
print("#Myx")
Myx = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)
# print(Myx)
print("#Nx")
Nx = sm.expand((E1 * h / (1-mu12 * mu21)) * (ex + mu21 * ey))
# print(Nx)

print("#Ny")
Ny = sm.expand((E2 * h / (1-mu12 * mu21)) * (ey + mu12 * ex))
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





#
#
#
# Epp=Nx * ex+Ny * ey+1 / 2 * (Nxy + Nyx) * gxy +Mx * varkappa1+My * varkappa2 +(Mxy + Myx) * varkappa12 +Qx * (Psix - Theta1)+Qy * (Psiy - Theta2)
# EPp=sm.expand(Epp)
#

Epp1 = Nx * ex+Ny * ey
print("Razbienie1")
# print(Epp1)
del(Nx , ex,Ny , ey)
Epp3 =Epp1+ 1 / 2 * (Nxy + Nyx) * gxy
del(Nxy , Nyx, gxy)
print("Razbienie3")
Epp4 =Epp3+ Mx * varkappa1+My * varkappa2
print("Razbienie4")
del(Mx ,varkappa1,My , varkappa2)
Epp6 =Epp4+ (Mxy + Myx) * varkappa12
print("Razbienie6")
del(Mxy ,Myx, varkappa12)
Epp7 =Epp6+ Qx * (Psix - Theta1)
print("Razbienie7")
del(Qx ,Psix , Theta1)
Epp8 =Epp7+ Qy * (Psiy - Theta2)
print("Razbienie8")
del(Qy ,Psiy,Theta2)

AllEpp=Epp8
# print(AllEpp)
del(Epp1,Epp3,Epp4,Epp6,Epp7,Epp8)

EPp=sm.expand(AllEpp)

print(integrate.quad(lambda xx: EPp.subs(x, xx), 0, 5.4)[0])

# print("Время раскрытия скобок")
# print(datetime.now() - start_time)
#
#
#
#
# print("Время раскрытия скобок")
# print(datetime.now() - start_time)
# print("Resul1")
#
#
#
#
# allin=0
# for el in EP:
#     allin=allin+1/2*sm.expand(el)
#
# print("Allin")
# allin=sm.expand(allin)
#
#
#
# AA = sp.integrate(sp.integrate((Px * U + Py * V + W * q) * A * B, (y, 0, 5.4)), (x, 0, 5.4))
# #print(AA)
#
# AA=sm.expand(AA)
# Es = allin - AA
#
# Es=sm.expand(Es)
# print("Es")
# print(Es)
