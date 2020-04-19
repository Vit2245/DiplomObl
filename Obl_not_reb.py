import math as mth
from math import *
import numpy as np
import sympy as sp
from sympy.integrals.quadrature import gauss_legendre
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

n =3
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



print("Время раскрытия скобок")
print(datetime.now() - start_time)
# print("EPp")
# print(EPp)
# print(type(EPp))
# print("split")
Epp=str(EPp).split('+')
# print("Epp")
# print(Epp)

del(Epp[0])


variable_x=[]
for i in range(1,n+1):
    for st in range(4, 0, -1):
        variable_x.append(str(sp.sin(2 * i * round(mth.pi, 5) * x / aa)**st))
        variable_x.append(str(sp.cos(2 * i * round(mth.pi, 5) * x / aa)**st))
        variable_x.append(str(sp.sin((2 * i - 1) * round(mth.pi, 5) * x / aa)**st))
        variable_x.append(str(sp.cos((2 * i - 1) * round(mth.pi, 5) * x / aa)**st))
variable_y = []
for i in range(1,n+1):
    for st in range(4, 0, -1):
        variable_y.append(str(sp.sin(2 * i * round(mth.pi, 5) * y / bb)**st))
        variable_y.append(str(sp.cos(2 * i * round(mth.pi, 5) * y / bb)**st))
        variable_y.append(str(sp.sin((2 * i - 1) * round(mth.pi, 5) * y / bb)**st))
        variable_y.append(str(sp.cos((2 * i - 1) * round(mth.pi, 5) * y / bb)**st))
number=0

EP=[]

for xc in Epp:
    Epp=xc.split('-')
    if len(Epp)>1:

        EP.append(Epp[0])
        del(Epp[0])
        for el in Epp:
            EP.append('-'+el)
    else:
        EP.append(Epp[0])

print(len(EP))
for i in EP:
    print(i)


my_dict={}
dict_x={}
dict_y={}



int_x=open('out_x.txt')

int_y=open('out_y.txt')
with int_x as inp:
    for i in inp.readlines():
        key,val = i.strip().split(':')
        dict_x[key] = val.strip(' ')
with int_y as inp:
    for i in inp.readlines():
        key,val = i.strip().split(':')
        dict_y[key] = val.strip(' ')
int_x.close()
int_y.close()
zamena_x=''
zamena_y=''
number=0

for elem in EP:
    # print(number)
    for zz in variable_x:
        if elem.find(zz):
            i = elem.find(zz)
            if i < 0:
                continue
            elem = elem[:i] + '' + elem[i + len(zz)+1:]
            zamena_x=zamena_x+'*'+zz
            EP[number]=elem
        continue
    if zamena_x[1:] in dict_x:
        # print("нашел")
        resul_x = dict_x.get(zamena_x[1:])
        # print(resul_x)
        EP[number] = elem+str(resul_x)
        zamena_x = ""
        number = number + 1
        int_x.close()
        continue
    else:
        # print("не нашел")
        zam_x=sm.expand(zamena_x[1:])
        # print(zam_x)
        # resul_x=Quadrature.simpson(lambda xx: (zam_x*A).subs(x,xx), 0, 5.4, rtol=1e-10)
        resul_x= integrate.quad(lambda xx: (zam_x*A).subs(x,xx), 0, 5.4)[0]
        dict_x.update({zamena_x[1:]:resul_x})

        EP[number] = elem+str(resul_x)
        zamena_x=""

        number=number+1
print("Время раскрытия скобок")
print(datetime.now() - start_time)

with open('out_x.txt','w') as out:
    for key,val in dict_x.items():
        out.write('{}:{}\n'.format(key,val))

number=0
for elem in EP:
    # print(number)
    for zz in variable_y:
        if elem.find(zz):
            i = elem.find(zz)
            if i < 0:
                continue
            elem = elem[:i] + '' + elem[i + len(zz)+1:]
            zamena_y=zamena_y+'*'+zz
            EP[number]=elem
        continue
    if zamena_y[1:] in dict_y:
        # print("нашел")
        resul_y = dict_y.get(zamena_y[1:])

        EP[number] = elem+'*'+str(resul_y)
        zamena_y= ""
        number = number + 1
        int_y.close()
        continue
    else:
        # print("не нашел")
        zam_y=sm.expand(zamena_y[1:])
        # resul_y=Quadrature.simpson(lambda yy: (zam_y*B).subs(y,yy), 0, 5.4, rtol=1e-10)
        resul_y= integrate.quad(lambda yy: (zam_y*B).subs(y,yy), 0, 5.4)[0]
        dict_y.update({zamena_y[1:]:resul_y})

        EP[number] = elem+'*'+str(resul_y)
        zamena_y=""
        number=number+1






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

with open('out_x.txt','w') as out:
    for key,val in dict_x.items():
        out.write('{}:{}\n'.format(key,val))
with open('out_y.txt','w') as out:
    for key,val in dict_y.items():
        out.write('{}:{}\n'.format(key,val))

int_x.close()
int_y.close()
number=0
print("Resul")
for el in EP:
     print(sm.expand(el))

print("Время раскрытия скобок")
print(datetime.now() - start_time)
print("Resul1")

with open('агт.txt','w') as out:
    for el in EP:
        out.write(str(sm.expand(el))+"\n")


allin=0
for el in EP:
    allin=allin+1/2*sm.expand(el)

print("Allin")
allin=sm.expand(allin)



AA = sp.integrate(sp.integrate((Px * U + Py * V + W * q) * A * B, (y, 0, 5.4)), (x, 0, 5.4))
#print(AA)

AA=sm.expand(AA)
Es = allin - AA

Es=sm.expand(Es)
print("Es")
print(Es)

#print(Coef)
Jacobi2=np.array([0] * 5 * N)
print(Jacobi2)


Jacobi = [0] * 5 * N


k = 0
# for i in range(0, n ):
#     for j in range(0, n ):
#         Jacobi[k] =         sm.expand(sp.diff(Es, u[i][j]))
#         Jacobi[k + N] =     sm.expand(sp.diff(Es, v[i][j]))
#         Jacobi[k + 2 * N] = sm.expand(sp.diff(Es, w[i][j]))
#         Jacobi[k + 3 * N] = sm.expand(sp.diff(Es, psix[i][j]))
#         Jacobi[k + 4 * N] = sm.expand(sp.diff(Es, psiy[i][j]))
#         k = k + 1
# print("Jacobi")
# print(Jacobi[0])
Jacobi=[]

for i in SN:
    Jacobi.append(sm.expand(Es.diff(i)))



Deter=[]

for i in range(0, len(Jacobi)):
    dpU = Jacobi[i]
    for columnsOfHessian in range(0, len(Jacobi)):
        lineOfHessian = []
        for symb in SN:
            lineOfHessian.append(sm.expand(dpU.diff(symb)))
    Deter.append(lineOfHessian)

Jacobi=sp.Matrix(Jacobi)
Deter=sp.Matrix(Deter)

print('Начальный нулевой вектор ... ', end='')


epsillon=1*10**(-5)


Coef=np.zeros(len(SN),dtype=np.float)
# print(Coef)
XkPred=np.array(Coef)

MasRes = []
res2=[]

WC=[]
WCC=[]
wcWW=[]
WC2=[]

BufV=np.zeros((5*N),dtype=float)
Buf = np.zeros((5*N),dtype=float)



delq = 0.01
MAX = 330
Q_y=[]



dict_coef=dict(zip(SN, list(Coef)))

# def Fun3(mat,dict_zam):
#     Deter1=mat.evalf(subs=dict_zam)
#     return Deter1
#
# def Fun4(mat,dict_zam):
#     Jacobi1=np.array(mat.evalf(subs=dict_zam))
#     return Jacobi1



dict_coef=dict(zip(SN, list(Coef)))
for qi in range(0,MAX+1):
        qq=round(delq*qi,2) # Увеличиваем нагрузку
        dict_coef.update({q: qq})
        print('Увеличиваем нагрузку qq={: f}'.format(qq)," коэффициенты: ", end="")
        delta=1
        kol_iter=0
        print(dict_coef)
        while delta>epsillon:

            dict_coef.update(zip(SN, list(Coef)))

            Deter1 = Deter.evalf(subs=dict_coef)
            Jacobi1 = np.array(Jacobi.evalf(subs=dict_coef))

            Rans=np.dot(np.array(Deter1.inv()), Jacobi1).reshape(Coef.shape)
            tmp=Coef-Rans
            Coef=np.array(tmp)	# Находим решение методом Ньютона

            delta= np.sum(np.abs(Coef - XkPred))/len(Coef)# ??

            XkPred=np.array(Coef)

            kol_iter=kol_iter+1
            if kol_iter>15:
                delta=0

        print("kol_iter=",kol_iter,"delta=",delta)
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

print(datetime.now() - start_time)
fig=plt.figure(num=1, figsize=(8, 6))
plt.plot(WC,Q_y,color='r', linestyle='--', marker='o', markersize=3, label='W((a+a1)/2,b/2)')
plt.plot(WC2,Q_y,color='b', linestyle='--', marker='o', markersize=3, label='W((a+a1)/4,b/4)')
plt.legend(loc='upper left')
grid1 = plt.grid(True)
plt.xlabel("W,м")
plt.ylabel("q,МПа")
plt.title('График прогиба W')
plt.show()
#
print(datetime.now() - start_time)