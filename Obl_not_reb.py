import math as mth
from math import *
import numpy as np
from sympy import *
import symengine as sm
from datetime import datetime
import time

from scipy.integrate import quad
start_time = datetime.now()

h = 0.09

n =1
N = np.power(n, 2)

a =60 * h
a1 = 0
b = 60 * h

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
x = symbols('x')
y = symbols('y')
q = symbols('q')
print(mth.pi)
X1 = lambda i: sin(2 * i * round(mth.pi, 5) * x / a)
X2 = lambda i: sin((2 * i - 1) * round(mth.pi, 5) * x / a)
X3 = lambda i: sin((2 * i - 1) * round(mth.pi, 5)*  x / a)
X4 = lambda i: cos((2 * i - 1) * round(mth.pi, 5) * x / a)
X5 = lambda i: sin((2 * i - 1) * round(mth.pi, 5) * x / a)
Y1 = lambda i: sin((2 * i - 1) * round(mth.pi, 5) * y / b)
Y2 = lambda i: sin(2 * i * round(mth.pi, 5) * y / b)
Y3 = lambda i: sin((2 * i - 1) * round(mth.pi, 5) * y / b)
Y4 = lambda i: sin((2 * i - 1) * round(mth.pi, 5) * y / b)
Y5 = lambda i: cos((2 * i - 1) * round(mth.pi, 5) * y / b)

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
        u[i-1][j-1]=(Symbol('u' + str(i) + '' + str(j)))
        v[i-1][j-1]=(Symbol('v' + str(i) + '' + str(j)))
        w[i-1][j-1]=(Symbol('w' + str(i) + '' + str(j)))
        psix[i-1][j-1]=(Symbol('psix' + str(i) + '' + str(j)))
        psiy[i-1][j-1]=(Symbol('psiy' + str(i) + '' + str(j)))



#print(u)
#print(v)
#print(w)
#print(psix)
#print(psiy)

SN=[]
SN.append(u)
SN.append(v)
SN.append(w)
SN.append(psix)
SN.append(psiy)
#print(SN)




for i in range(1, int(np.sqrt(N)+1)):
    for j in range(1, int(np.sqrt(N)+1)):
        U = U + u[i-1][j-1] * X1(i) * Y1(j)
        V = V + v[i-1][j-1] * X2(i) * Y2(j)
        W = W + w[i-1][j-1] * X3(i) * Y3(j)
        Psix = Psix + psix[i-1][j-1] * X4(i) * Y4(j)
        Psiy = Psiy + psiy[i-1][j-1] * X5(i) * Y5(j)

print(U)
print(V)
print(W)
print(Psix)
print(Psiy)


Theta1 = sm.expand(-(diff(W, x)) / A - kx * U)
print("Theta1")
#print(Theta1)



Theta2 =sm.expand( -(diff(W, y)) / B - ky * V)
print("Theta2")
#print(Theta2)
ex =sm.expand( (diff(U, x)) / A + (diff(A, y)) * V / (A * B) - kx * W + (1 / 2) * Theta1 ** 2)
print("#ex")
#print(ex)
ey = sm.expand((diff(V, y)) / B + (diff(B, x)) * U / (A * B) - ky * W + (1 / 2) * Theta2 ** 2)
print("#ey")
#print(ey)
gxy =sm.expand( (diff(V, x)) / A + (diff(U, y)) / B - (diff(A, y)) * U / (A * B) - (diff(B, x)) * V / (A * B) + Theta1 * Theta2)
print("#gxy")
#print(gxy)
gxz =sm.expand( k * (f(z)) * (Psix - Theta1))
gyz =sm.expand( k * (f(z)) * (Psiy - Theta2))
print("#gxz")
#print(gxz)
print("#gyz")
#print(gyz)
varkappa1 =sm.expand( (diff(Psix, x)) / A + (diff(A, y)) * Psiy / (A * B))
varkappa2 =sm.expand( (diff(Psiy, y)) / B + (diff(B, x)) * Psix / (A * B))

varkappa12 =expand( 1 / 2 * ((diff(Psiy, x)) / A + (diff(Psix, y)) / B - ((diff(A, y)) * Psix + (diff(B, x)) * Psiy) / (A * B)))
print("#varkappa1")

#print(varkappa1)
print("#varkappa2")
#print(varkappa2)
print("#varkappa12")
#print(varkappa12)


Mx = sm.expand((1 / 12) * E1 * h ** 3 * (mu21 * varkappa2 + varkappa1) / (1-mu12 * mu21))
My = sm.expand((1 / 12) * E2 * h ** 3 * (mu12 * varkappa1 + varkappa2) / (1-mu12 * mu21))
Mxy = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)
Myx = sm.expand((1 / 6) * G12 * h ** 3 * varkappa12)
Nx = sm.expand((E1 * h / (1-mu12 * mu21)) * (ex + mu21 * ey))
Ny = sm.expand((E2 * h / (1-mu12 * mu21)) * (ey + mu12 * ex))

print( sm.expand((E1 * h / (1-mu12 * mu21)) * (ex + mu21 * ey)))



print("#Mx")

#print(Mx)
print("#My")
#print(My)
print("#Mxy")
#print(Mxy)
print("#")

#print(Myx)
print("#Nx")
print(Nx)
print("#Ny")
print(Ny)



print("#")
Nxy = sm.expand(G12 * h * gxy)
Nyx = sm.expand(G12 * h * gxy)
Px = 0
Py = 0
Qx = sm.expand(G13 * k * h * (Psix - Theta1))
Qy = sm.expand(G23 * k * h * (Psiy - Theta2))
print("#")


"""


Epp=Nx * ex+Ny * ey+1 / 2 * (Nxy + Nyx) * gxy +Mx * varkappa1+My * varkappa2 +(Mxy + Myx) * varkappa12 +Qx * (Psix - Theta1)+Qy * (Psiy - Theta2)
EPp=sm.expand(Epp)
"""

Epp1 = Nx * ex+Ny * ey
print("Razbienie1")
print(Epp1)
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
print(AllEpp)
del(Epp1,Epp3,Epp4,Epp6,Epp7,Epp8)

EPp=sm.expand(AllEpp)



print(datetime.now() - start_time)


print("EPp")

print(EPp)

##print(type(EPp))
#print("split")
#Epp=str(EPp).split('+')
#print("Epp")
#print(Epp)
#print("len(Epp)")
##print(len(Epp))
##print(Epp[0])
#print("Epp")
#
#
#EP=0
#zxc=0
"""
for i in Epp:
    print(zxc)
    print(i)
    EP=EP+1 / 2 * quad(quad(i*A*B,  0, b),  a1, a)

    print(EP)
    zxc=zxc+1


print("EP")
print(sm.expand(EP))



AA = integrate(integrate((Px * U + Py * V + W * q) * A * B, (y, 0, b)), (x, a1, a))
print("#")


Es = EP - AA
Es=expand(Es)
print("Es")
print(Es)



Coef=[]
for l in range(0, 5 * N):
    Coef.append(0)
#print(Coef)

Jacobi = [0] * 5 * N

k = 0
for i in range(0, n ):
    for j in range(0, n ):
        Jacobi[k] = diff(Es,u[i][j])
        Jacobi[k + N] = diff(Es, v[i][j])
        Jacobi[k + 2 * N] = diff(Es, w[i][j])
        Jacobi[k + 3 * N] = diff(Es, psix[i][j])
        Jacobi[k + 4 * N] = diff(Es, psiy[i][j])
        k = k + 1
print("Jacobi")
for i in Jacobi:
    print(i)
print("#+#+")
Deter=[]

for l in range(0, 5 * N):
    Deter.append([0]*5)

#for b in Deter:
#    print(b)

print("!@#$%")
for l in range(0, 5 * N):
    k = 0
    for i in range(0, n ):
        for j in range(0, n ):
            Deter[l][k] = diff(Jacobi[l], u[i][j])
            Deter[l][k + N] = diff(Jacobi[l],v[i][j])
            Deter[l][k + N * 2] = diff(Jacobi[l], w[i][j])
            Deter[l][k + N * 3] = diff(Jacobi[l], psix[i][j])
            Deter[l][k + N * 4] = diff(Jacobi[l], psiy[i][j])
            k = k + 1

print("3")
for b in Deter:
    print(b)
print(Deter)

Prob3 = []
for l in range(0, 5 * N):
    Prob3.append([0]*5)

MAX = 340
epsillon=10**-5
delq = 0.1e-1
qq = 0


AnsMatr=np.zeros((MAX+1, 5*N+5),dtype=float)
Jacobi1 = [0] * 5 * N
Deter1=[]

for l in range(0, 5 * N):
    Deter1.append([0]*5)
#print("@")

bb=np.zeros((5*N),dtype=float)
cc=np.zeros((5*N, 5*N),dtype=float)
Coef=np.zeros((5*N),dtype=float)
BufV=np.zeros((5*N),dtype=float)
Buf = np.zeros((5*N),dtype=float)
#print(Coef)
#print(Buf)
#print(BufV)
print("12345")
#print(AnsMatr)
Resul = [0] * (MAX+25)

print("Основной цикл")
for p in range(0, MAX+1):
    delt = 1
    mm=0
    while delt >epsillon:
        mm=mm+1
        for l in  range(0, 5*N):
            BufV[l]=Coef[l]
        k=0
        for l in range(0, len(Jacobi)):
            k = 0
            for i in range(0, n):
                for j in range(0, n):
                    bb[l] = Jacobi[l].subs([(Symbol('u' + str(i + 1) + '' + str(j + 1)), Coef[k]),
                                        (Symbol('v' + str(i + 1) + '' + str(j + 1)), Coef[k + N]),
                                        (Symbol('w' + str(i + 1) + '' + str(j + 1)), Coef[k + N * 2]),
                                        (Symbol('psix' + str(i + 1) + '' + str(j + 1)), Coef[k + N * 3]),
                                        (Symbol('psiy' + str(i + 1) + '' + str(j + 1)), Coef[k + N * 4]),(q, p/100)])
                    k = k + 1
        for ii in range(len(Deter)):
            for jj in range(len(Deter[ii])):
                k = 0
                for i in range(0, n):
                    for j in range(0, n):
                        cc[ii][jj] = Deter[ii][jj].subs([(Symbol('u' + str(i + 1) + '' + str(j + 1)), Coef[k]),
                                                    (Symbol('v' + str(i + 1) + '' + str(j + 1)), Coef[k + N]),
                                                    (Symbol('w' + str(i + 1) + '' + str(j + 1)), Coef[k + N * 2]),
                                                    (Symbol('psix' + str(i + 1) + '' + str(j + 1)), Coef[k + N * 3]),
                                                    (Symbol('psiy' + str(i + 1) + '' + str(j + 1)), Coef[k + N * 4]),
                                                    (q, p/100)])
                        k = k + 1

        Rans = np.linalg.inv(cc).dot(bb)
        #print(Rans)
        Buf = Coef
        Coef = Buf - Rans
        delt = abs(BufV[1] - Coef[1])

        for l in  range(0, 5*N):
            if abs(BufV[l]-Coef[l])>delt:
                delt=abs(BufV[l]-Coef[l])
    print("Jac")
    print(bb)
    print("Det")
    print(cc)
    print(Rans)
    for l in range(0, 5 * N):
        AnsMatr[p][l + 1]= Coef[l]
            
    AnsMatr[p][0] = p/100


    for i in range(0, sqrt(N)):
        for j in range(0, sqrt(N)):
            Resul[p] = W.subs([(x, a/2),(y, a/2),(Symbol('w' + str(i + 1) + '' + str(j + 1)), Coef[k + N * 2])])



print("123456789")
for i in Resul:
    print(i)

for i in AnsMatr:
    print(i)


#EPN = EPN.subs([(Symbol('u' + str(1) + '' + str(1)), 0),(Symbol('v' + str(i) + '' + str(j)), 0),(Symbol('w' + str(i) + '' + str(j)), 0),(Symbol('psix' + str(i) + '' + str(j)), 0),(Symbol('psiy' + str(i) + '' + str(j)), 0)])



            
            del := abs(evalf(BufV[1]-Coef[1]));
            for l to 5*N do
                if abs(evalf(BufV[l]-Coef[l])) > del then del := abs(evalf(BufV[l]-Coef[l]))
                end if
            end do
        end do;
        for l to 5*N do
            AnsMatr[p, l+1] := Coef[l]
        end do;
        AnsMatr[p, 1] := qq;
        AnsMatr[p, 2] := subs({x = (1/2)*a, y = (1/2)*b}, W);
        AnsMatr[p, 3] := subs({x = (1/4)*a, y = (1/4)*b}, W);
        qq := qq+delq;
        print(qq)
end do;

evalm(AnsMatr);
with(plots);

gr1 := pointplot([seq([AnsMatr[i, 2], AnsMatr[i, 1]], i = 1 .. MAX)], color = red, axis = [gridlines = [10, color = black]], labels = ["W", "q"]);
gr2 := pointplot([seq([AnsMatr[i, 3], AnsMatr[i, 1]], i = 1 .. MAX)], color = blue, axis = [gridlines = [10, color = black]], labels = ["W", "q"]);
print(display([gr1, gr2]));

NULL;
"""
