import numpy as np
import math as mth
from sympy import *
import symengine as sm
h = 0.09

n =1
N = np.power(n, 2)

a =round(60 * h,2)
a1 = 0
b = round(60 * h,2)


x = symbols('x')
y = symbols('y')

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

e=['0.0 ', ' 4.74499410004808*psix11**2*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 ', ' 2475.0*psix11**2*sin(0.581775925925926*y)**2*cos(0.581775925925926*x)**2 ', ' 0.678534156306875*psix11**2*cos(0.581775925925926*y)**2*cos(0.581775925925926*x)**2 ', ' 4.74499410004808*psiy11**2*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 ', ' 2475.0*psiy11**2*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 ', ' 0.678534156306875*psiy11**2*cos(0.581775925925926*y)**2*cos(0.581775925925926*x)**2 ', ' 6.03566529492454*u11**2*sin(1.16355185185185*x)**2*sin(0.581775925925926*y)**2 ', ' 1005.2357871213*u11**2*sin(1.16355185185185*x)**2*cos(0.581775925925926*y)**2 ', ' 28118.4835558404*u11**2*cos(1.16355185185185*x)**2*sin(0.581775925925926*y)**2 ', ' 0.0308787925851719*u11**4*sin(1.16355185185185*x)**4*sin(0.581775925925926*y)**4 ', ' 6.03566529492454*v11**2*sin(1.16355185185185*y)**2*sin(0.581775925925926*x)**2 ', ' 1005.2357871213*v11**2*sin(1.16355185185185*y)**2*cos(0.581775925925926*x)**2 ', ' 28118.4835558404*v11**2*cos(1.16355185185185*y)**2*sin(0.581775925925926*x)**2 ', ' 0.0308787925851719*v11**4*sin(1.16355185185185*y)**4*sin(0.581775925925926*x)**4 ', ' 131.687242798354*w11**2*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 ', ' 837.696489267749*w11**2*sin(0.581775925925926*y)**2*cos(0.581775925925926*x)**2 ', ' 837.696489267749*w11**2*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 ', ' 594.817044400514*w11**4*sin(0.581775925925926*y)**4*cos(0.581775925925926*x)**4 ', ' 594.817044400514*w11**4*cos(0.581775925925926*y)**4*sin(0.581775925925926*x)**4 ', ' 2.84699646002885*psiy11*psix11*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 ', ' 1.35706831261375*psiy11*psix11*cos(0.581775925925926*y)**2*cos(0.581775925925926*x)**2 ', ' 58.9326673935725*u11**3*sin(1.16355185185185*x)**2*cos(1.16355185185185*x)*sin(0.581775925925926*y)**3 ', ' 58.9326673935725*v11**3*sin(1.16355185185185*y)**2*cos(1.16355185185185*y)*sin(0.581775925925926*x)**3 ', ' 2879.79083333334*w11*psix11*sin(0.581775925925926*y)**2*cos(0.581775925925926*x)**2 ', ' 2879.79083333334*w11*psiy11*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 - 451.284303982624*w11**3*sin(0.581775925925926*y)*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**3 - 451.284303982624*w11**3*sin(0.581775925925926*y)**3*sin(0.581775925925926*x)*cos(0.581775925925926*x)**2 ', ' 244.444444444444*u11*psix11*sin(1.16355185185185*x)*sin(0.581775925925926*y)**2*cos(0.581775925925926*x) ', ' 244.444444444444*v11*psiy11*sin(1.16355185185185*y)*cos(0.581775925925926*y)*sin(0.581775925925926*x)**2 ', ' 142.211893004115*w11*u11*sin(1.16355185185185*x)*sin(0.581775925925926*y)**2*cos(0.581775925925926*x) - 3102.8049382716*w11*u11*cos(1.16355185185185*x)*sin(0.581775925925926*y)**2*sin(0.581775925925926*x) - 3.25153685921861*w11*u11**2*sin(1.16355185185185*x)**2*sin(0.581775925925926*y)**3*sin(0.581775925925926*x) ', ' 1.45512758996475*w11*u11**3*sin(1.16355185185185*x)**3*sin(0.581775925925926*y)**4*cos(0.581775925925926*x) ', ' 142.211893004115*w11*v11*sin(1.16355185185185*y)*cos(0.581775925925926*y)*sin(0.581775925925926*x)**2 - 3102.8049382716*w11*v11*cos(1.16355185185185*y)*sin(0.581775925925926*y)*sin(0.581775925925926*x)**2 - 3.25153685921861*w11*v11**2*sin(1.16355185185185*y)**2*sin(0.581775925925926*y)*sin(0.581775925925926*x)**3 ', ' 1.45512758996475*w11*v11**3*sin(1.16355185185185*y)**3*cos(0.581775925925926*y)*sin(0.581775925925926*x)**4 ', ' 8179.32840316598*w11**2*u11*cos(1.16355185185185*x)*sin(0.581775925925926*y)**3*cos(0.581775925925926*x)**2 ', ' 25.7142053551352*w11**2*u11**2*sin(1.16355185185185*x)**2*sin(0.581775925925926*y)**4*cos(0.581775925925926*x)**2 ', ' 8179.32840316598*w11**2*v11*cos(1.16355185185185*y)*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**3 ', ' 25.7142053551352*w11**2*v11**2*sin(1.16355185185185*y)**2*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**4 ', ' 201.958726004098*w11**3*u11*sin(1.16355185185185*x)*sin(0.581775925925926*y)**4*cos(0.581775925925926*x)**3 ', ' 201.958726004098*w11**3*v11*sin(1.16355185185185*y)*cos(0.581775925925926*y)**3*sin(0.581775925925926*x)**4 ', ' 697.125576037404*w11**4*sin(0.581775925925926*y)**2*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**2*cos(0.581775925925926*x)**2 ', ' 2010.47157424259*v11*u11*sin(1.16355185185185*x)*sin(1.16355185185185*y)*cos(0.581775925925926*y)*cos(0.581775925925926*x) ', ' 16871.0901335042*v11*u11*cos(1.16355185185185*x)*cos(1.16355185185185*y)*sin(0.581775925925926*y)*sin(0.581775925925926*x) ', ' 17.6798002180718*v11*u11**2*sin(1.16355185185185*x)**2*cos(1.16355185185185*y)*sin(0.581775925925926*y)**2*sin(0.581775925925926*x) ', ' 17.6798002180718*v11**2*u11*sin(1.16355185185185*y)**2*cos(1.16355185185185*x)*sin(0.581775925925926*y)*sin(0.581775925925926*x)**2 ', ' 0.0361899449098215*v11**2*u11**2*sin(1.16355185185185*x)**2*sin(1.16355185185185*y)**2*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 ', ' 1388.5670891773*w11*u11**2*sin(1.16355185185185*x)*cos(1.16355185185185*x)*sin(0.581775925925926*y)**3*cos(0.581775925925926*x) ', ' 99.2825468761774*w11*u11**2*sin(1.16355185185185*x)**2*sin(0.581775925925926*y)*cos(0.581775925925926*y)**2*sin(0.581775925925926*x) ', ' 1388.5670891773*w11*v11**2*sin(1.16355185185185*y)*cos(1.16355185185185*y)*cos(0.581775925925926*y)*sin(0.581775925925926*x)**3 ', ' 99.2825468761774*w11*v11**2*sin(1.16355185185185*y)**2*sin(0.581775925925926*y)*sin(0.581775925925926*x)*cos(0.581775925925926*x)**2 - 76.6124676116444*w11**2*u11*sin(1.16355185185185*x)*sin(0.581775925925926*y)**3*sin(0.581775925925926*x)*cos(0.581775925925926*x) ', ' 2453.79852094979*w11**2*u11*cos(1.16355185185185*x)*sin(0.581775925925926*y)*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 ', ' 5.02284144603642*w11**2*u11**2*sin(1.16355185185185*x)**2*sin(0.581775925925926*y)**2*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**2 - 76.6124676116444*w11**2*v11*sin(1.16355185185185*y)*sin(0.581775925925926*y)*cos(0.581775925925926*y)*sin(0.581775925925926*x)**3 ', ' 2453.79852094979*w11**2*v11*cos(1.16355185185185*y)*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)*cos(0.581775925925926*x)**2 ', ' 5.02284144603642*w11**2*v11**2*sin(1.16355185185185*y)**2*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)**2*cos(0.581775925925926*x)**2 ', ' 8.4273714372809*v11*u11**2*sin(1.16355185185185*x)**2*sin(1.16355185185185*y)*sin(0.581775925925926*y)*cos(0.581775925925926*y)*sin(0.581775925925926*x) ', ' 8.4273714372809*v11**2*u11*sin(1.16355185185185*x)*sin(1.16355185185185*y)**2*sin(0.581775925925926*y)*sin(0.581775925925926*x)*cos(0.581775925925926*x) ', ' 1169.64396165274*w11**2*u11*sin(1.16355185185185*x)*sin(0.581775925925926*y)*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)*cos(0.581775925925926*x) ', ' 1169.64396165274*w11**2*v11*sin(1.16355185185185*y)*sin(0.581775925925926*y)*cos(0.581775925925926*y)*sin(0.581775925925926*x)*cos(0.581775925925926*x)**2 ', ' 118.347813438402*w11**3*u11*sin(1.16355185185185*x)*sin(0.581775925925926*y)**2*cos(0.581775925925926*y)**2*sin(0.581775925925926*x)**2*cos(0.581775925925926*x) ', ' 118.347813438402*w11**3*v11*sin(1.16355185185185*y)*sin(0.581775925925926*y)**2*cos(0.581775925925926*y)*sin(0.581775925925926*x)**2*cos(0.581775925925926*x)**2 ', ' 416.57012675319*w11*v11*u11*sin(1.16355185185185*x)*cos(1.16355185185185*y)*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)*cos(0.581775925925926*x) ', ' 416.57012675319*w11*v11*u11*sin(1.16355185185185*y)*cos(1.16355185185185*x)*sin(0.581775925925926*y)*cos(0.581775925925926*y)*sin(0.581775925925926*x)**2 ', ' 0.852704767719347*w11*v11*u11**2*sin(1.16355185185185*x)**2*sin(1.16355185185185*y)*sin(0.581775925925926*y)**2*cos(0.581775925925926*y)*sin(0.581775925925926*x)**2 ', ' 0.852704767719347*w11*v11**2*u11*sin(1.16355185185185*x)*sin(1.16355185185185*y)**2*sin(0.581775925925926*y)**2*sin(0.581775925925926*x)**2*cos(0.581775925925926*x) ', ' 198.565093752355*w11*v11*u11*sin(1.16355185185185*x)*sin(1.16355185185185*y)*sin(0.581775925925926*y)*cos(0.581775925925926*y)*sin(0.581775925925926*x)*cos(0.581775925925926*x) ', ' 20.0913657841457*w11**2*v11*u11*sin(1.16355185185185*x)*sin(1.16355185185185*y)*sin(0.581775925925926*y)**2*cos(0.581775925925926*y)*sin(0.581775925925926*x)**2*cos(0.581775925925926*x)']
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

    def trapezoid(func, x0, x1, rtol = 1e-10, nseg0 = 1):
        """Интегрирование методом трапеций с заданной точностью.
           rtol - относительная точность,
           nseg0 - число отрезков начального разбиения"""
        ans = Quadrature.__restart(func, x0, x1, nseg0)
        old_ans = 0.0
        err_est = max(1, abs(ans))
        while (err_est > abs(rtol * ans)):
            old_ans = ans
            ans = Quadrature.__double_nseg(func, x0, x1)
            err_est = abs(old_ans - ans)

        print("Total function calls: " + str(Quadrature.__ncalls))
        return ans

    def simpson(func, x0, x1, rtol = 1.0e-10, nseg0 = 1):
        """Интегрирование методом парабол с заданной точностью.
           rtol - относительная точность,
           nseg0 - число отрезков начального разбиения"""
        old_trapez_sum = Quadrature.__restart(func, x0, x1, nseg0)
        new_trapez_sum = Quadrature.__double_nseg(func, x0, x1)
        ans = (4 * new_trapez_sum - old_trapez_sum) / 3
        old_ans = 0.0
        err_est = max(1, abs(ans))
        while (err_est > abs(rtol * ans)):
            old_ans = ans
            old_trapez_sum = new_trapez_sum
            new_trapez_sum = Quadrature.__double_nseg(func, x0, x1)
            ans = (4 * new_trapez_sum - old_trapez_sum) / 3
            err_est = abs(old_ans - ans)

        print("Total function calls: " + str(Quadrature.__ncalls))
        return ans

# zxc=str(X4(1)**2)
# print(zxc)
# subStrNew=str(integrate(X2(1)**2,(x, 0, b)))
# lenStrOld=len(zxc)
#
#
variable=[]
variable.append(str(X1(1)))
variable.append(str(integrate(X1(1),(x, 0, b))))
variable.append(str(X2(1)))
variable.append(str(integrate(X2(1),(x, 0, b))))
variable.append(str(X3(1)))
variable.append(str(integrate(X3(1),(x, 0, b))))
variable.append(str(X4(1)))
variable.append(str(integrate(X4(1),(x, 0, b))))
variable.append(str(X5(1)))
variable.append(str(integrate(X5(1),(x, 0, b))))
variable.append(str(X1(1)**2))
variable.append(str(integrate(X1(1)**2,(x, 0, b))))
variable.append(str(X2(1)**2))
variable.append(str(integrate(X2(1)**2,(x, 0, b))))
variable.append(str(X3(1)**2))
variable.append(str(integrate(X3(1)**2,(x, 0, b))))
variable.append(str(X4(1)**2))
variable.append(str(integrate(X4(1)**2,(x, 0, b))))
variable.append(str(X5(1)**2))
variable.append(str(integrate(X5(1)**2,(x, 0, b))))
variable.append(str(X1(1)**3))
variable.append(str(integrate(X1(1)**3,(x, 0, b))))
variable.append(str(X2(1)**3))
variable.append(str(integrate(X2(1)**3,(x, 0, b))))
variable.append(str(X3(1)**3))
variable.append(str(integrate(X3(1)**3,(x, 0, b))))
variable.append(str(X4(1)**3))
variable.append(str(integrate(X4(1)**3,(x, 0, b))))
variable.append(str(X5(1)**3))
variable.append(str(integrate(X5(1)**3,(x, 0, b))))
variable.append(str(X1(1)**4))
variable.append(str(integrate(X1(1)**4,(x, 0, b))))
variable.append(str(X2(1)**4))
variable.append(str(integrate(X2(1)**4,(x, 0, b))))
variable.append(str(X3(1)**4))
variable.append(str(integrate(X3(1)**4,(x, 0, b))))
variable.append(str(X4(1)**4))
variable.append(str(integrate(X4(1)**4,(x, 0, b))))
variable.append(str(X5(1)**4))
variable.append(str(integrate(X5(1)**4,(x, 0, b))))



variable.append(str(Y1(1)))
variable.append(str(integrate(Y1(1),(y, 0, b))))
variable.append(str(Y2(1)))
variable.append(str(integrate(Y2(1),(y, 0, b))))
variable.append(str(Y3(1)))
variable.append(str(integrate(Y3(1),(y, 0, b))))
variable.append(str(Y4(1)))
variable.append(str(integrate(Y4(1),(y, 0, b))))
variable.append(str(Y5(1)))
variable.append(str(integrate(Y5(1),(y, 0, b))))
variable.append(str(Y1(1)**2))
variable.append(str(integrate(Y1(1)**2,(y, 0, b))))
variable.append(str(Y2(1)**2))
variable.append(str(integrate(Y2(1)**2,(y, 0, b))))
variable.append(str(Y3(1)**2))
variable.append(str(integrate(Y3(1)**2,(y, 0, b))))
variable.append(str(Y4(1)**2))
variable.append(str(integrate(Y4(1)**2,(y, 0, b))))
variable.append(str(Y5(1)**2))
variable.append(str(integrate(Y5(1)**2,(y, 0, b))))
variable.append(str(Y1(1)**3))
variable.append(str(integrate(Y1(1)**3,(y, 0, b))))
variable.append(str(Y2(1)**3))
variable.append(str(integrate(Y2(1)**3,(y, 0, b))))
variable.append(str(Y3(1)**3))
variable.append(str(integrate(Y3(1)**3,(y, 0, b))))
variable.append(str(Y4(1)**3))
variable.append(str(integrate(Y4(1)**3,(y, 0, b))))
variable.append(str(Y5(1)**3))
variable.append(str(integrate(Y5(1)**3,(y, 0, b))))
variable.append(str(Y1(1)**4))
variable.append(str(integrate(Y1(1)**4,(y, 0, b))))
variable.append(str(Y2(1)**4))
variable.append(str(integrate(Y2(1)**4,(y, 0, b))))
variable.append(str(Y3(1)**4))
variable.append(str(integrate(Y3(1)**4,(y, 0, b))))
variable.append(str(Y4(1)**4))
variable.append(str(integrate(Y4(1)**4,(y, 0, b))))
variable.append(str(Y5(1)**4))
variable.append(str(integrate(Y5(1)**4,(y, 0, b))))

print(e[9])

nn=[]
number=0
for xc in e:
    for z in range(0,len(variable),2):
        while xc.find(variable[z]) > 0:
            subStrNew=variable[z+1]
            i = xc.find(variable[z])
            xc = xc[:i] + subStrNew + xc[i+len(variable[z]):]
            e[number]=xc
    number=number+1

print(e)
#
# bv=0
# for cv in e:
#     bv=bv+sm.symbols(cv)
# print(sm.expand(bv))

