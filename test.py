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


x = sp.symbols('x')
y = sp.symbols('y')
q = sp.symbols('q')
def fun(a,b):
    return a+b

x=np.arange(0,5.4,5.4/30).astype(np.float64)
y=np.arange(0,5.4,5.4/30).astype(np.float64)
X,Y=np.meshgrid(x,y)
Z=fun(X,Y)

dat = np.random.random(200).reshape(20,10) # создаём матрицу значений


fig = plt.figure()
cf = plt.contourf(Z)
plt.xlabel("x,м")
plt.ylabel("y, ряд")

plt.colorbar(cf)
plt.title('Simple contourf plot')



