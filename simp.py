#https://habr.com/ru/post/238703/
#https://ru.wikipedia.org/wiki/Cython
#cython --embed simp.py -o simp.c
#gcc -g -O2 -c -o simp.o simp.c `python-config --includes --ldflags`
import os, sys
import sympy as mp
import numpy as np
#sys.path.append(os.getcwd())

def simp(f,a,b,nn):
    n=2
#    if True:
    try:
        dx=((b-a)/n).evalf()
    except AttributeError:
        dx=((b-a)/n)
    s1=(f(a)+f(b)).evalf()
    s4=(f(a+dx)).evalf()
    s2=0.0
    s0=(dx/3.0*(s1+4.0*s4))
    while True:
        n*=2
        dx/=2
        s2=(s2+s4).evalf()
        s4=0.0
        for ix in range(1,n,2):
            s4=s4+f(a+dx*ix)
        s4=s4.evalf()
        s=(dx/3.0*(s1+4.0*s4+2.0*s2))
#        ds=abs(((s-s0)/s).evalf())
#        s0=s
#        print(f"{n}: dx={dx} ds={ds}",file=sys.stderr)
        if n>nn:
            break
    return s


