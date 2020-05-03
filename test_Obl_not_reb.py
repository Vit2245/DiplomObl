from unittest import TestCase

import matplotlib
from sympy import Add, sin, cos, Mul, Pow, Symbol, Basic, symbols, Function, Derivative
from Obl_not_reb import recursive_cleaning, show_approximates, create_plot, create_functional, show_approximate
import matplotlib.pyplot as plt


class Test(TestCase):
    def test_recursive_cleaning(self):
        x = Symbol('x')
        y = Symbol('y')
        test_expr = recursive_cleaning(Add(3, sin(x)), x)
        test_expr2 = recursive_cleaning(Add(3, Mul(sin(x), 5)), x)
        test_expr3 = recursive_cleaning(Add(3, Mul(sin(x), x)), x)
        test_expr4 = recursive_cleaning(cos(Add(3, Mul(sin(x), x))), x)
        test_expr5 = recursive_cleaning(Pow(cos(Add(3, Mul(sin(x), x))), 3), x)
        test_expr6 = recursive_cleaning(Pow(cos(Add(3, Mul(sin(x), y))), 3), x, y)
        test_expr7 = recursive_cleaning(Pow(cos(Add(3, Mul(sin(x), y))), 3), *{1: x, 2: y}.values())

        f = Function('f')(x)
        test_expr8 = recursive_cleaning(Derivative(Pow(cos(Add(3, Mul(sin(x), f))), 3), x), x)

        self.assertEqual(test_expr, sin(x))
        self.assertEqual(test_expr2, sin(x))
        self.assertEqual(test_expr3, Mul(sin(x), x))
        self.assertEqual(test_expr4, cos(Mul(sin(x), x)))
        self.assertEqual(test_expr5, Pow(cos(Mul(sin(x), x)), 3))
        self.assertEqual(test_expr6, Pow(cos(Mul(sin(x), y)), 3))
        self.assertEqual(test_expr7, Pow(cos(Mul(sin(x), y)), 3))
        self.assertEqual(test_expr8, Derivative(cos(f*sin(x) + 3)**3, x))

    def test_show_approximates(self):
        show_approximates()

    def test_create_plot(self):
        plot = create_plot(lambda x: x ** 2, plt)
        # plot.show()

    def test_create_functional(self):
        Es, SN, W = create_functional(1)

        self.assertIsInstance(Es, Basic)
        self.assertIsInstance(W, Basic)
        self.assertEqual(SN, list(symbols('u11, v11, w11, psix11, psiy11')))

    def test_show_approximate(self):
        show_approximate(3)
