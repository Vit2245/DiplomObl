from unittest import TestCase
from sympy import Add, sin, cos, Mul, Pow, Symbol
from Obl_not_reb import recursive_cleaning


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

        self.assertEqual(test_expr, sin(x))
        self.assertEqual(test_expr2, sin(x))
        self.assertEqual(test_expr3, Mul(sin(x), x))
        self.assertEqual(test_expr4, cos(Mul(sin(x), x)))
        self.assertEqual(test_expr5, Pow(cos(Mul(sin(x), x)), 3))
        self.assertEqual(test_expr6, Pow(cos(Mul(sin(x), y)), 3))
        self.assertEqual(test_expr7, Pow(cos(Mul(sin(x), y)), 3))
