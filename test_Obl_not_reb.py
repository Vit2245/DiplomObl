from unittest import TestCase
from sympy import Add, sin, cos, Mul, Pow, Symbol
from Obl_not_reb import recursive_cleaning_internal


class Test(TestCase):
    def test_recursive_cleaning_internal(self):
        x = Symbol('x')
        test_expr = recursive_cleaning_internal(Add(3, sin(x)), x)
        test_expr2 = recursive_cleaning_internal(Add(3, Mul(sin(x), 5)), x)
        test_expr3 = recursive_cleaning_internal(Add(3, Mul(sin(x), x)), x)
        test_expr4 = recursive_cleaning_internal(cos(Add(3, Mul(sin(x), x))), x)
        test_expr5 = recursive_cleaning_internal(Pow(cos(Add(3, Mul(sin(x), x)))), x)

        self.assertEqual(test_expr, sin(x))
        self.assertEqual(test_expr2, sin(x))
        self.assertEqual(test_expr3, Mul(sin(x), x))
        self.assertEqual(test_expr4, cos(Mul(sin(x), x)))
        self.assertEqual(test_expr5, Pow(cos(Mul(sin(x), x))))
