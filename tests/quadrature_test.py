import unittest
import math
import numpy as np
import quadrature as quad 
from slater_test_f import slater_1s
import itertools
import operator

class QuadratureTest(unittest.TestCase):

    def test_euler_maclaurin(self):
        # Define a piecewise function such that 4 points should give exact result of 5
        def func(x):
            if(x <= 1):
                return 2
            elif(x <= 2):
                return 1
            elif(x <= 3):
                return 2
            else:
                raise ValueError("Bad x value")

        weights,roots = quad.euler_maclaurin_quad(0,3,4)
        integral = sum(weights[i]*func(roots[i]) for i in range(0,len(roots))) 
        self.assertEqual(integral,5)

    def test_gauss_legendre_polynomial(self):
        # Polynomials degree 0-17 inclusive
        polys = [lambda x: 1,
                 lambda x: 5*x - 3,
                 lambda x: 16*pow(x,2) - 4*x + 3,
                 lambda x: 7*pow(x,3) + 2*pow(x,2) + 5,
                 lambda x: -6*pow(x,4) + 9*pow(x,3) - 2*pow(x,2) -7,
                 lambda x: 3*pow(x,5) - 7*pow(x,4) - 5*pow(x,3) + 4*pow(x,2) - x + 7,
                 lambda x: -15*pow(x,6) + 5*pow(x,5) + 2*pow(x,4) - 9*pow(x,3) - pow(x,2) + 8*x - 15,
                 lambda x: 121*pow(x,7)-18*pow(x,6) - 201*pow(x,4) + pow(x,3) + 71*pow(x,2) + 15*x - 2,
                 lambda x: 9*pow(x,8) - 6*pow(x,2),
                 ]
        degrees = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0] 
        results = [2,-6,50.0/3,34.0/3,-266.0/15,208.0/15,-3586.0/105,-4432.0/105,-2,]
        for i in range(0,len(polys)):
            n_points = max(int(math.ceil((degrees[i]+1.0)/2)),2)
            weights,roots = quad.gauss_legendre_quad(n_points)
            integral = sum(weights[m] * polys[i](roots[m]) for m in range(0,len(roots)))
            self.assertLess(abs(integral-(results[i])),1e-14)
