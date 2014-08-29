import unittest
import math
import numpy as np
import quadrature as quad 
from slater_test_f import slater_1s
import pdb

class QuadratureTest(unittest.TestCase):

    def test_euler_maclaurin(self):
        # Becke used this instead of Bragg-Slater value of 0.25A
        h_becke_rad = 0.35
        weights,roots = quad.euler_maclaurin_quad(0,1,1000)
        weights = weights[1:-2]
        roots = roots[1:-2]
        #print "\n".join(map(str,roots))
        sf = slater_1s(1,np.array((0,0,0)))
        # Change of variables to improve quadrature accuracy from Murray/Handy/Laming
        m = 2
        r_of_q = lambda q: h_becke_rad*pow(q,m)*pow(1 - q,-m)
        jacob_of_q = lambda q: h_becke_rad*m*pow(q,m-1)*pow(1 - q,(-m)-1)
        r_roots = map(r_of_q,roots)
        #print "\n".join(map(str,r_roots))
        jacobian_vals = map(jacob_of_q,roots)
        quad_point_fn = lambda i: pow(r_roots[i],2)*jacobian_vals[i]*pow(sf(np.array((0,0,r_roots[i]))),2)
        integral = 4*math.pi*sum(weights[i] * quad_point_fn(i) for i in range(0,len(weights)))
        self.assertLess(abs(1-integral),1e-14)

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