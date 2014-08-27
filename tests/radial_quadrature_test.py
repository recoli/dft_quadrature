import unittest
import math
import numpy as np
from radial_quadrature import euler_maclaurin_quad
from slater_test_f import slater_1s
import pdb

class RadialQuadratureTest(unittest.TestCase):

    def test_euler_maclaurin(self):
        # Becke used this instead of Bragg-Slater value of 0.25A
        h_becke_rad = 0.35
        weights,roots = euler_maclaurin_quad(0,1,1000)
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
        self.assertLess(1-integral,1e-14)
